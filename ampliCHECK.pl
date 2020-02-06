#!/usr/bin/perl -w
#
################################################################
#
# Name: ampliCHECK.pl
#
# Version: 1.0
#
# Author: Alvaro Sebastian
#
# Support: Alvaro Sebastian (bioquimicas@yahoo.es)
#
# License: GPL
#
# Evolutionary Biology Group
# Faculty of Biology
# Adam Mickiewicz University
#
# Description:
#   Fast pairwise comparison among higher depth sequences in amplicon sequencing (AS) data
#   Analyzes AS data and gives as output an Excel file where are shown the similarities among the high depth sequences
#   This results can be used to optimize the parameters for a further analysis with AmpliSAS
#
# Requires as input a FASTA or FASTQ file with sequences/reads and a CSV format file with primer/tag data information
# Example: perl ampliCHECK.pl -d amplicon_data.csv -i reads.fq.gz -o results
#
# If the reads have been already demultiplexed into separate files (one file per sample), they can be packed into a single .zip or .tar.gz file and use it as input
# Example: perl ampliCHECK.pl -i reads.tar.gz -o results
#
# Alternatively, variants already extracted can be given in an Excel format file obtained with AmpliSAS/AmpliCHECK.
# Example: perl ampliCHECK.pl -i amplisas_results.xlsx -o results
#


my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Calculates lengths, frequencies and putative errors for the most frequent variants from an amplicon sequencing experiment.";

# Modules are in folder 'lib' in the path of the script
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options
# Default minimum per amplicon frequency
my $MIN_FREQ = 1; # 5%
# Default read demultiplexing alignment method
my $INP_align = 'match';
# Default allele references alignment
my %INP_allele_align = ('alignment' => 'dna blastn -evalue 1E-5 -ungapped',
			'aligned' => 1,
			'ident' => 1,
			);
# Maximum number of allowed sequences per amplicon
my $INP_nreads_amplicon; # = 100000;
# Clustering parameters
# my @ALL_CLUSTERING_PARAMS = ('substitution_threshold', 'indel_threshold', 'cluster_exact_length', 'cluster_inframe', 'cluster_nonstop', 'identity_threshold', 'min_dominant_frequency_threshold', 'min_amplicon_seq_frequency_threshold');
my @ALL_CLUSTERING_PARAMS = ('substitution_threshold', 'indel_threshold', 'identity_threshold');
# Sequence filtering parameters
# my @ALL_FILTERING_PARAMS = ('allowed_markers', 'allowed_samples', 'min_samples', 'min_samples_to_keep', 'min_amplicon_seq_frequency_to_keep', 'min_cluster_size', 'min_amplicon_depth', 'min_amplicon_seq_depth', 'min_amplicon_seq_frequency', 'min_amplicon_seq_identity', 'min_dominant_seq_frequency', 'discard_frameshifts', 'max_amplicon_length_error', 'min_chimera_length', 'degree_of_change');
my @ALL_FILTERING_PARAMS = ('allowed_markers', 'allowed_samples', 'min_samples', 'min_amplicon_depth', 'min_amplicon_seq_depth', 'min_amplicon_seq_frequency');
# Default recommended technology parameters
my $TECH_PARAMS = {
	'pyroseq' => { # 454 and IonTorrent
		'substitution_threshold' => { 'all' => [ '0.5' ] },
		'indel_threshold' => { 'all' => [ '1' ] },
	},
	'illumina' => { # Illumina
		'substitution_threshold' => { 'all' => [ '1' ] },
		'indel_threshold' => { 'all' => [ '0.001' ] },
	},
	'unknown' => { # Unknown technology
		'substitution_threshold' => { 'all' => [ '1' ] },
		'indel_threshold' => { 'all' => [ '1' ] },
	}
};
$TECH_PARAMS->{'454/iontorrent'} = $TECH_PARAMS->{'454'} = $TECH_PARAMS->{'iontorrent'} = $TECH_PARAMS->{'pyroseq'};
$TECH_PARAMS->{'solexa'} = $TECH_PARAMS->{'illumina'};

my ($INP_amplicons_file, $INP_reads_file, $INP_outpath, $INP_auto, $INP_tech, $INP_filter_freq, $INP_cluster_ident, $INP_threads, $INP_nreads, $INP_direct, $INP_shuffle, $INP_expand_results, $INP_test, $INP_allele_file, $INP_zip);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s' => \$INP_reads_file,
	'd|data=s' => \$INP_amplicons_file,
	'auto' => \$INP_auto,
	'o|output=s' => \$INP_outpath,
	't|tech=s' => \$INP_tech,
	'fr|frequency=f' => \$INP_filter_freq,
	'ci|clid=f' => \$INP_cluster_ident,
	'a|alleles=s' => \$INP_allele_file,
	'e|expand' => \$INP_expand_results,
# 	's|shuffle' => \$INP_shuffle,
	'n|number=i' => \$INP_nreads,
	'na|max=i' => \$INP_nreads_amplicon,
	'di|direct' => \$INP_direct,
	'thr|threads=i' => \$INP_threads,
	'test' => \$INP_test,
	'z|zip' => \$INP_zip,
	'<>' => \&usage,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -i <file>\tInput FASTQ or FASTA file (compressed or uncompressed) or set of files packed into a unique .ZIP or .TAR.GZ file or AmpliSAS format Excel file.\n";
	print "  -d <file>\tCSV file with primer/amplicon data.\n";
	print "  -o <path>\tOutput folder name.\n";
	print "  -t <tech>\tUse recommended technology parameters ('454', 'IonTorrent', 'Illumina', 'Unknown').\n";
	print "  -fr <freq>\tCheck only sequences with higher frequency. (default=$MIN_FREQ%).\n";
	print "  -ci <number>\tCluster together sequences with higher or equal identity.\n";
	print "  -a <file>\tFASTA file with allele names and sequences.\n";
	print "  -e\t\tExpand results in 3 Excel sheets (depths, freqs and errors).\n";
# 	print "  -s\t\tShuffle/randomize reads/sequences to analyze.\n";
	print "  -n <number>\tNumber of reads/sequences to analyze.\n";
	print "  -na <number>\tNumber of reads/sequences per amplicon to analyze.\n";
	print "  -di\t\tAnalyze reads only in direct sense.\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	#  print " -auto\t\tParameters are automatically adjusted by a preliminar analysis of amplicons (DEFAULT).\n";
	#  print " -test\t\t Run in 'test' mode and create a dump file of the alignments.\n";
	exit;
}

# Checks if an excel file or multifle are given as input with sequences previously processed
my ($INP_excel_file, $INP_multifile);
if (defined($INP_reads_file) && is_xlsx($INP_reads_file)){
	$INP_excel_file = $INP_reads_file;
# Checks if a set of demultiplexed files is given as input into a compressed file
} elsif (defined($INP_reads_file) && is_multifile($INP_reads_file)){
	$INP_multifile = $INP_reads_file;
# Prints usage help if no input file is specified
} elsif (!defined($INP_reads_file) || !-f $INP_reads_file){
	print "\nERROR: You must specify a sequence input file.\n\n";
	usage();
	exit;
} elsif (!defined($INP_amplicons_file) || !-f $INP_amplicons_file){
	print "\nERROR: You must specify an amplicon data input file.\n\n";
	usage();
	exit;
}
if (defined($INP_auto) && defined($INP_tech)){
	print "\nERROR: 'Auto mode' and 'Technology' options are incompatible.\n\n".
	"'Auto mode' will automatically detect the technology used and adjust\n".
	"the sequencing error rates and other analysis parameters.\n\n";
	usage();
	exit;
}
# Checks technology default parameters
if (defined($INP_tech)){
	$INP_tech =~ s/\s//g;
	$INP_tech = lc($INP_tech);
	if (!defined($TECH_PARAMS->{$INP_tech})){ 
		print "\nERROR: You must specify a valid sequencing technology ('454', 'IonTorrent', 'Illumina' or 'Unknown').\n\n";
		usage();
		exit;
	}
}
# Checks percentage threshold parameters, if defined
my %perc_params = (	'Filtering frequency percentage'=>$INP_filter_freq,
			'Clustering identity percentage'=>$INP_cluster_ident
			);
while (my ($param, $value) = each %perc_params) {
	if (defined($value)){
		if ($value =~ /([\d\.]+)%/) {
			$value = $1;
		}
		if (!is_numeric($value) || $value<0 || $value>100){
			printf("\nERROR: %s must be a number between 0 and 100.\n\n",$param);
			exit;
		}
	}
}
# Default output path
if (!defined($INP_outpath)){
	$INP_outpath  = lc((split('\.',$SCRIPT_NAME))[0]);
}

print "\nRunning '$COMMAND_LINE'\n";

# Assign dump file, only used in case of -test parameter chosen
my $dump_file = "$INP_outpath.ampli.dump";

# 2 PRIMERS => MARKER
# 1/2 TAGS => SAMPLE
# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

# my ($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads);
# Check and read sequences file
my ($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads,$ngs_tech);
if (!defined($INP_excel_file) && !defined($INP_multifile) && (!defined($INP_test) || !-e $dump_file)){
	($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads)
	= parse_sequence_file($INP_reads_file,$INP_nreads,['verbose']);
	# Detects the sequencing technology used, if reads file comes from the instrument
	$ngs_tech = ngs_tech($read_headers->[0]);
}
# # Randomizes/shuffles reads (if $INP_nreads is defined, seqs will be already random)
# if (defined($INP_shuffle) && !defined($INP_nreads)) {
# 	($read_headers,$read_seqs,$read_qualities) = shuffle_seqs($read_headers,$read_seqs,$read_qualities);
# }

# Check and read amplicons file
# Can be missing if an Excel file is given as input
my ($markerdata,$markers,$sampledata,$tags,$paramsdata,$alleledata);
if (!defined($INP_excel_file) && !defined($INP_multifile)){
	($markerdata,$markers,$sampledata,$tags,$paramsdata,$alleledata)
	= parse_amplicon_file($INP_amplicons_file,['verbose', 'skip samples']);
} elsif (defined($INP_amplicons_file)){
	($markerdata,$markers,$sampledata,$tags,$paramsdata,$alleledata)
	= parse_amplicon_file($INP_amplicons_file,['verbose','skip markers','skip samples']);
}

# Check and read alleles file (optional)
# Amplicon data file has preference over alleles in FASTA file
if (defined($INP_allele_file) && -e $INP_allele_file){
	print "\nReading allele sequences from '$INP_allele_file'.\n";
	$alleledata = read_allele_file($INP_allele_file);
}

# Retrieves clustering and filtering parameters from the CSV file
# If not defined and technology is specified, technogy defaults are used
# If there are not parameters in CSV file and no tech chosen, automatic mode is activated
my @CLUSTERING_PARAMS;
my @FILTERING_PARAMS;
if (defined($INP_auto) || (!defined($INP_tech) && !defined($paramsdata))){
	$INP_auto = 1;
	print "\n'AUTO MODE' activated.\n".
		"\tThe program will automatically detect the sequencing technology used\n".
		"\tand adjust the sequencing error rates and other analysis parameters.\n";
}

# Sets minimum per amplicon frequency if it's specified in command line
# Command line params have priority over CSV file and auto mode
if (defined($INP_filter_freq)){
	$paramsdata->{'min_amplicon_seq_frequency'}{'all'} = [ $INP_filter_freq ];
} elsif (!defined($paramsdata->{'min_amplicon_seq_frequency'})){
	print "\nNo minimum per amplicon frequency specified, setting it automatically to $MIN_FREQ%.\n";
	$paramsdata->{'min_amplicon_seq_frequency'}{'all'} = [ $MIN_FREQ ];
}
# Sets minimum identity clustering threshold if it's specified in command line
if (defined($INP_cluster_ident)){
	$paramsdata->{'identity_threshold'}{'all'} = [ $INP_cluster_ident ];
}
# Read always the params from CSV file, they will have priority over auto ones
foreach my $paramname (@ALL_CLUSTERING_PARAMS){
	if (defined($paramsdata->{$paramname})) {
		push(@CLUSTERING_PARAMS,$paramname);
	} elsif (defined($INP_tech) && defined($TECH_PARAMS->{$INP_tech}{$paramname})){
		$paramsdata->{$paramname} = $TECH_PARAMS->{$INP_tech}{$paramname};
		push(@CLUSTERING_PARAMS,$paramname);
	}
}
foreach my $paramname (@ALL_FILTERING_PARAMS){
	if (defined($paramsdata->{$paramname})) {
		push(@FILTERING_PARAMS,$paramname);
	} elsif (defined($INP_tech) && defined($TECH_PARAMS->{$INP_tech}{$paramname})){
		$paramsdata->{$paramname} = $TECH_PARAMS->{$INP_tech}{$paramname};
		push(@FILTERING_PARAMS,$paramname);
	}
}
if (!defined($INP_auto) && !@CLUSTERING_PARAMS) {
	print "\nERROR: No clustering parameters are specified, use auto mode to adjust them automatically. If you want to skip clustering use the option 'nocluster', otherwise use default parameters by specifying the sequencing technology used or include them into CSV data file.\n\n";
	usage();
	exit;
}
if (!defined($INP_auto) && !@FILTERING_PARAMS && !defined($INP_filter_freq)) {
	print "\nERROR: No filtering parameters are specified, use auto mode to adjust them automatically. If you want to skip filtering use the option 'nofilter', otherwise use default parameters by specifying the sequencing technology used or include them into CSV data file.\n\n";
	usage();
	exit;
}

# 	# Removes clustering and filtering params not used by AmpliCHECK (inherited from AmpliSAS)
# 	# Ex. IF WE INCLUDE CHIMERAS FILTER, MANY SEQUENCES WILL BE EXCLUDED FROM THE ANALYSIS
# 	foreach my $paramname (keys %{$paramsdata}){
# 		if (!in_array([@ALL_CLUSTERING_PARAMS, @ALL_FILTERING_PARAMS],$paramname)){
# 			delete($paramsdata->{$paramname});
# 		}
# 	}
# 	foreach my $paramname (@ALL_CLUSTERING_PARAMS){
# 		if (defined($paramsdata->{$paramname})) {
# 			push(@CLUSTERING_PARAMS,$paramname);
# 		} elsif (defined($INP_tech) && defined($TECH_PARAMS->{$INP_tech}{$paramname})){
# 			$paramsdata->{$paramname} = $TECH_PARAMS->{$INP_tech}{$paramname};
# 			push(@CLUSTERING_PARAMS,$paramname);
# 		}
# 	}
# 	foreach my $paramname (@ALL_FILTERING_PARAMS){
# 		# Sets minimum frequency if it's specified in command line
# 		if ($paramname eq 'min_amplicon_seq_frequency' && defined($INP_filter_freq)){
# 			$paramsdata->{$paramname} = { 'all' => [ $INP_filter_freq ] };
# 			push(@FILTERING_PARAMS,$paramname);
# 		} elsif (defined($paramsdata->{$paramname})) {
# 			push(@FILTERING_PARAMS,$paramname);
# 		# If minimum frequency and minimum depth are not specified, includes default min frequency
# 		} elsif ($paramname eq 'min_amplicon_seq_frequency' && !defined($paramsdata->{'min_amplicon_seq_depth'})){
# 			$paramsdata->{$paramname} = { 'all' => [ $MIN_FREQ ] };
# 			push(@FILTERING_PARAMS,$paramname);
# 		} elsif (defined($INP_tech) && defined($TECH_PARAMS->{$INP_tech}{$paramname})){
# 			$paramsdata->{$paramname} = $TECH_PARAMS->{$INP_tech}{$paramname};
# 			push(@FILTERING_PARAMS,$paramname);
# 		}
# 	}
# 	if (!@CLUSTERING_PARAMS || !@FILTERING_PARAMS) {
# 		print "\nERROR: Not enough comparison parameters are specified. Use default parameters by specifying the sequencing technology used or include them into CSV data file.\n\n";
# 		usage();
# 		exit;
# 	}
# 

# After previous checks to confirm that there are no errors in input data
# Creates output folders
my $INP_outpath_allseqs = "$INP_outpath/allseqs";
if (!-d $INP_outpath){
	mkdir($INP_outpath);
}
if (!defined($INP_excel_file) && !defined($INP_multifile) && !-d $INP_outpath_allseqs){
	mkdir($INP_outpath_allseqs);
}

# Print amplicon parameters into a file
print "\nPrinting amplicon data into '$INP_outpath/amplicon_data.csv'.\n";
# my $amplicon_data = ">command\n$COMMAND_LINE\n\n";
my $amplicon_data = '';
$amplicon_data .= print_amplicon_data($markerdata,'markers')."\n";
$amplicon_data .= print_amplicon_data($sampledata,'samples')."\n";
if (defined($paramsdata) && %$paramsdata){
	$amplicon_data .= print_amplicon_data($paramsdata,'params',\@CLUSTERING_PARAMS,\@FILTERING_PARAMS)."\n";
} # elsif (!defined($INP_nocluster) && defined($INP_nofilter)){
# 		$amplicon_data .= print_amplicon_data($paramsdata,'params',\@CLUSTERING_PARAMS)."\n";
# 	} elsif (defined($INP_nocluster) && !defined($INP_nofilter)){
# 		$amplicon_data .= print_amplicon_data($paramsdata,'params',\@FILTERING_PARAMS)."\n";
# 	}
if (defined($alleledata)){
	$amplicon_data .= print_amplicon_data($alleledata,'alleles')."\n";
}
write_to_file("$INP_outpath/amplicon_data.csv",$amplicon_data);

# Extracts amplicon sequences and depths
# $amplicon_raw_sequences stores the individual unique sequence depths
# $amplicon_raw_sequences->{$marker_name}{$sample_name}{$md5} = $depth;
# $amplicon_raw_depths stores the total depth of the sequences into an amplicon
# $amplicon_raw_depths->{$marker_name}{$sample_name}++;
my ($md5_to_sequence,$amplicon_raw_sequences,$amplicon_raw_depths);
my ($marker_seq_data, $amplicon_seq_data);
my $samples;
my $md5_to_name;
my ($marker_result_file,$marker_seq_files,$marker_matrix_files,$amplicon_seq_files);
if (!defined($INP_excel_file) && !defined($INP_multifile)){
	if (!defined($INP_test) || !-e $dump_file){
		if ($INP_align eq 'match') {
			print "\nDe-multiplexing amplicon sequences from reads.\n";
			# Creates a file with all sequences
			my $raw_seqs_file = write_to_file("/tmp/".random_file_name(),join("\n",@{$read_seqs}));
			my $match_options;
			if (defined($INP_direct)) {
				push(@$match_options, 'direct');
			}
			# Parse reads and primers+tags with GAWK to find matching amplicons
			# Is equivalent to 'align_amplicons+match_amplicons' with perfect matching, but very fast
			if (defined($INP_threads) && $INP_threads>1){
				($md5_to_sequence,$amplicon_raw_sequences,$amplicon_raw_depths)
				= match_amplicons_regex_with_threads($raw_seqs_file,$markerdata,$sampledata,$markers,$tags,$INP_nreads_amplicon,$match_options,$INP_threads);
			} else {
				($md5_to_sequence,$amplicon_raw_sequences,$amplicon_raw_depths)
				= match_amplicons_regex($raw_seqs_file,$markerdata,$sampledata,$markers,$tags,$INP_nreads_amplicon,$match_options);
			}
			if (defined($INP_test)) {
				print "\nDumping alignment data into '$dump_file'.\n";
				store_data_dump([$md5_to_sequence,$amplicon_raw_sequences,$amplicon_raw_depths], $dump_file, 'Storable');
			}
			`rm $raw_seqs_file`;
		} # else {
# 			# Aligns reads against primer/tag sequences
# 			print "\nAligning primer/tag sequences.\n";
# 			my $align_amplicon_data
# 			= align_amplicons($read_headers,$read_seqs,$primer_tag_headers,$primer_tag_seqs,$INP_align,$INP_revcomp,$INP_threads,$INP_test,$INP_outpath);
# 			print "\nDe-multiplexing amplicon sequences from reads.\n";
# 			# Loop reads with alignment results and extract amplicons (each amplicon is identified by a MD5 hash)
# 			($md5_to_sequence,$amplicon_raw_sequences,$amplicon_raw_depths)
# 			= match_amplicons($read_headers,$read_seqs,$markerdata,$sampledata,$align_amplicon_data,$INP_nreads_amplicon);
# 			# undef($align_amplicon_data);
# 		}
		# Free memory
		undef($read_headers);
		undef($read_seqs);
		undef($read_qualities);
	} else {
		print "\nRecovering de-multiplexed amplicon sequences from '$dump_file'.\n";
		($md5_to_sequence,$amplicon_raw_sequences,$amplicon_raw_depths) = @{recover_data_dump($dump_file, 'Storable')};
	}

	# Assigns samples/tags to markers

	foreach my $marker_name (@$markers){
		foreach my $sample_name (@$tags) {
			if (defined($amplicon_raw_depths->{$marker_name}{$sample_name})){
				push(@{$samples->{$marker_name}}, $sample_name);
			}
		}
	}

	# Align sequences to alleles and assign allele names to sequences
	if (defined($alleledata) && %$alleledata){
		print "\nMatching allele sequences.\n";
		$md5_to_name = match_alleles($alleledata,$md5_to_sequence,$md5_to_name,\%INP_allele_align,$INP_threads);
	}

	# Extracts marker/amplicon sequence data
	# $marker_seq_data stores the names and parameters of all unique sequences of a unique marker
	# $marker_seq_data->{$marker_name}{$md5} = { 'seq'=> $seq, 'name'=>$name, 'len'=>$len, 'depth'=>$unique_seq_depth, 'samples'=>$count_samples, 'mean_freq'=>$mean_freq, 'max_freq'=>$max_freq, 'min_freq'=>$min_freq,};
	# $amplicon_seq_data stores the names and parameters of all unique sequences of a unique amplicon
	# $amplicon_seq_data->{$marker_name}{$sample_name}{$md5} = { 'seq'=> $seq, 'name'=>$name, 'len'=>$len, 'depth'=>$unique_seq_depth, 'freq'=>$unique_seq_frequency, 'cluster_size'=>$cluster_size };
	print "\nExtracting de-multiplexed sequences into '$INP_outpath_allseqs'.\n";
	($marker_seq_data, $amplicon_seq_data, $md5_to_name)
	= retrieve_amplicon_data($markers,$samples,$amplicon_raw_sequences,$amplicon_raw_depths,$md5_to_sequence,$md5_to_name);
	# Prints marker sequences
	($marker_result_file,$marker_seq_files,$marker_matrix_files)
	= print_marker_sequences($markers,$samples,$marker_seq_data,$amplicon_seq_data,$amplicon_raw_depths,$INP_outpath_allseqs); # ,$STC_seq_data,$amplicon_raw_sequences);
	# Copy Excel results file to output parent folder
	# `cp $marker_result_file $INP_outpath/allseqs.xlsx`;
	# Prints amplicon sequences
	$amplicon_seq_files
	= print_amplicon_sequences($markers,$samples,$marker_seq_data,$amplicon_seq_data,$INP_outpath_allseqs);
	# Prints all amplicon sequences into an unique file
	# system("cat ".join(" ",@{$amplicon_seq_files}."> $INP_outpath/allseqs.fasta");

} elsif (defined($INP_multifile)) {

	printf("\nReading File '%s'.\n", $INP_multifile);
	my ($markers_, $samples_);
	($markers_, $samples_, $amplicon_raw_sequences, $amplicon_raw_depths, $md5_to_sequence)
	= read_compressed_file_amplicons($INP_multifile,[],$INP_nreads_amplicon);

	# Uses the markers specified by the multifile folders
	$markers = $markers_;
	
	# Assigns samples/tags to markers
	foreach my $marker_name (@$markers){
		foreach my $sample_name (@{$samples_->{$marker_name}}) {
			if (defined($amplicon_raw_depths->{$marker_name}{$sample_name})){
				push(@{$samples->{$marker_name}}, $sample_name);
			} elsif (is_numeric($sample_name) && defined($amplicon_raw_depths->{$marker_name}{sprintf('%d',$sample_name)})){
				# For example sample '001' could be writen into Excel as '1'
				push(@{$samples->{$marker_name}}, sprintf('%d',$sample_name));
			}
		}
	}

	# Recalculates marker and amplicon data
	($marker_seq_data, $amplicon_seq_data, $md5_to_name)
	= retrieve_amplicon_data($markers,$samples,$amplicon_raw_sequences,$amplicon_raw_depths,$md5_to_sequence,$md5_to_name);

	# Align sequences to alleles and assign allele names to sequences
	if (defined($alleledata) && %$alleledata){
		print "\nMatching allele sequences.\n";
		$md5_to_name = match_alleles($alleledata,$md5_to_sequence,$md5_to_name,\%INP_allele_align,$INP_threads);
	}

} elsif (defined($INP_excel_file)) {

	printf("\nReading File '%s'.\n", $INP_excel_file);
	my ($markers_, $samples_);
	($markers_, $samples_, $amplicon_raw_sequences, $amplicon_raw_depths, $md5_to_sequence)
	= read_amplisas_file_amplicons($INP_excel_file,[],$INP_nreads_amplicon);

	# Uses the markers specified in the Excel file
	$markers = $markers_;
	
	# Assigns samples/tags to markers
	foreach my $marker_name (@$markers){
		foreach my $sample_name (@{$samples_->{$marker_name}}) {
			if (defined($amplicon_raw_depths->{$marker_name}{$sample_name})){
				push(@{$samples->{$marker_name}}, $sample_name);
			} elsif (is_numeric($sample_name) && defined($amplicon_raw_depths->{$marker_name}{sprintf('%d',$sample_name)})){
				# For example sample '001' could be writen into Excel as '1'
				push(@{$samples->{$marker_name}}, sprintf('%d',$sample_name));
			}
		}
	}

	# Recalculates marker and amplicon data
	($marker_seq_data, $amplicon_seq_data, $md5_to_name)
	= retrieve_amplicon_data($markers,$samples,$amplicon_raw_sequences,$amplicon_raw_depths,$md5_to_sequence,$md5_to_name);

	# Align sequences to alleles and assign allele names to sequences
	if (defined($alleledata) && %$alleledata){
		print "\nMatching allele sequences.\n";
		$md5_to_name = match_alleles($alleledata,$md5_to_sequence,$md5_to_name,\%INP_allele_align,$INP_threads);
	}
}

# Detects automatically the sequencing technoloy
# And adjusts analysis parameters
if (defined($INP_auto)){
	my ($tech_auto, $lengths_auto, $params_auto);
	print "\nChecking data and setting automatically clustering parameters.\n";
	($tech_auto, $lengths_auto, $params_auto) = adjust_automatic_params($markers,$samples,$marker_seq_data,$amplicon_seq_data);
	if (defined($ngs_tech)){
		printf("\tTechnology: %s\n", ucfirst($ngs_tech));
	} elsif (defined($tech_auto)){
		printf("\tTechnology: %s\n", ucfirst($tech_auto));
	}
	printf("\tError rates: %s%% substitutions, %s%% indels\n",$params_auto->{'substitution_threshold'}{'all'}[0] ,$params_auto->{'indel_threshold'}{'all'}[0]);
	# Sets clustering parameters with the values retrieved by 'adjust_automatic_params'
	# 'indel_threshold', 'substitution_threshold', 'cluster_inframe', 'min_dominant_frequency_threshold'
	# Unless the parameter has been adjusted before on command line or CSV file
	foreach my $paramname (keys %$params_auto){
		push(@CLUSTERING_PARAMS,$paramname);
		foreach my $marker_name (keys %{$params_auto->{$paramname}}){
			if (!defined($paramsdata->{$paramname}{'all'}) && !defined($paramsdata->{$paramname}{$marker_name})){
				$paramsdata->{$paramname}{$marker_name} = $params_auto->{$paramname}{$marker_name};
			}
		}
	}
}

print "\nChecking data and setting marker lengths.\n";
# Checks if marker lengths are defined, if not they will be calculated automatically if they are required
foreach my $marker_name (@$markers){
	my $length_type = 'manual';
	# Calculates automatically marker length in case of not being specified
	#if (defined($paramsdata->{'cluster_exact_length'}{'all'}) || defined($paramsdata->{'cluster_inframe'}{'all'}) || defined($paramsdata->{'cluster_exact_length'}{$marker_name}) || defined($paramsdata->{'cluster_inframe'}{$marker_name})){
	if (!defined($markerdata->{$marker_name}{'length'})){
		$markerdata->{$marker_name}{'length'} = adjust_automatic_params([$marker_name],$samples,$marker_seq_data,$amplicon_seq_data,['only lengths'])->{$marker_name};
		$length_type = 'auto';
	}
	if ($#{$markerdata->{$marker_name}{'length'}} == 0){
		printf("\tMarker '%s' length: %s (%s)\n",$marker_name,$markerdata->{$marker_name}{'length'}[0],$length_type);
	} else {
		printf("\tMarker '%s' lengths: %s (%s)\n",$marker_name,join(',',@{$markerdata->{$marker_name}{'length'}}),$length_type);
	}
}

# Saves data after de-multiplexing for filtering purposes later
# Filter markers and samples with low coverage or not desired
foreach my $marker_name (keys %{$amplicon_seq_data}){
	# Skips markers without data
	if (!defined($amplicon_seq_data->{$marker_name})){
		delete($amplicon_seq_data->{$marker_name});
		next;
	}
	# Exclude amplicons not in the list
	if (defined($paramsdata->{'allowed_markers'}) && !defined($paramsdata->{'allowed_markers'}{'all'}) && !in_array($paramsdata->{'allowed_markers'},$marker_name)){
		delete($amplicon_seq_data->{$marker_name});
		next;
	}
	foreach my $sample_name (keys %{$amplicon_seq_data->{$marker_name}}){
		# Exclude samples not in the list
		if (defined($paramsdata->{'allowed_samples'}) && !defined($paramsdata->{'allowed_samples'}{'all'}) && !defined($paramsdata->{'allowed_samples'}{$marker_name})){
			delete($amplicon_seq_data->{$marker_name}{$sample_name});
			next;
		}
		if (defined($paramsdata->{'allowed_samples'}) && !defined($paramsdata->{'allowed_samples'}{'all'}) && ref($paramsdata->{'allowed_samples'}{$marker_name}) eq 'ARRAY' && !in_array($paramsdata->{'allowed_samples'}{$marker_name},$sample_name) ){
			delete($amplicon_seq_data->{$marker_name}{$sample_name});
			next;
		}
		# Exclude samples with low coverage
		my $total_seqs = $amplicon_raw_depths->{$marker_name}{$sample_name};
		if (defined($paramsdata->{'min_amplicon_depth'}) && defined($paramsdata->{'min_amplicon_depth'}{'all'}) && (!defined($total_seqs) || $total_seqs<$paramsdata->{'min_amplicon_depth'}{'all'}[0])){
			delete($amplicon_seq_data->{$marker_name}{$sample_name});
			next;
		} elsif (defined($paramsdata->{'min_amplicon_depth'}) && defined($paramsdata->{'min_amplicon_depth'}{$marker_name}) && (!defined($total_seqs) || $total_seqs<$paramsdata->{'min_amplicon_depth'}{$marker_name}[0])){
			delete($amplicon_seq_data->{$marker_name}{$sample_name});
			next;
		}
	}
}

# Apply filters to amplicon sequences
my ($amplicon_filtered_sequences, $amplicon_filtered_depths,$filters_output);
print "\nFiltering sequences with the following criteria ('filter' 'marker' 'values'):\n";
# In $paramsdata hash there are also clustering thresholds
foreach my $filtering_param (@FILTERING_PARAMS){
	if ($filtering_param ne 'allowed_markers' && defined($paramsdata->{$filtering_param})) {
		foreach my $marker_name (sort keys %{$paramsdata->{$filtering_param}}){
			print "\t$filtering_param\t$marker_name\t".join(',',@{$paramsdata->{$filtering_param}{$marker_name}})."\n";
		}
	} elsif (defined($paramsdata->{$filtering_param})) {
		print "\t$filtering_param\t".join(',',@{$paramsdata->{$filtering_param}})."\n";
	}
}
print "\n";
if (defined($INP_threads) && $INP_threads>1){
	($amplicon_filtered_sequences,$amplicon_filtered_depths,$filters_output)
	= filter_amplicon_sequences_with_threads($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$amplicon_raw_depths,$INP_threads);
} else {
	($amplicon_filtered_sequences,$amplicon_filtered_depths,$filters_output)
	= filter_amplicon_sequences($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$amplicon_raw_depths);
}

if (!defined($amplicon_filtered_sequences)){
	print "\nNo sequences passed the filters, please choose broader filtering values.\n\n";
	exit;
} else {
	# Free memory
	undef($marker_seq_data);
	undef($amplicon_seq_data);
	# Frequencies will be calculated with original depths, not after clustering+filtering ones
	($marker_seq_data, $amplicon_seq_data, $md5_to_name)
	= retrieve_amplicon_data($markers,$samples,$amplicon_filtered_sequences,$amplicon_raw_depths,$md5_to_sequence,$md5_to_name);
}

# Compares sequences with amplicon frequencies higher than threshold
# Retrieves similarities among sequences
# $amplicon_compared_sequences->{$marker_name}{$sample_name}{$dominant_md5}{$variant_md5} = { 'ident' => $identity, 'sub' => scalar @$substitutions, 'ins' => scalar @$insertions, 'del' => scalar @$deletions, 'homo_indel' => scalar @$homopolymer_indels };
print "\nComparing amplicon sequences that passed the filters with the following parameters ('threshold' 'marker' 'values'):\n";
# In $paramsdata hash there are also sequence filters
foreach my $clustering_param (@CLUSTERING_PARAMS){
	if (defined($paramsdata->{$clustering_param})) {
		foreach my $marker_name (sort keys %{$paramsdata->{$clustering_param}}){
			print "\t$clustering_param\t$marker_name\t".join(',',@{$paramsdata->{$clustering_param}{$marker_name}})."\n";
		}
	}
}
print "\n";
my $amplicon_seq_comparison_data;
if (defined($INP_threads) && $INP_threads>1){
	$amplicon_seq_comparison_data
	= compare_amplicon_sequences_with_threads($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$INP_threads);
} else {
	$amplicon_seq_comparison_data
	= compare_amplicon_sequences($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data);
}

my $print_options;
if (defined($INP_expand_results)){
	$print_options = ['expand results'];
}
# Prints marker sequences
my $comparison_result_file
= print_comparison_sequences($markers,$samples,$marker_seq_data,$amplicon_seq_comparison_data,$amplicon_raw_depths,"$INP_outpath/results.xlsx",$print_options);

# Print number of unique sequences per amplicon
my $summary_output = "Amplicon\tTotal\tUnique";
$summary_output .= "\n";
foreach my $marker_name (@$markers){
	foreach my $sample_name (@{$samples->{$marker_name}}) {
		$summary_output .= sprintf("%s-%s",$marker_name,$sample_name);
		if (defined($amplicon_raw_depths->{$marker_name}{$sample_name})){
			$summary_output .= sprintf("\t%d\t%d",$amplicon_raw_depths->{$marker_name}{$sample_name},scalar keys %{$amplicon_raw_sequences->{$marker_name}{$sample_name}});
		} else {
			$summary_output .= sprintf("\t%d\t%d",0,0);
		}
		$summary_output .= "\n";
	}
}
print "\nSequences per amplicon:\n$summary_output\n";

if (defined($INP_zip) && -d $INP_outpath && !is_folder_empty($INP_outpath)){
	my $cwd = getcwd;
	chdir($INP_outpath);
	my $outname = basename($INP_outpath);
	`zip -qrm $outname.zip *` ;
	`mv $outname.zip ..`;
	chdir($cwd);
	rmdir($INP_outpath);
	print "\nAnalysis results stored into '$INP_outpath.zip'.\n\n";
} elsif (-d $INP_outpath && !is_folder_empty($INP_outpath)){
	print "\nAnalysis results stored into '$INP_outpath'.\n\n";
} else {
	print "\nThere was some error in the analysis and no results were retrieved.\n\n";
}


#################################################################################

exit;
