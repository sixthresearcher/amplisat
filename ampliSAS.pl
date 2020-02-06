#!/usr/bin/perl -w
#
# Name: ampliSAS.pl
#
# Version: 1.0
#
# License: Creative Commons Attribution-ShareAlike 4.0 International
#          CC BY-NC-SA 4.0 (http://creativecommons.org/licenses/by-nc-sa/4.0/)
#
# Author: Alvaro Sebastian
#
# Support: Alvaro Sebastian (bioquimicas@yahoo.es)
#
# Evolutionary Biology Group
# Faculty of Biology
# Adam Mickiewicz University
#
# Description:
# Extracts and analyses amplicon sequencing (AS) data.
#
# Requires as input a FASTA or FASTQ file with sequences/reads and a CSV format file with primer/tag data information
# Example: perl ampliSAS.pl -d amplicon_data.csv -i reads.fq.gz -o results
#
# If the reads have been already demultiplexed into separate files (one file per sample), they can be packed into a single .zip or .tar.gz file and use it as input
# Example: perl ampliSAS.pl -i reads.tar.gz -o results
#
# Alternatively, variants already extracted can be given in an Excel format file obtained with AmpliSAS/AmpliCHECK.
# Example: perl ampliSAS.pl -i amplisas_results.xlsx -o results
#


my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Performs genotyping of amplicon sequencing data by clustering errors and filtering artefacts following a 3-step based analysis: de-multiplexing, clustering and filtering.";

# Modules are in folder 'lib' in the path of the script
use lib "lib";
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;
use Data::Dumper;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options
# Default read demultiplexing alignment method
my $INP_align = 'match';
# Default allele references alignment
my %INP_allele_align_params = ('
	alignment' => 'dna blastn -evalue 1E-5 -ungapped',
	'aligned' => 100,
	'ident' => 100,
);
# By default align direct and reverse complementary reads
my $INP_revcomp = 1;
# Maximum number of allowed sequences per amplicon
my $INP_nreads_amplicon;
# Amplicons with lower total coverage will be discarded
my $INP_min_amplicon_depth = 100;
# Clustering parameters
my @ALL_CLUSTERING_PARAMS = ('substitution_threshold', 'indel_threshold', 'cluster_exact_length', 'cluster_inframe', 'cluster_nonstop', 'identity_threshold', 'min_dominant_frequency_threshold', 'min_amplicon_seq_frequency_threshold', 'keep_singletons');
# Sequence filtering parameters
my @ALL_FILTERING_PARAMS = ('allowed_markers', 'allowed_samples', 'min_samples', 'min_samples_to_keep', 'min_amplicon_seq_frequency_to_keep', 'min_cluster_size', 'min_amplicon_depth', 'min_amplicon_seq_depth', 'min_amplicon_seq_frequency', 'min_amplicon_seq_identity', 'min_dominant_seq_frequency', 'discard_frameshifts', 'discard_noncoding', 'max_amplicon_length_error', 'min_chimera_length', 'degree_of_change', 'max_allele_number');
# Default recommended technology parameters
my $TECH_PARAMS = {
	'pyroseq' => { # 454 and IonTorrent
		'substitution_threshold' => { 'all' => [ '0.5' ] },
		'indel_threshold' => { 'all' => [ '1' ] },
		'min_dominant_frequency_threshold' => { 'all' => [ '25' ] },
		'cluster_inframe' => { 'all' => [ '1' ] },
		'min_amplicon_seq_frequency' => { 'all' => [ '3' ] },
		'min_chimera_length' => { 'all' => [ '10' ] },
	},
	'illumina' => { # Illumina
		'substitution_threshold' => { 'all' => [ '1' ] },
		'indel_threshold' => { 'all' => [ '0.001' ] },
		'min_dominant_frequency_threshold' => { 'all' => [ '25' ] },
		# 'cluster_inframe' => { 'all' => [ '1' ] },
		'min_amplicon_seq_frequency' => { 'all' => [ '3' ] },
		'min_chimera_length' => { 'all' => [ '10' ] },
	},
	'unknown' => { # Unknown technology
		'substitution_threshold' => { 'all' => [ '1' ] },
		'indel_threshold' => { 'all' => [ '1' ] },
		'min_dominant_frequency_threshold' => { 'all' => [ '25' ] },
		# 'cluster_inframe' => { 'all' => [ '1' ] },
		'min_amplicon_seq_frequency' => { 'all' => [ '3' ] },
		'min_chimera_length' => { 'all' => [ '10' ] },
	}
};
$TECH_PARAMS->{'454/iontorrent'} = $TECH_PARAMS->{'454'} = $TECH_PARAMS->{'iontorrent'} = $TECH_PARAMS->{'pyroseq'};
$TECH_PARAMS->{'solexa'} = $TECH_PARAMS->{'illumina'};

my ($INP_amplicons_file, $INP_reads_file, $INP_auto, $INP_tech, $INP_outpath, $INP_denovo, $INP_threads, $INP_shuffle, $INP_nreads, $INP_nocluster, $INP_nofilter, $INP_filter_freq, $INP_cluster_ident, $INP_keepsingletons, $INP_test, $INP_allele_file, $INP_clusref_file, $INP_allele_correction, $INP_verbose, $INP_zip);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s' => \$INP_reads_file,
	'd|data=s' => \$INP_amplicons_file,
	'auto' => \$INP_auto,
	't|tech=s' => \$INP_tech,
	'o|output=s' => \$INP_outpath,
	'dn|denovo' => \$INP_denovo,
	'di|direct' => \$INP_revcomp,
	's|shuffle' => \$INP_shuffle,
	'n|number=i' => \$INP_nreads,
	'na|max=i' => \$INP_nreads_amplicon,
	'fr|freq=f' => \$INP_filter_freq,
	'ci|clid=f' => \$INP_cluster_ident,
	'min=i' => \$INP_min_amplicon_depth,
	'ks|keepsingle' => \$INP_keepsingletons,
	'al|aligned=f' => \$INP_allele_align_params{'aligned'},
	'id|ident=f' => \$INP_allele_align_params{'ident'},
	'nc|nocluster' => \$INP_nocluster,
	'nf|nofilter' => \$INP_nofilter,
	'alcor|allelecorrect' => \$INP_allele_correction,
	'clusrefs=s' => \$INP_clusref_file,
	'a|alleles=s' => \$INP_allele_file,
	'thr|threads=i' => \$INP_threads,
	'test' => \$INP_test,
	'v|verbose' => \$INP_verbose,
	'z|zip' => \$INP_zip,
	'<>' => \&usage,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\n      If analysis parameters are not specified in the CSV file, they will be automatically adjusted by a preliminar analysis of amplicon sequences.\n";
	print "\nOptions:\n";
	print "  -i <file>\tInput FASTQ or FASTA file (compressed or uncompressed) or set of files packed into a unique .ZIP or .TAR.GZ file or AmpliSAS format Excel file.\n";
	print "  -d <file>\tCSV file with primer/amplicon data.\n";
	print "  -o <path>\tOutput folder name.\n";
# 	print "  -dn <tool>\tDe novo read clustering, without marker data.\n";
	print "  -a <file>\tFASTA file with allele names and sequences.\n";
	print "  -auto\t\tParameters are automatically adjusted by a preliminar analysis of amplicons (DEFAULT).\n";
	print "  -t <tech>\tUse recommended technology parameters ('454', 'IonTorrent', 'Illumina', 'Unknown').\n";
	print "  -fr <freq>\tFilter sequences with lower frequency after clustering.\n";
	print "  -ci <number>\tCluster together sequences with higher or equal identity.\n";
	print "  -min <number>\tAmplicons with lower total depth/coverage will be discarded (default=$INP_min_amplicon_depth).\n";
	print "  -ks\t\tKeep singletons after clustering.\n";
	print "  -al <perc>\tMinimum sequence percentage required to align to a reference to be annotated (default=".$INP_allele_align_params{'aligned'}."%).\n";
	print "  -id <perc>\tMinimum percentage of the aligned fragment required to be identical to a reference to be annotated (default=".$INP_allele_align_params{'ident'}."%).\n";
# 	print "  -e\t\tExpand results in 3 Excel sheets (depths, freqs and errors).\n";
	print "  -v\t\tVerbose output, prints additional pseudoFASTA files with details about clustered and filtered sequences.\n";
	print "  -s\t\tShuffle/randomize reads/sequences to analyze.\n";
	print "  -n <number>\tNumber of reads/sequences to analyze.\n";
	print "  -na <number>\tNumber of reads/sequences per amplicon to analyze.\n";
	print "  -nc\t\tDo not cluster sequences.\n";
	print "  -nf\t\tDo not filter sequences.\n";
	print "  -di\t\tAnalyze reads only in direct sense.\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	#   print "  -cluscorr\tCorrect clusters using sequences from alleles file.\n";
	#   print "  -clusrefs\tFASTA file with sequences to use as references in clustering process.\n";
	#   print "  -fast\t\tPerform fast full amplicon multiple alignments instead of more accurate sequence by sequence global alignments during clustering.\n";
	#   print "  -test\t\tRun in 'test' mode and create a dump file of the alignments.\n";
	exit;
}

amplisas();

exit;

#################################################################################

sub amplisas {

	# Downloads reads file if a link is provided
	my $INP_file_url;
	if ($INP_reads_file =~ /(ftp|http|https)\:/ ){
		$INP_file_url = $INP_reads_file;
		printf("\nDownloading file from '%s'\n",$INP_file_url);
		$INP_reads_file = download_file($INP_file_url);
		if (!defined($INP_reads_file)){
			printf("\nERROR: broken link: '%s'.\n\n",$INP_file_url);
			exit;
		}
	}

	# Checks if an excel file or multifle are given as input with sequences previously processed
	my ($INP_excel_file, $INP_multifile);
	if (defined($INP_reads_file) && is_xlsx($INP_reads_file)){
		$INP_excel_file = $INP_reads_file;
		# Next line commented because excel input could be de-multiplexed reads
		# $INP_nocluster = 1;
	# Checks if a set of demultiplexed files is given as input into a compressed file
	} elsif (defined($INP_reads_file) && is_multifile($INP_reads_file)){
		$INP_multifile = $INP_reads_file;
	# Prints usage help if no input file is specified
	} elsif (!defined($INP_reads_file) || !-f $INP_reads_file){
		print "\nERROR: You must specify a sequence input file.\n\n";
		usage();
		exit;
	# De novo clustering accepts single amplicon files:
	} elsif (defined($INP_denovo) && !defined($INP_amplicons_file)){
		$INP_multifile = $INP_reads_file;
	} elsif (!defined($INP_amplicons_file) || !-f $INP_amplicons_file){
		print "\nERROR: You must specify amplicon data input file.\n\n";
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
				'Clustering identity percentage'=>$INP_cluster_ident,
				'Aligned to reference percentage'=>$INP_allele_align_params{'aligned'},
				'Identical to reference percentage'=>$INP_allele_align_params{'ident'},
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
	if (defined($INP_allele_align_params{'aligned'})){
		$INP_allele_align_params{'aligned'} = sprintf("%.2f",$INP_allele_align_params{'aligned'}/100);
	}
	if (defined($INP_allele_align_params{'ident'})){
		$INP_allele_align_params{'ident'} = sprintf("%.2f",$INP_allele_align_params{'ident'}/100);
	}
	# Default output path
	if (!defined($INP_outpath)){
		$INP_outpath  = lc((split('\.',$SCRIPT_NAME))[0]);
	}
	# Tune parameters when AmpliSAS is called by external scripts to do de novo clustering (without primers and tags)
# 	if (defined($INP_denovo)){
# 		# Nothing to do
# 	}

	print "\nRunning '$COMMAND_LINE'\n";

	# Assign dump file, only used in case of -test parameter chosen
	my $dump_file = "$INP_outpath.ampli.dump";

	# 2 PRIMERS => MARKER
	# 1/2 TAGS => SAMPLE
	# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

	# Check and read sequences file
	my ($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads,$ngs_tech);
	if (!defined($INP_excel_file) && !defined($INP_multifile) && (!defined($INP_test) || !-e $dump_file)){
		my $options_ = ['verbose'];
		if (defined($INP_shuffle)){
			push(@$options_,'shuffle');
		}
		($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads)
		= parse_sequence_file($INP_reads_file,$INP_nreads,$options_);
		# Detects the sequencing technology used, if reads file comes from the instrument
		$ngs_tech = ngs_tech($read_headers->[0]);
	}

	# Checks and reads amplicons file
	# Can be missing if an Excel file is given as input
	my ($markerdata,$markers,$sampledata,$tags,$paramsdata,$alleledata);
	if (!defined($INP_excel_file) && !defined($INP_multifile)){
		($markerdata,$markers,$sampledata,$tags,$paramsdata,$alleledata)
		= parse_amplicon_file($INP_amplicons_file,['verbose', 'skip samples']);
	} elsif (defined($INP_amplicons_file)){
		($markerdata,$markers,$sampledata,$tags,$paramsdata,$alleledata)
		= parse_amplicon_file($INP_amplicons_file,['verbose','skip markers','skip samples']);
	}
	# All the files contained into a multifile will be assigned to the same marker, so only one marker is allowed
	if (defined($INP_multifile) && defined($markers) && scalar @{$markers}>1) {
		print "\nERROR: Only ONE marker can be specified because input file contains already demultiplexed sequence files.\n\n";
		exit;
	}
	# If there is no sample data, treat as a single sample named as the file
	if (!defined($tags) || !@$tags){
		my $filename = (reverse(split("\/",$INP_reads_file)))[0];
		# Gives the name of the internal file to the sample (without extension)
		if ($filename =~ /(.+?)\./){
			$tags = [$1];
		# If not gives the full name of the internal file to the sample
		} else {
			$tags = [$filename];
		}
	}

	# Retrieves clustering and filtering parameters from the CSV file
	# If not defined and technology is specified, technogy defaults are used
	# If there are not parameters in CSV file and no tech chosen, automatic mode is activated
	my (@CLUSTERING_PARAMS, @AUTO_CLUSTERING_PARAMS);
	my (@FILTERING_PARAMS, @AUTO_FILTERING_PARAMS);
	unless (defined($INP_nocluster) && defined($INP_nofilter)){
		if (defined($INP_auto) || (!defined($INP_tech) && !defined($paramsdata))){
			$INP_auto = 1;
			print "\n'AUTO MODE' activated.\n".
				"\tThe program will automatically detect the sequencing technology used\n".
				"\tand adjust the sequencing error rates and other analysis parameters.\n";
		}
		# Command line params have priority over CSV file and auto mode
		# Sets minimum per amplicon frequency if it's specified in command line
		if (defined($INP_filter_freq)){
			$paramsdata->{'min_amplicon_seq_frequency'}{'all'} = [ $INP_filter_freq ];
		}
		# Sets minimum identity clustering threshold if it's specified in command line
		if (defined($INP_cluster_ident)){
			$paramsdata->{'identity_threshold'}{'all'} = [ $INP_cluster_ident ];
		}
		# Sets minimum amplicon depth threshold if it's specified in command line
		if (!defined($INP_denovo) && defined($INP_min_amplicon_depth)){
			$paramsdata->{'min_amplicon_depth'}{'all'} = [ $INP_min_amplicon_depth ];
		}
		# Sets keeping singletons parameter if it's specified in command line
		if (defined($INP_keepsingletons)){
			$paramsdata->{'keep_singletons'}{'all'} = [ 'yes' ];
		}
		# Read always the params from CSV file, they will have priority over auto ones
		# Clustering parameters are mandatory
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
		if (!defined($INP_auto) && !@CLUSTERING_PARAMS && !defined($INP_nocluster) && !defined($INP_cluster_ident)) {
			print "\nERROR: No clustering parameters are specified, use auto mode to adjust them automatically. If you want to skip clustering use the option 'nocluster', otherwise use default parameters by specifying the sequencing technology used or include them into CSV data file.\n\n";
			usage();
			exit;
		}
		if (!defined($INP_auto) && !@FILTERING_PARAMS && !defined($INP_nofilter) && !defined($INP_filter_freq)) {
			print "\nERROR: No filtering parameters are specified, use auto mode to adjust them automatically. If you want to skip filtering use the option 'nofilter', otherwise use default parameters by specifying the sequencing technology used or include them into CSV data file.\n\n";
			usage();
			exit;
		}
	}

	# Check and read alleles file (optional)
	# Amplicon data file has preference over alleles in FASTA file
	if (defined($INP_allele_file) && -e $INP_allele_file){
		print "\nReading allele sequences from '$INP_allele_file'.\n";
		$alleledata = read_allele_file($INP_allele_file);
	}

	# Check and read clustering referece sequences file (optional)
	# Amplicon data file has preference over clusrefs in FASTA file
	my $clusrefdata;
	if (defined($INP_clusref_file) && -e $INP_clusref_file){
		print "\nReading clusref sequences from '$INP_clusref_file'.\n";
		my ($clusref_seqs, $clusref_names) = read_fasta_file($INP_clusref_file);
		#my @previous_clusref_seqs = map $clusrefdata->{$_}{'sequence'} , keys %$clusrefdata;
		for (my $i=0; $i<=$#{$clusref_names}; $i++) {
			if (defined($clusrefdata->{$clusref_names->[$i]}) ){
				print "\nERROR: Clustering reference name '".$clusref_names->[$i]."' is duplicated.\n\n";
	# 			exit;
			} else {
				$clusrefdata->{$clusref_names->[$i]}{'sequence'} = $clusref_seqs->[$i];
			}
		}
	}

	# After previous checks to confirm that there are no errors in input data
	# Creates output folders
	my $INP_outpath_allseqs = "$INP_outpath/allseqs";
	my $INP_outpath_clustered = "$INP_outpath/clustered";
	my $INP_outpath_filtered = "$INP_outpath/filtered";
	if (!-d $INP_outpath){
		mkdir($INP_outpath);
	}
	if (!-d $INP_outpath_allseqs){
		mkdir($INP_outpath_allseqs);
	}
	if (!defined($INP_nocluster) && !-d $INP_outpath_clustered){
		mkdir($INP_outpath_clustered);
	}
	if (!defined($INP_nofilter) && !-d $INP_outpath_filtered){
		mkdir($INP_outpath_filtered);
	}

	# Extracts amplicon sequences and depths
	# $amplicon_raw_sequences stores the individual unique sequence depths
	# $amplicon_raw_sequences->{$marker_name}{$sample_name}{$md5} = $depth;
	# $amplicon_raw_depths stores the total depth of the sequences into an amplicon
	# $amplicon_raw_depths->{$marker_name}{$sample_name} += $depth;;
	my ($md5_to_sequence,$amplicon_raw_sequences,$amplicon_raw_depths);
	my ($marker_seq_data, $amplicon_seq_data);
	my $samples;
	my $md5_to_name;
	my ($marker_result_file,$marker_seq_files,$marker_matrix_files,$amplicon_seq_files);
	if (!defined($INP_excel_file) && !defined($INP_multifile)){
		if (!defined($INP_test) || !-e $dump_file){
			if ($INP_align eq 'match') {
				# Free memory
				undef($read_headers);
				print "\nDe-multiplexing amplicon sequences from reads.\n";
				# Creates a file with all sequences
				my $raw_seqs_file = write_to_file("/tmp/".random_file_name(),join("\n",@{$read_seqs}));
				my $match_options;
				if (!$INP_revcomp) {
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
				# Free memory
				undef($read_seqs);
			} #else {
# 				# Aligns reads against primer/tag sequences
# 				print "\nAligning primer/tag sequences.\n";
# 				my $align_amplicon_data
# 				= align_amplicons($read_headers,$read_seqs,$primer_tag_headers,$primer_tag_seqs,$INP_align,$INP_revcomp,$INP_threads,$INP_test,$INP_outpath);
# 				print "\nDe-multiplexing amplicon sequences from reads.\n";
# 				# Loop reads with alignment results and extract amplicons (each amplicon is identified by a MD5 hash)
# 				($md5_to_sequence,$amplicon_raw_sequences,$amplicon_raw_depths)
# 				= match_amplicons($read_headers,$read_seqs,$markerdata,$sampledata,$align_amplicon_data,$INP_nreads_amplicon);
# 				# undef($align_amplicon_data);
# 				# Free memory
# 				undef($read_headers);
# 				undef($read_seqs);
# 
# 			}
		} else {
			print "\nRecovering de-multiplexed amplicon sequences from '$dump_file'.\n";
			($md5_to_sequence,$amplicon_raw_sequences,$amplicon_raw_depths) = @{recover_data_dump($dump_file, 'Storable')};
		}

		# Assigns samples/tags to markers
		foreach my $marker_name (@$markers){
			foreach my $sample_name (@$tags) {
				if (defined($amplicon_raw_depths->{$marker_name}{$sample_name})){
					push(@{$samples->{$marker_name}}, $sample_name);
				} else { # Removes amplicons without data
					delete($amplicon_raw_depths->{$marker_name}{$sample_name});
				}
			}
			# Removes markers without data
			if (!(%{$amplicon_raw_depths->{$marker_name}})){
				delete($amplicon_raw_depths->{$marker_name});
			}
		}
		
		if (!(%$amplicon_raw_depths)){
			print "\nERROR: No amplicons found, check if provided primer and tag sequences are correct.\n\n";
			exit;
		}


	} elsif (defined($INP_multifile)) {

		my $options_ = [];
		if (defined($INP_shuffle)){
			push(@$options_,'shuffle');
		}

		printf("\nReading File '%s'.\n", $INP_multifile);
		my ($markers_, $samples_);
		($markers_, $samples_, $amplicon_raw_sequences, $amplicon_raw_depths, $md5_to_sequence)
		= read_compressed_file_amplicons($INP_multifile,$options_,$INP_nreads_amplicon,$INP_shuffle);

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

	} elsif (defined($INP_excel_file)) {

		printf("\nReading File '%s'.\n", $INP_excel_file);
		my ($markers_, $samples_);
		($markers_, $samples_, $amplicon_raw_sequences, $amplicon_raw_depths, $md5_to_sequence)
		= read_amplisas_file_amplicons($INP_excel_file,[],$INP_nreads_amplicon);

		# Uses the markers specified in the Excel file
		$markers = $markers_;

		# Assigns samples/tags to markers
		foreach my $marker_name (@$markers_){
			foreach my $sample_name (@{$samples_->{$marker_name}}) {
				if (defined($amplicon_raw_depths->{$marker_name}{$sample_name})){
					push(@{$samples->{$marker_name}}, $sample_name);
				} elsif (is_numeric($sample_name) && defined($amplicon_raw_depths->{$marker_name}{sprintf('%d',$sample_name)})){
					# For example sample '001' could be writen into Excel as '1'
					push(@{$samples->{$marker_name}}, sprintf('%d',$sample_name));
				}
			}
		}

	}

	# Align sequences to alleles and assign allele names to sequences
	# Alleles are only assigned if demultiplexing is the last step
	if (defined($alleledata) && %$alleledata && (defined($INP_verbose) || (defined($INP_nocluster) && defined($INP_nofilter)))){
		print "\nMatching allele sequences.\n";
		$md5_to_name = match_alleles($alleledata,$md5_to_sequence,$md5_to_name,\%INP_allele_align_params,$INP_threads);
	}

	# Extracts (from reads) or recalculates (from excel or multifile) marker/amplicon sequence data
	# $marker_seq_data stores the names and parameters of all unique sequences of a unique marker
	# $marker_seq_data->{$marker_name}{$md5} = { 'seq'=> $seq, 'name'=>$name, 'len'=>$len, 'depth'=>$unique_seq_depth, 'samples'=>$count_samples, 'mean_freq'=>$mean_freq, 'max_freq'=>$max_freq, 'min_freq'=>$min_freq};
	# $amplicon_seq_data stores the names and parameters of all unique sequences of a unique amplicon
	# $amplicon_seq_data->{$marker_name}{$sample_name}{$md5} = { 'seq'=> $seq, 'name'=>$name, 'len'=>$len, 'depth'=>$unique_seq_depth, 'freq'=>$unique_seq_frequency, 'cluster_size'=>$cluster_size };
	print "\nExtracting de-multiplexed sequences into '$INP_outpath_allseqs'.\n";
	($marker_seq_data, $amplicon_seq_data, $md5_to_name)
	= retrieve_amplicon_data($markers,$samples,$amplicon_raw_sequences,$amplicon_raw_depths,$md5_to_sequence,$md5_to_name);
	# Prints marker sequences
	($marker_result_file,$marker_seq_files,$marker_matrix_files)
	= print_marker_sequences($markers,$samples,$marker_seq_data,$amplicon_seq_data,$amplicon_raw_depths,$INP_outpath_allseqs);
	# Copy Excel results file to output parent folder
	`cp $marker_result_file $INP_outpath/results.xlsx`;
	# Prints amplicon sequences
	$amplicon_seq_files
	= print_amplicon_sequences($markers,$samples,$marker_seq_data,$amplicon_seq_data,$INP_outpath_allseqs);
	# Prints all amplicon sequences into an unique file
	# system("cat ".join(" ",@{$amplicon_seq_files}."> $INP_outpath/allseqs.fasta");

	
	# Detects automatically the sequencing technoloy
	# And adjusts analysis parameters
	my ($tech_auto, $lengths_auto, $params_auto);
	if (defined($INP_auto) && !defined($INP_nocluster)){
		if (!defined($INP_cluster_ident)){
			print "\nChecking data and setting automatically clustering parameters.\n";
		} else {
			print "\nChecking automatically sequencing technology parameters.\n";
		}
		($tech_auto, $lengths_auto, $params_auto) = adjust_automatic_params($markers,$samples,$marker_seq_data,$amplicon_seq_data);
		if (defined($ngs_tech)){
			printf("\tTechnology: %s\n", ucfirst($ngs_tech));
		} elsif (defined($tech_auto)){
			printf("\tTechnology: %s\n", ucfirst($tech_auto));
		}
		printf("\tError rates: %s%% substitutions, %s%% indels\n",$params_auto->{'substitution_threshold'}{'all'}[0] ,$params_auto->{'indel_threshold'}{'all'}[0]);
		# Activates fast clustering mode if there are not indels
		if ($params_auto->{'indel_threshold'}{'all'}[0] <= 0.001){
			$params_auto->{'fast_clustering_alignment'}{'all'}[0] = 1;
		}
		# Sets clustering parameters with the values retrieved by 'adjust_automatic_params'
		# 'indel_threshold', 'substitution_threshold', 'cluster_inframe', 'min_dominant_frequency_threshold'
		# Unless the parameter has been adjusted before on command line or CSV file
		foreach my $paramname (keys %$params_auto){
			# De novo clustering mixes in the same data several amplicons so length-based thresholds are useless
			if (defined($INP_denovo) && $paramname !~ /indel_threshold|substitution_threshold|min_dominant_frequency_threshold|fast_clustering_alignment/){
				next;
			}
			# Cluster by identity doesn't require any other clustering thresholds
			if (defined($INP_cluster_ident) && $paramname !~ /fast_clustering_alignment/){
				next;
			}
			my $auto_param = 0;
			foreach my $marker_name (keys %{$params_auto->{$paramname}}){
				if (!defined($paramsdata->{$paramname}{'all'}) && !defined($paramsdata->{$paramname}{$marker_name})){
					$paramsdata->{$paramname}{$marker_name} = $params_auto->{$paramname}{$marker_name};
					$auto_param = 1;
				}
			}
			if ($auto_param){
				push(@AUTO_CLUSTERING_PARAMS,$paramname);
			}
		}
		# De novo clustering mixes in the same data several amplicons so automatic 'min_amplicon_depth' is not recommended
		if (!defined($paramsdata->{'min_amplicon_depth'}) && !defined($INP_denovo) && !defined($INP_nofilter) && defined($TECH_PARAMS->{lc($tech_auto)}) && defined($TECH_PARAMS->{lc($tech_auto)}{'min_amplicon_depth'}{'all'}[0])){
			$paramsdata->{'min_amplicon_depth'}{'all'}[0] = $TECH_PARAMS->{lc($tech_auto)}{'min_amplicon_depth'}{'all'}[0];
			push(@AUTO_FILTERING_PARAMS,'min_amplicon_depth');
			printf("\tMinimum amplicon depth: %s\n",$paramsdata->{'min_amplicon_depth'}{'all'}[0]);
		}
	}
	
	# De novo clustering mixes in the same data several amplicons
	if (!defined($INP_denovo)){
		print "\nChecking data and setting marker lengths.\n";
		# Checks if marker lengths are defined, if not they will be calculated automatically if they are required
		foreach my $marker_name (@$markers){
			# Skips markers without data
			if (!defined($amplicon_seq_data->{$marker_name})){
				delete($amplicon_seq_data->{$marker_name});
				next;
			}
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
		# Remove markers without data
		if (!(%{$amplicon_seq_data->{$marker_name}})){
			delete($amplicon_seq_data->{$marker_name});
		}
	}
	# If all amplicons have been excluded could be because a high amplicon depth threshold
	if (!(%$amplicon_seq_data)){
		print "\nERROR: All amplicons were discarded, try to decrease the 'minimum amplicon depth' threshold value.\n\n";
		exit;
	}
	my $amplicon_seq_raw_data = { %$amplicon_seq_data };

	# Cluster sequences if not specified the opposite
	# $clusters stores the clusters of a unique amplicon
	# $clusters->{$marker_name}{$sample_name} = array( array({'md5'=> $md5, 'seq'=>$allele_seq, 'depth'=>$allele_depth, 'errors'=>$errors }, { ... }, ...), array({ ... }, ... ) );
	my ($amplicon_clustered_sequences,$amplicon_clustered_depths,$clusters,$amplicon_clusters);
	if (!defined($INP_nocluster)){
		print "\nClustering amplicon sequences with the following parameters ('threshold' 'marker' 'values'):\n";
		# In $paramsdata hash there are also sequence filters
		foreach my $clustering_param (@AUTO_CLUSTERING_PARAMS){
			if (defined($paramsdata->{$clustering_param})) {
				foreach my $marker_name (sort keys %{$paramsdata->{$clustering_param}}){
					printf("\t*%s\t%s\t%s\n",$clustering_param,$marker_name,join(',',@{$paramsdata->{$clustering_param}{$marker_name}}));
				}
			}
		}
		foreach my $clustering_param (@CLUSTERING_PARAMS){
			if (defined($paramsdata->{$clustering_param})) {
				foreach my $marker_name (sort keys %{$paramsdata->{$clustering_param}}){
					printf("\t%s\t%s\t%s\n",$clustering_param,$marker_name,join(',',@{$paramsdata->{$clustering_param}{$marker_name}}));
				}
			}
		}
		if (@AUTO_CLUSTERING_PARAMS) {
			print "\n\t*WARNING: These parameters have been automatically adjusted and may require manual verification.\n";
		}
		print "\n";
		if (!defined($INP_denovo) && defined($INP_threads) && $INP_threads>1){
			($amplicon_clustered_sequences,$amplicon_clustered_depths,$clusters,$md5_to_sequence)
			= cluster_amplicon_sequences_with_threads($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$clusrefdata,$INP_threads);
		} elsif(!defined($INP_denovo)) {
			($amplicon_clustered_sequences,$amplicon_clustered_depths,$clusters,$md5_to_sequence)
			= cluster_amplicon_sequences($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$clusrefdata);
		} elsif (defined($INP_denovo) && defined($INP_threads) && $INP_threads>1){
			($amplicon_clustered_sequences,$amplicon_clustered_depths,$clusters,$md5_to_sequence)
			= cluster_denovo_amplicon_sequences_with_threads($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$INP_threads);
		} elsif(defined($INP_denovo)) {
			($amplicon_clustered_sequences,$amplicon_clustered_depths,$clusters,$md5_to_sequence)
			= cluster_denovo_amplicon_sequences($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data);
		}

		# If cluster correction with reference alleles is desired, cluster members exactly matching sequences in alleles file will be promoted as dominant in the cluster
		if (defined($INP_allele_correction) && defined($alleledata)){
			print "\nCorrecting clusters with allele sequences given in file '$INP_allele_file'.\n";
			my %allele_seqs;
			map $allele_seqs{$alleledata->{$_}{'sequence'}} = 1, keys %$alleledata;
			foreach my $marker_name (@$markers){
				if (!defined($clusters->{$marker_name})){
					next;
				}
				foreach my $sample_name (@{$samples->{$marker_name}}) {
					if (!defined($clusters->{$marker_name}{$sample_name})){
						next;
					}
					foreach my $cluster (@{$clusters->{$marker_name}{$sample_name}}) {
						my $first_member = $cluster->[0];
						for (my $i=0; $i<=$#{$cluster}; $i++){
							my $md5 = $cluster->[$i]{'md5'};
							if (defined($allele_seqs{$cluster->[$i]{'seq'}})){
								if ($i != 0){
									$amplicon_clustered_sequences->{$marker_name}{$sample_name}{$cluster->[$i]{'md5'}} = $amplicon_clustered_sequences->{$marker_name}{$sample_name}{$first_member->{'md5'}};
									delete($amplicon_clustered_sequences->{$marker_name}{$sample_name}{$first_member->{'md5'}});
									$md5_to_sequence->{$cluster->[$i]{'md5'}} = $cluster->[$i]{'seq'};
									$cluster->[0] = $cluster->[$i];
									$cluster->[$i] = $first_member;
									$cluster->[0]{'errors'} = 'CORRECTED';
									my $ref_var_aligned_seqs = align_seqs2one($cluster->[0]{'seq'},$cluster->[0]{'md5'},[map $_->{'seq'}, @$cluster],[map $_->{'md5'}, @$cluster],'needleall');
									for (my $j=1; $j<=$#{$cluster}; $j++){
										my ($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($ref_var_aligned_seqs->{$cluster->[$j]{'md5'}}[0],$ref_var_aligned_seqs->{$cluster->[$j]{'md5'}}[1]);
										$cluster->[$j]{'errors'} = sprintf("sub: %d, ins: %d, del: %d, homo_indel: %d", scalar @$substitutions, scalar @$insertions, scalar @$deletions, scalar @$homopolymer_indels);
										print '';
									}
								}
								last;
							}
						}
					}
				}
			}
		}

		# Writes FASTA files with real allele sequences and artifacts clustered together
		# And annotates in $amplicon_clusters the number of members per cluster (for filtering purposes)
		# $amplicon_clusters->{$marker_name}{$sample_name}{$md5} = $count_members;
		# if (defined($INP_verbose)){
			print "\nPrinting information about clustered and not clustered sequences into '$INP_outpath_clustered'.\n";
			foreach my $marker_name (@$markers){
				if (!defined($clusters->{$marker_name})){
					next;
				}
				foreach my $sample_name (@{$samples->{$marker_name}}) {
					if (!defined($clusters->{$marker_name}{$sample_name})){
						next;
					}
					# Prints debug sequences
					my @clusters_output;
					my $count_clusters = 0;
					foreach my $cluster (@{$clusters->{$marker_name}{$sample_name}}) {
						$count_clusters++;
						my $count_members = 0;
						foreach my $cluster_member (@$cluster) {
							$count_members++;
							my $md5 = $cluster_member->{'md5'};
							my $errors = $cluster_member->{'errors'};
							# Depth will be between parenthesis if the variant is high frequency
							# and has the same similarity with several dominant sequences
							my $depth = $cluster_member->{'depth'};
							if (defined($amplicon_seq_data->{$marker_name}{$sample_name}{$md5})){
								my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
								my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
								my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
								#my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
								my $frequency = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'};
								my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
								my $mean_freq = $marker_seq_data->{$marker_name}{$md5}{'mean_freq'};
								my $min_freq = $marker_seq_data->{$marker_name}{$md5}{'min_freq'};
								my $max_freq = $marker_seq_data->{$marker_name}{$md5}{'max_freq'};
								if ($#{$cluster}>0){
									no warnings; # No warning if any valued is undefined
									if ($count_members==1){
										push(@clusters_output, sprintf(">*%s | Cluster %d: dominant | hash=%s | len=%d | depth=%s | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f\n%s", $name, $count_clusters, $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq, $seq));
									} else {
										push(@clusters_output, sprintf(">#%s | Cluster %d: member, %s | hash=%s | len=%d | depth=%s | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f\n%s", $name, $count_clusters, $errors, $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq, $seq));
									}
								} else {
									no warnings; # No warning if any valued is undefined
									if (!$errors){
										push(@clusters_output, sprintf(">*%s | hash=%s | len=%d | depth=%d | freq=%.2f | samples=%s | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f\n%s", $name, $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq, $seq));
									} else {
										push(@clusters_output, sprintf(">#%s | %s | hash=%s | len=%d | depth=%d | freq=%.2f | samples=%s | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f\n%s", $name, $errors, $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq, $seq));
									}
								}
							} else { # For consensus sequences:
								my $seq = $cluster_member->{'seq'};
								my $len = length($seq);
								if ($count_members==1){
									push(@clusters_output, sprintf(">*%s | Cluster %d: dominant | hash=%s | len=%d\n%s", 'CONSENSUS', $count_clusters, $md5, $len, $seq));
								} else {
									push(@clusters_output, sprintf(">#%s | Cluster %d: member, %s | hash=%s | len=%d\n%s", 'CONSENSUS', $count_clusters, $errors, $md5, $len, $seq));
								}
								# Includes consensus sequences into $md5_to_sequence
# 								$md5_to_sequence->{$md5} = $seq;
							}
						}
						# Line break between clusters:
						# push(@clusters_output, '');
						# Do not save non clustered sequences
						if ($cluster->[0]{'errors'} ne 'no cluster'){
							$amplicon_clusters->{$marker_name}{$sample_name}{$cluster->[0]{'md5'}} = $count_members;
						}
					}
					# If no sequences are found and tries to write to file will exit the program with errors
					if (@clusters_output){
						write_to_file("$INP_outpath_clustered/$marker_name-$sample_name.verbose.fasta",join("\n",@clusters_output));
					}
				}
			}
# 		# Only annotates cluster sizes
# 		} else {
# 			foreach my $marker_name (@$markers){
# 				if (!defined($clusters->{$marker_name})){
# 					next;
# 				}
# 				foreach my $sample_name (@{$samples->{$marker_name}}) {
# 					if (!defined($clusters->{$marker_name}{$sample_name})){
# 						next;
# 					}
# 					foreach my $cluster (@{$clusters->{$marker_name}{$sample_name}}) {
# 						my $count_members = 0;
# 						foreach my $cluster_member (@$cluster) {
# 							$count_members++;
# 						}
# 						# Do not save non clustered sequences
# 						if ($cluster->[0]{'errors'} ne 'no cluster'){
# 							$amplicon_clusters->{$marker_name}{$sample_name}{$cluster->[0]{'md5'}} = $count_members;
# 						}
# 						# Includes consensus sequences into $md5_to_sequence and into RAW data
# # 						foreach my $cluster_member (@$cluster) {
# # 							my $md5 = $cluster_member->{'md5'};
# # 							my $seq = $cluster_member->{'seq'};
# # 							if (!defined($md5_to_sequence->{$md5})){
# # 								$md5_to_sequence->{$md5} = $seq;
# # 							}
# # 						}
# 					}
# 				}
# 			}
# 		}
		undef($clusters);
		
		# If no clusters are retrieved, it will be probably because wrong sequences lengths specified in the input
		if (!defined($amplicon_clusters)){
			print "\nERROR: No clusters retrieved, check if provided sequence lengths are correct or remove clustering length restrictions.\n\n";
			exit;
		}

		# Aligns sequences to alleles and assign allele names to sequences
		# Alleles are only assigned if clustering is the last step
		if (defined($alleledata) && %$alleledata && (defined($INP_verbose) || defined($INP_nofilter))){
			print "\nMatching allele sequences.\n";
			$md5_to_name = match_alleles($alleledata,$md5_to_sequence,$md5_to_name,\%INP_allele_align_params,$INP_threads);
		}

		# Extracts marker/amplicon sequence data
		print "\nExtracting clustered sequences into '$INP_outpath_clustered'.\n";
		# Frequencies will be calculated with original depths, not after clustering ones
		($marker_seq_data, $amplicon_seq_data, $md5_to_name)
		= retrieve_amplicon_data($markers,$samples,$amplicon_clustered_sequences,$amplicon_raw_depths,$md5_to_sequence,$md5_to_name,$amplicon_clusters);
		# Prints marker sequences
		($marker_result_file,$marker_seq_files,$marker_matrix_files)
		= print_marker_sequences($markers,$samples,$marker_seq_data,$amplicon_seq_data,$amplicon_raw_depths,$INP_outpath_clustered,$amplicon_raw_sequences);
		# Copy Excel results file to output parent folder
		`cp $marker_result_file $INP_outpath/results.xlsx`;
		# Prints amplicon sequences
		my $amplicon_seq_files
		= print_amplicon_sequences($markers,$samples,$marker_seq_data,$amplicon_seq_data,$INP_outpath_clustered);
	} else {
		($amplicon_clustered_sequences,$amplicon_clustered_depths)
		= ($amplicon_raw_sequences,$amplicon_raw_depths);

	}

	# Checks the clusters to decide which frequency threshold to use in filtering
	if (defined($INP_auto) && !defined($INP_nofilter)){

		if (!defined($INP_filter_freq) && !defined($paramsdata->{'min_amplicon_seq_frequency'}{'all'})){
			my $auto_param = 0;
			# Finds the minimum frequency threshold taking the median from an array with the frequencies of the last detected DOC alleles
			foreach my $marker_name (@$markers){
				if (defined($paramsdata->{'min_amplicon_seq_frequency'}{$marker_name}) || !defined($amplicon_clustered_sequences->{$marker_name})){
					next;
				}
				my $marker_allele_numbers;
				my @marker_min_amplicon_freqs_;
				foreach my $sample_name (@{$samples->{$marker_name}}) {
					if (!defined($amplicon_clustered_sequences->{$marker_name}{$sample_name})){
						next;
					}
					my $amplicon_total_depth = sum( map { $amplicon_raw_sequences->{$marker_name}{$sample_name}{$_} } keys %{$amplicon_raw_sequences->{$marker_name}{$sample_name}} );
					my @amplicon_sorted_depths = sort {$b<=>$a} values %{$amplicon_clustered_sequences->{$marker_name}{$sample_name}};
					my ($doc_allele_number,$DOCns) = degree_of_change(\@amplicon_sorted_depths, scalar @amplicon_sorted_depths);
					if (defined($amplicon_sorted_depths[$doc_allele_number-1])){
						my $min_amplicon_freq = sprintf("%.2f", 100 * $amplicon_sorted_depths[$doc_allele_number-1] / $amplicon_total_depth);
						push(@{$marker_allele_numbers->{$doc_allele_number}},$min_amplicon_freq);
					}
				}
				# Removes frequencies from homocygous loci if there are heterocygous
				if (scalar keys %$marker_allele_numbers > 1 && defined($marker_allele_numbers->{1})) {
					delete($marker_allele_numbers->{1});
				}
				# Sorts frequencies
				my @marker_min_amplicon_freqs;
				map push(@marker_min_amplicon_freqs,@{$marker_allele_numbers->{$_}}), keys %$marker_allele_numbers;
				$paramsdata->{'min_amplicon_seq_frequency'}{$marker_name}[0] = sprintf("%.2f", 0.9 * median(@marker_min_amplicon_freqs));
				$auto_param = 1;
			}
			# Includes the new filtering parameters in the list
			if ($auto_param) {
				push(@AUTO_FILTERING_PARAMS, 'min_amplicon_seq_frequency');
			}
		}

# 		elsif (defined($INP_filter_freq)) { # If a minimum per amplicon frequency is defined in the command line
# 
# 			$paramsdata->{'min_amplicon_seq_frequency'}{'all'} = [ $INP_filter_freq ];
# 
# 		}

		# Sets by default chimera filtering
		if (!defined($paramsdata->{'min_chimera_length'})){
			$paramsdata->{'min_chimera_length'}{'all'}[0] = 10;
			# Includes the new filtering parameters in the list
			push(@AUTO_FILTERING_PARAMS, 'min_chimera_length');
		}


	}

	# FILTERING:
	# if (defined($INP_excel_file)){
	# 	printf("\nReading File '%s'.\n", $INP_reads_file);
	# 	($markers, $samples, $amplicon_clustered_sequences, $amplicon_clustered_depths, $md5_to_sequence) = read_amplisas_file_amplicons($INP_excel_file,[],$INP_nreads_amplicon);
	# 	($amplicon_raw_sequences,$amplicon_raw_depths) = ($amplicon_clustered_sequences,$amplicon_clustered_depths);
	# 	$INP_outpath_filtered = "$INP_outpath/filtered";
	# 	if (!-d $INP_outpath){
	# 		mkdir($INP_outpath);
	# 	}
	# 	if (!defined($INP_nofilter) && !-d $INP_outpath_filtered){
	# 		mkdir($INP_outpath_filtered);
	# 	}
	# }

	# Apply filters to amplicon sequences (clustered or not)
	my ($amplicon_filtered_sequences, $amplicon_filtered_depths,$filters_output);
	if (!defined($INP_nofilter)){
		print "\nFiltering sequences with the following criteria ('filter' 'marker' 'values'):\n";
		# In $paramsdata hash there are also clustering thresholds
		foreach my $filtering_param (@AUTO_FILTERING_PARAMS){
			if (defined($paramsdata->{$filtering_param})) {
				foreach my $marker_name (sort keys %{$paramsdata->{$filtering_param}}){
					printf("\t*%s\t%s\t%s\n",$filtering_param,$marker_name,join(',',@{$paramsdata->{$filtering_param}{$marker_name}}));
				}
			}
		}
		foreach my $filtering_param (@FILTERING_PARAMS){
			if ($filtering_param ne 'allowed_markers' && defined($paramsdata->{$filtering_param})) {
				foreach my $marker_name (sort keys %{$paramsdata->{$filtering_param}}){
					printf("\t%s\t%s\t%s\n",$filtering_param,$marker_name,join(',',@{$paramsdata->{$filtering_param}{$marker_name}}));
				}
			} elsif (defined($paramsdata->{$filtering_param})) {
				printf("\t%s\t%s\n",$filtering_param,join(',',@{$paramsdata->{$filtering_param}}));
			}
		}
		if (@AUTO_FILTERING_PARAMS) {
			print "\n\t*WARNING: These parameters have been automatically adjusted and may require manual verification.\n";
		}
		print "\n";
		if (defined($INP_threads) && $INP_threads>1){
			($amplicon_filtered_sequences,$amplicon_filtered_depths,$filters_output,$md5_to_sequence)
			= filter_amplicon_sequences_with_threads($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$amplicon_raw_depths,$INP_threads);
		} else {
			($amplicon_filtered_sequences,$amplicon_filtered_depths,$filters_output,$md5_to_sequence)
			= filter_amplicon_sequences($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$amplicon_raw_depths);
		}
		# Clustering cross-checking to recover lost real alleles with low frequency or in small clusters
		# It cannot be included in 'filter_amplicon_sequences_with_threads' function because it needs the final results
		foreach my $marker_name (keys %{$amplicon_filtered_sequences}){
			my $min_samples_to_keep;
			if (defined($paramsdata->{'min_samples_to_keep'}{$marker_name})){
				$min_samples_to_keep = $paramsdata->{'min_samples_to_keep'}{$marker_name}[0];
			} elsif (defined($paramsdata->{'min_samples_to_keep'}{'all'})){
				$min_samples_to_keep = $paramsdata->{'min_samples_to_keep'}{'all'}[0];
			} else {
				next;
			}
			my $min_amplicon_seq_frequency_to_keep;
			if (defined($paramsdata->{'min_amplicon_seq_frequency_to_keep'}{$marker_name})){
				$min_amplicon_seq_frequency_to_keep = $paramsdata->{'min_amplicon_seq_frequency_to_keep'}{$marker_name}[0];
			} elsif (defined($paramsdata->{'min_amplicon_seq_frequency_to_keep'}{'all'})){
				$min_amplicon_seq_frequency_to_keep = $paramsdata->{'min_amplicon_seq_frequency_to_keep'}{'all'}[0];
			}
			# Annotates in how many sequences is present a sequence after filtering
			my $filtered_marker_md5s;
			foreach my $sample_name (keys %{$amplicon_filtered_sequences->{$marker_name}}){
				foreach my $md5 (keys %{$amplicon_filtered_sequences->{$marker_name}{$sample_name}}){
					push(@{$filtered_marker_md5s->{$md5}},$sample_name);
				}
			}
			foreach my $sample_name (keys %{$amplicon_filtered_sequences->{$marker_name}}){
				# Checks all the sequences before clustering and filtering for lost real alleles
	# 			foreach my $md5 (keys %{$amplicon_seq_raw_data->{$marker_name}{$sample_name}}){
				foreach my $md5 (keys %{$amplicon_seq_data->{$marker_name}{$sample_name}}){
	# if ($md5 eq '73ce3ddd162a2f04b061dfe658357ce6'){
	# print '';
	# }
					# Annotate sequence if it is present as clear allele after filtering in other samples
					if (!defined($amplicon_filtered_sequences->{$marker_name}{$sample_name}{$md5}) && defined($filtered_marker_md5s->{$md5}) && scalar @{$filtered_marker_md5s->{$md5}} >= $min_samples_to_keep ){
						my $seq_data;
						# If the sequence exists is a cluster ref, then use data from clustering
						if (defined($amplicon_seq_data->{$marker_name}{$sample_name}{$md5})) {
							if (!defined($min_amplicon_seq_frequency_to_keep) || $min_amplicon_seq_frequency_to_keep<=$amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'}) {
								$seq_data = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5};
							}
						# If the sequence has been clustered as a variant of another sequence, use the data before clustering
						} else {
							if (!defined($min_amplicon_seq_frequency_to_keep) || $min_amplicon_seq_frequency_to_keep<=$amplicon_seq_raw_data->{$marker_name}{$sample_name}{$md5}{'freq'}) {
								$seq_data = $amplicon_seq_raw_data->{$marker_name}{$sample_name}{$md5};
							}
						}
						if (!defined($seq_data)){
							next;
						}
						my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
						my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
						my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
						my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
						my $frequency = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'};
						my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
						my $mean_freq = $marker_seq_data->{$marker_name}{$md5}{'mean_freq'};
						my $max_freq = $marker_seq_data->{$marker_name}{$md5}{'max_freq'};
						my $min_freq = $marker_seq_data->{$marker_name}{$md5}{'min_freq'};
						my $header;
						if (defined($amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'cluster_size'})){
							my $cluster_size = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'cluster_size'};
							$header = sprintf("hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | cluster_size=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $md5, $len, $depth, $frequency, $count_samples, $cluster_size, $mean_freq, $max_freq, $min_freq);
						} else {
							$header = sprintf("hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq);
						}
						$amplicon_filtered_sequences->{$marker_name}{$sample_name}{$md5} = $seq_data->{'depth'};
						$amplicon_filtered_depths->{$marker_name}{$sample_name} += $seq_data->{'depth'};
						$filters_output->{$marker_name}{$sample_name} .= sprintf(">%s | Present in several samples: %s |%s\n%s\n", $name, join(',',@{$filtered_marker_md5s->{$md5}}), $header, $seq);
					}
				}
			}
		}
		# Writes FASTA files with real allele sequences and artifacts clustered together
# 		if (defined($INP_verbose)){
			print "\nPrinting information about filtered and non filtered sequences into '$INP_outpath_filtered'.\n";
			foreach my $marker_name (@$markers){
				if (!defined($filters_output->{$marker_name})){
					next;
				}
				foreach my $sample_name (@{$samples->{$marker_name}}) {
					if (!defined($filters_output->{$marker_name}{$sample_name})){
						next;
					}
					write_to_file("$INP_outpath_filtered/$marker_name-$sample_name.verbose.fasta",$filters_output->{$marker_name}{$sample_name});
				}
			}
# 		}
		if (!(%$amplicon_filtered_sequences)){
			`rm -rf $INP_outpath`;
			print "\nERROR: No sequences passed the filters, please choose broader filtering values.\n\n";
			exit;
		} else {

			# Aligns sequences to alleles and assign allele names to sequences
			# Alleles are only assigned if clustering is the last step
			if (defined($alleledata) && %$alleledata){
				print "\nMatching allele sequences.\n";
				$md5_to_name = match_alleles($alleledata,$md5_to_sequence,$md5_to_name,\%INP_allele_align_params,$INP_threads);
			}

			print "\nExtracting filtered sequences into '$INP_outpath_filtered'.\n";
			# Frequencies will be calculated with original depths, not after clustering+filtering ones
			($marker_seq_data, $amplicon_seq_data, $md5_to_name)
			= retrieve_amplicon_data($markers,$samples,$amplicon_filtered_sequences,$amplicon_raw_depths,$md5_to_sequence,$md5_to_name);
			($marker_result_file,$marker_seq_files,$marker_matrix_files)
			= print_marker_sequences($markers,$samples,$marker_seq_data,$amplicon_seq_data,$amplicon_raw_depths,$INP_outpath_filtered,$amplicon_raw_sequences);
			# Copy Excel results file to output parent folder
			`cp $marker_result_file $INP_outpath/results.xlsx`;
			$amplicon_seq_files
			= print_amplicon_sequences($markers,$samples,$marker_seq_data,$amplicon_seq_data,$INP_outpath_filtered);
		}
	} else {
		($amplicon_filtered_sequences, $amplicon_filtered_depths)
		= ($amplicon_clustered_sequences,$amplicon_clustered_depths)

	}

	# Print number of unique sequences per amplicon
	my $summary_output = "Amplicon\tTotal\tUnique";
	if (!defined($INP_nocluster)){
		$summary_output .= "\tReads-clustered\tVariants-clustered";
	}
	if (!defined($INP_nofilter)){
		$summary_output .= "\tReads-filtered\tVariants-filtered";
	}
	$summary_output .= "\n";
	foreach my $marker_name (@$markers){
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			$summary_output .= sprintf("%s-%s",$marker_name,$sample_name);
			if (defined($amplicon_raw_depths->{$marker_name}{$sample_name})){
				$summary_output .= sprintf("\t%d\t%d",$amplicon_raw_depths->{$marker_name}{$sample_name},scalar keys %{$amplicon_raw_sequences->{$marker_name}{$sample_name}});
			} else {
				$summary_output .= sprintf("\t%d\t%d",0,0);
			}
			if (!defined($INP_nocluster)){
				if (defined($amplicon_clustered_depths->{$marker_name}{$sample_name})){
					$summary_output .= sprintf("\t%d\t%d",$amplicon_clustered_depths->{$marker_name}{$sample_name},scalar keys %{$amplicon_clustered_sequences->{$marker_name}{$sample_name}});
				} else {
					$summary_output .= sprintf("\t%d\t%d",0,0);
				}
			}
			if (!defined($INP_nofilter)){
				if (defined($amplicon_filtered_depths->{$marker_name}{$sample_name})){
					$summary_output .= sprintf("\t%d\t%d",$amplicon_filtered_depths->{$marker_name}{$sample_name},scalar keys %{$amplicon_filtered_sequences->{$marker_name}{$sample_name}});
				} else {
					$summary_output .= sprintf("\t%d\t%d",0,0);
				}
			}
			$summary_output .= "\n";
		}
	}
	print "\nReads per amplicon:\n$summary_output\n";
	write_to_file("$INP_outpath/summary.txt",$summary_output);


	# Print amplicon parameters into a file
	print "\nPrinting amplicon data into '$INP_outpath/amplicon_data.csv'.\n";
	# my $amplicon_data = ">command\n$COMMAND_LINE\n\n";
	my $amplicon_data = '';
	$amplicon_data .= print_amplicon_data($markerdata,'markers')."\n";
	$amplicon_data .= print_amplicon_data($sampledata,'samples')."\n";
	if (defined($paramsdata) && %$paramsdata){
		$amplicon_data .= print_amplicon_data($paramsdata,'params',[unique(@CLUSTERING_PARAMS, @AUTO_CLUSTERING_PARAMS)],[unique(@FILTERING_PARAMS, @AUTO_FILTERING_PARAMS)])."\n";
	} # elsif (!defined($INP_nocluster) && defined($INP_nofilter)){
# 		$amplicon_data .= print_amplicon_data($paramsdata,'params',\@CLUSTERING_PARAMS)."\n";
# 	} elsif (defined($INP_nocluster) && !defined($INP_nofilter)){
# 		$amplicon_data .= print_amplicon_data($paramsdata,'params',\@FILTERING_PARAMS)."\n";
# 	}
	if (defined($alleledata)){
		$amplicon_data .= print_amplicon_data($alleledata,'alleles')."\n";
	}
	write_to_file("$INP_outpath/amplicon_data.csv",$amplicon_data);


	if (defined($INP_zip) && -d $INP_outpath && !is_folder_empty($INP_outpath)){
		my $cwd = getcwd;
		chdir($INP_outpath);
		my $outname = basename($INP_outpath);
		`zip -qrm $outname.zip *` ;
		`mv $outname.zip ..`;
		chdir($cwd);
		rmdir($INP_outpath);
		print "\nAnalysis results stored into '$INP_outpath.zip'.\n\n";
		return "$INP_outpath.zip";
	} elsif (-d $INP_outpath && !is_folder_empty($INP_outpath)){
		print "\nAnalysis results stored into '$INP_outpath'.\n\n";
		return $INP_outpath;
	} else {
		print "\nThere was some error in the analysis and no results were retrieved.\n\n";
		return 0;
	}

}

#################################################################################


