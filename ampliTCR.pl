#!/usr/bin/perl -w

# Looks for TCR primer sequences in read sequence file (FASTA or FASTQ compressed or uncompressed)
# Extracts TCR sequences trimmed on primer location
# Prints statistics of TCR sequence lengths
#
# Examples:
#  Clusters and extracts TCR variable segments based in a protein pattern and assigns names by similarity to references
#    perl ampliTCR.pl -i TCR-126_S1_merged.fq.gz -o TCR_126
#  Clusters and extracts TCR variable segments based in a protein pattern and assigns names by similarity to references
#    perl ampliTCR.pl -d -p -debug -thr 30 -vpat 'C.{60,}?C' -i TCR-126_S1_merged.fq.gz -o TCR_126 -vref unique_TCRBV.fa
#    perl ampliTCR.pl -d -p -regex -vpat 'C\w{10,12}Y\w[QK]\w{37,40}[LMI]\w{12}[YL]\wC' -debug -thr 30 -i TCR-126_S1_merged.fq.gz -o TCR_126 -vref unique_TCRBV.fa

my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Analyzes a set of genomic or transcriptomic TCR sequences recognizing and extracting their Variable, Joining, Diversity, CDR3 and/or Constant segments.";

# Modules are in folder '../' in the path of the script
use lib "lib";
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Bio::Sequences;
use Bio::Ampli;
use IO::Compress::Gzip qw(gzip $GzipError) ;
use Sort::Naturally;
# use Data::Dumper;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options:
# Default TCR seqment sequences clustering algorithm
my $DEFAULT_CLUSTERING = 'amplisas';
# TCR segment allele matching parameters
my $INP_allele_align = { 'alignment' => 'dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 100', 'aligned' => 1, 'ident' => 1 };

# Regular pattern to recognize TRV & TRJ segments
# TCR IMGT Numering (residue numbers are not correlative, because the are designed to number all IG and TCR chains)
# Variable: Glutamine (Q6), Proline (P8), Cys (C23), Trp (W41), Tyr(Y42), Leu/Met (L/M89), Cys (C104)
# CDR1: 27-38, CDR2: 56-65, CDR3: 105-117
# $INP_motif = '(.+C).+?(G\wG\w{2}L\w[VI])';
# Default RACE TCR primer in 5'->3' sense:
# my @DEFAULT_TCR_PRIMERS = ('GGACTCACCTTGCTCAGATCCT');

# Default patterns to recognize TCR beta segments in bank vole (valid for mouse and human):
my %DEFAULT_TCR_PATTERNS = (
	'tcrv' => 'Q\w[PS]\w{14}C\w{10,11}WY\w{39,42}[LM]\w{14}C',
# 	'tcrv' => 'C.{60,}?C';
	'tcrj' => 'G\wG\w{2}L\w[VI]',
	'tcrc' => 'EDL.+', # \wDLSKVS',
	'tcrd' => 'gggaca\w{6}|gggact\w{8}',
);

# Define thresholds for real sequences
my $CLUSTERING_THRESHOLDS;
# Clustering identity thresholds per segment:
$CLUSTERING_THRESHOLDS->{'tcrv'}{'identity'} = 95; # 95%
$CLUSTERING_THRESHOLDS->{'tcrj'}{'identity'} = 85; # 85%
$CLUSTERING_THRESHOLDS->{'tcrd'}{'identity'} = 80; # 80%
$CLUSTERING_THRESHOLDS->{'tcrc'}{'identity'} = 95; # 95%
# Minimum depth ratio of an unique sequence respect to the total number of seqs:
$CLUSTERING_THRESHOLDS->{'tcrv'}{'frequency'} = 0.01; # 0.01% of total seqs  | possible values: 0.01% 0.02%
$CLUSTERING_THRESHOLDS->{'tcrj'}{'frequency'} = 0.1; # 0.1% of total seqs  | possible values: 0.1% 0.01% 0.02%
$CLUSTERING_THRESHOLDS->{'tcrd'}{'frequency'} = 5; # 5% of total seqs (very strict but finds ok 1 loci, to find 2nd we have to use 0.2 and we have dozens of false positives)
$CLUSTERING_THRESHOLDS->{'tcrc'}{'frequency'} = 10;
# Minimum depth ratio of a sequence with a single error compared to major sequences:
$CLUSTERING_THRESHOLDS->{'tcrv'}{'dominant_frequency'} = 20; # 20% of the major sequence depth
$CLUSTERING_THRESHOLDS->{'tcrj'}{'dominant_frequency'} = 20; # 20% of the major sequence depth
$CLUSTERING_THRESHOLDS->{'tcrd'}{'dominant_frequency'} = 20; # 20% of the major sequence depth
$CLUSTERING_THRESHOLDS->{'tcrc'}{'dominant_frequency'} = 20; # 20% of the major sequence depth

my ($INP_reads_file, $INP_outpath, $INP_nreads, $INP_direct, @INP_primers, $INP_discover, $INP_prot, $INP_align, $INP_threads, $INP_zip, $INP_test, %INP_patterns, %INP_ref_files, $INP_default, $INP_regex, $INP_debug);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s' => \$INP_reads_file,
	'o|output=s' => \$INP_outpath,
	'n|number=i' => \$INP_nreads,
	'd|discover:s' => \$INP_discover,
	'p|prot' => \$INP_prot,
# 	'a|align:s' => \$INP_align,
	'pri|primers=s{1,}' => \@INP_primers,
	'default' => \$INP_default,
	'regex' => \$INP_regex,
	'di|direct' => \$INP_direct,
	'vpat=s' => \$INP_patterns{'tcrv'},
	'jpat=s' => \$INP_patterns{'tcrj'},
	'dpat=s' => \$INP_patterns{'tcrd'},
	'cpat=s' => \$INP_patterns{'tcrc'},
	'vref=s' => \$INP_ref_files{'tcrv'},
	'jref=s' => \$INP_ref_files{'tcrj'},
	'dref=s' => \$INP_ref_files{'tcrd'},
	'cref=s' => \$INP_ref_files{'tcrc'},
	'thr|threads:i' => \$INP_threads,
	'z|zip' => \$INP_zip,
	'debug' => \$INP_debug,
	'test' => \$INP_test,
	'<>' => \&usage,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> [-default] [options]\n";
	print "\nOptions:\n";
	print "  -i <file>\tInput FASTQ or FASTA file (compressed or uncompressed).\n";
	print "  -o <name>\tOutput files prefix.\n";
	print "  -n <number>\tNumber of reads/sequences to analyze.\n";
	print "  -d <algorithm>\tDiscover allele sequences for TCR segments by clustering sequencing errors with the desired algorithm: 'amplisas' (default) or 'cdhit'.\n";
	print "  -p\t\tPrint protein sequences.\n";
# 	print "  -a <type>\tAlignment type (match (default), gassst, blast).\n";
	print "  -pri <seq1> [ <seq2> ... <seqN> ]\n\t\tPrimer sequence/s in constant segment\n\t\t(only sequences containing the primer/s will be analyzed).\n";
	print "  -default\tUse default TCR patterns.\n";
	print "  -regex\tUse Perl REGEX for patterns.\n";
	print "  -di\t\tAnalyze reads only in direct sense.\n";
# 	print "  -debug\t\tPrints additional info for debugging.\n";
# 	print "  -old\t\tPerforms old clustering.\n";
	print "  -vpat <pat>\tProtein pattern to recognize the TCR V segment (Prosite format).\n"; #      (default='".$DEFAULT_TCR_\tProtein patternS{'tcrv'}."').\n";
	print "  -jpat <pat>\tProtein pattern to recognize the TCR J segment (Prosite format).\n"; #      (default='".$DEFAULT_TCR_\tProtein patternS{'tcrj'}."').\n";
	print "  -dpat <pat>\tNucleotide pattern to recognize the TCR D segment (IUPAC format).\n"; #      (default='".$DEFAULT_TCR_\tProtein patternS{'tcrd'}."').\n";
	print "  -cpat <pat>\tProtein pattern to recognize the TCR J segment (Prosite format).\n"; #      (default='".$DEFAULT_TCR_PATTERNS{'tcrc'}."').\n";
	print "  -vref <file>\tFASTA file with TRBV reference sequences.\n";
	print "  -jref <file>\tFASTA file with TRBJ reference sequences.\n";
	print "  -dref <file>\tFASTA file with TRBD reference sequences.\n";
	print "  -cref <file>\tFASTA file with TRBC reference sequences.\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}

# Checks if a multifile with files from several samples is given as input
my $INP_multifile;
# Checks if a set of files from several samples is given as input into a compressed file
if (defined($INP_reads_file) && is_multifile($INP_reads_file)){
	$INP_multifile = $INP_reads_file;
# Prints usage help if no input file is specified
} elsif (!defined($INP_reads_file)){
	print "\nERROR: You must specify input file.\n";
	usage();
	exit;
}

# Defines default patterns to recognize TCR beta segments in bank vole (also valid for human and mouse) if no patterns are specified
my @extract_tcr_options;
if (defined($INP_default) && (defined($INP_patterns{'tcrv'}) || defined($INP_patterns{'tcrj'}) || defined($INP_patterns{'tcrd'}) || defined($INP_patterns{'tcrc'}))){
	print "\nWARNING: Default patterns will override specified ones.\n";
}
if (defined($INP_default)){
	%INP_patterns = %DEFAULT_TCR_PATTERNS;
	push(@extract_tcr_options,'default');
} elsif (!defined($INP_patterns{'tcrv'}) && !defined($INP_patterns{'tcrj'}) && !defined($INP_patterns{'tcrd'}) && !defined($INP_patterns{'tcrc'})){
	print "\nERROR: You must specify at least one pattern to recognize the TCR segment/s or specify 'default' option.\n";
	usage();
	exit;
}
# Defines TCR segment sequence clustering method
if (defined($INP_discover)){
	if ($INP_discover =~ /cd-?hit/i){
		$INP_discover = 'cdhit';
	} elsif ($INP_discover =~ /amplisas/i){
		$INP_discover = 'amplisas';
	} else {
		$INP_discover = $DEFAULT_CLUSTERING;
	}
}

print "\nRunning '$COMMAND_LINE'\n";

# Translates PROSITE patterns to REGEX
if (!defined($INP_regex) && !defined($INP_default)){
	while (my ($tcr_segment, $pattern) = each %INP_patterns) {
		if (defined($pattern)){
			$INP_patterns{$tcr_segment} = prosite_to_regex($pattern);
		}
	}
}

# Align reads in direct and reverse complementary senses
if (defined($INP_direct)){
	push(@extract_tcr_options,'noreverse');
}
# Translates TCR sequences and prints proteins
if (defined($INP_prot)){
	push(@extract_tcr_options,'protein');
}

# Defines TCR segments
# TRBD and CDR3 segments must go after TRBV and TRBJ segments, to be extracted correctly by 'extract_tcr_sequences'
print "\nTCR patterns:\n";
my @tcr_segments; #  = ('tcrv','tcrj','tcrd','tcrc', 'cdr3');
foreach my $tcr_segment ( ('tcrv','tcrj','tcrd','tcrc') ){
	if (defined($INP_patterns{$tcr_segment}) && $INP_patterns{$tcr_segment}){
		printf("\t%s: '%s'\n", uc($tcr_segment), $INP_patterns{$tcr_segment});
		push(@tcr_segments,$tcr_segment);
	}
}
if (defined($INP_patterns{'tcrv'}) && defined($INP_patterns{'tcrj'})){
	push(@tcr_segments,'cdr3');
}


# Reads reference sequences from files
my ($tcr_ref_data, $tcr_ref_seqs, $tcr_ref_headers, $tcr_ref_names, $md5_to_tcr_ref_name);
foreach my $tcr_segment (keys %INP_ref_files) {
	if (!defined($INP_ref_files{$tcr_segment})){
		next;
	}
	printf("\nReading %s references from '%s'.\n", uc($tcr_segment), $INP_ref_files{$tcr_segment});
	$tcr_ref_data->{$tcr_segment} = read_allele_file($INP_ref_files{$tcr_segment});
# 	if ($INP_discover eq 'cdhit'){
# 		($tcr_ref_seqs->{$tcr_segment}, $tcr_ref_headers->{$tcr_segment}, $tcr_ref_names->{$tcr_segment}, $md5_to_tcr_ref_name->{$tcr_segment})
# 		= parse_reference_file($INP_ref_files{$tcr_segment});
# 	}
}

# Sets output files format
my @output_files;
my $output_format;
if (!defined($INP_debug) && !defined($INP_zip)){
	$output_format = 'gzip';
}

# Saves into an array the files to process
my @file_list = ($INP_reads_file);
my ($tmp_dir_name,$tmp_dir,@outfiles);
if (defined($INP_multifile)){
	# Uncompress the files in temporal folder
	$tmp_dir_name = random_file_name();
	$tmp_dir = "/tmp/".$tmp_dir_name;
	@file_list = decompress($INP_multifile,undef,$tmp_dir);
}

while (my $file = shift @file_list){

	my $file_name = $file;

	if (!-f $file || (!is_fasta($file) && !is_fastq($file))){
		next;
	}

	# Sets name and format for output files
	my $outpath;
	if ($file =~ /.+\/([^\.]+)/ || $file =~ /([^\.]+)/){
		$outpath = $1;
	} else {
		$outpath = $file;
	}
	if (defined($INP_outpath)){
		$outpath = $INP_outpath."/".$outpath;
	}

	# Checks format and number or reads
	my ($seqs_file_format,$total_reads)
	= parse_sequence_file($file,undef,['verbose','stats']);

	# Creates output dir if it doesn't exist
	if ($outpath =~ /(.+)\//){
		if (!-d $1){
			mkdir($1);
		}
	}
	# Extracts full TCR sequences and segments
	if (defined($INP_debug)){
		push(@extract_tcr_options,'debug','verbose');
	} else {
		push(@extract_tcr_options,'verbose');
	}
	my ($md5_to_seq, $md5_to_depth, $md5_to_md5p, $md5p_to_md5s, $md5p_to_seq, $md5p_to_depth, $outfiles)
	= extract_tcr_sequences($file, \@tcr_segments, \%INP_patterns, \@INP_primers, \@extract_tcr_options, $INP_nreads, $outpath);
	push(@output_files, @$outfiles);

	# Prints the recognized TCR segment unique sequences into FASTA files
	if (!defined($INP_discover)) {

		my @outfiles = print_tcr_sequences($outpath, \@tcr_segments, $md5_to_seq, $md5_to_depth, $md5_to_md5p, $md5p_to_md5s, $md5p_to_seq, $tcr_ref_data, \@extract_tcr_options);
		push(@output_files,@outfiles);

	# Clusters with AmpliSAS and prints only the unique major TCR segment sequences
	} elsif ($INP_discover eq 'amplisas') {

		my  $md5_clustered_to_depth
		= cluster_tcr_sequences(\@tcr_segments, $md5_to_seq, $md5_to_depth, $CLUSTERING_THRESHOLDS, \@extract_tcr_options);

		my @outfiles = print_tcr_sequences($outpath, \@tcr_segments, $md5_to_seq, $md5_clustered_to_depth, $md5_to_md5p, $md5p_to_md5s, $md5p_to_seq, $tcr_ref_data, \@extract_tcr_options);
		push(@output_files,@outfiles);

	# Clusters with CD-HIT_EST and prints only the unique major TCR segment sequences
	} elsif ($INP_discover eq 'cdhit') {

		my  $md5_clustered_to_depth
		= cluster_tcr_sequences_cdhit(\@tcr_segments, $md5_to_seq, $md5_to_depth, $CLUSTERING_THRESHOLDS, \@extract_tcr_options, $outpath, $INP_threads);

		my @outfiles = print_tcr_sequences($outpath, \@tcr_segments, $md5_to_seq, $md5_clustered_to_depth, $md5_to_md5p, $md5p_to_md5s, $md5p_to_seq, $tcr_ref_data, \@extract_tcr_options);
		push(@output_files,@outfiles);

	}
}

if (defined($INP_multifile)){
	`rm -rf $tmp_dir`;
}

if (defined($INP_zip) && @output_files){
	my $outfiles = "'".join("' '",@output_files)."'";
	`zip -qjm $INP_outpath.zip $outfiles`;
	print "\nAnalysis results stored into '$INP_outpath.zip'.\n\n";
} elsif (@output_files){
	my $outfiles = "'".join("', '",@output_files)."'";
	print "\nAnalysis results stored into $outfiles.\n\n";
} else {
	print "\nThere was some error in the analysis and no results were retrieved.\n\n";
}




exit;


#################################################################################

# Reads reference sequences and names from file
sub parse_reference_file {

	my ($ref_file) = @_;
	
	if (!-e $ref_file){
		print "\t\tERROR: Reference file '$ref_file' doesn't exist.\n";
		return;
	}
	
	my ($ref_seqs,$ref_headers) = read_fasta_file($ref_file);
	
	my ($ref_names, $md5_to_ref_name);

	for (my $i=0; $i<=$#{$ref_headers}; $i++) {
		my $md5 = generate_md5($ref_seqs->[$i]);
		my $ref_name;
		if ($ref_headers->[$i] =~ /TR[AB][VJD][\d\?]+(-\d+)?\*\d+(\|[\w\s]+)?/i){
			$ref_name = $&;
		} elsif ($ref_headers->[$i] =~ /(.+?)\s?\|/){
			$ref_name = $1;
		} else {
			$ref_name = $ref_headers->[$i];
		}
		if (in_array($ref_names, $ref_name)){
			print "\t\tERROR: Reference name '$ref_name' is duplicated.\n";
		}
		push(@{$ref_names}, $ref_name);
		if (defined($md5_to_ref_name->{$md5})){
			print "\t\tERROR: Reference sequence '$ref_name' is duplicated.\n";
		}
		$md5_to_ref_name->{$md5} = $ref_name;
	}
	print '';
	
	return ($ref_seqs, $ref_headers, $ref_names, $md5_to_ref_name);

}

#################################################################################

# Extracts full TCR sequences and segments
sub extract_tcr_sequences {

	my ($seqs_file, $tcr_segments, $tcr_patterns, $tcr_primers, $options, $read_limit, $outpath) = @_;

	# Defines variablers to store TCR data
	my @output_files;
	my (@tcr_headers, @tcr_seqs, @tcr_prots);
	my (@notcr_headers, @notcr_seqs, @notcr_prots);
	my (@noframe_headers, @noframe_seqs);
	my ($md5_to_seq, $md5_to_depth, $md5_to_md5p, $md5p_to_md5s, $md5p_to_seq, $md5p_to_depth);
	my ($count_reads, $count_tcr_primer_seqs, $count_tcr_noframe_seqs, $count_tcr_valid_seqs, $count_nontcr_pattern_seqs) = (0, 0, 0, 0, 0);
	
	my ($default,$debug,$verbose,$prot,$revcomp) = (0,0,0,0,1);
	if (in_array($options,'default')){
		$default = 1;
	}
	if (in_array($options,'debug')){
		$debug = 1;
	}
	if (in_array($options,'verbose')){
		$verbose = 1;
	}
	if (in_array($options,'prot') || in_array($options,'protein')){
		$prot = 1;
	}
	if (in_array($options,'noreverse')){
		$revcomp = 0;
	}

	# Defines a unique pattern for the primers (if defined)
	my $tcr_primer_pattern;
	if (defined($tcr_primers) && @$tcr_primers){
		$tcr_primer_pattern = join('|', map regex($_), @$tcr_primers);
		if ($revcomp){
			$tcr_primer_pattern .= "|".join('|', map regex(iupac_reverse_complementary($_)), @$tcr_primers);
		}
	}
	
	# Defines the full TCR pattern
# 	my $tcr_full_pat = sprintf("(%s).+?(%s).+?(%s)", $tcr_patterns->{'tcrv'},$tcr_patterns->{'tcrj'},$tcr_patterns->{'tcrc'});
	my (@tcr_pattern_segments, %tcr_pattern_order);
	my $count_pattern = 0;
	foreach my $tcr_segment (@$tcr_segments) {
		if (defined($tcr_patterns->{$tcr_segment}) && in_array(['tcrv','tcrj','tcrc'], $tcr_segment) ){
			$count_pattern++;
			push(@tcr_pattern_segments, $tcr_patterns->{$tcr_segment});
			$tcr_pattern_order{$tcr_segment} = $count_pattern;
		}
	}
	my $tcr_full_pat = "(".join(").+?(", @tcr_pattern_segments).")";

	# Opens reads file
	if (is_gzip($seqs_file) || is_zip($seqs_file)) {
		open(READSFILE, "zcat $seqs_file |") || die "\nERROR: cannot open compressed '$seqs_file'\n\n";
	} else {
		open(READSFILE, $seqs_file) || die "\n ERROR: cannot open '$seqs_file'\n\n";
	}

	# Extracts full TCR sequences and TCR segment sequences
	print "\nExtracting full TCR sequences and segments.\n";
	my ($tcr_seq,$header);
	while (<READSFILE>) {

		# Annotates headers
		if (/^[@>](.+)\n/){
			$header = $1;
		}

		# Skip non sequence lines or TCR sequences with undefined nucleotides
		if (!defined($header) || !/^[ACGTU]+$/i){
			next;
		}

		# Counts sequences and removes line break at the end
		$count_reads++;
		$tcr_seq = $_;
		chomp $tcr_seq;

		# Checks if the sequence contains the primer sequence
		if (defined($tcr_primers) && @$tcr_primers){
			if ($tcr_seq =~ /$tcr_primer_pattern/){
				$count_tcr_primer_seqs++;
			} else {
				next;
			}
		}

		# Extracts TCR V and J segments
		# $tcr_seq = 'CACTCTATCCGACAAGCAGTGGTATCAACGCAGAGTACGCGGGGACATCAGATCTTGCCTTGGTCCTGAGATGATCTTCAGGCTCTTCTGTGTGGCCTTGAGTCTCCTGTGTGCAAAACCCATGGAGGCTGCAGTCACCCAAAGCCCGAGAAACAAGGTGACAGTAACAGGAGGAGAGGTGGAACTAAGCTGTCACCAGACGGACAACCACAACATTATGTACTGGTATCGGCAGGACCTGGGTCATGGACTGAGGCTGATCCATTACTCATATGGTGATGGCAGCACTGAGAACGGAGATGTCCCTTATGGGTACAAGGCCACCAGACCGAAGACAGAGGAGTACTCTCTCATTCTGGAGAAGGCTTCCCCCTCTCAGACAGCTGTGTACTTCTGTGCCAGCAGAGGGGACAGCCAAAATTCGCCCCTCTTCTTTGCAGCAGGCACCAGGCTCACTGTGACAGAGGATCTGAGCA';
		my $found = 0;
		my $reversed = 0;
		my $noframe = 0;
		# my $stop_codons;
		for (my $frame=0; $frame<=2; $frame++){

			my $tcr_prot = dna_to_prot(substr($tcr_seq,$frame));
			# my $tcr_prot_len = length($tcr_prot);

			# Recognises the TCR segments using a pattern of conserved residues			
			if ($tcr_prot =~ /\*/){
				$noframe++;
# 				while ($tcr_prot =~ /\*/g){
# 					push(@{$stop_codons->{$frame}}, $-[0]);
# 				}
			} elsif ($tcr_prot =~ /$tcr_full_pat/i){
			
				$found = 1;

				my (%first_aa, %last_aa, %first_nt, %last_nt);

				# If the default pattern is used, annotates specific positions before and after TRBV and TRBJ pattern matches
				# segment definitions:
				# Malissen et al (1984): TRBJ: Gly(G)-[6-7] & Val(V)+1
				# Yanagi et al (1985): TRBV: Gln(Q)+3 & Cys(C)+2(7nt) ; TRBJ: ? & Val(V)+1
				# Bolotin et al (2012): TRBV: ? & Cys(C)+2(7nt) ; TRBJ: Gly(G)-5(14nt) & ? ; CDR3: Cys(C) & Gly(G)-1(3nt)
				# IMGT: TRBV: Gln(Q)-5 & Cys(C)+[0-5] ; TRBJ: Gly(G)-[6-7] & Val(V)+1
				
				if ($default){
					# @- contains the position of the beginning of the regex matches starting from 0 (add 1 for human readable position)
					# @+ contains the position of the ending of the regex matches, the position of next char (do not add 1)
					# TRBV segment defined as the segment recognized by 1st pattern plus 5 aa before Gln(Q) and 1 aa after Cys(C)
					$first_aa{'tcrv'} = $-[$tcr_pattern_order{'tcrv'}]+1-5;
					$last_aa{'tcrv'} = $+[$tcr_pattern_order{'tcrv'}]+1;
					# TRBJ segment defined as the segment recognized by 2rd pattern plus 4 aa before Gly(G) and 1 aa after Val|Ile(V,I)
					$first_aa{'tcrj'} = $-[$tcr_pattern_order{'tcrj'}]+1-4;
					$last_aa{'tcrj'} = $+[$tcr_pattern_order{'tcrj'}]+1;
# 					# CDR3 segment defined as the segment recognized between V and J segments
# 					$first_aa{'cdr3'} = $last_aa{'tcrv'}+1;
# 					$last_aa{'cdr3'} = $first_aa{'tcrj'}-1;
					# CDR3 segment defined as the segment recognized between last Cys(C) in V segment and one res (Phe,F) before first Gly(G) in J segment pattern
					$first_aa{'cdr3'} = $+[$tcr_pattern_order{'tcrv'}];
					$last_aa{'cdr3'} = $-[$tcr_pattern_order{'tcrj'}]+1-1;
					# TRBC segment defined as the segment recognized by 3rd pattern
					$first_aa{'tcrc'} = $-[$tcr_pattern_order{'tcrc'}]+1;
					$last_aa{'tcrc'} = $+[$tcr_pattern_order{'tcrc'}];
				} else {
					if (defined($tcr_pattern_order{'tcrv'})){
						$first_aa{'tcrv'} = $-[$tcr_pattern_order{'tcrv'}]+1;
						$last_aa{'tcrv'} = $+[$tcr_pattern_order{'tcrv'}];
					}
					if (defined($tcr_pattern_order{'tcrj'})){
						$first_aa{'tcrj'} = $-[$tcr_pattern_order{'tcrj'}]+1;
						$last_aa{'tcrj'} = $+[$tcr_pattern_order{'tcrj'}];
					}
					if (defined($tcr_pattern_order{'tcrv'}) && defined($tcr_pattern_order{'tcrj'})){
						$first_aa{'cdr3'} = $+[$tcr_pattern_order{'tcrv'}];
						$last_aa{'cdr3'} = $-[$tcr_pattern_order{'tcrj'}]+1-1;
					}
					if (defined($tcr_pattern_order{'tcrc'})){
						$first_aa{'tcrc'} = $-[$tcr_pattern_order{'tcrc'}]+1;
						$last_aa{'tcrc'} = $+[$tcr_pattern_order{'tcrc'}];
					}
				}
				# If TRBV or TRBJ segment sequences are not complete
				#if ($first_aa{'tcrv'}<1 || $last_aa{'tcrj'}>length($tcr_prot)){
					#next;
				#}
				# If any segment sequence is not complete
				my $partial_seq = 0;
				foreach my $tcr_segment (keys %tcr_pattern_order){
					if ($first_aa{$tcr_segment}<1 || $last_aa{$tcr_segment}>length($tcr_prot)){
						$partial_seq = 1;
					}
				}
				if ($partial_seq) {
					next;
				}

				# Proccess all the TCR segments defined before (Variable, Diversity, Joining...)
				foreach my $tcr_segment (@$tcr_segments) {
					# Check if some segment is anormally short
					if ((defined($first_aa{$tcr_segment}) && $first_aa{$tcr_segment}<1) || (defined($last_aa{$tcr_segment}) && $last_aa{$tcr_segment}>length($tcr_prot))){
						next;
					}
					# Extracts TCR segment sequence
					my ($seq,$qual);
					if ($tcr_segment eq 'tcrv' || $tcr_segment eq 'tcrj' || $tcr_segment eq 'tcrc') {
						$first_nt{$tcr_segment} = $first_aa{$tcr_segment}*3-2+$frame;
						$last_nt{$tcr_segment} = $last_aa{$tcr_segment}*3+$frame;
						$seq = substr($tcr_seq,$first_nt{$tcr_segment}-1,$last_nt{$tcr_segment}-$first_nt{$tcr_segment}+1);
						#$seq_prot = substr($tcr_prot,$first_aa{$tcr_segment}-1,$last_aa{$tcr_segment}-$first_aa{$tcr_segment}+1);
					} elsif ($tcr_segment eq 'cdr3') {
						$first_nt{$tcr_segment} = $first_aa{$tcr_segment}*3-2+$frame;
						$last_nt{$tcr_segment} = $last_aa{$tcr_segment}*3+$frame;
						# Variable and Joining segment nucleotides will be in lowercase
						$seq = lc(substr($tcr_seq,$first_nt{'cdr3'}-1,$last_nt{'tcrv'}-$first_nt{'cdr3'}+1))
						.uc(substr($tcr_seq,$last_nt{'tcrv'},$first_nt{'tcrj'}-$last_nt{'tcrv'}-1))
						.lc(substr($tcr_seq,$first_nt{'tcrj'}-1,$last_nt{'cdr3'}-$first_nt{'tcrj'}+1));
					} elsif ($tcr_segment eq 'tcrd') {
						# TRBD segment will be between segments V and J
						my $seq_ = substr($tcr_seq,$last_nt{'tcrv'},$first_nt{'tcrj'}-$last_nt{'tcrv'});
						if ($seq_ =~ /$INP_patterns{'tcrd'}/i) {
							$seq = $&; # = substr($tcr_seq,$+[0]-$-[0]+1)
						} else {
							next;
						}
						# if ($seq_ =~ /gggac.+?(g{5,7}|gggaggg)\w/i) {
						# if ($seq_ =~ /gggac\w{7,9}/i) {
						# if ($seq_ =~ /gggacagggggc|gggactggggggg\w/i) {
						# TRBD segment defined as the segment recognized by the following pattern
# 						if ($seq_ =~ /gggaca\w{6}|gggact\w{8}/i) {
# 							$seq = $&; # = substr($tcr_seq,$+[0]-$-[0]+1)
# 						} else {
# 							next;
# 						}
					}
					# Store the sequence with a unique MD5 identifier and check that identifiers are trully unique per sequence
					my $md5 = generate_md5($seq);
# if ($md5 eq '853cd0d194baa5c035015f9b0a2e2f11'){
# print '';
# }
					if (!defined($md5_to_seq->{$tcr_segment}{$md5})){
						$md5_to_seq->{$tcr_segment}{$md5} = $seq;
					} elsif (defined($md5_to_seq->{$tcr_segment}{$md5}) && $seq ne $md5_to_seq->{$tcr_segment}{$md5}){
						print "\nERROR: Sequence MD5 signature '$md5' is not unique.\n";
						print "SEQ1: ".$md5_to_seq->{$tcr_segment}{$md5}."\n";
						print "SEQ2: ".$seq."\n\n";
						exit;
					}
					$md5_to_depth->{$tcr_segment}{$md5}++;
					if ($prot){
						my $seq_prot = dna_to_prot($seq);
						my $md5p = generate_md5($seq_prot);
						if (!defined($md5p_to_seq->{$tcr_segment}{$md5p})){
							$md5p_to_seq->{$tcr_segment}{$md5p} = $seq_prot;
						}
						if (!defined($md5p_to_md5s->{$tcr_segment}{$md5p}) || !in_array($md5p_to_md5s->{$tcr_segment}{$md5p},$md5)){
							push(@{$md5p_to_md5s->{$tcr_segment}{$md5p}},$md5);
						}
						if (!defined($md5_to_md5p->{$tcr_segment}{$md5})){
							$md5_to_md5p->{$tcr_segment}{$md5} = $md5p;
						}
						$md5p_to_depth->{$tcr_segment}{$md5p}++;
					}
				}
			}
			# If the TCR pattern is found, annotates the TCR sequence and goes to the next sequence
			if ($found) {
# # DEBUGGING
# if ($tcr_seq =~ /ggacccaacgtcttgcagttcccaagtcatcaagtgacacgtgaggggcagatggtgatcctcagttgtgaccctgtctctaatcaccagtccttctattggtataaactgatcttgggacaaaacatagagtttctggtttctttctacagtggtaaacctatggaaaagtctaagctattcgaggaggatcgattttcagtcgaaaggcaagaggattcatatttcactctgaagatacagcccacaacgccggaggactcagccgtgtacttctgt/i) {
# print "\ntcr: $tcr_seq\n";
# }
				$count_tcr_valid_seqs++;
				push(@tcr_seqs, $tcr_seq);
				push(@tcr_headers, $header);
				if ($prot){
					push(@tcr_prots, $tcr_prot);
				}
				last;
			# Checks also reverse complementary sequence
			} elsif ($revcomp && $frame == 2 && !$reversed && !$found) {
				$frame = -1;
				$reversed = 1;
				$noframe = 0;
				$tcr_seq = iupac_reverse_complementary($tcr_seq);
			# Annotates non regular TCR sequences
			} elsif ($frame == 2 && !$found) {
				if ($noframe == 3){
					$count_tcr_noframe_seqs++;
					if ($debug){
						push(@noframe_seqs, $tcr_seq);
						push(@noframe_headers, $header);
					}
				} else {
					$count_nontcr_pattern_seqs++;
					if ($debug){
						push(@notcr_seqs, $tcr_seq);
						push(@notcr_headers, $header);
					}
				}
	# 			push(@notcr_prots, $tcr_prot);
			}
		}

		# Prints the number of processed reads and TCR seqs found in debug mode
		if ($verbose && $count_reads % 100000 == 0) {
			if ($count_reads == 100000) {
				if (defined($tcr_primers) && @$tcr_primers){
					printf("\n\t%s\t%s\t%s\t%s\t%s\n", 'Reads', 'Valid TCRs', 'Primer TCRs', 'Non pattern TCRs', 'Non frame TCRs');
				} else {
					printf("\n\t%s\t%s\t%s\t%s\n", 'Reads', 'Valid TCRs', 'Non pattern TCRs', 'Non frame TCRs');
				}
			}
			if (defined($tcr_primers) && @$tcr_primers){
				printf("\t%d\t%d\t%d\t%d\t%d\n", $count_reads, $count_tcr_valid_seqs, $count_tcr_primer_seqs, $count_nontcr_pattern_seqs, $count_tcr_noframe_seqs);
			} else {
				printf("\t%d\t%d\t%d\t%d\n", $count_reads, $count_tcr_valid_seqs, $count_nontcr_pattern_seqs, $count_tcr_noframe_seqs);
			}
		}

		# Cleans header and sequence variables
		$header = undef;
		$tcr_seq = undef;
		
		# Finish if a desired number of reads is reached
		if (defined($read_limit) && $read_limit == $count_reads){
			last;
		}

	}
	if ($verbose) {
		print "\n";
	}

	# Creates a FASTA file with TCR found sequences
	printf("\t%d valid sequences with ORF matching TCR patterns.\n", $count_tcr_valid_seqs);
	if (defined($tcr_primers) && @$tcr_primers){
		printf("\t%d sequences matching TCR primer/s.\n", $count_tcr_primer_seqs);
	}
	if ($debug){
		printf("\t%d sequences without ORF.\n", $count_tcr_noframe_seqs);
		printf("\t%d sequences with ORF not matching TCR patterns.\n", $count_nontcr_pattern_seqs);
	}
	# printf("\t%d Non ORF TCR sequences.\n", $count_tcr_nonframe_seqs);
	print "\nStoring full TCR sequences and segments.\n";
	my $outfile = create_fasta_file(\@tcr_seqs,\@tcr_headers,"$outpath.tcr.fa",$output_format);
	push(@output_files, $outfile);
	printf("\n\t%d TCR sequences stored into '%s'.\n", scalar @tcr_seqs, $outfile);
	if ($prot){
		$outfile = create_fasta_file(\@tcr_prots,\@tcr_headers,"$outpath.tcr.prot.fa",$output_format);
		push(@output_files, $outfile);
		printf("\t%d sequences matching TCR patterns stored into '%s'.\n", scalar @tcr_prots, $outfile);
	}
	if ($debug){
		$outfile = create_fasta_file(\@noframe_seqs,\@noframe_headers,"$outpath.noframe.fa",$output_format);
		push(@output_files, $outfile);
		printf("\t%d sequences without ORF stored into '%s'.\n", scalar @noframe_seqs, $outfile);
		$outfile = create_fasta_file(\@notcr_seqs,\@notcr_headers,"$outpath.notcr.fa",$output_format);
		push(@output_files, $outfile);
		printf("\t%d sequences with ORF not matching TCR patterns stored into '%s'.\n", scalar @notcr_seqs, $outfile);
	}
# 	print "\n";
	
	return ($md5_to_seq, $md5_to_depth, $md5_to_md5p, $md5p_to_md5s, $md5p_to_seq, $md5p_to_depth, \@output_files);
}

#################################################################################

# Prints the recognized TCR segment unique sequences into FASTA files
sub print_tcr_sequences {

	my ($outfile_name, $tcr_segments, $md5_to_seq, $md5_to_depth, $md5_to_md5p, $md5p_to_md5s, $md5p_to_seq, $tcr_ref_data, $options) = @_;
	
	my (@output_files, $md5_to_name);

	my ($default,$debug,$verbose,$prot) = (0,0,0,0);
	if (in_array($options,'default')){
		$default = 1;
	}
	if (in_array($options,'debug')){
		$debug = 1;
	}
	if (in_array($options,'verbose')){
		$verbose = 1;
	}
	if (in_array($options,'prot') || in_array($options,'protein')){
		$prot = 1;
	}

	# Print unique TCR segment sequences
	print "\nPrinting TCR sequences by segment.\n";
	my $depths;
	my $md5_to_tcr_ref_name;
	foreach my $tcr_segment (@$tcr_segments) {
	
		my $output_format_;
 		if ($tcr_segment eq 'cdr3' && defined($output_format)) {
			$output_format_	= $output_format;
		}

		# Matches reference sequences
		if (defined($tcr_ref_data->{$tcr_segment})){
			$md5_to_tcr_ref_name->{$tcr_segment} = match_alleles($tcr_ref_data->{$tcr_segment},$md5_to_seq->{$tcr_segment},undef,$INP_allele_align,$INP_threads);
		}

	# 	printf("Saving %s sequences.\n", uc($tcr_segment));
		my (@names, @seqs, %depths);
		# Sort TRBV sequences by number of reads
		my @sorted = sort { $md5_to_depth->{$tcr_segment}{$b} <=> $md5_to_depth->{$tcr_segment}{$a} } keys %{$md5_to_depth->{$tcr_segment}};
		my $zero_places = length(scalar @sorted);
		my ($count, $depth) = (0,0);
		foreach my $md5 (@sorted) {
# if ($md5 eq '853cd0d194baa5c035015f9b0a2e2f11'){
# print '';
# }
			$count++;
			# Annotate alleles from input files, if the sequence is not in the file, skip
			if (!defined($md5_to_tcr_ref_name->{$tcr_segment}) || !defined($md5_to_tcr_ref_name->{$tcr_segment}{$md5})){
				$md5_to_tcr_ref_name->{$tcr_segment}{$md5} = sprintf("%s-%0${zero_places}d", uc($tcr_segment), $count);
			}
			my $name = $md5_to_tcr_ref_name->{$tcr_segment}{$md5};
			$depths->{$tcr_segment}{$name} = $md5_to_depth->{$tcr_segment}{$md5};
			$depth += $md5_to_depth->{$tcr_segment}{$md5};
			$md5_to_name->{$tcr_segment}{$md5} = $name;
			push(@names, sprintf("%s | hash=%s | reads=%d", $name, $md5, $md5_to_depth->{$tcr_segment}{$md5}));
			# push(@names, sprintf("%s-%d", $md5, $md5_to_depth->{$tcr_segment}{$md5}));
			push(@seqs, $md5_to_seq->{$tcr_segment}{$md5});
		}
		my $outfile = create_fasta_file(\@seqs,\@names,"$outfile_name.$tcr_segment.fa",$output_format_);
		push(@output_files,$outfile);
		printf("\t%d %s unique sequences (total=%d) stored into '%s'.\n", $count, uc($tcr_segment), $depth, $outfile);
		if ($prot) {
			undef(@names);
			undef(@seqs);
			undef(%depths);
			($count, $depth) = (0,0);
			my (@psorted, $md5p_to_names, $md5p_to_depth);
			foreach my $md5 (@sorted) {
				my $md5p = $md5_to_md5p->{$tcr_segment}{$md5};
				if (!in_array(\@psorted, $md5_to_md5p->{$tcr_segment}{$md5})){
					push(@psorted, $md5p);
				}
				# Annotate alleles from input files, if the sequence is not in the file, skip
				push(@{$md5p_to_names->{$md5p}}, $md5_to_tcr_ref_name->{$tcr_segment}{$md5});
				$md5p_to_depth->{$md5p} += $md5_to_depth->{$tcr_segment}{$md5};
			}
			foreach my $md5p (@psorted) {
				$count++;
				my $name = join('|',@{$md5p_to_names->{$md5p}});
				$depth += $md5p_to_depth->{$md5p};
				push(@names, sprintf("%s | hash=%s | reads=%d", $name, $md5p, $md5p_to_depth->{$md5p}));
				push(@seqs, $md5p_to_seq->{$tcr_segment}{$md5p});
			}
			$outfile = create_fasta_file(\@seqs,\@names,"$outfile_name.$tcr_segment.prot.fa",$output_format_);
			push(@output_files,$outfile);
			printf("\t%d %s unique proteins (total=%d) stored into '%s'.\n", $count, uc($tcr_segment), $depth, $outfile);
		}
	}
	print "\n";

	# Annotates how many TCR sequences match the reference sequences
	if (defined($tcr_ref_data)){
		print "\nMatching TCR references by segment.\n";
		foreach my $tcr_segment (@$tcr_segments) {
			# Print allele statistics
			if (defined($tcr_ref_data->{$tcr_segment})){
				my $tcr_ref_content = '';
				my ($count_matches, $count_unique_matches) = (0,0);
				foreach my $tcr_ref_name (sort {$a cmp $b} keys %{$tcr_ref_data->{$tcr_segment}}) {
					if (defined($depths->{$tcr_segment}{$tcr_ref_name})){
						$tcr_ref_content .= sprintf("%s\t%d\n", $tcr_ref_name, $depths->{$tcr_segment}{$tcr_ref_name});
						$count_unique_matches++;
						$count_matches += $depths->{$tcr_segment}{$tcr_ref_name};
					} # else {
						# $tcr_ref_content .= sprintf("%s\t%d\n", $tcr_ref_name, 0);
					# }
				}
				$tcr_ref_content .= "\n";
				my $outfile = write_to_file("$outfile_name.$tcr_segment.txt", $tcr_ref_content);
				push(@output_files,$outfile);
				printf("\t%d %s unique sequences (total=%d) matching references annotated into '%s'.\n", $count_unique_matches, uc($tcr_segment), $count_matches, $outfile);
			}
		}
		print "\n";
	}
	
	return @output_files;

# 	if ($prot){
# 		print "\nTranslating unique TCR sequences by segment.\n";
# 		foreach my $tcr_segment (@$tcr_segments) {
# 			my (@names, @seqs);
# 			# Sort TRBV prot sequences by number of reads
# 			my @sorted = sort { $md5p_to_depth->{$tcr_segment}{$b} <=> $md5p_to_depth->{$tcr_segment}{$a} } keys %{$md5p_to_depth->{$tcr_segment}};
# 			my ($count, $depth) = (0,0);
# 			foreach my $md5p (@sorted) {
# 				my $name = 'a';
# 				$count++;
# # 				my @names_ = nsort(map $md5_to_name->{$tcr_segment}{$_}, keys %{$md5p_to_md5s->{$tcr_segment}{$md5p}});
# # 				if ($#names_ > 0){
# # 					$name = join("|",@names_);
# # 				} else {
# # 					$name = '';
# # 				}
# 				$depth += $md5p_to_depth->{$tcr_segment}{$md5p};
# 				push(@names, sprintf("%s | hash=%s | reads=%d", $name, $md5p, $md5p_to_depth->{$tcr_segment}{$md5p}));
# 				# push(@names, sprintf("%s-%d", $md5p, $md5p_to_depth->{$tcr_segment}{$md5p}));
# 				push(@seqs, $md5p_to_seq->{$tcr_segment}{$md5p});
# 			}
# 			my $outfile = create_fasta_file(\@seqs,\@names,"$outfile_name.$tcr_segment.prot.fa",$output_format_);
# 			push(@output_files,$outfile);
# 			printf("\t%d %s unique protein sequences (total=%d) stored into '%s'.\n", $count, uc($tcr_segment), $depth, $outfile);
# 			if (defined($debug)){
# 				`mafft --quiet $outfile > $outfile_name.$tcr_segment.prot.mafft.fa`;
# 			}
# 		}
# 		print "\n";
# 	}

	
}


#################################################################################

# Clusters and prints unique major TCR segment sequences
sub cluster_tcr_sequences {

	my ($tcr_segments, $md5_to_seq, $md5_to_depth, $clustering_thresholds, $options) = @_;

	my ($md5_clustered_to_depth, $md5p_clustered_to_depth);
	
	my ($debug,$verbose) = (0,0);
	if (in_array($options,'debug')){
		$debug = 1;
	}
	if (in_array($options,'verbose')){
		$verbose = 1;
	}

	print "\nClustering TCR sequences by segment.\n";

	foreach my $tcr_segment (@$tcr_segments) {
		# print "$tcr_segment: ";

		# Skips CDR3 segment
		if (!in_array(['tcrv','tcrj','tcrd','tcrc'], $tcr_segment)){
			$md5_clustered_to_depth->{$tcr_segment} = $md5_to_depth->{$tcr_segment};
			next;
		}

		# Sorts TCR segment sequences by number of reads
		my @sorted_seqs = sort { $md5_to_depth->{$tcr_segment}{$b} <=> $md5_to_depth->{$tcr_segment}{$a} } keys %{$md5_to_depth->{$tcr_segment}};

		# Counts number of sequences
		my $count_total_seqs;
		map $count_total_seqs += $md5_to_depth->{$tcr_segment}{$_} , @sorted_seqs;

		# Removes low depth sequences
		my (@seqs, @md5s);
		my $clustered_seqs = 0;
		foreach my $md5 (@sorted_seqs) {
# if ($md5 eq '853cd0d194baa5c035015f9b0a2e2f11'){
# print '';
# }
			if (100*$md5_to_depth->{$tcr_segment}{$md5}/$count_total_seqs < $clustering_thresholds->{$tcr_segment}{'frequency'} ){
				last;
			}
			push(@seqs, $md5_to_seq->{$tcr_segment}{$md5});
			push(@md5s, $md5);
			$clustered_seqs++;
		}

		# Cluster sequences
		my $count_seqs = 0;
		my @clustered_md5s;
		while (my $seq = shift @seqs) {

			my $md5 = shift @md5s;
			my $length = length($seq);
			$md5_clustered_to_depth->{$tcr_segment}{$md5} = $md5_to_depth->{$tcr_segment}{$md5};

			# Finish if last sequence found
			if ($#md5s <= 0) {
				$md5_clustered_to_depth->{$tcr_segment}{$md5} = $md5_to_depth->{$tcr_segment}{$md5};
				last;
			}

			# Aligns against the sequence all the other ones
			my $aligned_seqs = align_seqs2one($seq,$md5,\@seqs,\@md5s,'needleall');
			
			$count_seqs++;

			my @dominant_seqs; # Stores sequences very similar and with high frequencies
			for (my $i=0; $i<=$#md5s; $i++){
				my $md5_ = $md5s[$i];
# if ($md5 eq '853cd0d194baa5c035015f9b0a2e2f11' || $md5_ eq '853cd0d194baa5c035015f9b0a2e2f11'){
# print '';
# }
				# Skips if the pairwise global alignment fails
				if (!defined($aligned_seqs->{$md5_}[0]) || !defined($aligned_seqs->{$md5_}[1])){
					next;
				}

				# Skips seqs with lower identity than threshold
				my ($identical,$total) = binary_score_nts($aligned_seqs->{$md5_}[0],$aligned_seqs->{$md5_}[1]); #time=10%
				if ($identical/$length*100 < $clustering_thresholds->{$tcr_segment}{'identity'}){
					next;
				}

				# Skips if the sequence is more similar to another dominant
				my $skip_seq = 0;
				foreach my $dominant_seq (@dominant_seqs){
					my ($identical_,$total_) = binary_score_nts(align_2seqs($dominant_seq,$md5_to_seq->{$tcr_segment}{$md5_}, 'needle'));
					if ($identical_>$identical){
						$skip_seq = 1;
						last;
					}
				}
				if ($skip_seq) {
					next;
				}

				# Skips high depth/frequency sequences (before clustering) compared or not to the dominant frequency and process them in next clustering round
				if (100*$md5_to_depth->{$tcr_segment}{$md5_}/$md5_to_depth->{$tcr_segment}{$md5} >= $clustering_thresholds->{$tcr_segment}{'dominant_frequency'}){
					push(@dominant_seqs,$md5_to_seq->{$tcr_segment}{$md5_});
					next;
				}

				# If the variant passes the previous conditions, cluster it with the dominant one
				$md5_clustered_to_depth->{$tcr_segment}{$md5} += $md5_to_depth->{$tcr_segment}{$md5_};
				splice(@md5s,$i,1);
				splice(@seqs,$i,1);
				$i--;

				$count_seqs++;
				if ($verbose && $count_seqs % 10000 == 0) {
					if ($count_seqs == 10000) {
						printf("\n\t%s:\n\t%s\t%s\n", uc($tcr_segment), 'Clusters', 'Total' );
					}
					printf("\t%d\t%d\n", scalar keys %{$md5_clustered_to_depth->{$tcr_segment}}, $count_seqs);
				}
			}
			print '';
		}
		if ($verbose) {
			printf("\t%d clusters from %d %s unique sequences (total=%d, low_depth=%d).\n", scalar keys %{$md5_clustered_to_depth->{$tcr_segment}}, $clustered_seqs, uc($tcr_segment), scalar @sorted_seqs, scalar @sorted_seqs - $clustered_seqs);
		}
	}

	return $md5_clustered_to_depth;


}

#################################################################################

# Cluster and print unique major TCR segment sequences
sub cluster_tcr_sequences_cdhit {

	my ($tcr_segments,$md5_to_seq,$md5_to_depth,$clustering_thresholds,$options,$outpath,$threads) = @_;

	my ($md5_clustered_to_depth, $md5p_clustered_to_depth);

	my ($debug,$verbose) = (0,0);
	if (in_array($options,'debug')){
		$debug = 1;
	}
	if (in_array($options,'verbose')){
		$verbose = 1;
	}

	print "\nClustering TCR sequences by segment.\n";

	foreach my $tcr_segment (@$tcr_segments) {
		# print "$tcr_segment: ";

		# Skips CDR3 segment
		if (!in_array(['tcrv','tcrj','tcrd','tcrc'], $tcr_segment)){
			$md5_clustered_to_depth->{$tcr_segment} = $md5_to_depth->{$tcr_segment};
			next;
		}

		# Sort TCR segment sequences by number of reads
		my @sorted_seqs = sort { $md5_to_depth->{$tcr_segment}{$b} <=> $md5_to_depth->{$tcr_segment}{$a} } keys %{$md5_to_depth->{$tcr_segment}};

		# Counts number of sequences
		my $count_total_seqs;
		map $count_total_seqs += $md5_to_depth->{$tcr_segment}{$_} , @sorted_seqs;

		# Removes low depth sequences
		my (@seqs, @md5s);
		my $clustered_seqs = 0;
		foreach my $md5 (@sorted_seqs) {
			if (100*$md5_to_depth->{$tcr_segment}{$md5}/$count_total_seqs < $clustering_thresholds->{$tcr_segment}{'frequency'} ){
				last;
			}
			push(@seqs, $md5_to_seq->{$tcr_segment}{$md5});
			push(@md5s, $md5);
			$clustered_seqs++;
		}

		# Clusters similar sequences
		my $clusters;
		my $dump_file = "$outpath.$tcr_segment.cluster.dump";
		if (defined($INP_test) && -e $dump_file){
			print "\tRecovering cluster data from '$dump_file'.\n";
			$clusters = recover_data_dump($dump_file, 'Storable');
		} else {
			my $clustering_params;
			if (defined($INP_threads) && $INP_threads>1){
				$clustering_params = "cd-hit-est -g 1 -M 0 -T $INP_threads";
			} else {
				$clustering_params = "cd-hit-est -g 1 -M 0";
			}
			$clusters = cluster_seqs(\@seqs,\@md5s,$clustering_params,$clustering_thresholds->{$tcr_segment}{'identity'}/100);
			if (defined($INP_test) && defined($clusters) && !-e $dump_file) {
				print "\tDumping cluster data into '$dump_file'.\n";
				store_data_dump($clusters, $dump_file, 'Storable');
			}
		}

		# Sort clusters by number of reads
		my %count_cluster_seqs;
		for (my $i=0; $i<=$#{$clusters}; $i++) {
			for (my $j=0; $j<=$#{$clusters->[$i]}; $j++) {
				# $clusters->[$i][$j]{'name'} =~ /hash=(\w{32})/;
				# my $md5 = $1;
				my $md5 = $clusters->[$i][$j]{'name'};
				$clusters->[$i][$j]{'md5'} = $md5;
				$clusters->[$i][$j]{'depth'} = $md5_to_depth->{$tcr_segment}{$md5};
				$count_cluster_seqs{$i} += $md5_to_depth->{$tcr_segment}{$md5};
			}
		}

		my @sorted_clusters = sort { $count_cluster_seqs{$b} <=> $count_cluster_seqs{$a}} keys %count_cluster_seqs;
		my $count_clusters = 0;
		my $zero_places = length(scalar @sorted_clusters);
		# Loop the clusters and annotate the cluster members with higher depth
		foreach my $i (@sorted_clusters) {
			$count_clusters++;
			# Sort cluster members by number of reads
			my @sorted_cluster_members = sort { $clusters->[$i][$b]{'depth'} <=> $clusters->[$i][$a]{'depth'}} 0 .. $#{$clusters->[$i]};
			# Threshold for sequences that differ by only 1 error with the major one
			my $min_single_error_seq_depth = $clusters->[$i][$sorted_cluster_members[0]]{'depth'}*$clustering_thresholds->{$tcr_segment}{'dominant_frequency'}/100;
			# Print only the major sequences in cluster
			foreach my $j (@sorted_cluster_members){
				# Skip sequences with very low depth or low depth compared to the major one
	# 			if ($clusters->[$i][$j]{'depth'} < $threshold_depth1){ # || $clusters->[$i][$j]{'depth'} < $clusters->[$i][$sorted_cluster_members[0]]{'depth'}/20){
	# 				last;
	# 			}
				my $single_error = 0;
				my $md5 = $clusters->[$i][$j]{'md5'};
				# Do not annotate single errors with not very high depth
				foreach my $seq (@seqs) {
					my ($error_type, $error_pos) = compare_sequences($md5_to_seq->{$tcr_segment}{$md5}, $seq);
					if (defined($error_type) && $clusters->[$i][$j]{'depth'} < $min_single_error_seq_depth){
						$single_error = 1;
						last;
					}
				}
				if ($single_error) {
					next;
				}
				# If the cluster member passes the previous conditions, annotates it as a real variant
				$md5_clustered_to_depth->{$tcr_segment}{$md5} += $md5_to_depth->{$tcr_segment}{$md5};
# 				$count_cluster_members++;
# 				push(@names, sprintf("%s-%0${zero_places}d*%02d | hash=%s | reads=%d", uc($tcr_segment), $count_clusters, $count_cluster_members, $md5, $md5_to_depth->{$tcr_segment}{$md5}));
# 				push(@seqs, $md5_to_seq->{$tcr_segment}{$md5});
			}
			print '';
		}
		if ($verbose) {
			printf("\t%d clusters from %d %s unique sequences (total=%d, low_depth=%d).\n", scalar keys %{$md5_clustered_to_depth->{$tcr_segment}}, $clustered_seqs, uc($tcr_segment), scalar @sorted_seqs, scalar @sorted_seqs - $clustered_seqs);
		}
# 		printf("\t%d %s clusters processed.\n", scalar @sorted_clusters, uc($tcr_segment));
	}

	return $md5_clustered_to_depth;

}

#################################################################################

# # OBSOLETE:
# 
# # Cluster and print unique major TCR segment sequences
# sub cluster_tcr_segment_seqs_old {
# 
# 	my ($tcr_segment,$md5_to_seq,$md5_to_depth,$clustering_thresholds,$options,$outpath) = @_;
# 
# 	# Sort TCR segment sequences by number of reads
# 	my @sorted_seqs = sort { $md5_to_depth->{$tcr_segment}{$b} <=> $md5_to_depth->{$tcr_segment}{$a} } keys %{$md5_to_depth->{$tcr_segment}};
# 
# 	my $count_total_seqs;
# 	map $count_total_seqs += $md5_to_depth->{$tcr_segment}{$_} , @sorted_seqs;
# 
# 	# To have sequences ordered by depth helps in the clustering (major seqs will be the firsts of cluster)
# 	my (@names_, @seqs_);
# 	# my ($seqs_, $names_) = sequences_hash_to_array($md5_to_seq->{$tcr_segment});
# 	my $count_low_depth_seqs = 0;
# 	foreach my $md5 (@sorted_seqs) {
# 		# Skip single sequences with very low depth
# 		if (100*$md5_to_depth->{$tcr_segment}{$md5}/$count_total_seqs <= $clustering_thresholds->{'frequency'}) {
# 			$count_low_depth_seqs ++;
# 			next;
# 		}
# 		push(@names_, $md5);
# 		# push(@names, sprintf("%s-%d", $md5, $md5_to_depth->{$tcr_segment}{$md5}));
# 		push(@seqs_, $md5_to_seq->{$tcr_segment}{$md5});
# 	}
# 
# 	# Arrays to store unique major TCR segment sequences
# 	my ($clusters, @names, @seqs, @quals);
# 
# 	printf("Clustering %d %s unique sequences (unique=%d, low_depth=%d).\n", scalar @names_, uc($tcr_segment), scalar @sorted_seqs, $count_low_depth_seqs);
# 	my $dump_file = "$outpath.$tcr_segment.cluster.dump";
# 	if (defined($INP_test) && -e $dump_file){
# 		print "\tRecovering cluster data from '$dump_file'.\n";
# 		$clusters = recover_data_dump($dump_file, 'Storable');
# 	} else{
# 		my $clustering_params;
# 		if (defined($INP_threads) && $INP_threads>1){
# 			$clustering_params = "cd-hit-est -g 1 -M 0 -T $INP_threads";
# 		} else {
# 	# BORRAR THREADS
# 			$clustering_params = "cd-hit-est -g 1 -M 0 -T 30";
# 		}
# 		$clusters = cluster_seqs(\@seqs_,\@names_,$clustering_params,$clustering_thresholds->{'identity'}/100);
# 		if (defined($INP_test) && defined($clusters) && !-e $dump_file) {
# 			print "\tDumping cluster data into '$dump_file'.\n";
# 			store_data_dump($clusters, $dump_file, 'Storable');
# 		}
# 	}
# 
# 
# 	# 	store_data_dump($clusters,'kk.dump');
# # 	my $clusters = recover_data_dump('kk.dump');
# 
# 	# Sort clusters by number of reads
# 	my %count_cluster_seqs;
# 	for (my $i=0; $i<=$#{$clusters}; $i++) {
# 		for (my $j=0; $j<=$#{$clusters->[$i]}; $j++) {
# 			# $clusters->[$i][$j]{'name'} =~ /hash=(\w{32})/;
# 			# my $md5 = $1;
# 			my $md5 = $clusters->[$i][$j]{'name'};
# 			$clusters->[$i][$j]{'md5'} = $md5;
# 			$clusters->[$i][$j]{'depth'} = $md5_to_depth->{$tcr_segment}{$md5};
# 			$count_cluster_seqs{$i} += $md5_to_depth->{$tcr_segment}{$md5};
# 		}
# 	}
# 
# 	my @sorted_clusters = sort { $count_cluster_seqs{$b} <=> $count_cluster_seqs{$a}} keys %count_cluster_seqs;
# 	my $count_clusters = 0;
# 	my $zero_places = length(scalar @sorted_clusters);
# 	# Loop the clusters and annotate the cluster members with higher depth
# 	foreach my $i (@sorted_clusters) {
# # if ($i == 122){
# # print '';
# # }
# 		$count_clusters++;
# 
# 		# Sort cluster members by number of reads
# 		my @sorted_cluster_members = sort { $clusters->[$i][$b]{'depth'} <=> $clusters->[$i][$a]{'depth'}} 0 .. $#{$clusters->[$i]};
# 		# Threshold for sequences that differ by only 1 error with the major one
# 		my $min_single_error_seq_depth = $clusters->[$i][$sorted_cluster_members[0]]{'depth'}*$clustering_thresholds->{'dominant_frequency'}/100;
# 		# Print only the major sequences in cluster
# 		my $count_cluster_members = 0;
# 		foreach my $j (@sorted_cluster_members){
# 			# Skip sequences with very low depth or low depth compared to the major one
# # 			if ($clusters->[$i][$j]{'depth'} < $threshold_depth1){ # || $clusters->[$i][$j]{'depth'} < $clusters->[$i][$sorted_cluster_members[0]]{'depth'}/20){
# # 				last;
# # 			}
# 			my $single_error = 0;
# 			my $md5 = $clusters->[$i][$j]{'md5'};
# 			# Do not annotate single errors with not very high depth
# 			foreach my $seq (@seqs) {
# 				my ($error_type, $error_pos) = compare_sequences($md5_to_seq->{$tcr_segment}{$md5}, $seq);
# 				if (defined($error_type) && $clusters->[$i][$j]{'depth'} < $min_single_error_seq_depth){
# 					$single_error = 1;
# 					last;
# 				}
# 			}
# 			if ($single_error) {
# 				next;
# 			}
# 			$count_cluster_members++;
# 			push(@names, sprintf("%s-%0${zero_places}d*%02d | hash=%s | reads=%d", uc($tcr_segment), $count_clusters, $count_cluster_members, $md5, $md5_to_depth->{$tcr_segment}{$md5}));
# 			push(@seqs, $md5_to_seq->{$tcr_segment}{$md5});
# 		}
# 		print '';
# 	}
# 	printf("\t%d %s clusters processed.\n", scalar @sorted_clusters, uc($tcr_segment));
# 	
# 	return(\@seqs,\@names,\@quals);
# 
# }

#################################################################################


























