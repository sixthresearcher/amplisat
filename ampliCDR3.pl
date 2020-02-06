#!/usr/bin/perl -w

# MODIFIED VERSION FOR RACE+UMIs PROTOCOL

# File R1 will contain: Spacer (0-6bp) + RACE primer + UMI (15bp) + RACE primer + Partial V region
# R1: {N}(0-6)AAGCAGTGGTATCAACGCAGAGT{N}(15)CTTGGGGG....
# File R2 will contain: Spacer (0-6bp) + Barcode of sample (8bp) + Primer C region + CDR3 + Small part V region
# R2: {N}(0-6){GGACTCCT}GGACTCACCTTGCTCAGATCCT...

# Extracts TCR CDR3 sequences trimmed based on primer location
# Prints statistics of CDR3 sequence lengths
#
# Examples:
# Analysis of CDR3 sequences without clustering:
#    perl ampliCDR3.pl -cl 0 -i TCR_concat.fastq.gz -d primers.csv -o outfolder


my $VERSION = "v2.0";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Extracts TCR CDR3 sequences trimmed based on primer location.";

# Modules are in folder '../' in the path of the script
use lib "lib";
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;
use Sort::Naturally;
use List::Util qw(shuffle);
use IO::Compress::Gzip qw(gzip $GzipError) ;
use Time::HiRes qw(gettimeofday);
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Statistics::Descriptive;
use Data::Dumper;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options
# By default analyzes the sequences only in the original orientation
my $INP_revcomp = 0;
# Default sequence output format
my $INP_output_format = 'fasta';
# Maximum number of allowed primer+barcode sequences
# my $MAX_PRIMER_BARCODE_SEQS = 10000;
# Extract additional nucleotides to both sides of CDR3 region (to check that they are conserved and CDR3 is true)
my $INP_extra_length = 3; # 3 = MiXCR results
# Minimum CDR3 length
my $INP_minlen = 15;
# Maximum CDR3 length
my $INP_maxlen = 63;
# Filter CDR3 sequences not in-frame
my $INP_noframe = 1;
# Number of substitutions to cluster variants within UMIs
# If $INP_umi_errors = undef, then all variants within UMI will be clustered in 1
my $INP_umi_errors = 2;
# Number of substitutions to cluster variants together between UMIs
my $INP_cluster_errors = 0;
# Keep singletons
my $INP_singletons = 0;
# To print the options
my %yes_no = ( 0 => 'no', 1 => 'yes' );
# Steps to print statistics in Excel file;
my @analysis_steps_to_print = ('clustered');
# my @analysis_steps_to_print = ('filtered', 'clustered');
# my $TCRB_CDR3_PATTERN = 'TG[TC]([GA][CG]\w+)'; # Bank vole article
# my $TCRB_CDR3_PATTERN = 'TA[TC]\w{3}TG[TC]([GA][CG]\w+)\w{31}AGGACCTGA'; # Human
# my $TCRB_CDR3_PATTERN = 'T[AG][TC]\w{3}TG[TC]([GA][CG]\w+)\w{31}AGGA[CT]CTGA'; # Human + Mouse + Vole improved
# my $TCRA_CDR3_PATTERN = 'TA[TC]\w[TA][CT]TG[TC]([GA]\w+)\w{34}ATATCCAGAA'; # Human
# my $TCRA_CDR3_PATTERN = 'TA[TC]\w[TA][CT]TG[TC]([GA]\w+)\w{34}A[TC]ATCCAGAA'; # Human + Mouse
# my $CDR3_PARTIAL_PATTERN ='T[AG][TC]\w{3}TG[TC]([GA]\w+)'; # Alpha+Beta, Human+Mouse+Vole
my $CDR3_GENERIC_PATTERN ='T[AG][TC]\w{3}TG[TC]([GA]\w+?)\w{34}A[CT][CA]T(GA|CC)[AG]\wAA'; # Alpha+Beta, Human+Mouse+Vole

# The patterns with the J region do not work correctly, the patterns contain repeated sequences in the genes
my $CDR3_PATTERNS = {
	'human' => {	'alpha' => 'TA[TC]\w[TA][CT]TG[TC](\w{5,})\w{34}ATATCCAGAA', # 'TA[TC]\w[TA][CT]TG[TC](\w+?)(TTT|TTC|TGG|TGC)G[GC]\w{4}G[GA][AT]'
			#'TRAV37' => 'TTCTTCTGC(\w{5,})\w{34}ATATCCAGAA', # 'TRAV37' is the unique allele that doesn't match the previous pattern (it's not included in MiXCR)
			'beta'  => 'TA[TC]\w{3}[TC][GA][TC]([AG][GC]\w{3,})\w{31}AGGACCTGA', # 'TA[TC]\w{3}[TC][GA][TC]([AG][GC]\w+?)(TTT|TTC|GTC|TGG)GG\w{4}GG'
			#'TRBV24OR9-2' => 'TACTTCAGT([AG][GC]\w{3,})\w{31}AGGACCTGA', # 'TRBV24OR9-2' is the unique allele that doesn't match the previous pattern
			'gamma' => '[TA]A[TC][TC]ACTG[TC](\w{5,})\w{34}AACAACTTGA', # '[TA]A[TC][TC]ACTG[TC](\w+?)TTTG[GC]\w{4}GG[AG]AC[AT]A'
			'delta' => 'TACT[AT][CT]TGT(GC\w{3,})\w{33}AGAAGTCAG' }, # 'TACT[AT][CT]TGT(GC\w+?)TT[TC]GG[AC]A[AC]\wGG[AC]A'

	# 'monkey' => {	'beta'  => 'T\w[TC]\wT\w[TC]T\wTG[TC]([GCA][CG]\w{3,})(TTT|TTC|CTG)GG\w{4}GG',

	'mouse' => {	'alpha' => 'TA[TC]\w[TA][CT]TG[TC](\w{5,})\w{34}ACATCCAGAA', # 'TA[TC]\w[TA][CT]TG[TC](\w+?)((TTT|TTC|TTG|TTA|TCT|GTT)GG|GTTGA|TTTGC|TTTAG)\w{4}(GG|GA|GT|TG|AG|AT)'
			'beta'  => 'T[AT]\w[TC][TA][TCG][TG]G[TCG]([GAT][CG]\w{3,})\w{31}AGGATCTGA', # 'T[AT]\w[TC][TA][TCG][TG]G[TCG]([GAT][CG]\w+?)((TTT|TTC)GG|TCTTG|TTTGC|GATTG)\w{4}(GG|CC)'
			'gamma' => 'TA[TC]TACTGT(\w{5,})\w{34}(AAAAGCCAG|AAAGGCTTG|ACAAAGCTC)', # 'TA[TC]TACTGT(\w+?)TTTGC\w{4}GG[AG]AC[AT]A'
			'delta' => 'TA[TC][TC][AT]CTGT(G\w{4,})\w{33}AAAAGCCAG' }, # 'TA[TC][TC][AT]CTGT(G\w+?)TTTGGA[AC][CA]\wGG[AC]A'

	'vole' => {	'beta'  => 'TG[TC]([GA][CG]\w{3,})\w{31}AGGA[CT]CTGA' }, # 'T[AG][TC]\w{3}TG[TC]([GA][CG]\w+)(TTT|TTC)GG\w{4}GG'
	
	'tcrex' => {	'alpha' => 'TA[TC]\w[TA][CT]TG[TC](\w{5,})\w{34}ACATCCAGAA', # 'TA[TC]\w[TA][CT]TG[TC](\w+?)((TTT|TTC|TTG|TTA|TCT|GTT)GG|GTTGA|TTTGC|TTTAG)\w{4}(GG|GA|GT|TG|AG|AT)'
			# Original mouse beta C region sequence is AGGATCTGA, it has been changed to A[GA]GATCTGA to support the TR-x vector
			'beta'  => 'T[AT]\w[TC][TA][TCG][TG]G[TCG]([GAT][CG]\w{3,})\w{31}AAGATCTGA', # 'T[AT]\w[TC][TA][TCG][TG]G[TCG]([GAT][CG]\w+?)((TTT|TTC)GG|TCTTG|TTTGC|GATTG)\w{4}(GG|CC)'
			'gamma' => 'TA[TC]TACTGT(\w{5,})\w{34}(AAAAGCCAG|AAAGGCTTG|ACAAAGCTC)', # 'TA[TC]TACTGT(\w+?)TTTGC\w{4}GG[AG]AC[AT]A'
			'delta' => 'TA[TC][TC][AT]CTGT(G\w{4,})\w{33}AAAAGCCAG' }, # 'TA[TC][TC][AT]CTGT(G\w+?)TTTGGA[AC][CA]\wGG[AC]A'

			# Library preparatin PCR finishes before the C region, specific patterns must be used:
	'tcrlib' => {	'alpha' => 'TA[TC]\w[TA][CT]TG[TC](\w{5,})\w{34}ACCGAAGAGCAAG',
			'beta'  => 'T[AT]\w[TC][TA][TCG][TG]G[TCG]([GAT][CG]\w{3,})\w{34}AACGAAGAGCAAG' },

};
# 	'human' => {	'alpha' => 'TA[TC]\w[TA][CT]TG[TC](\w+?)\w{34}ATATCCAGAA', # 'TA[TC]\w[TA][CT]TG[TC](\w+?)(TTT|TTC|TGG|TGC)G[GC]\w{4}G[GA][AT]'
# 			'beta'  => 'TA[TC]\w{3}[TC][GA][TC]([AG][GC]\w+?)\w{31}AGGACCTGA', # 'TA[TC]\w{3}[TC][GA][TC]([AG][GC]\w+?)(TTT|TTC|GTC|TGG)GG\w{4}GG'
# 			'gamma' => '[TA]A[TC][TC]ACTG[TC](\w+?)\w{34}AACAACTTGA', # '[TA]A[TC][TC]ACTG[TC](\w+?)TTTG[GC]\w{4}GG[AG]AC[AT]A'
# 			'delta' => 'TACT[AT][CT]TGT(GC\w+?)\w{33}AGAAGTCAG' }, # 'TACT[AT][CT]TGT(GC\w+?)TT[TC]GG[AC]A[AC]\wGG[AC]A'
# 
# 	# 'monkey' => {	'beta'  => 'T\w[TC]\wT\w[TC]T\wTG[TC]([GCA][CG]\w+?)(TTT|TTC|CTG)GG\w{4}GG',
# 
# 	'mouse' => {	'alpha' => 'TA[TC]\w[TA][CT]TG[TC](\w+?)\w{34}ACATCCAGAA', # 'TA[TC]\w[TA][CT]TG[TC](\w+?)((TTT|TTC|TTG|TTA|TCT|GTT)GG|GTTGA|TTTGC|TTTAG)\w{4}(GG|GA|GT|TG|AG|AT)'
# 			'beta'  => 'T[AT]\w[TC][TA][TCG][TG]G[TCG]([GAT][CG]\w+?)\w{31}AGGATCTGA', # 'T[AT]\w[TC][TA][TCG][TG]G[TCG]([GAT][CG]\w+?)((TTT|TTC)GG|TCTTG|TTTGC|GATTG)\w{4}(GG|CC)'
# 			'gamma' => 'TA[TC]TACTGT(\w+?)\w{34}(AAAAGCCAG|AAAGGCTTG|ACAAAGCTC)', # 'TA[TC]TACTGT(\w+?)TTTGC\w{4}GG[AG]AC[AT]A'
# 			'delta' => 'TA[TC][TC][AT]CTGT(G\w+?)\w{33}AAAAGCCAG' }, # 'TA[TC][TC][AT]CTGT(G\w+?)TTTGGA[AC][CA]\wGG[AC]A'
# 
# 	'vole' => {	'beta'  => 'TG[TC]([GA][CG]\w+)\w{31}AGGA[CT]CTGA' }, # 'T[AG][TC]\w{3}TG[TC]([GA][CG]\w+)(TTT|TTC)GG\w{4}GG'
# };

my $REFERENCE_FOLDER = dirname(__FILE__).'/tcrefs/';
my $REFERENCE_FILES = {};
foreach my $species (keys %$CDR3_PATTERNS){
	# TR-x needs human reference sequences
	if ($species eq 'tcrex' || $species eq 'tcrlib') {
		$REFERENCE_FILES->{$species}{'tcrv'} = $REFERENCE_FOLDER.'human.tcrv.fna';
		$REFERENCE_FILES->{$species}{'tcrj'} = $REFERENCE_FOLDER.'human.tcrj.fna';
	} else {
		$REFERENCE_FILES->{$species}{'tcrv'} = $REFERENCE_FOLDER.$species.'.tcrv.fna';
		$REFERENCE_FILES->{$species}{'tcrj'} = $REFERENCE_FOLDER.$species.'.tcrj.fna';
	}
}


my ($INP_amplicons_file, @INP_read_files, $INP_cdr3_path, $INP_outpath, $INP_species, $INP_pattern, $INP_nreads, $INP_nseqs, $INP_numis, $INP_readsense, $INP_umi_cluster, $INP_threads, $INP_noref, %INP_ref_files, $INP_nocdr3, $INP_chao, $INP_onlystats, $INP_zip, $INP_verbose);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s{,}' => \@INP_read_files,
	'd|data=s' => \$INP_amplicons_file,
	'o|output=s' => \$INP_outpath,
	's|species=s' => \$INP_species,
	'p|pattern=s' => \$INP_pattern,
# 	'f|format=s' => \$INP_output_format,
	'e|extra:i' => \$INP_extra_length,
	'nr|numreads=i' => \$INP_nreads,
	'ns|numseqs=i' => \$INP_nseqs,
	'nu|numumis=i' => \$INP_numis,
	'rs|readsense=s' => \$INP_readsense,
	'rc|revcomp' => \$INP_revcomp,
	'min|minlen=i' => \$INP_minlen,
	'max|maxlen=i' => \$INP_maxlen,
	'if|inframe' => \$INP_noframe,
	'si|single' => \$INP_singletons,
	'cl|cluster=i' => \$INP_cluster_errors,
	'umi:s' => \$INP_umi_cluster,
	'ucl|umicluster=i' => \$INP_umi_errors,
	'noref' => \$INP_noref,
	'vref=s' => \$INP_ref_files{'tcrv'},
	'jref=s' => \$INP_ref_files{'tcrj'},
	'nocdr3' => \$INP_nocdr3,
	'cdr3=s' => \$INP_cdr3_path,
	'chao2' => \$INP_chao,
	'os|onlystats' => \$INP_onlystats,
	'thr|threads=i' => \$INP_threads,
	'v|verbose' => \$INP_verbose,
	'z|zip' => \$INP_zip,
	'<>' => \&usage,
);

# Usage help
sub usage
{
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -i (<file1> <file2>|<path>)\n\t\tInput paired-end read files in FASTQ format (compressed or uncompressed) or path.\n";
	print "  -d <file>\tCSV file with primer/amplicon data.\n";
	print "  -o <path>\tOutput folder name.\n";
	print "  -s <species>\tSpecies to analyze ('".join("','",keys %$CDR3_PATTERNS)."').\n";
	print "  -p <pattern>\tCustom pattern in REGEX format to extract the CDR3 region from the TCR sequences.\n";
# 	print "  -f <format>\tOutput format (default=$INP_output_format).\n";
	print "  -e <number>\tNumber of additional nucleotides to extract from both sides of CDR3 region.\n";
	print "  -nr <number>\tMaximum number of total reads/sequences to analyze.\n";
	print "  -ns <number>\tMaximum number of CDR3 sequences per sample to analyze.\n";
	print "  -nu <number>\tMaximum number of Unique Molecular Identifiers per sample to analyze.\n";
	print "  -rs [auto|R1+|R1-|R2+|R2-]\tIndicates the sense of the read containing the CDR3 (default=auto).\n";
	print "  -rc\t\tChecks the reverse complementary sequences if it fails to recognize the CDR3 region in the original orientation.\n";
	print "  -min <len>\tMinimum CDR3 sequence length (default=$INP_minlen)\n";
	print "  -max <len>\tMaximum CDR3 sequence length (default=$INP_maxlen)\n";
	print "  -if\t\tFilter CDR3 sequences out of frame or with stop codons(default=".$yes_no{$INP_noframe}.")\n";
	print "  -si\t\tKeep singletons\n";
	print "  -cl <number>\tCluster similar CDR3 sequences with the same or lower number of substitutions (default=$INP_cluster_errors)\n";
	print "  -umi\t\tAnnotates Unique Molecular Identifiers.\n";
	print "  -ucl <number>\tCluster similar CDR3 sequences within the same UMI with the same or lower number of substitutions (default=$INP_umi_errors)\n";
	print "  -noref\tSkips V-J gene allele annotation.\n";
	print "  -vref <file>\tFASTA file with TRBV reference sequences.\n";
	print "  -jref <file>\tFASTA file with TRBJ reference sequences.\n";
	print "  -nocdr3\tExtracts only reads that doesn't contain CDR3 sequences.\n";
	print "  -chao2\tPrints results from Chao2 estimator (requires multiple replicates).\n";
	print "  -os\t\tPrints only analysis statistcs, no sequence files.\n";
	print "  -cdr3 (<file>|<path>)\tInput FASTA file or path with already extracted CDR3 sequences.\n";
	#   print "  -debug Prints additional info for debugging.\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
	print "  -v\t\tVerbose mode, generates extra TXT files with additional information about the analysis.\n";
	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}

# Uses as input a FASTA file with already extracted CDR3 sequences
# my $INP_read_previous_results;
if (defined($INP_cdr3_path)){
	$INP_read_files[0] = $INP_cdr3_path;
	# $INP_read_previous_results = 1;
}
# Reads folder with previously extracted CDR3s
my ($INP_filepath, $INP_multifile);
if (-d $INP_read_files[0]) {
	$INP_filepath = 1;
# 	if (!defined($INP_outpath)){
# 		$INP_outpath = $INP_read_files[0];
# 	}
# Checks if a set of demultiplexed files is given as input into a compressed file
} elsif (defined($INP_read_files[0]) && is_multifile($INP_read_files[0])){
	$INP_multifile = $INP_read_files[0];
# Prints usage help if no input file is specified
} elsif (!defined($INP_read_files[0]) || !-f $INP_read_files[0]){
	print "\nERROR: You must specify a sequence input file.\n\n";
	usage();
	exit;
}
# if (!defined($INP_amplicons_file) || !-f $INP_amplicons_file){
# 	print "\nERROR: You must specify an amplicon data input file.\n\n";
# 	usage();
# 	exit;
# }
if (!defined($INP_ref_files{'tcrv'})){
	delete($INP_ref_files{'tcrv'});
} elsif ( !-f $INP_ref_files{'tcrv'} ){
	print "\nERROR: You must specify a valid FASTA file with TRV reference sequences.\n\n";
	usage();
	exit;
}
if (!defined($INP_ref_files{'tcrj'})){
	delete($INP_ref_files{'tcrj'});
} elsif ( !-f $INP_ref_files{'tcrj'} ){
	print "\nERROR: You must specify a valid FASTA file with TRJ reference sequences.\n\n";
	usage();
	exit;
}
if (defined($INP_species) && !defined($CDR3_PATTERNS->{$INP_species})){
	print "\nERROR: You must specify a valid species name.\n\n";
	usage();
	exit;
} elsif (defined($INP_species) && defined($INP_pattern)){
	print "\nWARNING: The internal pre-defined species CDR3 patterns will have priority respect to the custom one.\n\n";
# 	print "\nERROR: You cannot use a custom CDR3 pattern (-p) and at the same time choosing a particular species (-s). Species have CDR3 patterns defined by default.\n\n";
# 	exit;
}
if (defined($INP_readsense) && $INP_readsense !~ /auto|R[12][+-]/){
	print "\nERROR: You must specify a valid CDR3 read sense: auto, R1+, R1-, R2+ or R2-.\n\n";
	usage();
	exit;
}


print "\nRunning '$COMMAND_LINE'\n";

# Checks and reads amplicons information
my ($markerdata,$markers,$sampledata,$samples,$primer_seqs,$primer_headers,$paramsdata,$alleledata);
if (defined($INP_amplicons_file)){
	($markerdata,$markers,$sampledata,$samples,$primer_seqs,$primer_headers,$paramsdata,$alleledata)
	= parse_amplicon_file($INP_amplicons_file,['skip errors']);
}

# If there is no marker/primer data, then define empty markers/TCR chains
# Important for 'extract_cdr3_seqs' to work properly
if (!defined($markerdata) || !%$markerdata){
	# print "\nWARNING: No primer/amplicon data provided, default TCR alpha, beta, gamma and delta patterns will be used.\n\n";
	$markerdata = { 'TCRA' => { 'primer_f_id'=>[''], 'primer_r_id'=>[''], 'primer_f' => [''], 'primer_r' => [''] },
			'TCRB' => { 'primer_f_id'=>[''], 'primer_r_id'=>[''], 'primer_f' => [''], 'primer_r' => [''] },
			'TCRG' => { 'primer_f_id'=>[''], 'primer_r_id'=>[''], 'primer_f' => [''], 'primer_r' => [''] },
			'TCRD' => { 'primer_f_id'=>[''], 'primer_r_id'=>[''], 'primer_f' => [''], 'primer_r' => [''] }
			};
	$markers = ['TCRA', 'TCRB', 'TCRG', 'TCRD' ]
}
# If there is no sample data, then define a unique sample with the name of the file
my $sample_to_files;
if (!defined($sampledata) || !%$sampledata){
	print "\nWARNING: No sample tags/barcodes provided, the file/s will be treated as unique samples.\n\n";
	my $samplename = $INP_read_files[0];
	if ($INP_read_files[0] =~ /(.+?)(_R[12]_.+?)?\.(fa|fq|fasta|fastq)?\.?(gz|gzip)?$/){
		$samplename = $1;
	}
	$samplename = basename($samplename);
	$sampledata = { $samplename => {} };
	$samples = [ $samplename ];

	# Incorporates UMI sequence is neccesary
	if (defined($INP_umi_cluster) && $INP_umi_cluster =~ /[ACTG]+/){
		# Incorporates the parenthesis if it is ommited:
		if ($INP_umi_cluster !~ /\(.+\)/){
			$INP_umi_cluster =~ s/N(.+)N/(N$1N)/;
		}
		$sampledata->{$samplename}{'tag_f'} = $INP_umi_cluster;
	}
	$sample_to_files->{$samplename} = [ @INP_read_files ] ;
}

# Creates sample information, each file will be a sample
# But primer and sample data will be taken from .csv file
my $tmp_dir;
# if (!defined($INP_cdr3_path) && (defined($INP_multifile) || defined($INP_filepath))){
if (!defined($INP_cdr3_path) && (!defined($sampledata) || defined($INP_multifile) || defined($INP_filepath)) ){

	# TCR marker/primer data is defined in the .csv file
	my $marker_name = $markers->[0];
	my $sampledata_ = {};
	if (@$samples && defined($sampledata->{$samples->[0]})){
		$sampledata_ = $sampledata->{$samples->[0]};
	}
	# Incorporates UMI sequence is neccesary
	if (defined($INP_umi_cluster) && $INP_umi_cluster =~ /[ACTG]+/){
		# Incorporates the parenthesis if it is ommited:
		if ($INP_umi_cluster !~ /\(.+\)/){
			$INP_umi_cluster =~ s/N(.+)N/(N$1N)/;
		}
		$sampledata_->{'tag_f'} = $INP_umi_cluster;
	}
	# Removes previous sample name from .csv file
	undef($samples);
	undef($sampledata);

	# Uncompress the files in temporal folder
	my $paired_read_files = [];
	if (defined($INP_multifile)){
		$tmp_dir = "/tmp/".random_file_name();
		$paired_read_files= extract_paired_read_files_from_multifiles(\@INP_read_files,$tmp_dir);
		if (!@$paired_read_files){
			$paired_read_files = [ read_files_from_path($tmp_dir) ];
		}
	} elsif (defined($INP_filepath)){
		$paired_read_files = extract_paired_read_files_from_path($INP_read_files[0]);
		if (!@$paired_read_files){
			$paired_read_files = [ read_files_from_path($INP_read_files[0]) ];
		}
	} elsif (!defined($sampledata)){
		$paired_read_files = [ \@INP_read_files ];
	}
		
	# Each file will be an amplicon where the folder will give the marker name and the filename the sample name
	foreach my $paired_read_file_pair (@$paired_read_files) {
		my ($file1, $file2);
		if (ref($paired_read_file_pair) eq 'ARRAY') {
			$file1 = $paired_read_file_pair->[0];
			if (defined($paired_read_file_pair->[1])){
				$file2 = $paired_read_file_pair->[1];
			}
		} else {
			$file1 = $paired_read_file_pair;
		}
# 		unless (-f $file1 && (is_fasta($file1) || is_fastq($file1))){
# 			next;
# 		}
		# Gives the name of the internal file to the sample (without extension)
		my $sample_name;
		if ($file1 =~ /.+\/(.+?)\./){
			$sample_name = $1;
		# If not gives the full name of the internal file to the sample
		} else {
			$sample_name = $file1;
		}
		push(@$samples,$sample_name);
		# Copies the sample data from the .csv file
		$sampledata->{$sample_name} = $sampledata_;
		# Stores the path of the file/s
		if (!defined($file2)){
			$sample_to_files->{$sample_name} = [$file1];
		} else {
			$sample_to_files->{$sample_name} = [$file1, $file2];
		}
	}
	# Sorts samples by name
	$samples = [ nsort(@$samples) ];
}

# Creates variables to store statistics
my ($stats_cdr3_global, $stats_cdr3_common, $stats_cdr3_lengths, $stats_cdr3_depths, $stats_tcr_regions);

# Variables to store the CDR3 sequences and names of output files
my $output_files;
my $excel_outputfile;
my $cdr3_headers = {};
my $cdr3_seqs = {};
my $cdr3_qualities = {};
my $tcr_segments = {};
my $total_sample_reads = {};
my $total_tcr_reads = {};

# Configures CDR3 extraction options
my $options = [$INP_output_format];
if (defined($INP_umi_cluster)){
	push(@$options,'umi');
}
if (defined($INP_nocdr3)){
	push(@$options,'nocdr3');
}
if (defined($INP_readsense)){
	push(@$options,$INP_readsense);
} else {
	push(@$options,'auto');
}
if ($INP_revcomp==1){
	push(@$options,'revcomp');
}
# Defines a generic pattern to use when there is no an specific one available
my %INP_ref_patterns;
if (defined($INP_species)){
	$INP_species = lc($INP_species);
	foreach my $chain ('alpha','beta','gamma','delta'){
		if (!defined($INP_ref_patterns{$chain}) && defined($CDR3_PATTERNS->{$INP_species}{$chain})){
			$INP_ref_patterns{$chain} = $CDR3_PATTERNS->{$INP_species}{$chain};
		}
	}
	if (!defined($INP_noref) && !defined($INP_ref_files{'tcrj'}) && defined($REFERENCE_FILES->{$INP_species}{'tcrj'}) && -f $REFERENCE_FILES->{$INP_species}{'tcrj'}){
		$INP_ref_files{'tcrj'} = $REFERENCE_FILES->{$INP_species}{'tcrj'};
	}
	if (!defined($INP_noref) && !defined($INP_ref_files{'tcrv'}) && defined($REFERENCE_FILES->{$INP_species}{'tcrv'}) && -f $REFERENCE_FILES->{$INP_species}{'tcrv'}){
		$INP_ref_files{'tcrv'} = $REFERENCE_FILES->{$INP_species}{'tcrv'};
	}
}
if (defined($INP_pattern)){
	$INP_ref_patterns{'custom'} = $INP_pattern;
} elsif (!%INP_ref_patterns) {
	$INP_ref_patterns{'custom'} = $CDR3_GENERIC_PATTERN;
}

# Extracts CDR3 sequences if they must not be read from previous extraction files
# Counts total number of CDR3 sequences, calculates statistics and stores sequences
my $count_cdr3 = 0;
print "\n";
if (!defined($INP_cdr3_path)) {

	# Default filename for result files
	my ($output_folder,$output_name);
	if ($INP_read_files[0] =~ /(.+\/)?(.+?)\./){
		$output_folder = $1;
		$output_name = $2;
	} else {
		$output_folder = './';
		$output_name = 'amplicdr3';
	}

	if (!defined($INP_outpath)){
		$INP_outpath = $output_folder.$output_name;
	}

	# Creates folder to store results
	if (!-d $INP_outpath) {
		mkdir($INP_outpath);
	}
	# Variable to store the path of the Excel result file
	$excel_outputfile = "$INP_outpath/$output_name.stats.xlsx";

	# If the demultiplexed files are given as a compressed file or into a folder,
	# then the CDR3 sequences are extracted file by file (sample by sample)
	if (defined($INP_multifile) || defined($INP_filepath)){
		foreach my $sample (@$samples){
			my $read_files = $sample_to_files->{$sample};
			my ($read_file_format, $total_reads) = parse_sequence_file($read_files->[0],undef,['stats']);
			$total_sample_reads->{$sample} = $total_reads;
			my $cdr3_headers_ = {};
			my $cdr3_seqs_ = {};
			my $cdr3_qualities_ = {};
			my $tcr_segments_ = {};
			my $total_tcr_reads_ = {};
			if (defined($INP_nreads)){
				printf("Extracting '%s' TCR/CDR3 sequences from %d random reads (total=%d).\n",$sample,$INP_nreads,$total_reads);
				my $read_file_name = basename($read_files->[0]);
				if ($read_file_name =~ /(.+?)(\.(fa|fq|fasta|fastq)?)(\.?(gz|gzip)?)$/) {
					$read_file_name = "$INP_outpath/$1.$INP_nreads$2$4";
				}
				if ($INP_nreads < $total_reads) {
					if (scalar @$read_files == 2) {
						my $read_file_name2 = basename($read_files->[1]);
						if ($read_file_name2 =~ /(.+?)(\.(fa|fq|fasta|fastq)?)(\.?(gz|gzip)?)$/) {
							$read_file_name2 = "$INP_outpath/$1.$INP_nreads$2$4";
						}
						$read_files = shuffle_fastq_files($read_files,$INP_nreads,[$read_file_name,$read_file_name2],'gzip');
					} else {
						$read_files = [shuffle_fastq_file($read_files->[0],$INP_nreads,$read_file_name,'gzip')];
					}
					$total_sample_reads->{$sample} = $INP_nreads;
				} elsif ($INP_nreads > $total_reads) {
					printf("*WARNING: Sample '%s' has only %d reads.\n", $sample, $total_reads);
				}
			} else {
				printf("Extracting '%s' TCR/CDR3 sequences from %d total reads.\n",$sample,$total_reads);
			}
			if (!defined($INP_threads) || $INP_threads <= 1){
				($cdr3_headers_, $cdr3_seqs_, $cdr3_qualities_,$total_tcr_reads_,$tcr_segments_) = 
				extract_cdr3_seqs($read_files,$markerdata,{ $sample => $sampledata->{$sample} },$total_reads,\%INP_ref_patterns,$options); #,\%INP_ref_files,$options);
			} else {
				($cdr3_headers_, $cdr3_seqs_, $cdr3_qualities_,$total_tcr_reads_,$tcr_segments_) = 
				extract_cdr3_seqs_with_threads($read_files,$markerdata,{ $sample => $sampledata->{$sample} },$total_reads,\%INP_ref_patterns,$options,$INP_threads); #,\%INP_ref_files,$options,$INP_threads);
			}
			if (defined($cdr3_headers_->{$sample})){
				$cdr3_headers->{$sample} = $cdr3_headers_->{$sample};
				$cdr3_seqs->{$sample} = $cdr3_seqs_->{$sample};
			}
			if (defined($cdr3_qualities_->{$sample})){
				$cdr3_qualities->{$sample} = $cdr3_qualities_->{$sample};
			}
			if (defined($tcr_segments_->{$sample})){
				$tcr_segments->{$sample} = $tcr_segments_->{$sample};
			}
			if (defined($total_tcr_reads_->{$sample})){
				$total_tcr_reads->{$sample} = $total_tcr_reads_->{$sample};
			}
		}
		#print "\n";
		# Removes temporal dir
		if (defined($tmp_dir)){
			`rm -rf $tmp_dir`;
		}
	# If the DNA tags for each individual sample are specified in the amplicon data,
	# then the CDR3 sequences are extracted from the input reads file or the paired-end files
	} else {
		# Checks input file format and number of reads
		my ($read_file_format, $total_reads) = parse_sequence_file($INP_read_files[0],undef,['stats']);
		foreach my $sample (@$samples){
			$total_sample_reads->{$sample} = $total_reads;
		}
		# Shuffle reads if we are not analyzing all
		# Creates a temporal file with the desired number of reads to analyze`
		my $read_files = [@INP_read_files];
		if (defined($INP_nreads)){
			printf("Extracting TCR/CDR3 sequences from %d random reads (total=%d).\n",$INP_nreads,$total_reads);
			my $read_file_name = basename($read_files->[0]);
			if ($read_file_name =~ /(.+?)(\.(fa|fq|fasta|fastq)?)(\.?(gz|gzip)?)$/) {
				$read_file_name = "$INP_outpath/$1.$INP_nreads$2$4";
			}
			if ($INP_nreads < $total_reads) {
				if (scalar @$read_files == 2) {
					my $read_file_name2 = basename($read_files->[1]);
					if ($read_file_name2 =~ /(.+?)(\.(fa|fq|fasta|fastq)?)(\.?(gz|gzip)?)$/) {
						$read_file_name2 = "$INP_outpath/$1.$INP_nreads$2$4";
					}
					$read_files = shuffle_fastq_files($read_files,$INP_nreads,[$read_file_name,$read_file_name2],'gzip');
				} else {
					$read_files = [shuffle_fastq_file($read_files->[0],$INP_nreads,$read_file_name,'gzip')];
				}
				foreach my $sample (@$samples){
					$total_sample_reads->{$sample} = $INP_nreads;
				}
			} elsif ($INP_nreads > $total_reads) {
				printf("*WARNING: File '%s' has only %d reads.\n", $INP_read_files[0], $total_reads);
			}
		} else {
			print "Extracting TCR/CDR3 sequences from $total_reads total reads.\n";
		}
		if (!defined($INP_threads) || $INP_threads <= 1){
			($cdr3_headers, $cdr3_seqs, $cdr3_qualities,$total_tcr_reads,$tcr_segments) = 
			extract_cdr3_seqs($read_files,$markerdata,$sampledata,$total_reads,\%INP_ref_patterns,$options); #,\%INP_ref_files,$options);
		} else {

			($cdr3_headers, $cdr3_seqs, $cdr3_qualities,$total_tcr_reads,$tcr_segments) = 
			extract_cdr3_seqs_with_threads($read_files,$markerdata,$sampledata,$total_reads,\%INP_ref_patterns,$options,$INP_threads); #,\%INP_ref_patterns,\%INP_ref_files,$options,$INP_threads);
		}
	}
	# 	# Prints number of reads parsed by sample
	# 	foreach my $sample (@$samples){
	# 		printf("Extracted '%s' CDR3 sequences from %d reads.\n",$sample,$total_sample_reads->{$sample});
	# 	}
	print "\n";

	# Counts total number of CDR3 sequences, calculates statistics and stores sequences
	foreach my $sample (@$samples){
		# Prints statistics
		if (defined($total_tcr_reads->{$sample})) {
			printf("Sample '%s': %8d TCR sequences found.\n", $sample, $total_tcr_reads->{$sample});
		}
		my $sample_count_cdr3s = 0;
		foreach my $marker (@$markers){
			my $filename;
			# Writes sequences into a file
			if (defined($cdr3_headers->{$sample}{$marker}) && @{$cdr3_headers->{$sample}{$marker}}){
				# Writes sequences ordered by UMIs
				if (defined($INP_umi_cluster)){
					my $umi_cdr3_seqs;
					my (@cdr3_seqs_,@cdr3_headers_);
					for (my $i=0; $i<=$#{$cdr3_headers->{$sample}{$marker}}; $i++){
						if ($cdr3_headers->{$sample}{$marker}[$i] =~ /umi=([ACGT]+)/){
							push(@{$umi_cdr3_seqs->{$1}},$i);
						}
					}
					#foreach my $umi (sort { $#{$umi_cdr3_seqs->{$b}} <=> $#{$umi_cdr3_seqs->{$a}} } keys %$umi_cdr3_seqs){
					foreach my $umi (keys %$umi_cdr3_seqs){
						foreach my $i (@{$umi_cdr3_seqs->{$umi}}){
							push(@cdr3_seqs_,$cdr3_seqs->{$sample}{$marker}[$i]);
							push(@cdr3_headers_,$cdr3_headers->{$sample}{$marker}[$i]);
						}
					}
					$cdr3_seqs->{$sample}{$marker} = \@cdr3_seqs_;
					$cdr3_headers->{$sample}{$marker} = \@cdr3_headers_;
				}
				if (!defined($INP_onlystats)){
					my $file_prefix = "$INP_outpath/$sample.$marker.total.cdr3";
					if (defined($INP_nocdr3)){
						$file_prefix = "$INP_outpath/$sample.$marker.nocdr3";
					}
					if ($INP_output_format eq 'fastq'){
						$filename = create_fastq_file($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}, $cdr3_qualities->{$sample}{$marker}, "$file_prefix.fq", 'gzip');
					} else {
						$filename = create_fasta_file($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}, "$file_prefix.fa", 'gzip');
					}
				}
				# Counts total CDR3 sequences and calculate statistics
				$output_files->{'total'}{$sample}{$marker} = $filename;
				$sample_count_cdr3s += scalar @{$cdr3_headers->{$sample}{$marker}};
				$count_cdr3 += scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_global->{$sample}{$marker}{'total'} = scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_lengths->{$sample}{$marker}{'total'} = count_lengths($cdr3_seqs->{$sample}{$marker});
				# $stats_tcr_regions->{$sample}{$marker}{'total'} = count_tcr_regions($cdr3_headers->{$sample}{$marker},$cdr3_seqs->{$sample}{$marker});
			} else {
	# 			$cdr3_headers->{$sample}{$marker} = [];
	# 			$cdr3_seqs->{$sample}{$marker} = [];
	# 			$cdr3_qualities->{$sample}{$marker} = [];
				$stats_cdr3_global->{$sample}{$marker}{'total'} = 0;
	# 			$stats_cdr3_lengths->{$sample}{$marker}{'total'} = {};
	# 			$stats_tcr_regions->{$sample}{$marker}{'total'} = {};
			}
		}
		# Prints statistics
		my @marker_count_cdr3s = map sprintf("%s: %d",$_, $stats_cdr3_global->{$sample}{$_}{'total'}) , @$markers;
		if ($sample_count_cdr3s && !defined($INP_onlystats)){
			if (!defined($INP_nocdr3)){
				printf("Sample '%s': %8d CDR3 sequences (%s) written into '.total.cdr3' files.\n", $sample, $sample_count_cdr3s, join(", ", @marker_count_cdr3s)); # , join(", ", values %{$output_files->{'total'}{$sample}})
			} else {
				printf("Sample '%s': %8d reads not matching CDR3 sequences (%s) written into '.nocdr3' files.\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s)); # , join(", ", values %{$output_files->{'total'}{$sample}})
			}
		} elsif ($sample_count_cdr3s) {
			if (!defined($INP_nocdr3)){
				printf("Sample '%s': %8d CDR3 sequences found (%s).\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s));
			} else {
				printf("Sample '%s': %8d reads not matching CDR3 sequences found (%s).\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s));
			}
		}
	}

} elsif (defined($INP_cdr3_path)) {

	# Default filename for result files
	my ($output_folder,$output_name);
	if ($INP_read_files[0] =~ /(.+\/)?(.+?)\./){
		$output_folder = $1;
		$output_name = $2;
	} else {
		$output_folder = './';
		$output_name = 'amplicdr3';
	}

	if (!defined($INP_outpath)){
		$INP_outpath = $output_folder.$output_name;
	}

	# Creates folder to store results
	if (!-d $INP_outpath) {
		mkdir($INP_outpath);
	}
	# Variable to store the path of the Excel result file
	$excel_outputfile = "$INP_outpath/$output_name.stats.xlsx";

	printf("\nReading previously extracted CDR3 sequences from '$INP_cdr3_path'.\n\n",);

	if (-d $INP_cdr3_path) {
		# Reads sequences from previous existing FASTA/FASTQ GZIPPED file
		foreach my $sample (@$samples){
			foreach my $marker (@$markers){
				my $filename;
				my $file_prefix = $INP_read_files[0]."/$sample.$marker.total.cdr3";
				if (-e "$file_prefix.fq.gz"){
					$filename = "$file_prefix.fq.gz";
				} elsif (-e "$file_prefix.fa.gz"){
					$filename = "$file_prefix.fa.gz";
				}
				if (defined($filename)){
					($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}) = read_sequence_file($filename);
					$tcr_segments->{$sample}{$marker} = read_tcr_segments($cdr3_headers->{$sample}{$marker},$cdr3_seqs->{$sample}{$marker});
					printf("Sample '%s': %8d %s CDR3 sequences read from '%s'.\n", $sample, scalar @{$cdr3_headers->{$sample}{$marker}}, $marker, $filename);
				} else {
					printf("Sample '%s': %8s %s CDR3 sequence file found.\n", $sample, 'no', $marker);
					next;
				}
				if (!defined($cdr3_headers->{$sample}{$marker})){
					if (!defined($INP_nocdr3)){
						printf("Sample '%s': %8d %s CDR3 sequences read from '%s'.\n", $sample, 0, $marker, $filename);
					} else {
						printf("Sample '%s': %8d reads not matching %s CDR3 sequences.\n", $sample, 0, $marker);
					}
					next;
				}
				# Counts total CDR3 sequences and calculate statistics
				$output_files->{'total'}{$sample}{$marker} = $filename;
				$count_cdr3 += scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_global->{$sample}{$marker}{'total'} = scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_lengths->{$sample}{$marker}{'total'} = count_lengths($cdr3_seqs->{$sample}{$marker});
				# $stats_tcr_regions->{$sample}{$marker}{'total'} = count_tcr_regions($cdr3_headers->{$sample}{$marker},$cdr3_seqs->{$sample}{$marker});
			}
		}
# 	} elsif (-f $INP_cdr3_path) {
# 		foreach my $sample (@$samples){
# 			my $read_file; 
# 			if (defined($INP_multifile)){
# 				$read_file = $sample_to_files->{$sample};
# 			} else {
# 				$read_file = $INP_cdr3_path;
# 			}
# 			($cdr3_seqs->{$sample}, $cdr3_headers->{$sample}) = read_sequence_file($read_file);
# 			printf("Sample '%s': %8d CDR3 sequences read from '%s'.\n", $sample, scalar @{$cdr3_headers->{$sample}}, $read_file);
# 			if (!defined($cdr3_headers->{$sample})){
# 				if (!defined($INP_nocdr3)){
# 					printf("Marker '%s', sample '%s': %8d CDR3 sequences read from '%s'.\n", $marker, $sample, 0, $read_file);
# 				} else {
# 					printf("Marker '%s', sample '%s': %8d reads not matching CDR3 sequences.\n", $marker, $sample, 0);
# 				}
# 				next;
# 			}
# 			# Counts total CDR3 sequences and calculate statistics
# 			$output_files->{'total'}{$sample} = $read_file;
# 			$count_cdr3 += scalar @{$cdr3_headers->{$sample}};
# 			$stats_cdr3_global->{$sample}{'total'} = scalar @{$cdr3_headers->{$sample}};
# 			$stats_cdr3_lengths->{$sample}{'total'} = count_lengths($cdr3_seqs->{$sample});
# 			$stats_tcr_regions->{$sample}{'total'} = count_tcr_regions($cdr3_headers->{$sample},$cdr3_seqs->{$sample});
# 		}
# 		# Removes temporal dir
# 		if (defined($INP_multifile)){
# 			`rm -rf $tmp_dir`;
# 		}
	} else {
		printf("\nCould not read previously extracted CDR3 sequences from '$INP_cdr3_path'.\n\n",);
		exit;
	}
}

if (defined($INP_umi_cluster)) {
	foreach my $sample (@$samples){
		my @marker_count_umis;
		my $total_sample_umis = 0;
		foreach my $marker (@$markers){
			my %count_umis = count_umis($cdr3_headers->{$sample}{$marker});
			my $count_umis = scalar keys %count_umis;
			$total_sample_umis += $count_umis;
			push(@marker_count_umis, sprintf("%s: %d",$marker,$count_umis));
		}
		printf("Sample '%s': %8d Unique Molecular Identifiers found (%s).\n", $sample, $total_sample_umis, join(", ",@marker_count_umis));
	}
}

if (!defined($INP_nocdr3)){
	printf("\n%d CDR3 total sequences found.\n", $count_cdr3);
} else {
	printf("\n%d reads not matching CDR3 sequences found.\n", $count_cdr3);
	exit;
}


# #DEBUGGING:
# # Stops reading previous results
# undef($INP_read_previous_results);

# Selects only the desired number of CDR3 sequences, calculates statistics and stores sequences
if ($INP_nseqs){
	printf("\nSelecting %d CDR3 sequences per sample.\n\n",$INP_nseqs);
	$count_cdr3 = 0;
	# Takes the desired number of random CDR3 sequences and headers
	foreach my $sample (@$samples){
		my $sample_count_cdr3s = 0;
		foreach my $marker (@$markers){
			my $filename;
			my @random = shuffle 0..$#{$cdr3_headers->{$sample}{$marker}};
			if ($INP_nseqs < scalar @random){
				splice(@random, $INP_nseqs);
			} elsif ($INP_nseqs > scalar @random) {
				printf("*WARNING: Sample '%s' has only %d CDR3 sequences for %s.\n", $sample, scalar @{$cdr3_headers->{$sample}{$marker}}, $marker);
			}
			$cdr3_headers->{$sample}{$marker} = [ map $cdr3_headers->{$sample}{$marker}[$_], @random ];
			$cdr3_seqs->{$sample}{$marker} = [ map $cdr3_seqs->{$sample}{$marker}[$_], @random ];
			if (defined($cdr3_qualities) && defined($cdr3_qualities->{$sample}{$marker})){
				$cdr3_qualities->{$sample}{$marker} = [ map $cdr3_qualities->{$sample}{$marker}[$_], @random ];
			}
			if (defined($cdr3_headers->{$sample}{$marker}) && @{$cdr3_headers->{$sample}{$marker}}){
				# Writes sequences ordered by UMIs
				if (defined($INP_umi_cluster)){
					my $umi_cdr3_seqs;
					my (@cdr3_seqs_,@cdr3_headers_);
					for (my $i=0; $i<=$#{$cdr3_headers->{$sample}{$marker}}; $i++){
						if ($cdr3_headers->{$sample}{$marker}[$i] =~ /umi=([ACGT]+)/){
							push(@{$umi_cdr3_seqs->{$1}},$i);
						}
					}
					#foreach my $umi (sort { $#{$umi_cdr3_seqs->{$b}} <=> $#{$umi_cdr3_seqs->{$a}} } keys %$umi_cdr3_seqs){
					foreach my $umi (keys %$umi_cdr3_seqs){
						foreach my $i (@{$umi_cdr3_seqs->{$umi}}){
							push(@cdr3_seqs_,$cdr3_seqs->{$sample}{$marker}[$i]);
							push(@cdr3_headers_,$cdr3_headers->{$sample}{$marker}[$i]);
						}
					}
					$cdr3_seqs->{$sample}{$marker} = \@cdr3_seqs_;
					$cdr3_headers->{$sample}{$marker} = \@cdr3_headers_;
				}
				if (!defined($INP_onlystats)){
					my $file_prefix = "$INP_outpath/$sample.$marker.$INP_nseqs.cdr3";
					if ($INP_output_format eq 'fastq'){
						$filename = create_fastq_file($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}, $cdr3_qualities->{$sample}{$marker}, "$file_prefix.fq", 'gzip');
					} else {
						$filename = create_fasta_file($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}, "$file_prefix.fa", 'gzip');
					}
				}
				# Counts selected CDR3 sequences and calculate statistics
				$output_files->{'selected'}{$sample}{$marker} = $filename;
				$count_cdr3 += scalar @{$cdr3_headers->{$sample}{$marker}};
				$sample_count_cdr3s += scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_global->{$sample}{$marker}{'selected'} = scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_lengths->{$sample}{$marker}{'selected'} = count_lengths($cdr3_seqs->{$sample}{$marker});
				# $stats_tcr_regions->{$sample}{$marker}{'selected'} = count_tcr_regions($cdr3_headers->{$sample}{$marker},$cdr3_seqs->{$sample}{$marker});
			} else {
				$stats_cdr3_global->{$sample}{$marker}{'selected'} = 0;
			}
		}
		# Prints statistics
		my @marker_count_cdr3s = map sprintf("%s: %d",$_, $stats_cdr3_global->{$sample}{$_}{'selected'}) , @$markers;
		if (!defined($INP_onlystats)){
			printf("Sample '%s': %8d CDR3 sequences (%s) written into '%s'.\n", $sample, $sample_count_cdr3s, join(", ", @marker_count_cdr3s), join(", ", values %{$output_files->{'selected'}{$sample}}));
		} else {
			printf("Sample '%s': %8d CDR3 sequences (%s).\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s));
		}
		if (defined($INP_umi_cluster)) {
			foreach my $sample (@$samples){
				my @marker_count_umis;
				my $total_sample_umis = 0;
				foreach my $marker (@$markers){
					my %count_umis = count_umis($cdr3_headers->{$sample}{$marker});
					my $count_umis = scalar keys %count_umis;
					$total_sample_umis += $count_umis;
					push(@marker_count_umis, sprintf("%s: %d",$marker,$count_umis));
				}
				printf("Sample '%s': %8d Unique Molecular Identifiers found (%s).\n", $sample, $total_sample_umis, join(", ",@marker_count_umis));
			}
		}
	}
	printf("\n%d CDR3 random sequences selected.\n", $count_cdr3);

}

# if (defined($INP_umi_cluster)){
# 	exit;
# }

# DEBUGGING:
# Stops reading previous results
# undef($INP_cdr3_path);

# Filters CDR3 sequences (and group them by lengths)
# Counts filtered CDR3 sequences and calculates statistics
my $cdr3_filtered;
if (!defined($INP_cdr3_path)) {
	printf("\nFiltering %d CDR3 sequences (min_len: %d, max_len: %d, filter off-frame and stop-codon: %s).\n\n",$count_cdr3, $INP_minlen,$INP_maxlen,$yes_no{$INP_noframe});
	$count_cdr3 = 0;
	foreach my $sample (@$samples){
		my $sample_count_cdr3s = 0;
		my %sample_total_discarded;
		foreach my $marker (@$markers){
			$cdr3_filtered->{$sample}{$marker}{'minlen'} = 0;
			$cdr3_filtered->{$sample}{$marker}{'maxlen'} = 0;
			$cdr3_filtered->{$sample}{$marker}{'noframe'} = 0;
			$cdr3_filtered->{$sample}{$marker}{'stopcodon'} = 0;
			for (my $i=0; $i<=$#{$cdr3_headers->{$sample}{$marker}}; $i++){
				my $len = length($cdr3_seqs->{$sample}{$marker}[$i]);
				# Checks CDR3 length and if it is divisible by 3 (in-frame)
				my ($minlen, $maxlen, $noframe, $stopcodon) = (0,0,0,0);
				if ($len < $INP_minlen){
					$minlen = 1;
					$cdr3_filtered->{$sample}{$marker}{'minlen'}++;
				}
				if ($len > $INP_maxlen){
					$maxlen = 1;
					$cdr3_filtered->{$sample}{$marker}{'maxlen'}++;
				}
				if ($INP_noframe) {
					if ($len % 3 != 0) {
						$noframe = 1;
						$cdr3_filtered->{$sample}{$marker}{'noframe'}++;
					} elsif (dna_to_prot($cdr3_seqs->{$sample}{$marker}[$i]) =~ /\*/) {
						$stopcodon = 1;
						$cdr3_filtered->{$sample}{$marker}{'stopcodon'}++;
					}
				}
				if ($minlen || $maxlen || $noframe || $stopcodon){
					$sample_total_discarded{$marker}++;
					splice(@{$cdr3_headers->{$sample}{$marker}},$i,1);
					splice(@{$cdr3_seqs->{$sample}{$marker}},$i,1);
					if (defined($cdr3_qualities) && defined($cdr3_qualities->{$sample}{$marker})){
						splice(@{$cdr3_qualities->{$sample}{$marker}},$i,1);
					}
					$i--;
				}
			}
			# Writes sequences into file
			my $filename;
			if (defined($cdr3_headers->{$sample}{$marker}) && @{$cdr3_headers->{$sample}{$marker}}){
				if (!defined($INP_onlystats)){
					my $file_prefix = "$INP_outpath/$sample.$marker.filtered.cdr3";
					if ($INP_output_format eq 'fastq'){
						$filename = create_fastq_file($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}, $cdr3_qualities->{$sample}{$marker}, "$file_prefix.fq", 'gzip');
					} else {
						$filename = create_fasta_file($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}, "$file_prefix.fa", 'gzip');
					}
				}
				# Counts filtered CDR3 sequences and calculate statistics
				$output_files->{'filtered'}{$sample}{$marker} = $filename;
				$count_cdr3 += scalar @{$cdr3_headers->{$sample}{$marker}};
				$sample_count_cdr3s += scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_global->{$sample}{$marker}{'filtered'} = scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_lengths->{$sample}{$marker}{'filtered'} = count_lengths($cdr3_seqs->{$sample}{$marker});
				# $stats_tcr_regions->{$sample}{$marker}{'filtered'} = count_tcr_regions($cdr3_headers->{$sample}{$marker},$cdr3_seqs->{$sample}{$marker});
			} else {
				$stats_cdr3_global->{$sample}{$marker}{'filtered'} = 0;
			}
		}
		# Prints statistics
		my @marker_count_cdr3s = map sprintf("%s: %d",$_, $stats_cdr3_global->{$sample}{$_}{'filtered'}) , @$markers;
		if (!defined($INP_onlystats)){
			printf("Sample '%s': %8d filtered valid CDR3 sequences (%s) written into '.filtered.cdr3' files.\n", $sample, $sample_count_cdr3s, join(", ", @marker_count_cdr3s)); # , join(", ", values %{$output_files->{'filtered'}{$sample}}));
		} else {
			printf("Sample '%s': %8d filtered valid CDR3 sequences (%s).\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s));
		}
		# Prints statistics about discarded sequences
# 		my @marker_discarded = map sprintf("%8d discarded %s CDR3 sequences: %d (min_len: %d, max_len: %d, off-frame: %d, stop-codon: %d)", $sample_total_discarded{$marker}, $_, $cdr3_filtered->{$sample}{$_}{'minlen'}+$cdr3_filtered->{$sample}{$_}{'maxlen'}+$cdr3_filtered->{$sample}{$_}{'noframe'}+$cdr3_filtered->{$sample}{$_}{'stopcodon'}, $cdr3_filtered->{$sample}{$_}{'minlen'}, $cdr3_filtered->{$sample}{$_}{'maxlen'}, $cdr3_filtered->{$sample}{$_}{'noframe'}, $cdr3_filtered->{$sample}{$_}{'stopcodon'}) , @$markers;
# 		printf("%s %s.\n", " "x(10+length($sample)), join("\n", @marker_discarded));
		# printf("%s %8d discarded CDR3 sequences: ", " "x(10+length($sample)), $sample_total_discarded);
		foreach my $marker (@$markers){
			if (defined($sample_total_discarded{$marker})){
				my $minlen = $cdr3_filtered->{$sample}{$marker}{'minlen'};
				my $maxlen = $cdr3_filtered->{$sample}{$marker}{'maxlen'};
				my $noframe = $cdr3_filtered->{$sample}{$marker}{'noframe'};
				my $stopcodon = $cdr3_filtered->{$sample}{$marker}{'stopcodon'};
				printf("%s %8d discarded %s CDR3 sequences (min_len: %d, max_len: %d, off-frame: %d, stop-codon: %d)\n", " "x(10+length($sample)), $sample_total_discarded{$marker}, $marker, $minlen, $maxlen, $noframe, $stopcodon);
			}
		}
	}
} else {
	printf("\nReading previously filtered CDR3 sequences.\n\n");
	$count_cdr3 = 0;
	# Reads sequences from previous existing FASTA/FASTQ GZIPPED file
	foreach my $sample (@$samples){
		foreach my $marker (@$markers){
			my $filename;
			my $file_prefix = $INP_read_files[0]."/$sample.$marker.filtered.cdr3";
			if (-e "$file_prefix.fq.gz"){
				$filename = "$file_prefix.fq.gz";
			} elsif (-e "$file_prefix.fa.gz"){
				$filename = "$file_prefix.fa.gz";
			}
			if (defined($filename)){
				($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}) = read_sequence_file($filename);
				printf("Sample '%s': %8d filtered %s CDR3 sequences read from '%s'.\n", $sample, scalar @{$cdr3_headers->{$sample}{$marker}}, $marker, $filename);
			} else {
				printf("Sample '%s': %8s filtered %s CDR3 sequence file found.\n", $sample, 'no', $marker);
				next;
			}
			if (!defined($cdr3_headers->{$sample}{$marker})){
				printf("Sample '%s': %8d filtered %s CDR3 sequences read from '%s'.\n", $sample, 0, $marker, $filename);
			}
			# Counts filtered CDR3 sequences and calculate statistics
			$output_files->{'filtered'}{$sample}{$marker} = $filename;
			$count_cdr3 += scalar @{$cdr3_headers->{$sample}{$marker}};
			$stats_cdr3_global->{$sample}{$marker}{'filtered'} = scalar @{$cdr3_headers->{$sample}{$marker}};
			$stats_cdr3_lengths->{$sample}{$marker}{'filtered'} = count_lengths($cdr3_seqs->{$sample}{$marker});
			# $stats_tcr_regions->{$sample}{$marker}{'filtered'} = count_tcr_regions($cdr3_headers->{$sample}{$marker},$cdr3_seqs->{$sample}{$marker});
		}
	}
}
if (defined($INP_umi_cluster)) {
	foreach my $sample (@$samples){
		my @marker_count_umis;
		my $total_sample_umis = 0;
		foreach my $marker (@$markers){
			my %count_umis = count_umis($cdr3_headers->{$sample}{$marker});
			my $count_umis = scalar keys %count_umis;
			$total_sample_umis += $count_umis;
			push(@marker_count_umis, sprintf("%s: %d",$marker,$count_umis));
		}
		printf("Sample '%s': %8d Unique Molecular Identifiers found (%s).\n", $sample, $total_sample_umis, join(", ",@marker_count_umis));
	}
}

printf("\n%d CDR3 sequences after filtering.\n", $count_cdr3);

# DEBUGGING:
# Stops reading previous results
# undef($INP_cdr3_path);

# # Selects only the desired number of UMI sequences, calculates statistics and stores sequences
# if (defined($INP_numis)){
# 	printf("\nSelecting %d UMI sequences per sample.\n\n",$INP_numis);
# 	$count_cdr3 = 0;
# 	print "\n";
# 	# Takes the desired number of UMIs sequences and headers
# 	foreach my $sample (@$samples){
# 		my $filename;
# 		my @random = shuffle 0..$#{$cdr3_headers->{$sample}};
# 		$cdr3_headers->{$sample} = [ map $cdr3_headers->{$sample}[$_], @random ];
# 		$cdr3_seqs->{$sample} = [ map $cdr3_seqs->{$sample}[$_], @random ];
# 		if (defined($cdr3_qualities) && defined($cdr3_qualities->{$sample})){
# 			$cdr3_qualities->{$sample} = [ map $cdr3_qualities->{$sample}[$_], @random ];
# 		}
# 		if (defined($cdr3_headers->{$sample}{$marker}) && @{$cdr3_headers->{$sample}{$marker}}){
# 			# Writes sequences ordered by UMIs
# 			if (defined($INP_umi_cluster)){
# 				my $umi_cdr3_seqs;
# 				my (@cdr3_seqs_,@cdr3_headers_);
# 				for (my $i=0; $i<=$#{$cdr3_headers->{$sample}}; $i++){
# 					if ($cdr3_headers->{$sample}[$i] =~ /umi=([ACGT]+)/){
# 						push(@{$umi_cdr3_seqs->{$1}},$i);
# 					}
# 					if ($INP_numis == scalar keys %$umi_cdr3_seqs) {
# 						last;
# 					}
# 				}
# 				#foreach my $umi (sort { $#{$umi_cdr3_seqs->{$b}} <=> $#{$umi_cdr3_seqs->{$a}} } keys %$umi_cdr3_seqs){
# 				foreach my $umi (keys %$umi_cdr3_seqs){
# 					foreach my $i (@{$umi_cdr3_seqs->{$umi}}){
# 						push(@cdr3_seqs_,$cdr3_seqs->{$sample}[$i]);
# 						push(@cdr3_headers_,$cdr3_headers->{$sample}[$i]);
# 					}
# 				}
# 				$cdr3_seqs->{$sample} = \@cdr3_seqs_;
# 				$cdr3_headers->{$sample} = \@cdr3_headers_;
# 			}
# 			my $file_prefix = "$INP_outpath/$sample.$INP_numis.umis";
# 			if ($INP_output_format eq 'fastq'){
# 				$filename = create_fastq_file($cdr3_seqs->{$sample}, $cdr3_headers->{$sample}, $cdr3_qualities->{$sample}, "$file_prefix.fq", 'gzip');
# 			} else {
# 				$filename = create_fasta_file($cdr3_seqs->{$sample}, $cdr3_headers->{$sample}, "$file_prefix.fa", 'gzip');
# 			}
# 			printf("Sample '%s': %8d CDR3 sequences written into '%s'.\n", $sample, scalar @{$cdr3_headers->{$sample}}, $filename);
# 		} else {
# 			printf("Sample '%s': %8d CDR3 sequences.\n", $sample, 0);
# 			next;
# 		}
# 		# Counts selected CDR3 sequences and calculate statistics
# 		$output_files->{'selected'}{$sample} = $filename;
# 		$count_cdr3 += scalar @{$cdr3_headers->{$sample}};
# 		$stats_cdr3_global->{$sample}{'selected'} = scalar @{$cdr3_headers->{$sample}};
# 		$stats_cdr3_lengths->{$sample}{'selected'} = count_lengths($cdr3_seqs->{$sample});
# 		$stats_tcr_regions->{$sample}{'selected'} = count_tcr_regions($cdr3_headers->{$sample},$cdr3_seqs->{$sample});
# 	}
# 	if (defined($INP_umi_cluster)) {
# 		foreach my $sample (@$samples){
# 			printf("Sample '%s': %8d Unique Molecular Identifiers selected.\n", $sample, scalar keys count_umis($cdr3_headers->{$sample}));
# 		}
# 	}
# }

# Calculates common CDR3 sequences between different samples from the same individual
# $stats_cdr3_common->{'filtered'} = common_cdr3_seqs($cdr3_seqs);

# Clusters CDR3 sequences
my %cdr3_prots;
my %count_singletons;
if (!defined($INP_cdr3_path)) {
	if (defined($INP_numis)){
		printf("\nClustering %d UMIs from %d CDR3 sequences (max. subs. within UMI: %d, max. subs.: %d, keep singletons: %s).\n",$INP_numis,$count_cdr3,$INP_umi_errors,$INP_cluster_errors,$yes_no{$INP_singletons});
	} elsif (defined($INP_umi_cluster)){
		printf("\nClustering UMIs from %d CDR3 sequences (max. subs. within UMI: %d, max. subs.: %d, keep singletons: %s).\n",$count_cdr3,$INP_umi_errors,$INP_cluster_errors,$yes_no{$INP_singletons});
	} else {
		printf("\nClustering %d CDR3 sequences (max. subs.: %d, keep singletons: %s).\n",$count_cdr3,$INP_cluster_errors,$yes_no{$INP_singletons});
	}
	print "\n";
	# my $start_run = gettimeofday;
	# Clustering options
	my %options = ( 'errors' => $INP_cluster_errors );
	if (defined($INP_umi_cluster)){
		$options{'umi_errors'} = $INP_umi_errors;
		if (defined($INP_numis)){
			$options{'umi_limit'} = $INP_numis;
		}
		if ($INP_singletons){
			$options{'keep_singletons'} = 1;
		}
	}

	my $seq_clusters;
	my $seq_stats;
	if (!defined($INP_threads) || $INP_threads <= 1){
		foreach my $sample (@$samples){
			foreach my $marker (@$markers){
				# Cluster sequences of the same length
				if (defined($cdr3_seqs->{$sample}{$marker})){
			# 		($cdr3_seqs->{$sample}{$marker},$cdr3_headers->{$sample}{$marker}) = cluster_cdr3s($cdr3_seqs->{$sample}{$marker},$cdr3_headers->{$sample}{$marker},$INP_cluster_errors,$INP_numis);
					($seq_clusters->{$sample}{$marker}, $seq_stats->{$sample}{$marker}) = cluster_cdr3s($cdr3_seqs->{$sample}{$marker},$cdr3_headers->{$sample}{$marker},\%options);
				}
			}
		}
	} else {
		($seq_clusters, $seq_stats) = cluster_sample_cdr3s_with_threads($cdr3_seqs,$cdr3_headers,\%options,$INP_threads);
	}
	# printf("Total time: %.2fs\n\n",gettimeofday-$start_run);
	if (defined($options{'umi_errors'})) {
		foreach my $sample (@$samples){
			my ($sample_correct_umis,$sample_incorrect_umis,$sample_singleton_umis) = (0,0,0);
			my @marker_count_umis;
			foreach my $marker (@$markers){
				my ($sample_correct_umis_,$sample_incorrect_umis_,$sample_singleton_umis_) = (0,0,0);
				if (defined($seq_stats->{$sample}{$marker}{'umis_correct'})){
					$sample_correct_umis_ = $seq_stats->{$sample}{$marker}{'umis_correct'};
				}
				if (defined($seq_stats->{$sample}{$marker}{'umis_incorrect'})){
					$sample_incorrect_umis_ = $seq_stats->{$sample}{$marker}{'umis_incorrect'};
				}
				if (defined($seq_stats->{$sample}{$marker}{'umis_singletons'})){
					$sample_singleton_umis_ = $seq_stats->{$sample}{$marker}{'umis_singletons'};
				}
				$sample_correct_umis += $sample_correct_umis_;
				$sample_incorrect_umis += $sample_incorrect_umis_;
				$sample_singleton_umis += $sample_singleton_umis_;
				push(@marker_count_umis, sprintf("%s: %d (incorrect: %d, singletons: %d)", $marker, $sample_correct_umis_, $sample_incorrect_umis_, $sample_singleton_umis_));
			}
			if (defined($options{'umi_limit'})){
				printf("Sample '%s': %8d correct UMIs selected, %d incorrect and %d singleton UMIs discarded - %s.\n", $sample, $sample_correct_umis, $sample_incorrect_umis, $sample_singleton_umis, join(", ", @marker_count_umis));
			} else {
				printf("Sample '%s': %8d correct UMIs found, %d incorrect and %d singleton UMIs discarded - %s.\n", $sample, $sample_correct_umis, $sample_incorrect_umis, $sample_singleton_umis, join(", ", @marker_count_umis));
			}
		}
		if (defined($INP_verbose)) {
			foreach my $sample (@$samples){
				foreach my $marker (@$markers){
					# Prints into a file the histogram of the number of variants per UMI
					if (defined($seq_stats->{$sample}{$marker}{'variants'})){
						my $output_variant_stats = sprintf("%s\t%s\n", "Variants", "UMIs_count");
						my @sorted_variants_per_umi = sort {$a<=>$b} keys %{$seq_stats->{$sample}{$marker}{'variants'}};
						foreach my $variants_number ($sorted_variants_per_umi[0]..$sorted_variants_per_umi[-1]){
							if (defined($seq_stats->{$sample}{$marker}{'variants'}{$variants_number})){
								$output_variant_stats .= sprintf("%d\t%d\n", $variants_number, $seq_stats->{$sample}{$marker}{'variants'}{$variants_number});
							} else {
								$output_variant_stats .= sprintf("%d\t%d\n", $variants_number, 0);
							}
						}
						write_to_file("$INP_outpath/$sample.$marker.variants_umis.txt",$output_variant_stats);
					}
					# Prints into a file the histogram of the number of reads per UMI
					if (defined($seq_stats->{$sample}{$marker}{'reads'})){
						my $output_read_stats = sprintf("%s\t%s\n", "Sequences", "UMIs_count");
						my @sorted_seqs_per_umi = sort {$a<=>$b} keys %{$seq_stats->{$sample}{$marker}{'reads'}};
						foreach my $variants_number ($sorted_seqs_per_umi[0]..$sorted_seqs_per_umi[-1]){
							if (defined($seq_stats->{$sample}{$marker}{'reads'}{$variants_number})){
								$output_read_stats .= sprintf("%d\t%d\n", $variants_number, $seq_stats->{$sample}{$marker}{'reads'}{$variants_number});
							} else {
								$output_read_stats .= sprintf("%d\t%d\n", $variants_number, 0);
							}
						}
						write_to_file("$INP_outpath/$sample.$marker.seqs_umis.txt",$output_read_stats);
					}
				}
			}
		}
	}

# 	# Annotates Variable and Joining regions
# 	my $tcr_segment_alleles;
# 	if (%INP_ref_files) {
# # 		if ($INP_verbose) {
# # 			print "\n";
# # 		}
# 		foreach my $sample (@$samples){
# 			my @tcr_segments_ = nsort(keys %INP_ref_files);
# # 			if ($INP_verbose) {
# # 			printf("Matching %s %s alleles for %s sample\n", join('/',@$markers), uc(join('/',@tcr_segments_)), $sample);
# # 			}
# 			my $matched_alleles;
# 			foreach my $marker (@$markers){
# 				if (!defined($tcr_segments->{$sample}{$marker})){
# 					next;
# 				}
# 				# Matches reference sequences of each TCR region
# 				foreach my $tcr_segment (@tcr_segments_) {
# 					# The CDR3 variants will be the indexes and the TCR segment sequences the values
# 					my %cdr3_to_segment_seqs;
# 					foreach my $cdr3_seq (keys %{$seq_clusters->{$sample}{$marker}}){
# 						if (!defined($tcr_segments->{$sample}{$marker}{$cdr3_seq}{$tcr_segment})){
# 							next;
# 						}
# 						$cdr3_to_segment_seqs{$cdr3_seq} = most_frequent(@{$tcr_segments->{$sample}{$marker}{$cdr3_seq}{$tcr_segment}});
# 					}
# 					# TCR region allele matching parameters
# 					my $INP_allele_align;
# 	# 				if (scalar @{$tcr_segments->{$sample}{$marker}{$cdr3_seq}{$tcr_segment}} < 100000) {
# 					# THESE PARAMETERS ARE OPTIMAL, IF THERE ARE TCRV LOCI MISSING IS BECAUSE V REFERENCE SEQUENCES ARE TOO SHORT (TAKEN FROM OLD VERSION OF MIXCR)
# 					$INP_allele_align = { 'alignment' => 'dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 90', 'aligned' => 0.1, 'ident' => 0.9 };
# 					my $cdr3_to_tcr_ref_name = match_alleles($INP_ref_files{$tcr_segment},\%cdr3_to_segment_seqs,undef,$INP_allele_align,$INP_threads);
# 	# 				} else {
# 	# 					if (!defined($INP_threads)){
# 	# 						$INP_allele_align = { 'alignment' => 'dna bowtie2 --sensitive-local -k 2', 'aligned' => 0.1, 'ident' => 0.9 };
# 	# 					} else {
# 	# 						$INP_allele_align = { 'alignment' => "dna bowtie2 --sensitive-local -k 2 --threads $INP_threads", 'aligned' => 0.1, 'ident' => 0.9 };
# 	# 					}
# 	# 					$cdr3_to_tcr_ref_name->{$tcr_segment} = match_alleles($INP_ref_files{$tcr_segment},\%cdr3_to_segment_seqs,undef,$INP_allele_align);
# 	# 				}
# 	# 				print '';
# 					foreach my $cdr3_seq (keys %{$seq_clusters->{$sample}{$marker}}){
# 						if (defined($cdr3_to_tcr_ref_name->{$cdr3_seq})){
# 							# 2 precision digit for 'tcrv' region, Ex. BV_TCRBV13-01*01 BV_TCRBV31*01
# 							# 2 precision digits for 'tcrj' region, Ex. BV_TCRBJ1-4*01
# 							my @alleles;
# 							foreach my $allele (split(/\s*\|\s*/,$cdr3_to_tcr_ref_name->{$cdr3_seq})){
# 								# if ($allele =~ /(\w+\d+)/){ # family
# 								# if ($allele =~ /(\w+[\d\-\*\/]+)/){; # allele (full precision)
# 								if ($allele =~ /([^\*]+)/){ # locus (before asterisk)
# 									push(@alleles, $1);
# 								} else {
# 									push(@alleles,$allele);
# 								}
# 							}
# 							if (@alleles) {
# 								$tcr_segment_alleles->{$sample}{$marker}{$cdr3_seq}{$tcr_segment} = join(',',unique(nsort(@alleles)));
# 								$matched_alleles->{$tcr_segment}++;
# 							}
# 						}
# 					}
# 				}
# 			}
# 			my $first_tcr_segment = shift @tcr_segments_;
# 			my @matched_alleles = ( sprintf("%8d %s",$matched_alleles->{$first_tcr_segment},uc($first_tcr_segment)) );
# 			push(@matched_alleles, map sprintf("%d %s",$matched_alleles->{$_},uc($_)) , @tcr_segments_);
# 			printf("Sample '%s': %s alleles.\n", $sample, join(", ",@matched_alleles));
# 
# 		}
# 
# 	}

	# Annotates Variable and Joining regions
	my ($tcr_segment_alleles, @tcr_segments_in_refs);
	if (%INP_ref_files) {
		@tcr_segments_in_refs = reverse(nsort(keys %INP_ref_files));
	}

	# Formats CDR3 clusters into sequences and headers
	$count_cdr3 = 0;
	undef($cdr3_seqs);
	undef($cdr3_headers);
	my %count_total;
	foreach my $sample (@$samples){
		my $sample_count_cdr3s = 0;
		my %count_singletons;
		my $matched_alleles;
		foreach my $marker (@$markers){
			# Counts the number of corrected variants per UMI
			my %variants_per_cluster;
			# Counts the number of sequences per cluster
			my %seqs_per_cluster;

			my @clustered_seqs =  sort { scalar @{$seq_clusters->{$sample}{$marker}{$b}} <=> scalar @{$seq_clusters->{$sample}{$marker}{$a}} } keys %{$seq_clusters->{$sample}{$marker}};

			# Annotates number of cluster members/corrected variants
			$variants_per_cluster{scalar @clustered_seqs}++;

			# Defines the max. length of the sequence IDs
			my $len_id = length(scalar @clustered_seqs);
			# Annotates the protein sequences for the clustered CDR3 sequences
			# Generates unique sequence IDs for the CDR3 sequences that translate into the same protein
			my %prot_seq_ids;
			my (%count_prot_seqs, %count_dna_seqs);
			my $count_seqs = 0;
			foreach my $seq (@clustered_seqs){
				if (!defined($cdr3_prots{$seq})){
					$cdr3_prots{$seq} = dna_to_prot($seq);
				}
				if (!defined($count_prot_seqs{$cdr3_prots{$seq}})){
					$count_seqs++;
					$prot_seq_ids{$cdr3_prots{$seq}} = sprintf("%s-%s%0".$len_id."d", $sample, $marker, $count_seqs);
				}
				$count_prot_seqs{$cdr3_prots{$seq}}++;
				$count_dna_seqs{$seq} = $count_prot_seqs{$cdr3_prots{$seq}};
			}

			# Matches reference sequences of each TCR region
			if (%INP_ref_files) {
				foreach my $tcr_segment (@tcr_segments_in_refs) {
					# The CDR3 variants will be the indexes and the TCR segment sequences the values
					my %cdr3_to_segment_seqs;
					foreach my $cdr3_seq (@clustered_seqs){
						if (!defined($tcr_segments->{$sample}{$marker}{$cdr3_seq}{$tcr_segment})){
							next;
						}
						if (!$INP_singletons && scalar @{$seq_clusters->{$sample}{$marker}{$cdr3_seq}} == 1) {
							next;
						}
						$cdr3_to_segment_seqs{$cdr3_seq} = most_frequent(@{$tcr_segments->{$sample}{$marker}{$cdr3_seq}{$tcr_segment}});
					}
					# TCR region allele matching parameters
					my $INP_allele_align;
	# 				if (scalar @{$tcr_segments->{$sample}{$marker}{$cdr3_seq}{$tcr_segment}} < 100000) {
					# THESE PARAMETERS ARE OPTIMAL, IF THERE ARE TCRV LOCI MISSING IS BECAUSE V REFERENCE SEQUENCES ARE TOO SHORT (TAKEN FROM OLD VERSION OF MIXCR)
					$INP_allele_align = { 'alignment' => 'dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 90', 'aligned' => 0.1, 'ident' => 0.9 };
					my $cdr3_to_tcr_ref_name = match_alleles($INP_ref_files{$tcr_segment},\%cdr3_to_segment_seqs,undef,$INP_allele_align,$INP_threads);
	# 				} else {
	# 					if (!defined($INP_threads)){
	# 						$INP_allele_align = { 'alignment' => 'dna bowtie2 --sensitive-local -k 2', 'aligned' => 0.1, 'ident' => 0.9 };
	# 					} else {
	# 						$INP_allele_align = { 'alignment' => "dna bowtie2 --sensitive-local -k 2 --threads $INP_threads", 'aligned' => 0.1, 'ident' => 0.9 };
	# 					}
	# 					$cdr3_to_tcr_ref_name->{$tcr_segment} = match_alleles($INP_ref_files{$tcr_segment},\%cdr3_to_segment_seqs,undef,$INP_allele_align);
	# 				}
	# 				print '';
					foreach my $cdr3_seq (@clustered_seqs){
						if (defined($cdr3_to_tcr_ref_name->{$cdr3_seq})){
							# 2 precision digit for 'tcrv' region, Ex. BV_TCRBV13-01*01 BV_TCRBV31*01
							# 2 precision digits for 'tcrj' region, Ex. BV_TCRBJ1-4*01
							my @alleles;
							foreach my $allele (split(/\s*\|\s*/,$cdr3_to_tcr_ref_name->{$cdr3_seq})){
								# if ($allele =~ /(\w+\d+)/){ # family
								# if ($allele =~ /(\w+[\d\-\*\/]+)/){; # allele (full precision)
								if ($allele =~ /([^\*]+)/){ # locus (before asterisk)
									push(@alleles, $1);
								} else {
									push(@alleles,$allele);
								}
							}
							if (@alleles) {
								$tcr_segment_alleles->{$sample}{$marker}{$cdr3_seq}{$tcr_segment} = join(',',unique(@alleles));
								$matched_alleles->{$tcr_segment}++;
							}
						}
					}
				}
			}

			# Annotates the clustered CDR3 variant information
			my $count_prot_seqs = 0;
			my %annotated_seqs;
			$count_singletons{$marker} = 0;
			foreach my $seq (@clustered_seqs){
				my %count_umis = count_umis($seq_clusters->{$sample}{$marker}{$seq});
				my $count_umis = scalar keys %count_umis;
				my $depth = scalar @{$seq_clusters->{$sample}{$marker}{$seq}};
				# Annotates number of sequences per cluster
				$seqs_per_cluster{$depth}++;
				if ($depth == 1){
					$count_singletons{$marker}++;
					if (!$INP_singletons) {
						next;
					}
				}
				$count_total{$sample}{$marker} += $depth;
				my $md5 = generate_md5($seq);
				my $len = length($seq);
				my $count_samples = 0;
				foreach my $sample2 (@$samples){
					if (defined($seq_clusters->{$sample2}{$seq})){ # $sample ne $sample2 && 
						$count_samples++;
					}
				}
				my $seq_id;
				if ($count_prot_seqs{$cdr3_prots{$seq}} == 1){
					$seq_id = $prot_seq_ids{$cdr3_prots{$seq}};
				} else {
					$seq_id = sprintf("%s.%d", $prot_seq_ids{$cdr3_prots{$seq}}, $count_dna_seqs{$seq});
				}
				my $header = $seq_clusters->{$sample}{$marker}{$seq}[0];
				$header =~ s/^[^\|]+/$seq_id /;
				$header =~ s/ \| umi=\S+//;
				my $alleles = '';
				if (defined($tcr_segment_alleles) && defined($tcr_segment_alleles->{$sample}{$marker}{$seq})){
					my @alleles;
					foreach my $tcr_segment (keys %{$tcr_segment_alleles->{$sample}{$marker}{$seq}}){
						# my $allele = $tcr_segment_alleles->{$sample}{$marker}{$seq}{$tcr_segment};
						# $header =~ s/ \| $tcr_segment=\S+/ | $tcr_segment=$allele/;
						$header =~ s/ \| $tcr_segment=\S+//;
						push(@alleles,sprintf("%s=%s",$tcr_segment,$tcr_segment_alleles->{$sample}{$marker}{$seq}{$tcr_segment}));
					}
					$alleles = "\| ".join(" \| ", @alleles)." ";
				}
				push(@{$cdr3_seqs->{$sample}{$marker}}, $seq);
				push(@{$cdr3_headers->{$sample}{$marker}}, sprintf("%s | hash=%s %s| len=%d | count=%d | depth=%d | samples=%d", $header, $md5, $alleles, $len, $count_umis, $depth, $count_samples));
				$annotated_seqs{$seq}++;
			}

			# Prints a file with variants per cluster statistics
			if (defined($INP_verbose) && %seqs_per_cluster){
				my $variants_per_cluster_stats = sprintf("%s\t%s\n", "Variants", "Clusters_count");
				my @sorted_variants_per_cluster = sort {$a<=>$b} keys %seqs_per_cluster;
				foreach my $variants_number ($sorted_variants_per_cluster[0]..$sorted_variants_per_cluster[-1]){
					if (defined($seqs_per_cluster{$variants_number})){
						$variants_per_cluster_stats .= sprintf("%d\t%d\n", $variants_number, $seqs_per_cluster{$variants_number});
					} else {
						$variants_per_cluster_stats .= sprintf("%d\t%d\n", $variants_number, 0);
					}
				}
				write_to_file("$INP_outpath/$sample.$marker.variants_clusters.txt",$variants_per_cluster_stats);

				# Prints a file with sequences per cluster statistics
				my $seqs_per_cluster_stats = sprintf("%s\t%s\n", "Sequences", "Clusters_count");
				my @sorted_seqs_per_cluster = sort {$a<=>$b} keys %seqs_per_cluster;
				foreach my $variants_number ($sorted_seqs_per_cluster[0]..$sorted_seqs_per_cluster[-1]){
					if (defined($seqs_per_cluster{$variants_number})){
						$seqs_per_cluster_stats .= sprintf("%d\t%d\n", $variants_number, $seqs_per_cluster{$variants_number});
					} else {
						$seqs_per_cluster_stats .= sprintf("%d\t%d\n", $variants_number, 0);
					}
				}
				write_to_file("$INP_outpath/$sample.$marker.seqs_clusters.txt",$seqs_per_cluster_stats);
			}

			# Write sequences into file
			my $filename;
			if (defined($cdr3_headers->{$sample}{$marker}) && @{$cdr3_headers->{$sample}{$marker}}){
				my $file_prefix = "$INP_outpath/$sample.$marker.clustered.cdr3";
				if (!defined($INP_onlystats)){
					if ($INP_output_format eq 'fastq'){
						$filename = create_fastq_file($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}, $cdr3_qualities->{$sample}{$marker}, "$file_prefix.fq", 'gzip');
					} else {
						$filename = create_fasta_file($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}, "$file_prefix.fa", 'gzip');
					}
				} else {
				}
				# Counts clustered CDR3 sequences and calculate statistics
				$output_files->{'clustered'}{$sample}{$marker} = $filename;
				$count_cdr3 += scalar @{$cdr3_headers->{$sample}{$marker}};
				$sample_count_cdr3s += scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_global->{$sample}{$marker}{'clustered'} = scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_global->{$sample}{$marker}{'singletons'} = $count_singletons{$marker};
				$stats_cdr3_lengths->{$sample}{$marker}{'clustered'} = count_lengths($cdr3_seqs->{$sample}{$marker});
				$stats_cdr3_depths->{$sample}{$marker}{'clustered'} = count_depths($cdr3_headers->{$sample}{$marker});
				$stats_tcr_regions->{$sample}{$marker}{'clustered'} = count_tcr_regions($cdr3_headers->{$sample}{$marker},$cdr3_seqs->{$sample}{$marker});
				print '';
			} else {
				$stats_cdr3_global->{$sample}{$marker}{'clustered'} = 0;
			}
		}
		# Prints statistics
		my @marker_count_cdr3s = map sprintf("%s: %d",$_, $stats_cdr3_global->{$sample}{$_}{'clustered'}) , @$markers;
		my $total_sample_singletons = 0;
		map $total_sample_singletons += $count_singletons{$_}, @$markers;
		my @marker_count_singletons = map sprintf("%s: %d",$_, $count_singletons{$_}) , @$markers;
		if (!defined($INP_onlystats)){
			if (defined($INP_umi_cluster)) {
				printf("Sample '%s': %8d clustered CDR3 sequences (%s) written into 'clustered.cdr3' files.\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s)); # , join(", ", values %{$output_files->{'clustered'}{$sample}})
			} elsif ($INP_singletons) {
				printf("Sample '%s': %8d clustered CDR3 sequences (%s) written into 'clustered.cdr3' files, kept %5d singletons (%s).\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s), $total_sample_singletons, join(", ",@marker_count_singletons)); # , join(", ", values %{$output_files->{'clustered'}{$sample}})
			} else {
				printf("Sample '%s': %8d clustered CDR3 sequences (%s) written into 'clustered.cdr3' files, removed %5d singletons (%s).\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s), $total_sample_singletons, join(", ",@marker_count_singletons)); # , join(", ", values %{$output_files->{'clustered'}{$sample}})
			}
		} else {
			if (defined($INP_umi_cluster)) {
				printf("Sample '%s': %8d clustered CDR3 sequences (%s).\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s));
			} elsif ($INP_singletons) {
				printf("Sample '%s': %8d clustered CDR3 sequences (%s), kept %5d singletons (%s).\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s), $total_sample_singletons, join(", ",@marker_count_singletons));
			} else {
				printf("Sample '%s': %8d clustered CDR3 sequences (%s), removed %5d singletons (%s).\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s), $total_sample_singletons, join(", ",@marker_count_singletons));
			}
		}
		# Prints allele stats
		if (%INP_ref_files) {
			my @matched_alleles;
			foreach my $tcr_segment (@tcr_segments_in_refs) {
				if (defined($matched_alleles->{$tcr_segment})){
					push(@matched_alleles, sprintf("%8d %s",$matched_alleles->{$tcr_segment},uc($tcr_segment)));
				} else {
					push(@matched_alleles, sprintf("%8d %s",0,uc($tcr_segment)));
				}
			}
			printf("Sample '%s': %s alleles.\n", $sample, join(", ",@matched_alleles));
		}
	}
} else {
	# Reads previously clustered CDR3 sequences
	printf("\nReading previously clustered CDR3 sequences.\n\n",);
	$count_cdr3 = 0;
	foreach my $sample (@$samples){
		foreach my $marker (@$markers){
			my $filename;
			# Read sequences from previous existing FASTA/FASTQ GZIPPED file
			my $file_prefix = $INP_read_files[0]."/$sample.$marker.clustered.cdr3";
			if (-e "$file_prefix.fq.gz"){
				$filename = "$file_prefix.fq.gz";
			} elsif (-e "$file_prefix.fa.gz"){
				$filename = "$file_prefix.fa.gz";
			}
			if (defined($filename)){
				($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}) = read_sequence_file($filename);
				printf("Sample '%s': %8d clustered %s CDR3 sequences read from '%s'.\n", $sample, scalar @{$cdr3_headers->{$sample}{$marker}}, $marker, $filename);
			} else {
				printf("Sample '%s': %8s clustered %s CDR3 sequence file found.\n", $sample, 'no', $marker);
				next;
			}
			if (!defined($cdr3_headers->{$sample}{$marker})){
				printf("Sample '%s': %8d clustered %s CDR3 sequences read from '%s'.\n", $sample, 0, $marker, $filename);
				next;
			}
			# Counts clustered CDR3 sequences and calculate statistics
			$output_files->{'clustered'}{$sample}{$marker} = $filename;
			$count_cdr3 += scalar @{$cdr3_headers->{$sample}{$marker}};
			$stats_cdr3_global->{$sample}{$marker}{'clustered'} = scalar @{$cdr3_headers->{$sample}{$marker}};
			$stats_cdr3_lengths->{$sample}{$marker}{'clustered'} = count_lengths($cdr3_seqs->{$sample}{$marker});
			$stats_cdr3_depths->{$sample}{$marker}{'clustered'} = count_depths($cdr3_headers->{$sample}{$marker});
			$stats_tcr_regions->{$sample}{$marker}{'clustered'} = count_tcr_regions($cdr3_headers->{$sample}{$marker},$cdr3_seqs->{$sample}{$marker});
			$stats_cdr3_global->{$sample}{$marker}{'singletons'} = $stats_cdr3_depths->{$sample}{$marker}{'clustered'}{'1'};
		}
	}
}
if (defined($INP_umi_cluster)) {
	foreach my $sample (@$samples){
		my @marker_count_umis;
		my $total_sample_umis = 0;
		foreach my $marker (@$markers){
			my %count_umis = count_umis($cdr3_headers->{$sample}{$marker});
			my $count_umis = 0;
			if (%count_umis) {
				$count_umis = sum(values %count_umis);
			}
			$total_sample_umis += $count_umis;
			push(@marker_count_umis, sprintf("%s: %d",$marker,$count_umis));
		}
		printf("Sample '%s': %8d Unique Molecular Identifiers found (%s).\n", $sample, $total_sample_umis, join(", ",@marker_count_umis));
	}
}

printf("\n%d CDR3 sequences after clustering.\n", $count_cdr3);

# Stores in MIGEC/MIXCR format the CDR3 clonotype information after clustering
my $results;
foreach my $sample (@$samples){
	# my $file_prefix = "$INP_outpath/$sample.clones.cdr3";
	my $filename = "$INP_outpath/$sample.clones.cdr3.txt";
	my @fields;
	if (!defined($INP_umi_cluster)) {
		@fields = ('cloneId', 'cloneCount', 'cloneFraction', 'nSeqCDR3', 'aaSeqCDR3', 'allVHitsWithScore', 'allJHitsWithScore');
	} else {
		@fields = ('cloneId', 'cloneCount', 'cloneFraction', 'totalReads', 'readFraction', 'nSeqCDR3', 'aaSeqCDR3', 'allVHitsWithScore', 'allJHitsWithScore');
	}
	my $output = join("\t",@fields)."\n";
# 	my $len_id = 0;
# 	foreach my $marker (@$markers){
# 		if (defined($cdr3_headers->{$sample}{$marker}) && $len_id < length(scalar @{$cdr3_headers->{$sample}{$marker}})){
# 			$len_id = length(scalar @{$cdr3_headers->{$sample}{$marker}});
# 		}
# 	}
	my ($total_reads, $total_umis) = (0,0);
	foreach my $marker (@$markers){
		if (!defined($cdr3_headers->{$sample}{$marker})){
			next;
		}
		for (my $i=0; $i<=$#{$cdr3_headers->{$sample}{$marker}}; $i++){
			my %result;
			#$result{'cloneId'} = sprintf("%s-%s%0".$len_id."d",$sample,$marker,$i+1);
			$result{'nSeqCDR3'} = $cdr3_seqs->{$sample}{$marker}[$i];
			if (defined($cdr3_prots{$result{'nSeqCDR3'}})){
				$result{'aaSeqCDR3'} = $cdr3_prots{$result{'nSeqCDR3'}};
			} else {
				$result{'aaSeqCDR3'} = dna_to_prot($cdr3_seqs->{$sample}{$marker}[$i]);
			}
			my $header = $cdr3_headers->{$sample}{$marker}[$i];
			# s092_A-TCRB0001 | ampli=TCRB-FWD_RACE-REV_RACE | tcrj=TRBJ2 | hash=d61bd6ea94be0874ff2c129e731e35e3 | len=39 | count=9 | depth=9 | samples=0
			my ($ampli, $umi, $tcrv, $tcrj, $md5, $len, $count, $depth, $samples) = ('','','','','','','','','');
			if  ($header =~ /^(\S+)/){ $result{'cloneId'} = $1 }
			if  ($header =~ /ampli=(\S+)/){ $result{'amplicon'} = $1 }
			if  ($header =~ /umi=(\S+)/){ $result{'umi'} = $1 }
			if  ($header =~ /tcrv=(\S+)/ && $1 !~ /^[ACGT\-]+$/){ $result{'allVHitsWithScore'} = $1 }
			if  ($header =~ /tcrj=(\S+)/ && $1 !~ /^[ACGT\-]+$/){ $result{'allJHitsWithScore'} = $1 }
			if  ($header =~ /(hash|md5)=(\S+)/){ $result{'md5'} = $2 }
			if  ($header =~ /len=(\S+)/){ $result{'length'} = $1 }
			if  ($header =~ /samples=(\S+)/){ $result{'samples'} = $1 }
			# MIXCR 'cloneCount' is the equivalent to MIGEC 'totalReads', same for 'CloneFraction' and 'readFraction'
			if (!defined($INP_umi_cluster) && $header =~ /depth=(\S+)/){
				$result{'cloneCount'} = $1;
				$total_reads += $result{'cloneCount'};
			} elsif (defined($INP_umi_cluster) && $header =~ /count=(\S+)/){
				$result{'cloneCount'} = $1;
				$total_umis += $result{'cloneCount'};
				if ($header =~ /depth=(\S+)/){
					$result{'totalReads'} = $1;
					$total_reads += $result{'totalReads'};
				}
			}
			push(@{$results->{$sample}},\%result);
		}
	}
	foreach my $result (sort {$b->{'cloneCount'}<=>$a->{'cloneCount'}} @{$results->{$sample}}){
		# MIXCR 'cloneCount' is the equivalent to MIGEC 'totalReads', same for 'CloneFraction' and 'readFraction'
		if (!defined($INP_umi_cluster)){
			$result->{'cloneFraction'} = sprintf("%.6f",$result->{'cloneCount'}/$total_reads);
		} else {
			$result->{'cloneFraction'} = sprintf("%.6f",$result->{'cloneCount'}/$total_umis);
			$result->{'readFraction'} = sprintf("%.6f",$result->{'totalReads'}/$total_reads);
		}
		my @values;
		foreach my $field (@fields) {
			if (defined($result->{$field})){
				push(@values,$result->{$field});
			} else {
				push(@values,'');
			}
		}
		$output .= join("\t",@values)."\n";
	}
	write_to_file($filename,$output);
}


# Calculates common CDR3 sequences between different samples from the same individual
printf("\nCalculating common CDR3 sequences between different samples from the same individual.\n");
$stats_cdr3_common->{'clustered'}{'intraindividual'} = common_cdr3_seqs($cdr3_seqs,$INP_cluster_errors,'intraindividual');
my $individuals = {};
foreach my $marker (@$markers) {
	if (!defined($stats_cdr3_common->{'clustered'}{'intraindividual'}{$marker})){
		next;
	}
	$individuals->{$marker} = [ nsort(keys %{$stats_cdr3_common->{'clustered'}{'intraindividual'}{$marker}}) ];
	foreach my $individual (@{$individuals->{$marker}}) {
		if (!defined($stats_cdr3_common->{'clustered'}{'intraindividual'}{$marker}{$individual}{'total'})){
			next;
		}
		foreach my $count_samples (sort {$b<=>$a} keys %{$stats_cdr3_common->{'clustered'}{'intraindividual'}{$marker}{$individual}{'total'}}) {
			printf("Marker '%s', individual '%s': %8d sequences are shared by %d samples.\n", $marker, $individual, $stats_cdr3_common->{'clustered'}{'intraindividual'}{$marker}{$individual}{'total'}{$count_samples}, $count_samples);
		}
	}
}
# # VERY SLOW:
# printf("\nCalculating common CDR3 sequences between all samples from different individuals.\n");
# $stats_cdr3_common->{'clustered'}{'interindividual'} = common_cdr3_seqs($cdr3_seqs,$INP_cluster_errors,'interindividual');
# foreach my $marker (@$markers) {
# 	if (!defined($stats_cdr3_common->{'clustered'}{'interindividual'}{$marker}) || scalar @{$individuals->{$marker}} < 2 || !defined($stats_cdr3_common->{'clustered'}{'interindividual'}{$marker}{'all'})){
# 		next;
# 	}
# 	foreach my $count_individuals (sort {$b<=>$a} keys %{$stats_cdr3_common->{'clustered'}{'interindividual'}{$marker}{'all'}}) {
# 		printf("Marker '%s': %8d amino acid sequences are shared by %d individuals.\n", $marker, $stats_cdr3_common->{'clustered'}{'interindividual'}{$marker}{'all'}{$count_individuals}, $count_individuals);
# 	}
# }

# Calculates common CDR3 sequences between different individuals
# printf("\nCalculating common CDR3 sequences between all samples from different individuals.\n\n");
# $stats_cdr3_common->{'clustered'}{'interindividual'} = common_cdr3_seqs($cdr3_seqs,undef,'interindividual');

# DEBUGGING:
# Stops reading previous results
# undef($INP_cdr3_path);

# Translates CDR3 sequences to amino acids (CDR3s are extracted to start the reading frame from the 1st nucleotide)
if (!defined($INP_cdr3_path)) {
	printf("\nTranslating %d CDR3 sequences.\n",$count_cdr3);
	# Stores translated CDR3 sequences
	foreach my $sample (@$samples){
		foreach my $marker (@$markers){
			if (defined($cdr3_seqs->{$sample}{$marker})){
				($cdr3_seqs->{$sample}{$marker},$cdr3_headers->{$sample}{$marker}) = translate_seqs($cdr3_seqs->{$sample}{$marker},$cdr3_headers->{$sample}{$marker},'1',['nostop','unique']);
				# Converts the DNA IDs (TCRB0001.1, TCRB0001.2, TCRB0001.3) into protein IDs (TCRB0001)
				# s092_A-TCRB0001.4 | ampli=TCRB-FWD_RACE-REV_RACE | tcrj=TRBJ2 | hash=d61bd6ea94be0874ff2c129e731e35e3 | len=39 | count=9 | depth=9 | samples=0
				for (my $i=0; $i<=$#{$cdr3_headers->{$sample}{$marker}}; $i++){
					my $md5 = generate_md5($cdr3_seqs->{$sample}{$marker}[$i]);
					$cdr3_headers->{$sample}{$marker}[$i] =~ s/^(\S+)\.\d+/$1/;
					$cdr3_headers->{$sample}{$marker}[$i] =~ s/(hash|md5)=(\S+)/$1=$md5/;
				}
			}
		}
	}
}
# Stores translated CDR3 sequences
print "\n";
$count_cdr3 = 0;
if (!defined($INP_cdr3_path)) {
	foreach my $sample (@$samples){
		my $sample_count_cdr3s = 0;
		foreach my $marker (@$markers){
			my $filename;
			# Write sequences into file
			if (defined($cdr3_headers->{$sample}{$marker}) && @{$cdr3_headers->{$sample}{$marker}}){
				if (!defined($INP_onlystats)){
					my $file_prefix = "$INP_outpath/$sample.$marker.prot.cdr3";
					if ($INP_output_format eq 'fastq'){
						$filename = create_fastq_file($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}, $cdr3_qualities->{$sample}{$marker}, "$file_prefix.fq", 'gzip');
					} else {
						$filename = create_fasta_file($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}, "$file_prefix.fa", 'gzip');
					}
				}
				# Count translated CDR3 sequences and calculate statistics
				$output_files->{'translated'}{$sample}{$marker} = $filename;
				$count_cdr3 += scalar @{$cdr3_headers->{$sample}{$marker}};
				$sample_count_cdr3s += scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_global->{$sample}{$marker}{'translated'} = scalar @{$cdr3_headers->{$sample}{$marker}};
				$stats_cdr3_lengths->{$sample}{$marker}{'translated'} = count_lengths($cdr3_seqs->{$sample}{$marker});
				$stats_tcr_regions->{$sample}{$marker}{'translated'} = count_tcr_regions($cdr3_headers->{$sample}{$marker},$cdr3_seqs->{$sample}{$marker});
			} else {
				$stats_cdr3_global->{$sample}{$marker}{'translated'} = 0;
			}
		}
		# Prints statistics
		my @marker_count_cdr3s = map sprintf("%s: %d",$_, $stats_cdr3_global->{$sample}{$_}{'translated'}) , @$markers;
		if (!defined($INP_onlystats)){
			printf("Sample '%s': %8d translated CDR3 sequences (%s) written into '.prot.cdr3' files.\n", $sample, $sample_count_cdr3s, join(", ", @marker_count_cdr3s)); # , join(", ", values %{$output_files->{'translated'}{$sample}})
		} else {
			printf("Sample '%s': %8d translated CDR3 sequences (%s).\n", $sample, $sample_count_cdr3s, join(", ",@marker_count_cdr3s));
		}
	}
} else {
	# Reads previously translated CDR3 sequences
	printf("\nReading previously clustered CDR3 sequences.\n\n",);
	foreach my $sample (@$samples){
		foreach my $marker (@$markers){
			my $filename;
			# Read sequences from previous existing FASTA/FASTQ GZIPPED file
			my $file_prefix = $INP_read_files[0]."/$sample.$marker.prot.cdr3";
			if (-e "$file_prefix.fq.gz"){
				$filename = "$file_prefix.fq.gz";
			} elsif (-e "$file_prefix.fa.gz"){
				$filename = "$file_prefix.fa.gz";
			}
			if (defined($filename)){
				($cdr3_seqs->{$sample}{$marker}, $cdr3_headers->{$sample}{$marker}) = read_sequence_file($filename);
				printf("Sample '%s': %8d translated %s CDR3 sequences read from '%s'.\n", $sample, scalar @{$cdr3_headers->{$sample}{$marker}}, $marker, $filename);
			} else {
				printf("Sample '%s': %8s translated %s CDR3 sequence file found.\n", $sample, 'no', $marker);
				next;
			}
			if (!defined($cdr3_headers->{$sample}{$marker})){
				printf("Sample '%s': %8d translated %s CDR3 sequences read from '%s'.\n", $sample, 0, $marker, $filename);
			}
			# Counts translated CDR3 sequences and calculate statistics
			$output_files->{'translated'}{$sample}{$marker} = $filename;
			$count_cdr3 += scalar @{$cdr3_headers->{$sample}{$marker}};
			$stats_cdr3_global->{$sample}{$marker}{'translated'} = scalar @{$cdr3_headers->{$sample}{$marker}};
			$stats_cdr3_lengths->{$sample}{$marker}{'translated'} = count_lengths($cdr3_seqs->{$sample}{$marker});
			$stats_cdr3_depths->{$sample}{$marker}{'translated'} = count_depths($cdr3_headers->{$sample}{$marker});
			$stats_tcr_regions->{$sample}{$marker}{'translated'} = count_tcr_regions($cdr3_headers->{$sample}{$marker},$cdr3_seqs->{$sample}{$marker});
		}
	}
}

printf("\n%d CDR3 sequences after translating.\n", $count_cdr3);

# Calculates common CDR3 sequences between different samples from the same individual
my @analysis_steps_to_print_ = ('clustered','translated');
if (in_array(\@analysis_steps_to_print_,'translated')){
	printf("\nCalculating common CDR3 amino acid sequences between different samples from the same individual.\n\n");
	$stats_cdr3_common->{'translated'}{'intraindividual'} = common_cdr3_seqs($cdr3_seqs,$INP_cluster_errors,'intraindividual');
	foreach my $marker (@$markers) {
		if (!defined($stats_cdr3_common->{'translated'}{'intraindividual'}{$marker})){
			next;
		}
		foreach my $individual (@{$individuals->{$marker}}) {
			foreach my $count_samples (sort {$b<=>$a} keys %{$stats_cdr3_common->{'translated'}{'intraindividual'}{$marker}{$individual}{'total'}}) {
				printf("Marker '%s', individual '%s': %8d amino acid sequences are shared by %d samples.\n", $marker, $individual, $stats_cdr3_common->{'translated'}{'intraindividual'}{$marker}{$individual}{'total'}{$count_samples}, $count_samples);
			}
		}
	}
# 	# VERY SLOW:
# 	printf("\nCalculating common CDR3 sequences between all samples from different individuals.\n\n");
# 	$stats_cdr3_common->{'translated'}{'interindividual'} = common_cdr3_seqs($cdr3_seqs,$INP_cluster_errors,'interindividual');
# 	foreach my $marker (@$markers) {
# 		if (!defined($stats_cdr3_common->{'translated'}{'interindividual'}{$marker}) || scalar @{$individuals->{$marker}} < 2 || !defined($stats_cdr3_common->{'translated'}{'interindividual'}{$marker}{'all'})){
# 			next;
# 		}
# 		foreach my $count_individuals (sort {$b<=>$a} keys %{$stats_cdr3_common->{'translated'}{'interindividual'}{$marker}{'all'}}) {
# 			printf("Marker '%s': %8d amino acid sequences are shared by %d individuals.\n", $marker, $stats_cdr3_common->{'translated'}{'interindividual'}{$marker}{'all'}{$count_individuals}, $count_individuals);
# 		}
# 	}
}


# Creates an Excel file to print the statistics

my $workbook  = Excel::Writer::XLSX->new($excel_outputfile);
$workbook->set_properties(
	title    => "TCR CDR3 analysis results",
	author   => "Alvaro Sebastian",
	comments => "TCR CDR3 analysis results",
	company  => "Evolutionary Biology Group, Adam Mickiewicz University",
);
$workbook->compatibility_mode();
my $excel_style_bold = $workbook->add_format(bold => 1);
my @excel_red_colors = (
	$workbook->set_custom_color(40, '#ffffff'), # Lighter
	$workbook->set_custom_color(41, '#ffe5e5'),
	$workbook->set_custom_color(42, '#ffb2b2'),
	$workbook->set_custom_color(43, '#ff7f7f'),
	$workbook->set_custom_color(44, '#ff4c4c'), # Darker
);
my @excel_green_colors = (
	$workbook->set_custom_color(50, '#ffffff'), # Lighter
	$workbook->set_custom_color(51, '#e5f2e5'),
	$workbook->set_custom_color(52, '#b2d8b2'),
	$workbook->set_custom_color(53, '#7fbf7f'),
	$workbook->set_custom_color(54, '#4ca64c'), # Darker
);
my @excel_styles_red = map $workbook->add_format(bg_color => $_) , @excel_red_colors;
my @excel_styles_green = map $workbook->add_format(bg_color => $_) , @excel_green_colors;

# Prints global CDR3 statistics into Excel file
print "\nPrinting CDR3 extraction, filtering and clustering statistics.\n";
my $wsu = $workbook->add_worksheet("unique");
# Writes worksheet legend
my $wsu_row = 0;
$wsu->write($wsu_row, 0, "TCR CDR3 extraction, filtering and clustering statistics", $excel_style_bold); $wsu_row+=2;
if (defined($total_sample_reads) && %$total_sample_reads){
	my $total_reads = 0;
	if (defined($INP_multifile) || defined($INP_filepath)){
		map $total_reads += $total_sample_reads->{$_} , keys %$total_sample_reads;
		$wsu->write($wsu_row, 0, sprintf("Initial reads: %d (%d samples)", $total_reads, scalar @$samples), $excel_style_bold); $wsu_row+=2;
	} else {
		$total_reads = $total_sample_reads->{$samples->[0]};
		$wsu->write($wsu_row, 0, sprintf("Initial reads: %d (%d samples)", $total_reads, scalar @$samples), $excel_style_bold); $wsu_row+=2;
	}
}
# $wsu->write($wsu_row, 0, "Total: extracted CDR3 sequences before filtering and clustering", $excel_style_bold); $wsu_row++;
$wsu->write($wsu_row, 0, sprintf("Filtered: min_len: %d, max_len: %d, filter off-frame and stop-codon: %s", $INP_minlen,$INP_maxlen,$yes_no{$INP_noframe}), $excel_style_bold); $wsu_row++;
if (defined($INP_umi_cluster) && length($INP_umi_cluster)>1){
	$wsu->write($wsu_row, 0, sprintf("Clustered: UMIs: %s, max. subs. within UMI: %d, max. subs.: %d, keep singletons: %s",$INP_umi_cluster,$INP_umi_errors,$INP_cluster_errors,$yes_no{$INP_singletons}), $excel_style_bold); $wsu_row++;
} elsif (defined($INP_umi_cluster)){
	$wsu->write($wsu_row, 0, sprintf("Clustered: UMIs: %s, max. subs. within UMI: %d, max. subs.: %d, keep singletons: %s",$yes_no{$INP_umi_cluster},$INP_umi_errors,$INP_cluster_errors,$yes_no{$INP_singletons}), $excel_style_bold); $wsu_row++;
} else {
	$wsu->write($wsu_row, 0, sprintf("Clustered: max_subs: %d, keep singletons: %s",$INP_cluster_errors,$yes_no{$INP_singletons}), $excel_style_bold); $wsu_row++;
}
if ($INP_singletons){
	$wsu->write($wsu_row, 0, "Singletons: sequences of depth 1 included in the clustered group", $excel_style_bold); $wsu_row++;
} else {
	$wsu->write($wsu_row, 0, "Singletons: sequences of depth 1 removed during clustering", $excel_style_bold); $wsu_row++;
}
$wsu->write($wsu_row, 0, "Translated: translated sequences (removing the ones with stop codons)", $excel_style_bold); $wsu_row+=3;

my @cdr3_stats_fields=('total');
if ($INP_nseqs){
	push(@cdr3_stats_fields,'selected');
}
push(@cdr3_stats_fields,'filtered', 'singletons', 'clustered', 'translated');
my %stat_fields_headers = (	'total'=>'CDR3 reads',
				'selected'=>'Selected CDR3 reads',
				'filtered'=>'Filtered CDR3 reads',
				'clustered'=>'Clustered CDR3 variants',
				'translated'=>'Translated CDR3 variants',
			);
if ($INP_singletons){
	$stat_fields_headers{'singletons'} = 'Kept singleton CDR3 variants'
} else {
	$stat_fields_headers{'singletons'} = 'Removed singleton CDR3 variants'
}
$wsu->write_row($wsu_row, 0, [ 'Sample', 'Marker', 'Total reads', 'TCR reads', map($stat_fields_headers{$_},@cdr3_stats_fields) ] , $excel_style_bold); $wsu_row++;
my $wsu_col = 0;
my $chart_series;
$chart_series->{'all'}{'categories'}{'first'} = xl_col_to_name($wsu_col).($wsu_row+1);
$chart_series->{'all'}{'categories'}{'last'} = xl_col_to_name($wsu_col).($wsu_row + (scalar @$samples * scalar @$markers));
foreach my $marker (@$markers){
	foreach my $sample (@$samples){
		$wsu_col = 0;
		$wsu->write($wsu_row, $wsu_col, $sample, $excel_style_bold); $wsu_col++;
		$wsu->write($wsu_row, $wsu_col, $marker, $excel_style_bold); $wsu_col++;
		# If there are no reads for the sample, sets the number to 0 to avoid further errors
		if (!defined($total_sample_reads->{$sample})){
			$total_sample_reads->{$sample} = 0;
		}
		$wsu->write($wsu_row, $wsu_col, $total_sample_reads->{$sample}); $wsu_col++;
		$wsu->write($wsu_row, $wsu_col, $total_tcr_reads->{$sample}); $wsu_col++;
		foreach my $field (@cdr3_stats_fields){
			if (defined($stats_cdr3_global->{$sample}{$marker}{lc($field)})){
				$wsu->write($wsu_row, $wsu_col, $stats_cdr3_global->{$sample}{$marker}{$field});
			} else {
				$wsu->write($wsu_row, $wsu_col, 0);
			}
			if ($sample eq $samples->[0]){
				$chart_series->{lc($field)}{$marker}{'values'}{'first'} = xl_col_to_name($wsu_col).($wsu_row+1);
			} 
			if ($sample eq $samples->[-1]){
				$chart_series->{lc($field)}{$marker}{'values'}{'last'} = xl_col_to_name($wsu_col).($wsu_row+1);
			}
			$wsu_col++;
		}
		$wsu_row++;
	}
}
# Creates a new chart object. In this case an embedded chart.
my ($chart_col, $chart_row) = ($wsu_col+2, 3);
foreach my $step (reverse @analysis_steps_to_print){
	my $chart = $workbook->add_chart( type => 'column', embedded => 1 );
	# Ads a chart title and some axis labels.
	$chart->set_title ( name => sprintf("%s CDR3 sequences per sample", ucfirst($step)) );
	$chart->set_x_axis( name => "Samples" );
	$chart->set_y_axis( name => "Number of sequences", min => 0 );
	# Sets an Excel chart style. Colors with white outline and shadow.
	$chart->set_style( 2 );
	foreach my $marker (@$markers){
		# Configures the series.
		$chart->add_series(
			name       => ucfirst($marker),
			categories => "=unique!".$chart_series->{'all'}{'categories'}{'first'}.":".$chart_series->{'all'}{'categories'}{'last'},
			values => "=unique!".$chart_series->{$step}{$marker}{'values'}{'first'}.":".$chart_series->{$step}{$marker}{'values'}{'last'},
		);
	}
	# Inserts the chart into the worksheet (with an offset).
	$wsu->insert_chart( xl_col_to_name($chart_col).$chart_row, $chart, 0, 0, 1.5, 1.5 );
	$chart_row+=25;
}


# Prints CDR3 common sequences statistics into Excel file
if (scalar @$samples > 1 && (
   (defined($stats_cdr3_common->{$analysis_steps_to_print_[0]}{'intraindividual'}) && %{$stats_cdr3_common->{$analysis_steps_to_print_[0]}{'intraindividual'}}) || 
   (defined($stats_cdr3_common->{$analysis_steps_to_print_[0]}{'interindividual'}) && %{$stats_cdr3_common->{$analysis_steps_to_print_[0]}{'interindividual'}})
   ) ){
	print "Printing common CDR3 sequence statistics among samples.\n";
	my $wsc = $workbook->add_worksheet("common");
	my $wsc_row = 0;
	my $wsc_col = 0;
	if (%{$stats_cdr3_common->{$analysis_steps_to_print_[0]}{'intraindividual'}}){
		$wsc->write($wsc_row, $wsc_col, "Common CDR3 sequences among samples", $excel_style_bold);
		$wsc_row++;
		undef($chart_series);
		foreach my $marker (@$markers) {
			if (!defined($stats_cdr3_common->{$analysis_steps_to_print_[0]}{'intraindividual'}{$marker})){
				next;
			}
			$wsc_col = 1;
			$wsc_row++;
			$wsc->write($wsc_row, 0, "Identical $marker CDR3 sequences:", $excel_style_bold); $wsc_row+=2;
			$wsc->write($wsc_row, $wsc_col, "1 sample", $excel_style_bold); $wsc_col++;
			$chart_series->{'all'}{'categories'}{'first'} = xl_col_to_name($wsc_col-1).($wsc_row+1);
			# my $max_samples_per_ind = max(map scalar @{$dup_samples->{$_}}, keys %$dup_samples);
			my $max_samples_per_ind = max(map max(keys %{$stats_cdr3_common->{'clustered'}{'intraindividual'}{$marker}{$_}{'total'}}) , keys %{$stats_cdr3_common->{'clustered'}{'intraindividual'}{$marker}});
			if (defined($max_samples_per_ind)){
				for (my $count_samples=2; $count_samples<=$max_samples_per_ind; $count_samples++) {
					$wsc->write($wsc_row, $wsc_col, "$count_samples samples", $excel_style_bold); $wsc_col++;
				}
				$chart_series->{'all'}{'categories'}{'last'} = xl_col_to_name($wsc_col-1).($wsc_row+1);
				$wsc_row++;
				foreach my $individual (@{$individuals->{$marker}}) {
					foreach my $step (@analysis_steps_to_print_){
						$wsc_col = 0;
						$wsc->write($wsc_row, $wsc_col, "$individual $step", $excel_style_bold); $wsc_col++;
						$chart_series->{"$individual $step"}{$marker}{'values'}{'first'} = xl_col_to_name($wsc_col).($wsc_row+1);
				# 		foreach my $count_samples (sort {$a<=>$b} keys %{$stats_cdr3_common->{$step}{'intraindividual'}{$marker}{$individual}}) {
						for (my $count_samples=1; $count_samples<=$max_samples_per_ind; $count_samples++) {
							if (defined($stats_cdr3_common->{$step}{'intraindividual'}{$marker}{$individual}{'total'}{$count_samples})){
								$wsc->write($wsc_row, $wsc_col, $stats_cdr3_common->{$step}{'intraindividual'}{$marker}{$individual}{'total'}{$count_samples});$wsc_col++;
							} else {
								$wsc->write($wsc_row, $wsc_col, 0);$wsc_col++;
							}
						}
						$chart_series->{"$individual $step"}{$marker}{'values'}{'last'} = xl_col_to_name($wsc_col-1).($wsc_row+1);$wsc_col++;
						$wsc_row++;
						($chart_col, $chart_row) = ($wsc_col+2, 3);
					}
				}
			}
		}
		# Creates a new chart object. In this case an embedded chart.
		foreach my $step (@analysis_steps_to_print_){
			my $chart = $workbook->add_chart( type => 'column', embedded => 1 );
			foreach my $marker (@$markers) {
				# Ads a chart title and some axis labels.
				$chart->set_title ( name => sprintf("Number of common %s CDR3 sequences", $step) );
				$chart->set_x_axis( name => "Samples in common" );
				$chart->set_y_axis( name => "Number of sequences", min => 0 );
				# Sets an Excel chart style. Colors with white outline and shadow.
				$chart->set_style( 2 );
				# Configures the series.
				foreach my $individual (@{$individuals->{$marker}}) {
					$chart->add_series(
						name       => ucfirst($marker)."-".$individual,
						categories => "=common!".$chart_series->{'all'}{'categories'}{'first'}.":".$chart_series->{'all'}{'categories'}{'last'},
						values => "=common!".$chart_series->{"$individual $step"}{$marker}{'values'}{'first'}.":".$chart_series->{"$individual $step"}{$marker}{'values'}{'last'},
					);
				}
			}
			# Inserts the chart into the worksheet (with an offset).
			$wsc->insert_chart( xl_col_to_name($chart_col).$chart_row, $chart, 0, 0, 1.5, 1.5 );
			$chart_row+=25;
		}
		# Prints common sequences based in the number of substitutions used for clusterings
		foreach my $marker (@$markers) {
			if (!defined($stats_cdr3_common->{$analysis_steps_to_print_[0]}{'intraindividual'}{$marker})){
				next;
			}
			if (defined($INP_cluster_errors) && $INP_cluster_errors>0){
				$wsc_col = 1;
				$wsc_row+=2;
				$wsc->write($wsc_row, 0, "Similar $marker CDR3 sequences ($INP_cluster_errors subs):", $excel_style_bold); $wsc_row+=2;
				$wsc->write($wsc_row, $wsc_col, "1 sample", $excel_style_bold); $wsc_col++;
				my $max_samples_per_ind = max(map max(keys %{$stats_cdr3_common->{'clustered'}{'intraindividual'}{$marker}{$_}{'total'}}) , keys %{$stats_cdr3_common->{'clustered'}{'intraindividual'}{$marker}});
				for (my $count_samples=2; $count_samples<=$max_samples_per_ind; $count_samples++) {
					$wsc->write($wsc_row, $wsc_col, "$count_samples samples", $excel_style_bold); $wsc_col++;
				}
				$wsc_row++;
				foreach my $individual (@{$individuals->{$marker}}) {
					foreach my $step (@analysis_steps_to_print){
						$wsc_col = 0;
						$wsc->write($wsc_row, $wsc_col, "$individual $step", $excel_style_bold); $wsc_col++;
						for (my $count_samples=1; $count_samples<=$max_samples_per_ind; $count_samples++) {
							$wsc->write($wsc_row, $wsc_col, $stats_cdr3_common->{$step}{'intraindividual'}{$marker}{$individual}{'total_sim'}{$count_samples});$wsc_col++;
						}
						$wsc_row++;
					}
				}
			}
		}
		# Prints common sequences pair by pair of samples of the same individual and Chao2 estimations
		if (defined($INP_chao)){
			foreach my $marker (@$markers) {
				if (!defined($stats_cdr3_common->{$analysis_steps_to_print_[0]}{'intraindividual'}{$marker})){
					next;
				}
				$wsc_row+=2;
				$wsc->write($wsc_row, 0, "Identical $marker CDR3 sequences in pairs:", $excel_style_bold); $wsc_row+=2;
				foreach my $step (@analysis_steps_to_print_){
					$wsc_col = 1;
					$wsc->write($wsc_row, $wsc_col, "1 sample", $excel_style_bold); $wsc_col++;
					$wsc->write($wsc_row, $wsc_col, "2 samples", $excel_style_bold); $wsc_col++;
					$wsc->write($wsc_row, $wsc_col, "Chao2", $excel_style_bold); $wsc_col++;
					$wsc_row++;
					foreach my $individual (@{$individuals->{$marker}}) {
						foreach my $sample_pair (sort {$a cmp $b} keys %{$stats_cdr3_common->{$step}{'intraindividual'}{$marker}{$individual}{'pairs'}}){
							my (%chao2_cells, %chao2_values);
							$wsc_col = 0;
							my ($sample1, $sample2) = split(",", $sample_pair);
							my $common_seqs = 0;
							if (defined($stats_cdr3_common->{$step}{'intraindividual'}{$marker}{$individual}{'pairs'}{$sample_pair})){
								$common_seqs = $stats_cdr3_common->{$step}{'intraindividual'}{$marker}{$individual}{'pairs'}{$sample_pair};
							}
							my $noncommon_seqs = $stats_cdr3_global->{$sample1}{$marker}{$step}+$stats_cdr3_global->{$sample2}{$marker}{$step}-2*$common_seqs;
							$wsc->write($wsc_row, $wsc_col, sprintf("%s vs. %s %s",$sample1,$sample2,$step), $excel_style_bold); $wsc_col++;
							$wsc->write($wsc_row, $wsc_col, $noncommon_seqs); $wsc_col++;
							$wsc->write($wsc_row, $wsc_col, $common_seqs); $wsc_col++;
							# Chao2 estimator:
							my $chao2 = sprintf('%.0f', $noncommon_seqs+$common_seqs+(($noncommon_seqs^2)/(2*$common_seqs)) );
							$wsc->write_formula($wsc_row, $wsc_col, sprintf('=%s+%s+((%s^2)/(2*%s))', xl_rowcol_to_cell($wsc_row,$wsc_col-2), xl_rowcol_to_cell($wsc_row,$wsc_col-1), xl_rowcol_to_cell($wsc_row,$wsc_col-2), xl_rowcol_to_cell($wsc_row,$wsc_col-1)), undef, $chao2);
							$wsc_col++;
							$wsc_row++;
						}
					}
					$wsc_row++;
				}
			}
		}
	}
	# Prints a matrix of common clonotypes between samples (clustered and translated CDR3 sequences)
	if (defined($stats_cdr3_common->{$analysis_steps_to_print_[0]}{'interindividual'}) && %{$stats_cdr3_common->{$analysis_steps_to_print_[0]}{'interindividual'}}){
		foreach my $marker (@$markers) {
			if (!defined($stats_cdr3_common->{$analysis_steps_to_print_[0]}{'interindividual'}{$marker})){
				next;
			}
			foreach my $step (('clustered','translated')){
		# 		if (!defined($stats_cdr3_common->{$step}{'interindividual'}{$marker}) || scalar keys %{$stats_cdr3_common->{$step}{'interindividual'}{$marker}} < 2){
		# 			next;
		# 		}
				$wsc_row+=2;
				$wsc->write($wsc_row, 0, "Matrix of common $marker CDR3 $step sequences between samples:", $excel_style_bold); $wsc_row+=2;
				$wsc->write_row($wsc_row, 1, [nsort(keys %{$stats_cdr3_common->{$step}{'interindividual'}{$marker}})] , $excel_style_bold); $wsc_row++;
				foreach my $individual1 (nsort(keys %{$stats_cdr3_common->{$step}{'interindividual'}{$marker}})) {
					$wsc_col=0;
					$wsc->write($wsc_row, $wsc_col, $individual1, $excel_style_bold);
					$wsc_col++;
					foreach my $individual2 (nsort(keys %{$stats_cdr3_common->{$step}{'interindividual'}{$marker}})) {
						if (defined($stats_cdr3_common->{$step}{'interindividual'}{$marker}{$individual1}{$individual2})){
							$wsc->write($wsc_row, $wsc_col, $stats_cdr3_common->{$step}{'interindividual'}{$marker}{$individual1}{$individual2});
						} else {
							$wsc->write($wsc_row, $wsc_col, 0);
						}
						$wsc_col++;
					}
					$wsc_row++;;
				}
			}
		}
	}
}

# Prints CDR3 length statistics
print "Printing CDR3 length statistics.\n";
my $wsl = $workbook->add_worksheet("lengths");
my $wsl_row = 0;
$wsl->write($wsl_row, 0, "Distribution of lengths in CDR3 sequences", $excel_style_bold); $wsl_row+=3;
# $wsl->write_row($wsl_row, 0, [ 'Length', @$samples ] , $excel_style_bold); $wsl_row++;
undef($chart_series);
my $wsl_col = 0;
foreach my $marker (@$markers) {
	$wsl->write($wsl_row-1, $wsl_col, $marker, $excel_style_bold);
	# my $min_cdr3_len = min(map min(keys %{$stats_cdr3_lengths->{$_}{$marker}{'filtered'}}), keys %$stats_cdr3_lengths);
	# my $max_cdr3_len = max(map max(keys %{$stats_cdr3_lengths->{$_}{$marker}{'filtered'}}), keys %$stats_cdr3_lengths);
	$chart_series->{'all'}{$marker}{'categories'}{'first'} = xl_col_to_name($wsl_col).($wsl_row+2);
	$chart_series->{'all'}{$marker}{'categories'}{'last'} = xl_col_to_name($wsl_col).($wsl_row+($INP_maxlen-$INP_minlen)+2);
	$wsl->write_col($wsl_row, $wsl_col, [ "Length", $INP_minlen .. $INP_maxlen ], $excel_style_bold); $wsl_col++;
	foreach my $sample (@$samples){
		foreach my $step (@analysis_steps_to_print){
			my $wsl_row_ = $wsl_row;
			$chart_series->{"$sample $step"}{$marker}{'values'}{'first'} = xl_col_to_name($wsl_col).($wsl_row_+2);
			$wsl->write($wsl_row, $wsl_col, "$sample $step", $excel_style_bold); $wsl_row_++;
			my $stat = Statistics::Descriptive::Full->new();
			$stat->add_data(values %{$stats_cdr3_lengths->{$sample}{$marker}{$step}});
			for (my $len=$INP_minlen; $len<=$INP_maxlen; $len++){
				if (defined($stats_cdr3_lengths->{$sample}{$marker}{$step}{$len})){
					my $length_count = $stats_cdr3_lengths->{$sample}{$marker}{$step}{$len};
					my $excel_style_ = $excel_styles_red[assign_percentile_group($length_count, $stat)];
					$wsl->write($wsl_row_, $wsl_col, $stats_cdr3_lengths->{$sample}{$marker}{$step}{$len}, $excel_style_);
				} else {
					$wsl->write($wsl_row_, $wsl_col, 0);
				}
				$wsl_row_++;
			}
			$chart_series->{"$sample $step"}{$marker}{'values'}{'last'} = xl_col_to_name($wsl_col).($wsl_row_);
			$wsl_col++;
		}
	}
	$wsl_col++;
}
# Creates a new chart object. In this case an embedded chart.
($chart_col, $chart_row) = ($wsl_col+1, 4);
foreach my $step (reverse @analysis_steps_to_print){
	foreach my $marker (@$markers) {
		my $chart = $workbook->add_chart( type => 'column', embedded => 1 );
		# Ads a chart title and some axis labels.
		$chart->set_title ( name => "$marker CDR3 $step sequences length histogram" );
		# $chart->set_x_axis( name => "Lengths" );
		# $chart->set_y_axis( name => "Number of sequences", min => 0 );
		# Sets an Excel chart style. Colors with white outline and shadow.
		$chart->set_style( 2 );
		# Configures the series.
		foreach my $sample (@$samples){
			$chart->add_series(
				name       => $sample,
				categories => "=lengths!".$chart_series->{'all'}{$marker}{'categories'}{'first'}.":".$chart_series->{'all'}{$marker}{'categories'}{'last'},
				values => "=lengths!".$chart_series->{"$sample $step"}{$marker}{'values'}{'first'}.":".$chart_series->{"$sample $step"}{$marker}{'values'}{'last'},
			);
		}
		# Inserts the chart into the worksheet (with an offset).
		$wsl->insert_chart( xl_col_to_name($chart_col).$chart_row, $chart, 0, 0, 1.5, 1.5 );
		$chart_row+=25;
	}
}


# Prints CDR3 depth statistics
print "Printing CDR3 depth statistics.\n";
my $wsd = $workbook->add_worksheet("depths");
my $wsd_row = 0;
$wsd->write($wsd_row, 0, "Distribution of depths in clustered CDR3 sequences", $excel_style_bold); $wsd_row+=3;
# my @depths_to_print = sort {$a<=>$b} unique(map keys %{$stats_cdr3_depths->{$_}{'clustered'}}, keys %$stats_cdr3_depths);
# my $min_cdr3_depth = min(map min(keys %{$stats_cdr3_depths->{$_}{'clustered'}}), keys %$stats_cdr3_depths);
# my $max_cdr3_depth = max(map max(keys %{$stats_cdr3_depths->{$_}{'clustered'}}), keys %$stats_cdr3_depths);
my @depths_to_print;
foreach my $sample (keys %$stats_cdr3_depths){
	push(@depths_to_print, map keys %{$stats_cdr3_depths->{$sample}{$_}{'clustered'}}, keys %{$stats_cdr3_depths->{$sample}});
}
@depths_to_print = sort {$a<=>$b} unique(@depths_to_print);
undef($chart_series);
my $wsd_col = 0;
foreach my $marker (@$markers) {
	$wsl->write($wsd_row-1, $wsd_col, $marker, $excel_style_bold);
	# $wsd->write_row($wsd_row, 0, [ 'Length', @$samples ] , $excel_style_bold); $wsd_row++;
	$chart_series->{'all'}{$marker}{'categories'}{'first'} = xl_col_to_name($wsd_col).($wsd_row+2);
	$chart_series->{'all'}{$marker}{'categories'}{'last'} = xl_col_to_name($wsd_col).($wsd_row+(scalar @depths_to_print)+2);
	$wsd->write_col($wsd_row, $wsd_col, [ 'Depth', @depths_to_print ], $excel_style_bold); $wsd_col++;
	foreach my $sample (@$samples){
		foreach my $step (('clustered')){
			my $wsd_row_ = $wsd_row;
			$chart_series->{"$sample $step"}{$marker}{'values'}{'first'} = xl_col_to_name($wsd_col).($wsd_row_+2);
			$wsd->write($wsd_row, $wsd_col, "$sample $step", $excel_style_bold); $wsd_row_++;
			my $stat = Statistics::Descriptive::Full->new();
			$stat->add_data(values %{$stats_cdr3_depths->{$sample}{$marker}{$step}});
			foreach my $depth (@depths_to_print) {
				if (defined($stats_cdr3_depths->{$sample}{$marker}{$step}{$depth})){
					my $depth_count = $stats_cdr3_depths->{$sample}{$marker}{$step}{$depth};
					my $excel_style_ = $excel_styles_red[assign_percentile_group($depth_count, $stat)];
					$wsd->write($wsd_row_, $wsd_col, $stats_cdr3_depths->{$sample}{$marker}{$step}{$depth}, $excel_style_);
				} else {
					$wsd->write($wsd_row_, $wsd_col, 0);
				}
				$wsd_row_++;
			}
			$chart_series->{"$sample $step"}{$marker}{'values'}{'last'} = xl_col_to_name($wsd_col).($wsd_row_);
			$wsd_col++;
		}
	}
	$wsd_col++;
}
# Creates a new chart object. In this case an embedded chart.
($chart_col, $chart_row) = ($wsd_col+1, 4);
foreach my $step (('clustered')){
	foreach my $marker (@$markers) {
		my $chart = $workbook->add_chart( type => 'column', embedded => 1 );
		# Ads a chart title and some axis labels.
		$chart->set_title ( name => "$marker CDR3 $step sequences depth histogram" );
		# $chart->set_x_axis( name => "Lengths" );
		# $chart->set_y_axis( name => "Number of sequences", min => 0 );
		# Sets an Excel chart style. Colors with white outline and shadow.
		$chart->set_style( 2 );
		# Configures the series.
		foreach my $sample (@$samples){
			$chart->add_series(
				name       => $sample,
				categories => "=depths!".$chart_series->{'all'}{$marker}{'categories'}{'first'}.":".$chart_series->{'all'}{$marker}{'categories'}{'last'},
				values => "=depths!".$chart_series->{"$sample $step"}{$marker}{'values'}{'first'}.":".$chart_series->{"$sample $step"}{$marker}{'values'}{'last'},
			);
		}
		# Inserts the chart into the worksheet (with an offset).
		$wsd->insert_chart( xl_col_to_name($chart_col).$chart_row, $chart, 0, 0, 1.5, 1.5 );
		$chart_row+=25;
	}
}

# Prints TCR regions statistics
if (%INP_ref_files) {

	print "Printing CDR3 associated TCR regions statistics.\n";
	my $wsr = $workbook->add_worksheet("regions");
	my $wsr_row = 0;
	my $wsr_col = 0;
	$wsr->write($wsr_row, $wsr_col, "Distribution of V and J genes per sample", $excel_style_bold); $wsr_row+=2;

	$chart_row = $wsr_row+1;
	$chart_col = 0;
	undef($chart_series);
	foreach my $marker (@$markers) {
		foreach my $step (@analysis_steps_to_print){
			$wsr_col = 1;
			$chart_series->{'all'}{'categories'}{'first'} = xl_col_to_name($wsr_col).($wsr_row+1);
			my $tcr_alleles;
			$tcr_alleles->{'tcrj'} = [ nsort(unique( map keys %{$stats_tcr_regions->{$_}{$marker}{$step}{'tcrj'}}, @$samples )) ];
			$tcr_alleles->{'tcrv'} = [ nsort(unique( map keys %{$stats_tcr_regions->{$_}{$marker}{$step}{'tcrv'}}, @$samples )) ];
			my $chart_col_ = scalar @{$tcr_alleles->{'tcrj'}} + scalar @{$tcr_alleles->{'tcrv'}};
			if ($chart_col_ > $chart_col){
				$chart_col = $chart_col_;
			}
			foreach my $region ( ('tcrj', 'tcrv') ){
				foreach my $tcr_allele (@{$tcr_alleles->{$region}}){
					$wsr->write($wsr_row, $wsr_col, uc($tcr_allele), $excel_style_bold);
					$wsr_col++;
				}
			}
			$chart_series->{'all'}{'categories'}{'last'} = xl_col_to_name($wsr_col-1).($wsr_row+1);
			$wsr_row++;
			foreach my $sample (@$samples){
				my $wsr_col = 0;
				$wsr->write($wsr_row, $wsr_col, "$sample $step", $excel_style_bold); $wsr_col++;
				$chart_series->{"$sample $step"}{'values'}{'first'} = xl_col_to_name($wsr_col).($wsr_row+1);
				foreach my $region ( ('tcrj', 'tcrv') ){
					my $stat = Statistics::Descriptive::Full->new();
					$stat->add_data(values %{$stats_tcr_regions->{$sample}{$marker}{$step}{'tcrv'}});
					foreach my $tcr_allele (@{$tcr_alleles->{$region}}){
						if (defined($stats_tcr_regions->{$sample}{$marker}{$step}{$region}{$tcr_allele})) {
							my $allele_count = $stats_tcr_regions->{$sample}{$marker}{$step}{$region}{$tcr_allele};
							my $excel_style_;
							if ($region eq 'tcrj'){
								$excel_style_ = $excel_styles_red[assign_percentile_group($allele_count, $stat)];
							} elsif ($region eq 'tcrv'){
								$excel_style_ = $excel_styles_green[assign_percentile_group($allele_count, $stat)];
							}
							$wsr->write($wsr_row, $wsr_col, $allele_count, $excel_style_);
						} else {
							$wsr->write($wsr_row, $wsr_col, 0);
						}
						$wsr_col++;
					}
				}
				$chart_series->{"$sample $step"}{'values'}{'last'} = xl_col_to_name($wsr_col-1).($wsr_row+1);
				$wsr_row++;
			}
			$wsr_row++;
		}
		# Creates a new chart object. In this case an embedded chart.
		foreach my $step (reverse @analysis_steps_to_print){
			my $chart = $workbook->add_chart( type => 'column', embedded => 1 );
			# Ads a chart title and some axis labels.
			$chart->set_title ( name => "V and J genes distribution in $step CDR3s") ;
			#$chart->set_x_axis( name => "TCR alleles" );
			$chart->set_y_axis( name => "Number of sequences", min => 0 );
			# Sets an Excel chart style. Colors with white outline and shadow.
			$chart->set_style( 2 );
			# Configures the series.
			foreach my $sample (@$samples){
				$chart->add_series(
					name       => $sample,
					categories => "=regions!".$chart_series->{'all'}{'categories'}{'first'}.":".$chart_series->{'all'}{'categories'}{'last'},
					values => "=regions!".$chart_series->{"$sample $step"}{'values'}{'first'}.":".$chart_series->{"$sample $step"}{'values'}{'last'},
				);
			}
			# Inserts the chart into the worksheet (with an offset).
			$chart->set_size( width => 1000, height => 288 );
			$wsr->insert_chart( xl_col_to_name($chart_col+2).$chart_row, $chart, 0, 0, 1.5, 1.5 );
			$chart_row += 25;
		}
	}

	$wsr = $workbook->add_worksheet("regions2");
	$wsr_row = 0;
	$wsr->write($wsr_row, 0, "V-J gene usage per sample", $excel_style_bold); $wsr_row+=2;

	foreach my $step (@analysis_steps_to_print){

		$chart_col = 0;
		$chart_row = 3;
		# Prints all samples together TCRJ vs. TCRV. usage matrix
		foreach my $marker (@$markers) {
			$wsr_col = 0;
			$wsr->write($wsr_row, $wsr_col, "$marker $step", $excel_style_bold); $wsr_row++;
			my $stat = Statistics::Descriptive::Full->new();
			$wsr_col = 1;
			my $tcr_alleles;
			$tcr_alleles->{'tcrj'} = [ nsort(unique( map keys %{$stats_tcr_regions->{$_}{$marker}{$step}{'tcrj'}}, @$samples )) ];
			$tcr_alleles->{'tcrv'} = [ nsort(unique( map keys %{$stats_tcr_regions->{$_}{$marker}{$step}{'tcrv'}}, @$samples )) ];
			my $chart_col_ = scalar @{$tcr_alleles->{'tcrv'}};
			if ($chart_col_ > $chart_col){
				$chart_col = $chart_col_;
			}
			undef($chart_series);
			$chart_series->{'all'}{'categories'}{'first'} = xl_col_to_name($wsr_col).($wsr_row+1);
			foreach my $tcrv_allele (@{$tcr_alleles->{'tcrv'}}){
				$wsr->write($wsr_row, $wsr_col, uc($tcrv_allele), $excel_style_bold);
				$wsr_col++;
				foreach my $sample (@$samples){
					$stat->add_data(values %{$stats_tcr_regions->{$sample}{$marker}{$step}{'pairs'}{$tcrv_allele}});
				}
			}
			$chart_series->{'all'}{'categories'}{'last'} = xl_col_to_name($wsr_col-1).($wsr_row+1);
			$wsr_row++;
			foreach my $tcrj_allele (@{$tcr_alleles->{'tcrj'}}){
				$wsr_col = 0;
				$wsr->write($wsr_row, $wsr_col, uc($tcrj_allele), $excel_style_bold);$wsr_col++;
				$chart_series->{$tcrj_allele}{'values'}{'first'} = xl_col_to_name($wsr_col).($wsr_row+1);
				foreach my $tcrv_allele (@{$tcr_alleles->{'tcrv'}}){
					my $allele_pair_count = 0;
					foreach my $sample (@$samples){
						if (defined($stats_tcr_regions->{$sample}{$marker}{$step}{'pairs'}{$tcrv_allele}{$tcrj_allele})){
							$allele_pair_count += $stats_tcr_regions->{$sample}{$marker}{$step}{'pairs'}{$tcrv_allele}{$tcrj_allele};
						}
					}
					my $excel_style_ = $excel_styles_green[assign_percentile_group($allele_pair_count, $stat)];
					$wsr->write($wsr_row, $wsr_col, $allele_pair_count, $excel_style_);
					$wsr_col++;
				}
				$chart_series->{$tcrj_allele}{'values'}{'last'} = xl_col_to_name($wsr_col-1).($wsr_row+1);
				$wsr_row++;
			}
			$wsr_row++;

			# Creates a new chart object. In this case an embedded chart.
			my $chart = $workbook->add_chart( type => 'column', embedded => 1 );
			# Ads a chart title and some axis labels.
			$chart->set_title ( name => "$marker V-J gene usage in $step CDR3s") ;
			#$chart->set_x_axis( name => "TCR alleles" );
			$chart->set_y_axis( name => "Number of sequences", min => 0 );
			# Sets an Excel chart style. Colors with white outline and shadow.
			$chart->set_style( 2 );
			# Configures the series.
			foreach my $tcrj_allele (@{$tcr_alleles->{'tcrj'}}){
				$chart->add_series(
					name       => $tcrj_allele,
					categories => "=regions2!".$chart_series->{'all'}{'categories'}{'first'}.":".$chart_series->{'all'}{'categories'}{'last'},
					values => "=regions2!".$chart_series->{$tcrj_allele}{'values'}{'first'}.":".$chart_series->{$tcrj_allele}{'values'}{'last'},
				);
			}
			# Inserts the chart into the worksheet (with an offset).
			$chart->set_size( width => 1000, height => 288 );
			$wsr->insert_chart( xl_col_to_name($chart_col+2).$chart_row, $chart, 0, 0, 1.5, 1.5 );
			$chart_row += 25

		}
		$wsr_row++;

		# Prints sample by sample TCRJ vs. TCRV. usage matrix
		foreach my $sample (@$samples){
			foreach my $marker (@$markers) {
				$wsr_col = 0;
				$wsr->write($wsr_row, $wsr_col, "$marker $sample $step", $excel_style_bold); $wsr_row++;
				my $stat = Statistics::Descriptive::Full->new();
				$wsr_col = 1;
				my $tcr_alleles;
				$tcr_alleles->{'tcrj'} = [ nsort(unique( map keys %{$stats_tcr_regions->{$_}{$marker}{$step}{'tcrj'}}, @$samples )) ];
				$tcr_alleles->{'tcrv'} = [ nsort(unique( map keys %{$stats_tcr_regions->{$_}{$marker}{$step}{'tcrv'}}, @$samples )) ];
				foreach my $tcrv_allele (@{$tcr_alleles->{'tcrv'}}){
					$wsr->write($wsr_row, $wsr_col, uc($tcrv_allele), $excel_style_bold);
					$wsr_col++;
					$stat->add_data(values %{$stats_tcr_regions->{$sample}{$marker}{$step}{'pairs'}{$tcrv_allele}});
				}
				$wsr_row++;
				foreach my $tcrj_allele (@{$tcr_alleles->{'tcrj'}}){
					$wsr_col = 0;
					$wsr->write($wsr_row, $wsr_col, uc($tcrj_allele), $excel_style_bold);$wsr_col++;
					foreach my $tcrv_allele (@{$tcr_alleles->{'tcrv'}}){
						if (defined($stats_tcr_regions->{$sample}{$marker}{$step}{'pairs'}{$tcrv_allele}{$tcrj_allele})){
							my $allele_pair_count = $stats_tcr_regions->{$sample}{$marker}{$step}{'pairs'}{$tcrv_allele}{$tcrj_allele};
							my $excel_style_ = $excel_styles_green[assign_percentile_group($allele_pair_count, $stat)];
							$wsr->write($wsr_row, $wsr_col, $allele_pair_count, $excel_style_);
						} else {
							$wsr->write($wsr_row, $wsr_col, 0);
						}
						$wsr_col++;
					}
					$wsr_row++;
				}
				$wsr_row++;
			}
			$wsr_row++;
		}
		$wsr_row++;
	}
# 	foreach my $step (@analysis_steps_to_print){
# 		foreach my $sample (@$samples){
# 			$wsr_col = 0;
# 			$wsr->write($wsr_row, $wsr_col, "$sample $step", $excel_style_bold); $wsr_row++;
# 			my $stat = Statistics::Descriptive::Full->new();
# 			$wsr_col = 1;
# 			my $tcr_alleles;
# 			foreach my $marker (@$markers) {
# 				$tcr_alleles->{$marker}{'tcrj'} = [ nsort(unique( map keys %{$stats_tcr_regions->{$_}{$marker}{$step}{'tcrj'}}, @$samples )) ];
# 				$tcr_alleles->{$marker}{'tcrv'} = [ nsort(unique( map keys %{$stats_tcr_regions->{$_}{$marker}{$step}{'tcrv'}}, @$samples )) ];
# 				foreach my $tcrv_allele (@{$tcr_alleles->{$marker}{'tcrv'}}){
# 					$wsr->write($wsr_row, $wsr_col, uc($tcrv_allele), $excel_style_bold);
# 					$wsr_col++;
# 					$stat->add_data(values %{$stats_tcr_regions->{$sample}{$marker}{$step}{'pairs'}{$tcrv_allele}});
# 				}
# 			}
# 			$wsr_row++;
# 			foreach my $marker1 (@$markers) {
# 				foreach my $tcrj_allele (@{$tcr_alleles->{$marker1}{'tcrj'}}){
# 					$wsr_col = 0;
# 					$wsr->write($wsr_row, $wsr_col, uc($tcrj_allele), $excel_style_bold);$wsr_col++;
# 					foreach my $marker2 (@$markers) {
# 						foreach my $tcrv_allele (@{$tcr_alleles->{$marker2}{'tcrv'}}){
# 							if (defined($stats_tcr_regions->{$sample}{$marker1}{$step}{'pairs'}{$tcrv_allele}{$tcrj_allele})){
# 								my $allele_pair_count = $stats_tcr_regions->{$sample}{$marker1}{$step}{'pairs'}{$tcrv_allele}{$tcrj_allele};
# 								my $excel_style_ = $excel_styles_green[assign_percentile_group($allele_pair_count, $stat)];
# 								$wsr->write($wsr_row, $wsr_col, $allele_pair_count, $excel_style_);
# 							} else {
# 								$wsr->write($wsr_row, $wsr_col, 0);
# 							}
# 							$wsr_col++;
# 						}
# 					}
# 				}
# 				$wsr_row++;
# 			}
# 			$wsr_row++;
# 		}
# 		$wsr_row++;
# 	}

	
	$wsr = $workbook->add_worksheet("regions3");
	$wsr_row = 0;
	$wsr->write($wsr_row, 0, "CDR3 lengths vs. TCRJ and TCRV genes", $excel_style_bold); $wsr_row+=2;
	foreach my $marker (@$markers) {
		foreach my $step (@analysis_steps_to_print){
			my $tcr_alleles;
			$tcr_alleles->{'tcrj'} = [ nsort(unique( map keys %{$stats_tcr_regions->{$_}{$marker}{$step}{'tcrj'}}, @$samples )) ];
			$tcr_alleles->{'tcrv'} = [ nsort(unique( map keys %{$stats_tcr_regions->{$_}{$marker}{$step}{'tcrv'}}, @$samples )) ];
			foreach my $sample (@$samples){
				$wsr_col = 0;
				$wsr->write($wsr_row, $wsr_col, "$marker $sample $step", $excel_style_bold); $wsr_row++;
				$wsr->write($wsr_row, $wsr_col, "Length", $excel_style_bold); $wsr_col++;
				my $stat_tcrj = Statistics::Descriptive::Full->new();
				my $stat_tcrv = Statistics::Descriptive::Full->new();
				my %cdr3_lengths;
				foreach my $tcrj_allele (@{$tcr_alleles->{'tcrj'}}){
					$wsr->write($wsr_row, $wsr_col, uc($tcrj_allele), $excel_style_bold);
					$wsr_col++;
					$stat_tcrj->add_data(values %{$stats_tcr_regions->{$sample}{$marker}{$step}{'lengths'}{'tcrj'}{$tcrj_allele}});
					foreach my $len (keys %{$stats_tcr_regions->{$sample}{$marker}{$step}{'lengths'}{'tcrj'}{$tcrj_allele}}){
						$cdr3_lengths{$len} = 1;
					}
				}
				foreach my $tcrv_allele (@{$tcr_alleles->{'tcrv'}}){
					$wsr->write($wsr_row, $wsr_col, uc($tcrv_allele), $excel_style_bold);
					$wsr_col++;
					$stat_tcrv->add_data(values %{$stats_tcr_regions->{$sample}{$marker}{$step}{'lengths'}{'tcrv'}{$tcrv_allele}});
					foreach my $len (keys %{$stats_tcr_regions->{$sample}{$marker}{$step}{'lengths'}{'tcrv'}{$tcrv_allele}}){
						$cdr3_lengths{$len} = 1;
					}
				}
				$wsr_row++;
				foreach my $len (nsort(keys %cdr3_lengths)){
					$wsr_col = 0;
					$wsr->write($wsr_row, $wsr_col, $len, $excel_style_bold);$wsr_col++;
					foreach my $tcrj_allele (@{$tcr_alleles->{'tcrj'}}){
						if (defined($stats_tcr_regions->{$sample}{$marker}{$step}{'lengths'}{'tcrj'}{$tcrj_allele}{$len})){
							my $cdr3_length_count = $stats_tcr_regions->{$sample}{$marker}{$step}{'lengths'}{'tcrj'}{$tcrj_allele}{$len};
							my $excel_style_ = $excel_styles_red[assign_percentile_group($cdr3_length_count, $stat_tcrj)];
							$wsr->write($wsr_row, $wsr_col, $cdr3_length_count, $excel_style_);
						} else {
							$wsr->write($wsr_row, $wsr_col, 0);
						}
						$wsr_col++;
					}
					foreach my $tcrv_allele (@{$tcr_alleles->{'tcrv'}}){
						if (defined($stats_tcr_regions->{$sample}{$marker}{$step}{'lengths'}{'tcrv'}{$tcrv_allele}{$len})){
							my $cdr3_length_count = $stats_tcr_regions->{$sample}{$marker}{$step}{'lengths'}{'tcrv'}{$tcrv_allele}{$len};
							my $excel_style_ = $excel_styles_green[assign_percentile_group($cdr3_length_count, $stat_tcrv)];
							$wsr->write($wsr_row, $wsr_col, $cdr3_length_count, $excel_style_);
						} else {
							$wsr->write($wsr_row, $wsr_col, 0);
						}
						$wsr_col++;
					}
					$wsr_row++;
				}
				$wsr_row++;
			}
		}
	}
}


# Prints major CDR3 clonotypes into Excel file
print "\nPrinting top CDR3 clonotype data.\n";

# Checks how many clonotypes to print into a single sheet of one sheet per sample
my $count_clones_per_sample = 0;
map $count_clones_per_sample += scalar @{$results->{$_}}, map @$samples;
$count_clones_per_sample = $count_clones_per_sample / scalar @$samples;
my @fields;
if (!defined($INP_umi_cluster)) {
	@fields = ('cloneId', 'cloneCount', 'cloneFraction', 'nSeqCDR3', 'aaSeqCDR3', 'allVHitsWithScore', 'allJHitsWithScore');
} else {
	@fields = ('cloneId', 'cloneCount', 'cloneFraction', 'totalReads', 'readFraction', 'nSeqCDR3', 'aaSeqCDR3', 'allVHitsWithScore', 'allJHitsWithScore');
}
my $wsf;
my $wsf_row = 0;
my $wsf_col = 0;
my $samples_per_sheet;
if (scalar @$samples <= 10){
	$samples_per_sheet = 1;
} elsif (scalar @$samples <= 20){
	$samples_per_sheet = 2;
} elsif (scalar @$samples <= 50){
	$samples_per_sheet = 5;
} else {
	$samples_per_sheet = 10;
}
my ($count_samples, $count_sheets) = (0,0);
foreach my $sample (@$samples){
	if ($count_samples == 0 || $count_samples == $samples_per_sheet){
		$count_samples = 0;
		$count_sheets++;
		if ($samples_per_sheet == 1){
			$wsf = $workbook->add_worksheet(substr("chains_$sample",0,31));
		} else {
			$wsf = $workbook->add_worksheet(sprintf("chains%d",$count_sheets));
		}
		$wsf_col = 0;
	}
	$count_samples++;
	my $top_clones = scalar @{$results->{$sample}};
	$wsf_row = 0;
	if (scalar @{$results->{$sample}} <= 10000){
		$wsf->write($wsf_row, $wsf_col, "$top_clones CDR3 clonotypes found in sample '$sample':", $excel_style_bold); $wsf_row+=2;
	} else{
		$top_clones = 10000;
		$wsf->write($wsf_row, $wsf_col, "Top $top_clones CDR3 clonotypes in sample '$sample':", $excel_style_bold); $wsf_row+=2;
	}
	$wsf->write_row($wsf_row, $wsf_col, \@fields , $excel_style_bold); $wsf_row++;
	my $count_results = 0;
	foreach my $result (sort {$b->{'cloneCount'}<=>$a->{'cloneCount'}} @{$results->{$sample}}){
		$count_results++;
		if ($count_results > $top_clones){
			last;
		}
		my @values = map $result->{$_} , @fields;
		$wsf->write_row($wsf_row, $wsf_col, \@values); $wsf_row++;
	}
	$wsf_col += 2 + scalar @fields;
}


# if ($count_clones_per_sample <= 10) {
# 	$wsf = $workbook->add_worksheet("clones");
# 	$wsf->write_row($wsf_row, 0, ['sample', @fields] , $excel_style_bold); $wsf_row++;
# }
# foreach my $sample (@$samples){
# 	my $top_clones = scalar @{$results->{$sample}};
# 	# When there are more than 10 clonotypes per sample, writes the results in separate spreadsheets
# 	if ($count_clones_per_sample > 10) {
# 		$wsf = $workbook->add_worksheet(substr("chains_$sample",0,31));
# 		$wsf_row = 0;
# 		# Writes worksheet legend
# 		if (scalar @{$results->{$sample}} <= 10000){
# 			$wsf->write($wsf_row, 0, "$top_clones CDR3 clonotypes found in sample '$sample':", $excel_style_bold); $wsf_row+=2;
# 		} else{
# 			$top_clones = 10000;
# 			$wsf->write($wsf_row, 0, "Top $top_clones CDR3 clonotypes in sample '$sample':", $excel_style_bold); $wsf_row+=2;
# 		}
# 		$wsf->write_row($wsf_row, 0, \@fields , $excel_style_bold); $wsf_row++;
# 	}
# 	my $count_results = 0;
# 	foreach my $result (sort {$b->{'cloneCount'}<=>$a->{'cloneCount'}} @{$results->{$sample}}){
# 		$count_results++;
# 		if ($count_results > $top_clones){
# 			last;
# 		}
# 		my @values = map $result->{$_} , @fields;
# 		if ($count_clones_per_sample > 10) {
# 			$wsf->write_row($wsf_row, 0, \@values); $wsf_row++;
# 		} else {
# 			$wsf->write_row($wsf_row, 0, [$sample, @values]); $wsf_row++;
# 		}
# 	}
# }


# Closes Excel file
$workbook->close();


printf("\n'%s' Excel file created with CDR3 analysis statistics.\n", $excel_outputfile);
print "\n";

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


exit;


#################################################################################

sub count_umis {

	my $seq_headers = shift @_;
	
	my %umis;
	foreach (@$seq_headers){
		my ($umi,$count) = ('',0);
		if (/umi=([ACGT]+)/){
			$umi = $1;
		}
		if (/count=(\d+)/){
			$count = $1;
		}
		$umis{$umi} += $count;
	}
	
	return %umis;
	
}

#################################################################################

sub extract_cdr3_seqs {

	my ($read_files,$markerdata,$sampledata,$total_reads,$tcr_ref_patterns,$options) = @_; # $tcr_ref_files,$options) = @_;

	my ($umi,$verbose)=(0,1);
	# Reads and annotates Unique Molecular Identifiers (UMIs) into the sequence headers
	if (in_array($options, 'umi')){
		$umi = 1;
	}
	# Does not print any additional analysis information
	if (in_array($options, 'quiet')){
		$verbose = 0;
	}

	# Create search patterns in reads
	my ($cdr3_patterns, $tcr_patterns, $umi_data);
	my ($id_warning,$dist_warning)=(0,0);
	foreach my $sample (keys %$sampledata){
		foreach my $marker (keys %$markerdata){
			# Recognizes the chain defined by the marker
			my $chain;
			if ($marker =~ /TA|TRA|TCRA|alpha/i){
				$chain = 'alpha';
			} elsif ($marker =~ /TB|TRB|TCRB|beta/i){
				$chain = 'beta';
			} elsif ($marker =~ /TG|TRG|TCRG|gamma/i){
				$chain = 'gamma';
			} elsif ($marker =~ /TD|TRD|TCRD|delta/i){
				$chain = 'delta';
			} else {
				$chain = 'unknown';
			}
			# Assigns the CDR3 pattern to use
			my $cdr3_chain_pattern;
			if (defined($tcr_ref_patterns->{$chain}) && defined($tcr_ref_patterns->{'custom'})){
				$cdr3_chain_pattern = $tcr_ref_patterns->{$chain}.'|'.$tcr_ref_patterns->{'custom'};
			} elsif (defined($tcr_ref_patterns->{$chain})){
				$cdr3_chain_pattern = $tcr_ref_patterns->{$chain};
			} else {
				$cdr3_chain_pattern = $tcr_ref_patterns->{'custom'};
			}
			# In case the pattern is not defined, ex. a species that has no pattern for a specific chain
			if (!defined($cdr3_chain_pattern)){
				next;
			}
# 			if (defined($tcr_ref_patterns->{$chain})){
# 				$cdr3_chain_pattern = $tcr_ref_patterns->{$chain};
# 			} else {
# 				$cdr3_chain_pattern = $tcr_ref_patterns->{'custom'};
# 			}
			my $max_i = defined($markerdata->{$marker}{'primer_f'}) ? $#{$markerdata->{$marker}{'primer_f'}} : 0;
			for (my $i=0; $i<=$max_i; $i++){
# 				if (!@{$markerdata->{$marker}{'primer_r_dist'}} || !@{$markerdata->{$marker}{'primer_r_id'}} ){
# 					print "\nERROR: You must define 'primer_r_dist' and 'primer_r_id' values into '$INP_amplicons_file'.\n\n";
# 					exit;
# 				}
	# 			if (!@{$markerdata->{$marker}{'primer_f_dist'}} || !@{$markerdata->{$marker}{'primer_f_id'}} ||
	# 			!@{$markerdata->{$marker}{'primer_r_dist'}} || !@{$markerdata->{$marker}{'primer_r_id'}} ){
	# 				print "\nERROR: You must define 'primer_f_dist', 'primer_f_id', 'primer_r_dist' and 'primer_r_id' values into '$INP_amplicons_file'.\n\n";
	# 				exit;
	# 			}
				my $max_j = defined($markerdata->{$marker}{'primer_r'}) ? $#{$markerdata->{$marker}{'primer_r'}} : 0;
				for (my $j=0; $j<=$max_j; $j++){
					# Tags and primers must be redefined in every iteration to remove later the non standard symbols
					my $fwd_tag = defined($sampledata->{$sample}{'tag_f'}) ? $sampledata->{$sample}{'tag_f'} : '';
					my $rev_tag = defined($sampledata->{$sample}{'tag_r'}) ? $sampledata->{$sample}{'tag_r'} : '';
					my $fwd_primer = defined($markerdata->{$marker}{'primer_f'}[$i]) ? $markerdata->{$marker}{'primer_f'}[$i] : '';
					my $fwd_primer_id = defined($markerdata->{$marker}{'primer_f_id'}[$i]) ? $markerdata->{$marker}{'primer_f_id'}[$i] : '';
					my $rev_primer = defined($markerdata->{$marker}{'primer_r'}[$j]) ? $markerdata->{$marker}{'primer_r'}[$j] : '';
					my $rev_primer_dist;
					if (defined($markerdata->{$marker}{'primer_r_dist'})){
						$rev_primer_dist = $markerdata->{$marker}{'primer_r_dist'}[$j];
					}
					my $rev_primer_id = defined($markerdata->{$marker}{'primer_r_id'}[$j]) ? $markerdata->{$marker}{'primer_r_id'}[$j] : '';
					if (!$dist_warning && (!defined($fwd_primer_id) || !defined($rev_primer_id))){
						print "\tWARNING: some primers have no identifiers defined, primer sequence will be used as identifier.\n";
						$dist_warning = 1;
					}
					my $amplicon_id;
					if ($fwd_primer_id ne '' && $rev_primer_id ne ''){
						$amplicon_id = $marker."-".$fwd_primer_id."-".$rev_primer_id;
					} elsif ($fwd_primer_id ne ''){
						$amplicon_id = $marker."-".$fwd_primer_id;
					} elsif ($rev_primer_id ne ''){
						$amplicon_id = $marker."-".$rev_primer_id;
					} else {
						$amplicon_id= $marker;
					}
					# Reads UMI info
					# MODIFICATION: UMI CAN BE AT THE MIDDLE OF THE PATTERN (AFTER PRIMER)
					my %umi_data_;
					if ($umi){
# 						my $umi_pattern = '([^N]?)(N{5,})';
						my $umi_pattern = '\(([ACGTUN]+)\)';
						if ($fwd_tag =~ /$umi_pattern/){
							$umi_data_{'len'} = length($1);
							$umi_data_{'pos'} = $-[0]+1;
						} elsif ($rev_tag =~ /$umi_pattern/){
							$umi_data_{'len'} = -(length($1));
							$umi_data_{'pos'} = length($rev_primer)+$-[0]+1;
						} elsif ($fwd_primer =~ /$umi_pattern/){
							$umi_data_{'len'} = length($1);
							$umi_data_{'pos'} = length($fwd_tag)+$-[0]+1;
						} elsif ($rev_primer =~ /$umi_pattern/){
							$umi_data_{'len'} = -(length($1));
							$umi_data_{'pos'} = $-[0]+1;
						} else {
							print "\nERROR: You must include Unique Molecular Identifiers between parenthesis into tag or primer sequences at '$INP_amplicons_file'.\n\n";
							exit;
						}
					}
					# Remove any parenthesis or not standard symbol from sequences
					$fwd_tag =~ s/[^\^ACGTUN]//g;
					$rev_tag =~ s/[^\$ACGTUN]//g;
					$fwd_primer =~ s/[^ACGTUN]//g;
					$rev_primer =~ s/[^ACGTUN]//g;
					# With the TCR pattern we will be able to identify samples
					if ($fwd_tag.$fwd_primer ne '' && $rev_tag.$rev_primer ne ''){
						my $tcr_pattern = sprintf('%s(\w*)%s', regex($fwd_tag.$fwd_primer), regex(iupac_reverse_complementary($rev_primer).iupac_reverse_complementary($rev_tag)) );
						$tcr_patterns->{$sample}{$marker}{$amplicon_id} = $tcr_pattern;
						# $tcr_patterns->{$sample}{$marker}{$amplicon_id} = sprintf('%s%s(\w+)%s%s', regex($fwd_tag), regex($fwd_primer), regex(iupac_reverse_complementary($rev_primer)), regex(iupac_reverse_complementary($rev_tag)) );
					} else {
						$tcr_patterns->{$sample}{$marker}{$amplicon_id} = undef;
					}
					# IMPORTANT MODIFICATION: UMI (tag_f) IS AFTER FORWARD PRIMER (primer_f)
					# The motif TGYRSN is in the border of V region and CDR3 (TGY is in V region and RSN is inside CDR3)
					my $cdr3_pattern = $cdr3_chain_pattern;
					if (defined($rev_primer_dist)){
						$cdr3_pattern =~ s/\).+/)/;
						$cdr3_pattern = sprintf('%s[\w-]+%s\w{%d}%s', regex($fwd_tag.$fwd_primer), $cdr3_pattern, $rev_primer_dist, regex(iupac_reverse_complementary($rev_primer).iupac_reverse_complementary($rev_tag)) );
					} else {
						$cdr3_pattern = sprintf('%s[\w-]+%s\w*%s', regex($fwd_tag.$fwd_primer), $cdr3_pattern, regex(iupac_reverse_complementary($rev_primer).iupac_reverse_complementary($rev_tag)) );
					}
					$cdr3_patterns->{$sample}{$marker}{$amplicon_id} = $cdr3_pattern;
	# 				$cdr3_patterns->{$sample}{$marker}{$amplicon_id} = sprintf('%s%s\w{%d}(\w+)\w{%d}%s%s', regex($fwd_tag), regex($fwd_primer), $fwd_primer_dist-$INP_extra_length, $rev_primer_dist-$INP_extra_length, regex(iupac_reverse_complementary($rev_primer)), regex(iupac_reverse_complementary($rev_tag)) );
					# UMI information: position and length
					if ($umi){
						$umi_data->{$sample}{$marker}{$cdr3_pattern} = \%umi_data_;
					}
				}
			}
		}
	}

# 	print Dumper($umi_data);
# 	print Dumper($cdr3_patterns);
# 	print Dumper($tcr_patterns);exit;

	# Tests 100 reads to find the correct sense of the reads (the one with higher number of CDR3s annotated
	if (in_array($options, 'auto') || in_array($options, 'auto detect')){
		my $reads_to_check = 1000;
		if ($total_reads < $reads_to_check) {
			$reads_to_check = $total_reads;
		}
		my %count_readsense_cdr3s;
		if (defined($read_files->[1])){
			%count_readsense_cdr3s = ('R1+'=>0, 'R1-'=>0, 'R2+'=>0, 'R2-'=>0);
		} else {
			%count_readsense_cdr3s = ('R1+'=>0, 'R1-'=>0);
		}
		#my %count_readsense_vgene_lengths = ('R1+'=>0, 'R1-'=>0, 'R2+'=>0, 'R2-'=>0);
		my %count_readsense_vgenes = ('R1+'=>0, 'R1-'=>0, 'R2+'=>0, 'R2-'=>0);
		my %count_readsense_unique_vgenes = ('R1+'=>0, 'R1-'=>0, 'R2+'=>0, 'R2-'=>0);
		foreach my $readsense (keys %count_readsense_cdr3s){
			my ($cdr3_headers_, $cdr3_seqs_, $cdr3_qualities_, $total_tcr_reads_, $tcr_segments_) = detect_cdr3_seqs_from_files($read_files,$cdr3_patterns,$tcr_patterns,$umi_data,[@$options,$readsense],$reads_to_check);
			if (defined($total_tcr_reads_)){
				map $count_readsense_cdr3s{$readsense} += $_, values %$total_tcr_reads_;
			}
			if (defined($tcr_segments_) && %$tcr_segments_ && %INP_ref_files && defined($INP_ref_files{'tcrv'})){
# 				my $sample_ = (keys %$tcr_segments_)[0];
# 				my @segment_seq_lengths;
# 				while ((my $marker_, my $segment_seqs_) = each(%{$tcr_segments_->{$sample_}})){
# 					while ((my $cdr3_, my $cdr3_segment_seqs_) = each(%$segment_seqs_)){
# 						push(@segment_seq_lengths,length(most_frequent(@{$cdr3_segment_seqs_->{'tcrv'}})));
# 					}
# 				}
# 				$count_readsense_vgene_lengths{$readsense} = mean(@segment_seq_lengths);
				my $sample_ = (keys %$tcr_segments_)[0];
				# map $count_readsense_segments{$readsense} += scalar keys %{$_}, values %{$tcr_segments_->{$sample_}};
				my $align_params = { 'alignment' => 'dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 90', 'aligned' => 0.1, 'ident' => 0.9 };
				my $cdr3_to_segment_seqs;
				while ((my $marker_, my $segment_seqs_) = each(%{$tcr_segments_->{$sample_}})){
					while ((my $cdr3_, my $cdr3_segment_seqs_) = each(%$segment_seqs_)){
						$cdr3_to_segment_seqs->{$cdr3_} = most_frequent(@{$cdr3_segment_seqs_->{'tcrv'}});
					}
				}
				my $cdr3_to_tcrv_gene = match_alleles($INP_ref_files{'tcrv'},$cdr3_to_segment_seqs,undef,$align_params);
				$count_readsense_vgenes{$readsense} = scalar keys %{$cdr3_to_tcrv_gene};
				foreach my $vgene (values %{$cdr3_to_tcrv_gene}){
					# my $ngenes = () = $vgene =~ /\|/gi;
					# $count_readsense_unique_vgenes{$readsense} += 1/($ngenes+1);
					if ($vgene !~ /\|/) {
						$count_readsense_unique_vgenes{$readsense}++;
					}
				}
			}
		}
		my @best_readsenses = sort {$count_readsense_cdr3s{$b}<=>$count_readsense_cdr3s{$a}} keys %count_readsense_cdr3s;
		my $best_readsense;
		if ($count_readsense_cdr3s{$best_readsenses[0]} < 0.5*$reads_to_check) {
			die "\n ERROR: cannot detect the position of the CDR3 in the reads. Specify manually the sense of the read containing the CDR3, could be also required the primer/amplicon data.\n\n";
		# Checks if any read sense retrieves more CDR3s than the others
		} elsif ($count_readsense_cdr3s{$best_readsenses[0]} > $count_readsense_cdr3s{$best_readsenses[1]}){
			$best_readsense = $best_readsenses[0];
# 		# If so, retrieves the read sense with longer V gene segment sequences
# 		} else {
# 			@best_readsenses = sort {$count_readsense_vgene_lengths{$b}<=>$count_readsense_vgene_lengths{$a}} ($best_readsenses[0], $best_readsenses[1]);
# 			$best_readsense = $best_readsenses[0];
# 		}
# 		$options = [@$options, $best_readsense];
		} else {
			@best_readsenses = sort {$count_readsense_vgenes{$b}<=>$count_readsense_vgenes{$a}} ($best_readsenses[0], $best_readsenses[1]);
			# Checks if one read sense has more V genes assigned than the other
			if ($count_readsense_vgenes{$best_readsenses[0]} > $count_readsense_vgenes{$best_readsenses[1]}){
				$best_readsense = $best_readsenses[0];
			# If both read senses retrieve the same number of CDR3s and have the same number of V genes assigned, 
			# then it checks the number of unambiguously assigned V genes
			} else {
				@best_readsenses = sort {$count_readsense_unique_vgenes{$b}<=>$count_readsense_unique_vgenes{$a}} ($best_readsenses[0], $best_readsenses[1]);
				$best_readsense = $best_readsenses[0];
			}
		}
		$options = [@$options, $best_readsense];
	}

	# Processes the read files and extracts the CDR3 sequences
	return detect_cdr3_seqs_from_files($read_files,$cdr3_patterns,$tcr_patterns,$umi_data,$options);


}


#################################################################################

sub extract_cdr3_seqs_with_threads {

	my ($read_files,$markerdata,$sampledata,$total_reads,$tcr_ref_patterns,$options,$threads_limit) = @_; # ,$tcr_ref_files,$options,$threads_limit) = @_;

	my ($cdr3_headers, $cdr3_seqs, $cdr3_qualities, $total_tcr_reads, $tcr_segments);

	push(@$options, 'quiet');
	
	# Splits reads file into smaller files to process them in parallel
	my @splitted_read_files1 = split_fastq_file($read_files->[0],$threads_limit,undef,'gzip');
	my @splitted_read_files2;
	if (scalar @$read_files == 2){
		@splitted_read_files2 = split_fastq_file($read_files->[1],$threads_limit,undef,'gzip');
	}

	my @threads;
	for (my $i=0; $i<=$#splitted_read_files1; $i++){
		my $read_files = [ $splitted_read_files1[$i], $splitted_read_files2[$i] ];
		push(@threads, threads->create(\&extract_cdr3_seqs,$read_files,$markerdata,$sampledata,$total_reads,$tcr_ref_patterns,$options)); # ,$tcr_ref_files,$options));
# 		push(@threads, [extract_cdr3_seqs($read_files,$markerdata,$sampledata,$total_reads,$tcr_ref_patterns,$options)]);
	}
	while (){
		for (my $i=0; $i<=$#threads; $i++){
			unless ($threads[$i]->is_running()){
				my ($cdr3_headers_, $cdr3_seqs_, $cdr3_qualities_, $total_tcr_reads_, $tcr_segments_) = $threads[$i]->join;
				if (defined($cdr3_headers_)){
					foreach my $sample (keys %$cdr3_headers_){
						foreach my $marker (keys %{$cdr3_headers_->{$sample}}){
							if (defined($cdr3_headers_->{$sample}{$marker}) && @{$cdr3_headers_->{$sample}{$marker}}){
								push(@{$cdr3_headers->{$sample}{$marker}}, @{$cdr3_headers_->{$sample}{$marker}});
								push(@{$cdr3_seqs->{$sample}{$marker}}, @{$cdr3_seqs_->{$sample}{$marker}});
								if (defined($cdr3_qualities_->{$sample}{$marker})){
									push(@{$cdr3_qualities->{$sample}{$marker}}, @{$cdr3_qualities_->{$sample}{$marker}});
								}
								if (defined($tcr_segments_) && defined($tcr_segments_->{$sample}{$marker})){
									if (!defined($tcr_segments->{$sample}{$marker})){
										$tcr_segments->{$sample}{$marker} = {};
									}
									$tcr_segments->{$sample}{$marker} = { %{$tcr_segments->{$sample}{$marker}}, %{$tcr_segments_->{$sample}{$marker}} };
								}
								$total_tcr_reads->{$sample} += $total_tcr_reads_->{$sample};
							}
						}
					}
				}
				undef($threads[$i]);
				splice(@threads,$i,1);
				$i--;
			}
		}
		if (@threads){
			sleep(10);
		} else {
			last;
		}
	}

# 	# For debugging:
# 	for (my $i=0; $i<=$#threads; $i++){
# 		print '';
# 		my ($cdr3_headers_, $cdr3_seqs_, $cdr3_qualities_, $total_tcr_reads_, $tcr_segments_) = @{$threads[$i]};
# 		if (defined($cdr3_headers_)){
# 			foreach my $sample (keys %$cdr3_headers_){
# 				foreach my $marker (keys %{$cdr3_headers_->{$sample}}){
# 					if (defined($cdr3_headers_->{$sample}{$marker}) && @{$cdr3_headers_->{$sample}{$marker}}){
# 						push(@{$cdr3_headers->{$sample}{$marker}}, @{$cdr3_headers_->{$sample}{$marker}});
# 						push(@{$cdr3_seqs->{$sample}{$marker}}, @{$cdr3_seqs_->{$sample}{$marker}});
# 						if (defined($cdr3_qualities_->{$sample}{$marker})){
# 							push(@{$cdr3_qualities->{$sample}{$marker}}, @{$cdr3_qualities_->{$sample}{$marker}});
# 						}
# 						if (defined($tcr_segments_) && defined($tcr_segments_->{$sample}{$marker})){
# 							if (!defined($tcr_segments->{$sample}{$marker})){
# 								$tcr_segments->{$sample}{$marker} = {};
# 							}
# 							$tcr_segments->{$sample}{$marker} = { %{$tcr_segments->{$sample}{$marker}}, %{$tcr_segments_->{$sample}{$marker}} };
# 						}
# 						$total_tcr_reads->{$sample} += $total_tcr_reads_->{$sample};
# 					}
# 				}
# 			}
# 		}
# 		delete $threads[$i];
# 	}

	# Removes splitted files
	foreach my $file (@splitted_read_files1, @splitted_read_files2) {
		`rm $file`;
	}

	return ($cdr3_headers, $cdr3_seqs, $cdr3_qualities, $total_tcr_reads, $tcr_segments);

}


#################################################################################

sub detect_cdr3_seqs_from_files {

	my ($read_files,$cdr3_patterns,$tcr_patterns,$umi_data,$options,$limit) = @_;

	my %total_tcr_reads;
	my ($cdr3_headers, $cdr3_seqs, $cdr3_qualities);
	my ($nocdr3_headers, $nocdr3_seqs, $nocdr3_qualities);
	my $tcr_segments;

	my ($output_format,$readsense,$nocdr3,$revcomp,$verbose)=('fasta','R1+',0,0,1);
	# Retrieves CDR3 sequences in FASTQ format, including base calling qualities
	if (in_array($options, 'fastq')){
		$output_format = 'fastq';
	}
	# Retrieves sequences that do not match any CDR3 pattern
	if (in_array($options, 'nocdr3')){
		$nocdr3 = 1;
	}
	# Processes the sequences in direct and reverse complementary senses
	if (in_array($options, 'revcomp')){
		$revcomp = 1;
	}
	# Indicates the sense of the read containing the CDR3 region
	if (in_array($options, 'R1-')){
		$readsense = 'R1-';
	} elsif (in_array($options, 'R2+')){
		$readsense = 'R2+';
	} elsif (in_array($options, 'R2-')){
		$readsense = 'R2-';
	} else {
		$readsense = 'R1+';
	}

	# Open reads files
	if (is_gzip($read_files->[0]) || is_zip($read_files->[0])) {
		open(READSFILE1, "zcat $read_files->[0] |") || die "\nERROR: cannot open compressed '$read_files->[0]'\n\n";
	} else {
		open(READSFILE1, $read_files->[0]) || die "\n ERROR: cannot open '$read_files->[0]'\n\n";
	}
	if (defined($read_files->[1])){
		if (is_gzip($read_files->[1]) || is_zip($read_files->[1])) {
			open(READSFILE2, "zcat $read_files->[1] |") || die "\nERROR: cannot open compressed '$read_files->[1]'\n\n";
		} else {
			open(READSFILE2, $read_files->[1]) || die "\n ERROR: cannot open '$read_files->[1]'\n\n";
		}
	}

# 	my @random_reads;
# 	if (defined($limits)){
# 		if (@$limits && defined($limits->[1])){
# 			my @nums = 1 .. $limits->[1];
# 			@random_reads = map splice(@nums, int(rand @nums),1) , 1 .. $limits->[0];
# 		} elsif (@$limits) {
# 			@random_reads = 1..$limits->[0];
# 		} else {
# 			@random_reads = 1..$limits;
# 		}
# 	}
# 	@random_reads = sort {$b<=>$a} @random_reads;
# 	my $limit = $random_reads[0];
# # 	my %random_reads;
# # 	map $random_reads{$_}=1 , @random_reads;

	my (@lines,@headers,@seqs,@qualities);
	my $count_reads = 0;
	my $count_cdr3 = 0;
	my %count_cdr3; # Counts CDR3 seqs per sample
	map $count_cdr3{$_}=0 , keys %$cdr3_patterns;
	my $count_line = 0;
	while ($lines[0] = <READSFILE1>) {
		if (defined($read_files->[1])){
			$lines[1] = <READSFILE2>;
		}
		# Counts lines
		if ($count_line % 4 == 0){
			$count_line = 1;
		} else {
			$count_line++;
		}
		# Reads headers in FASTA or FASTQ format
		if ($count_line == 1 && $lines[0] =~ /^[@>](.+)/){ # (/^[@>](.+)\n/){
			$headers[0] = $1;
			next;
		# Reads sequences
		} elsif ($count_line == 2){
			$lines[0] =~ s/\012\015?|\015\012?|^\s+|\s+$//; # trim spaces and line break
			$lines[0] =~ s/[^ACGT]/N/ig; # trim spaces and line break
			$seqs[0] = $lines[0];
			if (defined($lines[1])){
				$lines[1] =~ s/\012\015?|\015\012?|^\s+|\s+$//; # trim spaces and line break
				$lines[1] =~ s/[^ACGT]/N/ig; # trim spaces and line break
				$seqs[1] = $lines[1];
			}
			if ($output_format eq 'fastq'){
				next;
			}
		# Skips non sequence lines or TCR sequences with undefined nucleotides
	# 	} elsif ($count_line != 2 && $output_format eq 'fasta'){ # (!defined($header) || !/^[ACGTU]+$/i){
	# 		next;
		# Reads qualities if the output will be in FASTQ format
		} elsif ($count_line == 4 && $output_format eq 'fastq'){
			$lines[0] =~ s/\012\015?|\015\012?//; # trim line break
			$qualities[0] = $lines[0];
			if (defined($lines[1])){
				$lines[1] =~ s/\012\015?|\015\012?//; # trim line break
				$qualities[1] = $lines[1];
			}
		} else {
			next;
		}
		$count_reads++;
# 		if (defined($limits)){
# 			if ($count_reads != $random_reads[-1]){
# 				next;
# 			} else {
# 				pop(@random_reads);
# 			}
# 		}
# 		# Prints number of processed reads
# 		if ($count_reads % 100000 == 0) {
# 			if ($count_reads == 100000) {
# 				print "\n\tReads\tCDR3s\n";
# 			}
# 			print "\t$count_reads\t$count_cdr3\n";
# 		}

		# Loops all the samples, all the markers and all the primers/amplicons to find the CDR3 region of the sequence to annotate it
		my ($cdr3_found,$cdr3_sample,$cdr3_marker,$cdr3_seq,$cdr3_header,$cdr3_quality,$tcrv_seq,$tcrj_seq) = 
		detect_cdr3_seq(\@seqs,\@headers,\@qualities,$cdr3_patterns,$tcr_patterns,$umi_data,$readsense,$revcomp);

		if ($cdr3_found) {
			push(@{$cdr3_seqs->{$cdr3_sample}{$cdr3_marker}}, $cdr3_seq);
			push(@{$cdr3_headers->{$cdr3_sample}{$cdr3_marker}}, $cdr3_header);
			if (defined($cdr3_quality)){
				push(@{$cdr3_qualities->{$cdr3_sample}{$cdr3_marker}}, $cdr3_quality);
			}
			if (defined($tcrv_seq)){
				push(@{$tcr_segments->{$cdr3_sample}{$cdr3_marker}{$cdr3_seq}{'tcrv'}}, $tcrv_seq);
			}
			if (defined($tcrj_seq)){
				push(@{$tcr_segments->{$cdr3_sample}{$cdr3_marker}{$cdr3_seq}{'tcrj'}}, $tcrj_seq);
			}
			$total_tcr_reads{$cdr3_sample}++;
			$count_cdr3{$cdr3_sample}++;
			$count_cdr3++;
		} elsif ($cdr3_found == -1){
			$total_tcr_reads{$cdr3_sample}++;
			if ($nocdr3){
				push(@{$nocdr3_headers->{$cdr3_sample}{$cdr3_marker}}, $cdr3_header);
				push(@{$nocdr3_seqs->{$cdr3_sample}{$cdr3_marker}},$cdr3_seq);
				if (defined($cdr3_quality)){
					push(@{$nocdr3_qualities->{$cdr3_sample}{$cdr3_marker}},$cdr3_quality);
				}
			}
		}
		undef(@seqs);
		undef(@headers);
		undef(@qualities);
		if (defined($limit) && $limit == $count_reads){
			last;
		}
	}

	close READSFILE1;
	if (defined($read_files->[1])){
		close READSFILE2;
	}

	if (!$nocdr3){
		return ($cdr3_headers, $cdr3_seqs, $cdr3_qualities, \%total_tcr_reads, $tcr_segments);
	} else {
		return ($nocdr3_headers, $nocdr3_seqs, $nocdr3_qualities, \%total_tcr_reads, $tcr_segments);
	}


}

#################################################################################

sub detect_cdr3_seq {

	my ($seqs,$headers,$qualities,$cdr3_patterns,$tcr_patterns,$umi_data,$readsense,$revcomp) = @_;
	
	my $umi = 0;
	if (defined($umi_data)){
		$umi = 1;
	}

	my $header = $headers->[0];
	# Concatenates sequences from both reads in the desired orientation
	my ($seq,$quality);
	if ($readsense eq 'R1-'){
		if (defined($seqs->[1])) {
			$seq = $seqs->[1].'-----'.iupac_reverse_complementary($seqs->[0]);
			if (@$qualities){
				$quality = $qualities->[1].'-----'.reverse_sequence($qualities->[0]);
			}
		} else {
			$seq = iupac_reverse_complementary($seqs->[0]);
			if (@$qualities){
				$quality = reverse_sequence($qualities->[0]);
			}
		}
	} elsif ($readsense eq 'R2-'){
		$seq = $seqs->[0].'-----'.iupac_reverse_complementary($seqs->[1]);
		if (@$qualities){
			$quality = $qualities->[0].'-----'.reverse_sequence($qualities->[1]);
		}
	} elsif ($readsense eq 'R2+'){
		$seq = iupac_reverse_complementary($seqs->[0]).'-----'.$seqs->[1];
		if (@$qualities){
			$quality = reverse_sequence($qualities->[0]).'-----'.$qualities->[1];
		}
	} else {
		if (defined($seqs->[1])) {
			$seq = iupac_reverse_complementary($seqs->[1]).'-----'.$seqs->[0];
			if (@$qualities){
				$quality = reverse_sequence($qualities->[1]).'-----'.$qualities->[0];
			}
		} else {
			$seq = $seqs->[0];
			if (@$qualities){
				$quality = $qualities->[0];
			}
		}
	}

	my ($cdr3_sample,$cdr3_marker,$cdr3_seq,$cdr3_header,$cdr3_quality,$tcrv_seq,$tcrj_seq);

	my $cdr3_found = 0;
	foreach my $sample (keys %$cdr3_patterns){
		foreach my $marker (keys %{$cdr3_patterns->{$sample}}){
# 			if (defined($limit_per_sample) && $count_cdr3{$sample} >= $limit_per_sample){
# 				next;
# 			}
			foreach my $amplicon_id (keys %{$cdr3_patterns->{$sample}{$marker}}){
				my $tcr_pattern = $tcr_patterns->{$sample}{$marker}{$amplicon_id};
				my $cdr3_pattern = $cdr3_patterns->{$sample}{$marker}{$amplicon_id};
				my ($umi_seq, $umi_len, $umi_pos);
				if ($umi){
					$umi_len = $umi_data->{$sample}{$marker}{$cdr3_pattern}{'len'};
					$umi_pos = $umi_data->{$sample}{$marker}{$cdr3_pattern}{'pos'};
				}
# if ($marker eq 'TCRB' && $readsense eq 'R1-'){
# 	print "$sample\n$seq\n$cdr3_pattern\n\n";
# 	print "";
# 
# }
				if ($seq =~ m/$cdr3_pattern/){
					$cdr3_found = 1;
					$cdr3_seq = substr($seq,($-[1])-$INP_extra_length,($+[1])-($-[1])+2*$INP_extra_length);
# if ($marker eq 'TCRB' && $cdr3_seq eq 'TGTGCCAGCAGCTCGGGACTGCGGAGGGGCAATGAGCAGTTCTTC'){
# print "$sample\n$cdr3_seq\n$seq\n$cdr3_pattern\n\n";
# exit;
# }
					# @tcrv_seqs = ( substr($seq, 0, $tcrv_len), substr($seq, $tcrv_len, $-[1]) );
					$tcrv_seq = substr($seq, 0, $-[1]);
					$tcrj_seq = substr($seq, $+[1], 30);
					if ($quality){
						$cdr3_quality = substr($quality, $-[1], $+[1] - $-[1]);
					}
					if ($umi && $umi_len > 0){
						$umi_seq = substr($seq, $-[0]+$umi_pos-1, $umi_len);
						# @tcrv_seqs = ( substr($seq, $-[0]+$umi_pos-1+$umi_len, $tcrv_len-($umi_pos-1+$umi_len)), substr($seq, $tcrv_len, $-[1]) );
					} elsif ($umi && $umi_len < 0){
						$umi_seq = iupac_reverse_complementary(substr($seq, $+[0]+$umi_pos+$umi_len-1, -$umi_len));
					}
				} elsif ($revcomp == 1){
					my $revcomp_seq = iupac_reverse_complementary($seq);
					if ($revcomp_seq =~ /$cdr3_pattern/){
						$cdr3_found = 1;
						$cdr3_seq = substr($revcomp_seq,($-[1])-$INP_extra_length,($+[1])-($-[1])+2*$INP_extra_length);
						# @tcrv_seqs = ( substr($revcomp_seq, 0, $tcrv_len), substr($revcomp_seq, $tcrv_len, $-[1]) );
						$tcrv_seq = substr($revcomp_seq, 0, $-[1]);
						$tcrj_seq = substr($revcomp_seq, $+[1], 30);
						if ($quality){
							$cdr3_quality = substr(reverse_sequence($quality), $-[1], $+[1] - $-[1]);
						}
						if ($umi && $umi_len > 0){
							$umi_seq = substr($revcomp_seq, $-[0], $umi_len);
							# @tcrv_seqs = ( substr($revcomp_seq, $-[0]+$umi_pos-1+$umi_len, length($seq)-$tcrv_len-($umi_pos-1+$umi_len)), substr($revcomp_seq, length($seq)-$tcrv_len, $-[1]) );
						} elsif ($umi && $umi_len < 0){
							$umi_seq = iupac_reverse_complementary(substr($revcomp_seq, $+[0]+$umi_len, -$umi_len));
						}
					}
				}
				if ($cdr3_found) {
					$cdr3_sample = $sample;
					$cdr3_marker = $marker;
					$cdr3_header = sprintf("%s | ampli=%s", $header, $amplicon_id);
					if (defined($umi_seq)){
						$cdr3_header .= sprintf(" | umi=%s", $umi_seq);
					}
					if (defined($tcrv_seq)){
						$cdr3_header .= sprintf(" | tcrv=%s", $tcrv_seq);
					}
					if (defined($tcrj_seq)){
						$cdr3_header .= sprintf(" | tcrj=%s", $tcrj_seq);
					}
					last;
				} elsif (defined($tcr_pattern) && $seq =~ /$tcr_pattern/){
					$cdr3_sample = $sample;
					$cdr3_marker = $marker;
					$cdr3_header = "$header | ampli=$amplicon_id";
					$cdr3_seq = $1;
					if ($quality){
						$cdr3_quality = substr($quality, $-[1], $+[1] - $-[1]);
					}
					$cdr3_found = -1;
					last;
				}
			}
			if ($cdr3_found != 0) {
				last;
			}
		}
		if ($cdr3_found != 0) {
			last;
		}
	}

	return ($cdr3_found,$cdr3_sample,$cdr3_marker,$cdr3_seq,$cdr3_header,$cdr3_quality,$tcrv_seq,$tcrj_seq);
}


#################################################################################

# Clusters of similar sequences
sub cluster_cdr3s {

	my ($cdr3_seqs,$cdr3_headers,$options) = @_;

	my ($seq_clusters, $seq_stats) = ({},{});

	my (@cdr3_seqs_,@cdr3_headers_);

	# First group sequences with the same UMI to remove artefacts
	# All sequences with the same UMI will be assigned to the dominant one
	if (defined($options->{'umi_errors'})) {
		# Group sequences with the same UMI
		my $same_length_seqs;
		my $umi_cdr3_seqs;
		for (my $i=0; $i<=$#{$cdr3_headers}; $i++){
			if ($cdr3_headers->[$i] =~ /umi=([ACGT]+)/){
				push(@{$umi_cdr3_seqs->{$1}},$i);
			}
		}

		# Do not sort if we want to take a subset of random UMIs
		my @sorted_umis;
		if (!defined($options->{'umi_limit'})) {
			@sorted_umis = sort { $#{$umi_cdr3_seqs->{$b}} <=> $#{$umi_cdr3_seqs->{$a}} } keys %$umi_cdr3_seqs;
		} else {
			@sorted_umis = keys %$umi_cdr3_seqs;
		}

		# Clusters sequences within the same UMI based in a max. number of substitutions
		foreach my $umi (@sorted_umis){
			# Singleton UMIs
			if (scalar @{$umi_cdr3_seqs->{$umi}} == 1){
				# Annotates number of sequences per UMI
				$seq_stats->{'reads'}++;
				# Annotates number of variants per UMI
				$seq_stats->{'variants'}++;
				# Annotates singleton UMIs
				$seq_stats->{'umis_singletons'}++;
				if (!defined($options->{'keep_singletons'})){
					next;
				}
			} else {
				my (@umi_seqs, @umi_headers);
				foreach my $i (@{$umi_cdr3_seqs->{$umi}}){
					push(@umi_seqs, $cdr3_seqs->[$i]);
					push(@umi_headers, $cdr3_headers->[$i]);
				}
				# Annotates number of sequences per UMI
				$seq_stats->{'reads'}{scalar @umi_seqs}++;
				# Annotates number of variants per UMI
				$seq_stats->{'variants'}{scalar unique(@umi_seqs)}++;

				# Clusters UMI sequences, 1 cluster will be created with identical/similar sequences
				my $umi_clusters_ = cluster_umi_seqs(\@umi_seqs,\@umi_headers,$options->{'umi_errors'});
				
				# Checks if at least half of the UMI sequences are similar to the dominant one
				# If not the UMI is discarded
				my $count_clustered_seqs = scalar @{$umi_clusters_->{(keys %$umi_clusters_)[0]}};
				# Counts the number of incorrect UMIs (less than half sequences after clustering)
				# UMIs that are singletons are removed unless it is specified
				# if ($count_clustered_seqs==1 || 2*$count_clustered_seqs <= scalar @umi_seqs) {
				if (2*$count_clustered_seqs <= scalar @umi_seqs){
					$seq_stats->{'umis_incorrect'}++;
					if (!defined($options->{'keep_singletons'})){
						next;
					}
				} else {
					$seq_stats->{'umis_correct'}++;
				}
				
				# Annotates the sequences from the first UMI cluster assigning the sequence of the dominant one
				# Only 1 $dominant_seq, UMIs clustering generates a unique cluster
				my $dominant_seq = (keys %$umi_clusters_)[0];
				foreach my $header (@{$umi_clusters_->{$dominant_seq}}){
					push(@cdr3_headers_,$header);
					# Includes in sequences array the dominant sequence for further clustering
					push(@cdr3_seqs_,$dominant_seq);
					$count_clustered_seqs++;
				}
	# 			# Corrects all the UMI clustered sequences assigning the sequence of the dominant one
	# 			foreach my $dominant_seq (keys %$umi_clusters_) { # Only 1 $dominant_seq, UMIs clustering generates a unique cluster
	# 				foreach my $header (@{$umi_clusters_->{$dominant_seq}}){
	# 					push(@cdr3_headers_,$header);
	# 					# Includes in sequences array the dominant sequence for further clustering
	# 					push(@cdr3_seqs_,$dominant_seq);
	# 					$count_clustered_seqs++;
	# 				}
	# 			}
			}
			# Stops if a fix number of UMIs per sample is desired
			if (defined($options->{'umi_limit'}) && $seq_stats->{'umis_correct'} == $options->{'umi_limit'}) {
				last;
			}
		}
	} else {
		@cdr3_seqs_ = @$cdr3_seqs;
		@cdr3_headers_ = @$cdr3_headers;
	}

	# Group same length sequences for faster clustering later
	my $same_length_seqs;
	for (my $i=0; $i<=$#cdr3_seqs_; $i++){
		my $len = length($cdr3_seqs_[$i]);
		push(@{$same_length_seqs->{$len}},$i);
	}
# 	my ($clustered_cdr3_seqs, $clustered_cdr3_headers);
# 	foreach my $len (sort {$b<=>$a} keys %{$same_length_seqs}){
	foreach my $len (sort { $#{$same_length_seqs->{$b}}<=>$#{$same_length_seqs->{$a}} } keys %{$same_length_seqs}){
		my @same_length_seqs = map $cdr3_seqs_[$_], @{$same_length_seqs->{$len}};
		my @same_length_headers = map $cdr3_headers_[$_], @{$same_length_seqs->{$len}};
		my $seq_clusters_ = cluster_illumina_seqs(\@same_length_seqs,\@same_length_headers,$options->{'errors'});
		$seq_clusters = { %$seq_clusters, %$seq_clusters_ };
# 		my ($clustered_seqs,$clustered_headers) = cluster_illumina_seqs(\@seqs_,\@headers_,$options->{'errors'});
# 		push(@{$clustered_cdr3_seqs}, @$clustered_seqs);
# 		push(@{$clustered_cdr3_headers}, @$clustered_headers);
	}
# 	return ($clustered_cdr3_seqs, $clustered_cdr3_headers);

	return ($seq_clusters, $seq_stats);

}

#################################################################################

# Clusters of similar sequences
sub cluster_sample_cdr3s_with_threads {

	my ($cdr3_seqs,$cdr3_headers,$options,$threads_limit) = @_;

	my ($seq_clusters, $seq_stats) = ({},{});

	my (@threads, @threads_samples, @threads_markers);
	my @samples = keys %$cdr3_seqs;
	foreach my $sample (@samples){
		my @markers = keys %{$cdr3_seqs->{$sample}};
		foreach my $marker (@markers){
			push(@threads, threads->create(\&cluster_cdr3s,$cdr3_seqs->{$sample}{$marker},$cdr3_headers->{$sample}{$marker},$options));
			push(@threads_samples, $sample); # It's required to know which sample is analyzed by each thread
			push(@threads_markers, $marker); # It's required to know which marker is analyzed by each thread
	# 		# For debugging:
	# 		push(@threads, [cluster_cdr3s($cdr3_seqs->{$sample}{$marker},$cdr3_headers->{$sample}{$marker},$options)]);
			# If maximum number of threads is reached or last sample is processed
			if (scalar @threads >= $threads_limit  || $sample eq $samples[-1]){
				my $check_threads = 1;
				while ($check_threads){
					for (my $i=0; $i<=$#threads; $i++){
						unless ($threads[$i]->is_running()){

							($seq_clusters->{$threads_samples[$i]}{$threads_markers[$i]}, $seq_stats->{$threads_samples[$i]}{$threads_markers[$i]})= $threads[$i]->join;
	# 						($cdr3_seqs->{$sample}{$marker},$cdr3_headers->{$sample}{$marker}) = $threads[$i]->join;
		
							undef($threads[$i]);
							splice(@threads,$i,1);
							splice(@threads_samples,$i,1);
							splice(@threads_markers,$i,1);
							$i = $i - 1;
							unless (($sample eq $samples[-1] && $marker eq $markers[-1]) && @threads){
								$check_threads = 0;
							}
						}
					}
					if ($check_threads){
						sleep(1);
					}
				}
	# 			# For debugging:
	# 			for (my $i=0; $i<=$#threads; $i++){
	# 				($seq_clusters->{$threads_samples[$i]}{$threads_markers[$i]}, $seq_stats->{$threads_samples[$i]}{$threads_markers[$i]}) = @{$threads[$i]};
	# 				delete $threads[$i];
	# 			}
			}
		}
	}
# 	print "\nOpen threads: ".scalar @threads." \n";
	
	return ($seq_clusters,$seq_stats);
# 	return ($cdr3_seqs,$cdr3_headers);
}

#################################################################################

# Calculates common clustered CDR3 sequences
# between different samples from the same individual
# or between all the samples from different individuals
sub common_cdr3_seqs {

	my ($cdr3_seqs,$substitutions,$type) = @_;

	if (!defined($type)){
		$type = 'intraindividual';
	}
	
	my $common_cdr3_seqs;
	my %markers;

	# Annotates duplicated samples/replicates per individual
	my $dup_samples;
	foreach my $sample (keys %$cdr3_seqs){
		if ($sample =~ /(.+?)[-_](.+)/){
			push(@{$dup_samples->{$1}},$sample)
		}
		foreach my $marker (keys %{$cdr3_seqs->{$sample}}){
			$markers{$marker} = 1;
		}
	}
	# TRICK:
# 	if ($type eq 'interindividual' && scalar keys %$dup_samples == 1){
# 		undef($dup_samples);
# 		foreach my $sample (keys %$cdr3_seqs){
# 			push(@{$dup_samples->{$sample}},$sample)
# 		}	
# 	}

	# Each marker must be evaluated individually
	foreach my $marker (keys %markers){

		if ($type eq 'intraindividual') {
			# my $max_samples_per_ind = max(map scalar @{$dup_samples->{$_}}, keys %$dup_samples);
			foreach my $individual (nsort(keys %$dup_samples)) {
				# @all_seqs stores sequences from all the samples/replicates of the same individual
				# $sample_seqs->{$sample}{$marker} is a hash with sequences from a single sample/replicate
				my (@all_seqs, $sample_seqs);
				foreach my $sample (@{$dup_samples->{$individual}}) {
					if (defined($cdr3_seqs->{$sample}{$marker})){
						push(@all_seqs, @{$cdr3_seqs->{$sample}{$marker}});
						$sample_seqs->{$sample}{$marker} = { map { $_ => undef } @{$cdr3_seqs->{$sample}{$marker}} };
					}
				}
				@all_seqs = unique(@all_seqs);
				foreach my $seq (@all_seqs) {
					my @common_samples;
					foreach my $sample (@{$dup_samples->{$individual}}) {
						if (exists($sample_seqs->{$sample}{$marker}{$seq})){
							push(@common_samples, $sample);
						}
					}
					$common_cdr3_seqs->{$marker}{$individual}{'total'}{scalar @common_samples}++;
					my %common_annotated;
					foreach my $sample1 (nsort(@common_samples)) {
						foreach my $sample2 (nsort(@common_samples)) {
							if ($sample1 ne $sample2 && !defined($common_annotated{"$sample2-$sample1"})) {
								$common_cdr3_seqs->{$marker}{$individual}{'pairs'}{"$sample1,$sample2"}++;
								$common_annotated{"$sample1-$sample2"} = 1;
							}
						}
					}
				}
		# 		foreach my $count_samples (sort {$b<=>$a} keys %{$common_cdr3_seqs->{$marker}{$individual}{'total'}}) {
		# 			printf("Individual '%s': %8d sequences are shared by %d samples.\n", $individual, $common_cdr3_seqs->{$marker}{$individual}{'total'}{$count_samples}, $count_samples);
		# 		}

				# Checks for common similar sequences between samples (that differ few substitutions)
				if (defined($substitutions) && $substitutions>0){
					my $count_seqs = 0;
					foreach my $seq (@all_seqs) {
						$count_seqs++;
						if ($count_seqs % 1000 == 0) {
							print "$count_seqs ";
						}
						my $common_samples = 0;
						foreach my $sample (@{$dup_samples->{$individual}}) {
							# If the identical sequence is found, it is annotated and goes to check the next sample
							if (exists($sample_seqs->{$sample}{$marker}{$seq})){
								$common_samples++;
								delete($sample_seqs->{$sample}{$marker}{$seq});
								next;
							}
							# Checks for similar sequences from the sample
							foreach my $seq_ (keys %{$sample_seqs->{$sample}{$marker}}) {
								my $diff = count_errors($seq,$seq_,$substitutions+1);
								# If a similar one is found, it is annotated and goes to check the next sample
								if (defined($diff) && $diff <= $substitutions) {
									$common_samples++;
									delete($sample_seqs->{$sample}{$marker}{$seq_});
									last;
								}
							}

						}
						if ($common_samples) {
							$common_cdr3_seqs->{$marker}{$individual}{'total_sim'}{$common_samples}++;
						}
					}
					# print "\n";
				}
			}

		} elsif ($type = 'interindividual') {
			
			# my $max_samples_per_ind = max(map scalar @{$dup_samples->{$_}}, keys %$dup_samples);
			# $sample_seqs->{$individual} is a hash with all sequences from all the samples/replicates of the same individual
			my $sample_seqs;
			foreach my $individual (nsort(keys %$dup_samples)) {
				# @all_seqs stores sequences from all the samples/replicates of the same individual
				my @all_seqs;
				foreach my $sample (@{$dup_samples->{$individual}}) {
					if (defined($cdr3_seqs->{$sample}{$marker})){
						push(@all_seqs, @{$cdr3_seqs->{$sample}{$marker}});
					}
				}
				@all_seqs = unique(@all_seqs);
				$sample_seqs->{$individual} = { map { $_ => undef } @all_seqs };
			}
			my %individual_seqs;
			foreach my $individual1 (nsort(keys %$dup_samples)) {
				foreach my $seq1 (keys %{$sample_seqs->{$individual1}}) {
					$individual_seqs{$seq1}++;
				}
				foreach my $individual2 (nsort(keys %$dup_samples)) {
					if ($individual1 eq $individual2) {
						$common_cdr3_seqs->{$marker}{$individual1}{$individual2} = scalar keys %{$sample_seqs->{$individual1}};
						next;
					}
					$common_cdr3_seqs->{$marker}{$individual1}{$individual2} = 0;
					foreach my $seq1 (keys %{$sample_seqs->{$individual1}}) {
						foreach my $seq2 (keys %{$sample_seqs->{$individual2}}) {
							if ($seq1 eq $seq2) {
								$common_cdr3_seqs->{$marker}{$individual1}{$individual2}++;
								last;
							}
						}
					}
				}
			}
			if (scalar keys %$dup_samples > 1){
				foreach my $common_indvs (values %individual_seqs) {
					$common_cdr3_seqs->{$marker}{'all'}{$common_indvs}++;
				}
			}
		}

	}

# 		# Checks for common similar sequences between samples (that differ few substitutions)
# 		if (defined($substitutions)){
# 			my $similar_seqs;
# 			for (my $i=0; $i<=$#all_seqs; $i++){
# 				if ($i % 1000 == 0) {
# 					print "$i ";
# 				}
# 				if (defined($similar_seqs->{$i})){
# 					next;
# 				}
# 				# Finds similar seqs among all remaining seqs
# 				for (my $j=0; $j<=$#all_seqs; $j++){
# 					if (defined($similar_seqs->{$j}) || $i == $j){
# 						next;
# 					}
# 					my $diff = count_errors($all_seqs[$i],$all_seqs[$j],$substitutions+1);
# 					if (defined($diff) && $diff <= $substitutions) {
# 						push(@{$similar_seqs->{$i}}, $j);
# 						push(@{$similar_seqs->{$j}}, $i);
# 					}
# 				}
# 			}
# 			print "\n";
# 			my %checked_seqs;
# 			for (my $i=0; $i<=$#all_seqs; $i++){
# 				if (defined($checked_seqs{$i})){
# 					next;
# 				}
# 				$checked_seqs{$i} = 1;
# 				my $common_samples = 1;
# 				foreach my $sample (@{$dup_samples->{$individual}}) {
# 					foreach my $j (@{$similar_seqs->{$i}}){
# 						if (defined($checked_seqs{$j})){
# 							next;
# 						}
# 						my $seq = $all_seqs[$j];
# 						if (exists($sample_seqs->{$sample}{$marker}{$seq})){
# 							$common_samples++;
# 							$checked_seqs{$j}=1;
# 							last;
# 						}
# 					}
# 				}
# 				$common_cdr3_seqs->{$marker}{$individual}{'total_sim'}{$common_samples}++;
# 			}
# 		}
# 	}
	
	return $common_cdr3_seqs;
}

#################################################################################

# Retrieves a hash with the V and J segments assigned to each CDR3 sequence
sub read_tcr_segments {

	my ($headers,$seqs) = @_;

	my $tcr_segments;
	for (my $i=0; $i<=$#{$headers}; $i++) {
		my $header = $headers->[$i];
		my $cdr3_seq = $seqs->[$i];
		while ($header =~ /tc?r(.?)=(\S+)/gi){
			my $segment = 'tcr'.lc($1);
			my $segment_seq = $2;
			push(@{$tcr_segments->{$cdr3_seq}{$segment}}, $segment_seq)
		}
	}

	return $tcr_segments;

}

#################################################################################

# Retrieves a hash with the TCR regions as keys and the number of sequences containing them as values
sub count_tcr_regions {

	my ($headers,$seqs) = @_;

	my $count_regions;
	for (my $i=0; $i<=$#{$headers}; $i++) {
		my $header = $headers->[$i];
		my $seq_len = length($seqs->[$i]);
		my ($tcrv, @tcrv_pairs);
		while ($header =~ /tc?r(.?)=(\S+)/gi){
			if ($2 =~ /^[ACGT\-]+$/){
				next;
			}
			my $region = 'tcr'.lc($1);
			my $allele = $2;
# 			# Does not annotates ambiguos assignments
# 			if ($allele =~ /,/){
# 				next;
# 			}
			$count_regions->{$region}{$allele}++;
			$count_regions->{'lengths'}{$region}{$allele}{$seq_len}++;
			if ($region =~ /v$/){
				$tcrv = $allele;
			} else {
				push(@tcrv_pairs, $allele);
			}
		}
		if (defined($tcrv) && @tcrv_pairs){
			map $count_regions->{'pairs'}{$tcrv}{$_}++ , @tcrv_pairs;
		}
# 		if ($header =~ /(TCR[JV])\:([\d\-\*\,]+) \| (TCR[JV])\:([\d\-\*\,]+)/){
# 			$count_regions->{lc($1)}{$2}++;
# 			$count_regions->{lc($3)}{$4}++;
# 			$count_regions->{'lengths'}{lc($1)}{$2}{$seq_len}++;
# 			$count_regions->{'lengths'}{lc($3)}{$4}{$seq_len}++;
# 			if ($1 eq 'TCRV'){
# 				$count_regions->{'pairs'}{$2}{$4}++;
# 			} else {
# 				$count_regions->{'pairs'}{$4}{$2}++;
# 			}
# 		} elsif ($header =~ /(TCR[JV])\:([\d\-\*\,]+)/){
# 			$count_regions->{lc($1)}{$2}++;
# 			$count_regions->{'lengths'}{lc($1)}{$2}{$seq_len}++;
# 		}
# # # 		while ($header =~ /(TCR[:,\w\d]+)/g){
# # 		while ($header =~ /(TCR\w)\:([\d\,]+)/g){
# # 			foreach my $allele (split($2)){
# # 				$count_regions->{lc($1)}{$allele}++;
# # 			}
# # 		}
	}

	return $count_regions;

}

#################################################################################

# Retrieves a hash with depths as keys and the number of sequences associated with that depth as values
sub count_depths {

	my $headers = shift @_;

	my $count_depths = {'1'=>0};
# 	my $count_singletons = 0;
	foreach (@$headers) {
		if (/depth=(\d+)/){
			$count_depths->{$1}++;
# 			if ($1 == 1){
# 				$count_singletons++;
# 			}
		}
	}

	return $count_depths;

}

#################################################################################

sub assign_percentile_group {

	my ($value, $stat_data) = @_;

	if (!defined($value) || !defined($stat_data->percentile('20'))){
		return 0;
	}
	if ($value < $stat_data->percentile('20')) {
		return 0;
	} elsif ($value < $stat_data->percentile('40')) {
		return 1;
	} elsif ($value < $stat_data->percentile('60')) {
		return 2;
	} elsif ($value < $stat_data->percentile('80')) {
		return 3;
	} else {
		return 4;
	}
}

#################################################################################


1;

