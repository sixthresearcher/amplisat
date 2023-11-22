#!/usr/bin/perl -w
#
# Name: ampliHLA.pl
#
# Version: 1.0
#
# License: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
#          CC BY-NC-SA 4.0 (http://creativecommons.org/licenses/by-nc-sa/4.0/)
#
# Author: Alvaro Sebastian
#
# Support: Alvaro Sebastian (sixthresearcher@gmail.com)
#
# Sixth Researcher - www.sixthresearcher.com
#
# Description:
#   Performs HLA typing taking as input clustered variants created by AmpliSAS.
#
# Requires as input a FASTA or FASTQ file with sequences/reads and a CSV format file with primer/tag data information
# Example: perl ampliHLA.pl -d amplicon_data.csv -i reads.fq.gz -o results
#
# If the reads have been already demultiplexed into separate files (one file per sample), they can be packed into a single .zip or .tar.gz file and use it as input
# Example: perl ampliHLA.pl -i reads.tar.gz -o results
#
# Alternatively, variants already extracted can be given in an Excel format file obtained with AmpliSAS/AmpliCHECK.
# Example: perl ampliHLA.pl -i amplisas_results.xlsx -o results
#
# How to create small FASTQ example files from RNA-Seq or WXS big files with only the reads that map to HLA references:
#  1. Change in $DEFAULT_PARAMS 'allele_alignment' to 'bowtie -v3' to map only one reference as maximum per read an run AmpliHLA analysis:
#        ampliHLA.pl -thr 20 -ty mapping -r nuc -o amplihla_results -i SRR387401_1.fastq.gz SRR387401_2.fastq.gz
#  2. Use PICARD tools to generate 2 FASTQ files from the BAM file created by previous AmpliHLA analysis:
#        java -Xmx2g -jar /opt/picard-tools-1.109/SamToFastq.jar I=amplihla_results/SRR387401.hla_nuc.bam F=SRR387401_R1.fq F2=SRR387401_R2.fq
#  3. Compress the FASTQ files
#        gzip SRR387401_R1.fq SRR387401_R2.fq
#  4. Change in $DEFAULT_PARAMS 'allele_alignment' to the original value 'bowtie -a -v3', run AmpliHLA with the new FASTQ files and compare the results:
#        ampliHLA.pl -thr 20 -ty mapping -r nuc -o amplihla_results -i SRR387401_R1.fq.gz SRR387401_R2.fq.gz
 


my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Performs HLA typing by genotyping amplicon sequencing data with AmpliSAS clustering algorithm and assigning variants to alleles from the IMGT/HLA database.";


# Modules are in folder 'lib' in the path of the script
use lib "lib";
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;
use Bio::Onco;
use Spreadsheet::XLSX;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Sort::Naturally;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $SCRIPT_DIR = dirname(__FILE__);

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options
# Analysis type
my $INP_analysis_type = 'amplisas';
# Data type
my $INP_reference_type = 'all'; # 'nuc' 'gen'
# HLA alleles to genotype by mapping
my @HLA_LOCI = ('A','B','C','DQA1','DQB1','DRB1');
# HLA sequence patterns to extract variable regions: exon2+exon3 (class I) or exon 2 (class II)
my $HLA_ALLELE_PATTERNS = {
# 	'A'    => 'ATGGCCGTCATGGCGCCCCG.+GAG[G|A]TGCTGGGCCCTG[G|A]GCT',
# 	'B'    => 'ATGCTGGTCATGGCGCCCCG.+GAGGTGCTGGGCCCTGGGTT',
# 	'C'    => 'ATGCGGGTCATGGCGCCCCG.+GAGGTGCTGGGCCCTGGGCT',
# 	'DQA1' => 'TAAACAAAGCTCTGCTGCTG.+CA[C|T]CCT[C|T]AT[C|T]TGTCTTGTGG',
# 	'DQB1' => 'GCCTTCGGGCAGCAACTGTG.+CCTGCTG[G|A]TCTGCTC[G|A]GTGA',
# 	'DRB1' => 'GCTCCTGCATGACAGCGCTG.+CCT[C|G]CTGGTCTG[T|C]TCTGTGA',
	'A'    => { '5'=>['GCTC[C|T]CACTCCA\w{11}','\w{12}TGAGGTATTTC'] , '3'=>['GGAGACGCTG\w{10}','\w{10}CAGCGCACGG'] },
	'B'    => { '5'=>['GCTCCCACTCCA\w{11}','\w{12}TGAGGTATTTC'] , '3'=>['ACCTGGAGA\w{10}','\w{9}ACGGGAAGGA'] },
	'C'    => { '5'=>['GCTCCCACTCC\w{12}','\w{11}ATGA[G|A]GTATTTC'] , '3'=>['GGGAAG[G|A]AGACG\w{11}','\w{12}CTGCAGCGCGC'] },
	'DQA1' => { '5'=>['GGTGTAAACT\w{9}','\w{10}TGTACCAGT'] , '3'=>['AGCCTCTTCT\w{11}','\w{10}GAAACACTGGG'] },
	'DQB1' => { '5'=>['AGGATTTCGTG[T|C]\w{12}','\w{12}[A|T]CCAGTTTAAGG'] , '3'=>['CCGCACGACCTT\w{11}','CCGCGGGATC[C|T]T\w{11}','\w{12}GCAGCGGCGAG','\w{12}GCAGAGGAGAG'] },
	'DRB1' => { '5'=>['TGT[C|T]ATTTCTTCAA\w{14}','\w{14}[T|C]GGGACGGAGCGGG'] , '3'=>['GAGAGCTTCAC[A|G]\w{13}','\w{12}GTGCAG[C|A]GGCGAG'] },
};
# Common human alleles to discard rare ones
my $HLA_FREQUENT_ALLELES = {
  'A'    => ['A*01:01','A*01:02','A*01:03','A*01:06','A*02:01','A*02:02','A*02:03','A*02:04','A*02:05','A*02:06','A*02:07','A*02:08','A*02:10','A*02:11','A*02:12','A*02:13','A*02:14','A*02:16','A*02:17','A*02:19','A*02:20','A*02:22','A*02:24','A*02:25','A*02:26','A*02:28','A*02:33','A*02:34','A*02:36','A*02:40','A*02:44','A*02:45','A*02:46','A*03:01','A*03:02','A*11:01','A*11:02','A*11:04','A*11:06','A*23:01','A*23:03','A*23:04','A*23:05','A*24:02','A*24:03','A*24:04','A*24:05','A*24:06','A*24:07','A*24:08','A*24:09N','A*24:10','A*24:13','A*24:14','A*24:15','A*24:17','A*24:18','A*24:21','A*24:22','A*24:23','A*25:01','A*25:02','A*26:01','A*26:02','A*26:03','A*26:05','A*26:06','A*26:08','A*26:09','A*26:11N','A*26:12','A*26:18','A*29:01','A*29:02','A*29:03','A*30:01','A*30:02','A*30:04','A*30:07','A*30:08','A*30:10','A*31:01','A*31:02','A*31:03','A*31:04','A*31:05','A*32:01','A*32:02','A*32:03','A*33:01','A*33:03','A*34:01','A*34:02','A*36:01','A*43:01','A*66:01','A*66:02','A*66:03','A*68:01','A*68:02','A*68:03','A*68:04','A*68:05','A*68:06','A*68:12','A*68:15','A*69:01','A*74:01','A*74:03','A*80:01'],
  'B'    => ['B*07:02','B*07:04','B*07:05','B*07:06','B*07:07','B*07:12','B*07:17','B*07:20','B*07:26','B*08:01','B*08:03','B*08:04','B*08:05','B*13:01','B*13:02','B*13:03','B*13:04','B*14:01','B*14:02','B*14:03','B*14:04','B*14:05','B*14:06','B*15:01','B*15:02','B*15:03','B*15:04','B*15:05','B*15:06','B*15:07','B*15:08','B*15:09','B*15:10','B*15:11','B*15:12','B*15:13','B*15:15','B*15:16','B*15:17','B*15:18','B*15:20','B*15:21','B*15:23','B*15:24','B*15:25','B*15:27','B*15:28','B*15:29','B*15:30','B*15:31','B*15:32','B*15:34','B*15:35','B*15:37','B*15:38','B*15:39','B*15:40','B*15:42','B*15:55','B*15:58','B*15:64','B*15:67','B*18:01','B*18:02','B*18:03','B*18:07','B*18:12','B*27:02','B*27:03','B*27:04','B*27:05','B*27:06','B*27:07','B*27:08','B*27:09','B*27:11','B*27:14','B*35:01','B*35:02','B*35:03','B*35:04','B*35:05','B*35:06','B*35:08','B*35:09','B*35:10','B*35:11','B*35:12','B*35:13','B*35:14','B*35:16','B*35:17','B*35:20','B*35:21','B*35:28','B*35:34','B*37:01','B*37:02','B*38:01','B*38:02','B*39:01','B*39:02','B*39:03','B*39:05','B*39:06','B*39:07','B*39:08','B*39:09','B*39:10','B*39:11','B*39:13','B*39:15','B*40:01','B*40:02','B*40:03','B*40:04','B*40:05','B*40:06','B*40:07','B*40:08','B*40:10','B*40:11','B*40:12','B*40:14','B*40:15','B*40:16','B*40:19','B*41:01','B*41:02','B*41:03','B*42:01','B*42:02','B*44:02','B*44:03','B*44:04','B*44:05','B*44:06','B*44:07','B*44:08','B*44:10','B*44:15','B*44:26','B*45:01','B*45:02','B*46:01','B*47:01','B*47:02','B*47:03','B*48:01','B*48:02','B*48:03','B*48:04','B*49:01','B*50:01','B*50:02','B*51:01','B*51:02','B*51:04','B*51:05','B*51:06','B*51:07','B*51:08','B*51:09','B*51:10','B*51:12','B*51:14','B*51:18','B*51:27N','B*52:01','B*52:02','B*53:01','B*53:02','B*53:03','B*53:05','B*54:01','B*55:01','B*55:02','B*55:03','B*55:04','B*55:07','B*56:01','B*56:02','B*56:03','B*56:04','B*56:05','B*57:01','B*57:02','B*57:03','B*57:04','B*58:01','B*58:02','B*59:01','B*67:01','B*73:01','B*78:01','B*78:02','B*78:03','B*81:01','B*82:01','B*82:02'],
  'C'    => ['C*01:02','C*01:03','C*01:04','C*02:02','C*02:03','C*03:02','C*03:03','C*03:04','C*03:05','C*03:06','C*03:07','C*03:08','C*03:10','C*04:01','C*04:03','C*04:04','C*04:06','C*04:07','C*05:01','C*05:04','C*06:02','C*06:04','C*07:01','C*07:02','C*07:03','C*07:04','C*07:05','C*07:06','C*07:07','C*07:08','C*07:13','C*07:14','C*08:01','C*08:02','C*08:03','C*08:04','C*08:05','C*08:06','C*12:02','C*12:03','C*12:04','C*12:05','C*12:07','C*14:02','C*14:03','C*15:02','C*15:03','C*15:04','C*15:05','C*15:07','C*15:08','C*15:09','C*16:01','C*16:02','C*16:04','C*17:01','C*18:01'],
  'DQA1' => ['DQA1*01:01','DQA1*01:02','DQA1*01:03','DQA1*01:04','DQA1*02:01','DQA1*03:01','DQA1*03:02','DQA1*03:03','DQA1*04:01','DQA1*05:01','DQA1*05:02','DQA1*06:01'],
  'DQB1' => ['DQB1*02:01','DQB1*02:02','DQB1*02:03','DQB1*03:01','DQB1*03:02','DQB1*03:03','DQB1*03:04','DQB1*03:05','DQB1*03:09','DQB1*04:01','DQB1*04:02','DQB1*05:01','DQB1*05:02','DQB1*05:03','DQB1*05:04','DQB1*06:01','DQB1*06:02','DQB1*06:03','DQB1*06:04','DQB1*06:05','DQB1*06:06','DQB1*06:07','DQB1*06:08','DQB1*06:09','DQB1*06:11','DQB1*06:13','DQB1*06:15'],
  'DRB1' => ['DRB1*01:01','DRB1*01:02','DRB1*01:03','DRB1*03:01','DRB1*03:02','DRB1*03:03','DRB1*03:05','DRB1*03:08','DRB1*03:17','DRB1*04:01','DRB1*04:02','DRB1*04:03','DRB1*04:04','DRB1*04:05','DRB1*04:06','DRB1*04:07','DRB1*04:08','DRB1*04:10','DRB1*04:11','DRB1*04:12','DRB1*04:13','DRB1*04:15','DRB1*04:16','DRB1*04:36','DRB1*07:01','DRB1*08:01','DRB1*08:02','DRB1*08:03','DRB1*08:04','DRB1*08:05','DRB1*08:06','DRB1*08:07','DRB1*08:08','DRB1*08:09','DRB1*08:11','DRB1*08:18','DRB1*09:01','DRB1*10:01','DRB1*11:01','DRB1*11:02','DRB1*11:03','DRB1*11:04','DRB1*11:08','DRB1*11:09','DRB1*11:11','DRB1*11:13','DRB1*11:14','DRB1*11:19','DRB1*11:30','DRB1*12:01','DRB1*12:02','DRB1*12:03','DRB1*12:04','DRB1*12:05','DRB1*13:01','DRB1*13:02','DRB1*13:03','DRB1*13:04','DRB1*13:05','DRB1*13:06','DRB1*13:07','DRB1*13:09','DRB1*13:10','DRB1*13:12','DRB1*13:17','DRB1*13:20','DRB1*13:23','DRB1*13:25','DRB1*13:27','DRB1*13:31','DRB1*14:01','DRB1*14:02','DRB1*14:03','DRB1*14:04','DRB1*14:05','DRB1*14:06','DRB1*14:07','DRB1*14:08','DRB1*14:09','DRB1*14:10','DRB1*14:12','DRB1*14:13','DRB1*14:15','DRB1*14:17','DRB1*14:18','DRB1*14:19','DRB1*14:21','DRB1*14:43','DRB1*15:01','DRB1*15:02','DRB1*15:03','DRB1*15:04','DRB1*15:05','DRB1*15:06','DRB1*16:01','DRB1*16:02','DRB1*16:03','DRB1*16:05'],
};
# Default analysis parameters
my $DEFAULT_PARAMS = {
	# Genomic+CDS full sequences
	'all' => {
		# Allele location folder
		'allele_path' => $SCRIPT_DIR.'/imgt',
		# Allele file - GENOTYPING MAY FAIL BECAUSE EVERY ALLELE HAS VARIABLE NUMBER OF ENTRIES WITH DIFFERENT LENGTHS AND IT SHOULD BE COUNTERBALANCED
		'allele_file' => $SCRIPT_DIR.'/imgt/hla_all.fa',
	},
	# Genomic trimmed sequences (variable region)
	'gen' => {
		# Allele location folder
		'allele_path' => $SCRIPT_DIR.'/imgt',
		# Allele file
		'allele_file' => $SCRIPT_DIR.'/imgt/hla_gen.fa',
	},
	# CDS trimmed sequences (variable region)
	'nuc' => {
		# Allele location folder
		'allele_path' => $SCRIPT_DIR.'/imgt',
		# Allele file
		'allele_file' => $SCRIPT_DIR.'/imgt/seq2hla.fa', # BOEGEL ET AL. CURATED VARIABLE REGION SEQUENCES (LOWER NUMBER OF MAPPED READS BUT MORE ACCURATE)
		# 'allele_file' => $SCRIPT_DIR.'/imgt/hla_nuc.fa', # SHORTER THAN BOEGEL ET AL. SEQUENCES (ARTIFICIALLY EXTENDED)
	},
	'amplisas' => {
		# Allele location folder
		# 'allele_path' => $SCRIPT_DIR.'/imgt',
		# Allele matching parameters
		'allele_alignment' => 'dna blastn -evalue 1E-5 -ungapped',
		# Minimum % of sequence aligned to a reference
		'min_allele_align' => 100,
		# Minimum % of identity in the portion aligned to a reference
		'min_allele_ident' => 100,
		# Amplicons with lower total coverage will be discarded
		'min_amplicon_depth' => 100,
		# Minimum frequency of variants after clustering
		# 'min_amplicon_seq_frequency' => 'auto',
		# Variants with equal or higher identity will be clustered together
		'identity_threshold' => undef,
		# Keep singletons after clustering
		'keep_singletons' => undef,
	},
	'mapping' => {
		# Allele matching parameters
		'allele_alignment' => 'bowtie -a -v3', # ORIGINAL AND THE BEST, ALLOWS 3 MISMATCHES (-v3) AND EACH READ IS MAPS AS MANY REFERENCES AS POSSIBLE
		# 'allele_alignment' => 'bowtie -v3', # USEFUL TO CREATE SMALL EXAMPLE FILES: EACH READ MAPS ONE REFERENCE AS MAXIMUM (READ INSTRUCTIONS AT THE TOP OF THE SCRIPT)
		# 'allele_alignment' => 'bowtie2 --end-to-end -k 2', # VERY BAD, MAPS VERY FEW READS
		# 'allele_alignment' => 'bowtie2 --sensitive-local', # -k 10', # FAILS BECAUSE READS MAP LOCALLY/PARTIALLY MANY DIFFERENT ALLELES, NOT ALLWAYS THE TRUE ONES
		# Max. number of allowed substitutions in mapping
		'max_substitutions' => 1,
		# Max. number of allowed indels in mapping
		'max_indels' => 0,
		# Keep or discard rare alleles
		'rare_alleles' => 0,
		# Minimum frequency percentage to annotate/discard low frequency alleles
# 		'min_allele_frequency' => 10, # 10%
		# Minimum coverage to annotate/discard low coverage alleles
# 		'min_allele_coverage' => 10, # 10
		# Minimum distance in bps to the start/end of the read to annotate the alleles
		'min_allele_position' => 2,

# 		# Minimum % of sequence aligned to a reference
# 		'min_allele_align' => 100,
# 		# Minimum % of identity in the portion aligned to a reference
# 		'min_allele_ident' => 100,
# 		# Amplicons with lower total coverage will be discarded
# 		'min_amplicon_depth' => 100,
# 		# Minimum frequency of variants after clustering
# 		# 'min_amplicon_seq_frequency' => 'auto',
# 		# Variants with equal or higher identity will be clustered together
# 		'identity_threshold' => undef,
# 		# Keep singletons after clustering
# 		'keep_singletons' => undef,
	},
};

# Input parameters variables (will be filled with command line options or default values)
my %INP_params;

my (@INP_reads_files,$INP_amplicons_file,$INP_params_file,$INP_outpath,$INP_allele_file,$INP_shuffle,$INP_tech,$INP_nreads,$INP_nreads_amplicon,$INP_threads,$INP_zip,$INP_verbose,$INP_update,$INP_test);

GetOptions(
	'h|help|?' =>  \&usage,
	'ty|type=s' => \$INP_analysis_type,
	'r|refs=s' => \$INP_reference_type,
	'i|input=s{,}' => \@INP_reads_files,
	'd|data=s' => \$INP_amplicons_file,
	'p|params=s' => \$INP_params_file,
	'o|output=s' => \$INP_outpath,
	'a|alleles=s' => \$INP_allele_file,
	't|tech=s' => \$INP_tech,
	's|shuffle' => \$INP_shuffle,
	'n|number=i' => \$INP_nreads,
	'na|max=i' => \$INP_nreads_amplicon,
	'al|aligned=f' => \$INP_params{'min_allele_align'},
	'id|ident=f' => \$INP_params{'min_allele_ident'},
	'fr|freq=f' => \$INP_params{'min_amplicon_seq_frequency'},
	'ci|clid=f' => \$INP_params{'identity_threshold'},
	'min=i' => \$INP_params{'min_amplicon_depth'},
	'ks|keepsingle' => \$INP_params{'keep_singletons'},
	'thr|threads=i' => \$INP_threads,
	'z|zip' => \$INP_zip,
	'v|verbose' => \$INP_verbose,
	'u|update' => \$INP_update,
	'test' => \$INP_test,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\nGlobal options:\n";
	print "  -ty <type>\tTyping strategy (default='$INP_analysis_type', 'mapping').\n";
	print "  -r <type>\tReference type (default='$INP_reference_type', 'nuc', 'gen').\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -v\t\tPrint AmpliSAS output and verbose FASTA files.\n";
	print "  -u\t\tUpdates the HLA alleles reference file to the latest version from the IMGT/HLA database.\n";
	print "  -h\t\tHelp.\n";
	print "\nAmpliSAS clustering options:\n";
	print "  -i <file>\tInput FASTQ or FASTA file (compressed or uncompressed) or set of files packed into a unique .ZIP or .TAR.GZ file or AmpliSAS format Excel file.\n";
	print "  -d <file>\tCSV file with primer/amplicon data.\n";
	print "  -o <path>\tOutput folder name.\n";
	print "  -a <file>\tFASTA file with allele reference sequences and names (using WHO nomenclature).\n";
	print "  -t <tech>\tUse recommended technology parameters ('Illumina', 'IonTorrent', '454', 'Unknown').\n";
	print "  -s\t\tShuffle/randomize reads/sequences to analyze.\n";
	print "  -n <number>\tNumber of reads/sequences to analyze.\n";
	print "  -na <number>\tNumber of reads/sequences per amplicon to analyze.\n";
	if (defined($DEFAULT_PARAMS->{$INP_analysis_type}{'min_allele_align'})){
		print "  -al <perc>\tMinimum sequence percentage required to align to a reference to be annotated (default=".$DEFAULT_PARAMS->{$INP_analysis_type}{'min_allele_align'}."%).\n";
	}
	if (defined($DEFAULT_PARAMS->{$INP_analysis_type}{'min_allele_ident'})){
		print "  -id <perc>\tMinimum percentage of the aligned fragment required to be identical to a reference to be annotated (default=".$DEFAULT_PARAMS->{$INP_analysis_type}{'min_allele_ident'}."%).\n";
	}
	print "  -fr <freq>\tFilter sequences with lower frequency after clustering.\n";
	print "  -ci <number>\tCluster together sequences with higher or equal identity.\n";
	if (defined($DEFAULT_PARAMS->{$INP_analysis_type}{'min_amplicon_depth'})){
		print "  -min <number>\tAmplicons with lower total depth/coverage will be discarded (default=".$DEFAULT_PARAMS->{$INP_analysis_type}{'min_amplicon_depth'}.").\n";
	}
	if (defined($DEFAULT_PARAMS->{$INP_analysis_type}{'keep_singletons'})){
		print "  -ks\t\tKeep singletons after clustering (default=".$DEFAULT_PARAMS->{$INP_analysis_type}{'keep_singletons'}.").\n";
	}
	print "\nMapping options:\n";
	print "  -i <file1> <file2>\n\t\tInput one single-end or two paired-end FASTQ or FASTA files (compressed or uncompressed).\n";
	print "  -p <file>\tCSV file with params data.\n";
	print "\n";
	exit;
}

# Updates FASTA file with HLA allele reference sequences from IMGT/HLA database
if (defined($INP_update)){

	print "\nRunning '$COMMAND_LINE'\n";
	
	print "\nUpdating HLA alleles reference file to the latest version from the IMGT/HLA database.\n";

	# Downloads HLA alleles from IMGT/HLA database (genomic and CDS sequences)
	my (@allele_files, @hla_all_seqs, @hla_all_headers);
	foreach my $data_type (('gen','nuc')){
		my $alleles_file = $DEFAULT_PARAMS->{$data_type}{'allele_file'};
		if ($alleles_file =~ /seq2hla/){
			push(@allele_files,$alleles_file);
			printf("\nReading '%s' HLA alleles reference file.\n",$alleles_file);
			printf("\tCurrent version contains %d allele sequences.\n", scalar count_seqs_from_fasta($alleles_file));
			next;
		}
		my $alleles_dir = $DEFAULT_PARAMS->{$data_type}{'allele_path'};
		if (!-d $alleles_dir){
			mkdir($alleles_dir);
		}
		my $old_seqs_count;
		if (-e $alleles_file) {
			$old_seqs_count = count_seqs_from_fasta($alleles_file);
		}
		printf("\nUpdating '%s' HLA alleles reference file.\n",$alleles_file);
		# Downloads and extracts sequences from exon2+exon3 (class I) or exon 2 (class II) and 75 additional nts like in Boegel et al.
		my (@hla_trimmed_seqs, @hla_trimmed_headers);
		foreach my $hla_locus (@HLA_LOCI){
			printf("\tProcessing '%s' locus.\n",$hla_locus);
			open(FASTA_DATA,"wget -q -O - ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/$hla_locus\_$data_type.fasta |") || die "# $0 : HLA alleles for '$hla_locus' locus cannot be downloaded from IMGT/HLA database server.\n\n";
			# open(FASTA_DATA,"cat ~/evobio/scripts/imgt/loci_files/$hla_locus\_$data_type.fa |") || die "# $0 : HLA alleles for '$hla_locus' locus cannot be downloaded from IMGT/HLA database server.\n\n";
			my ($hla_locus_seqs,$hla_locus_headers) = read_fasta([<FASTA_DATA>]);
			close(FASTA_DATA);
			push(@hla_all_seqs,@$hla_locus_seqs);
			push(@hla_all_headers,@$hla_locus_headers);
			my $count = 0;
			my (@hla_trimmed_seqs_, @hla_trimmed_headers_);
			# Trims the allele sequences according to the patterns
			for (my $i=0; $i<=$#{$hla_locus_headers}; $i++) {
				my $match = 0;
				foreach my $pat5 (@{$HLA_ALLELE_PATTERNS->{$hla_locus}{'5'}}){
					foreach my $pat3 (@{$HLA_ALLELE_PATTERNS->{$hla_locus}{'3'}}){
						my $pattern = $pat5.'.+'.$pat3;
						if ($hla_locus_seqs->[$i]=~ /$pattern/) {
							push(@hla_trimmed_seqs_,$&);
							push(@hla_trimmed_headers_,$hla_locus_headers->[$i]);
							$match = 1;
							$count++;
							last;
						}
					}
					if ($match){
						last;
					}
				}
			}
			# Removes redundant sequences
			my $seq_clusters = cluster_identical_seqs(\@hla_trimmed_seqs_,\@hla_trimmed_headers_);
			for (my $i=0; $i<=$#hla_trimmed_seqs_; $i++){
				if (defined($seq_clusters->{$hla_trimmed_seqs_[$i]})){
					my (%alleles,$seq_ids);
					foreach my $header (@{$seq_clusters->{$hla_trimmed_seqs_[$i]}}){
						if ($header =~ /HLA:(HLA\d+) (\w+)\*(\d+)\:(\d+)\:?([\d\:]+[NLSCAQ]?)?/) { # Ex. HLA:HLA00132 B*07:02:01N 3323 bp
							$alleles{"$2*$3:$4"}++;
							push(@{$seq_ids->{"$2*$3:$4"}},$1);
						}
					}
					my @sorted_alelles = sort {$alleles{$b}<=>$alleles{$a}} keys %alleles;
					my $major_allele = $sorted_alelles[0];
					if ($alleles{$major_allele}*2 > scalar keys %alleles) {
						push(@hla_trimmed_seqs,$hla_trimmed_seqs_[$i]);
						push(@hla_trimmed_headers,sprintf("%s %s", $seq_ids->{$major_allele}[0], $major_allele));
					} else {
						foreach my $allele (@sorted_alelles){
							push(@hla_trimmed_seqs,$hla_trimmed_seqs_[$i]);
							push(@hla_trimmed_headers,sprintf("%s %s", $seq_ids->{$allele}[0], $allele));
						}
					}
					delete($seq_clusters->{$hla_trimmed_seqs_[$i]});
				}
			}

			if (!$count) {
				printf("\nERROR: HLA alleles for '%s' locus cannot be found in the IMGT/HLA data.\n\n",$hla_locus);
				exit;
			} else {
				printf("\tSelected %d allele sequences from %d.\n",$count,scalar @$hla_locus_headers);
			}
		}
		push(@allele_files,create_fasta_file(\@hla_trimmed_seqs,\@hla_trimmed_headers,$alleles_file));
		#printf("\nUpdated '%s' HLA alleles reference file.\n",$alleles_file);
		if (defined($old_seqs_count)){
			printf("\tPrevious version contained %d allele sequences.\n", $old_seqs_count);
		}
		printf("\tNew version contains %d allele sequences.\n", scalar count_seqs_from_fasta($alleles_file));
	}
	push(@allele_files,create_fasta_file(\@hla_all_seqs,\@hla_all_headers,$DEFAULT_PARAMS->{'all'}{'allele_file'}));
	printf("\nUpdating '%s' HLA alleles reference file.\n",$DEFAULT_PARAMS->{'all'}{'allele_file'});
	print "\nCreating BOWTIE INDEX files.\n";
	foreach my $allele_file (@allele_files){
		`bowtie-build '$allele_file' '$allele_file' 2>&1`;
		`bowtie2-build '$allele_file' '$allele_file' 2>&1`;
	}
	print "\nCreated BOWTIE and BOWTIE2 INDEX files.\n";
	print "\n";

	exit;

}

# Checks if data type is correct
$INP_reference_type = lc($INP_reference_type);
if (!defined($DEFAULT_PARAMS->{$INP_reference_type})){
	print "\nERROR: '$INP_reference_type' data type is not supported ('all', 'nuc', 'gen').\n\n";
	usage();
	exit;
}
# Checks if analysis type is correct
$INP_analysis_type = lc($INP_analysis_type);
if (!defined($DEFAULT_PARAMS->{$INP_analysis_type})){
	print "\nERROR: '$INP_analysis_type' typing strategy is not supported ('amplisas', 'mapping').\n\n";
	usage();
	exit;
}
# Retrieves the location of the allele database
my $INP_allele_dir = $DEFAULT_PARAMS->{$INP_reference_type}{'allele_path'};
if (!defined($INP_allele_file)){
	$INP_allele_file = $DEFAULT_PARAMS->{$INP_reference_type}{'allele_file'};
# 	print "\nERROR: FASTA file with HLA alleles cannot be found.\n\n";
# 	usage();
# 	exit;
}
# Checks if exists alleles file
if (!-e $INP_allele_file){
	print "\nERROR: '$INP_allele_file' FASTA file with HLA alleles cannot be found.\n\n";
	usage();
	exit;
}


# Downloads input files if links are provided
my @INP_file_urls;
if (@INP_reads_files && $INP_reads_files[0] =~ /(ftp|http|https)\:/ ){
	for (my $i=0; $i<=$#INP_reads_files; $i++) {
		$INP_file_urls[$i] = $INP_reads_files[$i];
		printf("\nDownloading file from '%s'\n",$INP_file_urls[$i]);
		$INP_reads_files[$i] = download_file($INP_file_urls[$i]);
		if (!defined($INP_reads_files[$i])){
			printf("\nERROR: broken link: '%s'.\n\n",$INP_file_urls[$i]);
			exit;
		}
	}
}

my ($INP_excel_file, $INP_multifile);
if ($INP_analysis_type eq 'amplisas'){
	# Checks if an excel file or multifle are given as input with sequences previously processed
	if (defined($INP_reads_files[0]) && is_xlsx($INP_reads_files[0])){
		$INP_excel_file = $INP_reads_files[0];
	# Checks if a set of demultiplexed files is given as input into a compressed file
	} elsif (defined($INP_reads_files[0]) && is_multifile($INP_reads_files[0])){
		$INP_multifile = $INP_reads_files[0];
	# Prints usage help if no input file is specified
	} elsif (!defined($INP_reads_files[0])){
		print "\nERROR: You must specify sequence input file.\n\n";
		usage();
		exit;
	} elsif (!defined($INP_amplicons_file)){
		print "\nERROR: You must specify amplicon data input file.\n\n";
		usage();
		exit;
	}
} elsif ($INP_analysis_type eq 'mapping') {
	# Checks if an excel file or multifle are given as input with sequences previously processed
	if (!@INP_reads_files){
		print "\nERROR: You must specify at least one input file.\n\n";
		usage();
		exit;
	}
} else {
	print "\nERROR: Please specify a valid typing strategy.\n\n";
	usage();
	exit;

}

# Default output path
if (!defined($INP_outpath)){
	$INP_outpath  = lc((split('\.',$SCRIPT_NAME))[0]);
}

print "\nRunning '$COMMAND_LINE'\n";

# Creates output folder
if (!-d $INP_outpath){
	mkdir($INP_outpath);
}

my ($amplisas_results, $md5_to_name);
if ($INP_analysis_type eq 'amplisas'){

	# AMPLIHLA ONLY ALLOWS GLOBAL PARAMETERS FOR ALL THE MARKERS ('all')
	# Reads parameters data from CVS input file
	my %paramsdata;
	if (defined($INP_amplicons_file)){
		my $paramsdata_ = read_amplicon_data($INP_amplicons_file, 'params');
		# Stores the params in a simple hash
		map $paramsdata{$_} = $paramsdata_->{$_}{'all'}[0] , keys %$paramsdata_;
	}
	# Command line params have priority over CSV file and auto mode
	foreach my $param (keys %{$DEFAULT_PARAMS->{$INP_analysis_type}}) {
		if (defined($INP_params{$param})){
			$paramsdata{$param} = lc($INP_params{$param});
		} elsif (!defined($paramsdata{$param}) && defined($DEFAULT_PARAMS->{$INP_analysis_type}{$param})){
			$paramsdata{$param} = lc($DEFAULT_PARAMS->{$INP_analysis_type}{$param});
		}
	}
	# Checks percentage threshold parameters, if defined
	foreach my $param (('min_allele_align','min_allele_ident','min_amplicon_seq_frequency','identity_threshold')){
		if (defined($paramsdata{$param})){
			my $value = $paramsdata{$param};
			if (defined($value) && $value !~ /auto/){
				if ($value =~ /([\d\.]+)%/) {
					$value = $1;
				}
				if (!is_numeric($value) || $value<0 || $value>100){
					printf("\nERROR: '%s' value must be a number between 0 and 100.\n\n",$param);
					exit;
				} else {
					$paramsdata{$param} = $value;
				}
			}
		}
	}
	# Checks and reads alleles file
	my ($alleledata, $alleledata_);
	if (defined($INP_allele_file) && -e $INP_allele_file){
		print "\nReading HLA allele sequences from '$INP_allele_file'.\n";
		$alleledata_ = read_allele_file($INP_allele_file,['no warnings']);
		# Replaces allele names
		foreach my $allele_name (keys %{$alleledata_}){
			$allele_name =~ /HLA:HLA\d+\s(.+?)\s/; # Ex. HLA:HLA00132 B*07:02:01 3323 bp
			$alleledata->{$1} = $alleledata_->{$allele_name};
		}
	}
	# Runs AmpliSAS in auto mode before HLA typing if reads are given as input file
	# my $amplisas_results;
	if (!defined($INP_excel_file)){
		print "\nCalling AmpliSAS for sequence de-multiplexing, clustering and filtering.\n";
		my $amplisas_options = '';
		if (defined($INP_tech)){
			$amplisas_options .= " -t $INP_tech";
		}
		if (defined($INP_shuffle)){
			$amplisas_options .= " -s";
		}
		if (defined($INP_nreads)){
			$amplisas_options .= " -n $INP_nreads";
		}
		if (defined($INP_nreads_amplicon)){
			$amplisas_options .= " -na $INP_nreads_amplicon";
		}
		if (defined($paramsdata{'min_amplicon_seq_frequency'})){
			$amplisas_options .= " -fr ".$paramsdata{'min_amplicon_seq_frequency'};
		}
		if (defined($paramsdata{'identity_threshold'})){
			$amplisas_options .= " -ci ".$paramsdata{'identity_threshold'};
		} elsif (!defined($INP_tech)) {
			$amplisas_options .= " -auto";
		}
		if (defined($paramsdata{'keep_singletons'})){
			$amplisas_options .= " -ks";
		}
		if (defined($INP_verbose)) {
			$amplisas_options .= " -v";
		}
		if (defined($INP_threads)){
			$amplisas_options .= " -thr $INP_threads";
		}
		my ($amplisas_command, $amplisas_output);
		if (defined($INP_amplicons_file)){
			$amplisas_command = "$SCRIPT_DIR/ampliSAS.pl -denovo -i $INP_reads_files[0] -d $INP_amplicons_file -o $INP_outpath $amplisas_options";
		} else {
			$amplisas_command = "$SCRIPT_DIR/ampliSAS.pl -denovo -i $INP_reads_files[0] -o $INP_outpath $amplisas_options";
		}
	# 	print "\n$amplisas_command\n";
	# 	exit;
		# Prints live output during AmpliSAS running
		open CMD,'-|',$amplisas_command or die $@;
		while (my $line=<CMD>) {
			print $line;
			$amplisas_output .= $line;
		}
		close CMD;
		if ($amplisas_output =~ /results (stored|written) into '(.+)'/ && $2 eq $INP_outpath) {
			printf("\nReading AmpliSAS results.\n");
			$amplisas_results = read_amplisas_file_results("$INP_outpath/results.xlsx");
		} else {
			`rm -rf $INP_outpath`;
			if ($amplisas_output =~ /(ERROR:?\s*.+)/ ) {
				print "\n$1\n\n";
			} else {
	# 			if ($amplisas_output !~ /\nUsage:/){
	# 				printf("\nAmpliSAS output:\n\t%s\n",join("\n\t",split("\n",$amplisas_output)));
	# 			}
				print "\nERROR: AmpliSAS failed analyzing the data, try to run AmpliSAS independently and use the Excel results file as AmpliHLA input.\n\n";
			}
			exit;
		}
	} else { # Reads AmpliSAS results file
		printf("\nReading file '%s'.\n", $INP_excel_file);
		$amplisas_results = read_amplisas_file_results($INP_excel_file);
		# Creates output folder
		if (!-d $INP_outpath){
			mkdir($INP_outpath);
		}
	}

	# Rename variants with HLA allele names
	my %md5_to_sequence;
	foreach my $marker_name (keys %{$amplisas_results}){
		if (defined($amplisas_results->{$marker_name})){
			foreach my $md5 (keys %{$amplisas_results->{$marker_name}{'seq_data'}}){
				$md5_to_sequence{$md5} = $amplisas_results->{$marker_name}{'seq_data'}{$md5}{'sequence'};
			}
		}
	}
	print "\nMatching allele sequences.\n";
	my $allele_align_params = {
		'alignment' => $paramsdata{'allele_alignment'},
		'aligned' => sprintf("%.2f",$paramsdata{'min_allele_align'}/100),
		'ident' => sprintf("%.2f",$paramsdata{'min_allele_ident'}/100),
	};
	$md5_to_name = match_alleles($alleledata,\%md5_to_sequence,undef,$allele_align_params,$INP_threads);
	if (!defined($md5_to_name)){
		print "\nERROR: It was an error matching reference alleles.\n\n";
		exit;
	}

	# Perform genotyping by unifying all marker assignations

	# First discover which HLA locus is associated to each marker (A, B, DQ, DR...)
	print "\nAssigning HLA loci to markers.\n";
	my ($hla_loci, $samples, $sample_marker_alleles);
	foreach my $marker_name (sort {$a cmp $b} keys %{$amplisas_results}){
		if (defined($amplisas_results->{$marker_name})){
			my %count_hla_allele_loci;
			# Checks only the highest coverage variant of the marker (to avoid confusions)
			my $md5 = $amplisas_results->{$marker_name}{'seq_md5s'}[0];
			if (!defined($md5_to_name->{$md5})){
				next;
			}
			foreach my $hla_allele (split(' \| ',$md5_to_name->{$md5})) {
# 				if ($hla_allele =~ /HLA:HLA\d+ (\w+)\*.+/) { # Ex. HLA:HLA00132 B*07:02:01 3323 bp, exclude non expressed alleles or rare ones ex. HLA:HLA02369 A*03:21N 3095 bp
				if ($hla_allele =~ /^(\w+)\*/) { # Ex. HLA:HLA00132 B*07:02:01 3323 bp, exclude non expressed alleles or rare ones ex. HLA:HLA02369 A*03:21N 3095 bp
					$count_hla_allele_loci{$1}++;
				}
			}
			# The major allele locus will be the marker HLA locus
			my $hla_locus = (sort { $count_hla_allele_loci{$b} <=> $count_hla_allele_loci{$a} } keys %count_hla_allele_loci)[0];
			if (defined($hla_locus)){
				printf("\t%s locus assigned to marker '%s'\n", $hla_locus, $marker_name);
				push(@{$hla_loci->{$hla_locus}}, $marker_name);
				if (!defined($samples->{$hla_locus})){ $samples->{$hla_locus} = []; }
				$samples->{$hla_locus} = [ unique(@{$samples->{$hla_locus}},@{$amplisas_results->{$marker_name}{'samples'}}) ];
			} else {
				printf("\tCould not be assigned locus to '%s' marker\n", $marker_name);
				next;
			}
			# Annotates the alleles of the same locus for further proccesing
			foreach my $md5 (keys %{$amplisas_results->{$marker_name}{'seq_data'}}){
				if (!defined($md5_to_name->{$md5})){
					next;
				}
				foreach my $hla_allele (split(' \| ',$md5_to_name->{$md5})) {
					#if ($hla_allele =~ /HLA:HLA\d+ (\w+)\*(\d+)\:(\d+)\:?([\d\:]+[NLSCAQ]?)?/) { # Ex. HLA:HLA00132 B*07:02:01 3323 bp
					if ($hla_allele =~ /(\w+)\*(\d+)\:(\d+)\:?([\d\:]+[NLSCAQ]?)?/) { # Ex. B*07:02:01 3323 bp
						# Discards rare human alelles
						if ($hla_locus eq $1 && in_array($HLA_FREQUENT_ALLELES->{$hla_locus},"$1*$2:$3")) {
							foreach my $sample_name (sort {$a cmp $b} keys %{$amplisas_results->{$marker_name}{'assignments'}{$md5}}){
								push(@{$sample_marker_alleles->{$1}{$sample_name}{$marker_name}{$md5}},"$1*$2:$3");
							}
						}
					}
				}
			}
	#		# Checks all variants assigned to the marker
	# 		foreach my $md5 (keys %{$amplisas_results->{$marker_name}{'seq_data'}}){
	# # 			my @hla_alleles = split(' \| ',$amplisas_results->{$marker_name}{'seq_data'}{$md5}{'name'});
	# 			if (!defined($md5_to_name->{$md5})){
	# 				next;
	# 			}
	# 			my @hla_alleles = split(' \| ',$md5_to_name->{$md5});
	# 			foreach my $hla_allele (@hla_alleles) {
	# 				if ($hla_allele =~ /HLA:HLA\d+ (\w+)\*([\d\:]+[NLSCAQ]?)\s/) { # Ex. HLA:HLA00132 B*07:02:01 3323 bp, exclude non expressed alleles or rare ones ex. HLA:HLA02369 A*03:21N 3095 bp
	# 					$count_hla_allele_loci{$1}++;
	# 					foreach my $sample_name (sort {$a cmp $b} keys %{$amplisas_results->{$marker_name}{'assignments'}{$md5}}){
	# 						push(@{$sample_marker_alleles->{$1}{$sample_name}{$marker_name}{$md5}},$hla_allele);
	# 					}
	# 				}
	# 			}
	# 		}
		}
	}

	# AmpliHLA results will join several markers into a unique HLA locus
	my $amplihla_results;
	foreach my $allele_locus (keys %$hla_loci){
		# Find unique alleles in multiple markers
	# 	foreach my $sample_name (keys %{$sample_marker_alleles->{$allele_locus}}){
		foreach my $sample_name (sort {$a cmp $b} @{$samples->{$allele_locus}}){
	# if ($sample_name eq 'Daudi'){
	# print '';
	# }
			foreach my $marker_name_1 (sort {$a cmp $b} keys %{$sample_marker_alleles->{$allele_locus}{$sample_name}}){
				# Annotates common alleles if a variant is assigned to the same alleles in 2 or more different markers
				# Removes the rest of assignations and annotates it in genotyping results
				foreach my $md5_1 (keys %{$sample_marker_alleles->{$allele_locus}{$sample_name}{$marker_name_1}}){
					my @common_alleles = unique(@{$sample_marker_alleles->{$allele_locus}{$sample_name}{$marker_name_1}{$md5_1}});
					my @common_alleles_freqs = ($amplisas_results->{$marker_name_1}{'assignments'}{$md5_1}{$sample_name} / $amplisas_results->{$marker_name_1}{'sample_data'}{$sample_name}{'depth_amplicon'});
					my %common_alleles_md5s  = ($marker_name_1 => $md5_1);
					my %common_alleles_seqs  = ($marker_name_1 => $amplisas_results->{$marker_name_1}{'seq_data'}{$md5_1}{'sequence'});
					# Removes the variant from further assignations
					delete($sample_marker_alleles->{$allele_locus}{$sample_name}{$marker_name_1}{$md5_1});
					foreach my $marker_name_2 (sort {$a cmp $b} keys %{$sample_marker_alleles->{$allele_locus}{$sample_name}}){
	# 					if ($marker_name_1 eq $marker_name_2) {
	# 						next;
	# 					}
						foreach my $md5_2 (keys %{$sample_marker_alleles->{$allele_locus}{$sample_name}{$marker_name_2}}){
							my @common_alleles_ = intersect(\@common_alleles, $sample_marker_alleles->{$allele_locus}{$sample_name}{$marker_name_2}{$md5_2});
							if (@common_alleles_) {
								@common_alleles = @common_alleles_;
								push(@common_alleles_freqs, $amplisas_results->{$marker_name_2}{'assignments'}{$md5_2}{$sample_name} / $amplisas_results->{$marker_name_2}{'sample_data'}{$sample_name}{'depth_amplicon'});
								$common_alleles_md5s{$marker_name_2} = $md5_2;
								$common_alleles_seqs{$marker_name_2} = $amplisas_results->{$marker_name_2}{'seq_data'}{$md5_2}{'sequence'};
								# Remove the variant from further assignations
								delete($sample_marker_alleles->{$allele_locus}{$sample_name}{$marker_name_2}{$md5_2});
							}
						}
					}
					# Groups alleles with the highest common precision
					my $common_grouped_alleles = group_hla_alleles(@common_alleles);
					# push(@{$amplihla_results->{$allele_locus}{'sample_data'}{$sample_name}{'alleles'}},keys %$common_grouped_alleles);
	# 				foreach my $allele (keys %$common_grouped_alleles){
	# 					my $allele_freq = sprintf('%.2f',sum(@common_alleles_freqs)/(scalar keys %{$sample_marker_alleles->{$allele_locus}{$sample_name}}));
	# 					$amplihla_results->{$allele_locus}{'sample_data'}{$sample_name}{'freq'}{$allele} = $allele_freq;
	# 					$amplihla_results->{$allele_locus}{'sample_data'}{$sample_name}{'alleles'}{$allele} = $common_grouped_alleles->{$allele};
	# 				}
					my $allele = join(' | ', nsort(keys %$common_grouped_alleles));
					my $allele_freq = sprintf('%.2f',sum(@common_alleles_freqs)/(scalar keys %{$sample_marker_alleles->{$allele_locus}{$sample_name}}));
					$amplihla_results->{$allele_locus}{'sample_data'}{$sample_name}{'freq'}{$allele} = $allele_freq;
					$amplihla_results->{$allele_locus}{'allele_data'}{$allele}{'md5'} = \%common_alleles_md5s;
					$amplihla_results->{$allele_locus}{'allele_data'}{$allele}{'seq'} = \%common_alleles_seqs;
					my @common_grouped_alleles_;
					foreach my $allele_ (keys %$common_grouped_alleles){
						push(@common_grouped_alleles_,@{$common_grouped_alleles->{$allele_}});
					}
					$amplihla_results->{$allele_locus}{'sample_data'}{$sample_name}{'alleles'}{$allele} = [ nsort(@common_grouped_alleles_) ];
				}
			}
		}
	}

	# Print genotyping results
	my $excelfile_properties = { 
		'title' => "AmpliHLA genotyping results",
		'author' => "Alvaro Sebastian",
		'comments' => "AmpliHLA genotyping results",
		'company' => "Sixth Researcher - www.sixthresearcher.com",
	};
	my $amplihla_outfile = write_amplihla_file_results("$INP_outpath/results.xlsx",$amplihla_results,$excelfile_properties);

	#print "\nHLA typing results written into '$INP_output'.\n\n";

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


} elsif ($INP_analysis_type eq 'mapping') {

	# Default filename for result files
	my $infile_name;
	if ($INP_reads_files[0] =~ /(.+\/)?(.+?)(_R?1)?\./){
		$infile_name = $2;
	}

	# Reads parameters data from CVS input file if defined
	my %paramsdata;
	if (defined($INP_params_file)){
		my $paramsdata_ = read_amplicon_data($INP_params_file,'params',['only params']);
		# Stores the params in a simple hash
		map $paramsdata{$_} = $paramsdata_->{$_}[0] , keys %$paramsdata_;
	}
	# Command line params have priority over CSV file (if there is no command line neither CSV param, then use default value)
	foreach my $param (keys %{$DEFAULT_PARAMS->{$INP_analysis_type}}) {
		# If specified cancer type is not supported
		if (defined($INP_params{$param})){
			$paramsdata{$param} = lc($INP_params{$param});
		} elsif (!defined($paramsdata{$param})){
			$paramsdata{$param} = $DEFAULT_PARAMS->{$INP_analysis_type}{$param};
		}
	}
	write_to_file("$INP_outpath/amplihla_params.csv", print_amplicon_data(\%paramsdata,'params',[keys %{$DEFAULT_PARAMS->{$INP_analysis_type}}]));

	# Checks if exists reference sequences file
	my $reference_sequences_file = $INP_allele_file;
	if (!defined($reference_sequences_file) || !-e $reference_sequences_file){
		print "\nERROR: '$reference_sequences_file' FASTA file with reference sequences cannot be found.\n\n";
		usage();
		exit;
	} elsif (!is_fasta($reference_sequences_file)) {
		print "\nERROR: '$reference_sequences_file' must be in FASTA format.\n\n";
		exit;
	}

	my ($alleles, $reference_seqs, $align_data, $bamfile);
	if ( $paramsdata{'allele_alignment'} =~ /(bowtie|bowtie2) (.+)/i ){

		my $align_program = lc($1);
		my $align_params = $2;
		# No output messages:
		# $align_params .= " --quiet";
		# Maximum 2 alignments per read:
		#$align_params .= " -k 2";

		if (defined($INP_threads) && $INP_threads>1){
			$align_params .= " --threads $INP_threads";
		}
		if (is_bam($INP_reads_files[0])){
			$bamfile = $INP_reads_files[0];
		} elsif (defined($infile_name)){
			$bamfile = "$INP_outpath/$infile_name.hla_$INP_reference_type.bam";
		} else {
			$bamfile = "$INP_outpath/alleles.hla_$INP_reference_type.bam";
		}
		# Skip mapping if input file is a BAM file or test mode is on and the BAM file already exists
		unless ($bamfile eq $INP_reads_files[0] || (defined($INP_test) && is_bam($bamfile))){
			printf("\nMapping %d sequences to %d allele references with %s.\n",count_seqs_from_fastq($INP_reads_files[0]),count_seqs_from_fasta($reference_sequences_file), ucfirst($align_program));
			my $compression_type = is_compressed($reference_sequences_file);
			if ($compression_type){
				$reference_sequences_file = decompress($reference_sequences_file,$compression_type);
			}
			if ($align_program eq 'bowtie2'){
				execute_bowtie2(\@INP_reads_files,$reference_sequences_file,$align_params,$bamfile);
			} elsif ($align_program eq 'bowtie'){
				execute_bowtie(\@INP_reads_files,$reference_sequences_file,$align_params,$bamfile);
			}
			if ($compression_type){
				`rm $reference_sequences_file`;
			}
			if (!is_bam($bamfile)){
				print "\nERROR: It was an error in the read mapping and BAM file could not be generated.\n\n";
				exit;
			}
			# Indexes BAM file for allele retrieval from references and visualization
			execute_samtools($bamfile, 'index', "$bamfile.bai");
		}

		# Retrieves information about allele mutations from the alignment of reads/alleles to references
		print "\nRetrieving allele information from BAM data.\n";
		$alleles = retrieve_alleles_from_bam_data($bamfile,$reference_sequences_file,\%paramsdata);
		print '';

	} else {
		print "\nERROR: Alignment method not supported: '".$paramsdata{'allele_alignment'}."'.\n\n";
		exit;
	}

	if (!defined($alleles)) {
		print "\nERROR: It was an error retrieving alleles.\n\n";
		exit;
	}

# 	printf("\n%d alleles found.\n", scalar keys %$alleles);

	# Translates HLA sequence IDs (ex. HLA:HLA00005) into HLA allele names (ex. A*02:01:01:01)
	my $hla_alleles;
	# Reads references file
	my $references = read_fasta_file_hash($reference_sequences_file);
	# Parses reference names
	my %ref_name_to_allele;
	foreach my $ref_header (keys %$references){
		# Ex. >HLA:HLA00005 A*02:01:01:01 1098 bp
		my @header_fields = split('\s',$ref_header);
		# Annotates the classical HLA loci (and more)
		# if ($header_fields[1] =~ /(\w+\*\d+\:\d+)
		if ($header_fields[1] =~ /(A|B|C|DQA1|DQB1|DRB1)(\*\d+\:\d+)/){
			$ref_name_to_allele{$header_fields[0]} = $1.$2;
		}
	}
	# Reads allele data and annotate allele info per locus into $hla_alleles
	foreach my $allele_id (sort {$alleles->{$b}{'coverage'}<=>$alleles->{$a}{'coverage'}} keys %$alleles){
	# foreach my $allele_id (sort {$alleles->{$b}{'mapped_area'}<=>$alleles->{$a}{'mapped_area'}} keys %$alleles){
		# Discards non A, B, C, DQA1, DQB1 or DRB1 alleles
		if (!defined($ref_name_to_allele{$allele_id})){
			next;
		}
		my $allele_name = $ref_name_to_allele{$allele_id};
		my $allele_locus;
		if ($allele_name =~ /(A|B|C|DQA1|DQB1|DRB1)(\*\d+\:\d+)/){
			$allele_locus = $1;
		} else {
			next;
		}
		# Discards rare human alleles from A, B, C, DQA1, DQB1 or DRB1 loci
		if (!in_array($HLA_FREQUENT_ALLELES->{$allele_locus},$allele_name)){
			next;
		}
		# Keeps only the maximum coverage allele (as in Boegel et al. 2012)
		if (!defined($hla_alleles->{$allele_locus}{$allele_name}) || $alleles->{$allele_id}{'coverage'} > $hla_alleles->{$allele_locus}{$allele_name}{'coverage'}){
			$hla_alleles->{$allele_locus}{$allele_name} = { %{$alleles->{$allele_id}} };
		}
# 		# NOT RECOMMENDED: VERY SENSITIVE TO THE NUMBER OF ALLELES WITH THE SAME NAME IN THE REFERENCES
# 		# Keeps only the maximum coverage allele, incrementing its coverage with other alleles of the same type
# 		if (!defined($hla_alleles->{$allele_locus}{$allele_name})){
# 			$hla_alleles->{$allele_locus}{$allele_name} = $alleles->{$allele_id};
# 		} else {
# 			$hla_alleles->{$allele_locus}{$allele_name}{'coverage'} += $alleles->{$allele_id}{'coverage'};
# 			if ($hla_alleles->{$allele_locus}{$allele_name}{'mapped_length'} < $alleles->{$allele_id}{'mapped_length'}) {
# 				$hla_alleles->{$allele_locus}{$allele_name}{'mapped_length'} = $alleles->{$allele_id}{'mapped_length'};
# 			}
# 			$hla_alleles->{$allele_locus}{$allele_name}{'mapped_seqs'} = [ @{$hla_alleles->{$allele_locus}{$allele_name}{'mapped_seqs'}}, @{$alleles->{$allele_id}{'mapped_seqs'}} ];
# 		}
	}
	# Releases memory
	undef($alleles);
	
	# Annotates the major alleles per locus
	my ($grouped_hla_alleles, $sorted_alleles, $sorted_coverages);
	foreach my $allele_locus (@HLA_LOCI){
		# Counts the number of mapped reads for each allele
		# map $hla_alleles->{$allele_locus}{$_}{'count_mapped_seqs'} = scalar @{$hla_alleles->{$allele_locus}{$_}{'mapped_seqs'}}, keys %{$hla_alleles->{$allele_locus}};
		my (%corrected_allele_coverage, %all_previous_mapped_seqs);
		my ($previous_allele, $previous_allele_coverage, %previous_allele_mapped_seqs);
		# Re-calculates allele coverage, discarding from coverage reads previously mapped by other alleles (as in Boegel et al. 2012)
		foreach my $allele_name (sort {$hla_alleles->{$allele_locus}{$b}{'coverage'}<=>$hla_alleles->{$allele_locus}{$a}{'coverage'}} keys %{$hla_alleles->{$allele_locus}}){
		# foreach my $allele_name (sort {$hla_alleles->{$allele_locus}{$b}{'mapped_area'}<=>$hla_alleles->{$allele_locus}{$a}{'mapped_area'}} keys %{$hla_alleles->{$allele_locus}}){
			my @mapped_seqs = unique(@{$hla_alleles->{$allele_locus}{$allele_name}{'mapped_seqs'}});
# 			my $coverage = ${$hla_alleles->{$allele_locus}{$allele_name}}{'coverage'};
			my $coverage = scalar @mapped_seqs;
			# Checks if the allele has lower coverage than former one, or lower mapped length, or more mapping errors
			if (!defined($previous_allele)
			|| $coverage < $previous_allele_coverage
			|| $hla_alleles->{$allele_locus}{$allele_name}{'mapped_length'} < $hla_alleles->{$allele_locus}{$previous_allele}{'mapped_length'}
			|| $hla_alleles->{$allele_locus}{$allele_name}{'errors'} > $hla_alleles->{$allele_locus}{$previous_allele}{'mapped_length'}
			){
				# Incorporates reads mapped by previous allele for checking
				%all_previous_mapped_seqs = (%all_previous_mapped_seqs, %previous_allele_mapped_seqs);
				undef(%previous_allele_mapped_seqs);
				foreach my $mapped_seq (@mapped_seqs){
					if (defined($all_previous_mapped_seqs{$mapped_seq})){
						$coverage--;
					}
				}
				$corrected_allele_coverage{$allele_name} = $coverage;
			# If the allele has the same coverage than previous one
			} else {
				# Checks reads mapped by previous alleles with lower coverage (but not the mapped by the previous one)
				foreach my $mapped_seq (@mapped_seqs){
					if (defined($all_previous_mapped_seqs{$mapped_seq})){
						$coverage--;
					}
				}
				# If still coverage is still the same, groups the allele names
				# but does not annotate the coverage to skip the allele from DOC calculation
				if ($coverage == $previous_allele_coverage){
					$grouped_hla_alleles->{$allele_locus}{$previous_allele}{$allele_name} = 1;
					$grouped_hla_alleles->{$allele_locus}{$previous_allele}{$previous_allele} = 1;
				} else {
					$corrected_allele_coverage{$allele_name} = $coverage;
				}
			}
			$previous_allele_coverage = $coverage;
			$previous_allele = $allele_name;
			map $previous_allele_mapped_seqs{$_} = 1 , @mapped_seqs;
		}
		# Skip locus without mapped reads to alleles
		if (!%corrected_allele_coverage){
			next;
		}
		# Sorts alleles according to re-calculated coverages
		$sorted_alleles->{$allele_locus} = [ sort {$corrected_allele_coverage{$b}<=>$corrected_allele_coverage{$a}} keys %corrected_allele_coverage ];
		$sorted_coverages->{$allele_locus} = [ map $corrected_allele_coverage{$_} , @{$sorted_alleles->{$allele_locus}} ];
# 		# DEBUGGING INFO:
# 		print "\n$allele_locus\n";
# 		printf("%12s\t%s\n",'ALLELE',join("\t",('LEN','COV','COV2','LENCOV','LENCOV2','SCORE','SCORE2')));
# 		my ($total_coverage, $total_coverage2);
# 		map $total_coverage += $_ , @{$sorted_coverages->{$allele_locus}};
# 		map $total_coverage2 += $_ , values %corrected_allele_coverage;
# 		# foreach my $allele_name (sort {$hla_alleles->{$allele_locus}{$b}{'mapped_area'}<=>$hla_alleles->{$allele_locus}{$a}{'mapped_area'}} keys %{$hla_alleles->{$allele_locus}}){
# 		foreach my $allele_name (@{$sorted_alleles->{$allele_locus}}){
# 			my $mapped_len = $hla_alleles->{$allele_locus}{$allele_name}{'mapped_length'};
# 			my $coverage = $hla_alleles->{$allele_locus}{$allele_name}{'coverage'};
# 			my $mapped_area = $mapped_len*$coverage;
# 			my $score= 100*$coverage/$total_coverage;
# 			my ($coverage2, $mapped_area2, $score2) = (0,0,0);
# 			if ($corrected_allele_coverage{$allele_name}) {
# 				$coverage2 = $corrected_allele_coverage{$allele_name};
# 				$mapped_area2 = $mapped_len*$coverage2;
# 				$score2 = 100*$coverage2/$total_coverage2;
# 			}
# 			if (!defined($grouped_hla_alleles->{$allele_locus}{$allele_name})){
# 				printf("%12s\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f\n", $allele_name, $mapped_len, $coverage, $coverage2, $mapped_area, $mapped_area2, $score, $score2);
# 			} else {
# 				foreach my $allele_name_ (nsort(keys %{$grouped_hla_alleles->{$allele_locus}{$allele_name}})){
# 					my $mapped_len_ = $hla_alleles->{$allele_locus}{$allele_name_}{'mapped_length'};
# 					my $coverage_ = $hla_alleles->{$allele_locus}{$allele_name_}{'coverage'};
# 					my $mapped_area_ = $mapped_len*$coverage;
# 					my $score_ = 100*$coverage_/$total_coverage;
# 					my ($coverage2_, $mapped_area2_, $score2_) = (0,0,0);
# 					if ($corrected_allele_coverage{$allele_name}) {
# 						$coverage2_ = $corrected_allele_coverage{$allele_name};
# 						$mapped_area2_ = $mapped_len*$coverage2;
# 						$score2_ = 100*$coverage2_/$total_coverage2;
# 					}
# 					printf("%12s\t%d\t%d\t%d\t%d\t%d\t%.1f\t%.1f\n", $allele_name_, $mapped_len_, $coverage_, $coverage2_, $mapped_area_, $mapped_area2_, $score_, $score2_);
# 				}
# 			}
# 		}
# 		print "\n";
	}
	# Releases memory
	undef($hla_alleles);
	
	# Prints the results
	print "\nRESULTS:\n";
	printf("\n%s\t%s\t\t%s\t%s\t%s\n", 'LOCI', 'ALLELE', 'SCORE', 'READS', 'AMBIGUITIES');
	foreach my $allele_locus (@HLA_LOCI){
		if (!defined($sorted_coverages->{$allele_locus}) || !@{$sorted_coverages->{$allele_locus}}){
			printf("\n%s\tNO RESULTS\n",$allele_locus);
			next;
		}
		# If coverage is zero in 2nd or 3rd position, replaces 0 by 1
		if (defined($sorted_coverages->{$allele_locus}[1]) && $sorted_coverages->{$allele_locus}[1] == 0){
			$sorted_coverages->{$allele_locus}[1] = 1;
		} elsif (defined($sorted_coverages->{$allele_locus}[2]) && $sorted_coverages->{$allele_locus}[2] == 0){
			$sorted_coverages->{$allele_locus}[2] = 1;
		}
		my $total_coverage;
		map $total_coverage += $_ , @{$sorted_coverages->{$allele_locus}};
		# Assigns 1 or 2 alleles according to re-calculated coverages
		my $allele_number;
		if ($sorted_coverages->{$allele_locus}[0] > 5*$sorted_coverages->{$allele_locus}[1]){
			$allele_number = 1;
		} else {
			my $DOCns;
			($allele_number,$DOCns) = degree_of_change($sorted_coverages->{$allele_locus},2);
		}
		for (my $i=0; $i<$allele_number; $i++){
			if ($i == 0){
				printf("\n%s\t",$allele_locus);
			} else {
				printf("%s\t",'');
			}
			my $score = 100*$sorted_coverages->{$allele_locus}[$i]/$total_coverage;
			my $coverage = $sorted_coverages->{$allele_locus}[$i];
			my $allele = $sorted_alleles->{$allele_locus}[$i];
			my $ambiguities = '';
			if (defined($grouped_hla_alleles->{$allele_locus}{$sorted_alleles->{$allele_locus}[$i]})){
				# $allele = join('|',nsort(keys %{$grouped_hla_alleles->{$allele_locus}{$sorted_alleles->{$allele_locus}[$i]}}));
				my $common_grouped_alleles = group_hla_alleles(keys %{$grouped_hla_alleles->{$allele_locus}{$sorted_alleles->{$allele_locus}[$i]}});
				$allele = join(' | ', nsort(keys %$common_grouped_alleles));
				$ambiguities = join(',',nsort(keys %{$grouped_hla_alleles->{$allele_locus}{$sorted_alleles->{$allele_locus}[$i]}}));
			}
			if (length($allele)<10){
				printf("%s\t\t%.1f\t%d\t%s\n", $allele, $score, $coverage, $ambiguities);
			} else {
				printf("%s\t%.1f\t%d\t%s\n", $allele, $score, $coverage, $ambiguities);
			}
		}
	}
	print "\n";
}

# exit;


#########################################################################################

# Retrieves a common HLA genotype from a list of different (or similar) alleles
sub group_hla_alleles {

	my @input_alleles = @_;

	my $grouped_alleles;

# 	for (my $i=0; $i<=$#input_alleles; $i++) {
# 		my $allele1 = $input_alleles[$i];
	while (@input_alleles) {
		my $allele1 = shift @input_alleles;
		my ($type1, @fields1);
		if ($allele1 =~ /(\w+)\*([\d\:]+[NLSCAQ]?)/) { # Ex. HLA:HLA00132 B*07:02:01N 3323 bp
			$type1 = $1;
			@fields1 = split(':',$2);
			$allele1 = "$1*$2";
		}
		my $common_fields = 0;
		my @common_alleles;
		for (my $j=0; $j<=$#input_alleles; $j++) {
			my $allele2 = $input_alleles[$j];
			my ($type2, @fields2);
			if ($allele2 =~ /(\w+)\*([\d\:]+[NLSCAQ]?)/) {
				$type2 = $1;
				@fields2 = split(':',$2);
				$allele2 = "$1*$2";
			}
			if ($type1 ne $type2 || $allele1 eq $allele2) {
				next;
			}
			my $field_matchs = 0;
			for (my $n=0; $n<=$#fields1; $n++) {
				if ($fields1[$n] eq $fields2[$n]){
					$field_matchs++;
				} else {
					last;
				}
			}
			if ($field_matchs){ # Removes matched and annotates the lowest digit matched
				push(@common_alleles, $allele2);
				splice(@input_alleles,$j,1); $j--;
				if (!$common_fields || $field_matchs<$common_fields){
					$common_fields = $field_matchs;
				}
			}
		}
		if ($common_fields){ # If they have digits in common then simplify
			unshift(@common_alleles, $allele1);
			my $allele_consensus = $type1.'*'.join(':',splice(@fields1,0,$common_fields));
			push(@{$grouped_alleles->{$allele_consensus}}, @common_alleles);
		} else {
			my $allele_consensus = $type1.'*'.join(':',@fields1);
			push(@{$grouped_alleles->{$allele_consensus}}, $allele1);
		}
	}
	
	return ($grouped_alleles);

}

#########################################################################################

# Retrieves a common HLA genotype comparing a list of different (or similar) alleles to another list
sub group_hla_alleles_between {

	my ($query_alleles_,$sbjct_alleles_)= @_;

	my $query_alleles = [ @$query_alleles_ ];
	my $sbjct_alleles = [ @$sbjct_alleles_ ];

	my $grouped_alleles;

	while (@{$query_alleles}) {
		my $allele1 = shift @{$query_alleles};
		my ($type1, @fields1);
		if ($allele1 =~ /(\w+)\*([\d\:]+[NLSCAQ]?)/) { # Ex. HLA:HLA00132 B*07:02:01 3323 bp
			$type1 = $1;
			@fields1 = split(':',$2);
			$allele1 = "$1*$2";
		}
		my $common_fields = 0;
		my @common_alleles;
		for (my $j=0; $j<=$#{$sbjct_alleles}; $j++) {
			my $allele2 = $sbjct_alleles->[$j];
			my ($type2, @fields2);
			if ($allele2 =~ /(\w+)\*([\d\:]+[NLSCAQ]?)/) {
				$type2 = $1;
				@fields2 = split(':',$2);
				$allele2 = "$1*$2";
			}
			if ($type1 ne $type2 || $allele1 eq $allele2) {
				next;
			}
			my $field_matchs = 0;
			for (my $n=0; $n<=$#fields1; $n++) {
				if ($fields1[$n] eq $fields2[$n]){
					$field_matchs++;
				} else {
					last;
				}
			}
			if ($field_matchs){ # Removes matched and annotates the lowest digit matched
				push(@common_alleles, $allele2);
				splice(@{$sbjct_alleles},$j,1); $j--;
				if (!$common_fields || $field_matchs<$common_fields){
					$common_fields = $field_matchs;
				}
			}
		}
		if ($common_fields){ # If they have digits in common then simplify
			unshift(@common_alleles, $allele1);
			my $allele_consensus = $type1.'*'.join(':',splice(@fields1,0,$common_fields));
			push(@{$grouped_alleles->{$allele_consensus}}, @common_alleles);
		} else {
			my $allele_consensus = $type1.'*'.join(':',@fields1);
			push(@{$grouped_alleles->{$allele_consensus}}, $allele1);
		}
	}
	
	return ($grouped_alleles);

}

#########################################################################################



