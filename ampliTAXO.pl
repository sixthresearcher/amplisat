#!/usr/bin/perl -w
#
# Name: ampliTAXO.pl
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
#   Performs taxonomic classification of samples based in amplicon sequencing data.
#
# Requires as input a FASTA or FASTQ file with sequences/reads and a CSV format file with primer/tag data information
# Example: perl ampliTAXO.pl -d amplicon_data.csv -i reads.fq.gz -o results
#
# If the reads have been already demultiplexed into separate files (one file per sample), they can be packed into a single .zip or .tar.gz file and use it as input
# Example: perl ampliTAXO.pl -i reads.tar.gz -o results
#
# Alternatively, variants already extracted can be given in an Excel format file obtained with AmpliSAS/AmpliCHECK.
# Example: perl ampliTAXO.pl -i amplisas_results.xlsx -o results
#


my $VERSION = "1.0";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Performs taxonomic classification based in amplicon sequencing data.";


# Modules are in folder 'lib' in the path of the script
use lib "lib";
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;
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
my %DEFAULT_PARAMS = (
	# Allele matching parameters
	'allele_alignment' => 'dna blastn -evalue 1E-5 -ungapped',
	# Minimum % of sequence aligned to a reference
	'min_allele_align' => 70,
	# Minimum % of identity in the portion aligned to a reference
	'min_allele_ident' => 90,
	# Amplicons with lower total coverage will be discarded
	'min_amplicon_depth' => 100,
	# Minimum frequency of variants after clustering
	'min_amplicon_seq_frequency' => 0.002, # 0.005%
	# OTUs with equal or higher identity will be clustered together
	'identity_threshold' => undef,
	# Keep singletons after clustering
	'keep_singletons' => 'no',
	# Default reference database
	'taxo_database' => 'greengenes',
	# Lowest taxonomic level for grouping OTUs together
	'taxo_level' => 'genus',
	# Alpha diversity metric to perform rarefaction analysis
	'alpha_diversity' => 'otus',
	# Step size sample depth to perform rarefaction of reads
	'rarefraction_depth' => 'auto',
	# Number of rarefaction iterations/replicates for each sample size
	'rarefraction_replicates' => 10,
);

# Alpha diversity measures available
my %ALPHA_DIVERSITY = ( 'otus' => 'OTU number', 'chao1' => 'Chao1' );

# rRNA and ITSs default databases
# Created with 'reformat_amplitaxol.pl' script
my $TAXO_DATABASES = {
	'silva-ssu' => { 'name' => 'SILVA-SSU', 'description' => 'A comprehensive online resource for quality checked and aligned ribosomal RNA sequence data',
		'size' => 597607, 'file' => $SCRIPT_DIR.'/taxo/SILVA_123_SSURef_Nr99_tax_silva.utax.fa.gz', 'version' => 'SSU Ref NR99', 'url' => 'http://www.arb-silva.de', 'pubmed'=> '23193283' },
	'silva-lsu' => { 'name' => 'SILVA-LSU', 'description' => 'A comprehensive online resource for quality checked and aligned ribosomal RNA sequence data',
		'size' => 96642, 'file' => $SCRIPT_DIR.'/taxo/SILVA_123_LSURef_tax_silva.utax.fa.gz', 'version' => 'LSU Ref 123', 'url' => 'http://www.arb-silva.de', 'pubmed'=> '23193283' },
	'unite' => { 'name' => 'UNITE', 'description' => 'Unified system for the DNA based fungal species linked to the classification',
		'size' => 22770, 'file' => $SCRIPT_DIR.'/taxo/sh_general_release_dynamic_01.08.2015.utax.fa.gz', 'version' => 'v7', 'url' => 'https://unite.ut.ee', 'pubmed'=> '24112409' },
	'greengenes' => { 'name' => 'Greengenes', 'description' => 'Greengenes, a chimera-checked 16S rRNA gene database (only 90% identity core set)', # Greengenes clustered with cd-hit-est to 90% identity
		'size' => 4799, 'file' => $SCRIPT_DIR.'/taxo/greengenes_gg16S_90perc.utax.fa.gz', 'version' => '20110509', 'url' => 'http://greengenes.lbl.gov', 'pubmed'=> '16820507' },
	'core' => { 'name' => 'CORE',  'description' => 'A phylogenetically-curated 16S rDNA database of the core oral microbiome',
		'size' => 1262, 'file' => $SCRIPT_DIR.'/taxo/oralDB_012814.utax.fa.gz', 'version' => '20140128', 'url' => 'http://microbiome.osu.edu', 'pubmed' => '21544197' },
	'homd' => { 'name' => 'HOMD', 'description' => 'The Human Oral Microbiome Database',
		'size' => 831, 'file' => $SCRIPT_DIR.'/taxo/HOMD_16S_rRNA_RefSeq_V13.2.utax.fa.gz', 'version' => 'v13.2', 'url' => 'http://www.homd.org', 'pubmed' => '20624719' },
	'test' => { 'name' => 'TEST', 'description' => '',
		'size' => 1, 'file' => $SCRIPT_DIR.'/kk.utax.fa', 'version' => '', 'url' => '', 'pubmed'=> '' },
};
# Array with all available databases in lowercase
my @TAXO_DATABASES = map $TAXO_DATABASES->{$_}{'name'}, sort {$TAXO_DATABASES->{$b}{'size'} <=> $TAXO_DATABASES->{$a}{'size'}} keys %$TAXO_DATABASES;

# Input parameters variables (will be filled with command line options or default values)
my %INP_params;

my ($INP_reads_file,$INP_amplicons_file,$INP_outpath,$INP_allele_file,@INP_taxo_databases,$INP_shuffle,$INP_tech,$INP_nreads,$INP_nreads_amplicon,$INP_threads,$INP_zip,$INP_verbose);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s' => \$INP_reads_file,
	'd|data=s' => \$INP_amplicons_file,
	'o|output=s' => \$INP_outpath,
	'a|alleles=s' => \$INP_allele_file,
	'db|database=s{,}' => \@INP_taxo_databases,
	't|tech=s' => \$INP_tech,
	's|shuffle' => \$INP_shuffle,
	'n|number=i' => \$INP_nreads,
	'na|max=i' => \$INP_nreads_amplicon,
	'l|level=s' => \$INP_params{'taxo_level'},
	'al|aligned=f' => \$INP_params{'min_allele_align'},
	'id|ident=f' => \$INP_params{'min_allele_ident'},
	'fr|freq=f' => \$INP_params{'min_amplicon_seq_frequency'},
	'ci|clid=f' => \$INP_params{'identity_threshold'},
	'min=i' => \$INP_params{'min_amplicon_depth'},
	'ks|keepsingle' => \$INP_params{'keep_singletons'},
	'alpha:s' => \$INP_params{'alpha_diversity'},
	'rare:i' => \$INP_params{'rarefraction_depth'},
	'thr|threads=i' => \$INP_threads,
	'z|zip' => \$INP_zip,
	'v|verbose' => \$INP_verbose,
);

# Usage help
sub usage
{
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -i <file>\tInput FASTQ or FASTA file (compressed or uncompressed).\n";
	print "  -d <file>\tCSV file with primer/amplicon data.\n";
	print "  -o <path>\tOutput folder name.\n";
	print "  -db <db1> [ <db2> ... <dbN> ]\n\t\tDatabases to use in the analysis ('".join("','",@TAXO_DATABASES)."').\n\t\tDefault: ".$DEFAULT_PARAMS{'taxo_database'}."\n";
	print "  -a <file>\tFASTA file with sequences and taxonomic information in UTAX format.\n"; # (default=SILVA+UNITE).\n";
	print "  -t <tech>\tUse recommended technology parameters ('Illumina', 'IonTorrent', '454', 'Unknown').\n";
	print "  -s\t\tShuffle/randomize reads/sequences to analyze.\n";
	print "  -n <number>\tNumber of reads/sequences to analyze.\n";
	print "  -na <number>\tNumber of reads/sequences per amplicon to analyze.\n";
	print "  -l <level>\tLowest taxonomic level for grouping OTUs together (default='".$DEFAULT_PARAMS{'taxo_level'}."').\n";
	print "  -al <perc>\tMinimum sequence percentage required to align to a reference to be annotated (default=".$DEFAULT_PARAMS{'min_allele_align'}."%).\n";
	print "  -id <perc>\tMinimum percentage of the aligned fragment required to be identical to a reference to be annotated (default=".$DEFAULT_PARAMS{'min_allele_ident'}."%).\n";
	print "  -fr <freq>\tFilter sequences with lower frequency after clustering (default=".$DEFAULT_PARAMS{'min_amplicon_seq_frequency'}."%).\n";
	print "  -ci <number>\tCluster together sequences with higher or equal identity.\n";
	print "  -min <number>\tAmplicons with lower total depth/coverage will be discarded (default=".$DEFAULT_PARAMS{'min_amplicon_depth'}.").\n";
	print "  -ks\t\tKeep singletons after clustering (default=".$DEFAULT_PARAMS{'keep_singletons'}.").\n";
	print "  -alpha <type>\tAlpha diversity metric to perform rarefaction analysis (default=".$DEFAULT_PARAMS{'alpha_diversity'}.").\n";
	print "  -rare <number>Number of sample sequences to increment in each step of rarefaction analysis (default=".$DEFAULT_PARAMS{'rarefraction_depth'}.").\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -v\t\tPrint AmpliSAS output and verbose FASTA files.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
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
} elsif (!defined($INP_reads_file)){
	print "\nERROR: You must specify sequence input file.\n\n";
	usage();
	exit;
# Last version accepts single amplicon files:
} elsif (!defined($INP_amplicons_file)){
	$INP_multifile = $INP_reads_file;
}
# } elsif (!defined($INP_amplicons_file)){
# 	print "\nERROR: You must specify amplicon data input file.\n\n";
# 	usage();
# 	exit;
# }

# Checks if there is a references file
if (defined($INP_allele_file) && !-e $INP_allele_file){
	print "\nERROR: FASTA file with sequences and taxonomic information cannot be opened.\n\n";
	usage();
	exit;
}

# AMPLITAXO ONLY ALLOWS GLOBAL PARAMETERS FOR ALL THE MARKERS ('all')
# Reads parameters data from CVS input file
my %paramsdata;
if (defined($INP_amplicons_file)){
	my $paramsdata_ = read_amplicon_data($INP_amplicons_file, 'params');
	# Stores the params in a simple hash
	map $paramsdata{$_} = $paramsdata_->{$_}{'all'}[0] , keys %$paramsdata_;
}
# Command line params have priority over CSV file and auto mode
foreach my $param (keys %DEFAULT_PARAMS) {
	if (defined($INP_params{$param})){
		$paramsdata{$param} = lc($INP_params{$param});
	} elsif (!defined($paramsdata{$param}) && defined($DEFAULT_PARAMS{$param})){
		$paramsdata{$param} = lc($DEFAULT_PARAMS{$param});
	}
}
# Checks percentage threshold parameters, if defined
foreach my $param (('min_allele_align','min_allele_ident','min_amplicon_seq_frequency','identity_threshold')){
	if (defined($paramsdata{$param})){
		my $value = $paramsdata{$param};
		if (defined($value)){
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
# Checks if parameters have valid values
if (!in_array(['kingdom','domain','phylum','class','order','family','genus','species'],$paramsdata{'taxo_level'})){
	print "\nERROR: You must specify a valid taxonomic level (kingdom, domain, phylum, class, order, family, genus or species).\n\n";
	usage();
	exit;
}
if (!defined($ALPHA_DIVERSITY{$paramsdata{'alpha_diversity'}})){
	print "\nERROR: You must specify a valid alpha diversity metric : '".join("', '",keys %ALPHA_DIVERSITY)."'.\n\n";
	usage();
	exit;
}
if ($paramsdata{'rarefraction_depth'} ne 'auto' && (!is_numeric($paramsdata{'rarefraction_depth'}) || $paramsdata{'rarefraction_depth'}<10) ){
	print "\nERROR: You must specify a valid number of sample sequences higher than 10.\n\n";
	usage();
	exit;
}

# Default output path
if (!defined($INP_outpath)){
	$INP_outpath  = lc((split('\.',$SCRIPT_NAME))[0]);
}

print "\nRunning '$COMMAND_LINE'\n";

# # Checks and reads alleles file if defined
# # If any allele name or sequence is duplicated a message will be shown from 'read_alleles' function
# my $alleledata = {};
# my @allelesources;
# if (defined($INP_allele_file) && -e $INP_allele_file){
# 	push(@allelesources,$INP_allele_file);
# 	print "\nReading sequences and taxonomic information from '$INP_allele_file'.\n";
# 	$alleledata = read_fasta_file_hash_($INP_allele_file);
# 	if (!defined($alleledata) || !%$alleledata) {
# 		print "\nERROR: FASTA file with sequences and taxonomic information cannot be opened.\n\n";
# 		usage();
# 		exit;
# 	} else {
# 		printf("\t%d sequences read\n", scalar keys %$alleledata);
# 	}
# }
# Checks and reads taxonomic databases to be used in the analysis
if (@INP_taxo_databases) {
	@INP_taxo_databases = map lc($_) , @INP_taxo_databases;
	my @db_names;
	foreach my $db (@INP_taxo_databases) {
		if (!defined($TAXO_DATABASES->{$db})){
			print "\nERROR: '$db' is not a valid database.\n\n";
			usage();
			exit;
		} else {
			push(@db_names,$TAXO_DATABASES->{$db}{'name'});
		}
	}
	if ($#INP_taxo_databases>0) {
		print "\nDatabases '".join("', ",@db_names)."' will be used in the analysis.\n";
	} else {
		print "\nDatabase '".join("', ",@db_names)."' will be used in the analysis.\n";
	}
# If not defined allele file neither databases, uses default ones
} elsif (!defined($INP_allele_file)) {
	@INP_taxo_databases = ($paramsdata{'taxo_database'});
	if ($#INP_taxo_databases>0) {
		print "\nDefault databases '".join("', ",@INP_taxo_databases)."' will be used in the analysis.\n";
	} else {
		print "\nDefault database '".join("', ",@INP_taxo_databases)."' will be used in the analysis.\n";
	}
}
# # Reads sequences and taxonomic information from the databases
# if (@INP_taxo_databases) {
# 	foreach my $db (@INP_taxo_databases) {
# 		my $db_info = $TAXO_DATABASES->{lc($db)};
# 		push(@allelesources,$db_info->{'name'});
# 		printf("\nReading sequences and taxonomic information from '%s %s' database.\n", $db_info->{'name'}, $db_info->{'version'});
# 		my $alleledata_ = read_fasta_file_hash_($SCRIPT_DIR."/".$db_info->{'file'});
# 		if (!defined($alleledata_) || !%$alleledata_) {
# 			printf("\nERROR: '%s %s' database cannot be opened.\n\n", $db_info->{'name'}, $db_info->{'version'});
# 			usage();
# 			exit;
# 		} else {
# 			printf("\t%d sequences read\n", scalar keys %$alleledata_);
# 		}
# 		foreach my $allele_name (keys %$alleledata_) {
# 			if (defined($alleledata->{$allele_name}) ){
# 				print "\tERROR: Allele name '".$allele_name."' is duplicated.\n";
# 	# 			exit;
# 			} else {
# 				# Removes gaps to avoid errors in further alignments
# 				$alleledata_->{$allele_name} =~ s/-//g;
# 				$alleledata->{$allele_name} = $alleledata_->{$allele_name};
# 			}
# 		}
# 	}
# }
# if (!defined($alleledata) || !%$alleledata) {
# 	print "\nERROR: sequences and taxonomic information cannot be read.\n\n";
# 	usage();
# 	exit;
# }

# Runs AmpliSAS in auto mode before HLA typing if reads are given as input file
my $amplisas_results;
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
	if (defined($paramsdata{'keep_singletons'}) && $paramsdata{'keep_singletons'} !~ /no/i){
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
		$amplisas_command = "$SCRIPT_DIR/ampliSAS.pl -denovo -i $INP_reads_file -d $INP_amplicons_file -o $INP_outpath $amplisas_options";
	} else {
		$amplisas_command = "$SCRIPT_DIR/ampliSAS.pl -denovo -i $INP_reads_file -o $INP_outpath $amplisas_options";
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
		printf("Reading AmpliSAS results.\n");
		$amplisas_results = read_amplisas_file_results("$INP_outpath/results.xlsx");
	} else {
		`rm -rf $INP_outpath`;
		if ($amplisas_output =~ /(ERROR:?\s*.+)/ ) {
			print "\n$1\n\n";
		} else {
# 			if ($amplisas_output !~ /\nUsage:/){
# 				printf("\nAmpliSAS output:\n\t%s\n",join("\n\t",split("\n",$amplisas_output)));
# 			}
			print "\nERROR: AmpliSAS failed analyzing the data, try to run AmpliSAS independently and use the Excel results file as AmpliTAXO input.\n\n";
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


# Stores sequences into a hash associated to their MD5 signatures
# And creates a sorted array with unique samples from all the markers
# Because all the markers will be used to assign unique taxonomies per sample
my (%md5_to_sequence,@unique_samples);
foreach my $marker_name (keys %{$amplisas_results}){
	if (defined($amplisas_results->{$marker_name})){
		foreach my $md5 (keys %{$amplisas_results->{$marker_name}{'seq_data'}}){
			$md5_to_sequence{$md5} = $amplisas_results->{$marker_name}{'seq_data'}{$md5}{'sequence'};
		}
		@unique_samples = unique(@unique_samples, @{$amplisas_results->{$marker_name}{'samples'}});
	}
}
@unique_samples = nsort(@unique_samples);

# Assigns Taxonomy (names containing taxonomic info) to sequences based in the most similar sequences in the FASTA file with taxonomic data
# $INP_allele_align_params defines the similarity threshold to assign names
# print "\nMatching taxonomic sequences.\n";
my $allele_align_params = {
	'amplitaxo' => 1,
	'alignment' => $paramsdata{'allele_alignment'},
	'aligned' => sprintf("%.2f",$paramsdata{'min_allele_align'}/100),
	'ident' => sprintf("%.2f",$paramsdata{'min_allele_ident'}/100),
};
# my $md5_to_name = match_alleles($alleledata,\%md5_to_sequence,undef,$allele_align_params,$INP_threads);
# printf ("\nMatched %d alleles in %d unique sequences\n\n", scalar keys %$md5_to_name, scalar keys %md5_to_sequence);
# exit;

my $md5_to_name = {};
print "\nAssigning taxonomies by matching reference sequences.\n";
foreach my $db (@INP_taxo_databases) {
	my $db_info = $TAXO_DATABASES->{lc($db)};
	my $allele_count = count_seqs_from_file($db_info->{'file'});
	printf("\tReading %d sequences and taxonomic information from '%s %s' database.\n", $allele_count, $db_info->{'name'}, $db_info->{'version'});
# 	if ($allele_count < 100000) {
	if (!defined($INP_threads)){
		$allele_align_params->{ 'alignment'} = 'dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 90';
	} else {
		$allele_align_params->{ 'alignment'} = "dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 90 -num_threads $INP_threads";
	}
# 	} else {
# 		if (!defined($INP_threads)){
# 			$allele_align_params->{'alignment'} = 'dna bowtie2 --sensitive-local -k 2';
# 		} else {
# 			$allele_align_params->{'alignment'} = "dna bowtie2 --sensitive-local -k 2 --threads $INP_threads";
# 		}
# 	}
	my $matched_alleles = match_alleles($db_info->{'file'},\%md5_to_sequence,undef,$allele_align_params); # $INP_threads
	if (defined($matched_alleles)) {
		$md5_to_name = { %$md5_to_name, %{$matched_alleles} };
	}
	printf ("\tMatched %d alleles in %d unique sequences\n\n", scalar keys %$md5_to_name, scalar keys %md5_to_sequence);
# my $md5_to_name = match_alleles($alleledata,\%md5_to_sequence,undef,$allele_align_params,$INP_threads);
}
if (defined($INP_allele_file)) {
	my $allele_count = count_seqs_from_file($INP_allele_file);
	printf("\tReading %d sequences and taxonomic information from '%s'.\n", $allele_count, $INP_allele_file);
# 	if ($allele_count < 100000) {
	if (!defined($INP_threads)){
		$allele_align_params->{ 'alignment'} = 'dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 90';
	} else {
		$allele_align_params->{ 'alignment'} = "dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 90 -num_threads $INP_threads";
	}
# 	} else {
# 		if (!defined($INP_threads)){
# 			$allele_align_params->{'alignment'} = 'dna bowtie2 --sensitive-local -k 2';
# 		} else {
# 			$allele_align_params->{'alignment'} = "dna bowtie2 --sensitive-local -k 2 --threads $INP_threads";
# 		}
# 	}
	my $matched_alleles = match_alleles($INP_allele_file,\%md5_to_sequence,undef,$allele_align_params); # $INP_threads
	if (defined($matched_alleles)) {
		$md5_to_name = { %$md5_to_name, %$matched_alleles };
	}
	print '';
	printf ("\tMatched %d alleles in %d unique sequences\n\n", scalar keys %$md5_to_name, scalar keys %md5_to_sequence);
}

# If there are not taxonomy matches
if (!defined($md5_to_name)){
	print "\nERROR: No similar sequences were found in the reference database, please check the correctness of the analyzed sequences and the database ones.\n\n";
	exit;
}

# Annotates all the Taxonomic assignments and reference alignment details associated to the amplicon variants (also called Operational Taxonomic Units or OTUs)
# Also will unify the taxonomic assignments of several OTUs (Operational Taxonomic Units) and several markers
my @taxo_levels = ('class',$paramsdata{'taxo_level'});
my ($amplitaxo_results,$md5_to_taxo,%otu_depths,$otu_ambiguities,@otu_seqs,@otu_headers,@otu_names);
foreach my $marker_name (keys %{$amplisas_results}){
	foreach my $md5 (keys %{$amplisas_results->{$marker_name}{'seq_data'}}){
		$otu_depths{$md5} += $amplisas_results->{$marker_name}{'seq_data'}{$md5}{'depth'};
	}
}
my @otu_md5s = sort { $otu_depths{$b} <=> $otu_depths{$a} } keys %otu_depths;
foreach my $md5 (@otu_md5s) {
	foreach my $marker_name (keys %{$amplisas_results}){
		if (!defined($amplisas_results->{$marker_name}{'seq_data'}{$md5})){
			next;
		}

		my $otu_name = $amplisas_results->{$marker_name}{'seq_data'}{$md5}{'name'};
		my $otu_seq = $amplisas_results->{$marker_name}{'seq_data'}{$md5}{'sequence'};
		my ($otu_taxo,$otu_header);

		if (defined($md5_to_name->{$md5})) {

			my (@otu_taxo_assigns,@otu_taxo_aligns,@otu_taxo_idents,@otu_taxo_dbrefs);
			my @otu_db_matches = split(' \| ',$md5_to_name->{$md5});
			foreach my $otu_db_match (@otu_db_matches) {
				# Ex. AYMZ01000007;tax=k:Bacteria,p:Proteobacteria,c:Gammaproteobacteria,o:Legionellales,f:Coxiellaceae,g:Pseudomonas,s:moraviensis R28-S;aligned=215;ident=212;
				# my %otu_data = read_taxo_data($otu_db_match);
				my ($otu_align_data,$otu_align,$otu_ident);
				if ($otu_db_match =~ /(.+);aligned=(\d+);ident=(\d+);total=(\d+)/) {
					$otu_db_match = $1;
					$otu_align_data = sprintf("%.2f",($3/$4));
# 					$otu_align_data = sprintf("%d/%d",$3,$2);
					$otu_align = sprintf($2);
					$otu_ident = sprintf($3);
				}
				my ($otu_db_ref, $otu_taxo_assign) = split(/\s?[;\|]\s?/,$otu_db_match);
				#if ($otu_db_match =~ /.+;\s?(.+?)[;\s]?$/) {
				push(@otu_taxo_assigns, $otu_taxo_assign);
				push(@otu_taxo_aligns, $otu_align);
				push(@otu_taxo_idents, $otu_ident);
				push(@otu_taxo_dbrefs, $otu_db_ref);
			}
			# Adds taxonomic assignments to variant data
			push(@{$amplisas_results->{$marker_name}{'seq_data'}{$md5}{'taxo_assigns'}}, @otu_taxo_assigns);
			push(@{$amplisas_results->{$marker_name}{'seq_data'}{$md5}{'taxo_aligns'}}, @otu_taxo_aligns);
			push(@{$amplisas_results->{$marker_name}{'seq_data'}{$md5}{'taxo_idents'}}, @otu_taxo_idents);
			push(@{$amplisas_results->{$marker_name}{'seq_data'}{$md5}{'taxo_dbrefs'}}, @otu_taxo_dbrefs);

			
			# Stores info of unified Taxonomies till the desired taxo levels
			foreach my $taxo_level (@taxo_levels){

				# Groups taxonomic assignments into a unique taxonomy with the specified precision level
				$otu_taxo = group_taxo(\@otu_taxo_assigns,$taxo_level);

				push(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$otu_taxo}{'otu_md5s'}}, $md5);
				push(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$otu_taxo}{'otu_idents'}}, @otu_taxo_idents);
				push(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$otu_taxo}{'otu_aligns'}}, @otu_taxo_aligns);
				push(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$otu_taxo}{'otu_dbrefs'}}, @otu_taxo_dbrefs);
				push(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$otu_taxo}{'otu_names'}}, $otu_name);
				push(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$otu_taxo}{'otu_seqs'}}, $otu_seq);

				# Annotates putative ambiguities
				foreach my $otu_taxo_assign (@otu_taxo_assigns){
					$otu_ambiguities->{$taxo_level}{$otu_taxo}{$otu_taxo_assign}++;
				}

				# Stores in $md5_to_taxo associations between OTUs and Taxonomic assignments
				if (!defined($md5_to_taxo->{$md5}{$taxo_level}) || $md5_to_taxo->{$md5}{$taxo_level} eq $otu_taxo){
					$md5_to_taxo->{$md5}{$taxo_level} = $otu_taxo;
				} else {
					printf("\nERROR: OTU '%s' cannot be assigned 2 different taxonomies: '%s' and '%s'.\n\n",$md5,$md5_to_taxo->{$md5}{$taxo_level},$otu_taxo);
					exit;
				}

			}

			# Defines sequence header
			$otu_header = sprintf("%s | %s",$otu_name,$md5_to_name->{$md5});

		} else {
			# All OTUs without matches in the reference DB, will be cathegorized as 'Unassigned'
			push(@{$amplisas_results->{$marker_name}{'seq_data'}{$md5}{'taxo_assigns'}},'Unassigned');
			$otu_taxo = 'tax=Unassigned;';

			foreach my $taxo_level (@taxo_levels){
				push(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$otu_taxo}{'otu_names'}}, $otu_name);
				push(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$otu_taxo}{'otu_seqs'}}, $otu_seq);
				push(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$otu_taxo}{'otu_md5s'}}, $md5);
				# Stores in $md5_to_taxo associations between OTUs and Taxonomic assignments
				if (!defined($md5_to_taxo->{$md5}{$taxo_level}) || $md5_to_taxo->{$md5}{$taxo_level} eq $otu_taxo){
					$md5_to_taxo->{$md5}{$taxo_level} = $otu_taxo;
				} else {
					printf("\nERROR: OTU '%s' cannot be assigned 2 different taxonomies: '%s' and '%s'.\n\n",$md5,$md5_to_taxo->{$md5}{$taxo_level},$otu_taxo);
					exit;
				}
			}

			# Defines sequence header
			$otu_header = sprintf("%s | %s",$otu_name,'Unassigned');

		}
		
		# Stores and formats OTU sequences_hash
		push(@otu_seqs,$otu_seq);
		push(@otu_headers,$otu_header);
		push(@otu_names,$otu_name);

	}
}
# Annotates ambiguities
print "\nAnnotating ambiguous taxonomy assignments.\n";
foreach my $taxo_level (@taxo_levels){
	foreach my $otu_taxo (keys %{$otu_ambiguities->{$taxo_level}}){
		if (scalar keys %{$otu_ambiguities->{$taxo_level}{$otu_taxo}} > 1){
			$amplitaxo_results->{$taxo_level}{'taxo_data'}{$otu_taxo}{'ambiguities'} = [ keys %{$otu_ambiguities->{$taxo_level}{$otu_taxo}} ];
		}
	}
}
# Annotates the Taxonomic abundance sample by sample and prints statistics
print "\nCalculating taxonomy abundance by sample and other statistics.\n";
my (%min_sample_depth,%max_sample_depth);
my $summary_output = "Sample\tTotal depth\tOTUs depth\tUnique OTUs\tTaxonomies\n";
foreach my $taxo_level (@taxo_levels){
	foreach my $sample_name (@unique_samples){
		foreach my $marker_name (keys %{$amplisas_results}){
			foreach my $md5 (@{$amplisas_results->{$marker_name}{'seq_md5s'}}) {
				if (!defined($amplisas_results->{$marker_name}{'assignments'}{$md5}{$sample_name})){
					next;
				}
				my $taxo = $md5_to_taxo->{$md5}{$taxo_level};
				if (!defined($amplitaxo_results->{$taxo_level}{'assignments'}{$taxo}{$sample_name})){
					$amplitaxo_results->{$taxo_level}{'sample_data'}{$sample_name}{'count_taxa'}++;
				}
				$amplitaxo_results->{$taxo_level}{'assignments'}{$taxo}{$sample_name} += $amplisas_results->{$marker_name}{'assignments'}{$md5}{$sample_name};
				$amplitaxo_results->{$taxo_level}{'sample_data'}{$sample_name}{'depth_otus'} += $amplisas_results->{$marker_name}{'assignments'}{$md5}{$sample_name};
				$amplitaxo_results->{$taxo_level}{'sample_data'}{$sample_name}{'count_otus'}++;
			}
			$amplitaxo_results->{$taxo_level}{'sample_data'}{$sample_name}{'depth_total'} += $amplisas_results->{$marker_name}{'sample_data'}{$sample_name}{'depth_amplicon'};
# 			if (!defined($max_sample_depth{$marker_name})){
# 				$max_sample_depth{$marker_name} = $min_sample_depth{$marker_name} = $amplisas_results->{$marker_name}{'sample_data'}{$sample_name}{'depth_amplicon'};
# 			} elsif ($amplisas_results->{$marker_name}{'sample_data'}{$sample_name}{'depth_amplicon'} < $min_sample_depth{$marker_name}){
# 				$min_sample_depth{$marker_name} = $amplisas_results->{$marker_name}{'sample_data'}{$sample_name}{'depth_amplicon'};
# 			} elsif ($amplisas_results->{$marker_name}{'sample_data'}{$sample_name}{'depth_amplicon'} > $max_sample_depth{$marker_name}){
# 				$max_sample_depth{$marker_name} = $amplisas_results->{$marker_name}{'sample_data'}{$sample_name}{'depth_amplicon'};
# 			}
		}
		push(@{$amplitaxo_results->{$taxo_level}{'samples'}}, $sample_name);
		if ($taxo_level eq $paramsdata{'taxo_level'}) {
			$summary_output .= sprintf ("%s\t%d\t%d\t%d\t%d\n",$sample_name,
									$amplitaxo_results->{$taxo_level}{'sample_data'}{$sample_name}{'depth_total'},
									$amplitaxo_results->{$taxo_level}{'sample_data'}{$sample_name}{'depth_otus'},
									$amplitaxo_results->{$taxo_level}{'sample_data'}{$sample_name}{'count_otus'},
									$amplitaxo_results->{$taxo_level}{'sample_data'}{$sample_name}{'count_taxa'},
									);
		}
	}
}

# Decides wich rarefaction depth is optimal for alpha diversity analysis, it's 10 times the average sample depth
if (!defined($paramsdata{'rarefraction_depth'}) || $paramsdata{'rarefraction_depth'} eq 'auto'){
	my @sample_depths;
	foreach my $marker_name (keys %{$amplisas_results}){
		foreach my $sample_name (@unique_samples){
			push(@sample_depths,$amplisas_results->{$marker_name}{'sample_data'}{$sample_name}{'depth_alleles'});
		}
	}
	sprintf("%.0f", mean(@sample_depths)/10) =~ /(\d)(\d+)?/;
	$paramsdata{'rarefraction_depth'} = sprintf("%d%s", $1, 0 x length($2) );
}

# Decides wich rarefaction depth is optimal for beta diversity analysis, it's the closest to the lowest sample depth
my %beta_diversity_depths;
if (defined($paramsdata{'rarefraction_depth'}) && $paramsdata{'rarefraction_depth'} > 1){
	foreach my $marker_name (keys %{$amplisas_results}){
		my $min_sample_depth;
		foreach my $sample_name (@unique_samples){
			if (!defined($min_sample_depth) || $min_sample_depth > $amplisas_results->{$marker_name}{'sample_data'}{$sample_name}{'depth_alleles'}) {
				$min_sample_depth = $amplisas_results->{$marker_name}{'sample_data'}{$sample_name}{'depth_alleles'};
			}
		}
		my $depth = 0;
		while ($depth <= $min_sample_depth){
			$depth += $paramsdata{'rarefraction_depth'};
		}
		$beta_diversity_depths{$marker_name} = $depth;
	}
}

# Performs rarefactions of reads
# A rarefaction is a random collection of sequences from a sample, with a specified depth (number of sequences). 
# The step-size depth is specified by $paramsdata{'rarefraction_depth'}
printf("\nPerforming alpha rarefaction using '%s' estimator in %d steps and %d reads per step.\n", uc($paramsdata{'alpha_diversity'}), $paramsdata{'rarefraction_replicates'}, $paramsdata{'rarefraction_depth'});
my ($rarefaction_results,$beta_diversity_otus);
if (defined($paramsdata{'rarefraction_depth'}) && $paramsdata{'rarefraction_depth'} > 1){
	foreach my $marker_name (keys %{$amplisas_results}){
		foreach my $sample_name (@unique_samples){
			my (@otu_md5s,@cum_depths,$cum_depth);
			foreach my $md5 (@{$amplisas_results->{$marker_name}{'seq_md5s'}}){
				if (defined($amplisas_results->{$marker_name}{'assignments'}{$md5}{$sample_name})){
					push(@otu_md5s, $md5);
					$cum_depth += $amplisas_results->{$marker_name}{'assignments'}{$md5}{$sample_name};
					push(@cum_depths, $cum_depth);
				}
			}
			my $depth = 0;
			while ($depth <= $amplisas_results->{$marker_name}{'sample_data'}{$sample_name}{'depth_alleles'}) {
				$depth += $paramsdata{'rarefraction_depth'};
				# 10 replicates
				for (my $r=0; $r<$paramsdata{'rarefraction_replicates'}; $r++) {
					my %otus;
					for (my $i=0; $i<$depth; $i++){
						my $rand_depth = rand($cum_depth);
						my $prev_depth = 0;
						for (my $j=0; $j<=$#cum_depths; $j++){
							if ($rand_depth>=$prev_depth && $rand_depth<$cum_depths[$j]) {
								$otus{$otu_md5s[$j]}++;
								last;
							}
							$prev_depth = $cum_depths[$j];
						}
					}
					if ($paramsdata{'alpha_diversity'} eq 'chao1') {
						$rarefaction_results->{$paramsdata{'alpha_diversity'}}{$marker_name}{$depth}{$sample_name}[$r] = chao1(%otus);
					} else { # if ($paramsdata{'alpha_diversity'} eq 'otus') {
						$rarefaction_results->{$paramsdata{'alpha_diversity'}}{$marker_name}{$depth}{$sample_name}[$r] = scalar keys %otus;
					}
					if ($beta_diversity_depths{$marker_name} == $depth){
						$beta_diversity_otus->{$marker_name}{$sample_name}[$r] = { %otus };
					}
				}
			}
		}
	}
}

# Print number of sequences and OTUs per sample
print "\nSequences assigned to OTUs:\n$summary_output";

# Print OTU representative sequences into FASTA file:
my $amplitaxo_otuseqs_file = create_fasta_file(\@otu_seqs,\@otu_headers,"$INP_outpath/otu_sequences.fasta");
print "\nOTU sequences stored into '$amplitaxo_otuseqs_file'.\n";

# Print OTU table in BIOM and TXT formats:
my $amplitaxo_otutable_biom = write_amplitaxo_otutable("$INP_outpath/otu_table_depths.biom",$amplisas_results,['biom','depth']);
my $amplitaxo_otutable_txt = write_amplitaxo_otutable("$INP_outpath/otu_table_depths.txt",$amplisas_results,['txt','depth']);
print "\nOTU depth tables stored into '$amplitaxo_otutable_biom' and '$amplitaxo_otutable_txt'.\n";
$amplitaxo_otutable_biom = write_amplitaxo_otutable("$INP_outpath/otu_table_freqs.biom",$amplisas_results,['biom','frequency']);
$amplitaxo_otutable_txt = write_amplitaxo_otutable("$INP_outpath/otu_table_freqs.txt",$amplisas_results,['txt','frequency']);
print "\nOTU frequency tables stored into '$amplitaxo_otutable_biom' and '$amplitaxo_otutable_txt'.\n";

# Cluster by UPGMA OTU sequences according to NEEDLEALL global alignment scores
# my $otu_clusters = cluster_seqs(\@otu_seqs,\@otu_names,'upgma');


# Print genotyping results
my $excelfile_properties = { 
	'title' => "AmpliTAXO results",
	'author' => "Alvaro Sebastian",
	'comments' => "AmpliTAXO results",
	'company' => "Sixth Researcher - www.sixthresearcher.com",
};
my $amplitaxo_outfile = write_amplitaxo_file_results("$INP_outpath/results.xlsx", $amplisas_results, $amplitaxo_results, $rarefaction_results, $excelfile_properties);

# print "\nAnalysis results written into '$INP_output'.\n\n";

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


################################################################################

# Creates a hash from a FASTA file with names as keys and sequences as values
# When headers  are duplicated retrieves not warnings
sub read_fasta_file_hash_ {

	my($fasta_file) = @_;

	my ($sequences,$headers) = read_fasta_file($fasta_file);

	my %sequences_hash;
	for (my $i=0; $i<=$#{$sequences}; $i++) {
		if (!defined($sequences_hash{$headers->[$i]})){
			$sequences_hash{$headers->[$i]} = $sequences->[$i];
		} else {
			my ($header_name, $header_tax);
			if ($headers->[$i] =~ /(.+)(;tax=.+)/){
				$header_name = $1;
				$header_tax = $2;
			} else {
				$header_name = $headers->[$i];
				$header_tax = '';
			}
			my $counter = 0;
			while (defined($sequences_hash{$headers->[$i]})){
				$counter++;
				$headers->[$i] = sprintf("%s\_%d%s", $header_name, $counter, $header_tax);
			}
			$sequences_hash{$headers->[$i]} = $sequences->[$i];
		}
	}

	return \%sequences_hash;
}

#################################################################################

# Writes Ampli results into an Excel file
sub write_amplitaxo_file_results {

	my ($file,$amplisas_results,$amplitaxo_results,$rarefaction_results,$properties) = @_;

	my $workbook  = Excel::Writer::XLSX->new($file);
	$workbook->set_properties(
		title    => $properties->{'title'},
		author   => $properties->{'author'},
		comments => $properties->{'comments'},
		company  => $properties->{'company'}
	);
	$workbook->compatibility_mode();
	my $bold = $workbook->add_format( bold => 1 );
	my $red = $workbook->add_format(bg_color => 'red');
	my $green = $workbook->add_format(bg_color => 'green');
	my $blue = $workbook->add_format(bg_color => 'blue');
	my $yellow = $workbook->add_format(bg_color => 'yellow');
	my $magenta = $workbook->add_format(bg_color => 'magenta');
	my $cyan = $workbook->add_format(bg_color => 'cyan');
	my $decimal2 = $workbook->add_format( num_format => '[=0]0;0.##' );
	my $decimal4 = $workbook->add_format( num_format => '[=0]0;0.####' );

	# Obtains an array with the total unique samples
	my @unique_samples;
	foreach my $marker_name (keys %{$amplisas_results}){
		@unique_samples = unique(@unique_samples,@{$amplisas_results->{$marker_name}{'samples'}});
	}
	@unique_samples = nsort(@unique_samples);

	# Creates a worksheet per each marker and writes amplicon information
	# Like AmpliSAS format with additional Taxonomic assignments
	foreach my $marker_name (keys %{$amplisas_results}){

		my @seq_md5s = @{$amplisas_results->{$marker_name}{'seq_md5s'}};

		# NEW WORKSHEET

		# Creates worksheet and writes data headers
		# Sheetname must be <= 31 chars
		my $worksheet = $workbook->add_worksheet("OTUs_".substr($marker_name, 0, 26));
		$worksheet->set_column('A:A', 12, $bold); # OTU names column
		$worksheet->set_column('I:I', 75, $bold); # Taxonomy column
		my $ws_row = 0;
		my @seq_data_headers = ('OTU_ID', 'SEQUENCE', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'DB_REFERENCES', 'IDENTITIES', 'TAXONOMY');
		$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_AMPLICON'); $ws_row++;
		$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_OTUS'); $ws_row++;
		$worksheet->write($ws_row, $#seq_data_headers, 'COUNT_OTUS'); $ws_row++;
		$worksheet->write_row($ws_row, 0, \@seq_data_headers, $bold); $ws_row++;

		my $ws_row_first = $ws_row;
		my $ws_row_last = $ws_row_first+$#seq_md5s;
		my $ws_col = $#seq_data_headers;
		my $ws_col_first = $ws_col+1;

		my (%final_seq_depths,%final_seq_samples);
		foreach my $sample (@unique_samples){
			$ws_row = 0; $ws_col++;
			$worksheet->write($ws_row, $ws_col, $amplisas_results->{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'}); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $amplisas_results->{$marker_name}{'sample_data'}{$sample}{'depth_alleles'}); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $amplisas_results->{$marker_name}{'sample_data'}{$sample}{'count_alleles'}); $ws_row++;
			$worksheet->write($ws_row, $ws_col, $sample, $bold); $ws_row++;
			foreach my $seq_md5 (@seq_md5s) {
				if (defined($amplisas_results->{$marker_name}{'assignments'}{$seq_md5}{$sample})){
					$worksheet->write($ws_row, $ws_col, $amplisas_results->{$marker_name}{'assignments'}{$seq_md5}{$sample}); $ws_row++;
					$final_seq_depths{$seq_md5} += $amplisas_results->{$marker_name}{'assignments'}{$seq_md5}{$sample};
					$final_seq_samples{$seq_md5}++;
				} else {
					$worksheet->write($ws_row, $ws_col, ' '); $ws_row++;
				}
			}
		}

		my $ws_col_last = $ws_col;

		# Writes seq information
		$ws_row = $ws_row_first;
		foreach my $seq_md5 (@seq_md5s) {
			$ws_col = 0;
			$worksheet->write($ws_row, $ws_col,$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'name'}); $ws_col++;
			$worksheet->write($ws_row, $ws_col,$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'sequence'}); $ws_col++;
			$worksheet->write($ws_row, $ws_col,$seq_md5); $ws_col++;
			$worksheet->write($ws_row, $ws_col,length($amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'sequence'})); $ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $final_seq_depths{$seq_md5}); $ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $final_seq_samples{$seq_md5}); $ws_col++;
			if (defined($amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_assigns'}) && $amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_assigns'}[0] ne 'Unassigned'){
				$worksheet->write($ws_row, $ws_col,join(';',@{$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_dbrefs'}})); $ws_col++;
				my @sum_idents = map sprintf("%d/%d",$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_idents'}[$_],$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_aligns'}[$_]) , 0..$#{$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_aligns'}};
				my $sum_idents = sprintf("(%s)/%d",join('+',@sum_idents),scalar @{$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_aligns'}});
				$worksheet->write_formula($ws_row, $ws_col, $sum_idents, $decimal2, sprintf("%.4f",eval($sum_idents))); $ws_col++;
# 				$worksheet->write_formula($ws_row, $ws_col,sprintf("(%s)/%d",join('+',@{$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_idents'}}),scalar @{$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_idents'}}),undef,mean(@{$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_idents'}})); $ws_col++;
				$worksheet->write($ws_row, $ws_col,join(' | ',unique(@{$amplisas_results->{$marker_name}{'seq_data'}{$seq_md5}{'taxo_assigns'}}))); $ws_col++;
			} else {
				$worksheet->write($ws_row, $ws_col, ' '); $ws_col++;
				$worksheet->write($ws_row, $ws_col, ' '); $ws_col++;
				$worksheet->write($ws_row, $ws_col, 'Unassigned'); $ws_col++;
			}
			$ws_row++;
		}

	}

	# Creates a worksheet for each desired Taxonomy level with Taxonomy summary (all markers and OTUs merged in one page)
	foreach my $taxo_level (keys %$amplitaxo_results) {

		# NEW WORKSHEET

		# Creates worksheet and writes data headers
		my $worksheet_name = "Taxo_$taxo_level";
		my $worksheet = $workbook->add_worksheet($worksheet_name);
		if ($taxo_level eq 'class') {
			$worksheet->set_column('F:F', 50, $bold); # Taxonomy column
		} else {
			$worksheet->set_column('F:F', 75, $bold); # Taxonomy column
		}
		my $ws_row = 30;
		my @seq_data_headers = ('OTU_IDs', 'MEAN_IDENT', 'DEPTH', 'SAMPLES', 'MEAN_FREQ', 'TAXONOMY'); #('SEQUENCES', 'MD5S', 'NAMES', 'MEAN_IDENT', 'DEPTH', 'SAMPLES', 'MEAN_FREQ', 'TAXONOMY');
		$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH', $bold);
		$ws_row++;
		$worksheet->write($ws_row, $#seq_data_headers, 'COUNT', $bold);
		$ws_row++;
		$worksheet->write_row($ws_row, 0, \@seq_data_headers, $bold);
		$ws_row++;

		# Obtains all the taxa from all the samples to sorts them by depth
		my (%taxo_depths, %taxo_freqs, %taxo_samples);
		foreach my $taxo (keys %{$amplitaxo_results->{$taxo_level}{'assignments'}}) {
			foreach my $sample (keys %{$amplitaxo_results->{$taxo_level}{'assignments'}{$taxo}}) {
				$taxo_samples{$taxo}++;
				$taxo_depths{$taxo} += $amplitaxo_results->{$taxo_level}{'assignments'}{$taxo}{$sample};
				$taxo_freqs{$taxo} += $amplitaxo_results->{$taxo_level}{'assignments'}{$taxo}{$sample}/$amplitaxo_results->{$taxo_level}{'sample_data'}{$sample}{'depth_total'};
			}
			$taxo_freqs{$taxo} = $taxo_freqs{$taxo}/$taxo_samples{$taxo};
		}
		my @taxa = sort { $taxo_depths{$b} <=> $taxo_depths{$a} } keys %taxo_depths;

		my $ws_row_first = $ws_row;
		my $ws_row_last = $ws_row_first+$#taxa;
		my $ws_col = $#seq_data_headers;
		my $ws_col_first = $ws_col+1;
		my $ws_col_last = $ws_col+scalar(@unique_samples);
		my $chart_series;
		$chart_series->{'all'}{'categories'}{'first'} = xl_col_to_name($ws_col_first).($ws_row_first);
		$chart_series->{'all'}{'categories'}{'last'} = xl_col_to_name($ws_col_last).($ws_row_first);

		foreach my $sample (@unique_samples){
			$ws_row = 30; $ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $amplitaxo_results->{$taxo_level}{'sample_data'}{'depth_otus'});
			$ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $amplitaxo_results->{$taxo_level}{'sample_data'}{'count_otus'});
			$ws_row++;
			$worksheet->write($ws_row, $ws_col, $sample, $bold);
			$ws_row++;
			foreach my $taxo (@taxa) {
				if (defined($amplitaxo_results->{$taxo_level}{'assignments'}{$taxo}{$sample})){
					$worksheet->write($ws_row, $ws_col, $amplitaxo_results->{$taxo_level}{'assignments'}{$taxo}{$sample});
	# 				$worksheet->write($ws_row, $ws_col, sprintf("%.4f (%d)", $amplitaxo_results->{$taxo_level}{'assignments'}{$taxo}{$sample}/$amplitaxo_results->{$taxo_level}{'sample_data'}{$sample}{'depth_total'}, $amplitaxo_results->{$taxo_level}{'assignments'}{$taxo}{$sample}));
					$ws_row++;
				} else {
					$worksheet->write($ws_row, $ws_col, ' ');
					$ws_row++;
				}
				if ($sample eq $unique_samples[0]){
					$chart_series->{$taxo}{'values'}{'first'} = xl_col_to_name($ws_col).($ws_row);
				} 
				if ($sample eq $unique_samples[-1]){
					$chart_series->{$taxo}{'values'}{'last'} = xl_col_to_name($ws_col).($ws_row);
				}
			}
		}

	# 	my $ws_col_last = $ws_col;

		# Writes seq information
		$ws_row = $ws_row_first;
		foreach my $taxo (@taxa) {
			$ws_col = 0;
			# Each Taxonomy can have several OTUs/representative sequences
			$worksheet->write($ws_row, $ws_col, join(';',@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'otu_names'}}));
			$ws_col++;
			if (defined($amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'otu_idents'})){
				my @sum_idents = map sprintf("%d/%d",$amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'otu_idents'}[$_],$amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'otu_aligns'}[$_]) , 0..$#{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'otu_aligns'}};
				my $sum_idents = sprintf("(%s)/%d",join('+',@sum_idents),scalar @{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'otu_aligns'}});
				$worksheet->write($ws_row, $ws_col, sprintf("%.4f",eval($sum_idents)), $decimal2);
				# $worksheet->write_formula($ws_row, $ws_col,sprintf("(%s)/%d",join('+',@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'otu_idents'}}),scalar @{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'otu_idents'}}),$decimal2,mean(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'otu_idents'}}));
	# 			$worksheet->write($ws_row, $ws_col,sprintf("%.2f",mean(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'otu_idents'}})),$decimal2);
			} else {
				$worksheet->write($ws_row, $ws_col,' ',$decimal2);
			}
			$ws_col++;
	# 		$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $taxo_depths{$taxo});
			$worksheet->write($ws_row, $ws_col, $taxo_depths{$taxo});
			$ws_col++;
	# 		$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $taxo_samples{$taxo});
			$worksheet->write($ws_row, $ws_col, $taxo_samples{$taxo});
			$ws_col++;
			$worksheet->write($ws_row, $ws_col, sprintf("%.4f",$taxo_freqs{$taxo}), $decimal4);
			$ws_col++;
			$worksheet->write($ws_row, $ws_col,$taxo, $bold);
			$chart_series->{$taxo}{'name'} = xl_col_to_name($ws_col+1).($ws_row+1);
			$ws_col++;
			$ws_row++;
		}
		$ws_col--;
		$ws_row++;
		$worksheet->write($ws_row, $ws_col,"AMBIGUOUS TAXONOMY", $bold);
		$worksheet->write($ws_row, $ws_col+1,"DISAMBIGUATION", $bold);
		$ws_row++;
		foreach my $taxo (@taxa) {
			if (defined($amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'ambiguities'})){
				$worksheet->write($ws_row, $ws_col, $taxo, $bold);
				$worksheet->write($ws_row, $ws_col+1, join('; ', unique(@{$amplitaxo_results->{$taxo_level}{'taxo_data'}{$taxo}{'ambiguities'}})));
				$ws_row++;
			}
		}
		
		# Creates a new chart object. In this case an embedded chart.
		my ($chart_col, $chart_row) = (0, 1);
		my $chart = $workbook->add_chart( type => 'column', subtype => 'percent_stacked', embedded => 1 );
		# Ads a chart title and some axis labels.
		$chart->set_title ( name => "Taxonomy abundance" , name_font => { name => 'Arial', size => 16 } , num_font => { bold => 1 } );
		$chart->set_x_axis( name => "Samples" , name_font => { name => 'Arial', size => 12 } );
		$chart->set_y_axis( name => "Relative abundance" , name_font => { name => 'Arial', size => 12 } );
		$chart->set_legend( position => 'right' , font => { name => 'Arial', size => 8 } );
		# Sets an Excel chart style. Colors with white outline and shadow.
		$chart->set_style( 2 );
		# Configures the series.
		foreach my $taxo (reverse @taxa) {
			$chart->add_series(
				name       => $taxo, #"=Taxonomies!".$chart_series->{$taxo}{'name'},
				categories => "=$worksheet_name!".$chart_series->{'all'}{'categories'}{'first'}.":".$chart_series->{'all'}{'categories'}{'last'},
				values => "=$worksheet_name!".$chart_series->{$taxo}{'values'}{'first'}.":".$chart_series->{$taxo}{'values'}{'last'},
			);
		}
		# Inserts the chart into the worksheet (with an offset).
		$worksheet->insert_chart( xl_col_to_name($chart_col).$chart_row, $chart, 0, 0, 2, 2 );
	}
	
	# NEW WORKSHEET

	# Performs alpha diversity analysis and draws rarefaction plots
	# $rarefaction_results->{$alpha_diversity}{$marker_name}{$depth}{$sample_name}[$r]{$md5} = $depth;
	if (defined($rarefaction_results)){
		foreach my $alpha_diversity (keys %{$rarefaction_results}){
			foreach my $marker_name (keys %{$rarefaction_results->{$alpha_diversity}}){
				# Creates worksheet and writes data headers
				my $worksheet_name = "Alpha_$marker_name\_$alpha_diversity";
				# Worksheet names longer than 31 chars give error
				if (length($worksheet_name)>31){
					$worksheet_name = "Alpha_$alpha_diversity";
				}
				my $worksheet = $workbook->add_worksheet($worksheet_name);
				my $ws_row = 30;
				my $ws_col = 0;
				my $ws_col_first = $ws_col;
				my $chart_series;
				my $sd_row_first;
				# Prints mean and standard deviation values
				foreach my $measure (('MEAN', 'SD')){
					$ws_col = $ws_col_first;
					if ($measure eq 'MEAN'){
						$worksheet->write($ws_row, $ws_col, sprintf("Alpha diversity measure: %s average (%d replicates)", $ALPHA_DIVERSITY{$alpha_diversity}, $paramsdata{'rarefraction_replicates'}), $bold);
					} elsif ($measure eq 'SD'){
						$worksheet->write($ws_row, $ws_col, sprintf("Alpha diversity measure: %s standard deviation (%d replicates)", $ALPHA_DIVERSITY{$alpha_diversity}, $paramsdata{'rarefraction_replicates'}), $bold);
					}
					$ws_row+=2;
					my @seq_data_headers = ('DEPTH', @unique_samples);
					$worksheet->write_row($ws_row, $ws_col, \@seq_data_headers, $bold);
					$ws_row++;
					my $ws_row_first = $ws_row;
					foreach my $depth (sort {$a<=>$b} keys %{$rarefaction_results->{$alpha_diversity}{$marker_name}}){
						$ws_col = $ws_col_first;
						$worksheet->write($ws_row, $ws_col, $depth);$ws_col++;
						foreach my $sample_name (@unique_samples){
							if (defined($rarefaction_results->{$alpha_diversity}{$marker_name}{$depth}{$sample_name})){
								my @rarefaction_values;
								for (my $r=0; $r<=$#{$rarefaction_results->{$alpha_diversity}{$marker_name}{$depth}{$sample_name}}; $r++){
									push(@rarefaction_values,$rarefaction_results->{$alpha_diversity}{$marker_name}{$depth}{$sample_name}[$r]);
								}
								if ($measure eq 'MEAN'){
									$worksheet->write($ws_row, $ws_col, sprintf("%.2f",mean(@rarefaction_values)), $decimal2);
								} elsif ($measure eq 'SD'){
									$worksheet->write($ws_row, $ws_col, sprintf("%.4f",stdeviation(@rarefaction_values)), $decimal2);
								}
							}
							$ws_col++;
						}
						$ws_row++
					}
					my $ws_col_ = $ws_col_first + 1;
					if ($measure eq 'MEAN'){
						$chart_series->{'all'}{'categories'}{'first'} = xl_col_to_name($ws_col_first).($ws_row_first+1);
						$chart_series->{'all'}{'categories'}{'last'} = xl_col_to_name($ws_col_first).($ws_row);
						foreach my $sample_name (@unique_samples){
							$chart_series->{$sample_name}{'values'}{'first'} = xl_col_to_name($ws_col_).($ws_row_first+1);
							$chart_series->{$sample_name}{'values'}{'last'} = xl_col_to_name($ws_col_).($ws_row);
							$ws_col_++;
						}
					} elsif ($measure eq 'SD'){
						foreach my $sample_name (@unique_samples){
							$chart_series->{$sample_name}{'errors'}{'first'} = xl_col_to_name($ws_col_).($ws_row_first+1);
							$chart_series->{$sample_name}{'errors'}{'last'} = xl_col_to_name($ws_col_).($ws_row);
							$ws_col_++;
						}
					}
					$ws_row+=2;
				}
				# Creates a new chart object. In this case an embedded chart.
				my ($chart_col, $chart_row) = (0, 1);
				my $chart = $workbook->add_chart( type => 'line', embedded => 1 );
				# Ads a chart title and some axis labels.
				$chart->set_title ( name => "Alpha diversity" , name_font => { name => 'Arial', size => 16 } , num_font => { bold => 1 } );
				$chart->set_x_axis( name => "Sequences per sample" , name_font => { name => 'Arial', size => 12 } );
				$chart->set_y_axis( name => $ALPHA_DIVERSITY{$alpha_diversity} , name_font => { name => 'Arial', size => 12 } );
				$chart->set_legend( position => 'right' , font => { name => 'Arial', size => 8 } );
				# Sets an Excel chart style. Colors with white outline and shadow.
				$chart->set_style( 2 );
				# Configures the series.
				foreach my $sample_name (@unique_samples){
					$chart->add_series(
						name       => $sample_name,
						categories => "=$worksheet_name!".$chart_series->{'all'}{'categories'}{'first'}.":".$chart_series->{'all'}{'categories'}{'last'},
						values => "=$worksheet_name!".$chart_series->{$sample_name}{'values'}{'first'}.":".$chart_series->{$sample_name}{'values'}{'last'},
						smooth => 1 ,
						y_error_bars => {
							type         => 'custom',
							plus_values  => "=$worksheet_name!".$chart_series->{$sample_name}{'errors'}{'first'}.":".$chart_series->{$sample_name}{'errors'}{'last'},
							minus_values => "=$worksheet_name!".$chart_series->{$sample_name}{'errors'}{'first'}.":".$chart_series->{$sample_name}{'errors'}{'last'},
						},
					);
				}
				# Inserts the chart into the worksheet (with an offset).
				$worksheet->insert_chart( xl_col_to_name($chart_col).$chart_row, $chart, 0, 0, 2, 2 );
			}
		}
	}

	$workbook->close();
	
	return $file;

}

################################################################################


