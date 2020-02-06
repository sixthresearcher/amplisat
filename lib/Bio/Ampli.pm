package Bio::Ampli;

######################################################################################################
#
# Bio::Ampli - Package of subroutines to work with Ampplicon Sequencing data
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
# Package of subroutines to work with Ampplicon Sequencing data
#
######################################################################################################

# # Returns full path and file name
# sub dirname {
# 	my $path = shift;
# 	if ($path =~ m{^(.*/)?.*}s){
# 		if (defined($1)) {
# 			return $1;
# 		} else {
# 			return '';
# 		}
# 	} else {
# 		print "\nERROR: The path of the script couldn't be found, run it in a Linux system.\n\n";
# 		exit;
# 	}
# }
# # Libraries are in folder '../' in the path of the script
# use lib dirname(__FILE__).'../';

# Modules are in folder '../' in the path of the script
use File::FindLib '../';
# Perl modules necessaries for the correct working of the script
use List::Util qw(shuffle);
use Text::Iconv;
use Spreadsheet::XLSX;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Sort::Naturally;
use threads;
use threads::shared;
use Bio::Sequences;
use 5.004;
use strict;
use warnings;
no warnings ('uninitialized', 'substr');
use Time::HiRes qw(gettimeofday);
use File::Basename;
# use Data::Dumper;

# Turn autoflush on
local $| = 1;


use Exporter;
 
our @ISA       = qw(Exporter);
our $VERSION   = '1.0';
our @EXPORT = qw(
parse_sequence_file
parse_amplicon_file
read_amplicon_data
print_amplicon_data
extract_primer_tag_seqs
extract_primer_tag_iupac_seqs
align_amplicons
match_amplicons
match_amplicons_regex
match_amplicons_regex_with_threads
find_amplicon_reads
find_amplicon_reads_with_threads
read_allele_file
read_alleles
match_alleles
retrieve_amplihla_alleles
retrieve_amplicancer_refs
retrieve_amplicon_data
print_marker_sequences
print_amplicon_sequences
print_comparison_sequences
annotate_low_depth
filter_amplicon_sequences
filter_amplicon_sequences_with_threads
cluster_amplicon_sequences
cluster_amplicon_sequences_with_threads
cluster_denovo_amplicon_sequences
cluster_denovo_amplicon_sequences_with_threads
cluster_denovo_sequence_file
cluster_denovo_sequences
cluster_denovo_sequences_old
is_dominant
compare_amplicon_sequences
compare_amplicon_sequences_with_threads
genotype_amplicon_sequences
genotype_amplicon_sequences_with_threads
extract_alleles_freq_threshold
read_compressed_file_amplicons
read_amplisas_file_amplicons
read_amplisas_file_results
write_amplisas_file_results
write_amplicdr3_file_results
write_amplihla_file_results
write_amplicancer_file_results
write_amplitaxo_otutable
create_amplitaxo_otutable
degree_of_change
generate_barcode
adjust_automatic_params

);

######################################################################################################

# # Routes to binary files and databases
# my $TOOLSDIR = dirname (__FILE__).'/tools/';
# my $NEEDLEALLEXE = $TOOLSDIR.'needleall';
# my $NEEDLEMANWUNSCHEXE = $TOOLSDIR.'needleman_wunsch';
# my $SMITHWATERMANEXE = $TOOLSDIR.'smith_waterman';
# my $MAFFTEXE = $TOOLSDIR.'mafft';

#################################################################################

# Parses reads file from command input in amplicon sequencing scripts
sub parse_sequence_file {

	my ($seqs_file,$seq_number,$options) = @_;
	
	
	if (!-e $seqs_file){
		print "\nERROR: Sequences file '$seqs_file' doesn't exist.\n\n";
		exit;
	}

	my ($read_qualities,$only_stats,$verbose,$shuffle)=(0,0,0,0);
	if (in_array($options, 'qualities')){
		$read_qualities = 1;
	}
	if (in_array($options, 'verbose')){
		$verbose = 1;
	}
	if (in_array($options, 'stats')){
		$only_stats = 1;
	}
	if (in_array($options, 'shuffle')){
		$shuffle = 1;
	}
	# my $outdir = dirname($outfile);
	# my $outname = basename($outfile);
	# if (!defined(outdir)){
	# 	$outdir = '.';
	# }
	# #print "outdir, outname";

	if ($verbose) {
		print "\nChecking input sequence file '$seqs_file'.\n";
	}

	my $seqs_file_format;
	if (is_fastq($seqs_file)){
		$seqs_file_format = 'fastq';
		if ($verbose) {
			print "\tSequences are in FASTQ format.\n";
		}
	} elsif (is_fasta($seqs_file)){
		$seqs_file_format = 'fasta';
		if ($verbose) {
			print "\tSequences are in FASTA format.\n";
		}
	} else {
		print "\nERROR: Sequence file '$seqs_file' must be FASTA or FASTQ format (compressed or uncompressed).\n\n";
		exit;
	}

	my $total_reads;
	if ($seqs_file_format eq 'fastq'){
		$total_reads = count_seqs_from_fastq($seqs_file);
	} elsif ($seqs_file_format eq 'fasta'){
		$total_reads = count_seqs_from_fasta($seqs_file);
	}
	if ($verbose) {
		print "\tSequences number: $total_reads.\n";
	}
	
	if ($only_stats){
		return ($seqs_file_format,$total_reads);
	}
	
	# Reads input sequences
	my ($seqs,$headers,$qualities);
	if (!defined($seq_number)){
		if ($verbose) {
			if (!$shuffle){
				print "\nReading sequence data.\n";
			} else {
				print "\nReading shuffled sequence data.\n";
			}
		}
		if ($seqs_file_format eq 'fastq'){
# 			($seqs,$headers) = read_fastq_file($seqs_file); # 1=Read qualities
			($seqs,$headers,$qualities) = read_fastq_file($seqs_file,$read_qualities); # 1=Read qualities
		} elsif ($seqs_file_format eq 'fasta'){
			($seqs,$headers) = read_fasta_file($seqs_file);
		}
		if ($shuffle){
			($headers,$seqs,$qualities) = shuffle_seqs($headers,$seqs,$qualities);
		}
	} else {
		if ($verbose) {
			print "\nReading $seq_number sequences from randomized data.\n";
		}
		if ($seqs_file_format eq 'fastq'){
			($seqs,$headers,$qualities) = extract_random_seqs_from_fastq($seqs_file, $seq_number, $read_qualities);
		} elsif ($seqs_file_format eq 'fasta'){
			($seqs,$headers) = extract_random_seqs_from_fasta($seqs_file, $seq_number);
		}
	}
	
	return ($seqs_file_format,$seqs,$headers,$qualities,$total_reads);

}

#################################################################################

# Parses amplicons file from command input in amplicon sequencing scripts
sub parse_amplicon_file {

	# 2 PRIMERS => MARKER
	# 1/2 TAGS => SAMPLE
	# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

	my ($amplicons_file, $options) = @_;

	my ($skip_errors, $verbose, $skip_markers, $skip_samples) = (0,0,0,0);
	if (in_array($options, 'skip errors')){
		$skip_errors = 1;
	}
	if (in_array($options, 'verbose')){
		$verbose = 1;
	}
	if (in_array($options, 'skip markers')){
		$skip_markers = 1;
	}
	if (in_array($options, 'skip samples')){
		$skip_samples = 1;
	}

	# Reads markers/primers data from CVS input files
	if ($verbose) {
		print "\nReading amplicon data from file '$amplicons_file'.\n";
	}
	my ($markerdata, $markers) = read_amplicon_data($amplicons_file, 'markers', $options);
	if (!$skip_markers && !$skip_errors && !defined($markerdata)){
		print "\nERROR: Amplicon data file has not the correct format or marker data is missing.\n\n";
		exit;
	}
	# print Dumper($markerdata);
	# exit;
	if (!$skip_markers && $verbose) {
		print "\tNumber of markers: ".scalar @{$markers}.".\n";
	}

# # 	if (defined($align_type)){
# 	my ($primer_seqs, $primer_headers) = extract_primer_tag_seqs($markerdata, $markers);
# 	if ($verbose) {
# 		print "\tNumber of unique primer sequences: ".scalar @{$primer_headers}.".\n";
# 	}
# # 	} else {
# # 		($primer_seqs, $primer_headers) = extract_primer_tag_iupac_seqs($markerdata, $markers);
# # 	}
# 
# 	# # Creates a FASTA file with primer sequences
# 	# print "\tStoring primer sequences into file '$outfile.primers.fa'.\n";
# 	# create_fasta_file($primer_seqs,$primer_headers,"$outfile.primers.fa");
# 	# Stores primers into a hash and checks that there are not name duplications
# 	my %primer_seqs_hash;
# 	map $primer_seqs_hash{$primer_headers->[$_]} = $primer_seqs->[$_], 0 .. $#{$primer_headers};
# 	if (!$skip_errors && scalar @{$primer_seqs} != scalar keys %primer_seqs_hash ){
# 		print "\nERROR: Primers cannot have duplicated names.\n\n";
# 		exit;
# 	}

	# Reads only barcodes data from CVS input file
	my ($sampledata, $samples) = read_amplicon_data($amplicons_file, 'samples', $options);
	if (!$skip_samples && !$skip_errors && !defined($sampledata)) {
		print "\nERROR: Amplicon data file has not the correct format or sample barcodes are missing in '$amplicons_file' file.\n\n";
		exit;
	}
	# print Dumper($sampledata);
	# exit;
	if (!$skip_samples && $verbose) {
		print "\tNumber of samples: ".scalar @{$samples}.".\n";
	}

# 	# Reads amplicon/primers with barcodes data from CVS input files
# 	# my ($amplicondata, $amplicons) = read_amplicon_data($amplicons_file, 'amplicons', $options);
# 	# print Dumper($amplicondata);
# 	# exit;


# # 	if (defined($align_type)){
# 	my ($primer_tag_seqs, $primer_tag_headers) = extract_primer_tag_seqs($amplicondata, $amplicons);
# 	if ($verbose) {
# 		print "\tNumber of unique primer+barcode sequences: ".scalar @{$primer_tag_headers}.".\n";
# 	}
# # 	} else {
# # 		($primer_tag_seqs, $primer_tag_headers) = extract_primer_tag_iupac_seqs($amplicondata, $amplicons);
# # 	}
# 
# 	# # Creates a FASTA file with primer/barcode sequences
# 	# print "\tStoring primer+barcode sequences into file '$outfile.primers+barcodes.fa'.\n";
# 	# create_fasta_file($primer_tag_seqs,$primer_tag_headers,"$outfile.primers+barcodes.fa");
# 	# Stores barcodes into a hash and checks that there are not name duplications
# 	my %tag_seqs_hash;
# 	map $tag_seqs_hash{$_} = $sampledata->{$_}{'tag_f'}.'...'.$sampledata->{$_}{'tag_rc'}, @{$samples};
# 	if (!$skip_errors && scalar @{$samples} != scalar keys %tag_seqs_hash ){
# 		print "\nERROR: Tags cannot have duplicated names.\n\n";
# 		exit;
# 	}
	
	# Reads allele data from CVS input file
	my $alleledata = read_amplicon_data($amplicons_file, 'alleles', $options);
	if ($verbose && defined($alleledata) && %$alleledata){
		print "\tNumber of alleles: ".scalar(keys %{$alleledata}).".\n";
	}

	# Reads filters and clustering thresholds data from CVS input file
# 	my $paramsdata;
	my $paramsdata = read_amplicon_data($amplicons_file, 'params', $options);
# 	if (!defined($paramsdata_)) {
# 		# print "\nERROR: Amplicon data file has not the correct format or there are no clustering/filtering parameters in '$amplicons_file' file.\n\n";
# 	} else {
# 		foreach my $filtername (keys %{$paramsdata_}){
# 			if ($filtername ne 'allowed_markers') {
# 				foreach my $markername (keys %{$paramsdata_->{$filtername}}){
# 					foreach my $value (@{$paramsdata_->{$filtername}{$markername}}){
# 						if (!in_array($paramsdata->{$filtername}{$markername}, $value)) {
# 							push(@{$paramsdata->{$filtername}{$markername}}, $value);
# 						} elsif (!$skip_errors) {
# 							print "\nERROR: Parameter '$filtername/$markername' has duplicated value '$value'.\n\n";
# 							exit;
# 						}
# 					}
# 				}
# 			} else {
# 				push(@{$paramsdata->{$filtername}}, @{$paramsdata_->{$filtername}});
# 			}
# 		}
# 	}

	return ($markerdata,$markers,$sampledata,$samples,$paramsdata,$alleledata);

}

#################################################################################

# Reads amplicon data from .csv file
sub read_amplicon_data {

	# 2 PRIMERS => MARKER
	# 1/2 TAGS => SAMPLE
	# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

	my ($file, $data_type, $options) = @_;

	my $skip_errors = 0;
	if (in_array($options, 'skip errors')){
		$skip_errors = 1;
	}
	my $only_params = 0;
	if (in_array($options, 'only params')){
		$only_params = 1;
	}

	my ($sampledata, $markerdata, $amplicondata, $paramsdata, $alleledata);
	my @samples; # To keep the original sample order
	my @markers; # To keep the original marker/primer order
	my @amplicons;
	my @fields;
	my $cvsdataformat = '';
	my %sampleseqs; # To avoid repeated sample barcode sequences
	my %markerseqs; # To avoid repeated marker primer sequences
	
	# Example of CVS file:
	# marker,primer_f,primer_r,gene,feature,specie,length
	# MHC1E2,GsTGCTCCTrCTGCTGGC,CCTCGCTCTGGTTGTAGT,MHC class I,exon2,Myodes glareolus,317
	# MHC1E3,ACTACAACCAGAGCGAGG,TGTGCCTTTGGGsGAwCT,MHC class I,exon3,Myodes glareolus,313
	# MHC2DQA,rTCCTCGCCCTGACCACC,GGGTGTTGGGCTGACCCA,MHC class II DQA,exon2,Myodes glareolus,364
	# MHC2DQB,AGCTGTGGTGCTGATGGT,TCrAGCCGCCGCAGGGAA,MHC class II DQB,exon2,Myodes glareolus,330
	# MHC2DRB,TGGCAGCTGTGATCCTGA,AGCAGACCAGGAGGTTGT,MHC class II DRB,exon2,Myodes glareolus,405
	# sample_name,tag_f,tag_r
	# BV001a,AGAGAC,TGTACA
	# BV001b,CACAGT,CGTCAC
	# BV002a,AGAGAC,TGATCC
	#
	# OR:
	#
	# sample;tag_f;primer_f;tag_r;primer_r;marker;species;length
	# BV001a-MHC1E2;AGAGAC;GsTGCTCCTrCTGCTGGC;TGTACA;CCTCGCTCTGGTTGTAGT;MHC class I exon 2;Myodes glareolus;317
	# BV001a-MHC1E3;AGAGAC;ACTACAACCAGAGCGAGG;TGTACA;TGTGCCTTTGGGsGAwCT;MHC class I exon 3;Myodes glareolus;313
	# BV001b-MHC1E2;CACAGT;GsTGCTCCTrCTGCTGGC;CGTCAC;CCTCGCTCTGGTTGTAGT;MHC class I exon 2;Myodes glareolus;317
	# BV001b-MHC1E3;CACAGT;ACTACAACCAGAGCGAGG;CGTCAC;TGTGCCTTTGGGsGAwCT;MHC class I exon 3;Myodes glareolus;313
	# BV002a-MHC1E2;AGAGAC;GsTGCTCCTrCTGCTGGC;TGATCC;CCTCGCTCTGGTTGTAGT;MHC class I exon 2;Myodes glareolus;317
	# BV002a-MHC1E3;AGAGAC;ACTACAACCAGAGCGAGG;TGATCC;TGTGCCTTTGGGsGAwCT;MHC class I exon 3;Myodes glareolus;313

	
	open(INFILE,"$file")|| die "# $0 : cannot open $file\n";

	while (my $line = <INFILE>){
		if ($line =~ /^#/) { next; }
		if ($line =~ /^\s*$/) { next; }
		$line = trim($line);
		my @values = split(/[,|;|\t]/,$line,-1); # -1 forces to split also empty fields at the end (ex. MHCII,,)
		# Skip lines with less than 2 columns
		if (scalar @values < 2) { next; }
		# Checks the kind of info to annotate
		if ($line =~ /^>(.+)/){
			$line =~ s/\s+//g; # trim spaces
			@fields = split(/[,|;|\t]/,$1,-1);
			if ($fields[0] =~ /marker|primer/) {
				$cvsdataformat = 'markers';
			} elsif ($fields[0] =~ /sample|tag|barcode/){
				$cvsdataformat = 'samples';
			} elsif ($fields[0] =~ /param|filter|threshold/) {
				$cvsdataformat = 'params';
			} elsif ($fields[0] =~ /allele/) {
				$cvsdataformat = 'alleles';
			}
		# Annotates sample type data:
		# sample_name,tag_f,tag_r
		# BV001a,AGAGAC,TGTACA
		# BV001b,CACAGT,CGTCAC
		#
		# OR:
		#
		# sample_name;tag_f;primer_f;tag_r;primer_r;marker;species;length
		# BV001a-MHC1E2;AGAGAC;GsTGCTCCTrCTGCTGGC;TGTACA;CCTCGCTCTGGTTGTAGT;MHC class I exon 2;Myodes glareolus;317
		# BV001a-MHC1E3;AGAGAC;ACTACAACCAGAGCGAGG;TGTACA;TGTGCCTTTGGGsGAwCT;MHC class I exon 3;Myodes glareolus;313
		} elsif ($cvsdataformat eq 'samples' && (!defined($data_type) || $data_type eq 'samples' || $data_type eq 'amplicons') ) {
			my $samplename = convert_to_alphanum($values[0]);
			# Checks if sample names are duplicated
			if (!$skip_errors && defined($sampledata->{$samplename})) {
				#printf("\nWARNING: Sample '%s' is duplicated in CVS datafile.\n", $samplename);
				printf("\nERROR: Sample '%s' is duplicated in CVS datafile (samples must have different names).\n", $samplename);
				exit;
			}
			push(@samples, $samplename);
			# Reads fields
			for (my $i=0; $i<=$#fields; $i++){
				my $field = $fields[$i];
				if ($i == 0) { next; }
				# Finds undefined fields
				if (!$skip_errors && !defined($values[$i])) {
					print "\nERROR: Samples must have correct number of fields in CVS datafile (Sample '$samplename').\n\n";
					exit;
				# Skips empty fields
				} elsif ($values[$i] eq '' ){
					next;
				}
				# Old compatibility with 'tag' fields
				if ($field eq 'barcode_f'){ $field = 'tag_f'; }
				if ($field eq 'barcode_r'){ $field = 'tag_r'; }
				# Annotates sequences in uppercase
				if (in_array(['tag_f', 'tag_r'], $field)){
					$values[$i] = uc($values[$i]);
				}
				$values[$i] = trim($values[$i]);
				# Annotate all fields
				if ($values[$i] ne ''){
					$sampledata->{$samplename}{$field}=$values[$i];
				# Adds reverse-complementary reverse-barcode sequence (for easier processing of data later)
					if ($field eq 'tag_r'){
						$sampledata->{$samplename}{'tag_rc'}=iupac_reverse_complementary($values[$i]);
					}
				}
			}
			# Checks if any barcode sequence exists
			if (!$skip_errors && !defined($sampledata->{$samplename}{'tag_f'}) && !defined($sampledata->{$samplename}{'tag_r'})){
				printf("\nERROR: No barcode sequences found.\n\n", $sampleseqs{$sampledata->{$samplename}{'tag_f'}.$sampledata->{$samplename}{'tag_r'}}, $samplename);
				exit;
			} elsif (!defined($sampledata->{$samplename}{'tag_f'})){
				$sampledata->{$samplename}{'tag_f'} = '';
			} elsif (!defined($sampledata->{$samplename}{'tag_r'})){
				$sampledata->{$samplename}{'tag_r'} = '';
				$sampledata->{$samplename}{'tag_rc'} = '';
			}
			# Checks duplicated barcode pairs
			if (!$skip_errors && defined($sampleseqs{$sampledata->{$samplename}{'tag_f'}.$sampledata->{$samplename}{'tag_r'}} && $sampleseqs{$sampledata->{$samplename}{'tag_f'}.$sampledata->{$samplename}{'tag_r'}} ne $samplename) ){
# 				printf("\nWARNING: Samples '%s' & '%s' have identical barcodes in CVS datafile.\n", $sampleseqs{$sampledata->{$samplename}{'tag_f'}.$sampledata->{$samplename}{'tag_r'}}, $samplename);
				printf("\nERROR: Samples must have different sequence barcode pairs in CVS datafile (Samples: '%s' & '%s').\n\n", $sampleseqs{$sampledata->{$samplename}{'tag_f'}.$sampledata->{$samplename}{'tag_r'}}, $samplename);
				exit;
			} else {
				$sampleseqs{$sampledata->{$samplename}{'tag_f'}.$sampledata->{$samplename}{'tag_r'}} = $samplename;
			}
		# Annotates marker/primer type data:
		# marker,primer_f,primer_r,gene,feature,specie,length
		# MHC1E2,GsTGCTCCTrCTGCTGGC,CCTCGCTCTGGTTGTAGT,MHC class I,exon2,Myodes glareolus,317
		# MHC1E3,ACTACAACCAGAGCGAGG,TGTGCCTTTGGGsGAwCT,MHC class I,exon3,Myodes glareolus,313
		} elsif ($cvsdataformat eq 'markers' && (!defined($data_type) || $data_type eq 'markers' || $data_type eq 'amplicons') ) {
			my $markername = convert_to_alphanum($values[0]);
			# Checks if marker names are duplicated, can be several pairs of primers for the same gene
			if (!defined($markerdata->{$markername})){
				push(@markers,$markername);
			} #else {
				# die "\nMarkers must have different names in CVS datafile."
			#}
			# Reads fields
			for (my $i=0; $i<=$#fields; $i++){
				my $field = $fields[$i];
				if ($i == 0) { next; }
				# Finds undefined fields
				if (!defined($values[$i])) {
					print "\nERROR: Markers must have correct number of fields in CVS datafile (Marker '$markername').\n\n";
					exit;
				# Skips empty fields
				} elsif ($values[$i] eq '' ){
					next;
				}
				$values[$i] = trim($values[$i]);
				if ($values[$i] ne ''){
					if (in_array(['primer_f', 'primer_r'], $field)){
						# Annotates multiple primers and disambiguate them:
						#primer_set = disambiguate(values[i].strip())
						#for primer in primer_set:
						# Annotates all primer sequences in uppercase
						my $primer = uc($values[$i]);
						# Only annotates unique primers
						if (!in_array($markerdata->{$markername}{$field}, $primer)) {
							push(@{$markerdata->{$markername}{$field}},$primer);
							# Adds reverse-complementary reverse-primer sequence (for easier processing of data later)
							if ($field eq 'primer_r'){
								push(@{$markerdata->{$markername}{'primer_rc'}}, iupac_reverse_complementary($primer));
							}
						}
					# Annotates primer properties
					} elsif ($field =~ /^primer_/){
						push(@{$markerdata->{$markername}{$field}},$values[$i]);
					# Annotate multiple lengths
					} elsif ($field eq 'length') {
						my @lengths = split(/\s+/,trim($values[$i]));
						foreach my $length (@lengths) {
							if ($length =~ /(\d+)-(\d+)/){
								foreach my $len ($1..$2){
									if (!in_array($markerdata->{$markername}{'length'},$len)){
										push(@{$markerdata->{$markername}{'length'}},$len);
									}
								}
							} elsif (!in_array($markerdata->{$markername}{'length'},$length)) {
								push(@{$markerdata->{$markername}{'length'}}, $length);
							}
						}
					# Annotate other fields
					} else {
						$markerdata->{$markername}{$field} = $values[$i];
					}
				}
			}
		} elsif ($cvsdataformat eq 'params' && (!defined($data_type) || $data_type eq 'params') ) {
			if (!$only_params){
				my $paramname = trim($values[0]);
				my $markername = convert_to_alphanum($values[1]);
				if ($markername eq '') { $markername = 'all'; }
				my @filter_values = @values[2..$#values];
				if ($paramname ne 'allowed_markers') {
					# Checks duplicated filters
					if (defined($paramsdata->{$paramname}{$markername})) {
						foreach my $value (@filter_values){
							if (!in_array($paramsdata->{$paramname}{$markername}, $value)) {
								push(@{$paramsdata->{$paramname}{$markername}}, $value);
							} else {
								print "\nERROR: Parameter '$paramname/$markername' has duplicated value '$value'.\n\n";
								exit;
							}
						}
					} else {
						push(@{$paramsdata->{$paramname}{$markername}},@filter_values);
					}
				} else {
					@filter_values = @values[1..$#values];
					push(@{$paramsdata->{'allowed_markers'}}, @filter_values);
				}
			} else {
				my $paramname = trim($values[0]);
				my @filter_values = @values[1..$#values];
				# Checks duplicated params
				if (defined($paramsdata->{$paramname})) {
					print "\nERROR: Parameter '$paramname' has duplicated values '".join(',',@{$paramsdata->{$paramname}})."'.\n\n";
					exit;
				} else {
					push(@{$paramsdata->{$paramname}},@filter_values);
				}
			}
		} elsif ($cvsdataformat eq 'alleles' && (!defined($data_type) || $data_type eq 'alleles') ) {
			my $allelename = trim($values[0]);
			# Checks if allele names are duplicated
			if (defined($alleledata->{$allelename})) {
				print "\nERROR: Alleles must have different names in CVS datafile (Allele '$allelename').\n\n";
				exit;
			}
			# Reads fields
			for (my $i=0; $i<=$#fields; $i++){
				my $field = $fields[$i];
				if ($i == 0) { next; }
				$values[$i] = trim($values[$i]);
				# Annotate sequence
				if ($field =~ /seq/) {
					$alleledata->{$allelename}{'sequence'} = uc($values[$i]);
				# Annotate other fields
				} else {
					$alleledata->{$allelename}{$field} = $values[$i];
				}
			}
		}
	}
	close(INFILE);
	#pprint(sampledata)
	#pprint(markerdata)

	# Mix samples and primers into amplicons
	if ($data_type eq 'amplicons'){
		foreach my $sample (@samples) {
			my $sampledata_ = $sampledata->{$sample};
			foreach my $marker (@markers) {
				my $markerdata_ = $markerdata->{$marker};
				$amplicondata->{$sample."-".$marker} = { %$sampledata_, %$markerdata_ };
				$amplicondata->{$sample."-".$marker}{'sample'} = $sample;
				$amplicondata->{$sample."-".$marker}{'amplicon'} = $marker;
				push(@amplicons, $sample."-".$marker);
			}
		}
	}
	if (!defined($data_type)) {
		return ($markerdata,\@markers,$sampledata,\@samples,$paramsdata,$alleledata);
	} elsif (defined($data_type) && $data_type eq 'amplicons') {
		return ($amplicondata, \@amplicons);
	} elsif (defined($data_type) && $data_type eq 'samples') {
		return ($sampledata, \@samples);
	} elsif (defined($data_type) && $data_type eq 'markers') {
		return ($markerdata, \@markers);
	} elsif (defined($data_type) && $data_type eq 'params') {
		return $paramsdata;
	} elsif (defined($data_type) && $data_type eq 'alleles') {
		return $alleledata;
	} else {
		return undef;
	}

}

#################################################################################

# Prints amplicon data into a .csv file
sub print_amplicon_data {

	# 2 PRIMERS => MARKER
	# 1/2 TAGS => SAMPLE
	# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

	my ($data, $data_type, $clustering_params, $filtering_params) = @_;
	my $output;

	if ($data_type eq 'samples' || $data_type eq 'amplicons') {
		my @fields = ('tag_f','tag_r');
		$output .= '>sample,'.join(',',@fields)."\n";
		foreach my $samplename (nsort(keys %$data)) {
			my @values = map $data->{$samplename}{$_}, @fields;
			$output .= $samplename.','.join(',',@values)."\n";
		}
	} elsif ($data_type eq 'markers' || $data_type eq 'amplicons') {
		my @fields = ('length','primer_f','primer_r','gene','feature','specie');
		$output .= '>marker,'.join(',',@fields)."\n";
		foreach my $markername (nsort(keys %$data)) {
			my @values;
			foreach my $field (@fields) {
				if (defined($data->{$markername}{$field}) && ref($data->{$markername}{$field}) eq 'ARRAY') {
					push(@values, join(' ',@{$data->{$markername}{$field}}));
				} elsif (defined($data->{$markername}{$field})) {
					push(@values, $data->{$markername}{$field});
				} else {
					push(@values, '');
				}
			}
			$output .= $markername.','.join(',',@values)."\n";
		}
	} elsif ($data_type eq 'params') {
		my $output_ = '';
		my @fields;
		foreach my $paramname ( (@$clustering_params, @$filtering_params) ) {
			if (!defined($data->{$paramname})) {
				$output_ .= $paramname."\n";
			} elsif (ref($data->{$paramname}) eq 'HASH') {
# 			} elsif ($paramname ne 'allowed_markers') {
				foreach my $markername (nsort(keys %{$data->{$paramname}})) {
					$output_ .= $paramname.','.$markername.','.join(',',@{$data->{$paramname}{$markername}})."\n";
				}
				@fields = ('marker','value');
			} elsif (ref($data->{$paramname}) eq 'ARRAY') {
				$output_ .= $paramname.','.join(',',@{$data->{$paramname}})."\n";
				if (!@fields){
					@fields = ('value');
				}
			} else {
				$output_ .= $paramname.','.$data->{$paramname}."\n";
				if (!@fields){
					@fields = ('value');
				}
			}
		}
		$output .= ">param,".join(',',@fields)."\n".$output_;
	} elsif ($data_type eq 'alleles') {
		$output .= ">allele,sequence\n";
		foreach my $allelename (nsort(keys %$data)) {
			$output .= $allelename.','.$data->{$allelename}{'sequence'}."\n";
		}
	}

	return $output;


}

#################################################################################

# Creates a FASTA file with sequences of primers and barcodes from CSV file with amplicon data
sub extract_primer_tag_seqs {

	my ($sampledata, $samples) = @_;

	# Extracts unique primer/barcode combinations and create unique names for them
	my (@seqs, @headers, $unique_seqs, @ordered_unique_seqs);
	foreach my $samplename (@{$samples}){
		my $sampledata_ = $sampledata->{$samplename};
		my (@f_seqs,@r_seqs);
		my $count_seqs = 0;
		if (!defined($sampledata_->{'primer_f'})){
			$sampledata_->{'primer_f'} = [''];
		}
		foreach my $primer_f (@{$sampledata_->{'primer_f'}}){
			# if (!$primer_f) { next; }
			my $tag_primer_f;
			if (defined($sampledata_->{'tag_f'})){
				$tag_primer_f = uc($sampledata_->{'tag_f'}.$primer_f);
			} else {
				$tag_primer_f = uc($primer_f);
			}
			foreach my $seq (unambiguous_dna_sequences($tag_primer_f)){
				$count_seqs++;
				push(@seqs,$seq);
				push(@headers,sprintf('%s_F%03d', $samplename, $count_seqs));
				push(@f_seqs, $seq);
			}
		}
		$count_seqs = 0;
		if (!defined($sampledata_->{'primer_r'})){
			$sampledata_->{'primer_r'} = [''];
		}
		foreach my $primer_r (@{$sampledata_->{'primer_r'}}){
			# if (!$primer_r) { next; }
			my $tag_primer_r;
			if (defined($sampledata_->{'tag_r'})){
				$tag_primer_r = uc(iupac_reverse_complementary($sampledata_->{'tag_r'}.$primer_r));
			} else {
				$tag_primer_r = uc(iupac_reverse_complementary($primer_r));
			}
			foreach my $seq (unambiguous_dna_sequences($tag_primer_r)){
				$count_seqs++;
				push(@seqs,$seq);
				push(@headers,sprintf('%s_R%03d', $samplename, $count_seqs));
				push(@r_seqs, $seq);
			}
		}
		# Checks if sequences are repeated
		my $f_count = 0;
		foreach my $f_seq (@f_seqs) {
			$f_count++;
			my $r_count = 0;
			foreach my $r_seq (@r_seqs) {
				$r_count++;
# 				if (!defined($unique_seqs->{$f_seq.$r_seq})){
# 					unique_seqs[f_seq+r_seq] = []
				push(@{$unique_seqs->{$f_seq.$r_seq}}, sprintf('%s F%03d-R%03d', $samplename, $f_count, $r_count));
				if (!in_array(\@ordered_unique_seqs, $f_seq.$r_seq)){
					push(@ordered_unique_seqs, $f_seq.$r_seq);
				}
			}
		}
	}
	foreach my $seq (@ordered_unique_seqs){
		if ($#{$unique_seqs->{$seq}} > 0) {
			printf("\t'%s' have common forward and reverse sequences.\n", join("' and '", @{$unique_seqs->{$seq}}));
		}
	}

	return \@seqs, \@headers;

}


#################################################################################

# Creates a FASTA file with IUPAC sequences of primers and barcodes from CSV file with amplicon data
sub extract_primer_tag_iupac_seqs {

	my ($sampledata, $samples) = @_;

	# Extracts unique primer/barcode combinations and create unique names for them
	my (@seqs, @headers, $unique_seqs, @ordered_unique_seqs);
	foreach my $samplename (@{$samples}){
		my $sampledata_ = $sampledata->{$samplename};
		my (@f_seqs,@r_seqs);
		my $count_seqs = 0;
		for my $primer_f (@{$sampledata_->{'primer_f'}}){
			if (!$primer_f) { next; }
			my $tag_primer_f;
			if (defined($sampledata_->{'tag_f'})){
				$tag_primer_f = uc($sampledata_->{'tag_f'}.$primer_f);
			} else {
				$tag_primer_f = uc($primer_f);
			}
			$count_seqs++;
			push(@seqs,$tag_primer_f);
			push(@headers,sprintf('%s_F%03d', $samplename, $count_seqs));
			push(@f_seqs, $tag_primer_f);
		}
		$count_seqs = 0;
		foreach my $primer_r (@{$sampledata_->{'primer_r'}}){
			if (!$primer_r) { next; }
			my $tag_primer_r;
			if (defined($sampledata_->{'tag_r'})){
				$tag_primer_r = uc(iupac_reverse_complementary($sampledata_->{'tag_r'}.$primer_r));
			} else {
				$tag_primer_r = uc(iupac_reverse_complementary($primer_r));
			}
			$count_seqs++;
			push(@seqs,$tag_primer_r);
			push(@headers,sprintf('%s_R%03d', $samplename, $count_seqs));
			push(@r_seqs, $tag_primer_r);
		}
		# Checks if sequences are repeated
		my $f_count = 0;
		foreach my $f_seq (@f_seqs) {
			$f_count++;
			my $r_count = 0;
			foreach my $r_seq (@r_seqs) {
				$r_count++;
# 				if (!defined($unique_seqs->{$f_seq.$r_seq})){
# 					unique_seqs[f_seq+r_seq] = []
				push(@{$unique_seqs->{$f_seq.$r_seq}}, sprintf('%s F%03d-R%03d', $samplename, $f_count, $r_count));
				if (!in_array(\@ordered_unique_seqs, $f_seq.$r_seq)){
					push(@ordered_unique_seqs, $f_seq.$r_seq);
				}
			}
		}
	}
	foreach my $seq (@ordered_unique_seqs){
		if ($#{$unique_seqs->{$seq}} > 0) {
			printf("\t'%s' have common forward and reverse sequences.\n", join("' and '", @{$unique_seqs->{$seq}}));
		}
	}

	return \@seqs, \@headers;

}

#################################################################################

# Aligns reads against primer+barcode sequences
sub align_amplicons {

	# 2 PRIMERS => MARKER
	# 1/2 TAGS => SAMPLE
	# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

	my ($read_headers,$read_seqs,$primer_tag_headers,$primer_tag_seqs,$align_type,$revcomp,$threads) = @_;

	# Default alignment
	if (!defined($align_type)){
		$align_type = 'match';
	}

	# Perform alignment of all the reads against primer/barcode sequences
	# GASSST program allows to align the full primer sequence with the read (local+global alignment)
	my ($align_primer_data,$align_data);
	# Normal blastn finds less than half of the alignments because some high e-value problems
	my $align_options;
	if ($align_type =~ /match/){
		$align_options = "dna match";
		if (defined($revcomp)){
			$align_options = " dna match revcomp";
		}
	} elsif ($align_type =~ /gassst/){
		$align_options = "dna gassst -w 6 -p 90 -l 0 -s 5 -r 0";
		if (defined($revcomp)){
			$align_options = "dna gassst -w 6 -p 90 -l 0 -s 5 -r 1";
		}
	} elsif ($align_type =~ /short/){
		$align_options = "dna blastn-short -evalue 0.001 -strand plus";
		if (defined($revcomp)){
			$align_options = "dna blastn-short -evalue 0.001 -strand both";
		}
	} elsif ($align_type =~ /blast/){
		$align_options = "dna blastn -evalue 0.001 -strand plus";
		if (defined($revcomp)){
			$align_options = "dna blastn -evalue 0.001 -strand both";
		}
	}
# 	print "\nAligning primer/barcode sequences.\n";
	if (defined($threads)) {
		$align_data = align_seqs_with_threads($read_headers,$read_seqs,$primer_tag_headers,$primer_tag_seqs,0,$align_options,$threads);
	} else {
		$align_data = align_seqs($read_headers,$read_seqs,$primer_tag_headers,$primer_tag_seqs,0,$align_options);
	}

	return $align_data;

}

#################################################################################

# Finds amplicons in alignment results of reads and primers+barcodes
sub match_amplicons {

	# 2 PRIMERS => MARKER
	# 1/2 TAGS => SAMPLE
	# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

	my ($read_headers,$read_seqs,$markerdata,$sampledata,$align_data,$max_depth) = @_;

	# Initialize variables and counters to store analysis data
	my $amplicon_depths;
	my $amplicon_sequences;
	my $md5_to_marker_name;
	my $md5_to_sequence;
	my $md5_to_quality;

	# Stores matched reads
	# my ($read_headers_matched, $read_seqs_matched);

	# Loops reads with alignment results
	for (my $i=0; $i<=$#{$read_headers}; $i++){
		my $read_header = $read_headers->[$i];
		if (!defined($align_data->{$read_header})){
			next;
		}
		my $read_seq = $read_seqs->[$i];
# 		my $read_quality;
# 		if (defined($read_qualities) && @{$read_qualities}){
# 			$read_quality = $read_qualities->[$i];
# 		}
# 		my $read_length = length($read_seq);

		# Find between results matches with common forward and reverse sequences (primer+barcode)
		my ($forward_seqs, $reverse_seqs);
		my $amplicon_found;
		foreach my $result (@{$align_data->{$read_header}}) {
			my $primer_tag_header = $result->{'NAME'};
			$primer_tag_header =~ /(.+)_([F|R])\d+/;
			# If matched sequence is forward primer/barcode
			if ($2 eq 'F'){
				# Annotate forward sequences matched
				if (!defined($forward_seqs->{$1})){
					$forward_seqs->{$1} = $result;
				}
				# Stop checking results if the same amplicon has been detected in reverse sequences
				if (defined($reverse_seqs) && defined($reverse_seqs->{$1})){
					$amplicon_found = $1;
					last;
				}
			# If matched sequence is reverse primer/barcode
			} elsif ($2 eq 'R'){
				# Annotate reverse sequences matched
				if (!defined($reverse_seqs->{$1})){
					$reverse_seqs->{$1} = $result;
				}
				# Stop checking results if the same primer_barcode has been detected in forward sequences
				if (defined($forward_seqs) && defined($forward_seqs->{$1})){
					$amplicon_found = $1;
					last;
				}
			}
		}

		# If the two primers+barcodes forward and reverse of the same amplicon+sample are found
		if (defined($amplicon_found)){
			# Extracts names and sequences of the primers+barcodes, it should be the same name for forward and reverse
			$forward_seqs->{$amplicon_found}{'NAME'} =~ /(.+)-(.+)_([F|R])\d+/;
			my $sample_name = $1;
			my $marker_name = $2;
			if (!(defined($max_depth) && $max_depth>0 && $amplicon_depths->{$marker_name}{$sample_name}>$max_depth)) {
				# Annotates the amplicon sequence
				my ($forward_read_aligned_cols, $forward_amplicon_aligned_cols) = split("\n",$forward_seqs->{$amplicon_found}{'COLS'});
				my @forward_read_aligned_cols = split(",",$forward_read_aligned_cols);
				my @forward_amplicon_aligned_cols = split(",",$forward_amplicon_aligned_cols);
				my ($first_amplicon_pos,$last_amplicon_pos,$amplicon_seq);
				my $is_direct_seq = 0;
				if ($forward_amplicon_aligned_cols[-1]>$forward_amplicon_aligned_cols[0]){
					# One nucleotide after forward primer
					$first_amplicon_pos = $forward_read_aligned_cols[-1]+1;
					$is_direct_seq++;
				} else {
					# One nucleotide before forward primer aligned in reverse complementary position
					$last_amplicon_pos = $forward_read_aligned_cols[0]-1;
					$is_direct_seq--;
				}
				my ($reverse_read_aligned_cols, $reverse_amplicon_aligned_cols) = split("\n",$reverse_seqs->{$amplicon_found}{'COLS'});
				my @reverse_read_aligned_cols = split(",",$reverse_read_aligned_cols);
				my @reverse_amplicon_aligned_cols = split(",",$reverse_amplicon_aligned_cols);
				if ($reverse_amplicon_aligned_cols[-1]>$reverse_amplicon_aligned_cols[0]){
					# One nucleotide before reverse primer
					$last_amplicon_pos = $reverse_read_aligned_cols[0]-1;
					$is_direct_seq++;
				} else {
					# One nucleotide after forward primer
					$first_amplicon_pos = $reverse_read_aligned_cols[-1]+1;
					$is_direct_seq--;
				}
	# 			my $tag_fwd_length = length($sampledata->{$sample_name}{'tag_f'});
	# 			my $tag_rev_length = length($sampledata->{$sample_name}{'tag_r'});
				if ($is_direct_seq == 2){
					$amplicon_seq = substr($read_seq,$first_amplicon_pos-1, $last_amplicon_pos-$first_amplicon_pos+1);
				} elsif ($is_direct_seq == -2){
					$amplicon_seq = iupac_reverse_complementary(substr($read_seq,$first_amplicon_pos-1, $last_amplicon_pos-$first_amplicon_pos+1));
				} else {
					next;
				}
				my $md5 = generate_md5($amplicon_seq); # , 8, ['base64']
	# 			if (defined($read_quality)) {
	# 				my $amplicon_qual = substr($read_quality,$first_amplicon_pos-1, $last_amplicon_pos-$first_amplicon_pos+1);
	# 				push(@{$md5_to_quality->{$md5}},$amplicon_qual);
	# 			}
	# 			# Check if amplicon length is correct
	# 			if (defined($markerdata->{$marker_name}{'length'}) && defined($max_amplicon_length_error)){
	# 				my $amplicon_length = ($reverse_read_aligned_cols[-1]-$forward_read_aligned_cols[0]+1)-($tag_fwd_length+$tag_fwd_length);
	# 				if ($amplicon_length > $markerdata->{$marker_name}{'length'}+$max_amplicon_length_error || $amplicon_length < $markerdata->{$marker_name}{'length'}-$max_amplicon_length_error ){
	# 					next;
	# 				}
	# 			}
				if (!defined($md5_to_sequence->{$md5})){ 
					$md5_to_sequence->{$md5}=$amplicon_seq;
				} elsif ($md5_to_sequence->{$md5} ne $amplicon_seq){
					print "\nERROR: Sequence ID '$md5' is not unique.\n";
					print "SEQ1: ".$md5_to_sequence->{$md5}."\n";
					print "SEQ2: ".$amplicon_seq."\n\n";
					exit;
				}
				if (!defined($md5_to_marker_name->{$md5})){ 
					$md5_to_marker_name->{$md5}=$marker_name;
				} elsif ($md5_to_marker_name->{$md5} ne $marker_name){
					print "\nERROR: Sequence '$md5' has more than one amplicon assigned.\n\n";
					exit;
				}
				# print "\n$read_header $sample_name $marker_name $amplicon_seq\n\n";exit;
				$amplicon_sequences->{$marker_name}{$sample_name}{$md5}++;
				# Counts how many times a sample/amplicon assignment is found
				$amplicon_depths->{$marker_name}{$sample_name}++;
				
	# 			# Store reads matched for printing
	# 			push(@$read_headers_matched, $read_header);
	# 			push(@$read_seqs_matched, $read_seq);
			}
		}
	}

# 	foreach my $md5 (keys %{$md5_to_quality}){
# 		$md5_to_quality->{$md5} = mean_phred_quality($md5_to_quality->{$md5});
# 		print '';
# 	}
	
	# Create a file with reads matched
	# create_fasta_file($read_seqs_matched,$read_headers_matched,"matched_reads.fa.gz",1);

	return ($md5_to_sequence,$amplicon_sequences,$amplicon_depths);

}


#################################################################################

# Parse a file with one sequence per line and primers+barcodes with REGEX (AWK or PERL) to find matching amplicons
# Is equivalent to 'align_amplicons+match_amplicons' with perfect matching, but very fast
sub match_amplicons_regex {

	# 2 PRIMERS => MARKER
	# 1/2 TAGS => SAMPLE
	# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

	my ($seqs_file,$markerdata,$sampledata,$primers,$barcodes,$max_depth,$options) = @_;

	my ($md5_to_marker_name,$md5_to_sequence,$amplicon_sequences,$amplicon_depths);
	
	my ($verbose,$revcomp,$partial)=(1,1,0);
	if (in_array($options, 'quiet')){
		$verbose = 0;
	}
	if (in_array($options, 'direct')){
		$revcomp = 0;
	}
	if (in_array($options, 'partial')){
		$partial = 1;
	}

	# Assigns reads to primer/barcode sequences
	my $patterns_file = "/tmp/".random_file_name();

# 	# Stores matched reads
# 	my @read_seqs = @{read_from_file($seqs_file)};
# 	my ($read_headers_matched, $read_seqs_matched);

	my (%non_unique_md5s,%multimarker_md5s);
	foreach my $marker_name (@$primers){

		if (!defined($markerdata->{$marker_name}{'primer_f'}) && !defined($markerdata->{$marker_name}{'primer_r'}) && !defined($markerdata->{$marker_name}{'primer_rc'})){
			printf("\nERROR: Marker '%s' has no primers defined.\n\n", $marker_name);
			exit;
		}

		my @primers_f;
		if (defined($markerdata->{$marker_name}{'primer_f'})){
			@primers_f = @{$markerdata->{$marker_name}{'primer_f'}};
		} else {
			@primers_f = ('');
		}
		
		my @primers_rc;
		if (defined($markerdata->{$marker_name}{'primer_rc'})){
			@primers_rc = @{$markerdata->{$marker_name}{'primer_rc'}};
		} elsif (defined($markerdata->{$marker_name}{'primer_r'})){
			foreach my $primer_r (@{$markerdata->{$marker_name}{'primer_r'}}){
				push(@primers_rc, uc(iupac_reverse_complementary($primer_r)));
			}
		} else {
			@primers_rc = ('');
		}

		# If no barcodes are provided, match only primers
		if (!defined($barcodes) || !@$barcodes) {
			$sampledata->{''}{'tag_f'} = '';
			$sampledata->{''}{'tag_rc'} = '';
			@$barcodes = ('');
		}

		foreach my $sample_name (@$barcodes) {

			if ($verbose && $sample_name ne ''){
				print "\t$marker_name-$sample_name de-multiplexing\n";
			} elsif ($verbose) {
				print "\t$marker_name de-multiplexing\n";
			}

			my $tag_f = $sampledata->{$sample_name}{'tag_f'};
			my $tag_rc = $sampledata->{$sample_name}{'tag_rc'};
			my @patterns;
			foreach my $primer_f (@primers_f) {
				foreach my $primer_rc (@primers_rc) {
					# Very important the parenthesis in the regex for later locate the start position and the length of the sequence between primers and barcodes
					$patterns[0] .= regex($tag_f.$primer_f."(.+)".$primer_rc.$tag_rc)."\n";
					if ($revcomp) {
						$patterns[1] .= regex(iupac_reverse_complementary($primer_rc.$tag_rc)."(.+)".iupac_reverse_complementary($tag_f.$primer_f))."\n";
					}
				}
			}
			my %pos_seqs;
			my $count_total_seqs_amplicon = 0;
			for (my $i=0; $i<=$#patterns; $i++) {
				# If defined $max_depth, limit the number of sequences per amplicon
				if (defined($max_depth) && $max_depth>0 && $count_total_seqs_amplicon>=$max_depth) {
					last;
				}
				write_to_file($patterns_file, $patterns[$i]);
				# Match in $seqs_file the primer+barcode patterns from $patterns_file and returns the matched line, the start position and the length of the first parenthesis match
				# print "awk 'FNR==NR{a[\$0];next} {for (i in a) {if (p=match(\$0, i, arr)) print FNR, arr[1,\"start\"], arr[1,\"length\"], substr(\$0, arr[1,\"start\"], arr[1,\"length\"])} }' $patterns_file $seqs_file\n"; exit;
				# open(AWK_MATCHES, "awk 'FNR==NR{a[\$0];next} {for (i in a) {if (p=match(\$0, i, arr)) print FNR, arr[1,\"start\"], arr[1,\"length\"], substr(\$0, arr[1,\"start\"], arr[1,\"length\"])} }' $patterns_file $seqs_file|");
				# print "perl -e 'open(PRIMERFILE,\$ARGV[0]);while(<PRIMERFILE>){chomp;if(\$_){push(\@primers,\$_);}}close(PRIMERFILE);open(SEQSFILE,\$ARGV[1]);\$line=0;while(<SEQSFILE>){\$line++;chomp;if(\$_){foreach \$primer(\@primers){if(/\$primer/i){printf(\"\%d\\t\%d\\t\%d\\t\%s\\n\",\$line,\$-[1]+1,length(\$1),\$1);}}}}close(SEQSFILE);' $patterns_file $seqs_file\n"; exit;
				open(REGEX_MATCHES, "perl -e 'open(PRIMERFILE,\$ARGV[0]);while(<PRIMERFILE>){chomp;if(\$_){push(\@primers,\$_);}}close(PRIMERFILE);open(SEQSFILE,\$ARGV[1]);\$line=0;while(<SEQSFILE>){\$line++;chomp;if(\$_){foreach \$primer(\@primers){if(/\$primer/i){printf(\"\%d\\t\%d\\t\%d\\t\%s\\n\",\$line,\$-[1]+1,length(\$1),\$1);}}}}close(SEQSFILE);' $patterns_file $seqs_file|");
				while (<REGEX_MATCHES>) {
					chomp;
					my ($pos_seq, $first_amplicon_pos, $amplicon_length, $amplicon_seq) = split("\t");
					# If the read has already been matched by another primer in the same amplicon (primers could be duplicated in 2 degenerated ones)
					if (!defined($amplicon_seq) || defined($pos_seqs{$pos_seq})){
						next;
					}
					$pos_seqs{$pos_seq} = 1;

					if ($i == 1) {
						$amplicon_seq = iupac_reverse_complementary($amplicon_seq);
					}
					my $md5 = generate_md5($amplicon_seq);
# if ($md5 eq 'd82eba7823d1ae6dad22173f13c3a506' || $md5 eq '31214183d13676b38cd11cff4b9f4a15'){
# print '';
# }

					if (!defined($md5_to_sequence->{$md5})){ 
						$md5_to_sequence->{$md5}=$amplicon_seq;
					} elsif ($verbose && !defined($non_unique_md5s{$md5}) && $md5_to_sequence->{$md5} ne $amplicon_seq){
						printf("\tWARNING: Sequence '%s' is not unique:\n\t\tSEQ1: %s\n\t\tSEQ2: %s\n", $md5, $md5_to_sequence->{$md5}, $amplicon_seq);
						$non_unique_md5s{$md5} = 1;
# 						exit;
					}
					if (!defined($md5_to_marker_name->{$md5})){ 
						$md5_to_marker_name->{$md5}=$marker_name;
					} elsif ($verbose && !defined($multimarker_md5s{$md5}) && $md5_to_marker_name->{$md5} ne $marker_name){
						printf("\tWARNING: Sequence '%s' has more than one marker assigned: '%s' & '%s'.\n", $md5, $md5_to_marker_name->{$md5}, $marker_name);
						$multimarker_md5s{$md5} = 1;
# 						exit;
					}
			# 		print "\n$read_header $sample_name $marker_name $amplicon_seq\n\n";exit;
					$amplicon_sequences->{$marker_name}{$sample_name}{$md5}++;
					# Counts how many times a sample/amplicon assignment is found
					$amplicon_depths->{$marker_name}{$sample_name}++;
					# Count total number of sequences in an amplicon
					$count_total_seqs_amplicon++;
			
					# If defined $max_depth, limit the number of sequences per amplicon
					if (defined($max_depth) && $max_depth>0 && $count_total_seqs_amplicon>=$max_depth) {
						last;
					}

# 					if (defined($read_qualities) && @{$read_qualities}){
# 						my $amplicon_qual;
# 						if ($i == 0){
# 							$amplicon_qual = substr($read_qualities->[$pos_seq-1],$first_amplicon_pos-1, $amplicon_length);
# 						} else {
# 							$amplicon_qual = reverse_sequence(substr($read_qualities->[$pos_seq-1],$first_amplicon_pos-1, $amplicon_length));
# 						}
# 						push(@{$md5_to_quality->{$md5}},$amplicon_qual);
# 					}

# 					# Stores reads matched for printing
# 					push(@$read_headers_matched, $pos_seq);
# 					push(@$read_seqs_matched, $read_seqs[$pos_seq-1]);

				}
				close(REGEX_MATCHES);
			}

			if ($verbose && $sample_name ne ''){
				printf("\t%s-%s de-multiplexed (%d reads, %d variants)\n", $marker_name, $sample_name, $count_total_seqs_amplicon, scalar keys %{$amplicon_sequences->{$marker_name}{$sample_name}});
			} elsif ($verbose) {
				printf("\t%s de-multiplexed (%d reads, %d variants)\n", $marker_name, $count_total_seqs_amplicon, scalar keys %{$amplicon_sequences->{$marker_name}{$sample_name}});
			}

		}
	}

	`rm $patterns_file`;
	
# 	# Creates a file with reads matched
# 	print "\nMatched reads printed into 'matched_reads.fa.gz'.\n\n";
# 	create_fasta_file($read_seqs_matched,$read_headers_matched,"matched_reads.fa.gz",'gzip');
# 	exit;
	
	return ($md5_to_sequence,$amplicon_sequences,$amplicon_depths);

}

#################################################################################

# Parse reads and primers+barcodes with REGEX (AWK or PERL) to find matching amplicons
# Is equivalent to 'align_amplicons+match_amplicons' with perfect matching, but very fast
sub match_amplicons_regex_with_threads {

	my ($seqs_file,$markerdata,$sampledata,$primers,$barcodes,$max_depth,$options,$threads_limit) = @_;

	if (!defined($threads_limit)){
		$threads_limit = 4;
	}

	my ($md5_to_sequence,$amplicon_sequences,$amplicon_depths) = ({},{},{},{});

	my @amplicons;
	foreach my $marker_name (@{$primers}){
		# If no barcodes are provided, match only primers
		if (!defined($barcodes) || !@$barcodes) {
			@$barcodes = ('');
		}
		foreach my $sample_name (@$barcodes) {
			push(@amplicons,"$marker_name-$sample_name");
		}
	}

	my @threads;
	for (my $count_amplicon=0; $count_amplicon<=$#amplicons; $count_amplicon++){

		$amplicons[$count_amplicon] =~ /(.*)-(.*)/;
		my ($marker_name,$sample_name) = ($1, $2);

		push(@threads, threads->create(\&match_amplicons_regex,$seqs_file,$markerdata,$sampledata,[$marker_name],[$sample_name],$max_depth,$options));
# 		print "\n";

# 		# For debugging:
# 		push(@threads, [match_amplicons_regex($seqs_file,$markerdata,$sampledata,[$marker_name],[$sample_name],$max_depth,$options)]);

		# If maximum number of threads is reached or last sbjct of a query is processed
		if (scalar @threads >= $threads_limit  || $count_amplicon == $#amplicons){
			my $check_threads = 1;
			while ($check_threads){
				for (my $i=0; $i<=$#threads; $i++){
					unless ($threads[$i]->is_running()){
						my ($md5_to_sequence_,$amplicon_sequences_,$amplicon_depths_) = $threads[$i]->join;
						if (defined($amplicon_sequences_)){
							my $marker_name = (keys %$amplicon_sequences_)[0];
							my $sample_name = (keys %{$amplicon_sequences_->{$marker_name}})[0];
	# 						print "\t$marker_name-$sample_name finished\n";
							if (defined($amplicon_sequences_->{$marker_name}{$sample_name})){
								$amplicon_sequences->{$marker_name}{$sample_name} = $amplicon_sequences_->{$marker_name}{$sample_name};
								$amplicon_depths->{$marker_name}{$sample_name} = $amplicon_depths_->{$marker_name}{$sample_name};
							}
							if (defined($md5_to_sequence_)){
								$md5_to_sequence = { %$md5_to_sequence, %$md5_to_sequence_ };
							}
# 							if (defined($md5_to_quality_)){
# 								$md5_to_quality = { %$md5_to_quality, %$md5_to_quality_ };
# 							}
						}
						undef($threads[$i]);
						splice(@threads,$i,1);
						$i--;
						unless ($count_amplicon == $#amplicons && @threads){
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
# 				print '';
# 				my ($md5_to_sequence_,$amplicon_sequences_,$amplicon_depths_) = @{$threads[$i]};
# 				if (defined($amplicon_sequences_)){
# 					my $marker_name = (keys %$amplicon_sequences_)[0];
# 					my $sample_name = (keys %{$amplicon_sequences_->{$marker_name}})[0];
# 	# 				print "\t$marker_name-$sample_name finished\n";
# 					$amplicon_sequences->{$marker_name}{$sample_name} = $amplicon_sequences_->{$marker_name}{$sample_name};
# 					$amplicon_depths->{$marker_name}{$sample_name} = $amplicon_depths_->{$marker_name}{$sample_name};
# 					$md5_to_sequence = { %$md5_to_sequence, %$md5_to_sequence_ };
# # 					if (defined($md5_to_quality_)){
# # 						$md5_to_quality = { %$md5_to_quality, %$md5_to_quality_ };
# # 					}
# 				}
# 				delete $threads[$i];
# 			}
		}
	}
	
	return ($md5_to_sequence,$amplicon_sequences,$amplicon_depths);

}


#################################################################################

# Parse reads and primers+barcodes with REGEX (AWK or PERL) to find matching amplicons
# And retrieves full sequences and headers of matching reads
sub find_amplicon_reads {

	# 2 PRIMERS => MARKER
	# 1/2 TAGS => SAMPLE
	# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

	my ($seqs_file,$markerdata,$sampledata,$primers,$barcodes,$max_depth,$options) = @_;

	my ($verbose,$revcomp,$partial)=(1,1,0);
	if (in_array($options, 'quiet')){
		$verbose = 0;
	}
	if (in_array($options, 'direct')){
		$revcomp = 0;
	}
	if (in_array($options, 'partial')){
		$partial = 1;
	}

	# Assigns reads to primer/barcode sequences
	my $patterns_file = "/tmp/".random_file_name();

	# Saves the positions of matched sequences
	my %pos_seqs_matched;

	foreach my $marker_name (@$primers){

		if (!defined($markerdata->{$marker_name}{'primer_f'}) && !defined($markerdata->{$marker_name}{'primer_r'}) && !defined($markerdata->{$marker_name}{'primer_rc'})){
			printf("\nERROR: Marker '%s' has no primers defined.\n\n", $marker_name);
			exit;
		}

		my @primers_f;
		if (defined($markerdata->{$marker_name}{'primer_f'})){
			@primers_f = @{$markerdata->{$marker_name}{'primer_f'}};
		} else {
			@primers_f = ('');
		}
		
		my @primers_rc;
		if (defined($markerdata->{$marker_name}{'primer_rc'})){
			@primers_rc = @{$markerdata->{$marker_name}{'primer_rc'}};
		} elsif (defined($markerdata->{$marker_name}{'primer_r'})){
			foreach my $primer_r (@{$markerdata->{$marker_name}{'primer_r'}}){
				push(@primers_rc, uc(iupac_reverse_complementary($primer_r)));
			}
		} else {
			@primers_rc = ('')
		}

		# If no barcodes are provided, matches only primers
		if (!defined($barcodes) || !@$barcodes) {
			$sampledata->{''}{'tag_f'} = '';
			$sampledata->{''}{'tag_rc'} = '';
			@$barcodes = ('');
		}


		foreach my $sample_name (@$barcodes) {

			if ($verbose && $sample_name ne ''){
				print "\t$marker_name-$sample_name processing\n";
			} elsif ($verbose) {
				print "\t$marker_name processing\n";
			}

			my $tag_f = $sampledata->{$sample_name}{'tag_f'};
			my $tag_rc = $sampledata->{$sample_name}{'tag_rc'};
			my @patterns;
			foreach my $primer_f (@primers_f) {
				foreach my $primer_rc (@primers_rc) {
					# Very important the parenthesis in the regex for later locate the start position and the length of the sequence between primers and barcodes
					$patterns[0] .= regex($tag_f.$primer_f."(.+)".$primer_rc.$tag_rc)."\n";
					if ($revcomp) {
						$patterns[1] .= regex(iupac_reverse_complementary($primer_rc.$tag_rc)."(.+)".iupac_reverse_complementary($tag_f.$primer_f))."\n";
					}
				}
			}
			my $count_total_seqs_amplicon = 0;
			for (my $i=0; $i<=$#patterns; $i++) {
				# If defined $max_depth, limit the number of sequences per amplicon
				if (defined($max_depth) && $max_depth>0 && $count_total_seqs_amplicon>=$max_depth) {
					last;
				}
				# Annotates if the primer is found in direct or complementary DNA chain
				my $dna_dir = 'd';
				if ($i == 1){ $dna_dir = 'c' }
				write_to_file($patterns_file, $patterns[$i]);
				# Match in $seqs_file the primer+barcode patterns from $patterns_file and returns the matched line, the start position and the length of the first parenthesis match
				# print "awk 'FNR==NR{a[\$0];next} {for (i in a) {if (p=match(\$0, i, arr)) print FNR, arr[1,\"start\"], arr[1,\"length\"], substr(\$0, arr[1,\"start\"], arr[1,\"length\"])} }' $patterns_file $seqs_file\n"; exit;
				# open(AWK_MATCHES, "awk 'FNR==NR{a[\$0];next} {for (i in a) {if (p=match(\$0, i, arr)) print FNR, arr[1,\"start\"], arr[1,\"length\"], substr(\$0, arr[1,\"start\"], arr[1,\"length\"])} }' $patterns_file $seqs_file|");
				# print "perl -e 'open(PRIMERFILE,\$ARGV[0]);while(<PRIMERFILE>){chomp;if(\$_){push(\@primers,\$_);}}close(PRIMERFILE);open(SEQSFILE,\$ARGV[1]);\$line=0;while(<SEQSFILE>){\$line++;chomp;if(\$_){foreach \$primer(\@primers){if(/\$primer/){printf(\"\%d\\t\%d\\t\%d\\t\%s\\n\",\$line,\$-[1]+1,length(\$1),\$1);}}}}close(SEQSFILE);' $patterns_file $seqs_file\n"; exit;
				open(REGEX_MATCHES, "perl -e 'open(PRIMERFILE,\$ARGV[0]);while(<PRIMERFILE>){chomp;if(\$_){push(\@primers,\$_);}}close(PRIMERFILE);open(SEQSFILE,\$ARGV[1]);\$line=0;while(<SEQSFILE>){\$line++;chomp;if(\$_){foreach \$primer(\@primers){if(/\$primer/){printf(\"\%d\\t\%d\\t\%d\\t\%s\\n\",\$line,\$-[1]+1,length(\$1),\$1);}}}}close(SEQSFILE);' $patterns_file $seqs_file|");
				while (<REGEX_MATCHES>) {
					my $pos_seq = (split("\t"))[0];
					# If the read has already been matched by another primer in the same amplicon (primers could be duplicated in 2 degenerated ones)
					if (defined($pos_seqs_matched{$pos_seq})){
						next;
					}
					$pos_seqs_matched{$pos_seq} = $dna_dir;
					# Count total number of sequences in an amplicon
					$count_total_seqs_amplicon++;
					# If defined $max_depth, limit the number of sequences per amplicon
					if (defined($max_depth) && $max_depth>0 && $count_total_seqs_amplicon>=$max_depth) {
						last;
					}
				}
				close(REGEX_MATCHES);
			}
			if ($verbose && $sample_name ne ''){
				printf("\t%s-%s processed, found %d sequences\n", $marker_name, $sample_name, $count_total_seqs_amplicon);
			} elsif ($verbose) {
				printf("\t%s processed, found %d sequences\n", $marker_name, $count_total_seqs_amplicon);
			}
		}
	}
	`rm  $patterns_file`;

	return \%pos_seqs_matched;

}


#################################################################################

# Parse reads and primers+barcodes with REGEX (AWK or PERL) to find matching amplicons
# And retrieves full sequences and headers of matching reads
sub find_amplicon_reads_with_threads {

	my ($seqs_file,$markerdata,$sampledata,$primers,$barcodes,$max_depth,$options,$threads_limit) = @_;

	if (!defined($threads_limit)){
		$threads_limit = 4;
	}

	# Saves the positions of matched sequences
	my $pos_seqs_matched = {};

	my @amplicons;
	foreach my $marker_name (@{$primers}){
		# If no barcodes are provided, match only primers
		if (!defined($barcodes) || !@$barcodes) {
			@$barcodes = ('');
		}
		foreach my $sample_name (@$barcodes) {
			push(@amplicons,"$marker_name-$sample_name");
		}
	}

	my @threads;
	for (my $count_amplicon=0; $count_amplicon<=$#amplicons; $count_amplicon++){

		$amplicons[$count_amplicon] =~ /(.*)-(.*)/;
		my ($marker_name,$sample_name) = ($1, $2);

		push(@threads, threads->create(\&find_amplicon_reads,$seqs_file,$markerdata,$sampledata,[$marker_name],[$sample_name],$max_depth,$options));
# 		print "\n";

# 		# For debugging:
# 		push(@threads, find_amplicon_reads($seqs_file,$markerdata,$sampledata,[$marker_name],[$sample_name],$max_depth,$options));

		# If maximum number of threads is reached or last sbjct of a query is processed
		if (scalar @threads >= $threads_limit  || $count_amplicon == $#amplicons){
			my $check_threads = 1;
			while ($check_threads){
				for (my $i=0; $i<=$#threads; $i++){
					unless ($threads[$i]->is_running()){
						my $pos_seqs_matched_ = $threads[$i]->join;
						if (defined($pos_seqs_matched_)){
							$pos_seqs_matched = { %$pos_seqs_matched, %$pos_seqs_matched_ };
						}
						undef($threads[$i]);
						splice(@threads,$i,1);
						$i--;
						unless ($count_amplicon == $#amplicons && @threads){
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
# 				print '';
# 				my $pos_seqs_matched_ = $threads[$i];
# 				if (defined($pos_seqs_matched_)){
# 					$pos_seqs_matched = { %$pos_seqs_matched, %$pos_seqs_matched_ };
# 				}
# 				delete $threads[$i];
# 			}
		}
	}
	
	return $pos_seqs_matched;


}

#################################################################################

# Reads a sequence file to extract alleles in $alleledata format
sub read_allele_file {

	my ($seqs_file,$options) = @_;

	my ($seqs,$headers,$seqs_file_format);

	if (is_fastq($seqs_file)){
		$seqs_file_format = 'fastq';
	} elsif (is_fasta($seqs_file)){
		$seqs_file_format = 'fasta';
	}
	
	if ($seqs_file_format eq 'fastq'){
		($seqs,$headers) = read_fastq_file($seqs_file)
	} elsif ($seqs_file_format eq 'fasta'){
		($seqs,$headers) = read_fasta_file($seqs_file);
	}

	return read_alleles($seqs,$headers,$options);

}

#################################################################################

# Reads a sequence file to extract alleles in $alleledata format
sub read_alleles {

	my ($allele_seqs,$allele_names,$options) = @_;

	my $alleledata;

	my $show_warnings = 1;
	if (defined($options) && in_array($options, 'no warnings')){
		$show_warnings = 0;
	}
	
	my @previous_allele_seqs;
	for (my $i=0; $i<=$#{$allele_names}; $i++) {
		if ($show_warnings && defined($alleledata->{$allele_names->[$i]}) ){
			print "\tWARNING: Allele name '".$allele_names->[$i]."' is duplicated.\n";
# 			exit;
		} elsif ($show_warnings && in_array(\@previous_allele_seqs,$allele_seqs->[$i])) {
			print "\tWARNING: Allele '".$allele_names->[$i]."' sequence is duplicated.\n";
		} else {
			my $allele_data = extract_header_data($allele_names->[$i]);
			# Removes gaps to avoid errors in further alignments
			$allele_seqs->[$i] =~ s/-//g;
			$alleledata->{$allele_data->{'name'}}{'sequence'} = $allele_seqs->[$i];
		}
	}
	
	return $alleledata;
}


#################################################################################

# Assigns alleles to sequences ($alleledata can be a HASH ref with allele data, a HASH ref with sequences or a FASTA/FASTQ file)
sub match_alleles {

	my ($alleledata,$md5_to_sequence,$md5_to_name,$options,$threads) = @_;

	# Typical allele matching options
	# my $options = { 'alignment' => 'dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 100', 'aligned' => 1, 'ident' => 1 };
	# my $options = { 'alignment' => 'dna match minlen 10 revcomp' };

	if (!defined($options->{'alignment'})){
		# $options->{'alignment'} = 'dna match';
		$options->{'alignment'} = 'dna blastn -evalue 1E-5 -ungapped';
	}
	# Comma-separated values output format
	if ($options->{'alignment'} !~ /-outfmt \d+/){
		$options->{'alignment'} .= ' -outfmt 10';
	}
	# Default value for % of query length aligned
	if (!defined($options->{'aligned'})){
		$options->{'aligned'} = 0.9;
	# Converts thresholds given as percentages to ratios
	} elsif ($options->{'aligned'} =~ /([\d\.]+)%/) {
		$options->{'aligned'} = $1/100;
	}
	# Default value for % of identity in the aligned region
	if (!defined($options->{'ident'})){
		$options->{'ident'} = 1;
	# Converts thresholds given as percentages to ratios
	} elsif ($options->{'ident'} =~ /([\d\.]+)%/) {
		$options->{'ident'} = $1/100;
	}
	# Defines BLASTN extra parameters to run the alignment faster
	if ($options->{'alignment'} =~ /dna blastn/){
		if ($options->{'alignment'} !~ /word_size/){
			$options->{'alignment'} .= ' -word_size 10';
		}
# 		# Some HLA alleles have more than 10 identical matches
# 		if ($options->{'alignment'} !~ /num_alignments/){
# 			$options->{'alignment'} .= ' -num_alignments 10';
# 		}
		if ($options->{'alignment'} !~ /perc_identity/ && defined($options->{'ident'})){
			$options->{'alignment'} .= sprintf(' -perc_identity %d',100*$options->{'ident'});
		}
	}

	my ($unique_seqs, $unique_seq_md5s) = sequences_hash_to_array($md5_to_sequence);

	my ($allele_align_data, $allele_align_info);

	# my $query_file = create_fasta_file($unique_seqs,$unique_seq_md5s);

	# If the input is a sequence file
	if (ref($alleledata) ne "HASH" && -e $alleledata && (is_fasta($alleledata) || is_fastq($alleledata))) {
		#$allele_align_data = align_seqs_from_file($query_file,$alleledata,0,$options->{'alignment'});
		$alleledata = read_allele_file($alleledata,['no warnings']);
	} # else {
		my @allele_names = keys %$alleledata;
		my @allele_seqs;
		if (ref($alleledata->{$allele_names[0]}) eq "HASH" && defined($alleledata->{$allele_names[0]}{'sequence'})){
			@allele_seqs = map $alleledata->{$_}{'sequence'}, @allele_names;
		} else {
			@allele_seqs = map $alleledata->{$_}, @allele_names;
		}
		# my ($allele_seqs, $allele_names) = sequences_hash_to_array(\%alleles);
		if (defined($threads) && $threads>1){
			$allele_align_data = align_seqs_with_threads($unique_seq_md5s,$unique_seqs,\@allele_names,\@allele_seqs,0,$options->{'alignment'},$threads);
		} else {
			$allele_align_data = align_seqs($unique_seq_md5s,$unique_seqs,\@allele_names,\@allele_seqs,0,$options->{'alignment'});
		}
	# }
	# `rm $query_file`;

	# my @previous_assigned_alleles;
	for (my $i=0; $i<=$#{$unique_seq_md5s}; $i++){
		my $seq_md5 = $unique_seq_md5s->[$i];
		my $seq_len = length($unique_seqs->[$i]);
# if ($seq_md5 eq 'TGTGCCAGCAGCCACGAGCAGGGGGAGTCCTACGAGCAGTACTTC'){
# print '';
# }
		if (defined($allele_align_data->{$seq_md5})){
			# If there are several regions of the sequence matching the same allele, sums/merges all of them
			my $allele_align_data_merged = {};
			foreach my $align_data (@{$allele_align_data->{$seq_md5}}) {
				my $allele_name = $align_data->{'NAME'};
				if (!defined($allele_align_data_merged->{$allele_name})){
					$allele_align_data_merged->{$allele_name} = $align_data;
				} else {
					$allele_align_data_merged->{$allele_name}{'IDENT'} += $align_data->{'IDENT'};
					$allele_align_data_merged->{$allele_name}{'ALIGNED'} += $align_data->{'ALIGNED'};
					$allele_align_data_merged->{$allele_name}{'EVALUE'} *= $align_data->{'EVALUE'};
				}
			}
			# Reads the best allele candidates from the merged alignment results
			my (@allele_candidates);
			my $max_ident = 0;
			#for (my $j=0; $j<=$#{$allele_align_data->{$seq_md5}}; $j++) {
			my @sorted_allele_names = sort {$allele_align_data_merged->{$a}{'EVALUE'}<=>$allele_align_data_merged->{$b}{'EVALUE'}} keys %$allele_align_data_merged;
			foreach my $allele_name (@sorted_allele_names){
# if ($allele_name eq 'TRBV12-3*01'){
# print '';
# }
				#my $allele_name = $allele_align_data->{$seq_md5}[$j]{'NAME'};
				my $allele_len = length($alleledata->{$allele_name}{'sequence'});
				my $min_len;
				if ($seq_len <= $allele_len){
					$min_len = $seq_len;
				} else {
					$min_len = $allele_len;
				}

				if (defined($options->{'aligned'}) && $allele_align_data_merged->{$allele_name}{'ALIGNED'} < $options->{'aligned'}*$min_len) {
					next;
				}
				if (defined($options->{'ident'}) && $allele_align_data_merged->{$allele_name}{'IDENT'} < $options->{'ident'}*$allele_align_data_merged->{$allele_name}{'ALIGNED'}) {
					next;
				}
				# First blast alignment is not always the best
				if ($options->{'alignment'} =~ /blast/) {
					# if ($allele_align_data->{$seq_md5}[0]{'EVALUE'} > 1E-5) { last; }
					# Annotates other results with the same or higher identity and aligned length
					if ($max_ident < $allele_align_data_merged->{$allele_name}{'IDENT'}) {
						$max_ident = $allele_align_data_merged->{$allele_name}{'IDENT'};
						#@allele_candidates = ("ident: $max_ident | ".$allele_name);
						if (defined($options->{'amplitaxo'})){
							@allele_candidates = sprintf("%s".";aligned=%d;ident=%d;total=%d;",$allele_name,$allele_align_data_merged->{$allele_name}{'ALIGNED'},$allele_align_data_merged->{$allele_name}{'IDENT'},$seq_len);
						} else {
							@allele_candidates = ($allele_name);
						}
					} elsif ($max_ident == $allele_align_data_merged->{$allele_name}{'IDENT'}) {
						$max_ident = $allele_align_data_merged->{$allele_name}{'IDENT'};
						if (defined($options->{'amplitaxo'})){
							push(@allele_candidates, sprintf("%s".";aligned=%d;ident=%d;total=%d;",$allele_name,$allele_align_data_merged->{$allele_name}{'ALIGNED'},$allele_align_data_merged->{$allele_name}{'IDENT'},$seq_len));
						} else {
							push(@allele_candidates, $allele_name);
						}
					} else {
						next;
					}
				} else {
					# Annotates the best results and other results with the same or higher identity and aligned length
					if (!@allele_candidates || ( $allele_align_data_merged->{$allele_name}{'ALIGNED'} >= $allele_align_data->{$seq_md5}[0]{'ALIGNED'} 
					&& $allele_align_data_merged->{$allele_name}{'IDENT'} >= $allele_align_data->{$seq_md5}[0]{'IDENT'} )) {
# 					&& $allele_align_data_merged->{$allele_name}{'IDENT'} == $allele_align_data_merged->{$allele_name}{'ALIGNED'} ) ) {
						if (defined($options->{'amplitaxo'})){
							push(@allele_candidates, sprintf("%s".";aligned=%d;ident=%d;total=%d;",$allele_name,$allele_align_data_merged->{$allele_name}{'ALIGNED'},$allele_align_data_merged->{$allele_name}{'IDENT'},$seq_len));
						} else {
							push(@allele_candidates, $allele_name);
						}
# 						if ($output_aligned) {
# 							$allele_align_info->{$seq_md5}{$allele_name}{'aligned'} = $allele_align_data_merged->{$allele_name}{'ALIGNED'};
# 						}
# 						if ($output_ident) {
# 							$allele_align_info->{$seq_md5}{$allele_name}{'ident'} = $allele_align_data_merged->{$allele_name}{'IDENT'};
# 						}
					} else {
						last;
					}
				} 
			}
			if (@allele_candidates) {
				if (scalar @allele_candidates > 1){
					@allele_candidates = shuffle @allele_candidates;
				}
				my $allele_name = join(' | ', @allele_candidates);
				$md5_to_name->{$seq_md5} = $allele_name;
				# push(@previous_assigned_alleles,@allele_candidates);
			}
		}
	}
	
# 	if (defined($options->{'output'})){
# 		return ($md5_to_name,$allele_align_info);
# 	} else {
	return $md5_to_name;
# 	}

}

#################################################################################

# Extracts marker/amplicon sequence data
sub retrieve_amplicon_data {

	my ($markers,$samples,$amplicon_sequences,$amplicon_depths,$md5_to_sequence,$md5_to_name,$amplicon_clusters) = @_;

	my ($marker_seq_data,$amplicon_seq_data);
	
	my $unique_seq_number = 0;

	foreach my $marker_name (@$markers){

		if (!defined($amplicon_sequences->{$marker_name})){
			next;
		}

		my ($unique_seq_frequencies, $unique_seq_depths, $per_amplicon_frequencies);

		# Checks only samples with sequences (after filtering) in the original order
		my @sample_names;
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			if (defined($amplicon_sequences->{$marker_name}) && defined($amplicon_sequences->{$marker_name}{$sample_name})){
				push(@sample_names, $sample_name);
			}
		}

		# Loops samples and annotates frequencies
		foreach my $sample_name (@sample_names) {

			my @unique_seqs = keys %{$amplicon_sequences->{$marker_name}{$sample_name}};
			my $total_seqs = $amplicon_depths->{$marker_name}{$sample_name};

			# Loops all the unique sequences in the amplicon
			foreach my $md5 (@unique_seqs) {
				my $unique_seq_depth = $amplicon_sequences->{$marker_name}{$sample_name}{$md5};
				my $unique_seq_frequency = $unique_seq_depth/$total_seqs*100;
				push(@{$unique_seq_frequencies->{$md5}},$unique_seq_frequency);
				push(@{$unique_seq_depths->{$md5}},$unique_seq_depth);
				$per_amplicon_frequencies->{$sample_name}{$md5} = $unique_seq_frequency;
			}
		}

		# Calculates mean, max and min FREQs and also number of samples that contain the same unique sequence
		# These data is for all unique sequences of a single marker
		my %unique_seq_count_samples = map { $_ => scalar @{$unique_seq_frequencies->{$_}} } keys %{$unique_seq_frequencies};
		my %unique_seq_mean_frequencies = map { $_ => mean(@{$unique_seq_frequencies->{$_}}) } keys %{$unique_seq_frequencies};
		my %unique_seq_max_frequencies = map { $_ => max(@{$unique_seq_frequencies->{$_}}) } keys %{$unique_seq_frequencies};
		my %unique_seq_min_frequencies = map { $_ => min(@{$unique_seq_frequencies->{$_}}) } keys %{$unique_seq_frequencies};
		my %unique_seq_sum_depths = map { $_ => sum(@{$unique_seq_depths->{$_}}) } keys %{$unique_seq_depths};

		# Annotates unique sequences for a single marker
		# 1. Sorts unique sequences by depth (sum of all samples depths)
		# 2. Gives a name to the each sorted sequence
		# 3. Annotates sequence parameters (hash, length, depth, samples, FREQ...)
		my @sorted_depth_unique_seqs = sort { $unique_seq_sum_depths{$b} <=> $unique_seq_sum_depths{$a} } keys %unique_seq_sum_depths;
		my $zeros = length(scalar @sorted_depth_unique_seqs);
		foreach my $md5 (@sorted_depth_unique_seqs) {
			my $name;
			# Checks if the sequence already has a name (eg. allele name)
			if (defined($md5_to_name->{$md5})){
				$name = $md5_to_name->{$md5};
				# Leaves empty this sequence number
				$unique_seq_number++;
			# If not, creates a new name with an autoincrement number
			} else {
				$unique_seq_number++;
				$name = sprintf("%s-%0".$zeros."d", $marker_name, $unique_seq_number);
				$md5_to_name->{$md5} = $name;
			}
			my $seq = $md5_to_sequence->{$md5};
			my $len = length($seq);
			my $unique_seq_depth = $unique_seq_sum_depths{$md5};
			my $count_samples = $unique_seq_count_samples{$md5};
			my $mean_freq = $unique_seq_mean_frequencies{$md5};
			my $min_freq = $unique_seq_min_frequencies{$md5};
			my $max_freq = $unique_seq_max_frequencies{$md5};
			$marker_seq_data->{$marker_name}{$md5} = { 'seq'=> $seq, 'name'=>$name, 'len'=>$len, 'depth'=>$unique_seq_depth, 'samples'=>$count_samples, 'mean_freq'=>$mean_freq, 'max_freq'=>$max_freq, 'min_freq'=>$min_freq };
		}

		# Annotates unique sequences for a single amplicon (marker+sample)
		# 1. Sorts unique sequences by depth
		# 2. Gives a name to the each sorted sequence
		# 3. Annotates sequence parameters (hash, length, depth, samples, FREQ...)
		foreach my $sample_name (@sample_names) {
			# Sort by depth
			my @sorted_unique_seqs =  sort { $amplicon_sequences->{$marker_name}{$sample_name}{$b} <=> $amplicon_sequences->{$marker_name}{$sample_name}{$a} }  keys %{$amplicon_sequences->{$marker_name}{$sample_name}};
# 			if (defined($amplicon_depths->{$marker_name}{$sample_name})){
			# my $total_seqs = $amplicon_depths->{$marker_name}{$sample_name};
			foreach my $md5 (@sorted_unique_seqs) {
				# 'name', 'seq', 'len', 'count_samples', 'mean_freq', 'min_freq' and 'max_freq' are already annotated in $marker_seq_data
				# my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
				# my $seq = $md5_to_sequence->{$md5};
				# my $len = length($seq);
				my $unique_seq_depth = $amplicon_sequences->{$marker_name}{$sample_name}{$md5};
				my $unique_seq_frequency = $per_amplicon_frequencies->{$sample_name}{$md5};
				# my $count_samples = $unique_seq_count_samples{$md5};
				# my $mean_freq = $unique_seq_mean_frequencies{$md5};
				# my $min_freq = $unique_seq_min_frequencies{$md5};
				# my $max_freq = $unique_seq_max_frequencies{$md5};
				if (defined($amplicon_clusters->{$marker_name}{$sample_name}{$md5})){
					$amplicon_seq_data->{$marker_name}{$sample_name}{$md5} = { 'depth'=>$unique_seq_depth, 'freq'=>$unique_seq_frequency, 'cluster_size'=>$amplicon_clusters->{$marker_name}{$sample_name}{$md5} };
				} else {
					$amplicon_seq_data->{$marker_name}{$sample_name}{$md5} = { 'depth'=>$unique_seq_depth, 'freq'=>$unique_seq_frequency };
				}
			}
		}
	}
	return ($marker_seq_data, $amplicon_seq_data, $md5_to_name);

}

#################################################################################

# Prints an Excel file with genotyping results, with each amplicon in a different Spreadsheet
# Also prints the genotyping results into equivalent TXT files, one per amplicon
# And prints all amplicon sequences in one FASTA file
sub print_marker_sequences {

	my ($markers,$samples,$marker_seq_data,$amplicon_seq_data,$amplicon_depths,$outpath,$amplicon_raw_sequences) = @_; #,$STC_seq_data,$amplicon_raw_sequences) = @_;

	my ($marker_seq_files,$marker_matrix_files);

	my $marker_result_file = "$outpath/results.xlsx";
	my $workbook  = Excel::Writer::XLSX->new($marker_result_file);
	$workbook->set_properties(
		title    => "AmpliSAT results",
		author   => "Alvaro Sebastian",
		comments => "AmpliSAT results",
		company  => "SixthResearcher",
	);
	$workbook->compatibility_mode();
	my $bold = $workbook->add_format(bold => 1);
	my $red = $workbook->add_format(bg_color => 'red');
	my $green = $workbook->add_format(bg_color => 'green');
	my $blue = $workbook->add_format(bg_color => 'blue');
	my $yellow = $workbook->add_format(bg_color => 'yellow');
	my $magenta = $workbook->add_format(bg_color => 'magenta');
	my $cyan = $workbook->add_format(bg_color => 'cyan');

	foreach my $marker_name (@$markers){

		# Skips markers without data
		if (!defined($marker_seq_data->{$marker_name})){
			next;
		}

		# Checks only samples with sequences (after filtering) in the original order
		my (@sample_names, @sample_amplicon_depths, @sample_allele_depths, @sample_allele_counts);
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			if (defined($amplicon_seq_data->{$marker_name}{$sample_name})){
				push(@sample_names, $sample_name);
				push(@sample_amplicon_depths, $amplicon_depths->{$marker_name}{$sample_name});
				push(@sample_allele_depths, sum( map { $amplicon_seq_data->{$marker_name}{$sample_name}{$_}{'depth'} } keys %{$amplicon_seq_data->{$marker_name}{$sample_name}} ));
				push(@sample_allele_counts, scalar keys %{$amplicon_seq_data->{$marker_name}{$sample_name}});
			}
		}

		# Writes FASTA files with all unique sequences for a single marker
		my @sorted_depth_unique_seqs = sort { $marker_seq_data->{$marker_name}{$b}{'depth'} <=> $marker_seq_data->{$marker_name}{$a}{'depth'} } keys %{$marker_seq_data->{$marker_name}};
		my (@unique_seq_headers, @unique_seq_seqs);
		foreach my $md5 (@sorted_depth_unique_seqs) {
			my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
			my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
			my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
			my $depth = $marker_seq_data->{$marker_name}{$md5}{'depth'};
			my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
			my $mean_freq = $marker_seq_data->{$marker_name}{$md5}{'mean_freq'};
			my $max_freq = $marker_seq_data->{$marker_name}{$md5}{'max_freq'};
			my $min_freq = $marker_seq_data->{$marker_name}{$md5}{'min_freq'};
			push(@unique_seq_headers, sprintf("%s | hash=%s | len=%d | depth=%d | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $name, $md5, $len, $depth, $count_samples, $mean_freq, $max_freq, $min_freq));
			push(@unique_seq_seqs, $seq);
		}
		create_fasta_file(\@unique_seq_seqs,\@unique_seq_headers,"$outpath/$marker_name.fasta");
		$marker_seq_files->{$marker_name} = "$outpath/$marker_name.fasta";

		# Writes matrix of unique sequences vs. samples for a single marker
		my $worksheet = $workbook->add_worksheet("$marker_name");
		$worksheet->set_row(3, undef, $bold);
		$worksheet->set_column('I:I', undef, $bold);

		# Writes sample amplicon and allele depths, and alleles occurences
		$worksheet->write_row(0, 8, ['DEPTH_AMPLICON', @sample_amplicon_depths]);
		$worksheet->write(1, 8, 'DEPTH_ALLELES');
		$worksheet->write(2, 8, 'COUNT_ALLELES');
		$worksheet->write_row(3, 0, ['SEQUENCE', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'MEAN_FREQ', 'MAX_FREQ', 'MIN_FREQ', '',@sample_names]);
		my $ws_col = 9;
		my $ws_row_first = 4;
		my $ws_row_last = $ws_row_first+$#sorted_depth_unique_seqs;
		for (my $i=0; $i<=$#sample_names; $i++) {
			$worksheet->write(0, $ws_col, $sample_amplicon_depths[$i]);
			$worksheet->write_formula(1, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $sample_allele_depths[$i]);
			$worksheet->write_formula(2, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $sample_allele_counts[$i]);
			$worksheet->write(3, $ws_col, $sample_names[$i]);
			$ws_col++;
		}
		my $ws_row = 4;
# 		my @freq_matrix_headers = ('NAME', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'MEAN_FREQ', 'MAX_FREQ', 'MIN_FREQ', @sample_names);
		my @freq_matrix_headers = ('NAME', 'SEQ', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', @sample_names);
		my $freq_matrix = join("\t",@freq_matrix_headers)."\n";
# 		my $freq_matrix = sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t", 'NAME', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'MEAN_FREQ', 'MAX_FREQ', 'MIN_FREQ');
# 		$freq_matrix .= join("\t",@sample_names)."\n";
		# Annotate sequence depths for each unique sequence in each sample
		foreach my $md5 (@sorted_depth_unique_seqs) {
			my (@depths_formulas, @depths);
			foreach my $sample_name (@sample_names) {
				if (defined($amplicon_seq_data->{$marker_name}{$sample_name}{$md5}) && defined($amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'})){
					push(@depths,$amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'});
					if (defined($amplicon_raw_sequences) && defined($amplicon_raw_sequences->{$marker_name}{$sample_name}{$md5})){
						push(@depths_formulas,sprintf("=%d+%d", $amplicon_raw_sequences->{$marker_name}{$sample_name}{$md5}, $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'}-$amplicon_raw_sequences->{$marker_name}{$sample_name}{$md5}));
					} else {
						push(@depths_formulas,undef);
					}
				} else {
					push(@depths,'');
					push(@depths_formulas,undef);
				}
			}
			$ws_col = 0;
			my $ws_col_first = 9;
			my $ws_col_last = $ws_col_first+$#sample_names;
			my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
			my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
			my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
			my $depth = $marker_seq_data->{$marker_name}{$md5}{'depth'};
			my $depth_formula = sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last));
			my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
			my $count_samples_formula = sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last));
			my $mean_freq = sprintf("%.2f",$marker_seq_data->{$marker_name}{$md5}{'mean_freq'});
			my $max_freq = sprintf("%.2f",$marker_seq_data->{$marker_name}{$md5}{'max_freq'});
			my $min_freq = sprintf("%.2f",$marker_seq_data->{$marker_name}{$md5}{'min_freq'});
			$worksheet->write($ws_row, $ws_col, $seq);$ws_col++;
			$worksheet->write($ws_row, $ws_col, $md5);$ws_col++;
			$worksheet->write($ws_row, $ws_col, $len);$ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, $depth_formula, undef, $depth);$ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, $count_samples_formula, undef, $count_samples);$ws_col++;
			$worksheet->write($ws_row, $ws_col, $mean_freq);$ws_col++;
			$worksheet->write($ws_row, $ws_col, $max_freq);$ws_col++;
			$worksheet->write($ws_row, $ws_col, $min_freq);$ws_col++;
			$worksheet->write($ws_row, $ws_col, $name);$ws_col++;
			for (my $i=0; $i<=$#depths; $i++){
				if (defined($depths_formulas[$i])){
					$worksheet->write_formula($ws_row, $ws_col, $depths_formulas[$i], undef, $depths[$i]);$ws_col++;
				} else {
					$worksheet->write($ws_row, $ws_col, $depths[$i]);$ws_col++;
				}
			}
# 			my @freq_matrix_values = ($name, $md5, $len, $depth, $count_samples, sprintf("%.2f",$mean_freq), sprintf("%.2f",$max_freq), sprintf("%.2f",$min_freq), @depths);
			my @freq_matrix_values = ($name, $seq, $md5, $len, $depth, $count_samples, @depths);
# 			my @worksheet_values = ($seq, $md5, $len, $depth, $count_samples, sprintf("%.2f",$mean_freq), sprintf("%.2f",$max_freq), sprintf("%.2f",$min_freq), $name, @depths);
			$freq_matrix .= join("\t",@freq_matrix_values)."\n";
# 			$worksheet->write_row($ws_row, 0, \@worksheet_values);
			$ws_row++;
# 			$freq_matrix .= sprintf("%s\t%s\t%s\t%d\t%d\t%.2f\t%.2f\t%.2f\t%s\n", $name, $md5, $len, $depth, $count_samples, $mean_freq, $max_freq, $min_freq, join("\t",@depths));
			# push(@unique_seq_headers,sprintf("%s | hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | mean_freq=%.2f | min_freq=%.2f | max_freq=%.2f", $name, $md5, $len, $unique_seq_depth, $unique_seq_frequency, $count_samples, $mean_freq, $min_freq, $max_freq));
			
		}
		$freq_matrix .= "\n";
		write_to_file("$outpath/$marker_name.txt", $freq_matrix);
		$marker_matrix_files->{$marker_name} = "$outpath/$marker_name.txt";
	}
	$workbook->close();
	
	return ($marker_result_file,$marker_seq_files,$marker_matrix_files);

}

#################################################################################

# Prints amplicon sequences in individual FASTA files 
sub print_amplicon_sequences {

	my ($markers,$samples,$marker_seq_data,$amplicon_seq_data,$outpath) = @_;

	my ($amplicon_seq_files);

	foreach my $marker_name (@$markers){

		# Skips markers without data
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}
		
		# Checks only samples with sequences (after filtering) in the original order
		my @sample_names;
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			if (defined($amplicon_seq_data->{$marker_name}{$sample_name})){
				push(@sample_names, $sample_name);
			}
		}
		
		# Writes FASTA files with all unique sequences for a single amplicon (marker+sample)
		foreach my $sample_name (@sample_names) {
			my @sorted_depth_amplicon_unique_seqs =  sort { $amplicon_seq_data->{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->{$marker_name}{$sample_name}{$a}{'depth'} }  keys %{$amplicon_seq_data->{$marker_name}{$sample_name}};
			my (@amplicon_unique_seq_headers, @amplicon_unique_seq_seqs);
			foreach my $md5 (@sorted_depth_amplicon_unique_seqs) {
				my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
				my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
				my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
				my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
				my $frequency = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'};
				my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
				my $mean_freq = $marker_seq_data->{$marker_name}{$md5}{'mean_freq'};
				my $min_freq = $marker_seq_data->{$marker_name}{$md5}{'min_freq'};
				my $max_freq = $marker_seq_data->{$marker_name}{$md5}{'max_freq'};
				push(@amplicon_unique_seq_headers, sprintf("%s | hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $name, $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq));
				push(@amplicon_unique_seq_seqs, $seq);
			}
			create_fasta_file(\@amplicon_unique_seq_seqs,\@amplicon_unique_seq_headers,"$outpath/$marker_name-$sample_name.fasta");
			$amplicon_seq_files->{$marker_name}{$sample_name} = "$outpath/$marker_name-$sample_name.fasta";
		}
	}
	
	return $amplicon_seq_files;

}

#################################################################################

# Prints an excel file with a comparison of amplicon sequences retrieved by 'compare_amplicon_sequences'
# This function is a modified version of 'print_marker_sequences'
sub print_comparison_sequences {

	my ($markers,$samples,$marker_seq_data,$amplicon_seq_data,$amplicon_depths,$comparison_result_file,$options) = @_; #,$STC_seq_data,$amplicon_raw_sequences) = @_;

	my $expand_results = 0;
	if (in_array($options, 'expand results')){
		$expand_results = 1;
	}

	my $workbook  = Excel::Writer::XLSX->new($comparison_result_file);
	$workbook->set_properties(
		title    => "AmpliSAS results",
		author   => "Alvaro Sebastian",
		comments => "AmpliSAS results",
		company  => "SixthResearcher",
	);
	$workbook->compatibility_mode();
	my $bold = $workbook->add_format(bold => 1);
	my $red = $workbook->add_format(bg_color => 'red');
	my $red_bold = $workbook->add_format(bg_color => 'red', bold => 1);
	my $green = $workbook->add_format(bg_color => 'green');
	my $green_bold = $workbook->add_format(bg_color => 'green', bold => 1);
	my $blue = $workbook->add_format(bg_color => 'blue');
	my $yellow = $workbook->add_format(bg_color => 'yellow');
	my $yellow_bold = $workbook->add_format(bg_color => 'yellow', bold => 1);
	my $magenta = $workbook->add_format(bg_color => 'magenta');
	my $cyan = $workbook->add_format(bg_color => 'cyan');

	foreach my $marker_name (@$markers){

		# Skips markers without data
		if (!defined($marker_seq_data->{$marker_name})){
			next;
		}

		# Checks only samples with sequences (after filtering) in the original order
		my (@sample_names, @sample_amplicon_depths, @sample_allele_depths, @sample_allele_freqs, @sample_allele_counts);
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			if (defined($amplicon_seq_data->{$marker_name}{$sample_name})){
				push(@sample_names, $sample_name);
				push(@sample_amplicon_depths, $amplicon_depths->{$marker_name}{$sample_name});
				push(@sample_allele_depths, sum( map { $amplicon_seq_data->{$marker_name}{$sample_name}{$_}{'depth'} } keys %{$amplicon_seq_data->{$marker_name}{$sample_name}} ));
				push(@sample_allele_freqs, sum( map { $amplicon_seq_data->{$marker_name}{$sample_name}{$_}{'freq'} } keys %{$amplicon_seq_data->{$marker_name}{$sample_name}} ));
				push(@sample_allele_counts, scalar keys %{$amplicon_seq_data->{$marker_name}{$sample_name}});
			}
		}

		# Annotates sequences sample by sample
		# Stores comparisons to be annotated, when the sequence is similar to a previous one
		# '$mark_sequences' marks sequences as good (0) suspicious (1) and probably bad (2)
		my ($comparisons, $mark_sequences);
		foreach my $sample_name (@sample_names) {
			# Order sequences by depth
			my @sorted_depth_unique_sample_seqs = sort { $amplicon_seq_data->{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->{$marker_name}{$sample_name}{$a}{'depth'} } keys %{$amplicon_seq_data->{$marker_name}{$sample_name}};
			# If a new sequence is similar to any of the previous ones, then their comparisons with others will be annotated
			my %previous_seqs;
			foreach my $md5 (@sorted_depth_unique_sample_seqs) {
				if (defined($amplicon_seq_data->{$marker_name}{$sample_name}{$md5}) && defined($amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'similar_seq_comparisons'})){
					$mark_sequences->{$sample_name}{$md5} = 0;
					$comparisons->{$sample_name}{$md5} = [];
					# Annotates chimeras
					if (@{$amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'similar_seq_chimeras'}}){
						# Mark in red as bad sequence (2) if is chimera from previous annotated seqs
						$mark_sequences->{$sample_name}{$md5} = 2;
						my @similar_seq_chimeras = map "CH: $_", @{$amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'similar_seq_chimeras'}};
						# Mark in red as a bad sequence (2)
						$mark_sequences->{$sample_name}{$md5} = 2;
						push(@{$comparisons->{$sample_name}{$md5}}, @similar_seq_chimeras);
					}
					# Annotates sequence errors
					if (@{$amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'similar_seq_md5s'}}){
						my @similar_seq_md5s = @{$amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'similar_seq_md5s'}};
						my @similar_seq_names = @{$amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'similar_seq_names'}};
						my @similar_seq_comparisons = @{$amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'similar_seq_comparisons'}};
						for (my $i=0; $i<=$#similar_seq_md5s; $i++) {
							# Mark as a suspicious (!) or bad sequence (2) if the errors come from a previous sequence 
							if (defined($previous_seqs{$similar_seq_md5s[$i]})){
								# Mark in red as a bad sequence (2), if the errors come from a previous sequence with much higher depth and there are no errors already detected
								if (!defined($mark_sequences->{$sample_name}{$md5}) && $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'}<0.25*$amplicon_seq_data->{$marker_name}{$sample_name}{$similar_seq_md5s[$i]}{'depth'}){
									$mark_sequences->{$sample_name}{$md5} = 2;
								# Mark in yellow as a suspicious sequence (1), if the errors come from a previous sequence with not very much higher depth
								} elsif (!defined($mark_sequences->{$sample_name}{$md5})){
									$mark_sequences->{$sample_name}{$md5} = 1;
								}
								push(@{$comparisons->{$sample_name}{$md5}}, sprintf("%s: %s", $similar_seq_comparisons[$i], $similar_seq_names[$i]));
							}
# 							push(@{$comparisons->{$sample_name}{$md5}}, sprintf("%s: %s", $similar_seq_comparisons[$i], $similar_seq_names[$i]));
						}
					}
				}
				$previous_seqs{$md5} = 1;
			}
		}

		# Order sequences by global depth
		my @sorted_depth_unique_seqs = sort { $marker_seq_data->{$marker_name}{$b}{'depth'} <=> $marker_seq_data->{$marker_name}{$a}{'depth'} } keys %{$marker_seq_data->{$marker_name}};

		# Writes matrix of unique sequences vs. samples for a single marker
		my ($worksheet, $worksheet1, $worksheet2, $worksheet3);
		if (!$expand_results) {
			$worksheet = $workbook->add_worksheet("$marker_name");
		} else {
			$worksheet1 = $workbook->add_worksheet("$marker_name\_depths");
			$worksheet2 = $workbook->add_worksheet("$marker_name\_freqs");
			$worksheet3 = $workbook->add_worksheet("$marker_name\_errors");
		}

		my $ws_row = 3;
		if (!$expand_results) {
			$worksheet->set_row($ws_row, undef, $bold);
			$worksheet->set_column('I:I', undef, $bold);
		} else {
			$worksheet1->set_row($ws_row, undef, $bold);
			$worksheet1->set_column('I:I', undef, $bold);
			$worksheet2->set_row($ws_row, undef, $bold);
			$worksheet2->set_column('I:I', undef, $bold);
			$worksheet3->set_row($ws_row, undef, $bold);
			$worksheet3->set_column('I:I', undef, $bold);
		}

		# Writes color legend
		$ws_row = 1;
		if (!$expand_results) {
			$worksheet->write($ws_row, 0, "Good sequence", $green_bold);
			$worksheet->write($ws_row, 1, "Suspicious sequence", $yellow_bold);
			$worksheet->write($ws_row, 2, "Probably artifact", $red_bold);
		} else {
			$worksheet1->write($ws_row, 0, "Good sequence", $green_bold);
			$worksheet1->write($ws_row, 1, "Suspicious sequence", $yellow_bold);
			$worksheet1->write($ws_row, 2, "Probably artifact", $red_bold);
			$worksheet2->write($ws_row, 0, "Good sequence", $green_bold);
			$worksheet2->write($ws_row, 1, "Suspicious sequence", $yellow_bold);
			$worksheet2->write($ws_row, 2, "Probably artifact", $red_bold);
			$worksheet3->write($ws_row, 0, "Good sequence", $green_bold);
			$worksheet3->write($ws_row, 1, "Suspicious sequence", $yellow_bold);
			$worksheet3->write($ws_row, 2, "Probably artifact", $red_bold);
		}

		# Writes sample amplicon and allele depths, and alleles occurences
		my @worksheet_headers1 = ('', '', '', '', '', '', '', '','DEPTH_AMPLICON', @sample_amplicon_depths);
		my @worksheet_headers2 = ('', '', '', '', '', '', '', '','DEPTH_SEQUENCES'); # ,@sample_allele_depths);
		my @worksheet_headers3 = ('', '', '', '', '', '', '', '', 'COUNT_SEQUENCES'); # ,@sample_allele_counts);
		my @worksheet_headers4 = ('SEQUENCE', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'MEAN_FREQ', 'MAX_FREQ', 'MIN_FREQ', '',@sample_names);
		if (!$expand_results) {
			$worksheet->write_row(0, 0, \@worksheet_headers1);
			$worksheet->write_row(1, 0, \@worksheet_headers2);
			$worksheet->write_row(2, 0, \@worksheet_headers3);
			$worksheet->write_row(3, 0, \@worksheet_headers4);
		} else {
			my @worksheet1_headers2 = ('', '', '', '', '', '', '', '','DEPTH_SEQUENCES'); # ,@sample_allele_depths);
			my @worksheet2_headers2 = ('', '', '', '', '', '', '', '','FREQ_SEQUENCES');
			my @worksheet3_headers3 = ('', '', '', '', '', '', '', '', 'COUNT_SEQUENCES');
			$worksheet1->write_row(0, 0, \@worksheet_headers1);
			$worksheet1->write_row(1, 0, \@worksheet1_headers2);
			$worksheet1->write_row(2, 0, \@worksheet_headers3);
			$worksheet1->write_row(3, 0, \@worksheet_headers4);
			$worksheet2->write_row(0, 0, \@worksheet_headers1);
			$worksheet2->write_row(1, 0, \@worksheet2_headers2);
			$worksheet2->write_row(2, 0, \@worksheet_headers3);
			$worksheet2->write_row(3, 0, \@worksheet_headers4);
			$worksheet3->write_row(1, 0, \@worksheet_headers1);
			#$worksheet3->write_row(1, 0, \@worksheet_headers2_);
			$worksheet3->write_row(2, 0, \@worksheet3_headers3);
			$worksheet3->write_row(3, 0, \@worksheet_headers4);
		}
		my $ws_col = 9;
		$ws_row = 4;
		for (my $i=0; $i<=$#sample_names; $i++) {
			if (!$expand_results) {
				$worksheet->write(1, $ws_col, $sample_allele_depths[$i]);
				$worksheet->write(2, $ws_col, $sample_allele_counts[$i]);
			} else {
				$worksheet1->write_formula(1, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col), xl_rowcol_to_cell($ws_row+$#sorted_depth_unique_seqs,$ws_col)), undef, $sample_allele_depths[$i]);
				$worksheet1->write_formula(2, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col), xl_rowcol_to_cell($ws_row+$#sorted_depth_unique_seqs,$ws_col)), undef, $sample_allele_counts[$i]);
				$worksheet2->write_formula(1, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col), xl_rowcol_to_cell($ws_row+$#sorted_depth_unique_seqs,$ws_col)), undef, $sample_allele_freqs[$i]);
				$worksheet2->write_formula(2, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col), xl_rowcol_to_cell($ws_row+$#sorted_depth_unique_seqs,$ws_col)), undef, $sample_allele_counts[$i]);
				$worksheet3->write_formula(2, $ws_col, sprintf('=COUNTIF(%s:%s,"<>-")', xl_rowcol_to_cell($ws_row,$ws_col), xl_rowcol_to_cell($ws_row+$#sorted_depth_unique_seqs,$ws_col)), undef, $sample_allele_counts[$i]);
			}
			$ws_col++;
		}
		$ws_col = 9;
		$ws_row = 4;


		# Prints each ordered sequence in one row of Excel file
		foreach my $md5 (@sorted_depth_unique_seqs) {
# if ($md5 eq '1cce5ae72968c419e6f7d64cc1ba6104') {
# print '';
# }
			my $ws_col = 9;
			my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
			my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
			my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
			my $seq_depth = $marker_seq_data->{$marker_name}{$md5}{'depth'};
			my $depth_ws1 = sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col), xl_rowcol_to_cell($ws_row,$ws_col+$#sample_names));
			my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
			my $count_samples_ws1 = sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col), xl_rowcol_to_cell($ws_row,$ws_col+$#sample_names));
			my $count_samples_ws3 = sprintf('=COUNTIF(%s:%s,"<>-")', xl_rowcol_to_cell($ws_row,$ws_col), xl_rowcol_to_cell($ws_row,$ws_col+$#sample_names));
			my $mean_freq = $marker_seq_data->{$marker_name}{$md5}{'mean_freq'};
			my $max_freq = $marker_seq_data->{$marker_name}{$md5}{'max_freq'};
			my $min_freq = $marker_seq_data->{$marker_name}{$md5}{'min_freq'};
			if (!$expand_results) {
				$worksheet->write_row($ws_row, 0, [$seq, $md5, $len, $seq_depth, $count_samples]);
				$worksheet->write_row($ws_row, 5, [sprintf("%.2f",$mean_freq), sprintf("%.2f",$max_freq), sprintf("%.2f",$min_freq)]);
			} else {
				$worksheet1->write_row($ws_row, 0, [$seq, $md5, $len]);
				$worksheet2->write_row($ws_row, 0, [$seq, $md5, $len, $seq_depth]);
				$worksheet3->write_row($ws_row, 0, [$seq, $md5, $len, $seq_depth]);
				$worksheet1->write_formula($ws_row, 3, $depth_ws1, undef, $seq_depth);
				$worksheet1->write_formula($ws_row, 4, $count_samples_ws1, undef, $count_samples);
				$worksheet2->write_formula($ws_row, 4, $count_samples_ws1, undef, $count_samples);
				$worksheet3->write_formula($ws_row, 4, $count_samples_ws3, undef, $count_samples);
				$worksheet1->write_row($ws_row, 5, [sprintf("%.2f",$mean_freq), sprintf("%.2f",$max_freq), sprintf("%.2f",$min_freq)]);
				$worksheet2->write_row($ws_row, 5, [sprintf("%.2f",$mean_freq), sprintf("%.2f",$max_freq), sprintf("%.2f",$min_freq)]);
				$worksheet3->write_row($ws_row, 5, [sprintf("%.2f",$mean_freq), sprintf("%.2f",$max_freq), sprintf("%.2f",$min_freq)]);
			}
			my ($is_good_seq, $is_bad_seq, $is_suspicious_seq) = (0, 0, 0);
			# foreach my $sample_name (@sample_names) {
			for (my $i=0; $i<=$#sample_names; $i++) {
				my $sample_name = $sample_names[$i];
				if (defined($mark_sequences->{$sample_name}{$md5})){
					my $style;
					if ($mark_sequences->{$sample_name}{$md5} == 2) {
						$style = $red;
						$is_bad_seq++;
					} elsif ($mark_sequences->{$sample_name}{$md5} == 1) {
						$style = $yellow;
						$is_suspicious_seq++;
					} else {
						$is_good_seq++;
					}
					my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
					my $freq = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'};
					if (!$expand_results) {
						my $zeros = length($sample_allele_depths[$i]);
						if (@{$comparisons->{$sample_name}{$md5}}){
							$worksheet->write($ws_row, $ws_col, sprintf("%0".$zeros."d; %04.2f%%; %s" , $depth, $freq, join(", ", @{$comparisons->{$sample_name}{$md5}})), $style);
						} else {
	# 						$worksheet->write($ws_row, $ws_col, sprintf("%s; %s%%; %s" , $depth, $freq, 'NO ERRORS FOUND', $style);
							$worksheet->write($ws_row, $ws_col, sprintf("%0".$zeros."d; %04.2f%%" , $depth, $freq), $style);
						}
					} else {
						$worksheet1->write($ws_row, $ws_col, $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'}, $style);
						$worksheet2->write($ws_row, $ws_col, $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'}, $style);
						if (@{$comparisons->{$sample_name}{$md5}}){
							$worksheet3->write($ws_row, $ws_col, join(", ", @{$comparisons->{$sample_name}{$md5}}), $style);
						} else {
							$worksheet3->write($ws_row, $ws_col, 'NO ERRORS FOUND');
						}
					}
				} else {
					# Makes the spreadsheet easier to read the
					if (!$expand_results) {
						$worksheet->write($ws_row, $ws_col, ' ');
					} else {
						$worksheet3->write($ws_row, $ws_col, ' ');
					}
				}
				$ws_col++;
				# If the sequence looks good seq in most of the samples
				if ($is_good_seq > 5*($is_bad_seq+$is_suspicious_seq) ) {
					if (!$expand_results) {
						$worksheet->write($ws_row, 8, $name, $green);
					} else {
						$worksheet1->write($ws_row, 8, $name, $green);
						$worksheet2->write($ws_row, 8, $name, $green);
						$worksheet3->write($ws_row, 8, $name, $green);
					}
				# If the sequence looks good seq in some samples and artifact in others
				} elsif (($is_good_seq+$is_suspicious_seq) > $is_bad_seq) {
					if (!$expand_results) {
						$worksheet->write($ws_row, 8, $name, $yellow);
					} else {
						$worksheet1->write($ws_row, 8, $name, $yellow);
						$worksheet2->write($ws_row, 8, $name, $yellow);
						$worksheet3->write($ws_row, 8, $name, $yellow);
					}
				# If the sequence looks bad seq in all the samples
				} else {
					if (!$expand_results) {
						$worksheet->write($ws_row, 8, $name, $red);
					} else {
						$worksheet1->write($ws_row, 8, $name, $red);
						$worksheet2->write($ws_row, 8, $name, $red);
						$worksheet3->write($ws_row, 8, $name, $red);
					}
				}
			}
			# Print empty space in next column to the last one
			if (!$expand_results) {
				$worksheet->write($ws_row, $ws_col, ' ');
			} else {
				$worksheet3->write($ws_row, $ws_col, ' ');
			}
			$ws_row++;
		}
		# Writes error legend
		$ws_row++;
		if (!$expand_results) {
			$worksheet->write($ws_row, 0, "LEGEND:", $bold); $ws_row++;
			$worksheet->write($ws_row, 0, "Depth; Frequency; Putative errors", $bold); $ws_row++; $ws_row++;
			$worksheet->write($ws_row, 0, "ERRORS:", $bold); $ws_row++;
			$worksheet->write($ws_row, 0, "CH = Chimera", $red_bold); $ws_row++;
			$worksheet->write($ws_row, 0, "X = Substitutions", $red_bold); $ws_row++;
			$worksheet->write($ws_row, 0, "I = Insertions", $red_bold); $ws_row++;
			$worksheet->write($ws_row, 0, "D = Deletions", $red_bold); $ws_row++;
			$worksheet->write($ws_row, 0, "H = Homopolymer indels", $red_bold); $ws_row++;
		} else {
			$worksheet3->write($ws_row, 0, "LEGEND:", $bold); $ws_row++;
			$worksheet3->write($ws_row, 0, "CH = Chimera", $bold); $ws_row++;
			$worksheet3->write($ws_row, 0, "X = Substitutions", $bold); $ws_row++;
			$worksheet3->write($ws_row, 0, "I = Insertions", $bold); $ws_row++;
			$worksheet3->write($ws_row, 0, "D = Deletions", $bold); $ws_row++;
			$worksheet3->write($ws_row, 0, "H = Homopolymer indels", $bold); $ws_row++;
		}
	}

	$workbook->close();
	
	return $comparison_result_file;

}

#################################################################################

# Annotates low coverage sample/amplicons
sub annotate_low_depth {

	my ($markers,$samples,$amplicon_depths,$min_amplicon_depth) = @_;

	print "\nAnnotating low coverage sample/amplicon.\n";

	# Annotate low coverage sample/amplicon
	my $low_depth;
	foreach my $marker_name (@$markers){
		foreach my $sample_name (@{$samples->{$marker_name}}){
			if (defined($amplicon_depths->{$marker_name}{$sample_name})){
				if ($amplicon_depths->{$marker_name}{$sample_name} < $min_amplicon_depth){
					$low_depth->{$marker_name}{$sample_name} = 1;
				}
			} else {
				$low_depth->{$marker_name}{$sample_name} = 1;
			}
		}
	}
	return $low_depth;
}

#################################################################################

# # Filters alleles according to thresholds
sub filter_amplicon_sequences {

	my ($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$amplicon_depths) = @_;

	my ($filtered_amplicon_sequences, $filtered_amplicon_depths, $filters_output,$md5_to_sequence);
	foreach my $marker_name (@$markers){

		# Skips markers without data
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}

		# Exclude amplicons not in the list
		if (defined($paramsdata->{'allowed_markers'}) && !defined($paramsdata->{'allowed_markers'}{'all'}) && !in_array($paramsdata->{'allowed_markers'},$marker_name)){
			next;
		}

		# Extract correct marker lengths
		my @marker_lengths;
		if (defined($markerdata->{$marker_name}) && defined($markerdata->{$marker_name}{'length'})){
			@marker_lengths = @{$markerdata->{$marker_name}{'length'}};
		}

		# Extracts filtering parameters
		my $filtering_parameters;
		foreach my $param (keys %$paramsdata) {
			if (defined($paramsdata->{$param}{$marker_name})){
				$filtering_parameters->{$param} = $paramsdata->{$param}{$marker_name}[0];
			} elsif (defined($paramsdata->{$param}{'all'})){
				$filtering_parameters->{$param} = $paramsdata->{$param}{'all'}[0];
			}
		}
		# Threshold for parental sequences identity (respect to chimera seq) in chimera detection
		my $max_chimera_ident;
		if (defined($filtering_parameters->{'min_chimera_length'}) && defined($filtering_parameters->{'substitution_threshold'}) && defined($filtering_parameters->{'indel_threshold'})){
			if ($filtering_parameters->{'substitution_threshold'}>$filtering_parameters->{'indel_threshold'}){
				$max_chimera_ident = $filtering_parameters->{'substitution_threshold'}."%";
			} else {
				$max_chimera_ident = $filtering_parameters->{'indel_threshold'}."%";
			}
		}

# 		SAMPLE NUMBER DATA ALREADY INCLUDED IN VARIABLE '$amplicon_seq_data'
# 		# Loops samples and annotates in how many samples is each unique sequence
# 		my %unique_seq_count_samples;
# 		foreach my $sample_name (keys %{$amplicon_sequences->{$marker_name}}) {
# 			my @unique_seqs = keys %{$amplicon_sequences->{$marker_name}{$sample_name}};
# 			# Loops all the unique sequences in the amplicon
# 			foreach my $md5 (@unique_seqs) {
# 				$unique_seq_count_samples{$md5}++;
# 			}
# 		}

		# Annotates sequences that pass all the filters, to do later cross-checking to recover lost real alleles with low frequency or in small clusters
		my %filtered_marker_md5s;

		foreach my $sample_name (@{$samples->{$marker_name}}) {

			# Process only samples with sequences
			if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
# 				print "\t$marker_name-$sample_name doesn't have sequences to filter.\n";
				next;
			}
			print "\t$marker_name-$sample_name filtering\n";

# 			# Exclude samples not in the list
# 			if (defined($paramsdata->{'allowed_samples'}) && !defined($paramsdata->{'allowed_samples'}{'all'}) && (!defined($paramsdata->{'allowed_samples'}{$marker_name}) || !in_array($paramsdata->{'allowed_samples'}{$marker_name},$sample_name))){
# 				next;
# 			}
			# Exclude samples not in the list
			if (defined($filtering_parameters->{'allowed_samples'}) && ref($filtering_parameters->{'allowed_samples'}) eq 'ARRAY' && !in_array($filtering_parameters->{'allowed_samples'},$sample_name)){
				next;
			}

			# Exclude samples with low coverage
			my $total_seqs = $amplicon_depths->{$marker_name}{$sample_name};
			# map $total_seqs+=$amplicon_seq_data->{$marker_name}{$sample_name}{$_}{'depth'}, keys %{$amplicon_seq_data->{$marker_name}{$sample_name}};
			if (defined($filtering_parameters->{'min_amplicon_depth'}) && (!defined($total_seqs) || $total_seqs<$filtering_parameters->{'min_amplicon_depth'})){
				next;
			}

			# Order by coverage
# 			my @sorted_depth_amplicon_unique_md5s = sort { $amplicon_sequences->{$marker_name}{$sample_name}{$b} <=> $amplicon_sequences->{$marker_name}{$sample_name}{$a} } keys %{$amplicon_sequences->{$marker_name}{$sample_name}};
# 			my @sorted_depth_amplicon_unique_seqs = map $marker_seq_data->{$marker_name}{$_}{'seq'}, @sorted_depth_amplicon_unique_md5s;
			my @sorted_depth_amplicon_unique_md5s = sort { $amplicon_seq_data->{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->{$marker_name}{$sample_name}{$a}{'depth'} } keys %{$amplicon_seq_data->{$marker_name}{$sample_name}};
			my @sorted_depth_amplicon_unique_seqs = map $marker_seq_data->{$marker_name}{$_}{'seq'}, @sorted_depth_amplicon_unique_md5s;

			# Aligns sequences to filter by identity
			my ($aligned_seqs,$aligned_md5s,$highest_freq_aligned_seq);
			if (!defined($filtering_parameters->{'degree_of_change'}) && defined($filtering_parameters->{'min_amplicon_seq_identity'}) && $filtering_parameters->{'min_amplicon_seq_identity'}>0 ){
				($aligned_seqs,$aligned_md5s) = multiple_align_seqs(\@sorted_depth_amplicon_unique_seqs, \@sorted_depth_amplicon_unique_md5s, 'mafft --auto');
				$highest_freq_aligned_seq = $aligned_seqs->[0];
			}
			my $highest_amplicon_frequency = $amplicon_seq_data->{$marker_name}{$sample_name}{$sorted_depth_amplicon_unique_md5s[0]}{'freq'};

			# Calculates values for DOC filtering (ROCs = depths)
			my ($doc_allele_number,$DOCns);
			if (defined($filtering_parameters->{'degree_of_change'}) && $#sorted_depth_amplicon_unique_md5s>0){
				my $max_allele_number = 10;
				if (defined($filtering_parameters->{'max_allele_number'})){
					$max_allele_number = $filtering_parameters->{'max_allele_number'};
				}
				my @sorted_depths = map $amplicon_seq_data->{$marker_name}{$sample_name}{$_}{'depth'}, @sorted_depth_amplicon_unique_md5s;
				($doc_allele_number,$DOCns) = degree_of_change(\@sorted_depths,$max_allele_number);
			}

# 			my $doc_allele_number;
# 			if (defined($filtering_parameters->{'degree_of_change'}) && $#sorted_depth_amplicon_unique_md5s>0){
# 				# $filtering_parameters->{'degree_of_change'} stores the number of max. expected alleles
# 				# Calculates DOCs (ROC[$i]/ROC[$i+i] = depth[$i]/depth[$i+i])
# 				my $max_DOCn = 0;
# 				my (@DOCs,@DOCns);
# 				my $i_max = $filtering_parameters->{'degree_of_change'}-1; # Original paper 10-1
# 				if ($#sorted_depth_amplicon_unique_md5s<$i_max) {
# 					$i_max = $#sorted_depth_amplicon_unique_md5s;
# 				}
# 				for (my $i=0; $i<$i_max; $i++){
# 					push(@DOCs, $amplicon_seq_data->{$marker_name}{$sample_name}{$sorted_depth_amplicon_unique_md5s[$i]}{'depth'}/$amplicon_seq_data->{$marker_name}{$sample_name}{$sorted_depth_amplicon_unique_md5s[$i+1]}{'depth'});
# 				}
# 				my $sum_DOCs = sum(@DOCs);
# 				for (my $i=0; $i<=$#DOCs; $i++){
# 					push(@DOCns, $DOCs[$i]/$sum_DOCs*100);
# 				}
# 				# Finds the max DOCn position (max. number of alleles)
# 				$DOCns[$max_DOCn] > $DOCns[$_] or $max_DOCn = $_ for 1 .. $#DOCns;
# 				$doc_allele_number = $max_DOCn+1
# 			}

# IMPROVED VERSION: but not 100% like original Lighten method
# 			if (defined($filtering_parameters->{'degree_of_change'}) && $#sorted_depth_amplicon_unique_md5s>0){
# 				# Calculates DOCs (ROC[$i]/ROC[$i+i] = depth[$i]/depth[$i+i])
# 				my (@DOCns,$doc_allele_number);
# # 				my $i_max = 9;
# # 				if ($#sorted_depth_amplicon_unique_md5s<$i_max) {
# 					my $i_max = $#sorted_depth_amplicon_unique_md5s;
# # 				}
# # 				for (my $i=0; $i<$#sorted_depth_amplicon_unique_md5s; $i++){
# 				for (my $i=0; $i<$i_max; $i++){
# 					push(@DOCs, $amplicon_seq_data->{$marker_name}{$sample_name}{$sorted_depth_amplicon_unique_md5s[$i]}{'depth'}/$amplicon_seq_data->{$marker_name}{$sample_name}{$sorted_depth_amplicon_unique_md5s[$i+1]}{'depth'});
# 				}
# 				my $sum_DOCs = sum(@DOCs);
# 				# Includes an additional DOC of 1 (equivalent to include a copy of the last sequence, it will allow to identify the last sequence as allele)
# 				# ANYWAY IT FAILS WHEN THERE ARE NO ARTIFACTS AND ALLELES HAVE SIMILAR DEPTHS
# 				push(@DOCs, 1);
# 				for (my $i=0; $i<=$#DOCs; $i++){
# 					push(@DOCns, $DOCs[$i]/$sum_DOCs*100);
# 				}
# 				# Finds the min DOCn position
# 				my $min_DOCn = $#DOCns;
# 				$DOCns[$min_DOCn] < $DOCns[$_] or $min_DOCn = $_ for $#DOCns .. 1;
# 				# Finds the max DOCn position (max. number of alleles)
# 				$DOCns[$max_DOCn] > $DOCns[$_] or $max_DOCn = $_ for 1 .. $#DOCns;
# # 				# If the difference between $min_DOCn and $max_DOCn is small, then take all the alleles
# # 				if ($DOCns[$max_DOCn]-$DOCns[$min_DOCn]<5) {
# # 					$max_DOCn = $#DOCns;
# # 				}
# 				# If both $min_DOCn and $max_DOCn have high values, then take all the alleles
# 				# TO RECOVER ALLELES WITH SIMILAR DEPTHS AND WITHOUT ARTIFACTS
# 				if ($DOCns[$max_DOCn]>25 && $DOCns[$min_DOCn]>25) {
# 					$max_DOCn = $#DOCns;
# 				}
# 				$doc_allele_number = $max_DOCn+1
# 			}

			my (@filtered_amplicon_sequences, @filtered_amplicon_md5s);
			my $reading_frame;
			for (my $i=0; $i<=$#sorted_depth_amplicon_unique_md5s; $i++){
				my $md5 = $sorted_depth_amplicon_unique_md5s[$i];
				my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
				my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
				my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
				my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
				my $frequency = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'};
				my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
				my $mean_freq = $marker_seq_data->{$marker_name}{$md5}{'mean_freq'};
				my $min_freq = $marker_seq_data->{$marker_name}{$md5}{'min_freq'};
				my $max_freq = $marker_seq_data->{$marker_name}{$md5}{'max_freq'};
				my ($header,$cluster_size);
				if (defined($amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'cluster_size'})){
					$cluster_size = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'cluster_size'};
					$header = sprintf("hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | cluster_size=%s | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $md5, $len, $depth, $frequency, $count_samples, $cluster_size, $mean_freq, $max_freq, $min_freq);
				} else {
					$header = sprintf("hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq);
				}

# if ($md5 eq '1c725d3cbe7bef46efff738741f7643c'){
# print '';
# }
				# Record which filters are not passed
				my @filter_criteria;
				
# 				# Skip filters in case
# 				my $skip_filters = 0;
# 
# 				# Do not filter sequences present in several samples
# 				if (defined($filtering_parameters->{'min_samples_to_keep'}) && $count_samples>=$filtering_parameters->{'min_samples_to_keep'}){
# 					$skip_filters = 1;
# 					push(@filter_criteria, "samples: $count_samples");
# 				# Exclude sequences present in few samples
# 				} elsif (defined($filtering_parameters->{'min_samples'}) && $count_samples<$filtering_parameters->{'min_samples'}){
# 					push(@filter_criteria, "samples: $count_samples");
# 					# next;
# 				}

				# Stop filtering when there are as many clusters as expected alleles (saves time without altering results)
				if (defined($filtering_parameters->{'max_allele_number'}) && scalar @filtered_amplicon_md5s >= $filtering_parameters->{'max_allele_number'}){

					push(@filter_criteria, sprintf("maximum allele number reached (%d)", $filtering_parameters->{'max_allele_number'}));

				} else {

					# Exclude amplicon sequences after sudden depth change (DOC method)
					if (defined($filtering_parameters->{'degree_of_change'})){

						if ($i>=$doc_allele_number && defined($DOCns->[$i])){
							push(@filter_criteria, sprintf("DOC %.2f<%.2f", $DOCns->[$i], $DOCns->[$doc_allele_number-1]));
							# next;
						} elsif ($i>=$doc_allele_number) {
							push(@filter_criteria, "DOC");
						}

					} else {

						# Exclude small clusters
						if (defined($cluster_size) && defined($filtering_parameters->{'min_cluster_size'}) && $cluster_size<$filtering_parameters->{'min_cluster_size'}){
							push(@filter_criteria, "cluster_size: $cluster_size");
							# next;
						}
# if ($md5 eq '472074ed52cf30d4ce76e49392e0527e'){
# print '';
# }

						# Exclude amplicon sequences with low depth
						if (defined($filtering_parameters->{'min_amplicon_seq_depth'}) && $depth<$filtering_parameters->{'min_amplicon_seq_depth'} ){
							push(@filter_criteria, "depth: $depth");
							# next;
						}
						
						# Exclude amplicon sequences with low frequency
						if (defined($filtering_parameters->{'min_amplicon_seq_frequency'}) && $frequency<$filtering_parameters->{'min_amplicon_seq_frequency'} ){
							push(@filter_criteria, sprintf("frequency: %.2f", $frequency));
							# next;
						}

						# Exclude amplicon sequences with low frequency compared with previous ones
						if (defined($filtering_parameters->{'min_dominant_seq_frequency'}) && $frequency<$highest_amplicon_frequency*$filtering_parameters->{'min_dominant_seq_frequency'}/100 ){
							push(@filter_criteria, sprintf("dominant_freq: %.2f", $frequency));
							# next;
						}

						# Exclude amplicon sequences with erroneous length
						if (defined($filtering_parameters->{'max_amplicon_length_error'})){
							my $right_len = 0;
							foreach my $marker_len (@marker_lengths){
								my $max_len = $marker_len + $filtering_parameters->{'max_amplicon_length_error'};
								my $min_len = $marker_len - $filtering_parameters->{'max_amplicon_length_error'};
								if ($len>=$min_len && $len<=$max_len){
									$right_len = 1;
									last;
								}
							}
							if (!$right_len) {
								push(@filter_criteria, "length: $len");
								# next;
							}
						}

						# Exclude amplicon sequences with erroneous length
						if (defined($filtering_parameters->{'discard_frameshifts'})){
							my $right_len = 0;
							foreach my $marker_len (@marker_lengths){
								if (abs($marker_len-$len) % 3 == 0){
									$right_len = 1;
									last;
								}
							}
							if (!$right_len) {
								push(@filter_criteria, "frameshift");
								# next;
							}
						}

						# Exclude amplicon sequences with low identity
						if (defined($filtering_parameters->{'min_amplicon_seq_identity'}) && $filtering_parameters->{'min_amplicon_seq_identity'}>0 &&
						defined($highest_freq_aligned_seq) && defined($aligned_seqs->[$i])){
							my ($identity,$total) = binary_score_nts($highest_freq_aligned_seq,$aligned_seqs->[$i]);
							if ($identity < $filtering_parameters->{'min_amplicon_seq_identity'}){
								push(@filter_criteria, sprintf("identity: %.2f", $identity));
								# next;
							}
						}

						# Exclude amplicon non coding sequences (with stop codons in reading frame)
						if (defined($filtering_parameters->{'discard_noncoding'})){
							my $coding = 0;
							# Sets the reading frame with the first dominant sequence
							if (!defined($reading_frame)){
								for (my $frame=0; $frame<=2; $frame++){
									my $prot_seq = dna_to_prot(substr($seq,$frame));
									if ($prot_seq !~ /\*/) {
										$coding = 1;
										$reading_frame = $frame;
										last;
									}
								}
							} else {
								my $prot_seq = dna_to_prot(substr($seq,$reading_frame));
								if ($prot_seq !~ /\*/) {
									$coding = 1;
								}
							}
							if (!$coding) {
								push(@filter_criteria, "noncoding");
								# next;
							}
						}

						# Exclude chimeras
						if (defined($filtering_parameters->{'min_chimera_length'}) && $#filtered_amplicon_sequences>0){
							my ($is_chimera,$i_seq1,$i_seq2) = is_chimera($seq,\@filtered_amplicon_sequences,$filtering_parameters->{'min_chimera_length'},$max_chimera_ident);
							if ($is_chimera){
								push(@filter_criteria, sprintf("chimera: '%s'+'%s'",$filtered_amplicon_md5s[$i_seq1],$filtered_amplicon_md5s[$i_seq2]));
								# print "chimera: '$md5'='$seq1_md5'+'$seq2_md5'\n";
							}
						}

					}

				}

				if (!@filter_criteria) {
					$filtered_marker_md5s{$md5}++;
					push(@filtered_amplicon_sequences, $seq);
					push(@filtered_amplicon_md5s, $md5);
					$md5_to_sequence->{$md5} = $seq;
					$filtered_amplicon_sequences->{$marker_name}{$sample_name}{$md5} = $depth;
					$filtered_amplicon_depths->{$marker_name}{$sample_name} += $depth;
					$filters_output->{$marker_name}{$sample_name} .= sprintf(">*%s | %s\n%s\n", $name, $header, $seq);
				} else {
					$filters_output->{$marker_name}{$sample_name} .= sprintf(">#%s | Filtered: %s | %s\n%s\n", $name, join(', ',@filter_criteria), $header, $seq);
				}

			}
			if (defined($filtered_amplicon_depths->{$marker_name}{$sample_name})){
				if ($sample_name ne ''){
					printf("\t%s-%s filtered (%d reads, %d variants)\n", $marker_name, $sample_name, $filtered_amplicon_depths->{$marker_name}{$sample_name}, scalar keys %{$filtered_amplicon_sequences->{$marker_name}{$sample_name}});
				} else {
					printf("\t%s filtered (%d reads, %d variants)\n", $marker_name, $filtered_amplicon_depths->{$marker_name}{$sample_name}, scalar keys %{$filtered_amplicon_sequences->{$marker_name}{$sample_name}});
				}
			} else {
				if ($sample_name ne ''){
					printf("\t%s-%s filtered (0 sequences)\n", $marker_name, $sample_name);
				} else {
					printf("\t%s filtered (0 sequences)\n", $marker_name);
				}
			}

		}

# 		# Clustering cross-checking to recover lost real alleles with low frequency or in small clusters
# 		if (defined($filtering_parameters->{'min_samples_to_keep'})){
# 			foreach my $sample_name (keys %{$filtered_amplicon_sequences->{$marker_name}}){
# 				foreach my $md5 (keys %{$amplicon_seq_data->{$marker_name}{$sample_name}}){
# if ($md5 eq '73ce3ddd162a2f04b061dfe658357ce6'){
# print '';
# }
# 					# Annotate sequence if it is present as clear allele in other samples
# 					if (!defined($filtered_amplicon_sequences->{$marker_name}{$sample_name}{$md5}) && defined($filtered_marker_md5s{$md5}) && $filtered_marker_md5s{$md5} >= $filtering_parameters->{'min_samples_to_keep'}){
# 						my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
# 						my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
# 						my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
# 						my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
# 						my $frequency = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'};
# 						my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
# 						my $mean_freq = $marker_seq_data->{$marker_name}{$md5}{'mean_freq'};
# 						my $min_freq = $marker_seq_data->{$marker_name}{$md5}{'min_freq'};
# 						my $max_freq = $marker_seq_data->{$marker_name}{$md5}{'max_freq'};
# 						my $cluster_size = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'cluster_size'};
# 						my $header = sprintf("%s | hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $name, $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq);
# 
# 						$filtered_amplicon_sequences->{$marker_name}{$sample_name}{$md5} = $depth;
# 						$filtered_amplicon_depths->{$marker_name}{$sample_name} += $depth;
# 						$filters_output->{$marker_name}{$sample_name} .= sprintf(">%s\n%s\n%s\n\n", $header, "samples: $filtered_marker_md5s{$md5}", $seq);
# 					}
# 				}
# 			}
# 		}
	}

	return ($filtered_amplicon_sequences,$filtered_amplicon_depths,$filters_output,$md5_to_sequence);

}

#################################################################################

# # Filters alleles according to thresholds
sub filter_amplicon_sequences_with_threads {

	my ($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$amplicon_depths,$threads_limit) = @_;

	if (!defined($threads_limit)){
		$threads_limit = 4;
	}

	my ($filtered_amplicon_sequences, $filtered_amplicon_depths, $filters_output);
	my $md5_to_sequence = {};
	my $one_amplicon_seq_data;
	foreach my $marker_name (@$markers){
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}
		my $filtered_amplicon_depths->{$marker_name} = {};
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
				next;
			}
			my $one_amplicon_seq_data_;
			my $filtered_amplicon_sequences->{$marker_name}{$sample_name} = {};
			$one_amplicon_seq_data_->{$marker_name}{$sample_name} = $amplicon_seq_data->{$marker_name}{$sample_name};
			push(@{$one_amplicon_seq_data},$one_amplicon_seq_data_);
			# print '';
		}
	}

	my @threads;
	for (my $count_amplicon=0; $count_amplicon<=$#{$one_amplicon_seq_data}; $count_amplicon++){

		push(@threads, threads->create(\&filter_amplicon_sequences,$markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$one_amplicon_seq_data->[$count_amplicon],$amplicon_depths));
# 		print "\n";

# 		# For debugging:
# 		push(@threads, [filter_amplicon_sequences($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$one_amplicon_seq_data->[$count_amplicon],$amplicon_depths)]);

		# If maximum number of threads is reached or last sbjct of a query is processed
		if (scalar @threads >= $threads_limit  || $count_amplicon == $#{$one_amplicon_seq_data}){
			my $check_threads = 1;
			while ($check_threads){
				for (my $i=0; $i<=$#threads; $i++){
					unless ($threads[$i]->is_running()){
						my ($filtered_amplicon_sequences_,$filtered_amplicon_depths_,$filters_output_,$md5_to_sequence_) = $threads[$i]->join;
						if (defined($filtered_amplicon_sequences_)){
							my $marker_name = (keys %$filtered_amplicon_sequences_)[0];
							my $sample_name = (keys %{$filtered_amplicon_sequences_->{$marker_name}})[0];
	# 						print "\t$marker_name-$sample_name finished\n";
							if (defined($filtered_amplicon_sequences_->{$marker_name}{$sample_name})){
								$filtered_amplicon_sequences->{$marker_name}{$sample_name} = $filtered_amplicon_sequences_->{$marker_name}{$sample_name};
								$filtered_amplicon_depths->{$marker_name}{$sample_name} = $filtered_amplicon_depths_->{$marker_name}{$sample_name};
								$filters_output->{$marker_name}{$sample_name} .= $filters_output_->{$marker_name}{$sample_name};
								$md5_to_sequence = { %$md5_to_sequence, %$md5_to_sequence_};
							}
						}
						undef($threads[$i]);
						splice(@threads,$i,1);
						$i = $i - 1;
						unless ($count_amplicon == $#{$one_amplicon_seq_data} && @threads){
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
# 				print '';
# 				my ($filtered_amplicon_sequences_,$filtered_amplicon_depths_,$filters_output_,$md5_to_sequence_) = @{$threads[$i]};
# 				my $marker_name = (keys %$filtered_amplicon_sequences_)[0];
# 				my $sample_name = (keys %{$filtered_amplicon_sequences_->{$marker_name}})[0];
# 				$filtered_amplicon_sequences->{$marker_name}{$sample_name} = $filtered_amplicon_sequences_->{$marker_name}{$sample_name};
# 				$filtered_amplicon_depths->{$marker_name}{$sample_name} = $filtered_amplicon_depths_->{$marker_name}{$sample_name};
# 				$filters_output->{$marker_name}{$sample_name} .= $filters_output_->{$marker_name}{$sample_name};
# 				$md5_to_sequence = { %$md5_to_sequence, %$md5_to_sequence_};
# 				delete $threads[$i];
# 			}
		}
	}
	
	return ($filtered_amplicon_sequences,$filtered_amplicon_depths,$filters_output,$md5_to_sequence);

}

#################################################################################

# Clusters unique sequences to extract real sequences/alleles
sub cluster_amplicon_sequences {

	my ($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$referencedata) = @_;

	# Retrieves reference sequences to improve clustering results
	my (%reference_seqs,$reference_seqs,$reference_names);
	if (defined($referencedata)){
		%reference_seqs = map { $_ => $referencedata->{$_}{'sequence'} } keys %$referencedata;
		($reference_seqs, $reference_names) = sequences_hash_to_array(\%reference_seqs);
	}
	
	# REMOVED, TOO MANY ERRONEOUS ALLELES WITH FAST MULTIPLE ALIGNMENT
	# If fast multiple alignment is specified, the full amplicon unique sequences will be multiple aligned only once before clustering (fast)
	# If accurate alignment is desired (default), the sequences will be globally aligned one by one, not the full amplicon at once
	# Accurate alignment is 30 times slower
# 	if (!defined($fast_multiple_align)) {
# 		$fast_multiple_align = 0;
# 	}

	my ($amplicon_clustered_sequences, $amplicon_clustered_depths, $clusters, $md5_to_sequence);
	foreach my $marker_name (@$markers){

		# Skips markers without data
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}

		# Extracts clustering parameters
		my %clustering_thresholds;
		foreach my $clustering_threshold (keys %$paramsdata) {
			if (defined($paramsdata->{$clustering_threshold}{$marker_name})){
				$clustering_thresholds{$clustering_threshold} = $paramsdata->{$clustering_threshold}{$marker_name}[0];
			} elsif (defined($paramsdata->{$clustering_threshold}{'all'})){
				$clustering_thresholds{$clustering_threshold} = $paramsdata->{$clustering_threshold}{'all'}[0];
			}
		}

		# Loops samples/amplicons and performs clustering rounds
		foreach my $sample_name (@{$samples->{$marker_name}}) {

			# Process only samples with sequences (after filtering) in the original order
			if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
# 				print "\t$marker_name-$sample_name doesn't have sequences to cluster.\n";
				next;
			}
			if (defined($clustering_thresholds{'cluster_exact_length'}) || defined($clustering_thresholds{'cluster_inframe'})){
				if (!defined($markerdata->{$marker_name}{'length'})){
					print "\nERROR: Marker '$marker_name' has no correct length/s specified, please include 'length' value/s into CSV data file or run the analysis in 'Automatic' mode.\n\n";
					exit;
				}
			} 
			if ($sample_name ne ''){
				printf("\t%s-%s clustering\n",$marker_name,$sample_name);
			} else {
				printf("\t%s clustering\n",$marker_name);
			}

			# Reads amplicon unique sequences, full headers and depths and stores them into hashes with md5s as keys
			my (%amplicon_sequences, %amplicon_headers, %amplicon_names, %amplicon_depths, %amplicon_frequencies);
			foreach my $md5 (keys %{$amplicon_seq_data->{$marker_name}{$sample_name}}) {
				$amplicon_sequences{$md5} = $marker_seq_data->{$marker_name}{$md5}{'seq'};
				my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
				my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
				my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
				my $frequency = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'};
				my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
				my $mean_freq = $marker_seq_data->{$marker_name}{$md5}{'mean_freq'};
				my $min_freq = $marker_seq_data->{$marker_name}{$md5}{'min_freq'};
				my $max_freq = $marker_seq_data->{$marker_name}{$md5}{'max_freq'};
				$amplicon_headers{$md5} = sprintf("%s | hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $name, $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq);
				$amplicon_names{$md5} = $name;
				$amplicon_depths{$md5} = $depth;
				$amplicon_frequencies{$md5} = $frequency;
			}

# 			print "...";

			# Sorts amplicon sequences by number of reads
			# To have sequences ordered by depth helps in the clustering (major seqs will be the firsts of cluster)
			my $reference_seq;
			my (@unique_seq_seqs, @unique_seq_md5s);
			my @sorted_depth_amplicon_unique_seqs = sort { $amplicon_depths{$b} <=> $amplicon_depths{$a} } keys %amplicon_depths;
			my $count_unique_seqs = 0;
			foreach my $md5 (@sorted_depth_amplicon_unique_seqs) {

				push(@unique_seq_seqs, $amplicon_sequences{$md5});
				push(@unique_seq_md5s, $md5);
				# push(@{$marker_names}, sprintf("%s-%d", $md5, $amplicon_depths{$md5}));

				# CLUSTERING MORE THAN 500-1000 UNIQUE SEQUENCES IS EXTREMELY SLOW
				$count_unique_seqs++;
# 				if ($count_unique_seqs >= 1000) {
# 					last;
# 				}

				# If reference sequences are provided, finds the best reference for the amplicon
				# Align the unique sequences ordered by depth till one matches a reference with good Evalue
				# To use this reference sequence in 'is_dominant()' checking
				if (defined($referencedata) && !defined($reference_seq)){
					# Finds the most similar reference sequence for clustering
					my $ref_align_data = align_seqs([$md5],[$amplicon_sequences{$md5}],$reference_names,$reference_seqs,0,'dna blastn -evalue 0.001');
					if (defined($ref_align_data->{$md5})){
						my $reference_name = $ref_align_data->{$md5}[0]{'NAME'};
						$reference_seq = $reference_seqs{$reference_name};
# 						$reference_seq = (split("\n",$ref_align_data->{$md5}[0]{'ALIGN'}))[1];
					}
				}
			}
			# Doesn't use any reference if there are not suitable ones
			if (defined($referencedata) && !defined($reference_seq)){
				print "\tNo reference sequences match the amplicon '$marker_name-$sample_name'.\n";
				# Uses the first sequence if there are not similar ones
				# my $reference_name = $reference_names->[0];
				# $reference_seq = $reference_seqs->[0];
			}
			# Replaced by $dominant_depth
# 			my $highest_amplicon_frequency = $amplicon_depths{$unique_seq_md5s[0]};

			# Performs multiple alignment of all unique sequences into an amplicon (to use aligned sequences later, FAST METHOD)
			# Only when there are not indels (or very few) in the sequences (ex. Illumina reads), if not the alignment of hundreds of sequences will contain regions not correctly aligned
# 			my $unique_aligned_seqs;
# 			if ($clustering_thresholds{'indel_threshold'}/100*$amplicon_seq_data->{$marker_name}{$sample_name}{$unique_seq_md5s[0]}{'len'} < 1){
# 			if ($fast_multiple_align){
# 				$unique_aligned_seqs = sequences_array_to_hash(multiple_align_seqs(\@unique_seq_seqs, \@unique_seq_md5s, 'mafft --auto')); # 68 seg
# 			}

			# Cluster sequences
			my $clusters_amplicon; # Variable to store clusters
			my $cluster_count = 0; # Variable to number clusters
			my %variants_assignations; # Annotates md5s of sequences annotated as variants of others (included in clusters)
			my @variants_seqs = @unique_seq_seqs;
			my @variants_md5s = @unique_seq_md5s;
			# Stores variants that are equally similar to 2 or more high freq dominant sequences (to not to annotate their depth)
			my $high_freq_seq_variants;
			# Saves 'is_dominant' results to not to repeat the checking
			my %is_dominant;
			# Saves the first dominant in the sample to use it like a reference to check new dominant sequences
			# my ($first_dominant_md5, $first_dominant_seq);
# my @run_times;
			for (my $i=0; $i<=$#unique_seq_md5s; $i++) {
# if ($i+1 % 100 == 0) {
# 	printf("%d ", $i+1);
# }
				my $dominant_md5 = $unique_seq_md5s[$i];
# if ($dominant_md5 eq 'dfcf28d0734569a6a693bc8194de62bf'){
# print '';
# }
				# Checks if is already clustered
				if (!in_array(\@variants_md5s, $dominant_md5)){
					next;
				}
				# Checks if is not in the list of variants assigned to several clusters
				if (defined($high_freq_seq_variants) && defined($high_freq_seq_variants->{$dominant_md5})){
					delete($high_freq_seq_variants->{$dominant_md5});
					next;
				}
				
				my $dominant_seq = $unique_seq_seqs[$i];
				my $dominant_len = length($dominant_seq); # The length before alignment
				my $dominant_depth = $amplicon_depths{$dominant_md5};

				# Do not cluster remaining singletons (saves time without altering results)
				if ($dominant_depth == 1){ 
					if (!(defined($clustering_thresholds{'keep_singletons'}) && $clustering_thresholds{'keep_singletons'}) ){
						last;
					}
				}
				
				# Stop clustering when there are as many clusters as expected alleles (saves time without altering results)
				if (defined($clustering_thresholds{'max_allele_number'}) && $cluster_count >= $clustering_thresholds{'max_allele_number'}){
					last;
				}

				# Checks if the sequence pass the clutering thresholds
				if (!defined($is_dominant{$dominant_md5})){
					if (!is_dominant($amplicon_sequences{$dominant_md5},$markerdata->{$marker_name}{'length'},\%clustering_thresholds,$reference_seq)){
						$is_dominant{$dominant_md5} = 0;
					} else {
						$is_dominant{$dominant_md5} = 1;
					}
				}
				# Doesn't take non-dominant sequences as cluster references
				if (!$is_dominant{$dominant_md5}){
					next;
				}
				# 2 alleles can be very different, use as reference the dominant seq from the same cluster, not from the amplicon
# 				if (!defined($first_dominant_md5)){
# 					$first_dominant_md5 = $dominant_md5;
# 					$first_dominant_seq = $dominant_seq;
# 				}#time=0%

				# Variables to store data from cluster sequences
				my $cluster_data; #@cluster_md5s, @cluster_aligned_seqs, @cluster_depths, @cluster_identities, @cluster_errors);
# 				push(@{$cluster_data},{'md5' => $dominant_md5, 'aligned_seq' => $dominant_seq, 'depth' => $amplicon_depths{$dominant_md5}, 'identity' => 100, 'errors' => ''});
				push(@{$cluster_data},{'md5' => $dominant_md5, 'seq' => $dominant_seq, 'depth' => $amplicon_depths{$dominant_md5}, 'errors' => ''}); # 'identity' => 100,

				# Aligns the reference sequence with all the other sequences 1 by 1 (SLOWER than $fast_multiple_align option), the alignments will be accurate
				# If the multiple alignment is not previously calculated
# 				if (!$fast_multiple_align) {
				my $fast_align = 0;
				my $ref_var_aligned_seqs;
				if ((defined($clustering_thresholds{'indel_threshold'}) && $clustering_thresholds{'indel_threshold'}/100*$dominant_len < 1)
				|| (defined($clustering_thresholds{'fast_clustering_alignment'}) && $clustering_thresholds{'fast_clustering_alignment'})){
					$fast_align = 1;
				} else {
					$ref_var_aligned_seqs = align_seqs2one($dominant_seq,$dominant_md5,\@variants_seqs,\@variants_md5s,'needleall'); #time=20% # 104 seg
				}
# # 				my $ref_var_aligned_seqs = align_seqs2one($dominant_seq,$dominant_md5,\@variants_seqs,\@variants_md5s,'needleman-wunsch'); # 152 seg
# 					# SLOWER:
# # 					for (my $j=0; $j<=$#variants_md5s; $j++) {
# # 						$ref_var_aligned_seqs->{$variants_md5s[$j]} = [align_2seqs($dominant_seq,$variants_seqs[$j], 'needleman-wunsch')]; # 263 seg
# # 						$ref_var_aligned_seqs->{$variants_md5s[$j]} = [align_2seqs($dominant_seq,$variants_seqs[$j], 'fogsaa')]; # 355 seg
# # 					}
# 				}
			
				# Stores high frequency dominant sequences, to compare all variants also against them
				my (@high_freq_seqs, @high_freq_md5s);

				# Looks for variants
				for (my $j=0; $j<=$#variants_md5s; $j++) {

					my $variant_seq = $variants_seqs[$j];
					my $variant_md5 = $variants_md5s[$j];
					my $variant_len = length($variant_seq); # The length before alignment

# if ($dominant_md5 eq '6d68254ddf14fcec6c0cbd83448e00f4' && $variant_md5 eq '68a3a2777f6f1d25c248064bcaeadba7'){
# print '';
# }
					# If variant is the same that reference or if both sequences have been already compared and marked as dissimilar, skip
					if ($variant_md5 eq $dominant_md5){
# 						|| defined($dissimilar_pairs{"$dominant_md5-$variant_md5"}) || defined($dissimilar_pairs{"$variant_md5-$dominant_md5"})){
						next;
					}

					my ($identical, $total, $substitutions, $insertions, $deletions, $homopolymer_indels);
					if (!$fast_align){
						# Skip if the pairwise global alignment fails
						if (!defined($ref_var_aligned_seqs->{$variant_md5}[0]) || !defined($ref_var_aligned_seqs->{$variant_md5}[1])){
							#printf("ERROR:\n>%s\n%s\n>%s\n%s\n",$dominant_md5,$ref_var_aligned_seqs->{$variant_md5}[0],$variant_md5,$ref_var_aligned_seqs->{$variant_md5}[1]);
							next;
						}
	# 					my ($identical,$total);
	# 					if ($fast_multiple_align) {
	# 						($identical,$total) = binary_score_nts($unique_aligned_seqs->{$dominant_md5},$unique_aligned_seqs->{$variant_md5});
	# 					} else {
						# Skips seqs with lower identity than threshold
						if (defined($clustering_thresholds{'identity_threshold'})){
							($identical,$total) = binary_score_nts($ref_var_aligned_seqs->{$variant_md5}[0],$ref_var_aligned_seqs->{$variant_md5}[1]); #time=10%
							#my $identity = sprintf("%.2f", $identical/$dominant_len*100);
							if ($identical/$dominant_len*100 < $clustering_thresholds{'identity_threshold'}){
	# 							$dissimilar_pairs{"$variant_md5-$dominant_md5"} = 1;
								next;
							}
						}

	# 					my ($substitutions, $insertions, $deletions, $homopolymer_indels);
	# 					if ($fast_multiple_align) {
	# 						($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($unique_aligned_seqs->{$dominant_md5},$unique_aligned_seqs->{$variant_md5});
	# 					} else {
						($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($ref_var_aligned_seqs->{$variant_md5}[0],$ref_var_aligned_seqs->{$variant_md5}[1]); #time=7%

						# Skips seqs with more substitutions than threshold and process them in next clustering round
						# Ej. 3 subs > (2%*130=2.6 || 3)
						# Ej. 2 subst > (1%*130=1.3 || 1)
	# 					if (defined($clustering_thresholds{'substitution_threshold'}) && scalar @$substitutions > ($clustering_thresholds{'substitution_threshold'}/100*$dominant_len || sprintf("%.0f", $clustering_thresholds{'substitution_threshold'}/100*$dominant_len)) ){
# if (!defined($substitutions)){
# print '';
# }
						if (defined($clustering_thresholds{'substitution_threshold'}) && scalar @$substitutions > sprintf("%.0f", $clustering_thresholds{'substitution_threshold'}/100*$dominant_len) ){
	# 						$dissimilar_pairs{"$variant_md5-$dominant_md5"} = 1;
							next;
						}

						# Skips seqs with more non homopolymer indels than threshold and process them in next clustering round
						my $non_homopolymer_indels = scalar @$insertions + scalar @$deletions - scalar @$homopolymer_indels;
	# 					if (defined($clustering_thresholds{'indel_threshold'}) && $non_homopolymer_indels > ($clustering_thresholds{'indel_threshold'}/100*$dominant_len || sprintf("%.0f", $clustering_thresholds{'indel_threshold'}/100*$dominant_len)) ){
						if (defined($clustering_thresholds{'indel_threshold'}) && $non_homopolymer_indels > sprintf("%.0f", $clustering_thresholds{'indel_threshold'}/100*$dominant_len) ){
	# 						$dissimilar_pairs{"$variant_md5-$dominant_md5"} = 1;
							next;
						}
					} else { # If there are not indels, perform fast alignment, 4 times faster (5 times faster with C code)
						if ($dominant_len != $variant_len) {
							next;
						}
						my ($max_errors,$errors);
						if (defined($clustering_thresholds{'identity_threshold'}) && (!defined($clustering_thresholds{'substitution_threshold'})
						|| $clustering_thresholds{'identity_threshold'} > $clustering_thresholds{'substitution_threshold'})){
							$max_errors = sprintf("%.0f", (1-$clustering_thresholds{'identity_threshold'}/100)*$dominant_len);
						} elsif (defined($clustering_thresholds{'substitution_threshold'})){
							$max_errors = sprintf("%.0f", $clustering_thresholds{'substitution_threshold'}/100*$dominant_len);
						} 
						$errors = count_errors($dominant_seq,$variant_seq,$max_errors+1); # Perl code
# 						$errors = count_errors_c($dominant_seq,$variant_seq,$max_errors); # 15% faster C code
						if ($errors > $max_errors) {
							next;
						}
						if ($dominant_len > $variant_len){
							$total = $dominant_len;
						} else {
							$total = $variant_len;
						}
						$identical = $total - $errors; 
						($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($dominant_seq,$variant_seq); #time=7%
						$ref_var_aligned_seqs->{$variant_md5} = [$dominant_seq, $variant_seq];
					}
# my $start_run = gettimeofday;
					# Check this after checking the rest of similarity conditions, if not all the variants will be checked every time
					# Skips high depth/frequency sequences (before clustering) compared or not to the dominant frequency and process them in next clustering round
					# Only when reference freq is higher
					if ( ( defined($clustering_thresholds{'min_amplicon_seq_frequency_threshold'}) && $amplicon_frequencies{$variant_md5} >= $clustering_thresholds{'min_amplicon_seq_frequency_threshold'} )
					|| ( defined($clustering_thresholds{'min_dominant_frequency_threshold'}) && 100*$amplicon_depths{$variant_md5}/$dominant_depth >= $clustering_thresholds{'min_dominant_frequency_threshold'} ) ){
						# Checks if the variant passes the dominant thresholds
						# The dominant seq is given as reference to avoid false positives with 3 insertions or deletions not consecutive
						# The dominant has already been checked  and compared against references (if they exist), so is a good reference
						if (!defined($is_dominant{$variant_md5})){
							if (!defined($reference_seq)) {
								if (is_dominant($amplicon_sequences{$variant_md5},$markerdata->{$marker_name}{'length'},\%clustering_thresholds,$ref_var_aligned_seqs->{$variant_md5})){
									$is_dominant{$variant_md5} = 1;
								} else {
									$is_dominant{$variant_md5} = 0;
								}
							# Check if the variant passes the length thresholds compared with reference (if reference is given as input)
							# my $variant_len = length($amplicon_sequences{$variant_md5}); # The length before alignment
							} elsif (!defined($is_dominant{$variant_md5})) { 
								if (is_dominant($amplicon_sequences{$variant_md5},$markerdata->{$marker_name}{'length'},\%clustering_thresholds,$reference_seq)){
									$is_dominant{$variant_md5} = 1;
								} else {
									$is_dominant{$variant_md5} = 0;
								}
							}
						}
						if ($is_dominant{$variant_md5}) {
							push(@high_freq_seqs, $variant_seq);
							push(@high_freq_md5s, $variant_md5);
							# Go to next variant
							next;
						}
					} #time=40%
# push(@run_times, gettimeofday-$start_run);

					# Check all the high frequency sequences and skip variant if it is more similar to any of them
					if (@high_freq_md5s){
						if (!defined($identical)){
							($identical,$total) = binary_score_nts($ref_var_aligned_seqs->{$variant_md5}[0],$ref_var_aligned_seqs->{$variant_md5}[1]);
						}
						my $skip_variant = 0;
						for (my $k=0; $k<=$#high_freq_md5s; $k++){
							my $high_freq_seq = $high_freq_seqs[$k];
							my $high_freq_md5 = $high_freq_md5s[$k];
# 							my ($identical_,$total_);
# 							if ($fast_multiple_align) {
# 								($identical_,$total_) = binary_score_nts($unique_aligned_seqs->{$high_freq_md5}, $unique_aligned_seqs->{$variant_md5});
# 							} else {
							my ($identical_,$total_) = binary_score_nts(align_2seqs($high_freq_seq,$variant_seq, 'needle'));
							# If alignment fails
							if (!defined($identical_) || !defined($total_)) {
								next;
							# If identity is higher
							} elsif ($identical_>$identical){
								$skip_variant = 1;
								last;
							# If identity is the same, we mark the variant, we don't know to which one to assign the depth, so it will not be assigned to any of them
							} elsif ($identical_==$identical) { 
								$high_freq_seq_variants->{$variant_md5}{$dominant_md5} = 1;
								$high_freq_seq_variants->{$variant_md5}{$high_freq_md5} = 1;
								# last;
# 								splice(@variants_md5s,$j,1);
# 								splice(@variants_seqs,$j,1);
# 								$j--;
# 								last;
							}
						}
						if ($skip_variant) {
							next;
						}
					}#time=0%

					# Adds depth and depth to reference sequence data (real allele, or chimera)
# 					$amplicon_clustered_sequences->{$marker_name}{$sample_name}{$dominant_md5} += $amplicon_depths{$variant_md5};
# 					$amplicon_clustered_depths->{$marker_name}{$sample_name} += $amplicon_depths{$variant_md5};
					# Prints artifact sequence and errors
					my $variant_errors = sprintf("sub: %d, ins: %d, del: %d, homo_indel: %d", scalar @$substitutions, scalar @$insertions, scalar @$deletions, scalar @$homopolymer_indels);
# 					$variant_errors = print_sequence_errors($substitutions, $insertions, $deletions);
					my $depth;
					# Annotates the depth only if the variant is only assignated to this reference sequence
					if (!defined($variants_assignations{$variant_md5}) && (!defined($high_freq_seq_variants) || !defined($high_freq_seq_variants->{$variant_md5}))){
						$depth = $amplicon_depths{$variant_md5};
					# If the variant is assignated to several reference sequences, depth is annotated between parenthesis
					} else {
						$depth = '('.$amplicon_depths{$variant_md5}.')';
					}
# 					} elsif (defined($high_freq_seq_variants) && defined($high_freq_seq_variants->{$variant_md5})){
# 						$depth = '('.$amplicon_depths{$variant_md5}.')';
# 					} else {
# 						$depth = $amplicon_depths{$variant_md5};
# 					}
# 					# Removes the dominant sequence from the list of high freq seqs assigned to a variant
# 					# If there are not sequences assigned (empty hash), the variant will be removed
# 					if (defined($high_freq_seq_variants) && defined($high_freq_seq_variants->{$variant_md5}) && defined($high_freq_seq_variants->{$variant_md5}{$dominant_md5})){
# 						delete($high_freq_seq_variants->{$variant_md5}{$dominant_md5});
# 						if (!%{$high_freq_seq_variants->{$variant_md5}}){
# 							delete($high_freq_seq_variants->{$variant_md5});
# 						}
# 					}
					# Annotates cluster member data
					push(@{$cluster_data},{'md5' => $variant_md5, 'seq' => $variant_seq, 'depth' => $depth, 'errors' => $variant_errors}); #  'identity' => $identity,

				}#time=52%

				# Creates a new cluster
				$cluster_count++;
				my $cluster_name = $cluster_count;

				# Checks one by one all sequences in cluster and creates consensus sequence
				my ($consensus_seq, $consensus_md5);
				if ($#{$cluster_data} > 1){
					my @cluster_seqs = map $_->{'seq'} , @{$cluster_data};
					my @cluster_md5s = map $_->{'md5'} , @{$cluster_data};
					my ($cluster_aligned_seqs,$cluster_aligned_md5s) = multiple_align_seqs(\@cluster_seqs, \@cluster_md5s, 'mafft'); #time=9% # 'mafft --auto' #time=17%
					my $position_depth_data;
					for (my $j=0; $j<=$#{$cluster_data}; $j++){
						my @seq = split('', $cluster_aligned_seqs->[$j]);
						my $depth = $cluster_data->[$j]{'depth'};
						for (my $k=0; $k<=$#seq; $k++) {
							if (is_numeric($depth)){
								$position_depth_data->[$k]{$seq[$k]} += $depth;
							}
						}
					}
					my $consensus_aligned_seq;
					foreach my $position_depths (@{$position_depth_data}){
						my $max_nt;
						my $max_depth = -1;
						while ((my $nt, my $depth) = each %$position_depths) {
							if ($depth > $max_depth) {
								$max_depth = $depth;
								$max_nt = $nt;
							}
						}
						$consensus_aligned_seq .= $max_nt;
					}
					$consensus_seq = uc($consensus_aligned_seq);
					$consensus_seq =~ s/-//g;
					$consensus_md5 = generate_md5($consensus_seq);
# if ($consensus_md5 eq '58b33b5aeafa325afd426405ceaad436'){
# print '';
# }
					# If the consensus is not the reference sequence
					# and the consensus sequence passes dominant thresholds
					# then it's added as the reference of the cluster
					if ($consensus_md5 ne $dominant_md5) {
						if (!defined($variants_assignations{$consensus_md5})) {
							# Check if the consensus seq passes the dominant thresholds using the dominant seq as reference
							# The dominant seq is given as reference to avoid false positives with 3 insertions or deletions not consecutive
							# The dominant has already been checked  and compared against references (if they exist), so is a good reference
							#my $consensus_len = length($consensus_seq_); # The length before alignment
							if (!defined($is_dominant{$consensus_md5})){
								if (!defined($reference_seq)){
									if (is_dominant($consensus_seq,$markerdata->{$marker_name}{'length'},\%clustering_thresholds,$consensus_seq)){
										$is_dominant{$consensus_md5} = 1;
									} else {
										$is_dominant{$consensus_md5} = 0;
									}
								# Check if the conesnsus passes the length thresholds compared with reference (if reference is given as input)
								} elsif (!defined($is_dominant{$consensus_md5})) { 
									if (is_dominant($consensus_seq,$markerdata->{$marker_name}{'length'},\%clustering_thresholds,$reference_seq)){
										$is_dominant{$consensus_md5} = 1;
									} else {
										$is_dominant{$consensus_md5} = 0;
									}
								}
							}
							# If the consensus seq passes the dominant thresholds and it has not been previously clustered, then it will be the cluster ref
							if ($is_dominant{$consensus_md5}) {
								$dominant_md5 = $consensus_md5;
								$dominant_seq = $consensus_seq;
								# If the consensus sequence is in the same cluster, move to the first position
								my ($in_same_cluster,$pos) = in_array(\@cluster_md5s,$consensus_md5,1);
								if ($in_same_cluster) {
									unshift(@{$cluster_data},$cluster_data->[$pos->[0]]);
									splice(@{$cluster_data},$pos->[0]+1,1);
								# If the consensus sequence doesn't exist within the cluster
								} else {
									unshift(@{$cluster_data},{'md5' => $consensus_md5, 'seq' => $consensus_seq, 'depth' => 0, 'errors' => 'CONSENSUS'}); # 'identity' => '100',
								}
								# It could happen for a consensus sequence that it doesn't exist in the original data
								if (!defined($amplicon_depths{$consensus_md5})){
									# Add the sequence to sequence data variables
									#$md5_to_sequence->{$consensus_md5} = $consensus_seq;
									$amplicon_sequences{$consensus_md5} = $consensus_seq;
									$amplicon_depths{$consensus_md5} = 0;
									$amplicon_headers{$consensus_md5} = sprintf("%s | hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", 'CONSENSUS', $consensus_md5, length($consensus_seq), 0, 0, 0, 0, 0, 0);
									# printf("New sequence: %s\n", $consensus_md5);
								}
							}
						# If the consensus sequence has been already clustered
						} else {
							my $previous_cluster = $variants_assignations{$consensus_md5};
							# If the consensus is already a reference sequence from a cluster, then add the sequences to the previous cluster
							$cluster_name = $previous_cluster;
							$cluster_count--;
							# Reannotate variant errors respect to the reference sequence of the existing cluster
							# if ($cluster_name eq $previous_cluster) {
							my $previous_dominant_seq = $clusters_amplicon->{$previous_cluster}[0]{'seq'};
							my $previous_dominant_md5 = $clusters_amplicon->{$previous_cluster}[0]{'md5'};
							for (my $j=0; $j<=$#{$cluster_data}; $j++){
								my $variant_seq_ = $cluster_data->[$j]{'seq'};
								my $variant_md5_ = $cluster_data->[$j]{'md5'};
# 								my ($substitutions, $insertions, $deletions, $homopolymer_indels);
# 								if ($fast_multiple_align && defined($unique_aligned_seqs->{$previous_dominant_md5})) {
# 									($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($unique_aligned_seqs->{$previous_dominant_md5},$unique_aligned_seqs->{$variant_md5_});
# 								} else {
								my ($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors(align_2seqs($previous_dominant_seq,$variant_seq_, 'needle'));
								$cluster_data->[$j]{'errors'} = sprintf("sub: %d, ins: %d, del: %d, homo_indel: %d", scalar @$substitutions, scalar @$insertions, scalar @$deletions, scalar @$homopolymer_indels);
							}
							#}
						}#time=0%
						# Reannotate variant errors if the reference sequence has been replaced by the consensus
						if ($consensus_md5 eq $dominant_md5) {
							for (my $j=1; $j<=$#{$cluster_data}; $j++){
								my $variant_seq_ = $cluster_data->[$j]{'seq'};
								my $variant_md5_ = $cluster_data->[$j]{'md5'};
	# 							my ($substitutions, $insertions, $deletions, $homopolymer_indels);
	# 							if ($fast_multiple_align && defined($unique_aligned_seqs->{$consensus_md5})) {
	# 								($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($unique_aligned_seqs->{$consensus_md5},$unique_aligned_seqs->{$variant_md5_});
	# 							} else {
								my ($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors(align_2seqs($consensus_seq,$variant_seq_, 'needle'));
								$cluster_data->[$j]{'errors'} = sprintf("sub: %d, ins: %d, del: %d, homo_indel: %d", scalar @$substitutions, scalar @$insertions, scalar @$deletions, scalar @$homopolymer_indels);
							}
						}#time=0%
					}
				}#time=12% # End consensus section

				# Annotates cluster members
				# $clusters_amplicon->{$cluster_name} = $cluster_data;
				for (my $j=0; $j<=$#{$cluster_data}; $j++) {
					my $md5 = $cluster_data->[$j]{'md5'};
# if ($md5 eq '68a3a2777f6f1d25c248064bcaeadba7'){
# print '';
# }
# 					push(@{$clusters_amplicon->{$cluster_name}}, { 'md5'=>$md5, 'identity'=>$cluster_data->[$j]{'identity'}, 'depth'=>$cluster_data->[$j]{'depth'} });
					push(@{$clusters_amplicon->{$cluster_name}}, $cluster_data->[$j]);
					$variants_assignations{$md5} = $cluster_name;
				}
				# Prints a blank line between clusters in debug mode
				# $clusters_output->{$marker_name}{$sample_name}{$cluster_name} .= "\n";

				# MOVED TO THE BEGGINING FOR THE LOOP: if the dominant seq is in the list of $high_freq_seq_variants (is already a member of some clusters) it will not form a new cluster
				## Removes the dominant sequence from the list of high freq seqs
				#if (defined($high_freq_seq_variants) && defined($high_freq_seq_variants->{$dominant_md5})){
					#delete($high_freq_seq_variants->{$dominant_md5});
				#}

				# Removes clustered variants for next clustering
				my @cluster_md5s = map $_->{'md5'}, @{$cluster_data};
				for (my $j=0; $j<=$#variants_md5s; $j++) {
					my $variant_md5 = $variants_md5s[$j];
					if (in_array(\@cluster_md5s, $variant_md5) && (!defined($high_freq_seq_variants) || !defined($high_freq_seq_variants->{$variant_md5}))){
						splice(@variants_md5s,$j,1);
						splice(@variants_seqs,$j,1);
						$j--;
					}
				}
			}
# printf("Total time: %.2f\n",sum(@run_times));

			# Annotates sequences not clustered
			# Because they are singletons or they are not artifacts from any dominant one and they don't fulfil the dominant conditions
# 			printf("\nTotal not clustered: %d\n",scalar @variants_md5s);
# 			my $non_clustered_depth = 0;
			for (my $i=0; $i<=$#variants_md5s; $i++) {
				my $variant_md5 = $variants_md5s[$i];
# if ($variant_md5 eq 'e8507ff5f2b5775229db67ff0a4d3551'){
# print '';
# }
				if (defined($clustering_thresholds{'max_allele_number'}) && $cluster_count >= $clustering_thresholds{'max_allele_number'}){
					last;
				}
				if (!defined($variants_assignations{$variant_md5}) && (!defined($high_freq_seq_variants) || !defined($high_freq_seq_variants->{$variant_md5}))){
					$cluster_count++;
					my $cluster_name = $cluster_count;
					# Annotates at least one artefact for DOC filtering (if not it will not work properly)
					if ($i == 0 && defined($clustering_thresholds{'degree_of_change'})){
						push(@{$clusters_amplicon->{$cluster_name}},{'md5' => $variant_md5, 'seq' => $variants_seqs[$i], 'depth' => $amplicon_depths{$variant_md5}, 'errors' => 'DOC extra variant'}); # 'identity' => 100,
					} else {
						push(@{$clusters_amplicon->{$cluster_name}},{'md5' => $variant_md5, 'seq' => $variants_seqs[$i], 'depth' => $amplicon_depths{$variant_md5}, 'errors' => 'no cluster'}); # 'identity' => 100,
					}
				}
			}


			# Annotates allele depth and amplicon depth
			my $count_total_clustered_seqs_amplicon = 0;
			foreach my $cluster_name (sort {$a<=>$b} keys %$clusters_amplicon){
				# Annotates all sequences for further printing, even not clustered
				push(@{$clusters->{$marker_name}{$sample_name}}, $clusters_amplicon->{$cluster_name});
				# Finishes when non clustered seqs start to appear
				if ($#{$clusters_amplicon->{$cluster_name}} == 0 && $clusters_amplicon->{$cluster_name}[0]{'errors'} eq 'no cluster'){ # && !defined($clustering_thresholds{'degree_of_change'})){
					next;
				}
				my $allele_name = $clusters_amplicon->{$cluster_name}[0]{'md5'};
				$md5_to_sequence->{$allele_name} = $clusters_amplicon->{$cluster_name}[0]{'seq'};
				my $cluster_depth = 0;
				foreach my $cluster_member (@{$clusters_amplicon->{$cluster_name}}){
					if (is_numeric($cluster_member->{'depth'})){
						$cluster_depth += $cluster_member->{'depth'};
					}
				}
				$amplicon_clustered_sequences->{$marker_name}{$sample_name}{$allele_name} += $cluster_depth;
				$amplicon_clustered_depths->{$marker_name}{$sample_name} += $cluster_depth;
				$count_total_clustered_seqs_amplicon += $cluster_depth;
			}

			if ($sample_name ne ''){
				printf("\t%s-%s clustered (%d reads, %d variants)\n", $marker_name, $sample_name, $count_total_clustered_seqs_amplicon, scalar keys %{$amplicon_clustered_sequences->{$marker_name}{$sample_name}});
			} else {
				printf("\t%s clustered (%d reads, %d variants)\n", $marker_name, $count_total_clustered_seqs_amplicon, scalar keys %{$amplicon_clustered_sequences->{$marker_name}{$sample_name}});
			}

		}

	}


	return ($amplicon_clustered_sequences,$amplicon_clustered_depths,$clusters,$md5_to_sequence);

}

#################################################################################

# Clusters unique sequences to extract real sequences/alleles
sub cluster_amplicon_sequences_with_threads {

	my ($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$referencedata,$threads_limit) = @_;

	if (!defined($threads_limit)){
		$threads_limit = 4;
	}
	
	my ($clustered_amplicon_sequences, $clustered_amplicon_depths, $clusters);
	my $md5_to_sequence = {};
	my $one_amplicon_seq_data;
	foreach my $marker_name (@$markers){
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}
		my $clustered_amplicon_depths->{$marker_name} = {};
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
				next;
			}
			my $one_amplicon_seq_data_;
			my $clustered_amplicon_sequences->{$marker_name}{$sample_name} = {};
			$one_amplicon_seq_data_->{$marker_name}{$sample_name} = $amplicon_seq_data->{$marker_name}{$sample_name};
			push(@{$one_amplicon_seq_data},$one_amplicon_seq_data_);
			# print '';
		}
	}

	my @threads;
	for (my $count_amplicon=0; $count_amplicon<=$#{$one_amplicon_seq_data}; $count_amplicon++){

		push(@threads, threads->create(\&cluster_amplicon_sequences,$markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$one_amplicon_seq_data->[$count_amplicon],$referencedata));
# # 		push(@threads, 1);
# # 		print "\n";

# 		# For debugging:
# 		push(@threads, [cluster_amplicon_sequences($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$one_amplicon_seq_data->[$count_amplicon],$referencedata)]);

# 		# If maximum number of threads is reached or last sbjct of a query is processed
		if (scalar @threads >= $threads_limit  || $count_amplicon == $#{$one_amplicon_seq_data}){
# 			print scalar @threads."\n"; exit;
			my $check_threads = 1;
			while ($check_threads){
				for (my $i=0; $i<=$#threads; $i++){
					unless ($threads[$i]->is_running()){
						my ($clustered_amplicon_sequences_,$clustered_amplicon_depths_,$clusters_,$md5_to_sequence_);
						($clustered_amplicon_sequences_,$clustered_amplicon_depths_,$clusters_,$md5_to_sequence_) = $threads[$i]->join;
						if (defined($clustered_amplicon_sequences_)){
							my $marker_name = (keys %$clustered_amplicon_sequences_)[0];
							my $sample_name = (keys %{$clustered_amplicon_sequences_->{$marker_name}})[0];
	# 						print "\t$marker_name-$sample_name finished\n";
							if (defined($clustered_amplicon_sequences_->{$marker_name}{$sample_name})){
								$clustered_amplicon_sequences->{$marker_name}{$sample_name} = $clustered_amplicon_sequences_->{$marker_name}{$sample_name};
								$clustered_amplicon_depths->{$marker_name}{$sample_name} = $clustered_amplicon_depths_->{$marker_name}{$sample_name};
								$clusters->{$marker_name}{$sample_name} = $clusters_->{$marker_name}{$sample_name};
								if (defined($md5_to_sequence_) && %$md5_to_sequence_){
									$md5_to_sequence = { %$md5_to_sequence, %$md5_to_sequence_ };
								}
							}
						}
						undef($threads[$i]);
						splice(@threads,$i,1);
						$i = $i - 1;
						unless ($count_amplicon == $#{$one_amplicon_seq_data} && @threads){
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
# 				print '';
# 				my ($clustered_amplicon_sequences_,$clustered_amplicon_depths_,$clusters_,$md5_to_sequence_) = @{$threads[$i]};
# 				my $marker_name = (keys %$clustered_amplicon_sequences_)[0];
# # 				my $sample_name = (keys %{$clustered_amplicon_sequences_->{$marker_name}})[0];
# 				$clustered_amplicon_sequences->{$marker_name}{$sample_name} = $clustered_amplicon_sequences_->{$marker_name}{$sample_name};
# 				$clustered_amplicon_depths->{$marker_name}{$sample_name} = $clustered_amplicon_depths_->{$marker_name}{$sample_name};
# 				$clusters->{$marker_name}{$sample_name} = $clusters_->{$marker_name}{$sample_name};
# 				$md5_to_sequence = { %$md5_to_sequence, %$md5_to_sequence_ };
# 				delete $threads[$i];
# 			}
		}
	}
	
	return ($clustered_amplicon_sequences,$clustered_amplicon_depths,$clusters,$md5_to_sequence);
}

#################################################################################

# Clusters sequences de novo, without any primer/tag information
sub cluster_denovo_amplicon_sequences {

	my ($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data) = @_;

	# Output variables
	my ($amplicon_clustered_sequences, $amplicon_clustered_depths, $clusters, $md5_to_sequence);

	# Loops markers
	foreach my $marker_name (@$markers){

		# Skips markers without data
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}

		# Extracts clustering parameters
		my %clustering_thresholds;
		foreach my $clustering_threshold (keys %$paramsdata) {
			if (defined($paramsdata->{$clustering_threshold}{$marker_name})){
				$clustering_thresholds{$clustering_threshold} = $paramsdata->{$clustering_threshold}{$marker_name}[0];
			} elsif (defined($paramsdata->{$clustering_threshold}{'all'})){
				$clustering_thresholds{$clustering_threshold} = $paramsdata->{$clustering_threshold}{'all'}[0];
			}
		}

		# Loops samples/amplicons and performs clustering rounds
		foreach my $sample_name (@{$samples->{$marker_name}}) {

			# Process only samples with sequences (after filtering) in the original order
			if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
# 				print "\t$marker_name-$sample_name doesn't have sequences to cluster.\n";
				next;
			}

			if ($sample_name ne ''){
				printf("\t%s-%s clustering\n",$marker_name,$sample_name);
			} else {
				printf("\t%s clustering\n",$marker_name);
			}

			# Includes the sequences into $amplicon_seq_data
			foreach my $md5 (keys %{$amplicon_seq_data->{$marker_name}{$sample_name}}){
				$amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'seq'} = $marker_seq_data->{$marker_name}{$md5}{'seq'};
			}

			# Clusters the amplicons
			my $clusters_amplicon = cluster_denovo_sequences($amplicon_seq_data->{$marker_name}{$sample_name},\%clustering_thresholds);

			# Annotates allele depth and amplicon depth
			my $count_total_clustered_seqs_amplicon = 0;
			foreach my $cluster_name (sort {$a<=>$b} keys %$clusters_amplicon){
				# Annotates all sequences for further printing, even not clustered
				push(@{$clusters->{$marker_name}{$sample_name}}, $clusters_amplicon->{$cluster_name});
				# Finishes when non clustered seqs start to appear
				if ($#{$clusters_amplicon->{$cluster_name}} == 0 && $clusters_amplicon->{$cluster_name}[0]{'errors'} eq 'no cluster'){ # && !defined($clustering_thresholds{'degree_of_change'})){
					next;
				}
				my $allele_name = $clusters_amplicon->{$cluster_name}[0]{'md5'};
				$md5_to_sequence->{$allele_name} = $clusters_amplicon->{$cluster_name}[0]{'seq'};
				my $cluster_depth = 0;
				foreach my $cluster_member (@{$clusters_amplicon->{$cluster_name}}){
					if (is_numeric($cluster_member->{'depth'})){
						$cluster_depth += $cluster_member->{'depth'};
					}
				}
				$amplicon_clustered_sequences->{$marker_name}{$sample_name}{$allele_name} += $cluster_depth;
				$amplicon_clustered_depths->{$marker_name}{$sample_name} += $cluster_depth;
				$count_total_clustered_seqs_amplicon += $cluster_depth;
			}

			if ($sample_name ne ''){
				printf("\t%s-%s clustered (%d reads, %d variants)\n", $marker_name, $sample_name, $count_total_clustered_seqs_amplicon, scalar keys %{$amplicon_clustered_sequences->{$marker_name}{$sample_name}});
			} else {
				printf("\t%s clustered (%d reads, %d variants)\n", $marker_name, $count_total_clustered_seqs_amplicon, scalar keys %{$amplicon_clustered_sequences->{$marker_name}{$sample_name}});
			}
		}
	}

	return ($amplicon_clustered_sequences,$amplicon_clustered_depths,$clusters,$md5_to_sequence);

}


#################################################################################

# Clusters sequences de novo, without any primer/tag information
sub cluster_denovo_amplicon_sequences_with_threads {

	my ($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$threads_limit) = @_;

	if (!defined($threads_limit)){
		$threads_limit = 4;
	}
	
	my ($clustered_amplicon_sequences, $clustered_amplicon_depths, $clusters);
	my $md5_to_sequence = {};
	my $one_amplicon_seq_data;
	foreach my $marker_name (@$markers){
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}
		my $clustered_amplicon_depths->{$marker_name} = {};
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
				next;
			}
			my $one_amplicon_seq_data_;
			my $clustered_amplicon_sequences->{$marker_name}{$sample_name} = {};
			$one_amplicon_seq_data_->{$marker_name}{$sample_name} = $amplicon_seq_data->{$marker_name}{$sample_name};
			push(@{$one_amplicon_seq_data},$one_amplicon_seq_data_);
			# print '';
		}
	}

	my @threads;
	for (my $count_amplicon=0; $count_amplicon<=$#{$one_amplicon_seq_data}; $count_amplicon++){

		push(@threads, threads->create(\&cluster_denovo_amplicon_sequences,$markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$one_amplicon_seq_data->[$count_amplicon]));
# # 		push(@threads, 1);
# # 		print "\n";

# 		# For debugging:
# 		push(@threads, [cluster_amplicon_sequences($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$one_amplicon_seq_data->[$count_amplicon])]);

# 		# If maximum number of threads is reached or last sbjct of a query is processed
		if (scalar @threads >= $threads_limit  || $count_amplicon == $#{$one_amplicon_seq_data}){
# 			print scalar @threads."\n"; exit;
			my $check_threads = 1;
			while ($check_threads){
				for (my $i=0; $i<=$#threads; $i++){
					unless ($threads[$i]->is_running()){
						my ($clustered_amplicon_sequences_,$clustered_amplicon_depths_,$clusters_,$md5_to_sequence_);
						($clustered_amplicon_sequences_,$clustered_amplicon_depths_,$clusters_,$md5_to_sequence_) = $threads[$i]->join;
						if (defined($clustered_amplicon_sequences_)){
							my $marker_name = (keys %$clustered_amplicon_sequences_)[0];
							my $sample_name = (keys %{$clustered_amplicon_sequences_->{$marker_name}})[0];
	# 						print "\t$marker_name-$sample_name finished\n";
							if (defined($clustered_amplicon_sequences_->{$marker_name}{$sample_name})){
								$clustered_amplicon_sequences->{$marker_name}{$sample_name} = $clustered_amplicon_sequences_->{$marker_name}{$sample_name};
								$clustered_amplicon_depths->{$marker_name}{$sample_name} = $clustered_amplicon_depths_->{$marker_name}{$sample_name};
								$clusters->{$marker_name}{$sample_name} = $clusters_->{$marker_name}{$sample_name};
								if (defined($md5_to_sequence_) && %$md5_to_sequence_){
									$md5_to_sequence = { %$md5_to_sequence, %$md5_to_sequence_ };
								}
							}
						}
						undef($threads[$i]);
						splice(@threads,$i,1);
						$i = $i - 1;
						unless ($count_amplicon == $#{$one_amplicon_seq_data} && @threads){
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
# 				print '';
# 				my ($clustered_amplicon_sequences_,$clustered_amplicon_depths_,$clusters_,$md5_to_sequence_) = @{$threads[$i]};
# 				my $marker_name = (keys %$clustered_amplicon_sequences_)[0];
# # 				my $sample_name = (keys %{$clustered_amplicon_sequences_->{$marker_name}})[0];
# 				$clustered_amplicon_sequences->{$marker_name}{$sample_name} = $clustered_amplicon_sequences_->{$marker_name}{$sample_name};
# 				$clustered_amplicon_depths->{$marker_name}{$sample_name} = $clustered_amplicon_depths_->{$marker_name}{$sample_name};
# 				$clusters->{$marker_name}{$sample_name} = $clusters_->{$marker_name}{$sample_name};
# 				$md5_to_sequence = { %$md5_to_sequence, %$md5_to_sequence_ };
# 				delete $threads[$i];
# 			}
		}
	}
	
	return ($clustered_amplicon_sequences,$clustered_amplicon_depths,$clusters,$md5_to_sequence);
}

#################################################################################

# Clusters sequences from a file de novo, without any primer/tag information
sub cluster_denovo_sequence_file {

	my ($seqs_file,$paramsdata) = (@_);

	my ($sequences,$headers) = read_sequence_file($seqs_file);

	return cluster_denovo_sequences($sequences,$paramsdata);

}

#################################################################################

# Clusters sequences de novo, without any primer/tag information
sub cluster_denovo_sequences {

	my ($sequences,$paramsdata) = @_;
	
	# $sequences may be an ARRAY with sequences or a HASH with amplicon data

	# Extracts clustering parameters
	my %clustering_thresholds = %$paramsdata;

	# Loops all the sequences and annotates the md5s, depths and lengths
	my (%md5_to_seq, %md5_to_depth, %md5_to_length);
	my $total_depth = 0;
	if (ref($sequences) eq 'ARRAY'){
		foreach my $seq (@$sequences) {
			my $md5 = generate_md5($seq);
			if (defined($md5_to_depth{$md5})){
				$md5_to_depth{$md5}++;
				$total_depth++;
			} else {
				$md5_to_seq{$md5} = $seq;
				$md5_to_depth{$md5} = 1;
				$md5_to_length{$md5} = length($seq);
				$total_depth++;
			}
		}
	} elsif (ref($sequences) eq 'HASH'){
		foreach my $md5 (keys %$sequences) {
			$md5_to_seq{$md5} = $sequences->{$md5}{'seq'};
			$md5_to_depth{$md5} = $sequences->{$md5}{'depth'};
			$md5_to_length{$md5} = length($sequences->{$md5}{'seq'});
			$total_depth += $sequences->{$md5}{'depth'};
		}
	} else {
		print "\nERROR 'cluster_denovo_sequences': Invalid input data.\n";
		exit;
	}


# FIRST CLUSTERING GROUPS ONLY SAME LENGTH SEQUENCES BASED IN IDENTITY AND SUBSTITUTIONS - FAST

	my $similar_md5s;

	# Sorts sequences by depth
	my @sorted_md5s = sort { $md5_to_depth{$b} <=> $md5_to_depth{$a} } keys %md5_to_depth;

	# Variable to update the cluster depth
	my %cluster_depths;
	
	my ($count_clusters, $count_total_clustered_variants, $count_total_clustered_reads) = (0, 0, 0);
	while (@sorted_md5s) {

		# Data from the highest frequency variant
		my $md5 = shift @sorted_md5s;
		my $seq = $md5_to_seq{$md5};
		my $seq_len = $md5_to_length{$md5};
		my $seq_depth = $md5_to_depth{$md5};

# if ($md5 eq '6152c9d8beda2efcb458b18e9b9b03fd'){
# print '';
# }
		# Do not cluster singletons and stop clustering at very low coverages (0.002%) - saves time without altering results
		if (!defined($clustering_thresholds{'keep_singletons'}) || !$clustering_thresholds{'keep_singletons'}){
			if ($seq_depth == 1 || $seq_depth/$total_depth < 0.00002) {
				last;
			}
		}

		# Annotates cluster member data
		$count_clusters++;
		push(@{$similar_md5s->{$md5}},{'md5' => $md5, 'seq' => $seq, 'depth' => $seq_depth, 'errors' => ''});
		$cluster_depths{$md5} = $seq_depth;
		$count_total_clustered_variants++;
		$count_total_clustered_reads += $seq_depth;

		my ($max_errors,$max_substs,$max_indels,$min_ident) = (undef,$seq_len,$seq_len,0);
		# Sets maximum number of total errors allowed, minimum identity
		if (defined($clustering_thresholds{'identity_threshold'})){
			$min_ident = sprintf("%.0f", $clustering_thresholds{'identity_threshold'}/100*$seq_len);
			if (!defined($max_errors) || $seq_len-$min_ident>$max_errors) {
				$max_errors = $seq_len-$min_ident;
			}
		# Sets maximum number of substitution errors allowed
		} elsif (defined($clustering_thresholds{'substitution_threshold'})){
			$max_substs = sprintf("%.0f", $clustering_thresholds{'substitution_threshold'}/100*$seq_len);
			if (!defined($max_errors) || $max_substs>$max_errors) {
				$max_errors = $max_substs;
			}
		# Sets maximum number of substitution errors allowed
		} elsif (defined($clustering_thresholds{'max_substitutions'})){
			$max_substs = $clustering_thresholds{'max_substitutions'};
			if (!defined($max_errors) || $max_substs>$max_errors) {
				$max_errors = $max_substs;
			}
		}

		# Selects variants with the same length and similar to make clustering faster
		my @similar_md5s;
		for (my $i=0; $i<=$#sorted_md5s; $i++) {

# 			# Do not cluster more variants if the cluster depth is more than double than the highest frequency variant depth
# 			if ($cluster_depths{$md5} > 2*$seq_depth) {
# 				last;
# 			}

			my $md5_ = $sorted_md5s[$i];
			my $seq_ = $md5_to_seq{$md5_};
			my $seq_len_ = $md5_to_length{$md5_};
			my $seq_depth_ = $md5_to_depth{$md5_};

# if ($md5 eq 'dc4d714c37154137dd34e1282f42476a' && $md5_ =~ /6152c9d8beda2efcb458b18e9b9b03fd|a19d39b9ed247d6472cd1469fc96f3fa|2e628e48997f015ca57cbc9deab64a0b/){
# print '';
# }

			# Skips different length variants
			if ($seq_len_ != $seq_len){
				next;
			}
			# Skips non similar variants
			my $errors = count_errors($seq,$seq_,$max_errors+1);
			# if (count_errors_c($dominant_seq,$variant_seq,$max_errors+1) > $max_errors) { # 15% faster C code
			if ($errors > $max_errors) { # Perl code
				next;
			# Skips same depth variants
			} elsif ($seq_depth_ == $seq_depth){
				next;
			# Skips high amplicon depth/frequency variants
			} elsif (defined($clustering_thresholds{'min_amplicon_seq_frequency_threshold'}) && 100*$seq_depth_/$total_depth >= $clustering_thresholds{'min_amplicon_seq_frequency_threshold'}){
				next;
			# Skips high depth/frequency variants respect to the reference
			} elsif (defined($clustering_thresholds{'min_dominant_frequency_threshold'}) && 100*$seq_depth_/$seq_depth >= $clustering_thresholds{'min_dominant_frequency_threshold'}){
				next;
			# Clusters variant
			} else {
				# Annotates cluster member data
				push(@{$similar_md5s->{$md5}},{'md5' => $md5_, 'seq' => $seq_, 'depth' => $seq_depth_, 'errors' => "sub: $errors"}); #  'identity' => $identity,
				$cluster_depths{$md5} += $seq_depth_;
				$count_total_clustered_variants++;
				$count_total_clustered_reads += $seq_depth_;
				splice(@sorted_md5s,$i,1);
				$i--;
			}
		}
	}
	# DEBUGGING:
	# printf("First clustering: %d clusters, %d variants, %d reads\n", $count_clusters, $count_total_clustered_variants, $count_total_clustered_reads);

# SECOND CLUSTERING GROUPS PREVIOUSLY CLUSTERED SEQUENCES CONSIDERING INDELS - SLOW

	my $clusters;

	# Sorts clustered sequences by depth
	@sorted_md5s = sort { $cluster_depths{$b} <=> $cluster_depths{$a} } keys %$similar_md5s;
	# Array to order and rename the clusters later
	my @sorted_cluster_md5s = @sorted_md5s;

	($count_clusters, $count_total_clustered_variants, $count_total_clustered_reads) = (0, 0, 0);
	while (@sorted_md5s) {

		# Data from the highest frequency variant
		my $md5 = shift @sorted_md5s;
		my $seq = $md5_to_seq{$md5};
		my $seq_len = $md5_to_length{$md5};
		my $seq_depth = $cluster_depths{$md5}; # Depth of the cluster created in the first clustering round
		# Variable to store data from cluster sequences
		my $cluster_data;
		# Annotates cluster member data
		$count_clusters++;
		$clusters->{$count_clusters} = $similar_md5s->{$md5};
		$count_total_clustered_variants += scalar @{$similar_md5s->{$md5}};
		$count_total_clustered_reads += $seq_depth;

		# Stop clustering when there are as many clusters as expected alleles (saves time without altering results)
		if (!@sorted_md5s || defined($clustering_thresholds{'max_allele_number'}) && $count_clusters >= $clustering_thresholds{'max_allele_number'}){
			last;
		}

		# Sets maximum number of total errors allowed, minimum identity
		my ($max_substs,$max_indels,$min_ident) = ($seq_len,$seq_len,0);
		if (defined($clustering_thresholds{'identity_threshold'})){
			$min_ident = sprintf("%.0f", $clustering_thresholds{'identity_threshold'}/100*$seq_len);
		} else {
			# Sets maximum number of substitution errors allowed
			if (defined($clustering_thresholds{'substitution_threshold'})){
				$max_substs = sprintf("%.0f", $clustering_thresholds{'substitution_threshold'}/100*$seq_len);
			# Sets maximum number of substitution errors allowed
			} elsif (defined($clustering_thresholds{'max_substitutions'})){
				$max_substs = $clustering_thresholds{'max_substitutions'};
			}
			# Sets maximum number of indel errors allowed
			if (defined($clustering_thresholds{'indel_threshold'})){
				$max_indels = sprintf("%.0f", $clustering_thresholds{'indel_threshold'}/100*$seq_len);
			}
		}

		# Aligns seqs to reference if fast mode is not active
		my @variant_seqs = map $md5_to_seq{$_}, @sorted_md5s;
		my $aligned_seqs = align_seqs2one($seq,$md5,\@variant_seqs,\@sorted_md5s,'needleall');

		# Compares similar sequences with the high frequency ones
		for (my $i=0; $i<=$#sorted_md5s; $i++) {

# 			# Do not cluster more variants if the cluster depth is more than double than the highest frequency variant depth
# 			if ($cluster_depths{$md5} > 2*$seq_depth) {
# 				last;
# 			}

			my $md5_ = $sorted_md5s[$i];
			my $seq_len_ = $md5_to_length{$md5_};
			my $seq_depth_ = $cluster_depths{$md5_}; # Depth of the cluster created in the first clustering round

# if ($md5 eq 'dc4d714c37154137dd34e1282f42476a' && $md5_ =~ /6152c9d8beda2efcb458b18e9b9b03fd|a19d39b9ed247d6472cd1469fc96f3fa|2e628e48997f015ca57cbc9deab64a0b/){
# print '';
# }
			# Skip sequences with higher coverages
			if ($seq_depth_>=$seq_depth){
				next;
			# Skips high amplicon depth/frequency variants
			} elsif (defined($clustering_thresholds{'min_amplicon_seq_frequency_threshold'}) && 100*$seq_depth_/$total_depth >= $clustering_thresholds{'min_amplicon_seq_frequency_threshold'}){
				next;
			# Skips high depth/frequency variants respect to the reference
			} elsif (defined($clustering_thresholds{'min_dominant_frequency_threshold'}) && 100*$seq_depth_/$seq_depth >= $clustering_thresholds{'min_dominant_frequency_threshold'}){
				next;
			# Skip if the pairwise global alignment fails
			} elsif (!defined($aligned_seqs->{$md5_}[0]) || !defined($aligned_seqs->{$md5_}[1])){
				next;
			}

			# Skips seqs with lower identity than threshold
			my $error_data = '';
			if (defined($clustering_thresholds{'identity_threshold'})){
				my ($identical,$total) = binary_score_nts($aligned_seqs->{$md5_}[0],$aligned_seqs->{$md5_}[1]);
				if ($identical < $min_ident){
					next;
				} else {
					$error_data = sprintf("total: %d, ident: %d", $total, $identical);
				}
			# Skips seqs with more errors than substitution and indel thresholds
			} elsif (defined($clustering_thresholds{'substitution_threshold'}) || defined($clustering_thresholds{'indel_threshold'}) || defined($clustering_thresholds{'max_substitutions'})){
				my ($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($aligned_seqs->{$md5_}[0],$aligned_seqs->{$md5_}[1]);
				# Skips seqs with more substitutions than threshold
				if (scalar @$substitutions > $max_substs){
					next;
				}
				# Skips seqs with more non homopolymer indels than threshold and process them in next clustering round
				my $non_homopolymer_indels = scalar @$insertions + scalar @$deletions - scalar @$homopolymer_indels;
				if ($non_homopolymer_indels > $max_indels){
					next;
				}
				# Prints artifact sequence and errors
				$error_data = sprintf("sub: %d, ins: %d, del: %d, homo_indel: %d", scalar @$substitutions, scalar @$insertions, scalar @$deletions, scalar @$homopolymer_indels);
			}
			# Annotates cluster member data
			push(@{$clusters->{$count_clusters}},@{$similar_md5s->{$md5_}}); #  'identity' => $identity,
			$cluster_depths{$md5} += $seq_depth_;
			$count_total_clustered_reads += $seq_depth_;
			$count_total_clustered_variants += scalar @{$similar_md5s->{$md5_}};
			splice(@sorted_md5s,$i,1);
			$i--;
		}

	}

	return $clusters;


}

#################################################################################

# OBSOLETE:
# Clusters sequences de novo, without any primer/tag information
sub cluster_denovo_sequences_old {

	my ($sequences,$paramsdata) = @_;
	
	# $sequences may be an ARRAY with sequences or a HASH with amplicon data

	# Extracts clustering parameters
	my %clustering_thresholds = %$paramsdata;

	# Loops all the sequences and annotates the md5s, depths, lengths and 20 first nts
	my (@md5s, %md5_to_seq, %md5_to_depth, %md5_to_length, %md5_to_20nts);
	my $total_depth = 0;
	if (ref($sequences) eq 'ARRAY'){
		foreach my $seq (@$sequences) {
			my $md5 = generate_md5($seq);
			if (defined($md5_to_depth{$md5})){
				$md5_to_depth{$md5}++;
				$total_depth++;
			} else {
				push(@md5s,$md5);
				$md5_to_seq{$md5} = $seq;
				$md5_to_depth{$md5} = 1;
				$md5_to_length{$md5} = length($seq);
				$md5_to_20nts{$md5} = substr($seq,0,20);
				$total_depth++;
			}
		}
	} elsif (ref($sequences) eq 'HASH'){
		foreach my $md5 (keys %$sequences) {
			push(@md5s,$md5);
			$md5_to_seq{$md5} = $sequences->{$md5}{'seq'};
			$md5_to_depth{$md5} = $sequences->{$md5}{'depth'};
			$md5_to_length{$md5} = length($sequences->{$md5}{'seq'});
			$md5_to_20nts{$md5} = substr($sequences->{$md5}{'seq'},0,20);
			$total_depth += $sequences->{$md5}{'depth'};
		}
	} else {
		print "\nERROR 'cluster_denovo_sequences': Invalid input data.\n";
		exit;
	}

	# Sorts sequences by depth
	my @sorted_md5s = sort { $md5_to_depth{$b} <=> $md5_to_depth{$a} } @md5s;
	my @sorted_cluster_md5s = @sorted_md5s; # Array to order and rename the clusters later

	# Generates clusters of similar sequences
	my $clusters; # Variable to store clusters
	# Variables to count clusters, clustered, variants, reads...
	my ($count_rounds, $count_clusters, $count_total_clustered_variants, $count_total_clustered_reads) = (0, 0, 0, 0);
	while (@sorted_md5s) {

		$count_rounds++;

		# Data from the highest frequency variant
		my $md5 = shift @sorted_md5s;
		my $seq = $md5_to_seq{$md5};
		my $seq_len = $md5_to_length{$md5};
		my $seq_depth = $md5_to_depth{$md5};


		# Stop clustering when there are as many clusters as expected alleles (saves time without altering results)
		if (defined($clustering_thresholds{'max_allele_number'}) && $count_clusters >= $clustering_thresholds{'max_allele_number'}){
			last;
		}

		# Do not cluster singletons and stop clustering at very low coverages (0.002%) - saves time without altering results
		if (!(defined($clustering_thresholds{'keep_singletons'}) && $clustering_thresholds{'keep_singletons'}) ){
			if ($seq_depth == 1 || $seq_depth/$total_depth < 0.00002) {
				last;
			}
		}

		# Sets maximum number of substitution errors allowed
		my ($max_errors,$max_substs,$max_indels,$min_ident) = (undef,$seq_len,$seq_len,0);
		if (defined($clustering_thresholds{'substitution_threshold'})){
			$max_substs = sprintf("%.0f", $clustering_thresholds{'substitution_threshold'}/100*$seq_len);
			if (!defined($max_errors) || $max_substs>$max_errors) {
				$max_errors = $max_substs;
			}
		}
		# Sets maximum number of indel errors allowed
		if (defined($clustering_thresholds{'indel_threshold'})){
			$max_indels = sprintf("%.0f", $clustering_thresholds{'indel_threshold'}/100*$seq_len);
			if (!defined($max_errors) || $max_indels>$max_errors) {
				$max_errors = $max_indels;
			}
		}
		# Sets maximum number of total errors allowed, minimum identity
		if (defined($clustering_thresholds{'identity_threshold'})){
			$min_ident = sprintf("%.0f", $clustering_thresholds{'identity_threshold'}/100*$seq_len);
			if (!defined($max_errors) || $seq_len-$min_ident>$max_errors) {
				$max_errors = $seq_len-$min_ident;
			}
		}

		# Activates fast mode if indel threshold is minimum, no global alignment is required
		my $fast_align = 0;
		if ((defined($max_indels) && $max_indels < 1) || (defined($clustering_thresholds{'fast_clustering_alignment'}) && $clustering_thresholds{'fast_clustering_alignment'})){
			$fast_align = 1;
		}

		# Selects similar sequences to make clustering faster
		my @similar_md5s;
		for (my $i=0; $i<=$#sorted_md5s; $i++) {
			my $seq_len_ = $md5_to_length{$sorted_md5s[$i]};
			if ($fast_align){
				if ($seq_len_ != $seq_len){
					next;
# 				} elsif (count_errors_c($dominant_seq,$variant_seq,$max_errors+1) > $max_errors) { # 15% faster C code
				} elsif (count_errors($seq,$md5_to_seq{$sorted_md5s[$i]},$max_errors+1) > $max_errors) { # Perl code
					next;
				} else {
					push(@similar_md5s,$sorted_md5s[$i]);
				}
			} else {
				if (abs($seq_len_-$seq_len) > $max_errors){
					next;
				} elsif ($seq_len_ == $seq_len && count_errors($seq,$md5_to_seq{$sorted_md5s[$i]},$max_errors+1) > $max_errors){
					next;
				# Annotates high depth/frequency variants (before clustering) into different clusters
				} elsif (defined($clustering_thresholds{'min_amplicon_seq_frequency_threshold'}) && 100*$md5_to_depth{$sorted_md5s[$i]}/$total_depth >= $clustering_thresholds{'min_amplicon_seq_frequency_threshold'}){
					next;
				# Annotates high depth/frequency variants respect to the reference (before clustering) into different clusters
				} elsif (defined($clustering_thresholds{'min_dominant_frequency_threshold'}) && 100*$md5_to_depth{$sorted_md5s[$i]}/$seq_depth >= $clustering_thresholds{'min_dominant_frequency_threshold'}){
					next;
				} else {
					push(@similar_md5s,$sorted_md5s[$i]);
				}
			}
		}
		my @similar_seqs = map $md5_to_seq{$_}, @similar_md5s;

		# Detects high frequency variants among similar sequences
		# Fist is always the highest frequency variant
		my @high_freq_md5s = ($md5); # The same: @high_freq_md5s = ($similar_md5s[0]); 

		# Annotates clusters
		my $cluster_data; # Variable to store data from cluster sequences
		foreach my $high_freq_md5 (@high_freq_md5s){

			# Aligns seqs to reference if fast mode is not active
			my $aligned_seqs;
			if (!$fast_align && @similar_md5s) {
				$aligned_seqs = align_seqs2one($md5_to_seq{$high_freq_md5},$high_freq_md5,\@similar_seqs,\@similar_md5s,'needleall');
			}

			# Compares similar sequences with high frequency ones
			for (my $i=0; $i<=$#similar_md5s; $i++) {

				my $md5_ = $similar_md5s[$i];

				# If fast mode is active, both sequences should have same length and no indels
				if ($fast_align) {
					$aligned_seqs->{$md5_} = [ $md5_to_seq{$high_freq_md5}, $md5_to_seq{$md5_} ]
				}
				
# 				# Skip if the pairwise global alignment fails
				if (!defined($aligned_seqs->{$md5_}[0]) || !defined($aligned_seqs->{$md5_}[1])){
					next;
				}

				# Skip the same sequence, high frequency ones or sequences with higher coverages
				if ($high_freq_md5 eq $md5_ || in_array(\@high_freq_md5s,$md5_) || $md5_to_depth{$md5_}>$md5_to_depth{$high_freq_md5}){
					next;
				}

				# Skips seqs with lower identity than threshold
				my ($identical,$total) = binary_score_nts($aligned_seqs->{$md5_}[0],$aligned_seqs->{$md5_}[1]);
				if ($identical < $min_ident){
					next;
				}

				# Skips seqs with more errors than substitution and indel thresholds
				my ($error_data, $substitutions, $insertions, $deletions, $homopolymer_indels);
				if (defined($clustering_thresholds{'identity_threshold'})){
					$error_data = sprintf("total: %d, ident: %d", $total, $identical);
				} elsif (defined($clustering_thresholds{'substitution_threshold'}) || defined($clustering_thresholds{'indel_threshold'})){
					($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($aligned_seqs->{$md5_}[0],$aligned_seqs->{$md5_}[1]);
					# Skips seqs with more substitutions than threshold
					if (scalar @$substitutions > $max_substs){
						next;
					}
					# Skips seqs with more non homopolymer indels than threshold and process them in next clustering round
					my $non_homopolymer_indels = scalar @$insertions + scalar @$deletions - scalar @$homopolymer_indels;
					if ($non_homopolymer_indels > $max_indels){
						next;
					}
					# Prints artifact sequence and errors
					$error_data = sprintf("sub: %d, ins: %d, del: %d, homo_indel: %d", scalar @$substitutions, scalar @$insertions, scalar @$deletions, scalar @$homopolymer_indels);
				}

				# Annotates cluster member data
				push(@{$cluster_data->{$high_freq_md5}},{'md5' => $md5_, 'seq' => $md5_to_seq{$md5_}, 'depth' => $md5_to_depth{$md5_}, 'errors' => $error_data}); #  'identity' => $identity,

			}
		}

		# Annotates cluster members
# 		$count_clusters++;
# 		for (my $i=0; $i<=$#{$cluster_data->{$md5}}; $i++) {
# 			push(@{$clusters->{$count_clusters}}, $cluster_data->{$md5}[$i]);
# 		}
# 		my %clustered_md5s = map { $_->{'md5'} => 1 } @{$cluster_data->{$md5}};
		my ($count_clustered_variants, $count_clustered_reads) = (0, 0);
		my %clustered_md5s; # Annotates clustered MD5s to remove from next stages
		foreach my $high_freq_md5 (@high_freq_md5s){
			# Creates a new cluster if it doesn't exist previously
			if (!defined($clusters->{$high_freq_md5})){
				$count_clusters++;
				$clusters->{$high_freq_md5} = [{'md5' => $high_freq_md5, 'seq' => $md5_to_seq{$high_freq_md5}, 'depth' => $md5_to_depth{$high_freq_md5}, 'errors' => ''}]; #  'identity' => $identity,
				$count_clustered_variants++;
				$count_clustered_reads += $md5_to_depth{$high_freq_md5};
			}
			# Annotates low frequency variants into the cluster
			for (my $i=0; $i<=$#{$cluster_data->{$high_freq_md5}}; $i++) {
				push(@{$clusters->{$high_freq_md5}}, $cluster_data->{$high_freq_md5}[$i]);
				$count_clustered_variants++;
				$count_clustered_reads += $cluster_data->{$high_freq_md5}[$i]{'depth'};
			}
			my %clustered_md5s_ = map { $_->{'md5'} => 1 } @{$cluster_data->{$high_freq_md5}};
			%clustered_md5s = ( %clustered_md5s, %clustered_md5s_);
		}
		$count_total_clustered_variants += $count_clustered_variants;
		$count_total_clustered_reads += $count_clustered_reads;

		# Removes clustered variants for next clustering
		# Except high freq variants (the highest one has been already removed at the beginning)
		my %keep_md5s = map { $_ => 1 } @high_freq_md5s;
		for (my $i=0; $i<=$#sorted_md5s; $i++) {
			if (defined($clustered_md5s{$sorted_md5s[$i]}) && !defined($keep_md5s{$sorted_md5s[$i]})){
				splice(@sorted_md5s,$i,1);
				$i--;
			}
		}
		# DEBUGGING:
# 		printf("Cluster %d, %d reads clustered in this round, %d remaining reads, %d variants clustered in this round, %d remaining variants\n", $count_clusters, $count_clustered_reads, $total_depth-$count_total_clustered_reads, $count_clustered_variants, scalar @sorted_md5s);
		# if ($count_rounds == 10){
		# last;
		# }
	}
	
	# Replaces MD5s as cluster names for ordered numbers
	for (my $i=0; $i<=$#sorted_cluster_md5s; $i++){
		if (defined($clusters->{$sorted_cluster_md5s[$i]})){
			$clusters->{$i+1} = delete($clusters->{$sorted_cluster_md5s[$i]});
		}
		# Remove extra clusters when there are as many clusters as expected alleles
		if (defined($clustering_thresholds{'max_allele_number'}) && $i+1 > $clustering_thresholds{'max_allele_number'}){
			delete($clusters->{$sorted_cluster_md5s[$i]});
		}
	}

	return $clusters;

}

#################################################################################

# Checks if a sequence pass clustering thresholds
sub is_dominant {

	my ($sequence,$correct_lengths,$clustering_thresholds,$reference_seq) = @_;
	
	my $length = length($sequence);

	# First checks if length is in the range of marker lengths
	if (defined($clustering_thresholds->{'cluster_exact_length'}) && !in_array($correct_lengths,$length)){
		return 0;
	}
	# Second checks if length is in frame or not
	if (defined($clustering_thresholds->{'cluster_inframe'})){
		my $inframe = 0;
		foreach my $marker_len (@$correct_lengths){
			# Checks if length is in the range of 10%-190% of the marker length (to exclude very short or very long unusual seqs)
			if (abs($marker_len-$length) > 0.9*$marker_len){
				next;
			}
			if (abs($marker_len-$length) % 3 == 0){
				$inframe = 1;
				last;
			}
		}
		if (!$inframe){
			return 0;
		}
	}
	# Third checks if the sequence has stop codons
	if (defined($clustering_thresholds->{'cluster_nonstop'})){
		if (index(dna_to_prot(substr($sequence,$clustering_thresholds->{'cluster_nonstop'}-1)),'*') != -1) {
			return 0;
		}
	}

	# If reference sequences are provided, checks if the sequence is compatible with the best reference
	if (defined($reference_seq)){
		my ($substitutions, $insertions, $deletions, $homopolymer_indels);
		if (ref($reference_seq) eq 'ARRAY'){
			($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($reference_seq->[0],$reference_seq->[1]);
		} else {
			($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors(align_2seqs($sequence,$reference_seq,'needle'));
		}
# 		if (scalar @$insertions + scalar @$deletions > 0){
		# Allow a compensaroty indel in a range of 9bps
		if (scalar @$insertions == 1 && scalar @$deletions ==1 && abs($insertions->[0]-$deletions->[0]) <= 9){
			;
		# Only insertions or deletions multiple of 3nts are allowed
		} elsif (abs(scalar @$insertions) % 3 != 0 || abs(scalar @$deletions) % 3 != 0) {
			return 0;
		} elsif (@$insertions || @$deletions){
			while (my @insertions_ = splice(@$insertions,0,3)) {
				if ($insertions_[-1]-$insertions_[0]>2){
					return 0;
				}
			}
			while (my @deletions_ = splice(@$deletions,0,3)) {
				if ($deletions_[-1]-$deletions_[0]>2){
					return 0;
				}
			}
		}
	}

	return 1;

}

#################################################################################

# Compares high depth unique sequences among them
sub compare_amplicon_sequences {

	my ($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data) = @_;

	my $amplicon_seq_comparison_data;

	foreach my $marker_name (@$markers){

		# Skips markers without data
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}

		# Extracts clustering parameters
		# Only used: 'substitution_threshold', 'indel_threshold' and 'identity_threshold'
		my %clustering_thresholds;
		foreach my $clustering_threshold (keys %$paramsdata) {
			if (defined($paramsdata->{$clustering_threshold}{$marker_name})){
				$clustering_thresholds{$clustering_threshold} = $paramsdata->{$clustering_threshold}{$marker_name}[0];
			} elsif (defined($paramsdata->{$clustering_threshold}{'all'})){
				$clustering_thresholds{$clustering_threshold} = $paramsdata->{$clustering_threshold}{'all'}[0];
			}
		}

		# Loops samples/amplicons and performs comparisons
		foreach my $sample_name (@{$samples->{$marker_name}}) {

			# Process only samples with sequences (after filtering) in the original order
			if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
# 				print "\t$marker_name-$sample_name doesn't have sequences to compare.\n";
				next;
			}
			if ($sample_name ne ''){
				print "\t$marker_name-$sample_name pairwise comparing sequences\n";
			} else {
				print "\t$marker_name pairwise comparing sequences\n";
			}

			# Reads amplicon unique sequences, full headers and depths and stores them into hashes with md5s as keys
			my (%amplicon_sequences, %amplicon_depths, %amplicon_frequencies, %amplicon_names);
			foreach my $md5 (keys %{$amplicon_seq_data->{$marker_name}{$sample_name}}) {
				$amplicon_names{$md5} = $marker_seq_data->{$marker_name}{$md5}{'name'};
				$amplicon_sequences{$md5} = $marker_seq_data->{$marker_name}{$md5}{'seq'};
				$amplicon_depths{$md5} = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
				$amplicon_frequencies{$md5} = sprintf("%.2f", $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'});
			}

# 			print "...";

			# Sorts amplicon sequences by number of reads
			# Major seqs will be the firsts to be compared)
			my (@unique_seq_seqs, @unique_seq_md5s, @unique_seq_names);
			my @sorted_depth_amplicon_unique_seqs = sort { $amplicon_depths{$b} <=> $amplicon_depths{$a} } keys %amplicon_depths;
			foreach my $md5 (@sorted_depth_amplicon_unique_seqs) {
				push(@unique_seq_md5s, $md5);
				push(@unique_seq_seqs, $amplicon_sequences{$md5});
				push(@unique_seq_names, $amplicon_names{$md5});
				# push(@{$marker_names}, sprintf("%s-%d", $md5, $amplicon_depths{$md5}));
			}	
			# Variants will be the sequences not compared before against all the others
			my @variants_seqs = @unique_seq_seqs;
			my @variants_md5s = @unique_seq_md5s;
			my @variants_names = @unique_seq_names;

			# Stores highest frequency in the amplicon
# 			my $highest_amplicon_frequency = $amplicon_depths{$unique_seq_md5s[0]};

			# Stores comparisons, to avoid repeating the alignments
			my $comparisons;

			# Compares sequences
			for (my $i=0; $i<=$#unique_seq_md5s; $i++) {

				# Dominant sequence will be compared against all the others
				my $dominant_md5 = $unique_seq_md5s[$i];
				my $dominant_seq = $unique_seq_seqs[$i];
				my $dominant_len = length($dominant_seq);

# if ($dominant_md5 eq 'd0fb8aaf615d9688e6ee9802bd7aa842'){
# print '';
# }

				# Variables to store data from compared sequences
				my (@similar_seq_names, @similar_seq_md5s, @similar_seq_comparisons, @similar_seq_chimeras);

				# Aligns the reference sequence with all the other sequences 1 by 1, the alignments will be accurate
				my $ref_var_aligned_seqs = align_seqs2one($dominant_seq,$dominant_md5,\@variants_seqs,\@variants_md5s,'needleall');

				# Loops sequences not previously compared
				for (my $j=0; $j<=$#variants_md5s; $j++) {

					my $variant_seq = $variants_seqs[$j];
					my $variant_md5 = $variants_md5s[$j];
					my $variant_name = $variants_names[$j];
# if ($variant_md5 eq '0197e5f5edb65ab603b8169dfb796b7e'){
# print '';
# }
					# Skips the same sequence
					if ($dominant_md5 eq $variant_md5){
						next;
					}
					
					my %comparison;
					if (!defined($comparisons->{$variant_md5}{$dominant_md5})){
						# Skip if the pairwise global alignment fails
						if (!defined($ref_var_aligned_seqs->{$variant_md5}[0]) || !defined($ref_var_aligned_seqs->{$variant_md5}[1])){
							next;
						}
						my ($identical,$total) = binary_score_nts($ref_var_aligned_seqs->{$variant_md5}[0],$ref_var_aligned_seqs->{$variant_md5}[1]);
						my ($substitutions_, $insertions_, $deletions_, $homopolymer_indels_) = detect_sequence_errors($ref_var_aligned_seqs->{$variant_md5}[0],$ref_var_aligned_seqs->{$variant_md5}[1]);
						# Annotates comparison
						%comparison = (
							'identity' => scalar sprintf("%.2f", $identical/$dominant_len*100),
							'substitutions' => scalar @$substitutions_,
							'insertions' => scalar @$insertions_,
							'deletions' => scalar @$deletions_,
							'homopoymer_indels' => scalar @$homopolymer_indels_,
							'non_homopolymer_indels' => scalar @$insertions_ + scalar @$deletions_ - scalar @$homopolymer_indels_,
						);
						my %comparison_comp = %comparison;
						$comparison_comp{'insertions'} = scalar @$deletions_;
						$comparison_comp{'deletions'} = scalar @$insertions_;
						#$comparison = sprintf("ident: %d, sub: %d, ins: %d, del: %d, homo: %d", $identity_, scalar @$substitutions_, scalar @$insertions_, scalar @$deletions_, scalar @$homopolymer_indels_);
						$comparisons->{$dominant_md5}{$variant_md5} = \%comparison;
						$comparisons->{$variant_md5}{$dominant_md5} = \%comparison_comp;
					} else {
						%comparison = %{$comparisons->{$dominant_md5}{$variant_md5}};
					}
					#$comparison =~ /ident: (\d+), sub: (\d+), ins: (\d+), del: (\d+), homo: (\d+)/;
					#my ($identity, $substitutions, $insertions, $deletions, $homopolymer_indels) = ($1, $2, $3, $4, $5);
					# my $non_homopolymer_indels = $insertions + $deletions - $homopolymer_indels;
					# my $total = length($ref_var_aligned_seqs->{$variant_md5}[0]);
				
					# Decides is sequences are similar enough to be annotated
					my $annotate = 1;
					# Doesn't nnotate seqs with lower identity than threshold
					if ($annotate && defined($clustering_thresholds{'identity_threshold'}) && $comparison{'identity'} < $clustering_thresholds{'identity_threshold'}){
						$annotate = 0;
					}
					# Doesn't annotate seqs with more substitutions than threshold
					# Ej. 3 subs > (2%*130=2.6 || 3)
					# Ej. 2 subst > (1%*130=1.3 || 1)
					if ($annotate && defined($clustering_thresholds{'substitution_threshold'}) && $comparison{'substitutions'} > ($clustering_thresholds{'substitution_threshold'}/100*$dominant_len || sprintf("%.0f", $clustering_thresholds{'substitution_threshold'}/100*$dominant_len)) ){
						$annotate = 0;
					}
					# Doesn't annotate seqs with more non homopolymer indels than threshold
					if ($annotate && defined($clustering_thresholds{'indel_threshold'}) && $comparison{'non_homopolymer_indels'} > ($clustering_thresholds{'indel_threshold'}/100*$dominant_len || sprintf("%.0f", $clustering_thresholds{'indel_threshold'}/100*$dominant_len)) ){
						$annotate = 0;
					}

					# Annotates comparison data if the sequence is similar enough
					if ($annotate) {
						push(@similar_seq_names, $variant_name);
						push(@similar_seq_md5s, $variant_md5);
						my $comparison = '';
						if ($comparison{'substitutions'}>0){
							$comparison .= sprintf("X%d", $comparison{'substitutions'});
						}
						if ($comparison{'insertions'}>0){
							$comparison .= sprintf("I%d", $comparison{'insertions'});
						}
						if ($comparison{'deletions'}>0){
							$comparison .= sprintf("D%d", $comparison{'deletions'});
						}
						if ($comparison{'homopoymer_indels'}>0){
							$comparison .= sprintf("H%d", $comparison{'homopoymer_indels'});
						}
						push(@similar_seq_comparisons, $comparison);
					}
# 					# Finishes when there are 2 annotations
# 					if ($#similar_seq_md5s == 1) {
# 						 last;
# 					}
				}

				# Annotates seqs that are chimeras from other sequences with higher depth
				if ($i>1) {
					my ($is_chimera,$i_seq1,$i_seq2) = is_chimera($dominant_seq,[@unique_seq_seqs[0 .. $i-1]]);
					if ($is_chimera){
						push(@similar_seq_chimeras, sprintf("%s+%s",$unique_seq_names[$i_seq1],$unique_seq_names[$i_seq2]));
						# print "CHIMERA ".sprintf("%s+%s",$unique_seq_names[$i_seq1],$unique_seq_names[$i_seq2])."\n";
					}
				}

				# Annotates comparison data of all similar sequences into the amplicon
				$amplicon_seq_comparison_data->{$marker_name}{$sample_name}{$dominant_md5}{'depth'} = $amplicon_depths{$dominant_md5};
				$amplicon_seq_comparison_data->{$marker_name}{$sample_name}{$dominant_md5}{'freq'} = $amplicon_frequencies{$dominant_md5};
				$amplicon_seq_comparison_data->{$marker_name}{$sample_name}{$dominant_md5}{'similar_seq_comparisons'} = \@similar_seq_comparisons;
				$amplicon_seq_comparison_data->{$marker_name}{$sample_name}{$dominant_md5}{'similar_seq_md5s'} = \@similar_seq_md5s;
				$amplicon_seq_comparison_data->{$marker_name}{$sample_name}{$dominant_md5}{'similar_seq_names'} = \@similar_seq_names;
				$amplicon_seq_comparison_data->{$marker_name}{$sample_name}{$dominant_md5}{'similar_seq_chimeras'} = \@similar_seq_chimeras;

			}

			if ($sample_name ne ''){
				printf("\t%s-%s pairwise compared %d sequences\n", $marker_name, $sample_name, scalar @unique_seq_md5s);
			} else {
				printf("\t%s pairwise compared %d sequences\n", $marker_name, scalar @unique_seq_md5s);
			}

		}

	}

	return $amplicon_seq_comparison_data;

}

#################################################################################

# Compares high depth unique sequences among them
sub compare_amplicon_sequences_with_threads {

	my ($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$threads_limit) = @_;

	if (!defined($threads_limit)){
		$threads_limit = 4;
	}
	
	my $amplicon_seq_comparison_data;

	my $one_amplicon_seq_data;
	foreach my $marker_name (@$markers){
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}
		my $clustered_amplicon_depths->{$marker_name} = {};
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
				next;
			}
			my $one_amplicon_seq_data_;
			my $clustered_amplicon_sequences->{$marker_name}{$sample_name} = {};
			$one_amplicon_seq_data_->{$marker_name}{$sample_name} = $amplicon_seq_data->{$marker_name}{$sample_name};
			push(@{$one_amplicon_seq_data},$one_amplicon_seq_data_);
			# print '';
		}
	}

	my @threads;
	for (my $count_amplicon=0; $count_amplicon<=$#{$one_amplicon_seq_data}; $count_amplicon++){

		push(@threads, threads->create(\&compare_amplicon_sequences,$markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$one_amplicon_seq_data->[$count_amplicon]));
# # 		push(@threads, 1);
# # 		print "\n";

# 		# For debugging:
# 		push(@threads, [compare_amplicon_sequences($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$one_amplicon_seq_data->[$count_amplicon])]);

# 		# If maximum number of threads is reached or last sbjct of a query is processed
		if (scalar @threads >= $threads_limit  || $count_amplicon == $#{$one_amplicon_seq_data}){
# 			print scalar @threads."\n"; exit;
			my $check_threads = 1;
			while ($check_threads){
				for (my $i=0; $i<=$#threads; $i++){
					unless ($threads[$i]->is_running()){
						my $amplicon_seq_comparison_data_ = $threads[$i]->join;
						if (defined($amplicon_seq_comparison_data_)){
							my $marker_name = (keys %$amplicon_seq_comparison_data_)[0];
							my $sample_name = (keys %{$amplicon_seq_comparison_data_->{$marker_name}})[0];
	# 						print "\t$marker_name-$sample_name finished\n";
							if (defined($amplicon_seq_comparison_data_->{$marker_name}{$sample_name})){
								$amplicon_seq_comparison_data->{$marker_name}{$sample_name} = $amplicon_seq_comparison_data_->{$marker_name}{$sample_name};
							}
						}
						undef($threads[$i]);
						splice(@threads,$i,1);
						$i = $i - 1;
						unless ($count_amplicon == $#{$one_amplicon_seq_data} && @threads){
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
# 				print '';
# 				my $amplicon_seq_comparison_data_ = $threads[$i];
# 				my $marker_name = (keys %$amplicon_seq_comparison_data_)[0];
# # 				my $sample_name = (keys %{$amplicon_seq_comparison_data_->{$marker_name}})[0];
# 				$amplicon_seq_comparison_data->{$marker_name}{$sample_name} = $amplicon_seq_comparison_data_->{$marker_name}{$sample_name};
# 				delete $threads[$i];
# 			}
		}
	}
	
	return $amplicon_seq_comparison_data;
}


#################################################################################

# Genotypes amplicon sequences with the desired method ('Sommer', 'Lighten' or 'Radwan')
sub genotype_amplicon_sequences {

	my ($method,$markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data) = @_;

	my ($genotyped_amplicon_sequences,$genotyped_amplicon_depths,$genotyping_output);
	
	if ($method eq 'radwan'){
	
		# Herdegen genotyping method implementation
		# $all_seqs = $radwan_alleles + $radwan_artifacts
		# my $algorithm_parameters->{'error_threshold'} = 2; # Artefacts are max. 2 substitutions from a more abundant variant
		# my $algorithm_parameters->{'min_amplicon_seq_frequency'} = 3; # All variants <3% freq. are rejected as artifacts
		# my $algorithm_parameters->{'max_amplicon_seq_frequency'} = 12; # All variants >12% freq. are accepted as true alleles
	
		foreach my $marker_name (@$markers){

			# Skips markers without data
			if (!defined($amplicon_seq_data->{$marker_name})){
				next;
			}

			# Excludes amplicons not in the list
			if (defined($paramsdata->{'allowed_markers'}) && !defined($paramsdata->{'allowed_markers'}{'all'}) && !in_array($paramsdata->{'allowed_markers'},$marker_name)){
				next;
			}

			# Extracts genotyping parameters
			my $genotyping_parameters;
			foreach my $param (keys %$paramsdata) {
				if (defined($paramsdata->{$param}{$marker_name})){
					$genotyping_parameters->{$param} = $paramsdata->{$param}{$marker_name}[0];
				} elsif (defined($paramsdata->{$param}{'all'})){
					$genotyping_parameters->{$param} = $paramsdata->{$param}{'all'}[0];
				}
			}

			my ($md5_to_sequence, $md5_to_name);

			# Loops samples/amplicons and performs genotyping
			foreach my $sample_name (@{$samples->{$marker_name}}) {

				# Process only samples with sequences in the original order
				if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
	# 				print "\t$marker_name-$sample_name doesn't have sequences to compare.\n";
					next;
				}
				
				# Exclude samples not in the list
				if (defined($genotyping_parameters->{'allowed_samples'}) && ref($genotyping_parameters->{'allowed_samples'}) eq 'ARRAY' && !in_array($genotyping_parameters->{'allowed_samples'},$sample_name)){
					next;
				}

				if ($sample_name ne ''){
					print "\t$marker_name-$sample_name genotyping\n";
				} else {
					print "\t$marker_name genotyping\n";
				}

				my @sorted_md5s =  sort { $amplicon_seq_data->{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->{$marker_name}{$sample_name}{$a}{'depth'} }  keys %{$amplicon_seq_data->{$marker_name}{$sample_name}};
				my (@previous_seqs, @previous_md5s);
				foreach my $md5 (@sorted_md5s){
					my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
					my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
					my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
					my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
					my $frequency = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'};
					my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
					my $mean_freq = $marker_seq_data->{$marker_name}{$md5}{'mean_freq'};
					my $min_freq = $marker_seq_data->{$marker_name}{$md5}{'min_freq'};
					my $max_freq = $marker_seq_data->{$marker_name}{$md5}{'max_freq'};
					my $header = sprintf("hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq);
					my @genotyp_criteria;
					my $is_allele = 1;
					# First checks if sequence has correct length
					if (!is_dominant($seq,$markerdata->{$marker_name}{'length'},$genotyping_parameters)){
						unshift(@genotyp_criteria, sprintf("wrong length (%d)", $len));
						$is_allele = 0;
					}
					if ($frequency < $genotyping_parameters->{'min_amplicon_seq_frequency'}){
						unshift(@genotyp_criteria, sprintf("low freq (<%.2f%%)", $genotyping_parameters->{'min_amplicon_seq_frequency'}));
						$is_allele = 0;
					} elsif (@previous_seqs && $frequency >= $genotyping_parameters->{'min_amplicon_seq_frequency'} && $frequency <= $genotyping_parameters->{'max_amplicon_seq_frequency'}){ # Grey zone sequences
						my $aligned_seqs = align_seqs2one($seq,$md5,\@previous_seqs,\@previous_md5s,'needleall');
						foreach my $previous_md5 (@previous_md5s) {
							# Skips if the pairwise global alignment fails
							if (!defined($aligned_seqs->{$previous_md5}[0]) || !defined($aligned_seqs->{$previous_md5}[1])){
								next;
							}
							my ($identical,$total) = binary_score_nts($aligned_seqs->{$previous_md5}[0],$aligned_seqs->{$previous_md5}[1]);
							# If any previous sequence is N mismatches
							if ($total-$identical<=$genotyping_parameters->{'error_threshold'} && $depth<$genotyping_parameters->{'min_dominant_frequency_threshold'}*$amplicon_seq_data->{$marker_name}{$sample_name}{$previous_md5}{'depth'}/100) {
								#my ($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($aligned_seqs->{$previous_md5}[0],$aligned_seqs->{$previous_md5}[1]);
								unshift(@genotyp_criteria, sprintf("1-%dbp diff (%s)", $genotyping_parameters->{'error_threshold'}, $md5_to_name->{$previous_md5}));
								#push(@genotyp_criteria, sprintf("error: %s (%s)", $md5_to_name->{$previous_md5}, print_sequence_errors($substitutions, $insertions, $deletions, $homopolymer_indels)));
								$is_allele = 0;
								last;
							}
						}
						my ($is_chimera, $i_seq1, $i_seq2) = is_chimera($seq,\@previous_seqs, undef, $genotyping_parameters->{'error_threshold'});
						if ($is_chimera) {
							unshift(@genotyp_criteria, sprintf("chimera (%s+%s)",$md5_to_name->{$previous_md5s[$i_seq1]},$md5_to_name->{$previous_md5s[$i_seq2]}));
							$is_allele = 0;
							# push(@genotyp_criteria, sprintf("chimera: %s+%s",$md5_to_name->{$previous_md5s[$i_seq1]},$md5_to_name->{$previous_md5s[$i_seq2]}));
						}
						unshift(@genotyp_criteria, sprintf("grey zone (%.2f%%<=freq>=%.2f%%)", $genotyping_parameters->{'min_amplicon_seq_frequency'}, $genotyping_parameters->{'max_amplicon_seq_frequency'}));
					} else {
						unshift(@genotyp_criteria, sprintf("high freq (>%.2f%%)", $genotyping_parameters->{'max_amplicon_seq_frequency'}));
					}
					if ($is_allele){
						$genotyped_amplicon_sequences->{$marker_name}{$sample_name}{$md5} = $depth;
						$genotyped_amplicon_depths->{$marker_name}{$sample_name} += $depth;
						$genotyping_output->{$marker_name}{$sample_name} .= sprintf(">*%s | Allele: %s | %s\n%s\n", $name, join(', ', @genotyp_criteria), $header, $seq);
					} else {
						$genotyping_output->{$marker_name}{$sample_name} .= sprintf(">#%s | Artifact: %s | %s\n%s\n", $name, join(', ', @genotyp_criteria), $header, $seq);
					}
					push(@previous_seqs,$seq);
					push(@previous_md5s,$md5);
					if (!defined($md5_to_sequence->{$md5})){
						$md5_to_sequence->{$md5} = $seq;
						$md5_to_name->{$md5} = $name;
					}
				}
				if ($sample_name ne ''){
					printf("\t%s-%s genotyped (%d reads, %d variants)\n", $marker_name, $sample_name, $genotyped_amplicon_depths->{$marker_name}{$sample_name}, scalar keys %{$genotyped_amplicon_sequences->{$marker_name}{$sample_name}});
				} else {
					printf("\t%s genotyped (%d reads, %d variants)\n", $marker_name, $genotyped_amplicon_depths->{$marker_name}{$sample_name}, scalar keys %{$genotyped_amplicon_sequences->{$marker_name}{$sample_name}});
				}
			}
		}
	# End Herdegen method

	} elsif ($method eq 'lighten'){
	
		# Lighten genotyping method implementation
		# $all_seqs = $lighten_variants + $lighten_rpes
		# $lighten_variants = $lighten_alleles + $lighten_artifacts
		# my $algorithm_parameters->{'error_threshold'} = 3; # Errors are 1-3bp mismatches from parental PAs
		# my $algorithm_parameters->{'min_dominant_frequency_threshold'} = 2; # RPEs are <2% depth respect to parental PAs
		# my $algorithm_parameters->{'max_allele_number'} = 10; # Maximum number of expected alleles to calculate DOCs

		foreach my $marker_name (@$markers){

			# Skips markers without data
			if (!defined($amplicon_seq_data->{$marker_name})){
				next;
			}

			# Excludes amplicons not in the list
			if (defined($paramsdata->{'allowed_markers'}) && !defined($paramsdata->{'allowed_markers'}{'all'}) && !in_array($paramsdata->{'allowed_markers'},$marker_name)){
				next;
			}

			# Extracts genotyping parameters
			my $genotyping_parameters;
			foreach my $param (keys %$paramsdata) {
				if (defined($paramsdata->{$param}{$marker_name})){
					$genotyping_parameters->{$param} = $paramsdata->{$param}{$marker_name}[0];
				} elsif (defined($paramsdata->{$param}{'all'})){
					$genotyping_parameters->{$param} = $paramsdata->{$param}{'all'}[0];
				}
			}

			my ($md5_to_sequence, $md5_to_name);

			# Loops samples/amplicons and performs genotyping
			foreach my $sample_name (@{$samples->{$marker_name}}) {

				# Process only samples with sequences in the original order
				if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
	# 				print "\t$marker_name-$sample_name doesn't have sequences to compare.\n";
					next;
				}
				
				# Exclude samples not in the list
				if (defined($genotyping_parameters->{'allowed_samples'}) && ref($genotyping_parameters->{'allowed_samples'}) eq 'ARRAY' && !in_array($genotyping_parameters->{'allowed_samples'},$sample_name)){
					next;
				}

				print "\t$marker_name-$sample_name genotyping\n";

				# Saves new depths after clustering
				my (%amplicon_clustered_depths, %amplicon_clustered_artifacts, %amplicon_wrong_length) ;
				my @sorted_md5s_ =  sort { $amplicon_seq_data->{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->{$marker_name}{$sample_name}{$a}{'depth'} }  keys %{$amplicon_seq_data->{$marker_name}{$sample_name}};
				my @sorted_seqs_ = map $marker_seq_data->{$marker_name}{$_}{'seq'}, @sorted_md5s_;
				while (@sorted_md5s_){
					my $md5 = shift @sorted_md5s_;
					my $seq = shift @sorted_seqs_;
					my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
					# Do not start cluster sequences with incorrect length
					if (!is_dominant($seq,$markerdata->{$marker_name}{'length'},$genotyping_parameters)){
						$amplicon_wrong_length{$md5} = 1;
						next;
					}
					$amplicon_clustered_depths{$md5} = $depth;
					if (!@sorted_md5s_){
						last;
					}
					my $aligned_seqs = align_seqs2one($seq,$md5,\@sorted_seqs_,\@sorted_md5s_,'needleall');
					for (my $j=0; $j<=$#sorted_md5s_; $j++){
						my $md5_ = $sorted_md5s_[$j];
						my $seq_ = $sorted_seqs_[$j];
						my $depth_ = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5_}{'depth'};
						if ($depth_ < $depth*$genotyping_parameters->{'min_dominant_frequency_threshold'}/100 ){
							# Skips if the pairwise global alignment fails
							if (!defined($aligned_seqs->{$md5_}[0]) || !defined($aligned_seqs->{$md5_}[1])){
								next;
							}
							my ($identical,$total) = binary_score_nts($aligned_seqs->{$md5_}[0],$aligned_seqs->{$md5_}[1]);
							if ($total-$identical <= $genotyping_parameters->{'error_threshold'}){
								$amplicon_clustered_depths{$md5} += $depth_;
								$amplicon_clustered_artifacts{$md5_} = $md5;
								splice(@sorted_md5s_,$j,1);
								splice(@sorted_seqs_,$j,1);
								$j--;
							}
						}
					}
					if (scalar keys %amplicon_clustered_depths > $genotyping_parameters->{'max_allele_number'}){
						last;
					}
				}
				my @sorted_clustered_md5s =  sort { $amplicon_clustered_depths{$b} <=> $amplicon_clustered_depths{$a} } keys %amplicon_clustered_depths;
				my @sorted_clustered_depths = map $amplicon_clustered_depths{$_}, @sorted_clustered_md5s;
				my ($doc_allele_number,$DOCns) = degree_of_change(\@sorted_clustered_depths,$genotyping_parameters->{'max_allele_number'});
				my %amplicon_clustered_alleles;
				map $amplicon_clustered_alleles{$sorted_clustered_md5s[$_]} = $_+1 , 0..$#{$DOCns};

				#my @sorted_md5s =  sort { $amplicon_seq_data->{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->{$marker_name}{$sample_name}{$a}{'depth'} }  keys %{$amplicon_seq_data->{$marker_name}{$sample_name}};
				foreach my $md5 (@sorted_clustered_md5s){
					my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
					my $name = $marker_seq_data->{$marker_name}{$md5}{'name'};
					my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
					my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
					my $frequency = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'};
					my $count_samples = $marker_seq_data->{$marker_name}{$md5}{'samples'};
					my $mean_freq = $marker_seq_data->{$marker_name}{$md5}{'mean_freq'};
					my $min_freq = $marker_seq_data->{$marker_name}{$md5}{'min_freq'};
					my $max_freq = $marker_seq_data->{$marker_name}{$md5}{'max_freq'};
					my $header = sprintf("hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq);
					my @genotyp_criteria;
					my $is_allele = 0;
					# First checks if sequence has correct length
					if (defined($amplicon_wrong_length{$md5})){
						unshift(@genotyp_criteria, sprintf("wrong length (%d)", $len));
					} elsif (defined($amplicon_clustered_alleles{$md5})){
						if ($amplicon_clustered_alleles{$md5}<=$doc_allele_number){
							unshift(@genotyp_criteria, sprintf("DOC=%.2f%%", $DOCns->[$amplicon_clustered_alleles{$md5}-1]));
							$is_allele = 1;
						} else {
							unshift(@genotyp_criteria, sprintf("DOC=%.2f%%", $DOCns->[$amplicon_clustered_alleles{$md5}-1]));
						}
					} elsif (defined($amplicon_clustered_artifacts{$md5})){
						unshift(@genotyp_criteria, sprintf("1-%dbp diff (%s)", $genotyping_parameters->{'error_threshold'}, $md5_to_name->{$amplicon_clustered_artifacts{$md5}}));
					} else {
						unshift(@genotyp_criteria, sprintf("maximum allele number reached (%d)", $genotyping_parameters->{'max_allele_number'}));
					}
					if (defined($amplicon_clustered_depths{$md5}) && $amplicon_clustered_depths{$md5} > $depth){
						unshift(@genotyp_criteria, sprintf("clustered_depth=%d", $amplicon_clustered_depths{$md5}));
					}
					if ($is_allele){
						$genotyped_amplicon_sequences->{$marker_name}{$sample_name}{$md5} = $amplicon_clustered_depths{$md5};
						$genotyped_amplicon_depths->{$marker_name}{$sample_name} += $amplicon_clustered_depths{$md5};
						$genotyping_output->{$marker_name}{$sample_name} .= sprintf(">*%s | Allele: %s | %s\n%s\n", $name, join(', ', @genotyp_criteria), $header, $seq);
					} else {
						$genotyping_output->{$marker_name}{$sample_name} .= sprintf(">#%s | Artifact: %s | %s\n%s\n", $name, join(', ', @genotyp_criteria), $header, $seq);
					}
					if (!defined($md5_to_sequence->{$md5})){
						$md5_to_sequence->{$md5} = $seq;
						$md5_to_name->{$md5} = $name;
					}
				}
				if ($sample_name ne ''){
					printf("\t%s-%s genotyped (%d reads, %d variants)\n", $marker_name, $sample_name, $genotyped_amplicon_depths->{$marker_name}{$sample_name}, scalar keys %{$genotyped_amplicon_sequences->{$marker_name}{$sample_name}});
				} else {
					printf("\t%s genotyped (%d reads, %d variants)\n", $marker_name, $genotyped_amplicon_depths->{$marker_name}{$sample_name}, scalar keys %{$genotyped_amplicon_sequences->{$marker_name}{$sample_name}});
				}
			}
		}
	# End Lighten method

	} elsif ($method eq 'sommer'){
	
		# Sommer genotyping method implementation
		# $all_seqs = $sommer_alleles + $sommer_artifacts + $sommer_chimeras + $sommer_unclassified
		# my $algorithm_parameters->{'error_threshold'} = 2; # Artifacts are 1-2bp diff from putative alleles
		
		# Only in Sommer method, $amplicon_seq_data and $marker_seq_data will be two array references containing the data of the different experiments or duplicates
		# my $amplicon_seq_data = [ $amplicon_seq_data_rep0, $amplicon_seq_data_rep1 ... $amplicon_seq_data_repN ];
		# my $marker_seq_data = [ $marker_seq_data_rep0, $marker_seq_data_rep1 ... $marker_seq_data_repN ];

		foreach my $marker_name (@$markers){

			# Skips markers without data
			if (!defined($amplicon_seq_data->[0]{$marker_name})){
				next;
			}

			# Excludes amplicons not in the list
			if (defined($paramsdata->{'allowed_markers'}) && !defined($paramsdata->{'allowed_markers'}{'all'}) && !in_array($paramsdata->{'allowed_markers'},$marker_name)){
				next;
			}

			# Extracts genotyping parameters
			my $genotyping_parameters;
			foreach my $param (keys %$paramsdata) {
				if (defined($paramsdata->{$param}{$marker_name})){
					$genotyping_parameters->{$param} = $paramsdata->{$param}{$marker_name}[0];
				} elsif (defined($paramsdata->{$param}{'all'})){
					$genotyping_parameters->{$param} = $paramsdata->{$param}{'all'}[0];
				}
			}

			my ($md5_to_sequence, $md5_to_name);

			my ($sommer_alleles, $sommer_artifacts, $sommer_chimeras, $sommer_less2bps, $sommer_more2bps, $sommer_more2bps_2, $sommer_unclassified, $sommer_wrong_length, $sommer_low_freq);

			# Loops samples/amplicons and performs genotyping
			foreach my $sample_name (@{$samples->{$marker_name}}) {
				# Process only samples with sequences in the original order
				if (!defined($amplicon_seq_data->[0]{$marker_name}{$sample_name})){
					next;
				}
				# Exclude samples not in the list
				if (defined($genotyping_parameters->{'allowed_samples'}) && ref($genotyping_parameters->{'allowed_samples'}) eq 'ARRAY' && !in_array($genotyping_parameters->{'allowed_samples'},$sample_name)){
					next;
				}
				print "\t$marker_name-$sample_name Step I genotyping.\n";
				for (my $i=0; $i<=$#{$amplicon_seq_data}; $i++){
					my @sorted_md5s =  sort { $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$a}{'depth'} }  keys %{$amplicon_seq_data->[$i]{$marker_name}{$sample_name}};
					my (@previous_seqs, @previous_md5s);
					foreach my $md5 (@sorted_md5s) {
# if ($md5 eq '5f47972bc2b7db84e0a9a850f962e99c'){
# print '';
# }
						my $seq = $marker_seq_data->[$i]{$marker_name}{$md5}{'seq'};
						my $name = $marker_seq_data->[$i]{$marker_name}{$md5}{'name'};
						my $len = $marker_seq_data->[$i]{$marker_name}{$md5}{'len'};
						my $depth = $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$md5}{'depth'};
						my $frequency = $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$md5}{'freq'};
						# First checks if sequence has correct length
						if (!is_dominant($seq,$markerdata->{$marker_name}{'length'},$genotyping_parameters)){
							$sommer_wrong_length->[$i]{$sample_name}{$md5} = sprintf("Artifact: Step 0, wrong length (%d)", $len);
							next;
						} elsif ($depth<=1){
							$sommer_artifacts->[$i]{$sample_name}{$md5} = sprintf("Artifact: Step I, singleton");
							next;
						# Frequency threshold to exclude low frequency variants (NOT IN ORIGINAL ARTICLE)
						} elsif (defined($genotyping_parameters->{'min_amplicon_seq_frequency'}) && $frequency<$genotyping_parameters->{'min_amplicon_seq_frequency'}){
							$sommer_artifacts->[$i]{$sample_name}{$md5} = sprintf("Artifact: Step I, low freq (<%.2f%%)", $genotyping_parameters->{'min_amplicon_seq_frequency'});
							next;
						} elsif (!defined($sommer_alleles->[$i]{$sample_name})){
							$sommer_alleles->[$i]{$sample_name}{$md5} = sprintf("Allele: Step I, first cluster");
						} else {
							my ($is_chimera, $i_seq1, $i_seq2) = is_chimera($seq,\@previous_seqs,undef,$genotyping_parameters->{'error_threshold'});
							if ($is_chimera) {
								$sommer_chimeras->[$i]{$sample_name}{$md5} = sprintf("%s+%s",$md5_to_name->{$previous_md5s[$i_seq1]},$md5_to_name->{$previous_md5s[$i_seq2]});
							} else {
								my $aligned_seqs = align_seqs2one($seq,$md5,\@previous_seqs,\@previous_md5s,'needleall');
								my $is_less2bp = 0;
								foreach my $previous_md5 (@previous_md5s) {
									# Skips if the pairwise global alignment fails
									if (!defined($aligned_seqs->{$previous_md5}[0]) || !defined($aligned_seqs->{$previous_md5}[1])){
										next;
									}
									my ($identical,$total) = binary_score_nts($aligned_seqs->{$previous_md5}[0],$aligned_seqs->{$previous_md5}[1]);
									# If any previous sequence is 1 or 2 mismatches, can be changed
									if ($total-$identical<=$genotyping_parameters->{'error_threshold'}) {
										$sommer_less2bps->[$i]{$sample_name}{$md5} = sprintf("%s", $md5_to_name->{$previous_md5});
										$is_less2bp = 1;
										last;
									}
								}
								if (!$is_less2bp) {
									$sommer_more2bps->[$i]{$sample_name}{$md5} = $depth;
								}
							}
						}
						push(@previous_seqs,$seq);
						push(@previous_md5s,$md5);
						if (!defined($md5_to_sequence->{$md5})){
							$md5_to_sequence->{$md5} = $seq;
							$md5_to_name->{$md5} = $name;
						}
					}
				}
			}
			# Sommer method STEP II
			foreach my $sample_name (@{$samples->{$marker_name}}){
				# Process only samples with sequences in the original order
				if (!defined($amplicon_seq_data->[0]{$marker_name}{$sample_name})){
					next;
				}
				# Exclude samples not in the list
				if (defined($genotyping_parameters->{'allowed_samples'}) && ref($genotyping_parameters->{'allowed_samples'}) eq 'ARRAY' && !in_array($genotyping_parameters->{'allowed_samples'},$sample_name)){
					next;
				}
				print "\t$marker_name-$sample_name Step II genotyping.\n";
				for (my $i=0; $i<=$#{$amplicon_seq_data}; $i++){
					my @sorted_md5s =  sort { $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$a}{'depth'} }  keys %{$amplicon_seq_data->[$i]{$marker_name}{$sample_name}};
					my (@previous_seqs, @previous_md5s);
					foreach my $md5 (@sorted_md5s) {
						if (defined($sommer_alleles->[$i]{$sample_name}{$md5}) || defined($sommer_artifacts->[$i]{$sample_name}{$md5}) || defined($sommer_wrong_length->[$i]{$sample_name}{$md5})){
							next;
						}
						my $seq = $marker_seq_data->[$i]{$marker_name}{$md5}{'seq'};
						my $name = $marker_seq_data->[$i]{$marker_name}{$md5}{'name'};
						#my $len = $marker_seq_data->[$i]{$marker_name}{$md5}{'len'};
						my $depth = $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$md5}{'depth'};
						my $has_replicates = 1;
						if (!defined($amplicon_seq_data->[$i-1]{$marker_name}{$sample_name}{$md5})){
							$has_replicates = 0;
						}
						if (defined($sommer_less2bps->[$i]{$sample_name}{$md5})){
							if (!$has_replicates){
								$sommer_artifacts->[$i]{$sample_name}{$md5} = sprintf("Artifact: Step II, 1-%dbp diff (%s), no replicate",$genotyping_parameters->{'error_threshold'},$sommer_less2bps->[$i]{$sample_name}{$md5});
								delete($sommer_less2bps->[$i]{$sample_name}{$md5});
							}
							next;
						}
						if (defined($sommer_more2bps->[$i]{$sample_name}{$md5})){
							if (!$has_replicates) {
								my $another_sample = 0;
								foreach my $sample_name_ (@{$samples->{$marker_name}}){
									if ($sample_name_ ne $sample_name && defined($amplicon_seq_data->[$i]{$marker_name}{$sample_name_}{$md5})){
										$another_sample = 1;
										last;
									}
									if ($another_sample) { last; }
								}
								if (!$another_sample){
									$sommer_artifacts->[$i]{$sample_name}{$md5} = sprintf("Artifact: Step II, >%dbp diff, no replicate, no other individuals",$genotyping_parameters->{'error_threshold'});
								} else {
									$sommer_more2bps_2->[$i]{$sample_name}{$md5} = $depth;
								}
								delete($sommer_more2bps->[$i]{$sample_name}{$md5});
							}
							next;
						}
						if (defined($sommer_chimeras->[$i]{$sample_name}{$md5})){
							if (!$has_replicates){
								$sommer_artifacts->[$i]{$sample_name}{$md5} = sprintf("Artifact:  Step II, chimera (%s), no replicate",$sommer_chimeras->[$i]{$sample_name}{$md5});
								delete($sommer_chimeras->[$i]{$sample_name}{$md5});
							} else {
								my $has_chimera_replicates = 0;
								if (defined($sommer_chimeras->[$i-1]{$sample_name}{$md5})){
									$has_chimera_replicates = 1;	
									$sommer_artifacts->[$i]{$sample_name}{$md5} = sprintf("Artifact: Step II, chimera in replicate (%s)",$sommer_chimeras->[$i]{$sample_name}{$md5});
									delete($sommer_chimeras->[$i]{$sample_name}{$md5});
									$sommer_artifacts->[$i-1]{$sample_name}{$md5} = sprintf("Artifact: Step II, chimera in replicate (%s)",$sommer_chimeras->[$i-1]{$sample_name}{$md5});
									delete($sommer_chimeras->[$i-1]{$sample_name}{$md5});
								}
							}
							next;
						}
					}
				}
			}
			# Sommer method STEP III
			foreach my $sample_name (@{$samples->{$marker_name}}){
				# Process only samples with sequences in the original order
				if (!defined($amplicon_seq_data->[0]{$marker_name}{$sample_name})){
					next;
				}
				# Exclude samples not in the list
				if (defined($genotyping_parameters->{'allowed_samples'}) && ref($genotyping_parameters->{'allowed_samples'}) eq 'ARRAY' && !in_array($genotyping_parameters->{'allowed_samples'},$sample_name)){
					next;
				}
				print "\t$marker_name-$sample_name Step III genotyping.\n";
				for (my $i=0; $i<=$#{$amplicon_seq_data}; $i++){
					my @sommer_artifacts_step2;
					if (defined($sommer_artifacts->[$i]{$sample_name})){
						@sommer_artifacts_step2 = keys %{$sommer_artifacts->[$i]{$sample_name}};
					} else { # If there are not artifacts (ex. clustered data), includes an artificial singleton
						@sommer_artifacts_step2 = ('00000000000000000000000000000000');
						$amplicon_seq_data->[$i]{$marker_name}{$sample_name}{'00000000000000000000000000000000'}{'depth'} = 1;
					}
					my @sommer_artifacts_step2_freqs = map $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$_}{'freq'}, @sommer_artifacts_step2;
					my $max_artifacts_step2_freq = max(@sommer_artifacts_step2_freqs);
					my @sorted_md5s =  sort { $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$a}{'depth'} }  keys %{$amplicon_seq_data->[$i]{$marker_name}{$sample_name}};
					my (@previous_seqs, @previous_md5s);
					foreach my $md5 (@sorted_md5s) {
						if (defined($sommer_alleles->[$i]{$sample_name}{$md5}) || defined($sommer_artifacts->[$i]{$sample_name}{$md5}) || defined($sommer_wrong_length->[$i]{$sample_name}{$md5})){
							next;
						}
						my $seq = $marker_seq_data->[$i]{$marker_name}{$md5}{'seq'};
						my $name = $marker_seq_data->[$i]{$marker_name}{$md5}{'name'};
						#my $len = $marker_seq_data->[$i]{$marker_name}{$md5}{'len'};
						my $frequency = $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$md5}{'freq'};
						my $is_low_freq = 0;
						if ($frequency<$max_artifacts_step2_freq){
							$is_low_freq = 1;
						}
						if (defined($sommer_less2bps->[$i]{$sample_name}{$md5})){
							if (!$is_low_freq){
								$sommer_alleles->[$i]{$sample_name}{$md5} = sprintf("Allele: Step III, 1-%dbp diff (%s), no low freq (>=%.2f%%)",$genotyping_parameters->{'error_threshold'},$sommer_less2bps->[$i]{$sample_name}{$md5},$max_artifacts_step2_freq);
							} else {
								$sommer_unclassified->[$i]{$sample_name}{$md5} = sprintf("Unclassified: Step III, 1-%dbp diff (%s), low freq (<%.2f%%)",$genotyping_parameters->{'error_threshold'},$sommer_less2bps->[$i]{$sample_name}{$md5},$max_artifacts_step2_freq);
							}
							delete($sommer_less2bps->[$i]{$sample_name}{$md5});
							next;
						}
						if (defined($sommer_more2bps->[$i]{$sample_name}{$md5})){
							if (!$is_low_freq){
								$sommer_alleles->[$i]{$sample_name}{$md5} = sprintf("Allele: Step III, >%dbp diff, no low freq (>=%.2f%%)",$genotyping_parameters->{'error_threshold'},$max_artifacts_step2_freq);
							} else {
								$sommer_unclassified->[$i]{$sample_name}{$md5} = sprintf("Unclassified: Step III, >%dbp diff, low freq (<%.2f%%)",$genotyping_parameters->{'error_threshold'},$max_artifacts_step2_freq);
							}
							delete($sommer_more2bps->[$i]{$sample_name}{$md5});
							next;
						}
					}
					foreach my $md5 (keys %{$sommer_more2bps_2->[$i]{$sample_name}}){
						my $depth = $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$md5}{'depth'};
						my $another_sample;
						foreach my $sample_name_ (@{$samples->{$marker_name}}){
							if ($sample_name_ ne $sample_name && defined($sommer_alleles->[$i]{$sample_name_}{$md5})){ 
								$another_sample = sprintf("allele in individual %s",$sample_name_);
								last;
							} elsif ($sample_name_ ne $sample_name && defined($sommer_unclassified->[$i]{$sample_name_}{$md5})){
								$another_sample = sprintf("unclassified in individual %s",$sample_name_);
								last;
							}
							if (defined($another_sample)) { last; }
						}
						if (defined($another_sample)){
							$sommer_unclassified->[$i]{$sample_name}{$md5} = sprintf("Unclassified: Step III, >%dbp diff, %s",$genotyping_parameters->{'error_threshold'},$another_sample);
						} else {
							$sommer_artifacts->[$i]{$sample_name}{$md5} = sprintf("Artifact: Step III, >%dbp diff, no allele or unclassified in other individuals",$genotyping_parameters->{'error_threshold'});
						}
						delete($sommer_more2bps_2->[$i]{$sample_name}{$md5});
					}
				}
			}
			foreach my $sample_name (@{$samples->{$marker_name}}){
				# Process only samples with sequences in the original order
				if (!defined($amplicon_seq_data->[0]{$marker_name}{$sample_name})){
					next;
				}
				# Exclude samples not in the list
				if (defined($genotyping_parameters->{'allowed_samples'}) && ref($genotyping_parameters->{'allowed_samples'}) eq 'ARRAY' && !in_array($genotyping_parameters->{'allowed_samples'},$sample_name)){
					next;
				}
				for (my $i=0; $i<=$#{$amplicon_seq_data}; $i++){
					foreach my $md5 (keys %{$sommer_chimeras->[$i]{$sample_name}}){
						my $depth = $amplicon_seq_data->[$i]{$marker_name}{$sample_name}{$md5}{'depth'};
						my $another_sample;
						foreach my $sample_name_ (@{$samples->{$marker_name}}){
							if ($sample_name_ ne $sample_name && defined($amplicon_seq_data->[$i]{$marker_name}{$sample_name_}{$md5})){
								$another_sample = sprintf("individual %s",$sample_name_);
								last;
							}
							if (defined($another_sample)) { last; }
						}
						if (!defined($another_sample)){
							$sommer_unclassified->[$i]{$sample_name}{$md5} = sprintf("Unclassified: Step III, chimera (%s), not in other individuals",$sommer_chimeras->[$i]{$sample_name}{$md5});
						} else {
							my $is_allele = 0;
							if (defined($sommer_alleles->[$i-1]{$sample_name}{$md5})){
								$is_allele = 1;
								$sommer_alleles->[$i]{$sample_name}{$md5} = sprintf("Allele: Step III, chimera (%s), %s, allele in replicate",$sommer_chimeras->[$i]{$sample_name}{$md5},$another_sample);
							} else {
								$sommer_unclassified->[$i]{$sample_name}{$md5} = sprintf("Unclassified: Step III, chimera (%s), %s, no allele in replicate", $sommer_chimeras->[$i]{$sample_name}{$md5},$another_sample);
							}
						}
						delete($sommer_chimeras->[$i]{$sample_name}{$md5});
					}
					#print "\n".join("\n",@{$sommer_alleles->[$i]{$sample_name}})."\n";
				}
			}

			# Annonates Sommer genotyping results
			foreach my $sample_name (@{$samples->{$marker_name}}){
				# Process only samples with sequences in the original order
				if (!defined($amplicon_seq_data->[0]{$marker_name}{$sample_name})){
					next;
				}
				# Exclude samples not in the list
				if (defined($genotyping_parameters->{'allowed_samples'}) && ref($genotyping_parameters->{'allowed_samples'}) eq 'ARRAY' && !in_array($genotyping_parameters->{'allowed_samples'},$sample_name)){
					next;
				}
				my @sorted_md5s =  sort { $amplicon_seq_data->[0]{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->[0]{$marker_name}{$sample_name}{$a}{'depth'} }  keys %{$amplicon_seq_data->[0]{$marker_name}{$sample_name}};
				foreach my $md5 (@sorted_md5s){
					my $seq = $marker_seq_data->[0]{$marker_name}{$md5}{'seq'};
					my $name = $marker_seq_data->[0]{$marker_name}{$md5}{'name'};
					my $len = $marker_seq_data->[0]{$marker_name}{$md5}{'len'};
					my $depth = $amplicon_seq_data->[0]{$marker_name}{$sample_name}{$md5}{'depth'};
					my $frequency = $amplicon_seq_data->[0]{$marker_name}{$sample_name}{$md5}{'freq'};
					my $count_samples = $marker_seq_data->[0]{$marker_name}{$md5}{'samples'};
					my $mean_freq = $marker_seq_data->[0]{$marker_name}{$md5}{'mean_freq'};
					my $min_freq = $marker_seq_data->[0]{$marker_name}{$md5}{'min_freq'};
					my $max_freq = $marker_seq_data->[0]{$marker_name}{$md5}{'max_freq'};
	# 				my $genotyp_criteria;
	# 				if (defined($sommer_unclassified->[0]{$sample_name}{$md5})){
	# 					$genotyp_criteria = sprintf("unclassified variant");
	# 				} elsif (defined($sommer_artifacts->[0]{$sample_name}{$md5})){
	# 					$genotyp_criteria = sprintf("artifact");
	# 				}
					my $header = sprintf("hash=%s | len=%d | depth=%d | freq=%.2f | samples=%d | mean_freq=%.2f | max_freq=%.2f | min_freq=%.2f", $md5, $len, $depth, $frequency, $count_samples, $mean_freq, $max_freq, $min_freq);
					if (defined($sommer_artifacts->[0]{$sample_name}{$md5})){
						$genotyping_output->{$marker_name}{$sample_name} .= sprintf(">#%s | %s | %s\n%s\n", $name, $sommer_artifacts->[0]{$sample_name}{$md5}, $header, $seq);
					} elsif (defined($sommer_unclassified->[0]{$sample_name}{$md5})){
						$genotyping_output->{$marker_name}{$sample_name} .= sprintf(">#%s | %s | %s\n%s\n", $name, $sommer_unclassified->[0]{$sample_name}{$md5}, $header, $seq);
					} elsif (defined($sommer_wrong_length->[0]{$sample_name}{$md5})){
						$genotyping_output->{$marker_name}{$sample_name} .= sprintf(">#%s | %s | %s\n%s\n", $name, $sommer_wrong_length->[0]{$sample_name}{$md5}, $header, $seq);
					} elsif (defined($sommer_alleles->[0]{$sample_name}{$md5})){
						$genotyped_amplicon_sequences->{$marker_name}{$sample_name}{$md5} = $depth;
						$genotyped_amplicon_depths->{$marker_name}{$sample_name} += $depth;
						$genotyping_output->{$marker_name}{$sample_name} .= sprintf(">*%s | %s | %s\n%s\n", $name, $sommer_alleles->[0]{$sample_name}{$md5}, $header, $seq);
					} else {
						print sprintf("\nERROR: Variant '%s' has not been correctly classified.\n", $name);
					}
					if (!defined($md5_to_sequence->{$md5})){
						$md5_to_sequence->{$md5} = $seq;
						$md5_to_name->{$md5} = $name;
					}
				}
				if ($sample_name ne ''){
					printf("\t%s-%s genotyped (%d reads, %d variants)\n", $marker_name, $sample_name, $genotyped_amplicon_depths->{$marker_name}{$sample_name}, scalar keys %{$genotyped_amplicon_sequences->{$marker_name}{$sample_name}});
				} else {
					printf("\t%s genotyped (%d reads, %d variants)\n", $marker_name, $genotyped_amplicon_depths->{$marker_name}{$sample_name}, scalar keys %{$genotyped_amplicon_sequences->{$marker_name}{$sample_name}});
				}
			}
		}
		# End Sommer method

	}
	

	return ($genotyped_amplicon_sequences,$genotyped_amplicon_depths,$genotyping_output);

}


#################################################################################

# Genotypes amplicon sequences with the desired method ('Sommer', 'Lighten' or 'Herdegen')
sub genotype_amplicon_sequences_with_threads {

	my ($method,$markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data,$threads_limit) = @_;

	# Sommer method cannot be run in parallel
	if ($method eq 'sommer'){
		return genotype_amplicon_sequences($method,$markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$amplicon_seq_data);
	}

	if (!defined($threads_limit)){
		$threads_limit = 4;
	}

	my ($genotyped_amplicon_sequences, $genotyped_amplicon_depths, $genotyping_output);

	my $one_amplicon_seq_data;
	foreach my $marker_name (@$markers){
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}
		my $genotyped_amplicon_depths->{$marker_name} = {};
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
				next;
			}
			my $one_amplicon_seq_data_;
			my $genotyped_amplicon_sequences->{$marker_name}{$sample_name} = {};
			$one_amplicon_seq_data_->{$marker_name}{$sample_name} = $amplicon_seq_data->{$marker_name}{$sample_name};
			push(@{$one_amplicon_seq_data},$one_amplicon_seq_data_);
			# print '';
		}
	}

	my @threads;
	for (my $count_amplicon=0; $count_amplicon<=$#{$one_amplicon_seq_data}; $count_amplicon++){

		push(@threads, threads->create(\&genotype_amplicon_sequences,$method,$markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$one_amplicon_seq_data->[$count_amplicon]));
# 		print "\n";

# 		# For debugging:
# 		push(@threads, [genotyp_amplicon_sequences($markers,$samples,$markerdata,$paramsdata,$marker_seq_data,$one_amplicon_seq_data->[$count_amplicon],$amplicon_depths)]);

		# If maximum number of threads is reached or last sbjct of a query is processed
		if (scalar @threads >= $threads_limit  || $count_amplicon == $#{$one_amplicon_seq_data}){
			my $check_threads = 1;
			while ($check_threads){
				for (my $i=0; $i<=$#threads; $i++){
					unless ($threads[$i]->is_running()){
						my ($genotyped_amplicon_sequences_,$genotyped_amplicon_depths_,$genotyping_output_) = $threads[$i]->join;
						if (defined($genotyped_amplicon_sequences_)){
							my $marker_name = (keys %$genotyped_amplicon_sequences_)[0];
							my $sample_name = (keys %{$genotyped_amplicon_sequences_->{$marker_name}})[0];
	# 						print "\t$marker_name-$sample_name finished\n";
							if (defined($genotyped_amplicon_sequences_->{$marker_name}{$sample_name})){
								$genotyped_amplicon_sequences->{$marker_name}{$sample_name} = $genotyped_amplicon_sequences_->{$marker_name}{$sample_name};
								$genotyped_amplicon_depths->{$marker_name}{$sample_name} = $genotyped_amplicon_depths_->{$marker_name}{$sample_name};
								$genotyping_output->{$marker_name}{$sample_name} .= $genotyping_output_->{$marker_name}{$sample_name};
							}
						}
						undef($threads[$i]);
						splice(@threads,$i,1);
						$i = $i - 1;
						unless ($count_amplicon == $#{$one_amplicon_seq_data} && @threads){
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
# 				print '';
# 				my ($genotyped_amplicon_sequences_,$genotyped_amplicon_depths_,$genotyping_output_) = @{$threads[$i]};
# 				my $marker_name = (keys %$genotyped_amplicon_sequences_)[0];
# 				my $sample_name = (keys %{$genotyped_amplicon_sequences_->{$marker_name}})[0];
# 				$genotyped_amplicon_sequences->{$marker_name}{$sample_name} = $genotyped_amplicon_sequences_->{$marker_name}{$sample_name};
# 				$genotyped_amplicon_depths->{$marker_name}{$sample_name} = $genotyped_amplicon_depths_->{$marker_name}{$sample_name};
# 				$genotyping_output->{$marker_name}{$sample_name} .= $genotyping_output_->{$marker_name}{$sample_name};
# 				delete $threads[$i];
# 			}
		}
	}
	
	return ($genotyped_amplicon_sequences,$genotyped_amplicon_depths,$genotyping_output);

}

#################################################################################

# Analyzes sequence frequency per amplicon to extract real sequences/alleles
sub extract_alleles_freq_threshold {

	my ($markers, $samples, $markerdata, $low_depth, $amplicon_sample_assignment_seqs, $md5_to_amplicon_seq, $repair_errors, $min_amplicon_seq_depth, $max_amplicon_repairing_frequency, $min_amplicon_seq_frequency, $min_allele_mean_frequency) = @_;

	print "\nGenotyping data using 'Per Amplicon Frecuency Thresholds' method.\n";

	my ($real_alleles, %md5_to_real_alleles);
	my ($sorted_amplicon_md5s, $sorted_amplicon_seqs, $sorted_amplicon_names, $sorted_marker_md5s, $sorted_marker_seqs, $sorted_marker_names, $rank_thresholds);
	# Correct amplicon assignments if sequences are repaired
	my $amplicon_assignment_repaired_seqs;
	# Stores frequencies of each amplicon in each sample (individual values, not sumatory)
	my $amplicon_frequency_per_sample;

	foreach my $marker_name (@$markers){

		# Stores the sum of the frequencies of each amplicon in each sample
		my %amplicon_frequency_total;
		# Stores the highest frequency of each amplicon in each sample
		my %amplicon_frequency_highest;

		# Loop all the samples and annotate amplicon frequencies
		foreach my $sample_name (@{$samples->{$marker_name}}) {

			# Exclude from further analysis sequence with low coverage
			if (defined($low_depth->{$marker_name}{$sample_name})){
				next;
			}

			# Order the sequences by their coverage of the sample (number of reads with the sequence)
			$sorted_amplicon_md5s->{$marker_name}{$sample_name} = [ sort { $amplicon_sample_assignment_seqs->{$marker_name}{$sample_name}{$b} <=> $amplicon_sample_assignment_seqs->{$marker_name}{$sample_name}{$a} } keys %{$amplicon_sample_assignment_seqs->{$marker_name}{$sample_name}} ];
			$sorted_amplicon_seqs->{$marker_name}{$sample_name} = [ map $md5_to_amplicon_seq->{$_}, @{$sorted_amplicon_md5s->{$marker_name}{$sample_name}} ];
			my $total_seqs=0;
			map $total_seqs+=$amplicon_sample_assignment_seqs->{$marker_name}{$sample_name}{$_}, keys %{$amplicon_sample_assignment_seqs->{$marker_name}{$sample_name}};

			# Loop all the putative alleles with coverage higher than threshold
			my %amplicon_frequency;
			my $count_seqs = 1;
			for (my $i=0; $i<=$#{$sorted_amplicon_md5s->{$marker_name}{$sample_name}}; $i++) {
				my $md5 = $sorted_amplicon_md5s->{$marker_name}{$sample_name}[$i];
				# Annotate names
				my $name;
				if (defined($md5_to_real_alleles{$md5})){
					$name = $md5_to_real_alleles{$md5};
				} else {
					$name = sprintf('%s-%s-%03d', $marker_name, $sample_name, $count_seqs);
					$count_seqs++;
				}
				push(@{$sorted_amplicon_names->{$marker_name}{$sample_name}},$name);

				my $sample_amplicon_seq_depth = $amplicon_sample_assignment_seqs->{$marker_name}{$sample_name}{$md5};
				my $sample_amplicon_seq_frequency = $sample_amplicon_seq_depth/$total_seqs;
				if ($sample_amplicon_seq_depth >= $min_amplicon_seq_depth) {
					# Assign the frequency values to the reference sequence (without errors)
					if (defined($repair_errors) && $sample_amplicon_seq_frequency < $max_amplicon_repairing_frequency) {
						my $seq = $sorted_amplicon_seqs->{$marker_name}{$sample_name}[$i];
						for (my $j=0; $j<=$#{$sorted_amplicon_seqs->{$marker_name}{$sample_name}}; $j++){
							my $ref_seq =$sorted_amplicon_seqs->{$marker_name}{$sample_name}[$j];
							my ($error_type, $error_pos) = compare_sequences($seq, $ref_seq);
							# Repair only indels (substitution can be real alleles)
							# Read will be assigned to the major reference sequence
							if ($error_type ne 'substitution'){
								$md5 = $sorted_amplicon_md5s->{$marker_name}{$sample_name}[$j];
								last;
							}
						}
					}
					$amplicon_frequency{$md5} += $sample_amplicon_seq_frequency;
					$amplicon_assignment_repaired_seqs->{$marker_name}{$sample_name}{$md5} += $sample_amplicon_seq_depth;
				} else {
					# Sequences are ordered by coverage, so next sequences will have also low depth
					next;
				}
			}

			# Store the values of frequencies per sample
			foreach my $md5 (keys %amplicon_frequency) {
				# Annotates only frequencies higher than threshold
				if ($amplicon_frequency{$md5} >= $min_amplicon_seq_frequency){
					$amplicon_frequency_total{$md5} += $amplicon_frequency{$md5};
					push(@{$amplicon_frequency_per_sample->{$marker_name}{$md5}}, $amplicon_frequency{$md5});
					if (!defined($amplicon_frequency_highest{$md5}) || $amplicon_frequency_highest{$md5} < $amplicon_frequency{$md5}){
						$amplicon_frequency_highest{$md5} = $amplicon_frequency{$md5};
					}
				}
			}
		}

		# Order marker sequences by the total sum of frequencies per sample
	# 	$sorted_marker_md5s->{$marker_name} = [ sort { $amplicon_frequency_total{$b} <=> $amplicon_frequency_total{$a} } keys %amplicon_frequency_total ];
		$sorted_marker_md5s->{$marker_name} = [ sort { $amplicon_frequency_highest{$b} <=> $amplicon_frequency_highest{$a} } keys %amplicon_frequency_highest ];
	# 	$sorted_marker_md5s->{$marker_name} = [ sort { mean(@{$amplicon_frequency_per_sample->{$marker_name}{$b}}) <=> mean(@{$amplicon_frequency_per_sample->{$marker_name}{$a}}) } keys %amplicon_frequency_total ];
		if (!@{$sorted_marker_md5s->{$marker_name}}){
			print "\n\tThere are not samples with enough coverage to analyze '$marker_name' amplicon.\n";
			next;
		}
		$sorted_marker_seqs->{$marker_name} = [ map $md5_to_amplicon_seq->{$_}, @{$sorted_marker_md5s->{$marker_name}} ];
		$sorted_marker_names->{$marker_name} = [ map sprintf('%s-%03d', $marker_name, $_), 1 .. scalar @{$sorted_marker_md5s->{$marker_name}} ];


		# Will store the threshold value to annotate real alleles
		my $rank_threshold=0;
		my $freq_threshold=0;

		# Stores the frequency values to check where are located the real alleles
		my @frequencies_above_rank;
		# Calculates the mean of the frequencies of the alleles in lower rank positions
		my @frequencies_below_rank;
		# At the beginning all the frequencies are below rank
		for (my $i=0; $i<=$#{$sorted_marker_md5s->{$marker_name}}; $i++) {
			foreach my $frequency_per_sample (@{$amplicon_frequency_per_sample->{$marker_name}{$sorted_marker_md5s->{$marker_name}[$i]}}){
				push(@frequencies_below_rank, $frequency_per_sample);
			}
		}

		# Calculates frequency threshold of real alleles
		for (my $i=0; $i<=$#{$sorted_marker_md5s->{$marker_name}}; $i++) {
			my $md5 = $sorted_marker_md5s->{$marker_name}[$i];
			my $name = $sorted_marker_names->{$marker_name}[$i];
			my $highest_amplicon_frequency = 0;
			# Include frequencies in above rank array and remove them from below rank array
			foreach my $frequency_per_sample (@{$amplicon_frequency_per_sample->{$marker_name}{$md5}}){
				push(@frequencies_above_rank, $frequency_per_sample);
				if ($frequency_per_sample > $highest_amplicon_frequency){
					$highest_amplicon_frequency = $frequency_per_sample;
				}
				shift @frequencies_below_rank;
			}

			my $mean_frequencies_below_rank = mean(@frequencies_below_rank);
			my $mean_frequencies_above_rank = mean(@frequencies_above_rank);
			my $mean_amplicon_frequency = mean(@{$amplicon_frequency_per_sample->{$marker_name}{$md5}});

	# if ($mean_amplicon_frequency < $min_allele_mean_frequency){
	# print '';
	# }
			# Defines the threshold value to annotate real alleles
			# We can also use FREQ thresholds included in the amplicon CSV data file
			if (!defined($mean_frequencies_below_rank)){
				$rank_threshold = $i;
				$freq_threshold = $mean_amplicon_frequency;
			} elsif (!defined($markerdata->{$marker_name}{'freq'}) && !$rank_threshold && $mean_frequencies_below_rank < $min_allele_mean_frequency && $highest_amplicon_frequency < 4*$mean_frequencies_below_rank){
				$rank_threshold = $i-1;
				$freq_threshold = 4*$mean_frequencies_below_rank;
				last;
			} elsif (defined($markerdata->{$marker_name}{'freq'}) && !$rank_threshold && $mean_frequencies_below_rank < $markerdata->{$marker_name}{'freq'}){
				$rank_threshold = $i-1;
				$freq_threshold = $markerdata->{$marker_name}{'freq'};
				last;
			}
			my %allele_data = ('md5' => $md5, 'name' => $name);
			push(@{$real_alleles->{$marker_name}},\%allele_data);
			$md5_to_real_alleles{$md5} = $name;
		}
		$rank_thresholds->{$marker_name} = $rank_threshold;
	}

	return ($real_alleles, $sorted_amplicon_md5s, $sorted_amplicon_seqs, $sorted_amplicon_names, $sorted_marker_md5s, $sorted_marker_seqs, $sorted_marker_names, $rank_thresholds, $amplicon_assignment_repaired_seqs, $amplicon_frequency_per_sample);

}

#################################################################################

# Reads amplicons from a multifile containing demultiplexed sequence files
sub read_compressed_file_amplicons {

	my ($infile,$options,$max_depth) = @_;

	my ($type, $verbose, $shuffle) = (undef, 1, 0);
	if (in_array($options, 'zip')){
		$type = 'zip';
	} elsif (in_array($options, 'gzip')){
		$type = 'gzip';
	}
	if (in_array($options, 'quiet')){
		$verbose = 0;
	}
	if (in_array($options, 'shuffle')){
		$shuffle = 1;
	}

	# Reads amplicon sequences and depths
	# $amplicon_sequences stores the individual unique sequence depths
	# $amplicon_sequences->{$marker_name}{$sample_name}{$md5} = $depth;
	# $amplicon_depths stores the total depth of the sequences into an amplicon
	# $amplicon_depths->{$marker_name}{$sample_name} += $depth;
	# Reads marker/amplicon sequence data
	# $marker_seq_data stores the names and parameters of all unique sequences of a unique marker
	# $marker_seq_data->{$marker_name}{$md5} = { 'seq'=> $seq, 'name'=>$name, 'len'=>$len, 'depth'=>$unique_seq_depth, 'samples'=>$count_samples, 'mean_freq'=>$mean_freq, 'max_freq'=>$max_freq, 'min_freq'=>$min_freq };
	# $amplicon_seq_data stores the names and parameters of all unique sequences of a unique amplicon
	# $amplicon_seq_data->{$marker_name}{$sample_name}{$md5} = { 'seq'=> $seq, 'name'=>$name, 'len'=>$len, 'depth'=>$unique_seq_depth, 'freq'=>$unique_seq_frequency, 'cluster_size'=>$cluster_size };

	my ($markers, $samples, $amplicon_sequences, $amplicon_depths, $md5_to_sequence);
 
	# Uncompress the files in temporal folder
	my $tmp_dir_name = random_file_name();
	my $tmp_dir = "/tmp/".$tmp_dir_name;
	my @file_list = decompress($infile,$type,$tmp_dir);

	# Stores the first nucleotides from the first sequence of each marker
	# Is used to try to extract all the sequences in the same orientation
	my %first_seq_start;

	# Each file will be an amplicon where the folder will give the marker name and the filename the sample name
	foreach my $file (@file_list) {

# 		# Stores matched reads
# 		my ($read_headers_matched, $read_seqs_matched);

		my ($seqs,$headers);

		if (-f $file && is_fasta($file)){
			($seqs,$headers) = read_fasta_file($file);
		} elsif (-f $file && is_fastq($file)){
			($seqs,$headers) = read_fastq_file($file);
		} else {
			next;
		}
		my @path_names;
		# If the compressed file contains several files (ordered or not in folders)
		if ($#file_list){
			@path_names = reverse(split("\/",$file));
		# If the compressed file contains only 1 file
		} else {
			@path_names = reverse(split("\/",$infile));
		}
		my $marker_name;
		# Gives the name of the internal folder to the marker (if it is not the temporal dir)
		if ($#file_list && $#path_names && $path_names[1] ne $tmp_dir_name){
			$marker_name = convert_to_alphanum($path_names[1]);
		# If not gives the name of the compressed file without extension
		} elsif ($infile =~ /(.+\/)?(.+?)\./) {
			$marker_name = $2;
		# If not gives the full name of the compressed file
		} else {
			$marker_name = $infile;
		}
		if (!in_array($markers, $marker_name)){
			push(@$markers, $marker_name);
		}
		my $sample_name;
		# Gives the name of the internal file to the sample (without extension)
		if ($path_names[0] =~ /(.+?)\./){
			$sample_name = $1;
		# If not gives the full name of the internal file to the sample
		} else {
			$sample_name = $path_names[0];
		}
		push(@{$samples->{$marker_name}},$sample_name);
		if ($verbose){
			print "\t$marker_name-$sample_name de-multiplexing\n";
			# printf("\tReading file '%s'\n", $file);
		}
		# Extracts the first nucleotides from the first sequence
		if (!defined($first_seq_start{$marker_name})){
			$first_seq_start{$marker_name} = uc(substr($seqs->[0],0,5));
		
		}
		# Annotates the sequences and their depths
		my $count_total_seqs_amplicon = 0;
		my @seq_indexes = 0..$#{$headers};
		if ($shuffle) {
			@seq_indexes = shuffle @seq_indexes;
		}
		foreach my $i (@seq_indexes){
			my $amplicon_header = $headers->[$i];
			my $amplicon_seq = uc($seqs->[$i]);

			my $md5 = generate_md5($amplicon_seq);

			if (!defined($md5_to_sequence->{$md5})){ 
				# Checks if the sequence is the reverse complementary of another one already stored
				my $amplicon_seq_revcomp = iupac_reverse_complementary($amplicon_seq);
				my $md5_revcomp = generate_md5($amplicon_seq_revcomp);
				if (defined($md5_to_sequence->{$md5_revcomp})){
					$amplicon_seq = $amplicon_seq_revcomp;
					$md5 = $md5_revcomp;
				# Checks if the beginning of the sequence is the same that the beginning of the first sequence
				# If not, checks if the beginning of the reverse complementary one is similar to the beginning of the first sequence
				} elsif ($first_seq_start{$marker_name} ne substr($amplicon_seq,0,5) && $first_seq_start{$marker_name} eq iupac_reverse_complementary(substr($amplicon_seq,-5))) {
					$amplicon_seq = iupac_reverse_complementary($amplicon_seq);
					$md5 = generate_md5($amplicon_seq);
				}
			}

			# Stores new sequences and checks that MD5 signature is unique
			if (!defined($md5_to_sequence->{$md5})){
				$md5_to_sequence->{$md5}=$amplicon_seq;
			} elsif ($verbose && $md5_to_sequence->{$md5} ne $amplicon_seq){
				printf("\tERROR: Sequence '%s' is not unique:\n\t\tSEQ1: %s\n\t\tSEQ2: %s\n", $md5, $md5_to_sequence->{$md5}, $amplicon_seq);
# 				exit;
			}

			if ($amplicon_header !~ /\| depth=(\d+)/){
		# 		print "\n$read_header $sample_name $marker_name $amplicon_seq\n\n";exit;
				$amplicon_sequences->{$marker_name}{$sample_name}{$md5}++;
				# Counts how many times a sample/amplicon assignment is found
				$amplicon_depths->{$marker_name}{$sample_name}++;
				# Count total number of sequences in an amplicon
				$count_total_seqs_amplicon++;
			} else {
		# 		print "\n$read_header $sample_name $marker_name $amplicon_seq\n\n";exit;
				$amplicon_sequences->{$marker_name}{$sample_name}{$md5} += $1;
				# Counts how many times a sample/amplicon assignment is found
				$amplicon_depths->{$marker_name}{$sample_name} += $1;
				# Count total number of sequences in an amplicon
				$count_total_seqs_amplicon += $1;
			}

# 			# Stores reads matched for printing
# 			push(@$read_headers_matched, $amplicon_header);
# 			push(@$read_seqs_matched, $amplicon_seq);

			
			# If defined $max_depth, limit the number of sequences per amplicon
			if (defined($max_depth) && $max_depth>0 && $count_total_seqs_amplicon>=$max_depth) {
				last;
			}
		}

		if ($verbose){
			printf("\t%s-%s de-multiplexed (%d reads, %d variants)\n", $marker_name, $sample_name, $count_total_seqs_amplicon, scalar keys %{$amplicon_sequences->{$marker_name}{$sample_name}});
		}

# 		# Creates a file with reads matched
# 		if ($file =~ /(.+\/)?(.+?)\./){
# 			print "\nMatched reads printed into '$2.fasta'.\n\n";
# 			create_fasta_file($read_seqs_matched,$read_headers_matched,"$2.fasta");
# 		}
	}

	`rm -rf $tmp_dir`;
# 	exit;

	# Sort markers and samples by name
	$markers = [ nsort(@$markers) ];
	foreach my $marker_name (@$markers){
		$samples->{$marker_name} = [ nsort(@{$samples->{$marker_name}}) ];
	}

	return ($markers,$samples,$amplicon_sequences, $amplicon_depths, $md5_to_sequence);
}

#################################################################################

# Reads AmpliSAS results from Excel file
sub read_amplisas_file_amplicons {

	my ($file,$options,$max_depth) = @_;

	my ($type, $verbose) = ('amplisas', 1);
	if (in_array($options, 'amplicheck')){
		$type = 'amplicheck';
	}
	if (in_array($options, 'quiet')){
		$verbose = 0;
	}

	my $converter = Text::Iconv -> new ("utf-8", "windows-1251");
	my $excel = Spreadsheet::XLSX -> new ($file, $converter);
 
	# Reads amplicon sequences and depths
	# $amplicon_sequences stores the individual unique sequence depths
	# $amplicon_sequences->{$marker_name}{$sample_name}{$md5} = $depth;
	# $amplicon_depths stores the total depth of the sequences into an amplicon
	# $amplicon_depths->{$marker_name}{$sample_name} += $depth;
	# Reads marker/amplicon sequence data
	# $marker_seq_data stores the names and parameters of all unique sequences of a unique marker
	# $marker_seq_data->{$marker_name}{$md5} = { 'seq'=> $seq, 'name'=>$name, 'len'=>$len, 'depth'=>$unique_seq_depth, 'samples'=>$count_samples, 'mean_freq'=>$mean_freq, 'max_freq'=>$max_freq, 'min_freq'=>$min_freq };
	# $amplicon_seq_data stores the names and parameters of all unique sequences of a unique amplicon
	# $amplicon_seq_data->{$marker_name}{$sample_name}{$md5} = { 'seq'=> $seq, 'name'=>$name, 'len'=>$len, 'depth'=>$unique_seq_depth, 'freq'=>$unique_seq_frequency, 'cluster_size'=>$cluster_size };

	my ($markers, $samples, $amplicon_sequences, $amplicon_depths, %md5_to_sequence);
 
	foreach my $sheet (@{$excel -> {Worksheet}}) {
		my $marker = $sheet->{Name};
		if ($marker =~ /(.+)_depths$/)  {
			$type = 'amplicheck';
			$marker = $1;
		} elsif ($type eq 'amplicheck' && $marker !~ /_depths$/)  {
			next;
		}
		push(@$markers, $marker);
# 		if ($verbose){
# 			printf("\tReading Sheet '%s'\n", $sheet->{Name});
# 		}
		my (@seq_md5s, %seq_headers, %sample_headers);
		my (%seqs_per_sample, %samples_per_seq);
		my %md5_to_name;
		my ($samples_col,$headers_row);
		my $sample_data;
		$sheet -> {MaxRow} ||= $sheet -> {MinRow};
		foreach my $row ($sheet -> {MinRow} .. $sheet -> {MaxRow}) {
			my ($md5,$read_sample_data,%seq_data,%sample_assigns);
			$sheet -> {MaxCol} ||= $sheet -> {MinCol};
			foreach my $col ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {
				my $cell_data = $sheet -> {Cells} [$row] [$col];
				my $cell = $cell_data->{Val};
				if (!defined($cell) && (!defined($headers_row) || $row != $headers_row)) {
					next;
				# Skips AmpliCHECK legend
				} elsif (defined($cell) && $cell =~ /Good sequence|Suspicious sequence|Probably artifact/i){
					next;
				} elsif (!%seq_headers && defined($cell) && $cell =~/DEPTH_AMPLICON|DEPTH_ALLELES|DEPTH_SEQUENCES|COUNT_ALLELES|COUNT_SEQUENCES/i) {
					$read_sample_data = lc($cell);
					if (!defined($samples_col)){
						$samples_col = $col+1;
					}
				} elsif (defined($read_sample_data)){
					push(@{$sample_data->{$read_sample_data}}, $cell);
				} elsif (defined($samples_col) && $col<$samples_col && defined($cell) && $cell =~/SEQ|MD5|HASH|LEN|DEPTH|COV|SAMP|FREQ|NAME/i) {
					if ($cell =~/SEQ/i) {
						$seq_headers{$col}='seq';
					}elsif ($cell =~/LEN/i) {
						$seq_headers{$col}='len';
					}elsif ($cell =~/COV/i) {
						$seq_headers{$col}='depth';
					}elsif ($cell =~/SAMP/i) {
						$seq_headers{$col}='samples';
					} else {
						$seq_headers{$col}=lc($cell);
					}
					$read_sample_data = undef;
					if (!defined($headers_row)){
						$headers_row = $row;
					}
				} elsif (defined($headers_row) && $row==$headers_row && defined($samples_col) && $col<$samples_col && (!defined($cell) || $cell eq '')){
					$seq_headers{$col}='name';
				} elsif (defined($headers_row) && $row==$headers_row && defined($samples_col) && $col>=$samples_col){
					$sample_headers{$col}=$cell;
					push(@{$samples->{$marker}},$cell);
				} elsif (defined($headers_row) && $row>$headers_row && %seq_headers && %sample_headers){
					if (defined($seq_headers{$col})){
						$seq_data{lc($seq_headers{$col})}=$cell;
						if (lc($seq_headers{$col}) =~ /seq/i){
							$md5 = generate_md5($cell);
							if (!defined($md5_to_sequence{$md5})){
								$md5_to_sequence{$md5} = $cell;
							}
						}
						if (lc($seq_headers{$col}) =~ /name/i && !defined($md5_to_name{$md5})){
							$md5_to_name{$md5} = $cell;
						}
# 						if (lc($seq_headers{$col}) =~ /md5/i){
# 							$md5 = $cell;
# 							push(@seq_md5s,$md5);
# 						}
					} elsif (defined($md5) && defined($sample_headers{$col}) && (is_numeric($cell) || $cell =~ /^(\d+);/)) {
						# If the result comes from AmpliCHECK
						if ($1) { 
							$sample_assigns{$sample_headers{$col}} = $1 + 0; # to obtain a number without leading zeros
						} else {
							$sample_assigns{$sample_headers{$col}}=$cell;
						}
						$seqs_per_sample{$sample_headers{$col}}++;
					}
				}
# 				if ($cell) {
# 					printf("( %s , %s ) => %s\n", $row, $col, $cell_data -> {Val});
# 				}
			}
			if (scalar keys %seq_headers < 3){
				undef(%seq_headers);
			}
			# Annotate only if it's specified the sequence data and the sequence has sample assignments
			if (%seq_data && %sample_assigns) {
				foreach my $sample (keys %sample_assigns) {
					$amplicon_sequences->{$marker}{$sample}{$md5} = $sample_assigns{$sample};
# 					$marker_seq_data->{$marker}{$md5} = \%seq_data;
				}
			}
		}
		# Annotates full amplicon depth with original Excel value
		for (my $i=0; $i<=$#{$samples->{$marker}}; $i++) {
			my $sample = $samples->{$marker}[$i];
			if ($verbose){
				print "\t$marker-$sample de-multiplexing\n";
			}
# if ($sample eq 'Ion013'){
# print '';
# }
			if (defined($sample_data->{'depth_amplicon'})){
				if (!defined($max_depth) || $max_depth >= $sample_data->{'depth_amplicon'}[$i]){
					$amplicon_depths->{$marker}{$sample} = $sample_data->{'depth_amplicon'}[$i];
# 					# Calculates amplicon frequencies based in amplicon depth value
# 					foreach my $md5 (keys %{$amplicon_sequences->{$marker}{$sample}}){
# 						$amplicon_seq_data->{$marker}{$sample}{$md5}{'name'} = $md5_to_name{$md5};
# 						$amplicon_seq_data->{$marker}{$sample}{$md5}{'seq'} = $md5_to_sequence{$md5};
# 						$amplicon_seq_data->{$marker}{$sample}{$md5}{'depth'} = $amplicon_sequences->{$marker}{$sample}{$md5};
# 						if (is_numeric($amplicon_sequences->{$marker}{$sample}{$md5}) && $amplicon_depths->{$marker}{$sample}>0){
# 							$amplicon_seq_data->{$marker}{$sample}{$md5}{'freq'} = sprintf("%.2f",100*$amplicon_seq_data->{$marker}{$sample}{$md5}{'depth'}/$amplicon_depths->{$marker}{$sample});
# 						}
# 					}
				} else {
					my $ratio = 0;
					if ($amplicon_depths->{$marker}{$sample}){
						$ratio = 0.5*($max_depth/$amplicon_depths->{$marker}{$sample});
					}
					$amplicon_depths->{$marker}{$sample} = 0;
					# Recalculates amplicon frequencies based in maximum amplicon depth value
					my %amplicon_sequences_ = %{$amplicon_sequences->{$marker}{$sample}};
					while ($max_depth > $amplicon_depths->{$marker}{$sample}){
						delete($amplicon_sequences->{$marker}{$sample});
# 						delete($amplicon_seq_data->{$marker}{$sample});
						$amplicon_depths->{$marker}{$sample} = 0;
						foreach my $md5 (sort {$amplicon_sequences_{$b}<=>$amplicon_sequences_{$a}} keys %amplicon_sequences_){
							my $new_depth = sprintf("%.0f",$ratio*$amplicon_sequences_{$md5});
							$amplicon_depths->{$marker}{$sample} += $new_depth;
							if ($new_depth) {
								$amplicon_sequences->{$marker}{$sample}{$md5} = $new_depth;
# 								$amplicon_seq_data->{$marker}{$sample}{$md5}{'name'} = $md5_to_name{$md5};
# 								$amplicon_seq_data->{$marker}{$sample}{$md5}{'seq'} = $md5_to_sequence{$md5};
# 								$amplicon_seq_data->{$marker}{$sample}{$md5}{'depth'} = $amplicon_sequences_{$md5};
# 								if (is_numeric($amplicon_seq_data->{$marker}{$sample}{$md5}{'depth'}) && $amplicon_depths->{$marker}{$sample}>0){
# 									$amplicon_seq_data->{$marker}{$sample}{$md5}{'freq'} = sprintf("%.2f",100*$amplicon_seq_data->{$marker}{$sample}{$md5}{'depth'}/$amplicon_depths->{$marker}{$sample});
# 								}
							}
							if ($max_depth <= $amplicon_depths->{$marker}{$sample}){
								last;
							}
						}
						$ratio  = $ratio*1.1;
					}
				}
# print '';
			}
			if ($verbose){
				printf("\t%s-%s de-multiplexed (%d reads, %d variants)\n", $marker, $sample, $amplicon_depths->{$marker}{$sample}, scalar keys %{$amplicon_sequences->{$marker}{$sample}});
			}
		}
	}
	
	return ($markers,$samples,$amplicon_sequences, $amplicon_depths, \%md5_to_sequence);
}

#################################################################################

# Reads AmpliSAS results from Excel file
sub read_amplisas_file_results {

	my($file,$options) = @_;

	my ($type, $verbose) = ('amplisas', 1);
	if (in_array($options, 'amplicheck')){
		$type = 'amplicheck';
	}
	if (in_array($options, 'quiet')){
		$verbose = 0;
	}

	my $amplisas_results;
	
	my $converter = Text::Iconv -> new ("utf-8", "windows-1251");
	my $excel = Spreadsheet::XLSX -> new ($file, $converter);

	if (!defined($verbose)){
		$verbose = 0;
	}

	foreach my $sheet (@{$excel -> {Worksheet}}) {
		my $marker = $sheet->{Name};
		if ($marker =~ /(.+)_depths$/)  {
			$type = 'amplicheck';
			$marker = $1;
		} elsif ($type eq 'amplicheck' && $marker !~ /_depths$/)  {
			next;
		}
		if ($verbose){
			printf("\tReading Sheet '%s'\n", $sheet->{Name});
		}
		$sheet -> {MaxRow} ||= $sheet -> {MinRow};
		my (@seq_md5s, @samples, %seq_headers, %sample_headers);
		my (%seqs_per_sample, %samples_per_seq);
		my ($samples_col,$headers_row);
		my $sample_data;
		foreach my $row ($sheet -> {MinRow} .. $sheet -> {MaxRow}) {
			$sheet -> {MaxCol} ||= $sheet -> {MinCol};
			my ($md5,$read_sample_data,%seq_data,%sample_assigns);
			foreach my $col ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {
				my $cell_data = $sheet -> {Cells} [$row] [$col];
				my $cell = $cell_data->{Val};
				if (!defined($cell) && (!defined($headers_row) || $row != $headers_row)) {
					next;
				# Skips AmpliCHECK legend
				} elsif (defined($cell) && $cell =~ /Good sequence|Suspicious sequence|Probably artifact/i){
					next;
				} elsif (!%seq_headers && defined($cell) && $cell =~/DEPTH_AMPLICON|DEPTH_ALLELES|DEPTH_SEQUENCES|COUNT_ALLELES|COUNT_SEQUENCES/i) {
					$read_sample_data = lc($cell);
					if (!defined($samples_col)){
						$samples_col = $col+1;
					}
				} elsif (defined($read_sample_data)){
					push(@{$sample_data->{$read_sample_data}}, $cell);
				} elsif (defined($samples_col) && $col<$samples_col && defined($cell) && $cell =~/SEQ|MD5|HASH|LEN|DEPTH|COV|SAMP|FREQ|NAME|EFF/i) {
					if ($cell =~/SEQ/i) {
						$seq_headers{$col}='sequence';
					}elsif ($cell =~/LEN/i) {
						$seq_headers{$col}='length';
					}elsif ($cell =~/COV/i) {
						$seq_headers{$col}='depth';
					}elsif ($cell =~/SAMP/i) {
						$seq_headers{$col}='samples';
					} else {
						$seq_headers{$col}=lc($cell);
					}
					$read_sample_data = undef;
					if (!defined($headers_row)){
						$headers_row = $row;
					}
				} elsif (defined($headers_row) && $row==$headers_row && defined($samples_col) && $col<$samples_col && (!defined($cell) || $cell eq '')){
					$seq_headers{$col}='name';
				} elsif (defined($headers_row) && $row==$headers_row && defined($samples_col) && $col>=$samples_col){
					$sample_headers{$col}=$cell;
					push(@samples,$cell);
				} elsif (defined($headers_row) && $row>$headers_row && %seq_headers && %sample_headers){
					if (defined($seq_headers{$col})){
						$seq_data{lc($seq_headers{$col})}=$cell;
						if (lc($seq_headers{$col}) =~ /seq/i){
							$md5 = generate_md5($cell);
							push(@seq_md5s,$md5);
						}
# 						if (lc($seq_headers{$col}) =~ /md5/i){
# 							$md5 = $cell;
# 							push(@seq_md5s,$md5);
# 						}
					} elsif (defined($md5) && defined($sample_headers{$col})){
						# If the result comes from AmpliCHECK
						if ($cell =~ /^(\d+);/ && is_numeric($1)) {
							$sample_assigns{$sample_headers{$col}} = $1 + 0; # to obtain a number without leading zeros
							$seqs_per_sample{$sample_headers{$col}}++;
						} elsif (is_numeric($cell)) {
							$sample_assigns{$sample_headers{$col}}=$cell;
							$seqs_per_sample{$sample_headers{$col}}++;
						}
					}
				}
# 				if ($cell) {
# 					printf("( %s , %s ) => %s\n", $row, $col, $cell_data -> {Val});
# 				}
			}
			if (scalar keys %seq_headers < 3){
				undef(%seq_headers);
			}
			# Annotate only if it's specified the sequence data and the sequence has sample assignments
			if (%seq_data && %sample_assigns) {
				$amplisas_results->{$marker}{'seq_data'}{$md5} = \%seq_data;
				$amplisas_results->{$marker}{'assignments'}{$md5} = \%sample_assigns;
			} elsif (defined($md5)) {
				pop(@seq_md5s);
			}
		}
		# Removes samples without sequence assignments
		foreach my $sample (@samples) {
			if (defined($seqs_per_sample{$sample})){
				push(@{$amplisas_results->{$marker}{'samples'}}, $sample);
			}
			# Annotates amplicon depth, allele depth and number of alleles
			if ($sample_data){
				foreach my $data_param (keys %{$sample_data}){
					$amplisas_results->{$marker}{'sample_data'}{$sample}{$data_param} = shift @{$sample_data->{$data_param}};
				}
			}
		}
		$amplisas_results->{$marker}{'seq_md5s'} = \@seq_md5s;
	}
	
	return $amplisas_results;
}

#################################################################################

# Writes AmpliSAS results into an Excel file
sub write_amplisas_file_results {

	my ($file,$results,$properties) = @_;

	if (!defined($properties)){
		$properties = {
			'title'    => "AmpliSAS results file",
			'author'   => "Alvaro Sebastian",
			'comments' => "AmpliSAS results file",
			'company'  => "SixthResearcher",
		};
	}

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

	foreach my $marker_name (keys %{$results}){

		# Obtains all the variants in all the files
		my (%seq_depths,%final_seq_depths,%final_seq_samples);
		my (%md5_to_sequence, %md5_to_name);
		if (defined($results->{$marker_name})){
			foreach my $md5 (keys %{$results->{$marker_name}{'seq_data'}}){
				$seq_depths{$md5} += $results->{$marker_name}{'seq_data'}{$md5}{'depth'};
				$md5_to_sequence{$md5} = $results->{$marker_name}{'seq_data'}{$md5}{'sequence'};
				$md5_to_name{$md5} = $results->{$marker_name}{'seq_data'}{$md5}{'name'};
			}
		}
		if (!%seq_depths) { next; }

		my @seq_md5s = sort { $seq_depths{$b} <=> $seq_depths{$a} } keys %seq_depths;
		
		# Creates worksheet and writes data headers
		my $worksheet = $workbook->add_worksheet("$marker_name");
		$worksheet->set_column('F:F', undef, $bold);
		my $ws_row = 0;
		my @seq_data_headers = ('SEQUENCE', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'NAME');
		$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_AMPLICON'); $ws_row++;
		$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_ALLELES'); $ws_row++;
		$worksheet->write($ws_row, $#seq_data_headers, 'COUNT_ALLELES'); $ws_row++;
		$worksheet->write_row($ws_row, 0, \@seq_data_headers, $bold); $ws_row++;

		my $ws_row_first = $ws_row;
		my $ws_row_last = $ws_row_first+$#seq_md5s;
		my $ws_col = $#seq_data_headers;
		my $ws_col_first = $ws_col+1;

		my @samples = @{$results->{$marker_name}{'samples'}};
		foreach my $sample (@samples){

			$ws_row = 0; $ws_col++;
			$worksheet->write($ws_row, $ws_col, $results->{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'}); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $results->{$marker_name}{'sample_data'}{$sample}{'depth_alleles'}); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $results->{$marker_name}{'sample_data'}{$sample}{'count_alleles'}); $ws_row++;
			$worksheet->write($ws_row, $ws_col, $sample, $bold); $ws_row++;

			foreach my $seq_md5 (@seq_md5s) {

				if (defined($results->{$marker_name}{'assignments'}{$seq_md5}{$sample})){
					$worksheet->write($ws_row, $ws_col, $results->{$marker_name}{'assignments'}{$seq_md5}{$sample}); $ws_row++;
					$final_seq_depths{$seq_md5} += $results->{$marker_name}{'assignments'}{$seq_md5}{$sample};
					$final_seq_samples{$seq_md5}++;
				} else {
					$worksheet->write($ws_row, $ws_col, ''); $ws_row++;
				}
			}
		}

		my $ws_col_last = $ws_col;

		# Writes seq information
		$ws_row = $ws_row_first;
		foreach my $seq_md5 (@seq_md5s) {
			$ws_col = 0;
			$worksheet->write($ws_row, $ws_col,$md5_to_sequence{$seq_md5}); $ws_col++;
			$worksheet->write($ws_row, $ws_col,$seq_md5); $ws_col++;
			$worksheet->write($ws_row, $ws_col,length($md5_to_sequence{$seq_md5})); $ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $final_seq_depths{$seq_md5}); $ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $final_seq_samples{$seq_md5}); $ws_col++;
			$worksheet->write($ws_row, $ws_col,$md5_to_name{$seq_md5}); $ws_col++;
			$ws_row++;
		}

	}

	$workbook->close();
	
	return $file;

}

#################################################################################

# Writes AmpliCDR3 results into an Excel file
sub write_amplicdr3_file_results {

	my ($file,$results,$properties) = @_;

	if (!defined($properties)){
		$properties = {
			'title'    => "AmpliCDR3 results file",
			'author'   => "Alvaro Sebastian",
			'comments' => "AmpliSAS results file",
			'company'  => "SixthResearcher",
		};
	}

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

	foreach my $marker_name (keys %{$results}){

		# Obtains all the variants in all the files
		my (%seq_depths,%final_seq_depths,%final_seq_samples);
		my (%md5_to_sequence, %md5_to_name);
		if (defined($results->{$marker_name})){
			foreach my $md5 (keys %{$results->{$marker_name}{'seq_data'}}){
				$seq_depths{$md5} += $results->{$marker_name}{'seq_data'}{$md5}{'depth'};
				$md5_to_sequence{$md5} = $results->{$marker_name}{'seq_data'}{$md5}{'sequence'};
				$md5_to_name{$md5} = $results->{$marker_name}{'seq_data'}{$md5}{'name'};
			}
		}
		if (!%seq_depths) { next; }

		my @seq_md5s = sort { $seq_depths{$b} <=> $seq_depths{$a} } keys %seq_depths;
		
		# Creates worksheet and writes data headers
		my $worksheet = $workbook->add_worksheet("$marker_name");
		$worksheet->set_column('F:F', undef, $bold);
		my $ws_row = 0;
		my @seq_data_headers = ('SEQUENCE', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'NAME');
		$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_AMPLICON'); $ws_row++;
		$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_ALLELES'); $ws_row++;
		$worksheet->write($ws_row, $#seq_data_headers, 'COUNT_ALLELES'); $ws_row++;
		$worksheet->write_row($ws_row, 0, \@seq_data_headers, $bold); $ws_row++;

		my $ws_row_first = $ws_row;
		my $ws_row_last = $ws_row_first+$#seq_md5s;
		my $ws_col = $#seq_data_headers;
		my $ws_col_first = $ws_col+1;

		my @samples = @{$results->{$marker_name}{'samples'}};
		foreach my $sample (@samples){

			$ws_row = 0; $ws_col++;
			$worksheet->write($ws_row, $ws_col, $results->{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'}); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $results->{$marker_name}{'sample_data'}{$sample}{'depth_alleles'}); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $results->{$marker_name}{'sample_data'}{$sample}{'count_alleles'}); $ws_row++;
			$worksheet->write($ws_row, $ws_col, $sample, $bold); $ws_row++;

			foreach my $seq_md5 (@seq_md5s) {

				if (defined($results->{$marker_name}{'assignments'}{$seq_md5}{$sample})){
					$worksheet->write($ws_row, $ws_col, $results->{$marker_name}{'assignments'}{$seq_md5}{$sample}); $ws_row++;
					$final_seq_depths{$seq_md5} += $results->{$marker_name}{'assignments'}{$seq_md5}{$sample};
					$final_seq_samples{$seq_md5}++;
				} else {
					$worksheet->write($ws_row, $ws_col, ''); $ws_row++;
				}
			}
		}

		my $ws_col_last = $ws_col;

		# Writes seq information
		$ws_row = $ws_row_first;
		foreach my $seq_md5 (@seq_md5s) {
			$ws_col = 0;
			$worksheet->write($ws_row, $ws_col,$md5_to_sequence{$seq_md5}); $ws_col++;
			$worksheet->write($ws_row, $ws_col,$seq_md5); $ws_col++;
			$worksheet->write($ws_row, $ws_col,length($md5_to_sequence{$seq_md5})); $ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $final_seq_depths{$seq_md5}); $ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $final_seq_samples{$seq_md5}); $ws_col++;
			$worksheet->write($ws_row, $ws_col,$md5_to_name{$seq_md5}); $ws_col++;
			$ws_row++;
		}

	}

	$workbook->close();
	
	return $file;

}

#################################################################################

# Writes AmpliHLA results into an Excel file
sub write_amplihla_file_results {

	my ($file,$results,$properties) = @_;
	
	if (!defined($properties)){
		$properties = {
			'title'    => "AmpliSAS format file",
			'author'   => "Alvaro Sebastian",
			'comments' => "AmpliSAS format file",
			'company'  => "SixthResearcher",
		};
	}

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
	my $decimal = $workbook->add_format( num_format => '[=0]0;0.##' );
	foreach my $hla_type (keys %{$results}){

		my @samples = sort {$a cmp $b} keys %{$results->{$hla_type}{'sample_data'}};

		# Obtains all the alleles from all the samples
		my (%allele_freqs, %allele_samples, $allele_freqs);
		foreach my $sample (@samples){
			foreach my $allele (keys %{$results->{$hla_type}{'sample_data'}{$sample}{'freq'}}) {
				push(@{$allele_freqs->{$allele}}, $results->{$hla_type}{'sample_data'}{$sample}{'freq'}{$allele});
				$allele_samples{$allele}++;
			}
		}
		if (!%$allele_freqs) { next; }
		foreach my $allele (keys %$allele_freqs) {
			$allele_freqs{$allele} = mean(@{$allele_freqs->{$allele}});
		}

# 		my @alleles = sort { $allele_samples{$b} <=> $allele_samples{$a} } keys %allele_samples;
		my @alleles = sort { $allele_samples{$b}*$allele_freqs{$b} <=> $allele_samples{$a}*$allele_freqs{$a} } keys %allele_freqs;

		# Creates worksheet and writes data headers
		my $worksheet = $workbook->add_worksheet("HLA-$hla_type");
		# my $worksheet1 = $workbook->add_worksheet("$hla_type\_freqs");
		# my $worksheet2 = $workbook->add_worksheet("$hla_type\_alleles");
		my $ws_row = 0;
		my @seq_data_headers = ('SEQUENCES', 'MD5S', 'MEAN_FREQ', 'SAMPLES', 'ALLELE');
		$worksheet->write($ws_row, $#seq_data_headers, 'COUNT_ALLELES', $bold);
		# $worksheet1->write($ws_row, $#seq_data_headers, 'COUNT_ALLELES', $bold);
		# $worksheet2->write($ws_row, $#seq_data_headers, 'COUNT_ALLELES', $bold);
		$ws_row++;
		$worksheet->write_row($ws_row, 0, \@seq_data_headers, $bold);
		# $worksheet1->write_row($ws_row, 0, \@seq_data_headers, $bold);
		# $worksheet2->write_row($ws_row, 0, \@seq_data_headers, $bold);
		$ws_row++;

		my $ws_row_first = $ws_row;
		my $ws_row_last = $ws_row_first+$#alleles;
		my $ws_col = $#seq_data_headers;
		my $ws_col_first = $ws_col+1;

		my %allele_ambiguities;
		foreach my $sample (@samples){
					
			$ws_row = 0; $ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, scalar keys %{$results->{$hla_type}{'sample_data'}{$sample}{'freq'}});
			# $worksheet1->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, scalar keys %{$results->{$hla_type}{'sample_data'}{$sample}{'freq'}});
			# $worksheet2->write_formula($ws_row, $ws_col, sprintf('=COUNTA(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, scalar keys %{$results->{$hla_type}{'sample_data'}{$sample}{'freq'}});
			$ws_row++;
			$worksheet->write($ws_row, $ws_col, $sample, $bold);
			# $worksheet1->write($ws_row, $ws_col, $sample, $bold);
			# $worksheet2->write($ws_row, $ws_col, $sample, $bold);
			$ws_row++;

			foreach my $allele (@alleles) {
				if (defined($results->{$hla_type}{'sample_data'}{$sample}{'freq'}{$allele})){
					$worksheet->write($ws_row, $ws_col, $results->{$hla_type}{'sample_data'}{$sample}{'freq'}{$allele},$decimal);
					# $worksheet1->write($ws_row, $ws_col, $results->{$hla_type}{'sample_data'}{$sample}{'freq'}{$allele},$decimal);
					# $worksheet2->write($ws_row, $ws_col, join(' | ', @{$results->{$hla_type}{'sample_data'}{$sample}{'alleles'}{$allele}}));
					$ws_row++;
					# Checks if individual allele ambiguities are the same for all the samples
					my @ambiguities = sort {$a cmp $b} @{$results->{$hla_type}{'sample_data'}{$sample}{'alleles'}{$allele}};
					if ($#ambiguities > 0) {
						my $ambiguities = join (", ",@ambiguities);
						if (!defined ($allele_ambiguities{$allele})){
							$allele_ambiguities{$allele} = $ambiguities;
						# If two grouped alleles have diferent ambiguities, ej. A*01 -> A*0101, A*0104 & A*01 -> A*0101, A*0107
						} elsif ($allele_ambiguities{$allele} ne $ambiguities){
							$allele_ambiguities{"$allele ($sample)"} = $ambiguities;
						}
					}
# 					if ($#{$results->{$hla_type}{'sample_data'}{$sample}{'alleles'}{$allele}} > 0 ) {
# 						$allele_ambiguities->{$allele} = [ unique(@{$allele_ambiguities->{$allele}},@{$results->{$hla_type}{'sample_data'}{$sample}{'alleles'}{$allele}}) ];
# 					}
				} else {
					$worksheet->write($ws_row, $ws_col, ' ');
					# $worksheet1->write($ws_row, $ws_col, '');
					# $worksheet2->write($ws_row, $ws_col, '');
					$ws_row++;
				}
			}
		}

		my $ws_col_last = $ws_col;

		# Writes seq information
		$ws_row = $ws_row_first;
		foreach my $allele (@alleles) {
			$ws_col = 0;
			my @allele_md5s = map "$_: ".$results->{$hla_type}{'allele_data'}{$allele}{'md5'}{$_} , keys %{$results->{$hla_type}{'allele_data'}{$allele}{'md5'}};
			my @allele_seqs = map "$_: ".$results->{$hla_type}{'allele_data'}{$allele}{'seq'}{$_} , keys %{$results->{$hla_type}{'allele_data'}{$allele}{'seq'}};
			$worksheet->write($ws_row, $ws_col, join("\n",@allele_seqs));
			# $worksheet1->write($ws_row, $ws_col, join("\n",@allele_seqs));
			# $worksheet2->write($ws_row, $ws_col, join("\n",@allele_seqs));
			$ws_col++;
			$worksheet->write($ws_row, $ws_col, join("\n",@allele_md5s));
			# $worksheet1->write($ws_row, $ws_col, join("\n",@allele_md5s));
			# $worksheet2->write($ws_row, $ws_col, join("\n",@allele_md5s));
			$ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=AVERAGE(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), $decimal, $allele_freqs{$allele});
			# $worksheet1->write_formula($ws_row, $ws_col, sprintf('=AVERAGE(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), $decimal, $allele_freqs{$allele});
			# $worksheet2->write($ws_row, $ws_col, $allele_freqs{$allele}, $decimal);
			$ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $allele_samples{$allele});
			# $worksheet1->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $allele_samples{$allele});
			# $worksheet2->write_formula($ws_row, $ws_col, sprintf('=COUNTA(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $allele_samples{$allele});
			$ws_col++;
			$worksheet->write($ws_row, $ws_col,$allele, $bold);
			# $worksheet1->write($ws_row, $ws_col,$allele, $bold);
			# $worksheet2->write($ws_row, $ws_col,$allele, $bold);
			$ws_col++;
			$ws_row++;
		}
		$ws_col = 4;
		$ws_row++;
		$worksheet->write($ws_row, $ws_col,"ALLELE", $bold);
		$worksheet->write($ws_row, $ws_col+1,"AMBIGUITIES", $bold);
		$ws_row++;
		foreach my $allele (nsort(keys %allele_ambiguities)){
			$worksheet->write($ws_row, $ws_col,"$allele", $bold);
			$worksheet->write($ws_row, $ws_col+1,join(', ',$allele_ambiguities{$allele}));
			$ws_row++;
		}
	}

	$workbook->close();
	
	return $file;

}

#################################################################################

# Writes AmpliCANCER results into an Excel file
sub write_amplicancer_file_results {

	my ($file,$results,$properties) = @_;
	
	if (!defined($properties)){
		$properties = {
			'title'    => "AmpliSAS format file",
			'author'   => "Alvaro Sebastian",
			'comments' => "AmpliSAS format file",
			'company'  => "SixthResearcher",
		};
	}

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

	# SAME THAN AMPLISAS:

	# Should be only 1 marker name (the name of the input file)
	foreach my $marker_name (keys %{$results}){

		# Obtains all the variants in all the files
		my (%seq_depths,%final_seq_depths,%final_seq_samples);
		my (%md5_to_sequence, %md5_to_name);
		if (defined($results->{$marker_name})){
			foreach my $md5 (keys %{$results->{$marker_name}{'seq_data'}}){
				$seq_depths{$md5} += $results->{$marker_name}{'seq_data'}{$md5}{'depth'};
				$md5_to_sequence{$md5} = $results->{$marker_name}{'seq_data'}{$md5}{'sequence'};
				$md5_to_name{$md5} = $results->{$marker_name}{'seq_data'}{$md5}{'name'};
			}
		}
		if (!%seq_depths) { next; }

		my @seq_md5s = sort { $seq_depths{$b} <=> $seq_depths{$a} } keys %seq_depths;
		
		# Creates worksheet and writes data headers
		my $worksheet = $workbook->add_worksheet("variants");
		$worksheet->set_column('F:F', undef, $bold);
		my $ws_row = 0;
		my @seq_data_headers = ('SEQUENCE', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'NAME');
		$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_AMPLICON'); $ws_row++;
		$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_ALLELES'); $ws_row++;
		$worksheet->write($ws_row, $#seq_data_headers, 'COUNT_ALLELES'); $ws_row++;
		$worksheet->write_row($ws_row, 0, \@seq_data_headers, $bold); $ws_row++;

		my $ws_row_first = $ws_row;
		my $ws_row_last = $ws_row_first+$#seq_md5s;
		my $ws_col = $#seq_data_headers;
		my $ws_col_first = $ws_col+1;

		my @samples = @{$results->{$marker_name}{'samples'}};
		foreach my $sample (@samples){

			$ws_row = 0; $ws_col++;
			$worksheet->write($ws_row, $ws_col, $results->{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'}); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $results->{$marker_name}{'sample_data'}{$sample}{'depth_alleles'}); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $results->{$marker_name}{'sample_data'}{$sample}{'count_alleles'}); $ws_row++;
			$worksheet->write($ws_row, $ws_col, $sample, $bold); $ws_row++;

			foreach my $seq_md5 (@seq_md5s) {

				if (defined($results->{$marker_name}{'assignments'}{$seq_md5}{$sample})){
					$worksheet->write($ws_row, $ws_col, $results->{$marker_name}{'assignments'}{$seq_md5}{$sample}); $ws_row++;
					$final_seq_depths{$seq_md5} += $results->{$marker_name}{'assignments'}{$seq_md5}{$sample};
					$final_seq_samples{$seq_md5}++;
				} else {
					$worksheet->write($ws_row, $ws_col, ''); $ws_row++;
				}
			}
		}

		my $ws_col_last = $ws_col;

		# Writes seq information
		$ws_row = $ws_row_first;
		foreach my $seq_md5 (@seq_md5s) {
			$ws_col = 0;
			$worksheet->write($ws_row, $ws_col,$md5_to_sequence{$seq_md5}); $ws_col++;
			$worksheet->write($ws_row, $ws_col,$seq_md5); $ws_col++;
			$worksheet->write($ws_row, $ws_col,length($md5_to_sequence{$seq_md5})); $ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $final_seq_depths{$seq_md5}); $ws_col++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $final_seq_samples{$seq_md5}); $ws_col++;
			$worksheet->write($ws_row, $ws_col,$md5_to_name{$seq_md5}); $ws_col++;
			$ws_row++;
		}

	}


	# MUTATION ANNOTATION (ONE SHEET PER SAMPLE):
	
	# Should be only 1 marker name (the name of the input file)
	foreach my $marker_name (keys %{$results}){

		# Creates a list with all the samples
		my @samples = @{$results->{$marker_name}{'samples'}};
		# Creates a list  with all the UniProt IDs ordered by number of mutations
		my @uniprot_ids = sort { scalar keys %{$results->{$marker_name}{'prot_data'}{$b}{'mutations'}} <=> scalar keys %{$results->{$marker_name}{'prot_data'}{$a}{'mutations'}} } keys %{$results->{$marker_name}{'prot_data'}};

		foreach my $sample (@samples){

			# Creates a worksheet and writes protein errors
			my $worksheet = $workbook->add_worksheet($sample);
			$worksheet->set_column('E:E', undef, $bold);
			my $ws_row = 0;
			my @seq_data_headers = ('UNIPROT', 'DESCRIPTION', 'NAME');
			$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_AMPLICON'); $ws_row++;
			$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_REFERENCES'); $ws_row++;
			$worksheet->write($ws_row, $#seq_data_headers, 'COUNT_REFERENCES'); $ws_row++;
			$worksheet->write_row($ws_row, 0, \@seq_data_headers, $bold); $ws_row++;

			my $ws_row_first = $ws_row;
			my $ws_col = $#seq_data_headers+1;
			my $ws_col_first = $ws_col;

			my @sample_uniprot_ids;
			my $sample_uniprot_depth;
			foreach my $uniprot_id (@uniprot_ids) {
				my $depth = 0;
				my @mutations;
				foreach my $seq_md5 (@{$results->{$marker_name}{'prot_data'}{$uniprot_id}{'variants'}}){
					if (defined($results->{$marker_name}{'assignments'}{$seq_md5}{$sample})){
						$depth += $results->{$marker_name}{'assignments'}{$seq_md5}{$sample};
						push(@mutations, @{$results->{$marker_name}{'prot_data'}{$uniprot_id}{'mutations'}{$seq_md5}});
					}
				}
				if ($depth) {
					if (@mutations) {
						$worksheet->write($ws_row, $ws_col, sprintf("%s; %d reads", join(", ",@mutations), $depth), $red);
					} else {
						$worksheet->write($ws_row, $ws_col, sprintf("%d reads", $depth));
					}
					push(@sample_uniprot_ids, $uniprot_id);
					$sample_uniprot_depth += $depth;
				} else {
					# $worksheet->write($ws_row, $ws_col, '');
					next;
				}
				$ws_row++;
			}

			my $ws_col_last = $ws_col;
			my $ws_row_last = $ws_row_first+$#sample_uniprot_ids;

			# Writes sample names, and total coverage and number of matched protein references
			$ws_row = 0;
			$ws_col=$ws_col_first;
			$worksheet->write($ws_row, $ws_col, $results->{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'}); $ws_row++;
			$worksheet->write($ws_row, $ws_col, $sample_uniprot_depth); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, sprintf('%d',scalar @sample_uniprot_ids)); $ws_row++;
			$worksheet->write($ws_row, $ws_col, $sample, $bold); $ws_row++;

			# Writes seq information
			$ws_row = $ws_row_first;
			foreach my $uniprot_id (@sample_uniprot_ids) {
				$ws_col = 0;
				$worksheet->write($ws_row, $ws_col,$uniprot_id); $ws_col++;
				$worksheet->write($ws_row, $ws_col,$results->{$marker_name}{'prot_data'}{$uniprot_id}{'description'}); $ws_col++;
				$worksheet->write($ws_row, $ws_col,$results->{$marker_name}{'prot_data'}{$uniprot_id}{'name'}); $ws_col++;
				$ws_row++;
			}
		}

	}


	$workbook->close();
	
	return $file;

}

#################################################################################

# Writes AmpliTXO results into an OTU Table file
sub write_amplitaxo_otutable {

	my ($file,$results,$options) = @_;
	
	return write_to_file($file,create_amplitaxo_otutable($results,$options));
	
}

#################################################################################

# Creates an OTU Table from AmpliTAXO results
sub create_amplitaxo_otutable {

	my ($results,$options) = @_;

	my $otu_table = '';

	my ($format, $type) = ('txt','depth');
	if (in_array($options,'biom')){
		$format = 'biom';
	}
	if (in_array($options,'freq') || in_array($options,'frequency')){
		$type = 'freq';
	}

	# Obtains an array with the total unique samples
	my @unique_samples;
	foreach my $marker_name (keys %{$results}){
		@unique_samples = unique(@unique_samples,@{$results->{$marker_name}{'samples'}});
	}
	@unique_samples = nsort(@unique_samples);

	my $row = 0;
	my (@txt_matrix,@biom_taxa,@biom_samples,@biom_matrix);
	foreach my $marker_name (keys %{$results}){

		# Obtains all the variants in all the files
		my (%otu_depths,%md5_to_name);
		if (defined($results->{$marker_name})){
			foreach my $md5 (keys %{$results->{$marker_name}{'seq_data'}}){
				$otu_depths{$md5} += $results->{$marker_name}{'seq_data'}{$md5}{'depth'};
				$md5_to_name{$md5} = $results->{$marker_name}{'seq_data'}{$md5}{'name'};
			}
		}
		if (!%otu_depths) { next; }

		my @otu_md5s = sort { $otu_depths{$b} <=> $otu_depths{$a} } keys %otu_depths;

		foreach my $otu_md5 (@otu_md5s) {
			# tax=k:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Lactobacillaceae,g:Lactobacillus,s:reuteri;
			my $taxo = group_taxo($results->{$marker_name}{'seq_data'}{$otu_md5}{'taxo_assigns'});
			my @taxo_data;
			while ($taxo =~ /([kpcofgs])(\:|__)(\w+)/g){
				push(@taxo_data, sprintf("\"%s__%s\"",$1,$3));
			}
			push(@biom_taxa,sprintf("{\"id\": \"%s\", \"metadata\": {\"taxonomy\": [%s]}}", $md5_to_name{$otu_md5}, join(", ",@taxo_data)));

			my $col = 0;
			my @row_values;
			foreach my $sample (@unique_samples){
				my $value;
				if ($type eq 'freq'){
					$value = sprintf("%.4f",$results->{$marker_name}{'assignments'}{$otu_md5}{$sample} / $results->{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'});
				} else {
					$value = $results->{$marker_name}{'assignments'}{$otu_md5}{$sample};
				}
				if ($format eq 'biom'){
					push(@biom_matrix,sprintf("[%d,%d,%d]",$row,$col,$value));
				} elsif($format eq 'txt'){
					push(@row_values,$value);
				}
				$col++;
			}
			if ($format eq 'txt'){
				push(@txt_matrix, join("\t",$md5_to_name{$otu_md5},@row_values,$taxo));
			}
			$row++;
		}
	}

	if ($format eq 'biom'){
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
		my $date = sprintf("%d-%02d-%02d %02d:%02d:%02d",$year,$mon,$mday,$hour,$min,$sec);
		my $nrows = scalar @biom_taxa;
		my $ncols = scalar @unique_samples;
		@biom_samples = map '{"id": "'.$_.'", "metadata": null}', @unique_samples;
		$otu_table .=	'{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org",'.
				'"generated_by": "AmpliTAXO","date": "'.$date.'",'.
				'"matrix_element_type": "float","shape": ['.$nrows.', '.$ncols.'],"type": "OTU table","matrix_type": "sparse",';
		$otu_table .= '"data": ['.join(",",@biom_matrix).'],';
		$otu_table .= '"rows": ['.join(",", @biom_taxa).'],';
		$otu_table .= '"columns": ['.join(",", @biom_samples).']';
		$otu_table .= '}';
	} elsif($format eq 'txt'){
		$otu_table .= join("\t",'#OTU ID',@unique_samples,'Taxonomy')."\n";
		$otu_table .= join("\n",@txt_matrix)."\n";
	}
	
	return $otu_table;

}

#################################################################################

# # Creates an OTU Table from AmpliTAXO results
# sub create_amplitaxo_otutable {
# 
# 	my ($results,$options) = @_;
# 
# 	my $otu_table = '';
# 
# 	my ($format, $type) = ('txt','depth');
# 	if (in_array($options,'biom')){
# 		$format = 'biom';
# 	}
# 	if (in_array($options,'freq') || in_array($options,'frequency')){
# 		$type = 'freq';
# 	}
# 
# 	my @samples = @{$results->{'TAXO'}{'samples'}};
# 
# 	# Obtains all the taxa from all the samples
# 	my (%taxo_depths,%taxo_names);
# 	foreach my $taxo (keys %{$results->{'TAXO'}{'assignments'}}) {
# 		foreach my $sample (keys %{$results->{'TAXO'}{'assignments'}{$taxo}}) {
# 			$taxo_depths{$taxo} += $results->{'TAXO'}{'assignments'}{$taxo}{$sample};
# 			$taxo_names{$taxo} = join('|',@{$results->{'TAXO'}{'taxo_data'}{$taxo}{'otu_names'}});
# 		}
# 	}
# 
# 	my @taxa = sort { $taxo_depths{$b} <=> $taxo_depths{$a} } keys %taxo_depths;
# 
# 	if ($format eq 'biom'){
# 		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
# 		my $date = sprintf("%d-%02d-%02d %02d:%02d:%02d",$year,$mon,$mday,$hour,$min,$sec);
# 		my $nrows = scalar @taxa;
# 		my $ncols = scalar @samples;
# 		$otu_table .=	'{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org",'.
# 				'"generated_by": "AmpliTAXO","date": "'.$date.'",'.
# 				'"matrix_element_type": "float","shape": ['.$nrows.', '.$ncols.'],"type": "OTU table","matrix_type": "sparse",';
# 	} elsif($format eq 'txt'){
# 		$otu_table .= join("\t",'#OTU ID','Taxonomy',@samples)."\n";
# 	}
# 
# 	my $row = 0;
# 	my (@biom_taxa,@biom_samples,@biom_matrix);
# 	foreach my $taxo (@taxa) {
# 		my @row_values;
# 		# tax=k:Bacteria,p:Firmicutes,c:Bacilli,o:Lactobacillales,f:Lactobacillaceae,g:Lactobacillus,s:reuteri;
# 		my @taxa;
# 		while ($taxo =~ /([kpcofgs])(\:|__)(\w+)/g){
# 			push(@taxa, sprintf("\"%s__%s\"",$1,$3));
# 		}
# 		push(@biom_taxa,sprintf("{\"id\": \"%s\", \"metadata\": {\"taxonomy\": [%s]}}", $taxo_names{$taxo}, join(", ",@taxa)));
# 		my $col = 0;
# 		foreach my $sample (@samples){
# 			my $value = sprintf("%d",$results->{'TAXO'}{'assignments'}{$taxo}{$sample});
# 			if ($type eq 'freq'){
# 				$value = sprintf("%.4f",$results->{'TAXO'}{'assignments'}{$taxo}{$sample}/$results->{'TAXO'}{'sample_data'}{$sample}{'depth_total'});
# 			}
# 			if ($format eq 'biom'){
# 				push(@biom_matrix,sprintf("[%d,%d,%d]",$row,$col,$value));
# 			} elsif($format eq 'txt'){
# 				push(@row_values,$value);
# 			}
# 			$col++;
# 		}
# 		if ($format eq 'txt'){
# 			$otu_table .= join("\t",$taxo_names{$taxo},$taxo,@row_values)."\n";
# 		}
# 		$row++;
# 	}
# 	@biom_samples = map '{"id": "'.$_.'", "metadata": null}', @samples;
# 	if ($format eq 'biom'){
# 		$otu_table .= '"data": ['.join(",",@biom_matrix).'],';
# 		$otu_table .= '"rows": ['.join(",", @biom_taxa).'],';
# 		$otu_table .= '"columns": ['.join(",", @biom_samples).']';
# 		$otu_table .= '}';
# 	}
# 	
# 	return $otu_table;
# 
# }

#################################################################################

# Estimates the number of alleles by the Degree of Change method  (Lighten et al. 2014)
# Takes as input an array of unique sequence depths or frequencies and the max. number of expected alleles
sub degree_of_change {

	my ($depths, $max_alleles) = @_;

	# Default max. number of expected alleles
	if (!defined($max_alleles)){
		$max_alleles = 10
	}

	my $max_DOCn = 0;

	# Sorts non zero values from the array (includes only the first zero value)
	my @sorted_depths = sort {$b<=>$a} grep {$_ != 0} @$depths;

	my $i_max = $max_alleles;
	# If there are less depths than max. number of expected alleles
	if ($#sorted_depths<$i_max) {
		$i_max = $#sorted_depths;
	}

	# Calculates cumulative depths (only neccesary for some graphs)
	my @cum_depths = ($sorted_depths[0]);
	for (my $i=1; $i<=$i_max; $i++){
		$cum_depths[$i] = $cum_depths[$i-1]+$sorted_depths[$i];
	}

	my @DOCs;
	# Calculates DOCs = ROC[$i]/ROC[$i+1] = depth[$i]/depth[$i+1])
	for (my $i=0; $i<$i_max; $i++){
		push(@DOCs, $sorted_depths[$i]/$sorted_depths[$i+1]);
	}
	my $sum_DOCs = sum(@DOCs);

	# Calculates DOCns = DOC/sum_DOCs*100
	my @DOCns;
	for (my $i=0; $i<=$#DOCs; $i++){
		my $DOCn = $DOCs[$i]/$sum_DOCs*100;
		push(@DOCns, $DOCn);
		#printf("DOC-%03d = %.2f\n", $i+1, $DOCn);
	}

	# Finds the max_DOCn position (max. number of alleles)
	$DOCns[$max_DOCn] > $DOCns[$_] or $max_DOCn = $_ for 1 .. $#DOCns;

	return ($max_DOCn+1, \@DOCns, \@cum_depths);

}

#################################################################################

# Generates a random barcode sequence of the given length
sub generate_barcode {

	my ($length,$previous_barcodes) = @_;
	
	my @nts = ('A', 'C', 'G', 'T');
	
	my ($barcode,@barcode);
	my $right_barcode = 0;
	while (!$right_barcode) {
		@barcode = map $nts[rand @nts] , 1 .. $length;
		$barcode = join('',@barcode);
		$right_barcode = 1;
		 # No more than 2 consecutive identical nucleotides
		if ($barcode =~ /AA+|CC+|GG+|TT+/) {
			$right_barcode = 0;
			next;
		}
		foreach my $previous_barcode (@$previous_barcodes){
			my @previous_barcode = split('',$previous_barcode);
			my $simil = 0;
			for (my $i=0; $i<=$#barcode; $i++) {
				if ($barcode[$i] eq $previous_barcode[$i]) {
					$simil++;
				}
			}
			if ($length-$simil<2) { # At least 2nts of difference between barcodes
				$right_barcode = 0;
				last;
			}
		}
	}
	
	return $barcode;

}

#################################################################################

# Retrieves the sequencing technology used and sets automatically analysis parameters
sub adjust_automatic_params {

	my ($markers,$samples,$marker_seq_data,$amplicon_seq_data,$options) = @_;

	my ($technology, $marker_lengths, $custom_params, $major_errors);

	my ($only_lengths)=(0);
	if (in_array($options, 'only lengths')){
		$only_lengths = 1;
	}

	# Reads the major variants to annotates the error rates
	foreach my $marker_name (@$markers){
		if (!defined($amplicon_seq_data->{$marker_name})){
			next;
		}
		my $marker_total_depth = 0;
		my (%marker_major_len_depths,$marker_min_dominant_freq, %marker_all_len_depths);
		foreach my $sample_name (@{$samples->{$marker_name}}) {
			if (!defined($amplicon_seq_data->{$marker_name}{$sample_name})){
				next;
			}
			# Sorts amplicon sequences by depth
			my @sorted_depth_amplicon_unique_seqs = sort { $amplicon_seq_data->{$marker_name}{$sample_name}{$b}{'depth'} <=> $amplicon_seq_data->{$marker_name}{$sample_name}{$a}{'depth'} } keys %{$amplicon_seq_data->{$marker_name}{$sample_name}};
			# Reads the lengths of the first 10 variants in the amplicon (sorted by depth)
# 			my $seq_limit = scalar @sorted_depth_amplicon_unique_seqs;
# 			if ($seq_limit > 10){
# 				$seq_limit = 10;
# 			}
			my $amplicon_total_depth = sum( map { $amplicon_seq_data->{$marker_name}{$sample_name}{$_}{'depth'} } keys %{$amplicon_seq_data->{$marker_name}{$sample_name}} );
			$marker_total_depth += $amplicon_total_depth;
			#Analyzes only variants that accumulate the 25% of amplicon depth, with a limit of 10 variants
			my (@amplicon_major_variants, @amplicon_major_seqs);
			my $amplicon_cum_depth = 0;
			foreach my $md5 (@sorted_depth_amplicon_unique_seqs) {
				my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
				my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
				$marker_all_len_depths{$len} += $depth;
				if ($amplicon_cum_depth<0.25*$amplicon_total_depth && scalar @amplicon_major_variants < 10 ){
					push(@amplicon_major_variants,$md5);
					push(@amplicon_major_seqs,$marker_seq_data->{$marker_name}{$md5}{'seq'});
					$marker_major_len_depths{$len} += $depth;
					$amplicon_cum_depth += $depth;
				}
			}
			if ($only_lengths) {
				next;
			}
			foreach my $md5 (@amplicon_major_variants) {
# if ($md5 eq 'fd228dc8cac592f8dd00a2bc1a0c364d'){
# print '';
# }
				my $depth = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'depth'};
				my $len = $marker_seq_data->{$marker_name}{$md5}{'len'};
				my $seq = $marker_seq_data->{$marker_name}{$md5}{'seq'};
				my $freq = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5}{'freq'};
				# Aligns the most abundant variant against the lowest frequency ones to detect the class of sequencing errors
				if ($md5 eq $amplicon_major_variants[0]){
					# Reverse the order of variants (from lower to higher freq)
					foreach my $md5_2 (reverse(@sorted_depth_amplicon_unique_seqs)) {
# if ($sample_name eq 'S5'){
# print '';
# }
						my $len_2 = $marker_seq_data->{$marker_name}{$md5_2}{'len'};
						my $freq_2 = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5_2}{'freq'};
						# Check only if they are similar in length
						# if (abs($len-$len_2)/$len>0.01){
						if (abs($len-$len_2)/$len>0.02){ # Works for reads longer than 50bp
							next;
						}
						# Check only if they are different in frequency
						if (abs($freq-$freq_2)/$freq<0.9){
							next;
						}
						# Check only if they are different variants
						if ($md5 eq $md5_2) {
							next;
						}
						my $seq_2 = $marker_seq_data->{$marker_name}{$md5_2}{'seq'};
						my @aligned_seqs = align_2seqs($seq,$seq_2, 'needle');
						# Skips if the pairwise global alignment fails
						if (!defined($aligned_seqs[0]) || !defined($aligned_seqs[1])){
							next;
						}
						my ($identical,$total) = binary_score_nts(@aligned_seqs);
						# Checks only if they are similar enough to be sequencing errors
						# if ($identical/$total < 0.99 ) {
						if ($identical/$total < 0.98 ) { # Works for reads longer than 50bp
							next;
						}
						my ($substitutions, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors(@aligned_seqs);
# 						if (scalar @$insertions + scalar @$deletions > scalar @$substitutions) {
							push(@{$major_errors->{'indels'}},sprintf('%.2f', 100 * (scalar @$insertions + scalar @$deletions)/$total));
# 						} elsif (scalar @$substitutions > scalar @$insertions + scalar @$deletions) {
							push(@{$major_errors->{'substitutions'}},sprintf('%.2f', 100 * scalar @$substitutions/$total));
# 						}
						# Annotates only the first low freq variant that looks like an artifact of the major one
						last;
					}
				}
				# Aligns high freq variants to annotate the freq diffences between highly disimilar variants
				my $aligned_seqs = align_seqs2one($seq,$md5,\@amplicon_major_seqs,\@amplicon_major_variants,'needleall'); #time=20% # 104 seg
				foreach my $md5_2 (@amplicon_major_variants) {
					if ($md5 eq $md5_2) { next; }
					my $seq_2 = $marker_seq_data->{$marker_name}{$md5_2}{'seq'};
# 					my @aligned_seqs = align_2seqs($seq,$seq_2, 'needle');
					# Skips if the pairwise global alignment fails
					if (!defined($aligned_seqs->{$md5_2}[0]) || !defined($aligned_seqs->{$md5_2}[1])){
						next;
					}
					my ($identical,$total) = binary_score_nts(@{$aligned_seqs->{$md5_2}});
					# Checks only if they are disimilar enough not to be sequencing errors
					# if ($identical/$total >= 0.99 ) {
					if ($identical/$total >= 0.95 ) {
						next;
					}
					my $freq_2 = $amplicon_seq_data->{$marker_name}{$sample_name}{$md5_2}{'freq'};
					my $freq_ratio = sprintf('%.2f', 100 * $freq_2/$freq);
					# If frequency ratio respect to the major is lower than previous one, annotates it
					if (!defined($marker_min_dominant_freq) || $freq_ratio < $marker_min_dominant_freq){
						$marker_min_dominant_freq = $freq_ratio;
					}
				}
			}
		}
		# Takes the most abundant lengths as the correct amplicon lengths
		my @sorted_major_lengths = sort { $marker_major_len_depths{$b} <=> $marker_major_len_depths{$a} } keys %marker_major_len_depths;
		for (my $i=0; $i<=$#sorted_major_lengths; $i++){
			my $len = $sorted_major_lengths[$i];
			# If there are several lengths with high coverage (>75% of the previous one)
			if ($i==0 || $marker_major_len_depths{$len} > 0.75*$marker_major_len_depths{$marker_lengths->{$marker_name}[-1]}) {
				push(@{$marker_lengths->{$marker_name}}, $len);
			# Activates in-frame option if the major lengths are in frame
			} elsif ($only_lengths && abs($sorted_major_lengths[0] - $sorted_major_lengths[$i]) % 3 == 0) {
				$custom_params->{'cluster_inframe'}{$marker_name}[0] = 1;
			}
		}
		# Activates 'exact length' option if the major lengths have more than 95% of the total depth
		my $major_lengths_total_depth = 0;
		map $major_lengths_total_depth += $marker_all_len_depths{$_}, @{$marker_lengths->{$marker_name}};
		if ($major_lengths_total_depth >= 0.95*$marker_total_depth) {
			$custom_params->{'cluster_exact_length'}{$marker_name}[0] = 1;
		}

		# Annotates the minimum frequency respect to the dominant to consider new clusters of similar sequences
		if (defined($marker_min_dominant_freq)){
			if ($marker_min_dominant_freq<50) {
				$custom_params->{'min_dominant_frequency_threshold'}{$marker_name}[0] = $marker_min_dominant_freq;
			} else {
				$custom_params->{'min_dominant_frequency_threshold'}{$marker_name}[0] = 50;
			}
		}

	}
	

	# Chooses the sequencing technology used based on the kind of sequencing errors
	# And establishes indel and substitution error frequency parameters
	if (!$only_lengths) {
		my %error_stats;
		if (defined($major_errors->{'indels'})){
			$error_stats{'avg_indel_rate'} = sprintf('%.2f',mean(@{$major_errors->{'indels'}}));
			$error_stats{'max_indel_rate'} = max(@{$major_errors->{'indels'}});
			# $error_stats{'min_indel_rate'} = min(@{$major_errors->{'indels'}});
			$custom_params->{'indel_threshold'}{'all'}[0] = sprintf('%.2f', 1.1 * $error_stats{'max_indel_rate'});
		} else {
			$custom_params->{'indel_threshold'}{'all'}[0] = $error_stats{'avg_indel_rate'}
			= $error_stats{'max_indel_rate'} = $error_stats{'min_indel_rate'} = 0.001;
		}
		if (defined($major_errors->{'substitutions'})){
			$error_stats{'avg_subst_rate'} = sprintf('%.2f',mean(@{$major_errors->{'substitutions'}}));
			$error_stats{'max_subst_rate'} = max(@{$major_errors->{'substitutions'}});
			# $error_stats{'min_subst_rate'} = min(@{$major_errors->{'substitutions'}});
			$custom_params->{'substitution_threshold'}{'all'}[0] = sprintf('%.2f', 1.1 * $error_stats{'max_subst_rate'});
		} else {
			$custom_params->{'substitution_threshold'}{'all'}[0] = $error_stats{'avg_subst_rate'}
			= $error_stats{'max_subst_rate'} = $error_stats{'min_subst_rate'} = 0.001;
		}
		if ($error_stats{'avg_indel_rate'} > $error_stats{'avg_subst_rate'}){
			$technology = '454/IonTorrent';
		} elsif ($error_stats{'avg_subst_rate'} > 100*$error_stats{'avg_indel_rate'}){
			$technology = 'Illumina';
			# Do not use 'cluster_inframe' parameter with Illumina data
			delete($custom_params->{'cluster_inframe'});
		} else {
			$technology = 'Unknown';
		}
		return ($technology, $marker_lengths, $custom_params);
	} else {
		return $marker_lengths;
	}

}

#################################################################################



1;
