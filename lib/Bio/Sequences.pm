package Bio::Sequences;

######################################################################################################
#
# Bio::Sequences - Package of subroutines to work with DNA sequences, NGS data and genotyping
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
# Package of subroutines to work with DNA sequences, NGS data and genotyping
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
use threads;
use threads::shared;
use Scalar::Util 'looks_like_number';
use List::Util qw(shuffle);
use Digest::MD5 qw(md5);
use MIME::Base64;
use File::Basename;
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Compress::Zip qw(zip $ZipError);
use File::Find qw(find);
use Cwd qw(cwd);
# use Data::Dumper;
use Storable;
use 5.004;
use strict;
use warnings;

use Exporter;
 
our @ISA       = qw(Exporter);
our $VERSION   = '1.0';
our @EXPORT = qw(
trim
mean
median
max
min
sum
is_numeric
is_percentage
log10
variance
stdeviation
correlation
covariance
percentile
chao1
variations
binomial_probability
upgma
in_array
most_frequent
array_to_hash
extract_pattern_from_array
compare_arrays
compare_strings
unique
intersect
read_from_file
write_to_file
download_file
is_compressed
is_gzip
is_zip
is_bzip2
is_tar
is_tgz
is_tbz2
list_zip
list_gzip
list_tar
list_tgz
list_tbz2
is_xlsx
is_multifile
extract_paired_read_files_from_multifiles
extract_paired_read_files_from_path
extract_paired_read_files
read_files_from_path
compress
decompress
random_file_name
is_fasta
is_fastq
is_sam
is_bam
is_bcf
is_vcf
is_txt
is_dna
is_rna
is_prot
ngs_tech
read_fasta
read_fastq
read_sequence_file
read_fasta_file
read_fastq_file
read_gff_file
read_fasta_file_hash
read_fastq_file_hash
sequences_hash_to_array
sequences_array_to_hash
create_fasta_file
create_fastq_file
shuffle_fasta_file
shuffle_fastq_file
shuffle_fastq_files
split_fasta_file
split_fastq_file
split_file
convert_phred_ascii_to_numeric
convert_phred_numeric_to_ascii
mean_phred_ascii_quality
mean_multiple_phred_ascii_qualities
count_lengths
shuffle_seqs
splice_seqs
extract_random_seqs_from_file
extract_random_seqs_from_fasta
extract_random_seqs_from_fastq
extract_header_data
count_seqs_from_file
count_seqs_from_fasta
count_seqs_from_fastq
dna_to_prot
translate_seqs
retrieve_reading_frames
is_stop_codon
convert_nt_from_iupac
convert_aa_from_iupac
convert_nt_to_iupac
convert_aa_to_iupac
convert_aa_3_to_1
convert_aa_1_to_3
clean_sequence
complementary_sequence
reverse_sequence
iupac_reverse_complementary
consensus_sequence
random_sequence
unambiguous_dna_sequences
number_nt_combs
regex
prosite_to_regex
regex_to_prosite
binary_score_nts
detect_sequence_errors
insert_sequence_errors
print_sequence_errors
compare_sequences
count_errors
detect_chimeras
is_chimera
align_seqs_from_file
map_reads_from_files
align_seqs
align_seqs_with_threads
multiple_align_seqs
align_2seqs
align_seqs2one
align_seqs2seqs
cluster_seqs
cluster_identical_seqs
cluster_identical_seqs_with_threads
cluster_illumina_seqs
cluster_umi_seqs
execute_blastp
execute_blastpshort
execute_blastn
execute_blastnshort
execute_megablast
execute_gassst
execute_bwa
execute_bowtie
execute_bowtie2
execute_samtools
execute_bcftools
execute_regex_match
execute_cdhit_est
execute_cdhit_454
execute_mafft
execute_flash
extract_align_bam_data
extract_blast_data
extract_align_blast_data
extract_align_gassst_data
extract_align_regex_data
extract_cdhit_data
extract_vcf_data
extract_coverage_data
obtain_cols_from_alignment
alignment_to_single_line
store_data_dump
recover_data_dump
generate_md5
is_folder_empty
list_dir
convert_to_alphanum
group_taxo
read_taxo_data
print_taxo_data
matrices_to_transfac
transfac_to_matrices
normalize_transfac
matrix_to_sequences
execute_weblogo


);


######################################################################################################

# Routes to binary files and databases
my $TOOLSDIR = dirname (__FILE__).'/tools/';
my $BLASTPEXE = $TOOLSDIR.'blastp -task blastp'; # -num_threads 4";
my $BLASTPSHORTEXE = $TOOLSDIR.'blastp -task blastp-short'; # -num_threads 4";
my $BLASTNEXE = $TOOLSDIR.'blastn -task blastn'; # -num_threads 4";
my $BLASTNSHORTEXE = $TOOLSDIR.'blastn -task blastn-short'; # -num_threads 4";
my $MEGABLASTEXE = $TOOLSDIR.'blastn -task megablast'; # -num_threads 4";
my $MAKEBLASTDBEXE = $TOOLSDIR.'makeblastdb';
my $NCBI_NR_DATABASE = '/data/ncbi/nr'; # non-redundant protein
my $NCBI_NT_DATABASE = '/data/ncbi/nt'; # non-redundant nucleotide
my $GASSSTEXE = $TOOLSDIR.'gassst -m 0'; # Blast like format (human readable)
my $BWAMEMEXE = $TOOLSDIR.'bwa mem';
my $BWABUILDEXE = $TOOLSDIR.'bwa index';
my $BOWTIEEXE = $TOOLSDIR.'bowtie';
my $BOWTIEBUILDEXE = $TOOLSDIR.'bowtie-build';
my $BOWTIE2EXE = $TOOLSDIR.'bowtie2';
my $BOWTIE2BUILDEXE = $TOOLSDIR.'bowtie2-build';
my $SAMTOOLS = $TOOLSDIR.'samtools';
my $BCFTOOLS = $TOOLSDIR.'bcftools';
my $CDHITEXE = $TOOLSDIR.'cd-hit';
my $CDHITESTEXE = $TOOLSDIR.'cd-hit-est';
my $CDHIT454EXE = $TOOLSDIR.'cd-hit-454';
my $CLUSTALW2EXE = $TOOLSDIR.'clustalw2';
my $MAFFTEXE = $TOOLSDIR.'mafft';
my $NEEDLEMANWUNSCHEXE = $TOOLSDIR.'needleman_wunsch';
my $SMITHWATERMANEXE = $TOOLSDIR.'smith_waterman';
my $NEEDLEEXE = 'needle'; # It's necessary to install EMBOSS package because complex dependencies
my $NEEDLEALLEXE = 'needleall'; # It's necessary to install EMBOSS package because complex dependencies
my $FOGSAAEXE = $TOOLSDIR.'fogsaa'; # Ultrafast pairwise global alignment
my $FLASHEXE = $TOOLSDIR.'flash2'; # Paired-end read merging program
my $WEBLOGOEXE = $TOOLSDIR.'weblogo/weblogo'; # Generates sequence logos
# my $BLASTNEXE = dirname(__FILE__).'/../../bin/blastn -task blastn'; # -num_threads 4";
# my $BLASTNSHORTEXE = dirname(__FILE__).'/../../bin/blastn -task blastn-short'; # -num_threads 4";
# my $MEGABLASTEXE = dirname(__FILE__).'/../../bin/blastn -task megablast'; # -num_threads 4";
# my $MAKEBLASTDBEXE = dirname(__FILE__).'/../../bin/makeblastdb';
# my $NCBI_NT_DATABASE = '/data/nt/nt';
# my $GASSSTEXE = dirname(__FILE__).'/../../bin/gassst -m 0'; # Blast like format (human readable)
# my $CDHITEXE = dirname(__FILE__).'/../../bin/cd-hit';
# my $CDHITESTEXE = dirname(__FILE__).'/../../bin/cd-hit-est';
# my $CDHIT454EXE = dirname(__FILE__).'/../../bin/cd-hit-454';
# my $CLUSTALW2EXE = dirname(__FILE__).'/../../bin/clustalw2';
# my $MAFFTEXE = dirname(__FILE__).'/../../bin/mafft';

#################################################################################

# Trims spaces at start and end of a string
sub trim {
	return $_[0] =~ s/^\s+|\s+$//rg;
}

#################################################################################

# Calculates the mean of array values
sub mean {

	if (@_) {
		return sum(@_) / scalar @_;
	} else {
		return undef;
	}

}

#################################################################################

# Calculates the median of array values
sub median {

	my $count = scalar @_;
	
	@_ = sort {$a<=>$b} @_;

	if ($count == 1){
		return $_[0];
	## Odd or even
	} elsif ($count % 2 == 0){   
		return ( @_[$count/2] + @_[($count-2)/2] ) / 2;
	} else {
		return @_[($count-1)/2];
	}

}

#################################################################################

# Calculates the maximum number in an array
sub max {

	if (@_) {
		my $max;
		for (@_) {
			$max = $_ if !defined($max) || $_ > $max
		};
		return $max;
	} else {
		return undef;
	}

}

#################################################################################

# Calculates the minimum number in an array
sub min {

	if (@_) {
		my $min;
		for (@_) {
			$min = $_ if !defined($min) || $_ < $min;
		};
		return $min;
	} else {
		return undef;
	}

}

#################################################################################

# Calculates summatory of array values
sub sum {

	my $sum;

	if (@_) {
		map $sum += $_, @_;
		return $sum;
	} else {
		return undef;
	}

}

#########################################################################################

# Returns TRUE if an expression is numerical
sub is_numeric {

	my ($exp) = @_;

# 	if ($exp =~ /^(-?[\d.\-]*e[\d.\-]+|-?[\d.\-]\^[\d.\-]+|-?[\d.\-]+)$/){
	if (looks_like_number($exp)){
		return 1;
	} else {
		return 0;
	}

}

#########################################################################################

# Returns the percentage value if an expression is a percentage
sub is_percentage {

	my ($exp) = @_;

	if ($exp =~ /([\d\.]+)%/) {
		$exp = $1;
	}
	if (!is_numeric($exp) || $exp<0 || $exp>100){
		return -1;
	} else {
		return $exp;
	}
}

#################################################################################

# Calculates the decimal logarithm of a number
sub log10 {
	my $n = shift;
	return log($n)/log(10);
}

#################################################################################

# Calculates the variance of the array values (assuming they are a sample, N-1)
sub variance {


	my $mean = mean(@_);

	return sum( map {($_-$mean)**2} @_)/(scalar @_ - 1);

}

#################################################################################

# Calculates the standard deviation of the array values (assuming they are a sample, N-1)
sub stdeviation {

	return sqrt(variance(@_));

}

#################################################################################

# Calculates correlation between 2 vectors
sub correlation {

	my ($x,$y) = @_;

	my ($num,$denom)=(0,0);
	my ($diff_x, $diff_y);
	my ($sum_diff_x_sq, $sum_diff_y_sq)=(0,0);
	my ($mean_x,$mean_y)=(0,0);
	my $corr=0;

	if ($#{$x} != $#{$y}){
		die "\n# Cannot calculate correlation of two different length vectors\n";
	}

	for (my $i=0; $i<=$#{$x}; $i++){
		$mean_x += $x->[$i];
		$mean_y += $y->[$i];
	}

	$mean_x = $mean_x/(scalar @{$x});
	$mean_y = $mean_y/(scalar @{$y});

	if ($mean_x != 0 && $mean_y != 0){
		for (my $i=0; $i<=$#{$x}; $i++){
			$diff_x = $x->[$i] - $mean_x;
			$diff_y = $y->[$i] - $mean_y;
			$num += $diff_x * $diff_y;
			$sum_diff_x_sq += ($diff_x*$diff_x);
			$sum_diff_y_sq += ($diff_y*$diff_y);
		}
		if ($num != 0){
			$corr += $num / (sqrt($sum_diff_x_sq * $sum_diff_y_sq));
		} else {
			$corr += 0;
		}
	}

	return $corr;
}

#################################################################################

# Calculates covariance between 2 vectors
sub covariance {

	my ($x,$y) = @_;

	my ($num,$denom)=(0,0);
	my ($diff_x, $diff_y);
	my ($mean_x,$mean_y)=(0,0);
	my $var=0;

	if ($#{$x} != $#{$y}){ die "\n# Cannot calculate correlation of two different length vectors\n";}

	for (my $i=0; $i<=$#{$x}; $i++){
		$mean_x += $x->[$i];
		$mean_y += $y->[$i];
	}
	$mean_x = $mean_x/(scalar @{$x});
	$mean_y = $mean_y/(scalar @{$y});

	if ($mean_x != 0 && $mean_y != 0){
		for (my $i=0; $i<=$#{$x}; $i++){
			$diff_x = $x->[$i] - $mean_x;
			$diff_y = $y->[$i] - $mean_y;
			$num += $diff_x * $diff_y;
		}
		$var = $num / (scalar @{$x} - 1);
	}

	return $var;
}


#################################################################################

# Calculates percentiles
sub percentile {

	my ($values,$p) = @_;

	# Sort
	@$values = sort {$a <=> $b} @$values; 

	# Returns $p percentile
	return $values->[sprintf("%.0f",(($p/100)*($#{$values})))];

}

#################################################################################

# Calculate the bias-corrected Chao1 richness estimator
sub chao1 {

	my %obs = @_;
	
	my ($s,$f1,$f2) = (0,0,0);
	while (my ($event, $times) = each %obs) {
		if (defined($times) && $times) {
			$s++;
			if ($times > 1){
				$f2++;
			} else {
				$f1++;
			}
		}
	}

	return sprintf("%.0f", $s+(($f1*($f1-1))/(2*($f2+1))) );

}

#################################################################################

# Calculate the variations with repetition of a number of letters
sub variations {
	my ($letters,$num) = @_;
	my $last = [map $letters->[0] , 0 .. $num-1];
	my $result;
	while (join('',@{$last}) ne $letters->[$#{$letters}] x $num) {
		push(@{$result},[@{$last}]);
		$last = char_add($letters,$last,$num-1);
		print '';
	}
	push(@{$result}, $last);
	return $result;
}

sub char_add{
	my ($digits,$string,$char) = @_;
	if ($string->[$char] ne $digits->[$#{$digits}]){
		my ($match) = grep { $digits->[$_] eq $string->[$char]} 0 .. $#{$digits};
		$string->[$char] = $digits->[$match+1];
		return $string;
	} else {
		$string = changeall($string,$digits->[0],$char);
		return char_add($digits,$string,$char-1);
	}
}

sub changeall {
	my ($string,$char,$start,$end) = @_;
	if (!defined($start)){$start=0;}
	if (!defined($end)){$end=0;}
	if ($end == 0) {$end = $#{$string};}
	for(my $i=$start; $i<=$end; $i++){
		$string->[$i] = $char;
	}
	return $string;
}

#################################################################################

# Calculates probabilities conforming to the binomial distribution
# Returns the probability of an event occurring $k times
# in $n attempts, where the probability of it occurring
# in a single attempt is $p
sub binomial_probability {

	my ($n, $k, $p) = @_;
	
	return $k = 0 if $p ==0;
	return $k != $n if $p ==1;
	return binomial_coef($n, $k) * $p**$k * (1-$p)**($n-$k);

}

# Number of ways (combinations) to choose $k elements from a set of $n elements,
# when the order of selection is irrelevant
sub binomial_coef {

	my ($n, $k) = @_;
	my ($result, $j) = (1, 1);
	
	return 0 if $k > $n || $k < 0;
	$k = ($n - $k) if ($n - $k) < $k;
	
	while ( $j <= $k) {
		$result *= $n--;
		$result /= $j++;
	}
	
	return $result;

}


#################################################################################

# Performs WPGMA (Weighted Pair Group Method with Arithmetic Mean) clustering from a matrix of distances
sub wpgma {

	my ($distances_) = @_;
	
	return upgma($distances_, 'wpgma');
	
}

#################################################################################

# Performs UPGMA (Unweighted Pair Group Method with Arithmetic Mean) clustering from a matrix of distances
sub upgma {

	my ($distances_,$type) = @_;

	my $tree;

	if (!defined($type)){
		$type = 'wpgma';
	}
# 	$distances_ = { 
# 		'a' => { 'b' => 17, 'c' => 21, 'd' => 31, 'e' => 23 },
# 		'b' => { 'c' => 30, 'd' => 34, 'e' => 21 },
# 		'c' => { 'd' => 28, 'e' => 39 },
# 		'd' => { 'e' => 43 },
# 		};
# 	$distances_ = { 
# 		'a' => { 'b' => 20, 'c' => 60, 'd' => 100, 'e' => 90 },
# 		'b' => { 'c' => 50, 'd' => 90, 'e' => 80 },
# 		'c' => { 'd' => 40, 'e' => 50 },
# 		'd' => { 'e' => 30 },
# 		};
# 	$distances_ = { 
# 		'turtle' => { 'man' => 19, 'tuna' => 27, 'chicken' => 8, 'moth' => 33, 'monkey' => 18, 'dog' => 13 },
# 		'man' => { 'tuna' => 31, 'chicken' => 18, 'moth' => 36, 'monkey' => 1, 'dog' => 13 },
# 		'tuna' => { 'chicken' => 26, 'moth' => 41, 'monkey' => 32, 'dog' => 29 },
# 		'chicken' => { 'moth' => 31, 'monkey' => 17, 'dog' => 14 },
# 		'moth' => { 'monkey' => 35, 'dog' => 28 },
# 		'monkey' => { 'dog' => 12 },
# 		};

	my $distances = { %$distances_ };
	
	# Number of clustering steps
	my $step = scalar keys %$distances;

	my %node_size;
	foreach my $e1 (keys %$distances){
		$node_size{$e1} = 1;
		foreach my $e2 (keys %{$distances->{$e1}}){
			$node_size{$e2} = 1;
		}
	}
	my %tree_distances;
	for (my $i=1; $i<=$step; $i++) {
		$tree_distances{$step} = 0;
	}

	while ($step) {
		my @min_dist;
		my $min_dist;
		# Finds the minimum distance between members
		my %pairs;
		foreach my $e1 (keys %$distances){
			foreach my $e2 (keys %{$distances->{$e1}}){
				# Removes duplicates and distances to the same element
				if ($e1 eq $e2 || defined($pairs{"$e2-$e1"})){
					delete($distances->{$e1}{$e2});
					next;
				} elsif (!defined($distances->{$e1}{$e2})) {
					next;
				}
				if (!defined($min_dist) || $distances->{$e1}{$e2} < $min_dist){
					$min_dist = $distances->{$e1}{$e2};
					$min_dist[0] = $e1;
					$min_dist[1] = $e2;
				}
				$pairs{"$e1-$e2"} = 1;
			}
		}
		if (!defined($tree_distances{$min_dist[0]})){
			$tree_distances{$min_dist[0]} = 0;
		}
		if (!defined($tree_distances{$min_dist[1]})){
			$tree_distances{$min_dist[1]} = 0;
		}
		$tree ->{$step}{$min_dist[0]} = $min_dist/2 - $tree_distances{$min_dist[0]};
		$tree ->{$step}{$min_dist[1]} = $min_dist/2 - $tree_distances{$min_dist[1]};
		$tree_distances{$step} = max($tree ->{$step}{$min_dist[0]},$tree ->{$step}{$min_dist[1]});
		foreach my $e1 (keys %$distances){
			foreach my $e2 (keys %{$distances->{$e1}}){
				if ($e1 eq $min_dist[0] && $e2 eq $min_dist[1]
				|| $e2 eq $min_dist[0] && $e1 eq $min_dist[1]){
					delete($distances->{$e1}{$e2});
					delete($distances->{$e2}{$e1});
					next;
				}
				my ($dist1,$dist2,$size1,$size2);
				my @elems = ($e1, $e2);
				for (my $i=0; $i<2; $i++){
					$e1 = $elems[$i];
					$e2 = $elems[$i-1];
					if (($e1 eq $min_dist[0] || $e1 eq $min_dist[1])){
						if (defined($distances->{$min_dist[0]}{$e2})){
							$dist1 = $distances->{$min_dist[0]}{$e2};
							delete($distances->{$min_dist[0]}{$e2});
						} else {
							$dist1 = $distances->{$e2}{$min_dist[0]};
							delete($distances->{$e2}{$min_dist[0]});
						}
						if (defined($distances->{$min_dist[1]}{$e2})){
							$dist2 = $distances->{$min_dist[1]}{$e2};
							delete($distances->{$min_dist[1]}{$e2});
						} else {
							$dist2 = $distances->{$e2}{$min_dist[1]};
							delete($distances->{$e2}{$min_dist[1]});
						}
						if (defined($dist1) && defined($dist2) && !defined($distances->{$step}{$e2})){
							if ($type eq 'upgma'){
								$size1 = $node_size{$min_dist[0]}+$node_size{$e2}-1;
								$size2 = $node_size{$min_dist[1]}+$node_size{$e2}-1;
							} elsif ($type eq 'wpgma'){
								$size1 = 1;
								$size2 = 1;
							}
							$distances->{$step}{$e2} = ($dist1*$size1+$dist2*$size2)/($size1+$size2);
						}
					}
				}
			}
		}
		$node_size{$step} += $node_size{$min_dist[0]};
		$node_size{$step} += $node_size{$min_dist[1]};
		delete($distances->{$min_dist[0]});
		delete($distances->{$min_dist[1]});
		$step--;
	}
	
	return $tree;

}

#################################################################################

# Searches a string inside an array
# Example: in_array($sentence, "/$word/", 1);
sub in_array {

	my ($array,$search_for,$search_posic) = @_;

	my @posic;

	# If search a pattern
	if ($search_for =~ /^\/(.+)\/$/){
		for (my $i=0; $i<=$#{$array}; $i++) {
			if (defined($array->[$i]) && $array->[$i] =~ /$1/i){
				push(@posic, $i);
			}
		}
	} else {
		for (my $i=0; $i<=$#{$array}; $i++) {
			if (defined($array->[$i]) && $array->[$i] eq $search_for){
				push(@posic, $i);
			}
		}
	}

	if (defined($search_posic) && $search_posic){
		(scalar @posic)?return (1,\@posic):return (0);
	} else {
		(scalar @posic)?return 1:return 0;
	}
}


#################################################################################

# Searches the most frequent value into an array
sub most_frequent {

	$_[$_{$_}] = $_ for map{$_{$_}++; $_} @_;

	return $_[-1];

}

################################################################################

# Converts an array into a hash
sub array_to_hash {

	return map { $_ => 1 } @_;

}


#################################################################################

# Extracts a pattern from an array
sub extract_pattern_from_array {

	my ($array,$search_for) = @_;
	
	my $new_array;
	
	if ($search_for =~ /^\/(.+)\/$/){
		$search_for = $1;
	}
	foreach my $row (@$array) {
		if ($row =~ /($1)/i){
			push(@{$new_array}, $1);
		}
	}

	return $new_array;

}

#################################################################################

# Compares two arrays, returning 1 if they are identical and 0 and the first substitution position if not
sub compare_arrays {

	my ($array1,$array2) = @_;

	for (my $i=0; $i<=$#{$array1}; $i++) {
		if (!defined($array1->[$i]) || !defined($array2->[$i])){
			return (0,$i+1);
		}
		if ($array1->[$i] ne $array2->[$i]){
			return (0,$i+1);
		}
	}
	if ($#{$array1} == $#{$array2}) {
		return (1, undef);
	} else {
		return (0,$#{$array1}+2);
	}

}

#################################################################################

# Compares two strings, returning 1 if they are identical and 0 and the first substitution position if not
sub compare_strings {

	my ($string1,$string2) = @_;

	return compare_arrays([split('',$string1)],[split('',$string2)]);

}


#################################################################################

# Returns an array of unique items in the arguments list.
# Copyright (c) 2007 Sergei A. Fedorov.
sub unique {
	my @array = @_;
	my (%hash, @unique);
	foreach my $a (@array) {
		if (!defined($hash{$a})){
			push(@unique,$a);
			$hash{$a} = 1;
		}
	}
	return @unique ;
}

#################################################################################

# Returns an intersection of two arrays passed as arguments, keeping the order of the second parameter.
# Copyright (c) 2007 Sergei A. Fedorov.sub unique(@) {
sub intersect {
	my ($array1,$array2) = @_;
	my %e = map { $_ => undef } @{$array2};
	return grep { exists( $e{$_} ) } @{$array1};
}

#################################################################################

# Reads content from a file
sub read_from_file {
	my ($file) = @_;
	
	my @content;
	
	if (is_gzip($file)){
# 		# open FILE, "<:gzip", $file or die "# cannot read $file\n"; # Requires PerlIO::gzip that gives problems to copy the libraries (*.pm)
		open(FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_zip($file)){
		# open(FILE, "unzip -p '$file' |") or die "# cannot read $file\n";
		open(FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_bzip2($file)){
		# open(FILE, "unzip -p '$file' |") or die "# cannot read $file\n";
		open(FILE, "bzcat '$file' |") or die "# cannot read $file\n";
	} else {
		open(FILE,$file) || die "# cannot read $file\n";
	}

	@content = <FILE>;
	close FILE;

	return \@content;
}

#################################################################################

# Writes content to a file
sub write_to_file {

	my ($outfile,$content,$mode)=@_;
	if (!defined($content)){
		 die "\nERROR 'write_to_file': cannot create empty '$outfile'\n\n";
	}
	if (!defined($mode)){
		$mode=1;
	}
	if ($mode eq 1){
		open(FILE,">$outfile") || die "# cannot create $outfile\n";
		print FILE $content;
		close FILE;
	} elsif ($mode eq 2){
		open(FILE,">>$outfile") || die "# cannot create $outfile\n";
		print FILE $content;
		close FILE;
	} elsif ($mode eq 'zip'){
		$outfile .= ".zip";
		my $FILE = new IO::Compress::Zip("$outfile") || die "\nERROR 'write_to_file': cannot create '$outfile'\n\n";
		$FILE->print($content);
		$FILE->close();
	} elsif ($mode eq 'gzip'){
		$outfile .= ".gzip";
		my $FILE = new IO::Compress::Gzip("$outfile") || die "\nERROR 'write_to_file': cannot create '$outfile'\n\n";
		$FILE->print($content);
		$FILE->close();
	}
	return $outfile;
}

#################################################################################

# Downloads a file from a URL
sub download_file {

	my ($url,$outfile) = @_;

	if (!defined($outfile)){
		$outfile="/tmp/".random_file_name();
	}

	#print "wget -q -O $outfile $url\n";exit;
	`wget -q -O $outfile $url`;

	if (-e $outfile){
		return $outfile;
	} else {
		return undef;
	}
}

#################################################################################

# Checks if is a COMPRESSED file
sub is_compressed {

	my ($file) = @_;

	if (is_zip($file)){
		return 'zip';
	} elsif (is_gzip($file)){
		return 'gzip';
	} elsif (is_bzip2($file)){
		return 'bzip2';
	} elsif (is_tar($file)){
		return 'tar';
	} elsif (is_tgz($file)){
		return 'tgz';
	} elsif (is_tbz2($file)){
		return 'tbz2';
	} else {
		return 0;
	}

}

#################################################################################

# Checks if is GZIP format file
sub is_gzip {

	my ($infile) = @_;

# if (!defined($infile)){
# print '';
# }
	my $output = `file -L "$infile"`;

	if ($output =~ /gzip compressed data/g) {
		return  1;
	} else {
		return 0;
	}
}

#################################################################################

# Checks if is ZIP format file
sub is_zip {

	my ($infile) = @_;

	my $output = `file -L "$infile"`;

	if ($output =~ /Zip archive data/gi) {
		return  1;
	} else {
		return 0;
	}
}

#################################################################################

# Checks if is BZIP2 format file
sub is_bzip2 {

	my ($infile) = @_;

	my $output = `file -L "$infile"`;

	if ($output =~ /bzip2 compressed data/g) {
		return  1;
	} else {
		return 0;
	}
}

#################################################################################

# Checks if is TAR format file
sub is_tar {

	my ($infile) = @_;

	my $output = `file -L "$infile"`;

	if ($output =~ /POSIX tar archive/gi) {
		return  1;
	} else {
		return 0;
	}

}

#################################################################################

# Checks if is TGZ/TAR.GZ format file
sub is_tgz {

	my ($infile) = @_;

	if (is_gzip($infile)){
		my @list = list_gzip($infile);
		if ($#list == 0 && $list[0] =~ /\.tar$/) {
			return 1;
		}
	}
	
	return 0;

}

#################################################################################

# Checks if is TBZ2/TAR.BZ2 format file
sub is_tbz2 {

	my ($infile) = @_;

	if (is_bzip2($infile)){
		my @list = list_tbz2($infile);
		if ($#list == 0 && $list[0] =~ /\.tar$/) {
			return 1;
		}
	}
	
	return 0;

}

#################################################################################

# Retrieves an array with the names of the files/folders contained by a ZIP format file
sub list_zip {

	my ($infile) = @_;

	my @list;

	foreach (split("\n",`unzip -l "$infile"`)) {
		if (/\d+\:\d+\s+(.+)$/) {
			push(@list,trim($1));
		}
	}

	return @list;

}

#################################################################################

# Retrieves an array with the name of the file contained by a GZIP format file
# GZIP files only contain ONE file (TGZ files can contain several)
sub list_gzip {

	my ($infile) = @_;

	my @list;

	foreach (split("\n",`gunzip -l "$infile"`)) {
		if (/[\d\.]+\%\s+(.+)$/) {
			push(@list,trim($1));
		}
	}

	return @list;

}

#################################################################################

# Retrieves an array with the names of the files/folders contained by TAR or TGZ/TBZ2 format files
sub list_tar {

	my ($infile) = @_;

	my @list;

	foreach (split("\n",`tar tf "$infile" 2>&1`)) {
		if (!/^tar\:/) {
			push(@list,trim($_));
		}
	}

	return @list;

}

#################################################################################

# Retrieves an array with the names of the files/folders contained by TAR or TGZ/TBZ2 format files
sub list_tgz {

	return list_tar($_[0]);

}

#################################################################################

# Retrieves an array with the names of the files/folders contained by TAR or TGZ/TBZ2 format files
sub list_tbz2 {

	return list_tar($_[0]);

}

#################################################################################

# Checks if is Excel format file
sub is_xlsx {

	my ($infile) = @_;

	my $output = `file -L "$infile"`;

	if ($output =~ /Microsoft (OOXML|Excel)/gi) {
		return  1;
	} elsif ($infile =~ /\.xlsx$/ && $output =~ /Zip archive data/) {
		return 1;
	} else {
		return 0;
	}
}

#################################################################################

# Checks if is a compressed file with multiple files inside
sub is_multifile {

	my ($infile) = @_;

	if (is_tar($infile) && scalar list_tar($infile) > 1){
		return 1;
	} elsif (is_zip($infile) && scalar list_zip($infile) > 1){
		return 1;
	} elsif (is_tgz($infile) && scalar list_tgz($infile) > 1){
		return 1;
	} elsif (is_tbz2($infile) && scalar list_tbz2($infile) > 1){
		return 1;
	} else {
		return 0;
	}

}

#################################################################################

# Extracts paired end read files from one or two compressed files with multiple files inside (multifile)
sub extract_paired_read_files_from_multifiles {

	my ($multifiles,$outdir) = @_;
	
	# Uncompress the files in temporal folder
	if (!defined($outdir)){
		$outdir = "/tmp/".random_file_name();
	}

	# Checks if the input file is a set of packed files (.zip or .tar.gz/.tgz)
	my $multifile_lists = [];
	for (my $i=0; $i<=$#{$multifiles}; $i++) {
		if (is_multifile($multifiles->[$i])){
			$multifile_lists->[$i] = [ decompress($multifiles->[$i],undef,$outdir) ];
		} else {
			$multifile_lists->[$i] = [ $multifiles->[$i] ];
		}
	}
	
	return extract_paired_read_files($multifile_lists);

}

#################################################################################

# Extracts paired end read files from a path
sub extract_paired_read_files_from_path {

	my ($path) = @_;

	my $file_list = [];
	$file_list->[0] = [ read_files_from_path($path) ];

	return extract_paired_read_files($file_list);

}

#################################################################################

# Extracts paired end read files from one or two lists of files
sub extract_paired_read_files {

	my ($file_lists) = @_;

	my $paired_reads_files = [];

	my %previously_paired_tmp_files;
	foreach my $file1 (@{$file_lists->[0]}){

		if (defined($previously_paired_tmp_files{$file1})){ # Already paired file
			next;
		}

		my ($file1_name, $file1_prefix);
		my ($file2, $file2_name, $file2_prefix);
		# If only two files are specified, one for each paired-end reads
		if ($#{$file_lists->[0]}==0) {
			$file2 = $file_lists->[1][0];
		# If there are two multiple files, will search for the same file prefix in the files from the second mulifile
		} elsif (defined($file_lists->[1])) {
			if ($file1 =~ /(.+\/)?(.+?)(_?R1.+)/ || $file1 =~ /(.+\/)?(.+?)(_1\..+)/){
				$file1_prefix = $2;
				# $file1_name = $2.$3;
			} elsif ($file1 =~ /(.+\/)?(.+?)$/) {
				$file1_prefix = $file1_name = $2;
			} else {
				$file1_prefix = $file1_name = $file1;
			}
			foreach my $file2_ (@{$file_lists->[1]}){
				if (defined($previously_paired_tmp_files{$file2_})){ # Already paired file
					next;
				} elsif ($file2_ =~ /(.+\/)?(.+?)(_?R2.+)/ || $file2_ =~ /(.+\/)?(.+?)(_2\..+)/){
					$file2_prefix = $2;
					$file2_name = $2.$3;
				} elsif ($file2_ =~ /(.+\/)?(.+?)$/){
					$file2_prefix = $file2_name = $2;
				}
				if (defined($file2_prefix) && $file1_prefix eq $file2_prefix){
					$file2 = $file2_;
					last;
				}
			}
# 			if ($file1 =~ /(.+\/)?(.+?)$/){
# 				$file1_prefix = $file1_name = $2;
# 			} else {
# 				$file1_prefix = $file1_name = $file1;
# 			}
# 			foreach my $file2_ (@{$file_lists->[1]}){
# 				if (defined($previously_paired_tmp_files{$file2_})){ # Already paired file
# 					next;
# 				} elsif ($file2_ =~ /(.+\/)?(.+?)$/){
# 					$file2_prefix = $file2_name = $2;
# 				} else {
# 					$file2_prefix = $file2_name = $file2_;
# 				}
# 				if ($file1_prefix eq $file2_prefix){
# 					$file2 = $file2_;
# 					last;
# 				}
# 			}
		# If there is only one multiple file, will search among the other files of the single multifile
		} else{
			if ($file1 =~ /(.+\/)?(.+?)(_?R1.+)/ || $file1 =~ /(.+\/)?(.+?)(_1\..+)/){
				$file1_prefix = $2;
				# $file1_name = $2.$3;
				foreach my $file2_ (@{$file_lists->[0]}){
					if (defined($previously_paired_tmp_files{$file2_})){ # Already paired file
						next;
					} elsif ($file2_ =~ /(.+\/)?(.+?)(_?R2.+)/ || $file2_ =~ /(.+\/)?(.+?)(_2\..+)/){
						$file2_prefix = $2;
						$file2_name = $2.$3;
					}
					if (defined($file2_prefix) && $file1_prefix eq $file2_prefix){
						$file2 = $file2_;
						last;
					}
				}
			}
		}
		# Checks each paired-read file (R1 and R2 usually)
		if (!is_fastq($file1)){
			printf("\tWARNING: '%s' file has not FASTQ format.\n", $file1);
			next;
		}
		if (defined($file2) && !is_fastq($file2)){
			printf("\tWARNING: '%s' file has not FASTQ format.\n", $file1);
			next;
		}
		if (defined($file2) && count_seqs_from_fastq($file1) != count_seqs_from_fastq($file2)){
			printf("\tWARNING: '%s' and '%s' have different number of reads.\n", $file1, $file2);
			next;
		}
		if (defined($file2)){
			$previously_paired_tmp_files{$file1} = 1;
			$previously_paired_tmp_files{$file2} = 1;
			push(@$paired_reads_files, [$file1, $file2]);
		}
	}
	foreach my $file (@{$file_lists->[0]}){

		if (defined($previously_paired_tmp_files{$file})){ # Already paired file
			next;
		}
		
		printf("\tWARNING: Couldn't find the '%s' paired-end file.\n", $file);

		# Checks each paired-read file (R1 and R2 usually)
		if (!is_fastq($file)){
			printf("\tWARNING: '%s' file has not FASTQ format.\n", $file);
			next;
		}

		push(@$paired_reads_files, [$file]);

	}
	
	return $paired_reads_files;

}


#################################################################################

# Reads files from a path
sub read_files_from_path {

	my ($path) = @_;

	my @file_list;
	opendir(DIR, $path) or die "ERROR: Couldn't open $path: $!\n";
	while (readdir(DIR)) {
		next if /^\.+$/i;
		push(@file_list, $path."/".$_);
	}
	closedir DIR;

	return @file_list;

}

#################################################################################

# Compresses an input file
sub compress {
	my ($infile,$type,$outfile)= @_;
	
	if (!defined($type) || $type !~ /zip/) {
		$type = 'zip';
	}
	if (!defined($outfile)){
		$outfile = "/tmp/".random_file_name();
	}

	my $output;

	# If a single file is given in input
	if (ref($infile) ne "ARRAY"){
		if ($type eq 'gzip'){
			$output = `gzip -c "$infile" > "$outfile"`;
		} elsif ($type eq 'bzip2'){
			$output = `bzip2 -c "$infile" > "$outfile"`;
		} elsif ($type eq 'zip'){
			$output = `zip -qm "$outfile" "$infile"`;
		}
	# If not, creates a multifile
	} else {
		my $current_path = cwd();
		my $change_path;
		if ($infile->[0] =~ /(\/tmp\/.\d+|\/tmp.+)\//) {
			$change_path = $1;
			foreach my $infile_ (@$infile){
				if ($infile_ =~ /$change_path\/(.+)/) {
					$infile_ = $1;
				}
			}
		}
		chdir($change_path);
		my $infiles = '"'.join('" "',@$infile).'"';
		if ($type eq 'tar'){
			$output = `tar -cf "$current_path/$outfile" -C $change_path $infiles`;
		} elsif ($type eq 'gzip' || $type eq 'tgz'){
			# print "tar -czf \"$outfile\" -C $change_path $infiles\n";exit;
			$output = `tar -czf "$current_path/$outfile" -C $change_path $infiles`;
		} elsif ($type eq 'bzip2' || $type eq 'tbz2'){
			$output = `tar -cjf "$current_path/$outfile" -C $change_path $infiles`;
		} elsif ($type eq 'zip'){
			$output = `zip -qm "$current_path/$outfile" $infiles`;
		}
		chdir($current_path);
	}

	if (!defined($output)) {
		return undef;
	} else {
		return $outfile;
	}

}

#################################################################################

# Decompresses an input file
sub decompress {
	my ($infile,$type,$outpath)= @_;

	if (!defined($type) || $type !~ /^(zip|gzip|bzip2|tar|tgz|tbz2)$/) {
		if (is_zip($infile)){
			$type = 'zip';
		} elsif (is_gzip($infile)){
			$type = 'gzip';
		} elsif (is_bzip2($infile)){
			$type = 'bzip2';
		} elsif (is_tar($infile)){
			$type = 'tar';
		} elsif (is_tgz($infile)){
			$type = 'tgz';
		} elsif (is_tbz2($infile)){
			$type = 'tbz2';
		} else {
			print "\nERROR 'decompress': Unknown file type.\n\n";
			return undef;
		}
	}
	if (!defined($outpath)){
		$outpath = "/tmp/".random_file_name();
	}
	my $is_multifile = is_multifile($infile);
	
	if ($type eq 'gzip'){
		if (!$is_multifile){
			`gunzip -c "$infile" > "$outpath"`;
			return $outpath;
		} elsif (is_tgz($infile)) {
			my @output;
			mkdir($outpath);
			foreach (split("\n",`tar -xzvf "$infile" -C "$outpath"`)) {
				my $outfile = $outpath."/".trim($_);
				if (!/^tar\:/ && -e $outfile && !-d $outfile) {
					push(@output,$outfile);
				}
			}
			return @output;
		}
	}elsif ($type eq 'bzip2'){
		if (!$is_multifile){
			`bunzip2 -c "$infile" > "$outpath"`;
			return $outpath;
		} elsif (is_tgz($infile)) {
			my @output;
			mkdir($outpath);
			foreach (split("\n",`tar -xjvf "$infile" -C "$outpath"`)) {
				my $outfile = $outpath."/".trim($_);
				if (!/^tar\:/ && -e $outfile && !-d $outfile) {
					push(@output,$outfile);
				}
			}
			return @output;
		}
	} elsif ($type eq 'zip'){
		if (!$is_multifile){
			`unzip -p "$infile" > "$outpath"`;
			return $outpath;
		} else {
			my @output;
			mkdir($outpath);
			foreach (split("\n",`unzip "$infile" -d "$outpath"`)) {
				if (/(creating|inflating|extracting)\:\s+(.+)/) {
					push(@output,trim($2));
				}
			}
			return @output;
		}
	} elsif ($type eq 'tar'){
		my @output;
		mkdir($outpath);
		foreach (split("\n",`tar -xvf "$infile" -C "$outpath"`)) {
			my $outfile = $outpath."/".trim($_);
			if (!/^tar\:/ && -e $outfile && !-d $outfile) {
				push(@output,$outfile);
			}
		}
		return @output;
	} elsif ($type eq 'tgz'){
		my @output;
		mkdir($outpath);
		foreach (split("\n",`tar -zxvf "$infile" -C "$outpath"`)) {
			my $outfile = $outpath."/".trim($_);
			if (!/^tar\:/ && -e $outfile && !-d $outfile) {
				push(@output,$outfile);
			}
		}
		return @output;
	} elsif ($type eq 'tbz2'){
		my @output;
		mkdir($outpath);
		foreach (split("\n",`tar -jxvf "$infile" -C "$outpath"`)) {
			my $outfile = $outpath."/".trim($_);
			if (!/^tar\:/ && -e $outfile && !-d $outfile) {
				push(@output,$outfile);
			}
		}
		return @output;
	}
	
	return undef;
}

#################################################################################

# Creates a random file name
sub random_file_name {
	my ($length) = @_;

	if (!defined($length)){
		$length = 15;
	}

	my $max_number = "9"x$length;

	return sprintf("%0".$length."u",int(rand($max_number)));
}

#################################################################################

# Checks if is FASTA format file
sub is_fasta {
	my ($infile) = @_;
	
	my @lines4;
	if (is_gzip($infile) || is_zip($infile)) {
		@lines4 = split("\n", `zcat "$infile" | head -4`);
	} elsif (is_bzip2($infile)) {
		@lines4 = split("\n", `bzcat "$infile" | head -4`);
	} else {
		@lines4 = split("\n", `cat "$infile" | head -4`);
	}
	my $read_seq = 0;
	foreach my $line (@lines4){
		if ($line =~ /^>/){
			$read_seq = 1;
		} elsif ($read_seq) {
			if (is_dna($line) || is_rna($line) || is_prot($line)) {
				return 1;
			} else {
				return 0;
			}
		}
	}
	
	return 0;
}

#################################################################################

# Checks if is FASTQ format file
sub is_fastq {

	my ($infile) = @_;
	
	my @lines8;
	if (is_gzip($infile) || is_zip($infile)) {
		@lines8 = split("\n", `zcat "$infile" | head -8`);
	} elsif (is_bzip2($infile)) {
		@lines8 = split("\n", `bzcat "$infile" | head -8`);
	} else {
		@lines8 = split("\n", `cat "$infile" | head -8`);
	}
	if ( $#lines8 > 3 && $lines8[0] =~ /^@/ && is_dna($lines8[1]) && $lines8[2] =~ /^\+/ && $lines8[4] =~ /^@/ ){
		return 1;
	}
	if ( $#lines8 == 3 && $lines8[0] =~ /^@/ && is_dna($lines8[1]) && $lines8[2] =~ /^\+/ ){
		return 1;
	}
	
	return 0;
}

#################################################################################

# Checks if is SAM format file (Sequence Alignment/Map Format)
sub is_sam {

	my ($infile) = @_;
	
	my @lines2 = split("\n", `cat "$infile" | head -2`);

	my $read_seq = 0;
	foreach my $line (@lines2){
		if ($line =~ /^\@HD/){
			$read_seq = 1;
		} elsif ($read_seq) {
			if ($line =~ /^\@SQ/){
				return 1;
			} else {
				return 0;
			}
		}
	}
	
	return 0;
}

#################################################################################

# Checks if is BAM format file (Binary Sequence Alignment/Map)
sub is_bam {

	my ($infile) = @_;
	
	if (is_gzip($infile)){
		my $first_line = `zcat "$infile" | head -1`;
		if ($first_line =~ /^BAM/){
			return 1;
		} else {
			return 0;
		}
	}
	
	return 0;
}

#################################################################################

# Checks if is BCF format file (Binary Call Format)
sub is_bcf {

	my ($infile) = @_;
	
	if (is_gzip($infile)){
		my $first_line = `zcat "$infile" | head -1`;
		if ($first_line =~ /^BCF/){
			return 1;
		} else {
			return 0;
		}
	}
	
	return 0;
}

#################################################################################

# Checks if is VCF format file (Variant Call Format)
sub is_vcf {

	my ($infile) = @_;

	my $first_line = `cat "$infile" | head -1`;

	if ($first_line =~ /^##fileformat=VCF/){
		return 1;
	} else {
		return 0;
	}
}

#################################################################################

# Checks if is TXT format file
sub is_txt {

	my ($infile) = @_;

	my $output = `file -L "$infile"`;

	if ($output =~ /ASCII text/g) {
		return  1;
	} else {
		return 0;
	}
}

#################################################################################

# Checks if is a DNA sequence
sub is_dna {

	my ($seq) = @_;

	if ($seq =~ /^[ACGTMRWSYKVHDBXNI\-\s\?\[\]\(\)\.]+$/i){
		return 1;
	} else {
		return 0;
	}
}

#################################################################################

# Checks if is a RNA sequence
sub is_rna {

	my ($seq) = @_;

	if ($seq =~ /^[ACGUXN\-\s\?\[\]\(\)]+$/i){
		return 1;
	} else {
		return 0;
	}
}

#################################################################################

# Checks if is a protein sequence
sub is_prot {

	my ($seq) = @_;

	if ($seq =~ /^[ARNDCQEGHLIKMFPSTWYVOUJBZX\-\s\*\?\[\]\(\)\.]+$/i){
		return 1;
	} else {
		return 0;
	}
}


#################################################################################

# Recognizes the NGS technology used
sub ngs_tech {

	my ($data) = @_;

	foreach (split("\n",$data)){
		if (/^\w+:\d+:\w+:\d:\d{4}:\d{5}/) {
			return 'Illumina';
		} elsif (/^\w{5}:\d{5}:\d{5}/) {
			return 'IonTorrent';
		} elsif (/^\w{7}\d{2}\w{5}/) {
			return '454';
		}
	}
	
	return undef;

}

#################################################################################

# Creates 2 array references with sequences and headers from FASTA data
sub read_fasta {

	my ($fasta_data,$id_type) = @_;

	my (@sequences,@headers,$seq,$header);

	foreach (@$fasta_data){
		next if (/^#/); # allow comments
		next if (/^\s*$/); # empty lines
		s/\012\015?|\015\012?//; # trim line break
		if (/^>\s?(.+)/){
			if(defined($seq)){
				$seq =~ s/\d|\s//g;
				push(@sequences,$seq);
				push(@headers,$header);
				undef($header);
				undef($seq);
			}
			$header = $1;
			# If a specific ID type is required
			if (defined($id_type)){
				# /^\s*(gi\|\d+\||(gb|tpg|sp)\||jgi\|\w+\||tr\|)?\s*([^\t|^\||^;|^#|\:]+)/
				if ($id_type eq 'genbank' && $1 =~ /^\s*(gi\|\d+\||(gb|tpg|sp)\||jgi\|\w+\||tr\|)?\s*([^\t|^\||^;|^#|\:]+)/){
					$header = $3;
				} elsif ($id_type eq 'ensembl' && $1 =~ /(ENS.+?)\s/){
					$header = $1;
				} elsif ($id_type eq 'uniprot' && $1 =~ /sp\|([\w\-]+)\|/){
					$header = $1;
				} else {
					$header = (split(/[\|\s]/,$1))[0]
				}
			}
		} else {
			$seq .= $_;
		}
	}
	# take care of last FASTA sequence
	if(defined($seq)){
		$seq =~ s/\d|\s//g;
		push(@sequences,$seq);
		push(@headers,$header);
	}
	
	return (\@sequences,\@headers);
}

#################################################################################

# Create 3 array references with sequences, headers and qualities from FASTQ data
sub read_fastq {

	my ($fastq_data,$read_qualities) = @_;

	if (!defined($read_qualities)) {
		$read_qualities = 0;
	}

	my (@sequences, @headers, @qualities);

	my $field_count = 0;
	foreach my $row (@$fastq_data){
		next if ($field_count != 3 && $row =~ /^#/); # allow comments
		$row =~ s/\012\015?|\015\012?//; # trim line break
		if ($field_count) {
			$field_count++;
		}
		if ($field_count == 2){
			push(@sequences,$row);
		} elsif ($field_count == 4){ 
			$field_count = 0;
			if ($read_qualities) {
				push(@qualities,$row);
			}
		} elsif ($row =~ /^@/){
			push(@headers,substr($row,1));
			$field_count = 1;
		}
	}

	if ($read_qualities) {
		return (\@sequences,\@headers,\@qualities);
	} else {
		return (\@sequences,\@headers);
	}
}

#################################################################################

# Creates 2 array references with sequences and headers from a FASTA/FASTQ file
sub read_sequence_file {

	my ($file) = @_;
	
	my $format;

	if (is_fastq($file)){
		$format = 'fastq';
	} elsif (is_fasta($file)){
		$format = 'fasta';
	}

	if ($format eq 'fastq'){
		return read_fastq_file($file)
	} elsif ($format eq 'fasta'){
		return read_fasta_file($file);
	} else {
		print "\nERROR 'read_sequence_file': File format not recognized.\n\n";
		exit;
	}
}

#################################################################################

# Creates 2 array references with sequences and headers from a FASTA file
sub read_fasta_file {

	my($file,$id_type) = @_;
	
	my (@sequences,@headers);

	if (is_gzip($file)){
		open(FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_zip($file)){
		open(FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_bzip2($file)){
		open(FILE, "bzcat '$file' |") or die "# cannot read $file\n";
	} else {
		open(FILE,$file) || die "# cannot read $file\n";
	}

	my ($seq,$header);
	while (<FILE>){
		next if (/^#/); # allow comments
		next if (/^\s*$/); # empty lines
		s/\012\015?|\015\012?//; # trim line break
		if (/^>\s?(.+)/){
			if(defined($seq)){
				$seq =~ s/\d|\s//g;
				push(@sequences,$seq);
				push(@headers,$header);
				undef($header);
				undef($seq);
			}
			$header = $1;
			# If a specific ID type is required
			if (defined($id_type)){
				# /^\s*(gi\|\d+\||(gb|tpg|sp)\||jgi\|\w+\||tr\|)?\s*([^\t|^\||^;|^#|\:]+)/
				if ($id_type eq 'genbank' && $1 =~ /^\s*(gi\|\d+\||(gb|tpg|sp)\||jgi\|\w+\||tr\|)?\s*([^\t|^\||^;|^#|\:]+)/){
					$header = $3;
				} elsif ($id_type eq 'ensembl' && $1 =~ /(ENS.+?)\s/){
					$header = $1;
				} elsif ($id_type eq 'uniprot' && $1 =~ /sp\|([\w\-]+)\|/){
					$header = $1;
				} else {
					$header = (split(/[\|\s]/,$1))[0];
				}
			}
		} else {
			$seq .= $_;
		}
	}
	# take care of last FASTA sequence
	if(defined($seq)){
		$seq =~ s/\d|\s//g;
		push(@sequences,$seq);
		push(@headers,$header);
	}
	
	return (\@sequences,\@headers);

}


#################################################################################

# Creates 3 array references with sequences, headers and qualities from a FASTQ file
sub read_fastq_file {

	my($file,$read_qualities) = @_;

	if (!defined($read_qualities)) {
		$read_qualities = 0;
	}

	my (@sequences, @headers, @qualities);

	if (is_gzip($file)){
		open(FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_zip($file)){
		open(FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_bzip2($file)){
		open(FILE, "bzcat '$file' |") or die "# cannot read $file\n";
	} else {
		open(FILE,$file) || die "# cannot read $file\n";
	}

	my $field_count = 0;
	while (<FILE>){
		next if ($field_count != 3 && /^#/); # allow comments
		s/\012\015?|\015\012?//; # trim line break
		if ($field_count) {
			$field_count++;
		}
		if ($field_count == 2){
			push(@sequences,$_);
		} elsif ($field_count == 4){ 
			$field_count = 0;
			if ($read_qualities) {
				push(@qualities,$_);
			}
		} elsif (/^@/){
			push(@headers,substr($_,1));
			$field_count = 1;
		}
	}
	close FILE;

	if ($read_qualities) {
		return (\@sequences,\@headers,\@qualities);
	} else {
		return (\@sequences,\@headers);
	}

}

#################################################################################

# Reads a GFF/GFF3 annotation file
sub read_gff_file {

	my ($file,$source) = @_;

	my $gff_data;

	# print ("\nParsing GFF file '$file'.\n");

	if (is_gzip($file)){
		open(GFF_FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_zip($file)){
		open(GFF_FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_bzip2($file)){
		open(GFF_FILE, "bzcat '$file' |") or die "# cannot read $file\n";
	} else {
		open(GFF_FILE,$file) || die "# cannot read $file\n";
	}
	while (<GFF_FILE>) {
		if (/^#/){
			next;
		}
		s/\012\015?|\015\012?//; # trim line break
		my @gff_cols = split("\t");
		if ($#gff_cols < 8){
			next;
		}
		# If we desire only sequences from a specific source
		if (defined($source) && lc($gff_cols[1]) ne $source){
			next;
		}
		my $seq_id = $gff_cols[0];
		my $type = $gff_cols[2];
		if ($type eq 'gene' || $type eq 'mt_gene'){
			my @attrs = split( ";", $gff_cols[8] );
			# Only protein coding genes, exclude pseudogenes (biotype=processed_pseudogene)
			if ($attrs[2] !~ /protein_coding/){
				next;
			}
			if ($gff_cols[8] =~ /^ID=(gene:)?(\w+)/){
				push(@{$gff_data->{$seq_id}{'gene'}{$2}{'coords'}},[$gff_cols[3],$gff_cols[4]]);
				$gff_data->{$seq_id}{'gene'}{$2}{'strand'} = $gff_cols[6];
				# $gff_data->{$seq_id}{'gene'}{$2}{'start'} = $gff_cols[3];
				# $gff_data->{$seq_id}{'gene'}{$2}{'end'} = $gff_cols[4];
				# $gff_data->{$seq_id}{'gene'}{$2}{'source'} = $gff_cols[1];
				# $gff_data->{$seq_id}{'gene'}{$2}{'type'} = $gff_cols[2];
				# $gff_data->{$seq_id}{'gene'}{$2}{'score'} = $gff_cols[5];
				# $gff_data->{$seq_id}{'gene'}{$2}{'phase'} = $gff_cols[7];
				# $gff_data->{$seq_id}{'gene'}{$2}{'attributes'} = $gff_cols[8];
			}
		} elsif ($type eq 'mRNA' || $type eq 'transcript'){
			if ($gff_cols[8] =~ /^ID=(transcript:)?(\w+)/){
				push(@{$gff_data->{$seq_id}{'mRNA'}{$2}{'coords'}},[$gff_cols[3],$gff_cols[4]]);
				if (!defined($gff_data->{$seq_id}{'mRNA'}{$2}{'strand'})){
					$gff_data->{$seq_id}{'mRNA'}{$2}{'strand'} = $gff_cols[6];
				}
			}
		} elsif ($type eq 'CDS' ) {
			if ($gff_cols[8] =~ /Parent=(transcript:)?(\w+)/){
				push(@{$gff_data->{$seq_id}{'CDS'}{$2}{'coords'}},[$gff_cols[3],$gff_cols[4]]);
				if (!defined($gff_data->{$seq_id}{'CDS'}{$2}{'strand'})){
					$gff_data->{$seq_id}{'CDS'}{$2}{'strand'} = $gff_cols[6];
				}
				push(@{$gff_data->{$seq_id}{'cDNA'}{$2}{'coords'}},[$gff_cols[3],$gff_cols[4]]);
				if (!defined($gff_data->{$seq_id}{'cDNA'}{$2}{'strand'})){
					$gff_data->{$seq_id}{'cDNA'}{$2}{'strand'} = $gff_cols[6];
				}
			}
		} elsif ($type eq 'exon' ) {
			if ($gff_cols[8] =~ /Parent=(transcript:)?(\w+)/){
				push(@{$gff_data->{$seq_id}{'exon'}{$2}{'coords'}},[$gff_cols[3],$gff_cols[4]]);
				if (!defined($gff_data->{$seq_id}{'exon'}{$2}{'strand'})){
					$gff_data->{$seq_id}{'exon'}{$2}{'strand'} = $gff_cols[6];
				}
			}
		} elsif ($type eq 'five_prime_UTR' || $type eq 'three_prime_UTR') {
			if ($gff_cols[8] =~ /Parent=(transcript:)?(\w+)/){
				push(@{$gff_data->{$seq_id}{'cDNA'}{$2}{'coords'}},[$gff_cols[3],$gff_cols[4]]);
				if (!defined($gff_data->{$seq_id}{'cDNA'}{$2}{'strand'})){
					$gff_data->{$seq_id}{'cDNA'}{$2}{'strand'} = $gff_cols[6];
				}
			}
		}
	}
	close GFF_FILE;
	
	return $gff_data;

}

#################################################################################

# Creates a FASTA format file from a list of sequences
sub create_fasta_file {

	my ($seqs,$headers,$outfile,$compress_type) = @_;

	if (!defined($outfile)){
		$outfile = "/tmp/".random_file_name().".fasta";
	} elsif ($outfile =~ /\.zip$/){
		$compress_type = 'zip';
	} elsif ($outfile =~ /\.(gzip|gz)$/){
		$compress_type = 'gzip';
	}

	if (!defined($compress_type)){
		open(FASTAFILE,">$outfile") || die "\nERROR 'create_fasta_file': cannot create '$outfile'\n\n";
		for (my $i=0; $i<=$#{$seqs}; $i++){
			if (defined($headers) && defined($headers->[$i])){
				printf FASTAFILE ">%s\n%s\n" , $headers->[$i], $seqs->[$i];
			} else {
				printf FASTAFILE ">SEQ%d\n%s\n" , $i+1, $seqs->[$i];
			}
		}
		close FASTAFILE;
	} else {
		my $FASTAFILE;
		if ($compress_type eq 'zip'){
			if ($outfile !~ /\.zip$/){
				$outfile .= ".zip";
			}
			$FASTAFILE = new IO::Compress::Zip("$outfile") || die "\nERROR 'create_fasta_file': cannot create '$outfile'\n\n";
		} else {
			if ($outfile !~ /\.(gzip|gz)$/){
				$outfile .= ".gz";
			}
			$FASTAFILE = new IO::Compress::Gzip("$outfile") || die "\nERROR 'create_fasta_file': cannot create '$outfile'\n\n";
		}
		for (my $i=0; $i<=$#{$seqs}; $i++){
			if (defined($headers) && defined($headers->[$i])){
				$FASTAFILE->printf(">%s\n%s\n" , $headers->[$i], $seqs->[$i]);
			} else {
				$FASTAFILE->printf(">SEQ%d\n%s\n" , $i+1, $seqs->[$i]);
			}
		}
		$FASTAFILE->close();
	}

	return $outfile;

}

#################################################################################

# Creates a FASTA format file from a list of sequences
sub create_fastq_file {

	my ($seqs,$headers,$quals,$outfile,$compress_type) = @_;

	if (!defined($outfile)){
		$outfile = "/tmp/".random_file_name().".fastq";
	}

	if (!defined($compress_type)){
		open(FASTQFILE,">$outfile");
		for (my $i=0; $i<=$#{$seqs}; $i++){
			if (defined($quals) && defined($quals->[$i])){
				print FASTQFILE "@".$headers->[$i]."\n".$seqs->[$i]."\n+\n".$quals->[$i]."\n";
			} else {
				print FASTQFILE "@".$headers->[$i]."\n".$seqs->[$i]."\n+\n\n";
			}
		}
		close FASTQFILE;
	} else {
		my $FASTQFILE;
		if ($compress_type eq 'zip'){
			if ($outfile !~ /\.zip$/){
				$outfile .= ".zip";
			}
			$FASTQFILE = new IO::Compress::Zip("$outfile") || die "\nERROR 'create_fastq_file': cannot create '$outfile'\n\n";
		} else {
			if ($outfile !~ /\.(gzip|gz)$/){
				$outfile .= ".gz";
			}
			$FASTQFILE = new IO::Compress::Gzip("$outfile") || die "\nERROR 'create_fastq_file': cannot create '$outfile'\n\n";
		}
		for (my $i=0; $i<=$#{$seqs}; $i++){
			if (defined($quals) && defined($quals->[$i])){
				$FASTQFILE->printf("@%s\n%s\n+\n%s\n" , $headers->[$i], $seqs->[$i], $quals->[$i]);
			} else {
				$FASTQFILE->printf("@%s\n%s\n+\n\n" , $headers->[$i], $seqs->[$i]);
			}
		}
		$FASTQFILE->close();
	}

	return $outfile;

}

################################################################################

# Creates a hash from a FASTA file with names as keys and sequences as values
sub read_fasta_file_hash {

	my($fasta_file,$only_id) = @_;

	my ($sequences,$headers) = read_fasta_file($fasta_file,$only_id);

	my %sequences_hash;
	for (my $i=0; $i<=$#{$sequences}; $i++) {
		if (!defined($sequences_hash{$headers->[$i]})){
			$sequences_hash{$headers->[$i]} = $sequences->[$i];
		} else {
			print "\tWARNING 'read_fasta_file_hash': Header '".$headers->[$i]."' is duplicated.\n";
		}
	}

	return \%sequences_hash;
}

################################################################################

# Creates a hash from a FASTQ file with names as keys and sequences as values
sub read_fastq_file_hash {

	my($fastq_file) = @_;

	my ($sequences,$headers) = read_fastq_file($fastq_file);

	my %sequences_hash;
	for (my $i=0; $i<=$#{$sequences}; $i++) {
		if (!defined($sequences_hash{$headers->[$i]})){
			$sequences_hash{$headers->[$i]} = $sequences->[$i];
		} else {
			print "\tWARNING 'read_fastq_file_hash': Header '".$headers->[$i]."' is duplicated.\n";
		}
	}

	return \%sequences_hash;
}

################################################################################

# Converts a hash of sequences into two arrays of sequences and headers
sub sequences_hash_to_array {

	my $seqs_hash = shift @_;
	
	my (@sequences,@headers);
	foreach my $header (sort {$a cmp $b} keys %{$seqs_hash}){
		push(@headers,$header);
		push(@sequences,$seqs_hash->{$header});
	}

	return (\@sequences,\@headers);
}

################################################################################

# Converts two arrays of sequences and headers into a hash of sequences
sub sequences_array_to_hash {

	my ($seqs,$headers) = @_;
	
	my %sequences_hash;
	for (my $i=0; $i<=$#{$headers}; $i++) {
		if (!defined($sequences_hash{$headers->[$i]})) {
			$sequences_hash{$headers->[$i]} = $seqs->[$i];
		} else {
			print "\tWARNING 'sequences_array_to_hash': Header '".$headers->[$i]."' is duplicated.\n";
		}
	}

	return (\%sequences_hash);
}

#################################################################################

# Shuffles sequences from a FASTA format file
# Each sequence must be one line
sub shuffle_fasta_file {

	my ($file,$nseqs,$outfile,$compress_type) = @_;

	if (!defined($outfile)){
		$outfile = "/tmp/".random_file_name().'.fa';
	}

	my $cat_exe = 'cat';
	if (is_gzip($file) || is_zip($file)){
		$cat_exe = 'zcat';
	} elsif (is_bzip2($file)) {
		$cat_exe = 'bzcat';
	}

	my $command = "$cat_exe $file".' | awk \'{if ((NR%2)==0)print prev"X#&X"$0;prev=$0;}\' | shuf | ';

	if (defined($nseqs)){
		$command .= 'head -'.$nseqs.' | ';
	}
	
	$command .= 'sed \'s/X#&X/\n/g\' ';

	if (defined($compress_type)){
		$command .= '| '."$compress_type > $outfile";
	} else {
		$command .= '> '.$outfile;
	}
# 	print "$command\n";exit;
	`$command`;

	return $outfile;

}

#################################################################################

# Shuffles sequences from a FASTQ format file
sub shuffle_fastq_file {

	my ($file,$nseqs,$outfile,$compress_type) = @_;

	if (!defined($outfile)){
		$outfile = "/tmp/".random_file_name().'.fq';
	}

	my $cat_exe = 'cat';
	if (is_gzip($file) || is_zip($file)){
		$cat_exe = 'zcat';
	} elsif (is_bzip2($file)) {
		$cat_exe = 'bzcat';
	}

	my $command = "$cat_exe $file".' | awk \'{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("X#&X");} }\' | shuf | ';

	if (defined($nseqs)){
		$command .= 'head -'.$nseqs.' | ';
	}
	
	$command .= 'sed \'s/X#&X/\n/g\' ';

	if (defined($compress_type)){
		$command .= '| '."$compress_type > $outfile";
	} else {
		$command .= '> '.$outfile;
	}
# 	print "$command\n";exit;
	`$command`;

	return $outfile;

}


#################################################################################

# Shuffles sequences from two FASTQ format files
sub shuffle_fastq_files {

	my ($files,$nseqs,$outfiles,$compress_type) = @_;

	if (!defined($outfiles)){
		my $random_name = random_file_name();
		$outfiles = [ "/tmp/$random_name\_R1.fq", "/tmp/$random_name\_R2.fq" ];
	}

	my $subsample_exe = 'fq_subsample.py';

	my $command = "$subsample_exe -n $nseqs ".join(" ",@$files)." ".join(" ",@$outfiles);

# 	print "$command\n";exit;
	`$command`;
	
	if (defined($compress_type)){
		foreach my $outfile (@$outfiles) {
			if (!is_compressed($outfile)){
				$outfile = compress($outfile,$compress_type,$outfile.".$compress_type");
			}
		}
	}

	return $outfiles;

}

#################################################################################

# Splits a FASTA format file into several files
# Each sequence must be one line
sub split_fasta_file {

	my ($file,$parts,$outfile,$compress_type) = @_;

	if (!defined($outfile)){
		$outfile = "/tmp/".random_file_name();
	}
	
	my $total_lines;
	if (is_gzip($file) || is_zip($file)){
		$total_lines = `zcat $file | wc -l`;
	} elsif (is_bzip2($file)) {
		$total_lines = `bzcat $file | wc -l`;
	} else {
		$total_lines = `cat $file | wc -l`;
	}
	my $part_lines = int($total_lines/$parts);
	while ($part_lines % 2 !=0){
		$part_lines++;
	}

	return split_file($file,$part_lines,$outfile,$compress_type);

}

#################################################################################

# Splits a FASTQ format file into several files
sub split_fastq_file {

	my ($file,$parts,$outfile,$compress_type) = @_;

	if (!defined($outfile)){
		$outfile = "/tmp/".random_file_name();
	}
	
	my $total_lines;
	if (is_gzip($file) || is_zip($file)){
		$total_lines = `zcat $file | wc -l`;
	} elsif (is_bzip2($file)) {
		$total_lines = `bzcat $file | wc -l`;
	} else {
		$total_lines = `cat $file | wc -l`;
	}
	my $part_lines = int($total_lines/$parts);
	while ($part_lines % 4 !=0){
		$part_lines++;
	}

	return split_file($file,$part_lines,$outfile,$compress_type);

}

#################################################################################

# Splits a file into several files with the number of specified lines
sub split_file {

	my ($file,$lines,$outfile,$compress_type) = @_;

	my @splitted_files;

	if (!defined($outfile)){
		$outfile = "/tmp/".random_file_name();
	}

	my $cat_exe = 'cat';
	if (is_gzip($file) || is_zip($file)){
		$cat_exe = 'zcat';
	} elsif (is_bzip2($file)) {
		$cat_exe = 'bzcat';
	}

	my $split_output;
	if (defined($compress_type)){
		$split_output= `$cat_exe $file | split -l $lines -d --verbose --filter='$compress_type > \$FILE' - $outfile.`;
# 		print "$cat_exe $file | split -l $lines -d --verbose --filter='$compress_type > \$FILE' - $outfile.\n"; exit;
	} else {
		$split_output= `$cat_exe $file | split -l $lines -d --verbose - $outfile.`;
# 		print "$cat_exe $file | split -l $lines -d --verbose - $outfile.\n"; exit;
	}

	while ($split_output =~ /($outfile\.\d+)/g){
		push(@splitted_files,$1);
	}

	return @splitted_files;

}

#################################################################################

# Returns an array with Phred quality score numeric values
sub convert_phred_ascii_to_numeric {

	my ($quality_ascii, $format) = @_;
	
	if (!defined($format)){
		$format = 'sanger';
	}
	
	my %formats = (
			'sanger' => '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI',
			'solexa' => '@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh',
			'illumina-old' => '@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh',
			'illumina' => '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ',
	);

	my @quality_numeric;
	foreach my $char (split('',$quality_ascii)){
		push(@quality_numeric,index($formats{$format}, $char));
	}

	return \@quality_numeric;

}

#################################################################################

# Returns an array with Phred quality score ASCII values
sub convert_phred_numeric_to_ascii {

	my ($quality_numeric, $format) = @_;
	
	if (!defined($format)){
		$format = 'sanger';
	}
	
	my %formats = (
			'sanger' => '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI',
			'solexa' => '@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh',
			'illumina-old' => '@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh',
			'illumina' => '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ',
	);

	my $quality_ascii;
	foreach my $qual_value (@$quality_numeric){
		$quality_ascii .= substr($formats{$format},$qual_value,1);
	}

	return $quality_ascii;

}


#################################################################################

# Takes a unique sequence of ascii Phred qualities and returns a number with the average
sub mean_phred_ascii_quality {

	my ($quality_ascii, $format) = @_;
	
	if (!defined($format)){
		$format = 'sanger';
	}	

	my $quality_numeric = convert_phred_ascii_to_numeric($quality_ascii, $format);
	
	return mean(@$quality_numeric);

}

#################################################################################

# Takes an array of Phred qualities from several sequences and returns a unique sequence with the average
sub mean_multiple_phred_ascii_qualities {

	my ($seqs_qualities_ascii, $format) = @_;
	
	my @sum_qualities;
	foreach my $quality_ascii (@$seqs_qualities_ascii) {
		my $quality_numeric = convert_phred_ascii_to_numeric($quality_ascii);
		map $sum_qualities[$_]+=$quality_numeric->[$_], 0..$#{$quality_numeric};
	}
	my $count = scalar @$seqs_qualities_ascii;
	@sum_qualities = map sprintf("%.0f", $_/$count), @sum_qualities;

	return convert_phred_numeric_to_ascii(\@sum_qualities, $format);

}

#################################################################################

# Retrieves a hash with the sequence lenghts as keys and the number of sequences of the same length as values
sub count_lengths {

	my $seqs = shift @_;

	my %count_lengths;
	foreach my $len (map length($_) , @$seqs) {
		$count_lengths{$len}++;
	}

	return \%count_lengths;

}

#################################################################################

# Randomizes/shuffles sequences and headers (and optional qualities)
sub shuffle_seqs {

	my ($seqs,$headers,$quals) = @_;

	my @order = shuffle 0..$#{$seqs};
	
	if (defined($headers) && defined($quals)){
		return ([map $seqs->[$_], @order], [map $headers->[$_], @order], [map $quals->[$_], @order]);
	} elsif (defined($headers)){
		return ([map $seqs->[$_], @order], [map $headers->[$_], @order]);
	} else {
		return [map $seqs->[$_], @order];
	}

}

#################################################################################

# Takes a subset from a bigger set of sequences and headers (and optional qualities)
sub splice_seqs {

	my ($size,$seqs,$headers,$quals) = @_;

	if ($size > scalar @{$seqs}) {
		$size = scalar @{$seqs};
	}
	if (defined($headers) && defined($quals)){
		return ([splice(@{$seqs},0,$size)],[splice(@{$headers},0,$size)],[splice(@{$quals},0,$size)]);
	} elsif (defined($headers)){
		return ([splice(@{$seqs},0,$size)],[splice(@{$headers},0,$size)]);
	} else {
		return ([splice(@{$seqs},0,$size)]);
	}

}

#################################################################################

# Extracts 'n' random sequences from a sequence file
sub extract_random_seqs_from_file {

	my ($file, $number) = @_;
	
	my ($format);

	if (is_fastq($file)){
		$format = 'fastq';
	} elsif (is_fasta($file)){
		$format = 'fasta';
	}
	
	if ($format eq 'fastq'){
		return extract_random_seqs_from_fastq($file,$number)
	} elsif ($format eq 'fasta'){
		return extract_random_seqs_from_fasta($file,$number);
	} else {
		print "\nERROR 'extract_random_seqs_from_file': File format not recognized.\n\n";
		exit;
	}

}

#################################################################################

# Extracts 'n' random sequences from big FASTA file
sub extract_random_seqs_from_fasta {

	my ($file, $number) = @_;
	
	my $outfile = "/tmp/".random_file_name().".fa";

	my $cat_exe = 'cat';
	if (is_gzip($file) || is_zip($file)){
		$cat_exe = 'zcat';
	} elsif (is_bzip2($file)) {
		$cat_exe = 'bzcat';
	}

	system($cat_exe.' "'.$file.'" | awk \'{if ((NR%2)==0)print prev"X#&X"$0;prev=$0;}\' | shuf | head -'.$number.' | sed \'s/X#&X/\\n/g\' > "'.$outfile.'"');
	
	my ($sequences,$headers) = read_fasta_file($outfile);

	`rm "$outfile"`;
	
	return ($sequences,$headers);

}

#################################################################################

# Extracts 'n' random sequences from big FASTQ file
sub extract_random_seqs_from_fastq {

	my ($file, $number, $read_qualities) = @_;
	
	my $outfile = "/tmp/".random_file_name().".fq";

	my $cat_exe = 'cat';
	if (is_gzip($file) || is_zip($file)){
		$cat_exe = 'zcat';
	} elsif (is_bzip2($file)) {
		$cat_exe = 'bzcat';
	}

# 	print "\n".'cat '.$file.' | awk \'{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("X#&X");} }\' | shuf | head -'.$number.' | sed \'s/X#&X/\\n/g\' > '.$outfile."\n\n";exit;
	system($cat_exe.' "'.$file.'" | awk \'{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("X#&X");} }\' | shuf | head -'.$number.' | sed \'s/X#&X/\\n/g\' > "'.$outfile.'"');
	
	my ($sequences,$headers,$qualities) = read_fastq_file($outfile,$read_qualities);

	`rm "$outfile"`;
	
	if (defined($read_qualities) && $read_qualities) {
		return ($sequences,$headers,$qualities);
	} else {
		return ($sequences,$headers);
	}

}

#################################################################################

# Reads annotations from a sequence header
sub extract_header_data {

	my $header = shift @_;
	
	my %header_data;

	my @header_cols = split(/\s*[\|]\s*/,$header);

	for (my $i=0; $i<=$#header_cols; $i ++) {
		my $col = $header_cols[$i];
		if ($i == 0){
			$header_data{'name'} = $col;
		} elsif ($col =~ /(\w+)=(.+)/) {
			$header_data{$1} = $2;
		}
	}

	return \%header_data;

}

#################################################################################

# Counts sequences from sequences file
sub count_seqs_from_file {

	my ($file) = @_;
	
	my ($format);

	if (is_fastq($file)){
		$format = 'fastq';
	} elsif (is_fasta($file)){
		$format = 'fasta';
	}
	
	if ($format eq 'fastq'){
		return count_seqs_from_fastq($file)
	} elsif ($format eq 'fasta'){
		return count_seqs_from_fasta($file);
	} else {
		print "\nERROR 'count_seqs_from_file': File format not recognized.\n\n";
		exit;
	}

	
}

#################################################################################

# Counts sequences from FASTA file
sub count_seqs_from_fasta {

	my ($file) = @_;

	my $output;

	if (is_gzip($file) || is_zip($file)) {
		$output = `zgrep -c '^>' "$file"`;
	} elsif (is_bzip2($file)) {
		$output = `bzgrep -c '^>' "$file"`;
	} else {
		$output = `grep -c '^>' "$file"`;
	}
	
	chomp $output;
	
	return $output;
	
}

#################################################################################

# Counts sequences from FASTQ file
sub count_seqs_from_fastq {

	my ($file) = @_;

	my $output;

# 	print "\n\ncat $file | echo \$((\`wc -l\`/4))\n";exit;
	if (is_gzip($file) || is_zip($file)) {
		$output = `zcat "$file" | echo \$((\`wc -l\`/4))`;
	} elsif (is_bzip2($file)) {
		$output = `bzcat "$file" | echo \$((\`wc -l\`/4))`;
	} else {
		$output = `cat "$file" | echo \$((\`wc -l\`/4))`;
	}
	
	chomp $output;
	
	return $output;
	
}


#################################################################################

# Translates a cDNA/RNA sequence to protein
sub dna_to_prot {

	my ($seq) = @_;

	my $seq2;

	my(%genetic_code) = (
		'AAA' => 'K', 'AAG' => 'K', # Lysine
		'AAC' => 'N', 'AAT' => 'N', # Asparagine
		'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T', # Threonine
		'AGA' => 'R', 'AGG' => 'R', # Arginine
		'AGC' => 'S', 'AGT' => 'S', # Serine
		'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I', # Isoleucine
		'ATG' => 'M', # Methionine
		'CAA' => 'Q', 'CAG' => 'Q', # Glutamine
		'CAC' => 'H', 'CAT' => 'H', # Histidine
		'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P', # Proline
		'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', # Arginine
		'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', # Leucine
		'GAA' => 'E', 'GAG' => 'E', # Glutamic Acid
		'GAC' => 'D', 'GAT' => 'D', # Aspartic Acid
		'GCA' => 'A', 'GCC' => 'A', 'GCG' => 'A', 'GCT' => 'A', # Alanine
		'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',  # Glycine
		'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V', # Valine
		'TAA' => '*', 'TAG' => '*',
		'TAC' => 'Y', 'TAT' => 'Y', # Tyrosine
		'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', # Serine
		'TGA' => '*', # Stop
		'TGC' => 'C', 'TGT' => 'C', # Cysteine
		'TGG' => 'W', # Tryptophan
		'TTA' => 'L', 'TTG' => 'L', # Leucine
		'TTC' => 'F', 'TTT' => 'F', # Phenylalanine
	);
	
	$seq = uc($seq);
	$seq =~ s/U/T/g;

	while ( $seq =~ /\w{3}/g ) {
		my $codon = $&;
		if (defined($genetic_code{$codon})){
			$seq2 .= $genetic_code{$codon};
		} elsif ($codon =~ /[ACGTMRWSYKVHDBNX]{3}/){
			$seq2 .= '['.join('', unique(map $genetic_code{$_} , unambiguous_dna_sequences($codon))).']';
		} else {
			 print "# dna_to_prot: ERROR, non defined codon '$&'\n"
		}
	}

	return $seq2;

}


#################################################################################

# Translates multiple cDNA/RNA to protein sequences
sub translate_seqs {

	my ($dna_seqs,$dna_headers,$translation_frame,$options) = @_;

	my ($nostop,$dna_out,$revcomp,$unique) = (0,0,0,0);
	my ($auto_frame,$all_frames,$one_frame,$dna_trim) = (0,0,0,0);
	if ($translation_frame eq 'auto'){
		$auto_frame = 1; # detects automatically ORFs
	} elsif ($translation_frame eq 'all'){
		$all_frames = 1; # translates DNA into the 3 reading frames
	} elsif (is_numeric($translation_frame) && $translation_frame>=1 && $translation_frame<=3){
		$one_frame = $translation_frame; # translates into the desired reading frame
	} elsif ($translation_frame eq 'dna trim'){
		$dna_trim = 1;  # outputs trimmed DNA sequences to the nearest ORF
		$dna_out = 1;
	} else {
		$auto_frame = 1;
	}

	if (in_array($options, 'nostop')){
		$nostop = 1; # excludes translations with stop codons
	}
	if (in_array($options, 'dna out')){
		$dna_out = 1; # retrieves DNA sequences of translated sequences
	}
	if (in_array($options, 'revcomp')){
		$revcomp = 1; # retrieves DNA sequences of translated sequences
	}
	if (in_array($options, 'unique')){
		$unique = 1; # retrieves only unique translated sequences (removes redundant)
	}

	my ($prot_seqs,$prot_headers,%annotated_seqs);
	my ($dna_out_seqs,$dna_out_headers);
	for (my $i=0; $i<=$#{$dna_seqs}; $i++){
# if ($dna_headers->[$i] =~ /f51d75dc20754ba80cb827de6658c3b8/){
# print '';
# }
		for (my $j=0; $j<=$revcomp; $j++){
			my $seq = $dna_seqs->[$i];
			if ($j==1){
				$seq = iupac_reverse_complementary($seq)
			}
			my $sense = '';
			if ($revcomp && $j==0){
				$sense = '+';
			} elsif ($revcomp && $j==1){
				$sense = '-';
			}
			if ($auto_frame || $dna_trim){
				my %prot_seqs;
				my %prot_stops;
				for (my $frame=0; $frame<=2; $frame++){
					my $dna_seq = substr($seq,$frame);
					if (!$dna_seq){
						last;
					}
					$prot_seqs{$frame} = dna_to_prot($dna_seq);
					# Counts the number of stop codons:
					if (defined($prot_seqs{$frame})){
						$prot_stops{$frame} = () = $prot_seqs{$frame} =~ /\*/gi;
					}
				}
				if (!%prot_stops){
					next;
				}
				my @ordered_frames = sort { $prot_stops{$a}<=>$prot_stops{$b} } keys %prot_stops;
				my $min_stops = $prot_stops{$ordered_frames[0]};
				if ($nostop){
					$min_stops = 0;
				}
				for (my $frame=0; $frame<=2; $frame++){
					my $dna_seq = substr($seq,$frame);
					# Annotates only ORFs with minimum number of stop codons
					if ($prot_stops{$frame} > $min_stops){
						next;
					}
					# Protein sequence in the correct ORFs
					unless ($unique && !defined($annotated_seqs{$prot_seqs{$frame}})){
						push(@$prot_seqs, $prot_seqs{$frame});
						push(@$prot_headers, sprintf("%s | RF=%s%d", $dna_headers->[$i], $sense, $frame+1));
						if ($dna_out){
							push(@$dna_out_seqs, $dna_seq);
						}
						$annotated_seqs{$prot_seqs{$frame}} = 1;
					}
				}
			} elsif ($all_frames){
				for (my $frame=0; $frame<=2; $frame++){
					if (!substr($seq,$frame)){
						last;
					}
					my $dna_seq = substr($seq,$frame);
					my $prot_seq = dna_to_prot($dna_seq);
					if (!defined($prot_seq)){
						next;
					}
					my $count_stops = () = $prot_seq =~ /\*/gi;
					if ($nostop && $count_stops){
						next;
					}
					unless ($unique && defined($annotated_seqs{$prot_seq})){
						push(@$prot_seqs, $prot_seq);
						push(@$prot_headers, sprintf("%s | RF=%s%d", $dna_headers->[$i], $sense, $frame+1));
						if ($dna_out){
							push(@$dna_out_seqs, $dna_seq);
						}
						$annotated_seqs{$prot_seq} = 1;
					}
				}
			} elsif ($one_frame){
				if (!substr($seq,$one_frame-1)){
					last;
				}
				my $dna_seq = substr($seq,$one_frame-1);
				my $prot_seq = dna_to_prot($dna_seq);
				if (!defined($prot_seq)){
					next;
				}
				my $count_stops = () = $prot_seq =~ /\*/gi;
				if ($nostop && $count_stops){
					next;
				}
				unless ($unique && defined($annotated_seqs{$prot_seq})){
					push(@$prot_seqs, $prot_seq);
					push(@$prot_headers, sprintf("%s | RF=%s%d", $dna_headers->[$i], $sense, $one_frame));
					if ($dna_out){
						push(@$dna_out_seqs, $dna_seq);
					}
					$annotated_seqs{$prot_seq} = 1;
				}
			} else {
				my $dna_seq = $seq;
				my $prot_seq = dna_to_prot($dna_seq);
				my $count_stops = () = $prot_seq =~ /\*/gi;
				if ($nostop && $count_stops){
					next;
				}
				unless ($unique && defined($annotated_seqs{$prot_seq})){
					push(@$prot_seqs, $prot_seq);
					push(@$prot_headers, $dna_headers->[$i]);
					if ($dna_out){
						push(@$dna_out_seqs, $dna_seq);
					}
					$annotated_seqs{$prot_seq} = 1;
				}
			}
		}
	}

	if (!$dna_out){
		return ($prot_seqs,$prot_headers);
	} else {
		return ($prot_seqs,$prot_headers,$dna_out_seqs);
	}

}

#################################################################################

# Retrieves most probable reading frames
sub retrieve_reading_frames {

	my ($seq,$revcomp) = @_;

	my @reading_frames;
	foreach my $frame (1..3){
		if (!defined($revcomp)){
			if (index(dna_to_prot(substr($seq,$frame-1)),'*') == -1) {
				push(@reading_frames,$frame);
			}
		} else {
			if (index(dna_to_prot(substr(iupac_reverse_complementary($seq),$frame-1)),'*') == -1) {
				push(@reading_frames,$frame);
			}
		}
	}
	
	return @reading_frames;

}

#################################################################################

# Checks if there are stop codons in a sequence
sub is_stop_codon {

	my ($seq,$frame) = @_;
	
	if (index(dna_to_prot(substr($seq,$frame-1)),'*') != -1) {
		return 1;
	} else {
		return 0;
	}
}

#################################################################################

# Retrieves all possible nucleotides for a IUPAC symbol
sub convert_nt_from_iupac {

	my $nt = uc($_[0]);

	my %nt_iupac = (
	'A' => ['A'],
	'C' => ['C'],
	'G' => ['G'],
	'T' => ['T'],
	'U' => ['T'],
	'M' => ['A', 'C'],
	'R' => ['A', 'G'],
	'W' => ['A', 'T'],
	'S' => ['C', 'G'],
	'Y' => ['C', 'T'],
	'K' => ['G', 'T'],
	'V' => ['A', 'C', 'G'],
	'H' => ['A', 'C', 'T'],
	'D' => ['A', 'G', 'T'],
	'B' => ['C', 'G', 'T'],
	'N' => ['A', 'C', 'G', 'T'],
	'X' => ['A', 'C', 'G', 'T']
	);

	if (defined($nt_iupac{$nt})){
		return @{$nt_iupac{$nt}};
	} else {
		return undef;
	}

}

#################################################################################

# Retrieves all possible aminoacids for a IUPAC symbol
sub convert_aa_from_iupac {

	my $aa = uc($_[0]);

	my %aa_iupac = (
	'A' => ['A'],
	'R' => ['R'],
	'N' => ['N'],
	'D' => ['D'],
	'C' => ['C'],
	'Q' => ['Q'],
	'E' => ['E'],
	'G' => ['G'],
	'H' => ['H'],
	'L' => ['L'],
	'I' => ['E'],
	'K' => ['K'],
	'M' => ['M'],
	'F' => ['F'],
	'P' => ['P'],
	'S' => ['S'],
	'T' => ['T'],
	'W' => ['W'],
	'Y' => ['Y'],
	'V' => ['V'],
	'O' => ['O'],
	'U' => ['U'],
	'J' => ['L', 'I'],
	'B' => ['D', 'N'],
	'Z' => ['E', 'Q'],
	'X' => ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'L', 'I', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
	);

	return @{$aa_iupac{$aa}};

}

#################################################################################

# Converts a set of nucleotides to a IUPAC symbol
sub convert_nt_to_iupac{

	my @nts = @_;

	my %nt_iupac = (
	'A' => ['A'],
	'C' => ['C'],
	'G' => ['G'],
	'T' => ['T'],
# 	'U' => ['T'],
	'M' => ['A', 'C'],
	'R' => ['A', 'G'],
	'W' => ['A', 'T'],
	'S' => ['C', 'G'],
	'Y' => ['C', 'T'],
	'K' => ['G', 'T'],
	'V' => ['A', 'C', 'G'],
	'H' => ['A', 'C', 'T'],
	'D' => ['A', 'G', 'T'],
	'B' => ['C', 'G', 'T'],
	'N' => ['A', 'C', 'G', 'T']
	#'X' => ['A', 'C', 'G', 'T']
	);

	my $hits=0;
	foreach my $iupac (keys(%nt_iupac)){
		$hits=0;
		foreach my $nt_iupac (@{$nt_iupac{$iupac}}){
			foreach my $nt (@nts) {
				if (uc($nt) eq 'U') { $nt = 'T' }
				if ($nt_iupac eq uc($nt)){
					$hits++;
				}
				if ($hits == scalar(@nts) && $hits == scalar(@{$nt_iupac{$iupac}})){
					return $iupac;
				}
			}
		}
	}
	
	return undef;

}

#################################################################################

# Converts a set of aminoacids to a IUPAC symbol
sub convert_aa_to_iupac{

	my @aas = @_;

	my %aa_iupac = (
	'A' => ['A'],
	'R' => ['R'],
	'N' => ['N'],
	'D' => ['D'],
	'C' => ['C'],
	'Q' => ['Q'],
	'E' => ['E'],
	'G' => ['G'],
	'H' => ['H'],
	'L' => ['L'],
	'I' => ['E'],
	'K' => ['K'],
	'M' => ['M'],
	'F' => ['F'],
	'P' => ['P'],
	'S' => ['S'],
	'T' => ['T'],
	'W' => ['W'],
	'Y' => ['Y'],
	'V' => ['V'],
	'O' => ['O'],
	'U' => ['U'],
	'J' => ['L', 'I'],
	'B' => ['D', 'N'],
	'Z' => ['E', 'Q'],
	'X' => ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'L', 'I', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
	);

	my $hits=0;
	foreach my $iupac (keys(%aa_iupac)){
		$hits=0;
		foreach my $aa_iupac (@{$aa_iupac{$iupac}}){
			foreach my $aa (@aas) {
				if ($aa_iupac eq uc($aa)){
					$hits++;
				}
				if ($hits == scalar(@aas) && $hits == scalar(@{$aa_iupac{$iupac}})){
					return $iupac;
				} elsif ($hits>2) {
					return 'X';
				}
			}
		}
	}
	
	return undef;
}

#################################################################################

# Converts 3 letter aminoacid codes to 1 letter ones
sub convert_aa_3_to_1 {

	my $aa = ucfirst(lc($_[0]));

	my %aa_iupac = (
	'Ala' => 'A',
	'Arg' => 'R',
	'Asn' => 'N',
	'Asp' => 'D',
	'Cys' => 'C',
	'Gln' => 'Q',
	'Glu' => 'E',
	'Gly' => 'G',
	'His' => 'H',
	'Leu' => 'L',
	'Ile' => 'I',
	'Lys' => 'K',
	'Met' => 'M',
	'Mse' => 'M', # Seleno metionine
	'Phe' => 'F',
	'Pro' => 'P',
	'Ser' => 'S',
	'Thr' => 'T',
	'Trp' => 'W',
	'Tyr' => 'Y',
	'Val' => 'V',
	'Pyl' => 'O',
	'Sec' => 'U',
	'Xle' => 'J',
	'Asx' => 'B',
	'Glx' => 'Z',
	'Xaa' => 'X',
	'Ter' => '*',
	);

	return $aa_iupac{$aa};

}

#################################################################################

# Converts 1 letter aminoacid codes to 3 letter ones
sub convert_aa_1_to_3 {

	my $aa = uc($_[0]);

	my %aa_iupac = reverse (
	'Ala' => 'A',
	'Arg' => 'R',
	'Asn' => 'N',
	'Asp' => 'D',
	'Cys' => 'C',
	'Gln' => 'Q',
	'Glu' => 'E',
	'Gly' => 'G',
	'His' => 'H',
	'Leu' => 'L',
	'Ile' => 'I',
	'Lys' => 'K',
	'Met' => 'M',
	'Phe' => 'F',
	'Pro' => 'P',
	'Ser' => 'S',
	'Thr' => 'T',
	'Trp' => 'W',
	'Tyr' => 'Y',
	'Val' => 'V',
	'Pyl' => 'O',
	'Sec' => 'U',
	'Xle' => 'J',
	'Asx' => 'B',
	'Glx' => 'Z',
	'Xaa' => 'X',
	'Ter' => '*',
	);

	return $aa_iupac{$aa};

}

#################################################################################

# Clean dashes and spaces from sequence
sub clean_sequence {

	my ($seq) = @_;

	my $seq2 = $seq;
	
	$seq2 =~ s/[-\s,\.]//g;
	
	return $seq2;
	
}


#################################################################################

# Converts a IUPAC sequence in its complementary
sub complementary_sequence {

	my ($seq) = @_;

	my $seq2 = '';

	my %nt_iupac_comp = (
	'A' => 'T',
	'C' => 'G',
	'G' => 'C',
	'T' => 'A',
	'U' => 'A',
	'M' => 'K',
	'R' => 'Y',
	'W' => 'W',
	'S' => 'S',
	'Y' => 'R',
	'K' => 'M',
	'V' => 'B',
	'H' => 'D',
	'D' => 'H',
	'B' => 'V',
	'X' => 'X',
	'N' => 'N',
	'?' => '?',
	'-' => '-'
	);

	if ($seq){
		foreach my $nt (split(//, $seq)){
			if ($nt =~ /[\s\-_]/){
				$seq2 .= $nt;
			} elsif (defined($nt_iupac_comp{uc($nt)})){
				if ($nt =~ /[a-z]/){
					$seq2 .= lc($nt_iupac_comp{uc($nt)});
				} elsif ($nt =~ /[A-Z]/){
					$seq2 .= uc($nt_iupac_comp{uc($nt)});
				}
			} else {
				$seq2 .= $nt;
				print "# ERROR in \'$seq\',\n  nucleotide \'$nt\' doesn't have a complementary one.\n";
			}
		}
	}

	return $seq2;

}

#################################################################################

# Obtains the reverse sequence
sub reverse_sequence {

	my ($seq) = @_;

	my $seq2 = '';

	foreach my $nt (reverse(split(//, $seq))){
		$seq2 .= $nt;
	}

	return $seq2;
}

#################################################################################

# Obtains the reverse complementary IUPAC motif consensus
sub iupac_reverse_complementary {

	my ($iupac) = @_;

	return complementary_sequence(reverse_sequence($iupac));

}

#################################################################################

# Obtains the consensus from a group of sequences
sub consensus_sequence {

	my ($seqs,$depths,$type) = @_;

	my $consensus_seq;

	my ($is_dna, $is_prot) = (0, 0);
	if (!defined($type) && !defined($depths)){
		$is_dna = is_dna(join'',@$seqs);
		$is_prot = is_prot(join'',@$seqs);
		if ($is_dna + $is_prot == 0) {
			print "\nERROR 'consensus_sequence': Unrecognized sequence type.\n\n";
			return undef;
		} elsif ($is_dna + $is_prot == 2) {
			$is_dna = 1;
		}
	} elsif ($type =~ /dna|rna/){
		$is_dna = 1;
	} elsif ($type =~ /pro|aa|pep/){
		$is_prot = 1;
	}

	my ($aligned_seqs,$aligned_names) = multiple_align_seqs($seqs);

	my $position_depth_data;
	for (my $i=0; $i<=$#{$aligned_seqs}; $i++){
		my @seq = split('', $aligned_seqs->[$i]);
		my $depth = 1;
		if (defined($depths)){
			$depth = $depths->[$i];
		}
		for (my $j=0; $j<=$#seq; $j++) {
			my $nt = uc($seq[$j]);
			if (is_numeric($depth)){
				if ($nt =~ /[ACGTU]/i){
					$position_depth_data->[$j]{$nt} += $depth;
				} elsif ($nt =~ /[MRWSYKVHDBNX]/i){
					my @nts_ = convert_nt_from_iupac($nt);
					foreach my $nt_ (@nts_){
						$position_depth_data->[$j]{uc($nt_)} += $depth/scalar(@nts_);
					}
				}
			}
		}
	}

	my $consensus_aligned_seq;
	foreach my $position_depths (@{$position_depth_data}){
		if (defined($depths)){
			my $max_nt;
			my $max_depth = -1;
			while ((my $nt, my $depth) = each %$position_depths) {
				if ($depth > $max_depth) {
					$max_depth = $depth;
					$max_nt = $nt;
				}
			}
			$consensus_aligned_seq .= $max_nt;
		} elsif ($is_dna) {
			$consensus_aligned_seq .= convert_nt_to_iupac(keys %$position_depths);
		} elsif ($is_prot) {
			my @aas = keys %$position_depths;
			my $con = convert_aa_to_iupac(@aas);
			if (!defined($con) || $con eq 'X'){
				if (scalar @aas < 4){
					$con = '('.join(',',@aas).')';
				} else {
					$con = 'X';
				}
			}
			$consensus_aligned_seq .= $con;
		}
	}
	$consensus_seq = uc($consensus_aligned_seq);
	$consensus_seq =~ s/-//g;

	return $consensus_seq;

}

#################################################################################

# Creates a random sequence
sub random_sequence {
	my ($length,$type) = @_;

	if (!defined($type)){
		$type = 'dna';
	}

	my $letters = {
		'dna' => [ 'A', 'C', 'G', 'T' ],
		'rna' => [ 'A', 'C', 'G', 'U' ],
		'prot' => [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ],
	};
	

	my $seq = '';
	
	for (my $i=1; $i<=$length; $i++){
		$seq .= $letters->{$type}[rand @{$letters->{$type}}];
	}

	return $seq;
}

#################################################################################

# Extract all possible sequence combinations from ambiguous sequences
sub unambiguous_dna_sequences {

	my $seq=$_[0];

	my @nts = split(//, $seq);
	my @comb_seqs;
	my @seqs;
	my $count_combs=number_nt_combs($seq);
	my $index_combs=$count_combs;

	foreach my $nt (@nts){
		my @comb_nts = convert_nt_from_iupac($nt);
		my $count_nt_combs=0;
		my $index=0;
		$index_combs=$index_combs/scalar(@comb_nts);
		for (my $i=0;$i<$count_combs;$i++){
			$count_nt_combs++;
			push(@{$comb_seqs[$i]},$comb_nts[$index]);
			if ($count_nt_combs==$index_combs){
				$count_nt_combs=0;
				$index++;
				if ($index>$#comb_nts){$index=0;}
			}
		}
	}

	foreach my $comb_seq (@comb_seqs){
		push(@seqs,join('',@{$comb_seq}));
	}

	return @seqs;


}

#################################################################################

# Counts the number of possible sequence combinations from ambiguous sequences
sub number_nt_combs {

	my $seq=$_[0];

	my $count_combs = 1;

	my @nts = split(//, $seq);

	foreach my $nt (@nts){
		if ($nt =~ /[M|R|W|S|Y|K]/i) {$count_combs *= 2;}
		if ($nt =~ /[V|H|D|B]/i) {$count_combs *= 3;}
		if ($nt =~ /[X|N]/i) {$count_combs *= 4;}
	}

	return $count_combs;

}


#################################################################################

# Retrieves REGEX for IUPAC symbols
sub regex {

	my ($seq) = @_;

	my $mult_iupac = {
	'M' => ['A', 'C'], 'm' => ['a', 'c'],
	'R' => ['A', 'G'], 'r' => ['a', 'g'],
	'W' => ['A', 'T'], 'w' => ['a', 't'],
	'S' => ['C', 'G'], 's' => ['c', 'g'],
	'Y' => ['C', 'T'], 'y' => ['c', 't'],
	'K' => ['G', 'T'], 'k' => ['g', 't'],
	'V' => ['A', 'C', 'G'], 'v' => ['a', 'c', 'g'],
	'H' => ['A', 'C', 'T'], 'h' => ['a', 'c', 't'],
	'D' => ['A', 'G', 'T'], 'd' => ['a', 'g', 't'],
	'B' => ['C', 'G', 'T'], 'b' => ['c', 'g', 't'],
	'N' => ['A', 'C', 'G', 'T'], 'n' => ['a', 'c', 'g', 't'],
	'X' => ['A', 'C', 'G', 'T'], 'x' => ['a', 'c', 'g', 't']
	};
	
	my $regex = '';
	my @regex_ = ('');
	foreach my $nt (split('', $seq)) {
		my $regex_;
		if (!defined($mult_iupac->{$nt})){
			$regex_ = $nt;
		} else {
			$regex_ = '['.join('', @{$mult_iupac->{$nt}}).']';
		}
		if ($regex_ ne $regex_[0]){
			if (scalar @regex_ <= 2){
				$regex .= join('',@regex_);
			} else {
				$regex .= sprintf("%s{%d}",$regex_[0],scalar @regex_);
			}
			undef(@regex_);
		}
		push(@regex_,$regex_);
	}
	# Last positions:
	if (scalar @regex_ <= 2){
		$regex .= join('',@regex_);
	} else {
		$regex .= sprintf("%s{%d}",$regex_[0],scalar @regex_);
	}
	
	return $regex;
	
# 	# OBSOLETE:
# 	my $regex = '';
# 	foreach my $nt (split('', $seq)) {
# 		if (!defined($mult_iupac->{$nt})){
# 			$regex .= $nt;
# 		} else {
# 			$regex .= '['.join('', @{$mult_iupac->{$nt}}).']';
# 		}
# 	}

}

#################################################################################

# Retrieves REGEX from PROSITE pattern
sub prosite_to_regex {

	my ($prosite_pattern) = @_;
	
	# $pattern = '< A-x-[ST](2)-x(0,1)-V';
	
	my %prosite_to_regex = (
		'<' => '^',
		'>' => '$',
		'-' => '',
		' ' => '',
		'x' => '\w',
		'X' => '\w',
		'{' => '[^',
		'}' => ']',
		'(' => '{',
		')' => '}',
		'*' => '.+',
	);
	
	my $regex_pattern = '';
	my $ambiguous = 0;
	foreach my $pos (split('',$prosite_pattern)){
		if ($pos eq '['){
			$ambiguous = 1;
		} elsif ($pos eq ']'){
			$ambiguous = 0;
		} elsif ($pos =~ /[JBZ]/i) {
			if ($ambiguous){
				$pos = join('',convert_aa_from_iupac($pos));
			} else {
				$pos = '['.join('',convert_aa_from_iupac($pos)).']';
			}
		} elsif (defined($prosite_to_regex{$pos})){
			$pos = $prosite_to_regex{$pos};
		}
		$regex_pattern .= $pos;
	}

# 	while (my ($prosite, $regex) = each %prosite_to_regex) {
# 		$regex_pattern =~ s/quotemeta($prosite)/$regex/;
# 	}
	
	return $regex_pattern;
	
}

#################################################################################

# Retrieves PROSITE pattern from REGEX
sub regex_to_prosite {

	my ($regex_pattern) = @_;
	
	# $pattern = '^A\w[ST]{2}\w{0,1}V$';
	
	my %regex_to_prosite = (
		'^' => '<',
		'$' => '>',
		'-' => '',
		'\s' => '',
		'\w' => 'x',
		'{' => '(',
		'}' => ')',
		'[^' => '{',
		']' => '}',
		'.+' => '*',
		);

	my $prosite_pattern = '';
	my $except = 0;
	my $pos = '';
	foreach my $pos_ (split('',$regex_pattern)){
		if ($pos && defined($regex_to_prosite{$pos.$pos_})){
			if ($pos.$pos_ eq '[^'){
				$except = 1;
			}
			$prosite_pattern .= $regex_to_prosite{$pos.$pos_};
			$pos = '';
		} elsif (defined($regex_to_prosite{$pos_})){
			if ($pos_ eq ']' && !$except){
				$prosite_pattern .= $pos.$pos_;
			} else {
				$prosite_pattern .= $pos.$regex_to_prosite{$pos_};
				$except = 0;
			}
			$pos = '';
		} elsif ($pos) {
			$prosite_pattern .= $pos;
			$pos = $pos_;
		} else {
			$pos .= $pos_;
		}
	}
	$prosite_pattern .= $pos;

# 	my $prosite_pattern = $regex_pattern;
# 	while (my ($prosite, $regex) = each %regex_to_prosite) {
# 		$prosite = quotemeta($prosite);
# 		$prosite_pattern =~ s/$prosite/$regex/g;
# 	}
	
	return $prosite_pattern;
	
}


#################################################################################

# Score a DNA sequence with one and zero
sub binary_score_nts {

	my ($seq1,$seq2,$options)=@_;

	if (!defined($seq1) || !defined($seq2)){
		return;
	}

	my ($score,$total) = (0,0);

	if (length($seq1)==length($seq2)){

		my @seq1 = split(//,$seq1);
		my @seq2 = split(//,$seq2);

		for (my $i=0; $i<=$#seq1; $i++) {

			# Skip non IUPAC letters
			if ($seq1[$i] !~ /[ACGTUMRWSYKVHDBNX]/i && $seq2[$i] !~ /[ACGTUMRWSYKVHDBNX]/i) {
				next;
			}
			
			if (defined($options) && $options eq 'fine'){
				my @nts1 = convert_nt_from_iupac($seq1[$i]);
				my @nts2 = convert_nt_from_iupac($seq2[$i]);
				my $count_nts = 0;
				my $score_nts = 0;
				foreach my $nt1 (@nts1) {
					foreach my $nt2 (@nts2) {
						if ($nt1 eq $nt2){
							$score_nts++;
							last;
						}
					}
				}
				if (scalar @nts1 >= scalar @nts2) { $count_nts = scalar @nts1; } else { $count_nts = scalar @nts2; }
				$score += $score_nts/$count_nts;
				$total++;

			} else {

				if (lc($seq1[$i]) eq lc($seq2[$i])){
					$score++;
				}
				$total++;
			} 
		}

	} else  {
		return;
# 		die "\nERROR 'binary_score_nts': Both sequences must have the same length to score the alignment\n\n";
	}

	return ($score,$total);

}

#################################################################################

# Detect errors in a sequence compared to a reference one
sub detect_sequence_errors {

	my ($seq,$seq_ref,$options)=@_;

	if (!defined($seq) || !defined($seq_ref)){
		return;
	}

	my ($substitutions, $insertions, $deletions, $homopolymer_indels, $insertions_ref) = ([],[],[],[],[]);

	my (@homopolymers, @homopolymers_ref);

	if (length($seq)==length($seq_ref)){

		my @seq = split(//,$seq);
		my @seq_ref = split(//,$seq_ref);

		my ($pos, $pos_ref) = (0, 0);
		my ($last, $last_ref) = ('', '');
		my (@homopolymer, @homopolymer_ref);
		for (my $i=0; $i<=$#seq; $i++) {
			if ($seq[$i] ne '-'){
				$pos++;
				if ($seq[$i] ne $last){
					if (scalar @homopolymer >= 3){
						push(@homopolymers,@homopolymer);
					}
					undef @homopolymer;
				}
				push(@homopolymer,$pos);
				$last = $seq[$i];
			}
			if ($seq_ref[$i] ne '-'){
				$pos_ref++;
				if ($seq_ref[$i] ne $last_ref){
					if (scalar @homopolymer_ref >= 3){
						push(@homopolymers_ref,@homopolymer_ref);
					}
					undef @homopolymer_ref;
				}
				push(@homopolymer_ref,$pos_ref);
				$last_ref = $seq_ref[$i];
			}
			if (lc($seq[$i]) ne lc($seq_ref[$i])){
				if ($seq[$i] eq '-'){
					push(@$deletions, $pos+1);
					push(@$insertions_ref, $pos_ref);
				} elsif ($seq_ref[$i] eq '-'){
					push(@$insertions, $pos);
				} else {
					push(@$substitutions, $pos);
				}
			}
		}

	} else  {
		die "\nERROR 'detect_sequence_errors': Both sequences must have the same length to be compared.\n\n";
	}

	foreach my $pos (@$insertions){
		if (in_array(\@homopolymers,$pos)){
			push(@$homopolymer_indels, $pos);
		}
	}
	for (my $i=0; $i<=$#{$deletions}; $i++){
		my $pos = $deletions->[$i];
		my $pos_ref = $insertions_ref->[$i];
		if (in_array(\@homopolymers_ref,$pos_ref)){
			push(@$homopolymer_indels, $pos);
		}
	}

	return ($substitutions, $insertions, $deletions, $homopolymer_indels);

}

#################################################################################

# Inserts errors into a sequence
sub insert_sequence_errors {

	my ($seq,$substitutions,$insertions,$deletions)=@_;

	if (!defined($seq)){
		return;
	}

	my @nts = ( 'A', 'C', 'G', 'T' );
	
	my @seq = split(//,$seq);

	if (defined($substitutions) && @$substitutions) {
		foreach my $pos (@$substitutions){
			my $nt = $seq[$pos-1];
			while ($nt eq $seq[$pos-1]){
				$nt = $nts[rand @nts];
			}
			$seq[$pos-1] = $nt;
		}
	}

	if (defined($insertions) && @$insertions) {
		foreach my $pos (@$insertions){
			my $nt = $nts[rand @nts];
			splice(@seq, $pos-1, 0, $nt);
		}
	}

	if (defined($deletions) && @$deletions) {
		foreach my $pos (@$deletions){
			splice(@seq, $pos-1, 1);
		}
	}

	return join('',@seq);

}

#################################################################################

# Prints in an abbreviatted SAM like format sequence errors
sub print_sequence_errors {

	my ($substitutions, $insertions, $deletions, $homopolymer_indels) = @_;

	my $seq_errors = '';

	if (defined($substitutions) && scalar @$substitutions > 0){
		$seq_errors .= sprintf("X%d", scalar @$substitutions);
	}
	if (defined($insertions) && scalar @$insertions > 0){
		$seq_errors .= sprintf("I%d", scalar @$insertions);
	}
	if (defined($deletions) && scalar @$deletions > 0){
		$seq_errors .= sprintf("D%d", scalar @$deletions);
	}
	if (defined($homopolymer_indels) && scalar @$homopolymer_indels > 0){
		$seq_errors .= sprintf("H%d", scalar @$homopolymer_indels);
	}

	return $seq_errors;

# 	if (@$substitutions && @$insertions && @$deletions){
# 		return scalar @$substitutions."X".scalar @$insertions."I".scalar @$deletions."D";
# 	} elsif (@$substitutions && @$insertions){
# 		return scalar @$substitutions."X".scalar @$insertions."I";
# 	} elsif (@$substitutions && @$deletions){
# 		return scalar @$substitutions."X".scalar @$deletions."D";
# 	} elsif (@$insertions && @$deletions){
# 		return scalar @$insertions."I".scalar @$deletions."D";
# 	} elsif (@$substitutions){
# 		return scalar @$substitutions."X";
# 	} elsif (@$insertions){
# 		return scalar @$insertions."I";
# 	} elsif (@$deletions){
# 		return scalar @$deletions."D";
# 	}

}

#################################################################################

# Compares a set of reads-sequences looking for 1 substitution/indel differences
sub compare_sequences {

	my ($seq,$ref_seq)=@_;

	$seq = lc($seq);
	$ref_seq = lc($ref_seq);
	my @seq_nts = split('',$seq);
	my @ref_seq_nts = split('',$ref_seq);
	my %seq_errors;
	my $error = -1;
	for (my $j=0; $j<=$#seq_nts; $j++){
		if ($j<length($ref_seq) && $ref_seq_nts[$j] ne $seq_nts[$j]){
			$error=$j;
			last;
		}
	}
	if ($error>=0) {
		if ($ref_seq eq substr($seq,0,$error).substr($ref_seq,$error,1).substr($seq,$error) ) {
			return ('deletion', $error+1);
		} elsif ($ref_seq eq substr($seq,0,$error).substr($seq,$error+1) ) {
			return ('insertion',$error+1);
		} elsif ($ref_seq eq substr($seq,0,$error).substr($ref_seq,$error,1).substr($seq,$error+1) ) {
			return ('substitution',$error+1);
		}
	} else {
		return ('identical', undef);
	}
	
	return (undef);

}

#################################################################################

# Compares a set of reads/sequences looking for N substitutions
sub count_errors {

	my ($seq1,$seq2,$limit) = @_;

	if (!defined($seq1) || !defined($seq2) || length($seq1) != length($seq2)){
		return undef;
	}

	my $errors = 0;
# 	if (length($seq1)==length($seq2)){
	my @seq1 = split(//,$seq1);
	my @seq2 = split(//,$seq2);
	for (my $i=0; $i<=$#seq1; $i++) {
		if ($seq1[$i] ne $seq2[$i]){
			$errors++;
			if ($errors == $limit){
				last;
			}
		}
	}
	return $errors;
# 	} else  {
# 		die "\nERROR 'find_substitutions': Both sequences must have the same length to score the alignment\n\n";
# 	}

}

#################################################################################

# # Embedded C code that compares a set of reads/sequences looking for N substitutions
# # 15% faster than Perl implementation
# use Inline C => << 'EOC';
# 	int count_errors_c(char* x, char* y, int limit) {
# 		int i;
# 		int errors = 0;
# 		for(i=0; x[i] && y[i]; ++i) {
# 		if(x[i] != y[i]) {
# 			errors++;
# 			if (errors == limit){
# 				break;
# 			}
# 		}
# 		}
# 		return errors;
# 	}
# EOC

#################################################################################

# Takes a set of sequences and detects chimeras
sub detect_chimeras {

	my ($seqs,$min_chimera) = @_;

	# Detect only chimeras of more than a threshold
	# Threshold can be a number of nts or a percentage of sequence length
	if (!defined($min_chimera)){
		$min_chimera = 0; # 10%, 15
	}

	my $chimeras = {};
	my $partial_seqs;
	my $partial_seqs_rev;

	# Loops all the sequences
	for (my $i=0; $i<=$#{$seqs}; $i++){
		my $seq1 = $seqs->[$i];
		my $len_seq1 = length($seq1);
		# Define minimum length of a partial match in the chimera
		my $min_len;
		if ($min_chimera =~ /([\d\.]+)%/) {
			$min_len = int($len_seq1*$1/100);
		} else {
			$min_len = $min_chimera;
		}
		for (my $j=0; $j<=$#{$seqs}; $j++){
			if ($i == $j) { next; } # Skips identical sequences
			my $seq2 = $seqs->[$j];
			if ($seq1 eq $seq2) { next; } # Skips identical sequences
			my $len_seq2 = length($seq2);
			# Compare both sequences to annotate how many nts are in common in both sides
			# $diff_pos is the position of the first substitution (starting from 1...)
			my ($identical,$diff_pos) = compare_strings($seq1,$seq2);
			if ($identical == 0 && $diff_pos-1 >= $min_len){
				$partial_seqs->{$i}{$j} = $diff_pos-1;
			}
			($identical,$diff_pos) = compare_strings(scalar reverse($seq1), scalar reverse($seq2));
			if ($identical == 0 && $diff_pos-1 >= $min_len){
				$partial_seqs_rev->{$i}{$j} = $diff_pos-1;
			}
			# print '';
		}
	}

	# Loops all the sequences
	for (my $i=0; $i<=$#{$seqs}; $i++){
		my $len_seq = length($seqs->[$i]);
		# Checks sequences with common left side
		foreach my $j1 (keys %{$partial_seqs->{$i}}) {
			if (!defined($partial_seqs->{$i}{$j1})) { next; }
			# Checks sequences with common right side
			foreach my $j2 (keys %{$partial_seqs_rev->{$i}}) {
				if ($j1 == $j2 || !defined($partial_seqs_rev->{$i}{$j2})) { next; }
				# If both sequences can build the template one, annote positions and ranges
				if ( $partial_seqs->{$i}{$j1}+$partial_seqs_rev->{$i}{$j2} > $len_seq ){
					my %chimera = ('seq1'=>$j1, 'seq1_start'=>'1', 'seq1_end'=>sprintf("%d...%d", $len_seq-$partial_seqs_rev->{$i}{$j2},$partial_seqs->{$i}{$j1}), 
					'seq2'=>$j2, 'seq2_start'=>sprintf("%d...%d", $len_seq-$partial_seqs_rev->{$i}{$j2}+1, $partial_seqs->{$i}{$j1}+1), 'seq2_end'=>$len_seq);
					push(@{$chimeras->{$i}},\%chimera);
				# If there is only one combination of both sequences to give the template
				} elsif ( $partial_seqs->{$i}{$j1}+$partial_seqs_rev->{$i}{$j2} == $len_seq ){
					my %chimera = ('seq1'=>$j1, 'seq1_start'=>'1', 'seq1_end'=>$partial_seqs->{$i}{$j1}, 
					'seq2'=>$j2, 'seq2_start'=>$partial_seqs->{$i}{$j1}+1, 'seq2_end'=>$len_seq);
					push(@{$chimeras->{$i}},\%chimera);
				}
			}
		}
	}
	return $chimeras;

# WITH ERRORS AND SLOW:
# 	my $chimeras = {};
# 	for (my $i=0; $i<=$#{$seqs}; $i++){
# 		my $test_seq = $seqs->[$i];
# 		my @test_seq_nts = split('',$test_seq);
# 		my $min_len = int(length($test_seq)*$min_chimera/100);
# 		for (my $j=0; $j<=$#{$seqs}; $j++){
# 			my $seq1 = $seqs->[$j];
# 			if ($test_seq eq $seq1) { next;}
# 			my @seq1_nts = split('',$seq1);
# 			for (my $p1=0; $p1<=$#test_seq_nts; $p1++){
# 				if ($test_seq_nts[$p1] ne $seq1_nts[$p1]){
# 					for (my $k=0; $k<=$#{$seqs}; $k++){
# 						my $seq2 = $seqs->[$k];
# 						if ($test_seq eq $seq2 || $seq1 eq $seq2) { next;}
# 						my @seq2_nts = split('',$seq2);
# 						for (my $p2=$p1; $p2<=$#test_seq_nts; $p2++){
# 							if ($p2<=$#seq2_nts && $test_seq_nts[$p2] ne $seq2_nts[$p2]){
# 								last;
# 							} elsif ($p2==$#test_seq_nts) {
# 								my %chimera = ('seq1'=>$j, 'seq1_start'=>1, 'seq1_end'=>$p1, 'seq2'=>$k, 'seq2_start'=>$p1+1, 'seq2_end'=>$p2+1);
# 								push(@{$chimeras->{$i}},\%chimera);
# 							}
# 						}
# 					}
# 					last;
# 				}
# 			}
# 		}
# 	}

# 	return $chimeras;

}

#################################################################################

# Checks if a sequence is chimera from another two from a set of sequences
sub is_chimera {

	my ($test_seq,$seqs,$min_len,$max_ident) = @_;

	# If there are not enough parental sequences to look for chimeras
	if ($#{$seqs}<1){
		return 0;
	}

	my $test_seq_len = length($test_seq);

	# Detects only chimeras of more than a threshold
	# Threshold can be a number of nts or a percentage of sequence length
	# Defines minimum length of a partial match in the chimera (Default=10bp)
	if (!defined($min_len)){
		$min_len = 10; # 10%, 15
	} elsif ($min_len =~ /([\d\.]+)%/) {
		$min_len = int($test_seq_len*$1/100);
	}

	# Skip parental sequences very similar to the tested one
	my %skip_seqs;
	# $max_ident is used to check if parental sequence is disimilar enough
	# Ex.:
	# >Seq1 (false chimera)
	# CGTACAGTGGCTTTACGGCTGTGAGCTCGTGTCTGACGGGAGCATCCATGGTTCCGATCGGTTCGGCTACGATGGTCGGGATTTCATCTCCTTTGAGGGATCCGGGAGATTCGTGGCGGCCGACAGCGCTGCTGAGATCACCAGGAGGC
	# >Seq2 (1 substition in pos 8)
	# CGTACAGCGGCTTTACGGCTGTGAGCTCGTGTCTGACGGGAGCATCCATGGTTCCGATCGGTTCGGCTACGATGGTCGGGATTTCATCTCCTTTGAGGGATCCGGGAGATTCGTGGCGGCCGACAGCGCTGCTGAGATCACCAGGAGGC
	# >Seq3 ()
	# CGTACAGTGGCTTTACGGCTGTGAGCTCCTGTCTGACGGGAGCATCCATGGTTCCGATCGGTTCGGCTACGATGGTCGGGATTTCATCTCCTTTGAGCTGGGATCCGGGAGATTCGTGGCGGCCGACAGCGCTGCTGAGATCACCAGGA
	if (defined($max_ident)){
		# Identity threshold can be a number of nts or a percentage of sequence length
		if ($max_ident =~ /([\d\.]+)%/) {
			$max_ident = int($test_seq_len*$1/100);
		}
		my $aligned_seqs = align_seqs2one($test_seq,"test",$seqs,[0..$#{$seqs}],'needleall');
		if (defined($aligned_seqs)){
			foreach (my $i=0; $i<=$#{$seqs}; $i++){
				# But then true chimeras can be wrongly skipped because artificial indels in the alignment
				# >Seq1 (true chimera)
				# AGCTCAAAGACATCCAGTACATCAGGTCC-ATCTATTACAACAAGCTGGAGTTCATCAGGTTTGACAGCAACGTGGGGGAGTTTGTTGGATACACGGAGCTGGGAGTGAAGAACGCTGAGCGGCTCAACAACGACCCGTCACAGATCGCTGGGTTGAACGCTCAGAGGGAGGCCTACTGCATGAACCATGTTACTGCTTTCTACCCAAACGCTCTGGA
				# >Seq3 (parental sequence)    *  * 2bp diff
				# AGCTCAAAGACATCCAGTACATCAGGTCCTAT-TATTACAACAAGCTGGAGTTCATCAGGTTTGACAGCAACGTGGGGGAGTTTGTTGGATACACGGAGCTGGGAGTGAAGAACGCTGAGCGGCTCAACAACGACCCGTCACAGATCGCTGGGTTGAACGCTCAGAGGGAGGCCTACTGCATGAACCATGTTACTGCTTTCTACCCAAACGCTCTGGA
				# If the test sequence and the parental have the same length, removes all the indels from the alignment
				if ($test_seq_len == length($seqs->[$i])){
	# 			# If the number of alignment gaps is lower than identity threshold,
	# 			if (@{[$aligned_seqs->{$i}[0]=~/\-/g]} <= $max_ident){
					$aligned_seqs->{$i}[0] =~ s/\-//g;
					$aligned_seqs->{$i}[1] =~ s/\-//g;
				}
				my ($ident,$total) = binary_score_nts($aligned_seqs->{$i}[0],$aligned_seqs->{$i}[1]);
				if (defined($ident) && defined($total) && $total-$ident <= $max_ident){
					$skip_seqs{$i} = 1;
				}
			}
		}
	}
	
	# Saves previously calculated sequence differences
	my %matches_right;

	# Loops all the sequences
	for (my $j1=0; $j1<=$#{$seqs}; $j1++){
		my $seq1 = $seqs->[$j1];
		if ($seq1 eq $test_seq) { next; } # Skips identical sequences
		# Compare both sequences to annotate how many nts are in common in both sides
		# $diff_pos is the position of the first substitution (starting from 1...)
		my ($identical1,$diff_pos1) = compare_strings($test_seq,$seq1);
		# If the sequence match the test_seq in its right side, look for another sequence matching the left side
		if ($identical1 == 0 && $diff_pos1-1 >= $min_len){
			# If parental seq is not enough disimilar looks for next parental
			if (defined($skip_seqs{$j1})){
				next;
			}
			for (my $j2=0; $j2<=$#{$seqs}; $j2++){
				my $seq2 = $seqs->[$j2];
				if ($seq2 eq $test_seq || $seq1 eq $seq2) { next; } # Skips identical sequences
				my ($identical2,$diff_pos2);
				if (defined($matches_right{$j2})){
					($identical2,$diff_pos2) = (0,$matches_right{$j2});
				} else {
					($identical2,$diff_pos2) = compare_strings(scalar reverse($test_seq), scalar reverse($seq2));
				}
				if ($identical2 == 0 && $diff_pos2-1 >= $min_len) {
					# If parental seq is not enough disimilar looks for next parental
					if (defined($skip_seqs{$j2})){
						next;
					}
					if ($diff_pos1+$diff_pos2-2>=$test_seq_len){
						return (1,$j1,$j2);
					}
					if (!defined($matches_right{$j2})){
						$matches_right{$j2}=$diff_pos2;
					}
				} elsif (!defined($matches_right{$j2})){
					$matches_right{$j2}=0;
				}
			}
		}
	}

	return 0;

}


#################################################################################

# Align seqs from two FASTA format files
sub align_seqs_from_file {

	my ($query_file,$sbjct_file,$maxhits,$align_options,$keep_align_file) = @_;

	my ($align_data, $align_file);
	
	my ($align_type,$align_params) = ('','');
	if ($align_options =~ /(.+?)\s+(.+)/i){
		$align_type = $1;
		$align_params = $2;
	}
	$align_type =~ s/\s+//g;

	# If not Sbjct file defined, uses GeneBank or UniProt
	if (!defined($sbjct_file) && $align_type eq 'dna'){
		$sbjct_file = $NCBI_NT_DATABASE;
	} elsif (!defined($sbjct_file) && $align_type eq 'prot'){
		$sbjct_file = $NCBI_NR_DATABASE;
	}

	# Annotates desired blastn output format (Ex. 0 = Default pairwise format, 10 = Comma-separated values)
	my $format;
	if ($align_params =~ /-outfmt (\d+)/ && $1 != 11){
		$format = $1;
	}

	# If $align_options = 'dna blastn-short -evalue 0.001'
	if ($align_type eq 'dna' && $align_params =~ /short/i){
		if ($align_params =~ /short\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		} 
		$align_file = execute_blastnshort($query_file,$sbjct_file,$maxhits,$align_params);
		my $blast_data = extract_blast_data($align_file,$maxhits,$format);
		$align_data = extract_align_blast_data($blast_data);

	# If $align_options = 'dna blastn -evalue 0.001'
	} elsif ($align_type eq 'dna' && $align_params =~ /blastn/i){
		if ($align_params =~ /blastn\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		} 
		$align_file = execute_blastn($query_file,$sbjct_file,$maxhits,$align_params);
		my $blast_data = extract_blast_data($align_file,$maxhits,$format);
		$align_data = extract_align_blast_data($blast_data);

	# If $align_options = 'prot blastp-short -evalue 0.001'
	} elsif ($align_type eq 'prot' && $align_params =~ /short/i){
		if ($align_params =~ /short\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		} 
		$align_file = execute_blastpshort($query_file,$sbjct_file,$maxhits,$align_params);
		my $blast_data = extract_blast_data($align_file,$maxhits,$format);
		$align_data = extract_align_blast_data($blast_data);

	# If $align_options = 'prot blastp -evalue 0.001'
	} elsif ($align_type eq 'prot' && $align_params =~ /blastp/i){
		if ($align_params =~ /blastp\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		} 
		$align_file = execute_blastp($query_file,$sbjct_file,$maxhits,$align_params);
		my $blast_data = extract_blast_data($align_file,$maxhits,$format);
		$align_data = extract_align_blast_data($blast_data);

	# If $align_options = 'dna megablast -evalue 0.001'
	} elsif ($align_type eq 'dna' && $align_params =~ /megablast/i){
		if ($align_params =~ /megablast\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		} 
		$align_file = execute_megablast($query_file,$sbjct_file,$maxhits,$align_params);
		my $blast_data = extract_blast_data($align_file,$maxhits,$format);
		$align_data = extract_align_blast_data($blast_data);
	
	# If $align_options = 'dna bwa'
	} elsif ($align_type eq 'dna' && $align_params =~ /bwa/i){
		if ($align_params =~ /bwa\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		} 
		$align_file = execute_bwa($query_file,$sbjct_file,$align_params);
		$align_data = extract_align_bam_data($align_file,$sbjct_file);

	# If $align_options = 'dna bowtie'
	} elsif ($align_type eq 'dna' && $align_params =~ /bowtie/i){
		if ($align_params =~ /bowtie\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		} 
		$align_file = execute_bowtie($query_file,$sbjct_file,$align_params);
		$align_data = extract_align_bam_data($align_file,$sbjct_file);

	# If $align_options = 'dna bowtie2'
	} elsif ($align_type eq 'dna' && $align_params =~ /bowtie2/i){
		if ($align_params =~ /bowtie2\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		} 
		$align_file = execute_bowtie2($query_file,$sbjct_file,$align_params);
		$align_data = extract_align_bam_data($align_file,$sbjct_file);

	# If $align_options = 'dna gassst -r 0 -p 80'
	} elsif ($align_type eq 'dna' && $align_params =~ /gassst/i){
		if ($align_params =~ /gassst\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		} 
		$align_file = execute_gassst($query_file,$sbjct_file,0,$align_params);
# 		my $gassst_file = execute_gassst($query_file,$sbjct_file,$maxhits,$align_params);
		$align_data = extract_align_gassst_data($align_file,$maxhits);

	# If $align_options = 'dna regex -r 0 -p 80'
	} elsif ($align_type eq 'dna' && $align_params =~ /regex/i){
		if ($align_params =~ /regex\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		} 
		$align_file = execute_regex_match($query_file,$sbjct_file,$align_params);
		$align_data = extract_align_regex_data($align_file,$maxhits);

	# If $align_options = 'dna match revcomp'
	} elsif ($align_type eq 'dna' && $align_params =~ /match/i){
		if ($align_params =~ /match\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		}
		my $revcomp = 0;
		if ($align_params =~ /revcomp/i){
				$revcomp = 1;
		}
		my $minlen = 0;
		if ($align_params =~ /minlen (\d+)/i){
				$minlen = $1;
		}

		my ($query_seqs,$query_names) = read_FASTA_file($query_file);
		my ($sbjct_seqs,$sbjct_names) = read_FASTA_file($sbjct_file);
		
		$align_data =  align_seqs($query_names,$query_seqs,$sbjct_names,$sbjct_seqs,$maxhits,$align_options,$keep_align_file)

	}

	if (!defined($keep_align_file)){
		if (defined($align_file)){
			`rm "$align_file"`;
		}
		return $align_data;
	} else {
		return ($align_data, $align_file);
	}


}

#################################################################################

# UNDER DEVELOPMENT
# Maps reads from one single-end FASTQ file or 2 paired-end FASTQ files ($query_files) to references from another FASTA file ($sbjct_file)
sub map_reads_from_files {

	my ($query_files,$sbjct_file,$align_options,$keep_align_file) = @_;

	my ($align_data, $align_file);
	
	my ($align_type,$align_params) = ('','');
	if ($align_options =~ /(.+?)\s+(.+)/i){
		$align_type = $1;
		$align_params = $2;
	}
	$align_type =~ s/\s+//g;

	if ($align_type eq 'bwa'){
		$align_file = execute_bwa($query_files,$sbjct_file,$align_params);
	} elsif ($align_type eq 'bowtie'){
		$align_file = execute_bowtie($query_files,$sbjct_file,$align_params);
	} elsif ($align_type eq 'bowtie2'){
		$align_file = execute_bowtie2($query_files,$sbjct_file,$align_params);
	}
	if (!defined($keep_align_file)){
		if (defined($align_file)){
			`rm "$align_file"`;
		}
		return $align_data;
	} else {
		return ($align_data, $align_file);
	}

}

#################################################################################

# Align protein or DNA sequences
sub align_seqs {

	my ($query_names,$query_seqs,$sbjct_names,$sbjct_seqs,$maxhits,$align_options,$keep_align_file) = @_;

	my ($align_data, $align_file);

	my ($align_type,$align_params) = ('','');
	if ($align_options =~ /(.+?)\s+(.+)/i){
		$align_type = $1;
		$align_params = $2;
	}
	$align_type =~ s/\s+//g;
	
	# Creates FASTA files and runs the alignment software with 'align_seqs_from_file'
	unless ($align_type eq 'dna' && $align_params =~ /match/i){

		my $query_file = create_fasta_file($query_seqs,$query_names);
		my $sbjct_file;
		if (defined($sbjct_seqs) && defined($sbjct_names)){
			$sbjct_file = create_fasta_file($sbjct_seqs,$sbjct_names);
		} else {
			if ($align_type eq 'dna'){
				$sbjct_file = $NCBI_NT_DATABASE;
			} elsif ($align_type eq 'prot'){
				$sbjct_file = $NCBI_NR_DATABASE;
			}
		}
		if (!defined($keep_align_file)){
			$align_data = align_seqs_from_file($query_file,$sbjct_file,$maxhits,$align_options);
		} else {
			($align_data, $align_file) = align_seqs_from_file($query_file,$sbjct_file,$maxhits,$align_options,$keep_align_file);
		}
		`rm $query_file $sbjct_file`;

	# If $align_type = 'dna match revcomp'
	# FASTA files doesn't need to be created, there is no external software
	} else {

		if ($align_params =~ /match\s+(.+)/i){
			$align_params = $1;
		} else {
			$align_params = '';
		}
		my $revcomp = 0;
		if ($align_params =~ /revcomp/i){
				$revcomp = 1;
		}
		my $minlen = 0;
		if ($align_params =~ /minlen (\d+)/i){
				$minlen = $1;
		}
		for (my $i=0; $i<=$#{$query_seqs}; $i++) {
			my @query_seqs = ($query_seqs->[$i]);
			if ($revcomp){
				push(@query_seqs, iupac_reverse_complementary($query_seqs->[$i]));
			}
			for (my $j=0; $j<=$#{$sbjct_seqs}; $j++) {
				my %align_data_;
				for (my $k=0; $k<=$#query_seqs; $k++) {
					if ($query_seqs[$k] =~ /$sbjct_seqs->[$j]/i){
						if ($minlen && length($&)<$minlen) { next; }
						$align_data_{'NAME'} = $sbjct_names->[$j];
	# 					$align_data_{'ALIGNED'} = length($sbjct_seqs->[$j]);
	# 					$align_data_{'SIMIL'} = length($sbjct_seqs->[$j]);
						$align_data_{'ALIGN'} = $&."\n".$&; # Text that matchs
						$align_data_{'IDENT'} = $align_data_{'ALIGNED'} = length($&);
						if ($k == 0) {
							$align_data_{'COLS'} = join(",", $-[0]+1 .. $+[0])."\n".join(",", 1 .. length($sbjct_seqs->[$j]));
						} else {
							$align_data_{'COLS'} = join(",", length($query_seqs->[$i])-$+[0]+1 .. length($query_seqs->[$i])-$-[0])."\n".join(",", reverse 1 .. length($sbjct_seqs->[$j]) );
						}
# 						print "\n".$query_names->[$i]." ".$sbjct_names->[$j]."\n";
# 						print $query_seqs->[$i]."\n".$sbjct_seqs->[$j]."\n";
						push (@{$align_data->{$query_names->[$i]}},\%align_data_);
					} elsif ($sbjct_seqs->[$j] =~ /$query_seqs[$k]/i){ # NOT CHECKED
						if ($minlen && length($&)<$minlen) { next; }
						$align_data_{'NAME'} = $sbjct_names->[$j];
						$align_data_{'ALIGN'} = $&."\n".$&; # Text that matchs
						$align_data_{'IDENT'} = $align_data_{'ALIGNED'} = length($&);
						if ($k == 0) {
							$align_data_{'COLS'} = join(",", 1 .. length($query_seqs[$k]) ) . "\n" . join(",", $-[0]+1 .. $+[0]);
						} else {
							$align_data_{'COLS'} = join(",", reverse 1 .. length($query_seqs[$k]) ) . "\n" . join(",", length($sbjct_seqs->[$j])-$+[0]+1 .. length($sbjct_seqs->[$j])-$-[0] );
						}
						push (@{$align_data->{$query_names->[$i]}},\%align_data_);
					}
				}
			}
		}
	}

	if (!defined($keep_align_file)){
		return $align_data;
	} else {
		return ($align_data, $align_file);
	}

}

#################################################################################

# Use threads to align seqs quickly, doing alignments in different threats
# Is useless when less than 200-500 alignments per query

# THREADS MODULE HAVE A MEMORY LEAK IN PERL 5.10 WHEN MANY THREADS ARE CREATED, THEY ARE NOT COMPLETELY DELETED FROM MEMORY WHEN JOINED

sub align_seqs_with_threads {

	my ($query_names,$query_seqs,$sbjct_names,$sbjct_seqs,$maxhits,$align_type,$threads_limit,$keep_align_file) = @_;

	if (!defined($threads_limit)){
		$threads_limit = 4;
	}

	my $align_data = {};
	
	my $align_file;

	my (@threads, $thread_seq_data);

	my $thread_seq_data_number = abs(int(scalar @{$query_names} / $threads_limit))+1;

	my ($thread_query_names,$thread_query_seqs);

	for (my $count_query=0; $count_query<=$#{$query_seqs}; $count_query++){

		push(@{$thread_query_names}, $query_names->[$count_query]);
		push(@{$thread_query_seqs}, $query_seqs->[$count_query]);

		if (($count_query+1) % $thread_seq_data_number == 0  || $count_query == $#{$query_seqs}) {

			push(@threads, threads->create(\&align_seqs,$thread_query_names,$thread_query_seqs,$sbjct_names,$sbjct_seqs,$maxhits,$align_type,$keep_align_file));

# 			# For debugging:
# 			push(@threads, [align_seqs($thread_query_names,$thread_query_seqs,$sbjct_names,$sbjct_seqs,$maxhits,$align_type,$keep_align_file)]);
# 
			undef($thread_query_names);
			undef($thread_query_seqs);

		}

		
		# If maximum number of threads is reached or last sbjct of a query is processed
		if (scalar @threads >= $threads_limit  || $count_query == $#{$query_seqs}){
			my $check_threads = 1;
			while ($check_threads){
				for (my $i=0; $i<=$#threads; $i++){
					unless ($threads[$i]->is_running()){
						my ($align_data_, $align_file_);
						if (!defined($keep_align_file)){
							($align_data_) = $threads[$i]->join;
						} else {
							($align_data_,$align_file_) = $threads[$i]->join;
						}
						if (defined($align_data_)) {
							$align_data = { %$align_data, %$align_data_ };
						}
						if (!defined($align_file)){
							$align_file = $align_file_;
						} else {
							`cat "$align_file_" >> "$align_file"`;
						}
						undef($threads[$i]);
						splice(@threads,$i,1);
						$i = $i - 1;
						unless ($count_query == $#{$query_seqs} && @threads){
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
# 				my ($align_data_, $align_file_);
# 				if (!defined($keep_align_file)){
# 					($align_data_)  = @{$threads[$i]};
# 				} else {
# 					($align_data_,$align_file_)  = @{$threads[$i]};
# 				}
# 				if (defined($align_data_)) {
# 					$align_data = { %$align_data, %$align_data_ };
# 				}
# 				if (!defined($align_file)){
# 					$align_file = $align_file_;
# 				} else {
# 					`cat "$align_file_" >> "$align_file"`;
# 				}
# 				delete $threads[$i];
# 			}

		}


	}

	if (!defined($keep_align_file)){
		return $align_data;
	} else {
		return ($align_data, $align_file);
	}

}

#################################################################################

# Multiple align sequences
sub multiple_align_seqs {

	my ($seqs,$names,$align_type) = @_;

	if (!defined($align_type)){
		$align_type='mafft';
	}
	my $align_options = '';
	if ($align_type =~ /(.+?)\s+(.+)/i){
		$align_type = $1;
		$align_options = $2;
	}
	$align_type =~ s/\s+//g;

	# If $align_type = 'mafft'
	if ($align_type =~ /mafft/){
		my $query_file = create_fasta_file($seqs,$names);
		my $alignment_file = execute_mafft($query_file,$align_options);
		my ($sequences,$headers) = read_fasta_file($alignment_file);
		`rm "$query_file" "$alignment_file"`;
		return ($sequences,$headers);
	}
}

#################################################################################

# Global align 2 sequences
sub align_2seqs {

	my ($seq1,$seq2,$align_type) = @_;

	if (!defined($align_type)){
		$align_type='global';
	}
	my $align_options = '';
	if ($align_type =~ /(.+?)\s+(.+)/i){
		$align_type = $1;
		$align_options = $2;
	}
	$align_type =~ s/\s+//g;

	if ($align_type =~ /needle|wunsch|global|nw/){
# 		print "$NEEDLEEXE $align_options -auto -aformat3 fasta -stdout 'echo $seq1 |' 'echo $seq2 |' \n"; exit;
		open(ALIGN,"$NEEDLEEXE $align_options -auto -aformat3 fasta -stdout 'echo $seq1 |' 'echo $seq2 |' |") || die "# $0 : cannot run $NEEDLEEXE\n";
		my ($sequences,$headers) = read_fasta([<ALIGN>]);
		close(ALIGN);
		return ($sequences->[0],$sequences->[1]);
	} elsif ($align_type =~ /needleman-wunsch/){
		open(ALIGN,"$NEEDLEMANWUNSCHEXE $align_options $seq1 $seq2 |") || die "# $0 : cannot run $NEEDLEMANWUNSCHEXE\n";
		chomp(my @sequences = <ALIGN>);
		close(ALIGN);
		return ($sequences[0],$sequences[1]);
	} elsif ($align_type =~ /local|sw|smith|waterman/){
		open(ALIGN,"$SMITHWATERMANEXE --maxhits 1 $align_options $seq1 $seq2 |") || die "# $0 : cannot run $SMITHWATERMANEXE\n";
		my @sequences;
		while (<ALIGN>){
			if (/^  (.*)  \[pos: \d+; len: \d+\]$/i){
				push(@sequences, $1);
			}
		}
		close(ALIGN);
		return ($sequences[0],$sequences[1]);
	} elsif ($align_type =~ /mafft/){
		my $query_file = create_fasta_file([$seq1,$seq2],['1','2']);
		open(ALIGN,"$MAFFTEXE --quiet $align_options '$query_file' |")|| die "# $0 : cannot run $MAFFTEXE\n";
		`rm "$query_file"`;
		my ($sequences,$headers) = read_fasta([<ALIGN>]);
		close(ALIGN);
# 		my $alignment_file = execute_mafft($query_file,$align_options);
# 		my ($sequences,$headers) = read_fasta_file($alignment_file);
# 		`rm "$query_file" "$alignment_file"`;
		return ($sequences->[0],$sequences->[1]);
	} elsif ($align_type =~ /fogsaa/){
		open(ALIGN,"$FOGSAAEXE $seq1 $seq2 1 0 1 -1 -2 |")|| die "# $0 : cannot run $FOGSAAEXE\n";
		chomp(my @sequences = <ALIGN>);
# 		if (!@sequences) {
# 			print "\n\n$FOGSAAEXE $seq1 $seq2 1 0 1 -1 -2\n\n";exit;
# 		}
		close(ALIGN);
		return ($sequences[0],$sequences[1]);
	}

}

#################################################################################

# Globally align sbjct sequences, one by one, to a query sequence
sub align_seqs2one {

	my ($query_seq,$query_name,$sbjct_seqs,$sbjct_names,$align_type) = @_;

	my $results;

	if (!defined($align_type)){
		$align_type='global';
	}
	my $align_options = '';
	if ($align_type =~ /(.+?)\s+(.+)/i){
		$align_type = $1;
		$align_options = $2;
	}
	$align_type =~ s/\s+//g;

	my $outfile = sprintf("/tmp/%s.align",random_file_name());
	unlink($outfile);

	if ($align_type =~ /global|needleall/){
		my $query_file = create_fasta_file([$query_seq], [$query_name]);
		my $sbjct_file = create_fasta_file($sbjct_seqs, $sbjct_names);
# 		print  "$NEEDLEALLEXE -auto -aformat3 fasta -asequence '$query_file' -bsequence '$sbjct_file' -outfile '$outfile' -errfile '$outfile.error'\n";exit;
		system("$NEEDLEALLEXE -auto -aformat3 fasta -asequence '$query_file' -bsequence '$sbjct_file' -outfile '$outfile' -errfile '$outfile.error'");
		my ($aligned_seqs,$aligned_headers) = read_fasta_file($outfile);
		`rm '$query_file' '$sbjct_file' '$outfile' '$outfile.error'`;
# 		# DOESN'T WORK WHEN THERE ARE MANY SEQUENCES:
# 		# "/bin/sh: Argument list too long"
# 		my $needleall_query = ">$query_name\n$query_seq";
# 		my $needleall_sbjct;
# 		map $needleall_sbjct .= ">".$sbjct_names->[$_]."\n".$sbjct_seqs->[$_]."\n" , 0 .. $#{$sbjct_names};
# 		open(ALIGN, "$NEEDLEALLEXE -auto -aformat3 fasta -stdout 'echo $needleall_query |' 'echo $needleall_sbjct |' |") || die "# $0 : cannot run $NEEDLEALLEXE\n";
# 		my ($aligned_seqs,$aligned_headers) = read_fasta([<ALIGN>]);
# 		close(ALIGN);
		for (my $i=0; $i<=$#{$aligned_headers}; $i+=2) {
			if ($aligned_headers->[$i] eq $query_name){
				$results->{$aligned_headers->[$i+1]}=[$aligned_seqs->[$i], $aligned_seqs->[$i+1]];
			}
		}
	} elsif ($align_type =~ /needleman-wunsch/){
		my $query_file = create_fasta_file([($query_seq) x scalar @$sbjct_seqs]);
		my $sbjct_file = create_fasta_file($sbjct_seqs);
# 		print "$NEEDLEMANWUNSCHEXE $align_options --files '$query_file' '$sbjct_file' > '$outfile'\n";exit;
		open(ALIGN, "$NEEDLEMANWUNSCHEXE $align_options --files '$query_file' '$sbjct_file' |") || die "# $0 : cannot run $NEEDLEMANWUNSCHEXE\n";
		my ($i,$count_align,$seq1,$seq2) = (0,0,'','');
		while(<ALIGN>){
			chomp;
			if (/^$/) {
				$count_align=0;
			} elsif ($count_align==0) {
				$seq1 = $_;
				$count_align++;
			} elsif ($count_align==1) {
				$seq2 = $_;
				$results->{$sbjct_names->[$i]}=[$seq1, $seq2];
				$i++;
				$count_align=0;
			}
		}
		close(ALIGN);
		`rm '$query_file' '$sbjct_file'`;
		if (scalar @$sbjct_seqs != scalar keys %$results) {
			die "# $0 : # align_seqs2one: different number of alignments than sequences.\n";
		}
	} 
	return $results;
# 	elsif ($align_type =~ /local|sw|smith|waterman/){
# 		open(ALIGN,"$SMITHWATERMANEXE --maxhits 1 $align_options $seq1 $seq2 |")|| die "# $0 : cannot run $SMITHWATERMANEXE\n";
# 		my @sequences;
# 		while (<ALIGN>){
# 			if (/^  (.*)  \[pos: \d+; len: \d+\]$/i){
# 				push(@sequences, $1);
# 			}
# 		}
# 		return ($sequences[0],$sequences[1]);
# 		close(ALIGN);
# 	} elsif ($align_type =~ /mafft/){
# 		my $query_file = create_fasta_file([$seq1,$seq2],['1','2']);
# 		open(ALIGN,"$MAFFTEXE --quiet $align_options '$query_file' |")|| die "# $0 : cannot run $MAFFTEXE\n";
# 		`rm "$query_file"`;
# 		my ($sequences,$headers) = read_fasta([<ALIGN>]);
# 		close(ALIGN);
# # 		my $alignment_file = execute_mafft($query_file,$align_options);
# # 		my ($sequences,$headers) = read_fasta_file($alignment_file);
# # 		`rm "$query_file" "$alignment_file"`;
# 		return ($sequences->[0],$sequences->[1]);
# 	}
}

#################################################################################

# Globally align sbjct sequences to query sequences, one by one
sub align_seqs2seqs {

	my ($query_seqs,$query_names,$sbjct_seqs,$sbjct_names,$align_type) = @_;

	my $results;

	if (!defined($align_type)){
		$align_type='global';
	}

	my $outfile = sprintf("/tmp/%s.align",random_file_name());
	unlink($outfile);

	if ($align_type =~ /global|needleall/){
		my $query_file = create_fasta_file($query_seqs, $query_names);
		my $sbjct_file = create_fasta_file($sbjct_seqs, $sbjct_names);
# 		print  "$NEEDLEALLEXE -auto -aformat3 score -stdout -asequence '$query_file' -bsequence '$sbjct_file' -errfile '$outfile.error'\n";exit;
		open(ALIGN,"$NEEDLEALLEXE -auto -aformat3 score -stdout -asequence '$query_file' -bsequence '$sbjct_file' -errfile '$outfile.error' |");
		while (<ALIGN>){
			if (/(.+)\s(.+)\s(.+)\s\((.+)\)/){
				$results->{$2}{$1} = $4;
			}
		}
		close ALIGN;
		`rm '$query_file' '$sbjct_file' '$outfile.error'`;
	}

	return $results;
	
}

#################################################################################

# Cluster sequences
sub cluster_seqs {

	my ($query_seqs,$query_names,$cluster_type,$clustering_threshold) = @_;

	if (!defined($cluster_type)){
		$cluster_type='cd-hit';
	}
	my $cluster_options = '';
	if ($cluster_type =~ /(.+?)\s+(.+)/i){
		$cluster_type = $1;
		$cluster_options = $2;
	}
	$cluster_type =~ s/\s+//g;

	# If $cluster_type = 'cd-hit-454'
	if ($cluster_type =~ /cd-?hit-454/i){
		my $query_file = create_fasta_file($query_seqs,$query_names);
		my ($unique_seqs_file, $clusters_file) = execute_cdhit_454($query_file,$clustering_threshold,$cluster_options);
		my $clusters = extract_cdhit_data($clusters_file,$query_file);
		unlink($unique_seqs_file, $clusters_file, $query_file);
		return $clusters;
	}

	# If $cluster_type = 'cd-hit-est'
	if ($cluster_type =~ /cd-?hit-est/i){
		my $query_file = create_fasta_file($query_seqs,$query_names);
		my ($unique_seqs_file, $clusters_file) = execute_cdhit_est($query_file,$clustering_threshold,$cluster_options);
		my $clusters = extract_cdhit_data($clusters_file,$query_file);
		unlink($unique_seqs_file, $clusters_file, $query_file);
		return $clusters;
	}

	# If $cluster_type = 'cd-hit'
	if ($cluster_type =~ /cd-?hit/i){
		my $query_file = create_fasta_file($query_seqs,$query_names);
		my ($unique_seqs_file, $clusters_file) = execute_cdhit($query_file,$clustering_threshold,$cluster_options);
		my $clusters = extract_cdhit_data($clusters_file,$query_file);
		unlink($unique_seqs_file, $clusters_file, $query_file);
		#print $unique_seqs_file, $clusters_file, "\n";
		return $clusters;
	}
	
	# If $cluster_type = 'manual'
	if ($cluster_type =~ /manual/i){
		my ($aligned_seqs,$aligned_names) = multiple_align_seqs($query_seqs,$query_names,'mafft');
		my $clusters;
		my @clustered_seqs;
		my $cluster_count = 0;
		for (my $i=0; $i<=$#{$aligned_names}; $i++) {
			if (in_array(\@clustered_seqs,$i)) { next; }
			$cluster_count++;
			my $aligned_seq1 = $aligned_seqs->[$i];
			push(@{$clusters->{$cluster_count}}, { 'name'=>$aligned_names->[$i], 'identity'=>100 });
			push(@clustered_seqs,$i);
			for (my $j=0; $j<=$#{$aligned_names}; $j++) {
				if ($i==$j || in_array(\@clustered_seqs,$j)) { next; }
				my $aligned_seq2 = $aligned_seqs->[$j];
				my ($identity,$total) = binary_score_nts($aligned_seq1,$aligned_seq2);
				$identity = sprintf("%.2f",$identity/$total*100);
				# Skips seqs with lower identity than threshold
				if ($identity < $clustering_threshold){
					next;
				}
				push(@{$clusters->{$cluster_count}}, { 'name'=>$aligned_names->[$j], 'identity'=>$identity });
				push(@clustered_seqs,$j);
			}
		}
		return $clusters;
	}

	# If $cluster_type = 'manual'
	if ($cluster_type =~ /upgma/i){
		return upgma(align_seqs2seqs($query_seqs,$query_names,$query_seqs,$query_names))
	}


}

#################################################################################

# Returns a hash with unique sequences and their duplicates (headers)
sub cluster_identical_seqs {

	my ($seqs,$headers) = @_;

	my $seq_clusters;

	my $count_unique = 0;
	for (my $i=0; $i<=$#{$seqs}; $i++){
		# if (($i+1) % 10000 == 0) {
		# print "$i $count_unique\n";
		# }
		if (!defined($seq_clusters->{$seqs->[$i]})){
			$seq_clusters->{$seqs->[$i]} = [$headers->[$i]];
			$count_unique++;
		} else {
			push(@{$seq_clusters->{$seqs->[$i]}},$headers->[$i]);
		}
	}
	
	return $seq_clusters;

}

#################################################################################

# Returns a hash with unique sequences and their duplicates (headers)
sub cluster_identical_seqs_with_threads {

	my ($seqs,$headers,$threads) = @_;

	my $seq_clusters;

	my @threads;
	my $total_seqs = scalar @$seqs;
	my $part_seqs = int($total_seqs/$threads);
	while (@$seqs){
		my @seqs_ = splice(@$seqs,0,$part_seqs);
		my @headers_ = splice(@$headers,0,$part_seqs);
		push(@threads, threads->create(\&cluster_identical_seqs,\@seqs_,\@headers_));
# 		# For debugging:
# 		push(@threads, [cluster_identical_seqs(\@seqs_,\@headers_)]);
	}
	while (){
		for (my $i=0; $i<=$#threads; $i++){
			unless ($threads[$i]->is_running()){
				my ($seq_clusters_) = $threads[$i]->join;
				if (%$seq_clusters_){
					foreach my $unique_seq (keys %$seq_clusters_) {
						push(@{$seq_clusters->{$unique_seq}}, @{$seq_clusters_->{$unique_seq}});
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
# 		my ($seq_clusters_) = @{$threads[$i]};
# 		if (%$seq_clusters_){
# 			foreach my $unique_seq (keys %$seq_clusters_) {
# 				push(@{$seq_clusters->{$unique_seq}}, @{$seq_clusters_->{$unique_seq}});
# 			}
# 		}
# 		delete $threads[$i];
# 	}

	return $seq_clusters;

}

#################################################################################

# Returns a hash with clusters of similar sequences (headers)
sub cluster_illumina_seqs {

	my ($seqs,$headers,$substitutions,$depths) = @_;

	my $seq_clusters;

	# First clusters identical ones
	my $unique_seq_clusters = cluster_identical_seqs($seqs,$headers);
	if (!defined($substitutions) || $substitutions == 0){
		return $unique_seq_clusters;
	}

	# Generates a hash assigning depths to the sequence headers
	# Should not be duplicated headers
	my $header_to_depth;
	if (defined($depths)){
		if (ref($depths) eq "ARRAY") {
			map $header_to_depth->{$headers->[$_]} += $depths->[$_], 0..$#{$headers};
		} elsif (ref($depths) eq "HASH") {
			$header_to_depth = $depths;
		}
	}

	# Orders sequences by identical sequences cluster sizes
	my $unique_seq_clusters_depths;
	foreach my $unique_seq (keys %$unique_seq_clusters){
		if (defined($depths)){
			map $unique_seq_clusters_depths->{$unique_seq} += $header_to_depth->{$_} , @{$unique_seq_clusters->{$unique_seq}};
		} else {
			$unique_seq_clusters_depths->{$unique_seq} = scalar @{$unique_seq_clusters->{$unique_seq}};
		}
	}
	my @unique_seqs =  sort { $unique_seq_clusters_depths->{$b} <=> $unique_seq_clusters_depths->{$a} } keys %$unique_seq_clusters_depths;
	# my @unique_seqs =  sort { scalar @{$unique_seq_clusters->{$b}} <=> scalar @{$unique_seq_clusters->{$a}} } keys %$unique_seq_clusters;

	# Second clusters similar ones
	while (my $unique_seq = shift @unique_seqs){
		push(@{$seq_clusters->{$unique_seq}},@{$unique_seq_clusters->{$unique_seq}});
		# my $last_cluster_seq;
		for (my $i=0; $i<=$#unique_seqs; $i++){
			my $diff = count_errors($unique_seq,$unique_seqs[$i],$substitutions+1);
# 			my $diff = count_errors($unique_seq,$unique_seqs[$i],$substitutions+1);
			if (defined($diff) && $diff <= $substitutions) {
				push(@{$seq_clusters->{$unique_seq}},@{$unique_seq_clusters->{$unique_seqs[$i]}});
				# $last_cluster_seq = $unique_seqs[$i];
				splice(@unique_seqs,$i,1);
			}
		}
# 		# If the cluster is composed by only 2 reads, we will take the consensus sequence
# 		if (defined($last_cluster_seq) && $#{$seq_clusters->{$unique_seq}} == 1){
# 			my $consensus_seq = consensus_sequence([$unique_seq,$last_cluster_seq],undef,'dna');
# 			$seq_clusters->{$consensus_seq} = $seq_clusters->{$unique_seq};
# 			delete($seq_clusters->{$unique_seq});
# 		}
	}
# my $kk = 0;
# foreach my $cluster (keys %$unique_seq_clusters) {
# 	foreach my $seq (@{$seq_clusters->{$cluster}}){
# 		$kk += $header_to_depth->{$seq};
# 	}
# }
# print "$kk\n"; exit;
	return $seq_clusters;
# 	my @clustered_seqs =  sort { scalar @{$seq_clusters->{$b}} <=> scalar @{$seq_clusters->{$a}} } keys %$seq_clusters;
# 	my @clustered_headers;
# 	foreach my $seq (@clustered_seqs){
# # 		push(@clustered_headers, join(' | ', @{$seq_clusters->{$seq}}));
# # 		my $md5 = generate_md5($seq);
# # 		my $len = length($seq);
# 		
# 		push(@clustered_headers, join(' | ', @{$seq_clusters->{$seq}}));
# 	}
# 	return (\@clustered_seqs,\@clustered_headers);

}

#################################################################################

# Returns a hash with clusters of similar sequences (headers)
sub cluster_umi_seqs {

	my ($seqs,$headers,$substitutions) = @_;

	my $seq_clusters;

	# First clusters identical ones
	my $unique_seq_clusters = cluster_identical_seqs($seqs,$headers);
	if (!defined($substitutions) || $substitutions == 0){
		return $unique_seq_clusters;
	}

	# Second clusters all sequences with the most abundant one
	# And discards sequences with higher number of errors than $substitutions
	my @unique_seqs =  sort { scalar @{$unique_seq_clusters->{$b}} <=> scalar @{$unique_seq_clusters->{$a}} } keys %$unique_seq_clusters;
	my $unique_seq = shift @unique_seqs;
	push(@{$seq_clusters->{$unique_seq}},@{$unique_seq_clusters->{$unique_seq}});
	for (my $i=0; $i<=$#unique_seqs; $i++){
		my $diff = count_errors($unique_seq,$unique_seqs[$i],$substitutions+1);
		if (defined($diff) && $diff <= $substitutions) {
			push(@{$seq_clusters->{$unique_seq}},@{$unique_seq_clusters->{$unique_seqs[$i]}});
		}
	}
	return $seq_clusters;

}

#################################################################################

# Runs BLASTP against two files with FASTA sequences
sub execute_blastp {

	my ($query_file,$sbjct_file,$maxhits,$blast_options,$outfile) = @_;

	if($maxhits){$maxhits = "-num_descriptions $maxhits -num_alignments $maxhits";}else{$maxhits='';} # "-num_descriptions $maxhits -num_alignments $maxhits"

	if (!defined($blast_options)){$blast_options='';}

	my $random_name = random_file_name();

	if (!defined($outfile)){
		$outfile="/tmp/$random_name.blast";
	}
	unlink($outfile);

	my $blastdb_file;
	if ($sbjct_file eq $NCBI_NR_DATABASE) {
		$blastdb_file = $NCBI_NR_DATABASE;
	} elsif (-e "$sbjct_file.phr" && -e "$sbjct_file.pin" && -e "$sbjct_file.psq") {
		$blastdb_file = $sbjct_file;
	} elsif (is_compressed($sbjct_file)) {
		$blastdb_file = "/tmp/align_seqs_$random_name.db";
		`zcat '$sbjct_file' | $MAKEBLASTDBEXE -in - -title 'BLASTDB' -dbtype prot -out '$blastdb_file' 2>&1`;
	} else {
		$blastdb_file = "/tmp/align_seqs_$random_name.db";
		# print "\n$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype prot -out '$blastdb_file'\n"; exit;
# 		system("$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype prot -out '$blastdb_file' > /dev/null");
		`$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype prot -out '$blastdb_file' 2>&1`;
	}
	
	if (is_compressed($query_file)) {
		`zcat '$query_file' | $BLASTPEXE $blast_options $maxhits -query - -db '$blastdb_file' -out '$outfile' 2>&1`;
	} else{ 
	# 	print "\n$BLASTPEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile'\n";exit;
	# 	system("$BLASTPEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile' > /dev/null");
		`$BLASTPEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile' 2>&1`;
	}
	if ($blastdb_file =~ /\/tmp\/align_seqs_/) {
		`rm $blastdb_file*`;
	}

	return $outfile;

}

#################################################################################

# Runs BLASTP-SHORT against two files with FASTA sequences
sub execute_blastpshort {

	my ($query_file,$sbjct_file,$maxhits,$blast_options,$outfile) = @_;

	if($maxhits){$maxhits = "-num_descriptions $maxhits -num_alignments $maxhits";}else{$maxhits='';} # "-num_descriptions $maxhits -num_alignments $maxhits"

	if (!defined($blast_options)){$blast_options='';}

	my $random_name = random_file_name();

	if (!defined($outfile)){
		$outfile="/tmp/$random_name.blast";
	}
	unlink($outfile);

	my $blastdb_file;
	if ($sbjct_file eq $NCBI_NR_DATABASE) {
		$blastdb_file = $NCBI_NR_DATABASE;
	} elsif (-e "$sbjct_file.phr" && -e "$sbjct_file.pin" && -e "$sbjct_file.psq") {
		$blastdb_file = $sbjct_file;
	} elsif (is_compressed($sbjct_file)) {
		$blastdb_file = "/tmp/align_seqs_$random_name.db";
		`zcat '$sbjct_file' | $MAKEBLASTDBEXE -in - -title 'BLASTDB' -dbtype prot -out '$blastdb_file' 2>&1`;
	} else {
		$blastdb_file = "/tmp/align_seqs_$random_name.db";
# 		print "\n$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype prot -out '$blastdb_file'\n"; exit;
		`$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype prot -out '$blastdb_file' 2>&1`;
	}
	
	if (is_compressed($query_file)) {
		`zcat '$query_file' | $BLASTPSHORTEXE $blast_options $maxhits -query - -db '$blastdb_file' -out '$outfile' 2>&1`;
	} else {
	# 	print "\n$BLASTPSHORTEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile'\n";exit;
		`$BLASTPSHORTEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile' 2>&1`;
	}

	if ($blastdb_file =~ /\/tmp\/align_seqs_/) {
		`rm $blastdb_file*`;
	}

	return $outfile;

}

#################################################################################

# Runs BLASTN against two files with FASTA sequences
sub execute_blastn {

	my ($query_file,$sbjct_file,$maxhits,$blast_options,$outfile) = @_;

	if($maxhits){$maxhits = "-num_descriptions $maxhits -num_alignments $maxhits";}else{$maxhits='';} # "-num_descriptions $maxhits -num_alignments $maxhits"

	if (!defined($blast_options)){$blast_options='';}

	my $random_name = random_file_name();

	if (!defined($outfile)){
		$outfile="/tmp/$random_name.blast";
	}
	unlink($outfile);

	my $blastdb_file;
	if ($sbjct_file eq $NCBI_NT_DATABASE) {
		$blastdb_file = $NCBI_NT_DATABASE;
	} elsif (-e "$sbjct_file.nhr" && -e "$sbjct_file.nin" && -e "$sbjct_file.nsq") {
		$blastdb_file = $sbjct_file;
	} elsif (is_compressed($sbjct_file)) {
		$blastdb_file = "/tmp/align_seqs_$random_name.db";
		#print "\nzcat '$sbjct_file' | $MAKEBLASTDBEXE -in - -title 'BLASTDB' -dbtype nucl -out '$blastdb_file'\n";
		`zcat '$sbjct_file' | $MAKEBLASTDBEXE -in - -title 'BLASTDB' -dbtype nucl -out '$blastdb_file' 2>&1`;
	} else {
		$blastdb_file = "/tmp/align_seqs_$random_name.db";
		# print "\n$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype nucl -out '$blastdb_file'\n"; exit;
# 		system("$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype nucl -out '$blastdb_file' > /dev/null");
		`$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype nucl -out '$blastdb_file' 2>&1`;
	}
	
	if (is_compressed($query_file)) {
# 		print "\nzcat '$query_file' | $BLASTNEXE $blast_options $maxhits -query - -db '$blastdb_file' -out '$outfile'\n";exit;
		`zcat '$query_file' | $BLASTNEXE $blast_options $maxhits -query - -db '$blastdb_file' -out '$outfile' 2>&1`;
	} else{ 
# 		print "\n$BLASTNEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile'\n";exit;
	# 	system("$BLASTNEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile' > /dev/null");
		`$BLASTNEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile' 2>&1`;
	}

	if ($blastdb_file =~ /\/tmp\/align_seqs_/) {
		`rm $blastdb_file*`;
	}

	return $outfile;

}

#################################################################################

# Runs BLASTN-SHORT against two files with FASTA sequences
sub execute_blastnshort {

	my ($query_file,$sbjct_file,$maxhits,$blast_options,$outfile) = @_;

	if($maxhits){$maxhits = "-num_descriptions $maxhits -num_alignments $maxhits";}else{$maxhits='';} # "-num_descriptions $maxhits -num_alignments $maxhits"

	if (!defined($blast_options)){$blast_options='';}

	my $random_name = random_file_name();

	if (!defined($outfile)){
		$outfile="/tmp/$random_name.blast";
	}
	unlink($outfile);

	my $blastdb_file;
	if ($sbjct_file eq $NCBI_NT_DATABASE) {
		$blastdb_file = $NCBI_NT_DATABASE;
	} elsif (-e "$sbjct_file.nhr" && -e "$sbjct_file.nin" && -e "$sbjct_file.nsq") {
		$blastdb_file = $sbjct_file;
	} elsif (is_compressed($sbjct_file)) {
		$blastdb_file = "/tmp/align_seqs_$random_name.db";
		`zcat '$sbjct_file' | $MAKEBLASTDBEXE -in - -title 'BLASTDB' -dbtype nucl -out '$blastdb_file' 2>&1`;
	} else {
		$blastdb_file = "/tmp/align_seqs_$random_name.db";
# 		print "\n$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype nucl -out '$blastdb_file'\n"; exit;
		`$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype nucl -out '$blastdb_file' 2>&1`;
	}
	
	if (is_compressed($query_file)) {
		`zcat '$query_file' | $BLASTNSHORTEXE $blast_options $maxhits -query - -db '$blastdb_file' -out '$outfile' 2>&1`;
	} else {
	# 	print "\n$BLASTNSHORTEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile'\n";exit;
		`$BLASTNSHORTEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile' 2>&1`;
	}

	if ($blastdb_file =~ /\/tmp\/align_seqs_/) {
		`rm $blastdb_file*`;
	}

	return $outfile;

}

#################################################################################

# Runs MEGABLAST against two files with FASTA sequences
sub execute_megablast {

	my ($query_file,$sbjct_file,$maxhits,$blast_options,$outfile) = @_;

	if($maxhits){$maxhits = "-num_descriptions $maxhits -num_alignments $maxhits";}else{$maxhits='';} # "-num_descriptions $maxhits -num_alignments $maxhits"

	if (!defined($blast_options)){$blast_options='';}

	my $random_name = random_file_name();

	if (!defined($outfile)){
		$outfile="/tmp/$random_name.blast";
	}
	unlink($outfile);

	my $blastdb_file;
	if ($sbjct_file eq $NCBI_NT_DATABASE) {
		$blastdb_file = $NCBI_NT_DATABASE;
	} elsif (-e "$sbjct_file.nhr" && -e "$sbjct_file.nin" && -e "$sbjct_file.nsq") {
		$blastdb_file = $sbjct_file;
	} elsif (is_compressed($sbjct_file)) {
		$blastdb_file = "/tmp/align_seqs_$random_name.db";
		`zcat '$sbjct_file' | $MAKEBLASTDBEXE -in - -title 'BLASTDB' -dbtype nucl -out '$blastdb_file' 2>&1`;
	} else {
		$blastdb_file = "/tmp/align_seqs_$random_name.db";
# 		print "\n$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype nucl -out '$blastdb_file'\n"; exit;
		`$MAKEBLASTDBEXE -in '$sbjct_file' -dbtype nucl -out '$blastdb_file' 2>&1`;
	}
	
	if (is_compressed($query_file)) {
		`zcat '$query_file' | $MEGABLASTEXE $blast_options $maxhits -query - -db '$blastdb_file' -out '$outfile' 2>&1`;	
	} else {
	# 	print "\$MEGABLASTEXE $blast_options $maxhits -query '$query_file' -db '/tmp/align_seqs_$random_name.db' -out '$outfile'\n";exit;
		`$MEGABLASTEXE $blast_options $maxhits -query '$query_file' -db '$blastdb_file' -out '$outfile' 2>&1`;	
	}

	if ($blastdb_file =~ /\/tmp\/align_seqs_/) {
		`rm $blastdb_file*`;
	}

	return $outfile;

}

#################################################################################

# Runs GASSST (primer global alignment) against two files with FASTA sequences (Sbjct: short seqs)
sub execute_gassst {

	my ($query_file,$sbjct_file,$maxhits,$gassst_options,$outfile) = @_;

	if($maxhits){$maxhits = "-h $maxhits";}else{$maxhits='-h 0';} # "-num_descriptions $maxhits -num_alignments $maxhits"

	if (!defined($gassst_options)){$gassst_options='-p 90';}

	if (!defined($outfile)){
		my $random_name = random_file_name();
		$outfile="/tmp/$random_name.gassst";
	}
	if (-e $outfile){
		`rm "$outfile"`;
	}

# 	print "\n$GASSSTEXE $gassst_options $maxhits -i '$sbjct_file' -d '$query_file' -o '$outfile'\n"; exit;
	`$GASSSTEXE $gassst_options $maxhits -i '$sbjct_file' -d '$query_file' -o '$outfile' 2>&1`;
	
	return $outfile;

}

#################################################################################

# Runs BWA (fast read aligner) against two files with FASTA sequences (Query: reads or short seqs, Sbjct: reference genome)
sub execute_bwa {

	my ($query_files,$sbjct_file,$bwa_options,$outfile) = @_;

	# $query_files should be a reference array
	if (ref($query_files) ne "ARRAY"){
		$query_files = [$query_files];
	}

	if (!defined($bwa_options)){ $bwa_options=''; }

	my $random_name = random_file_name();

	if (defined($outfile)){
		$outfile =~ s/\.bam$//; # Bowtie adds the extension automatically
	} else {
		$outfile="/tmp/$random_name.bwa";
	}
	unlink($outfile);

	my $bwadb_file;
	if ($sbjct_file eq $NCBI_NT_DATABASE) {
		$bwadb_file = $NCBI_NT_DATABASE;
	} elsif (-e "$sbjct_file.bwt" && -e "$sbjct_file.pac" && -e "$sbjct_file.ann"
		&& -e "$sbjct_file.amb" && -e "$sbjct_file.sa") {
		$bwadb_file = $sbjct_file;
	} elsif ($sbjct_file =~ /(.+)\.(fa|fasta|fq|fastq)/ && -e "$1.bwt" && -e "$1.pac" && -e "$1.ann"
		&& -e "$1.amb" && -e "$1.sa") {
		$bwadb_file = $1;
	} elsif ($sbjct_file =~ /(.+\.(fa|fasta|fq|fastq))/ && -e "$1.bwt" && -e "$1.pac" && -e "$1.ann"
		&& -e "$1.amb" && -e "$1.sa") {
		$bwadb_file = $1;
	} elsif (-e "$sbjct_file.bwt" && -e "$sbjct_file.pac" && -e "$sbjct_file.ann"
		&& -e "$sbjct_file.amb" && -e "$sbjct_file.sa") {
		$bwadb_file = $1;
	} else {
		$bwadb_file = "/tmp/align_seqs_$random_name.db";
		# print "\n$BWABUILDEXE '$sbjct_file' '$bwadb_file'\n"; exit;
		`$BWABUILDEXE -p '$bwadb_file' '$sbjct_file' 2>&1`;
	}
	# BWA does not deal with GZIP files
	my $file_type1 = is_compressed($query_files->[0]);
	if ($file_type1){
		$query_files->[0] = "<(zcat '".$query_files->[0]."')";
	} else {
		$query_files->[0] = "'".$query_files->[0]."'";
	}
	if (scalar @$query_files == 2){
		my $file_type2 = is_compressed($query_files->[1]);
		if ($file_type2){
			$query_files->[1] = "<(zcat '".$query_files->[1]."')";
		} else {
			$query_files->[1] = "'".$query_files->[1]."'";
		}
	}
	# print "bash -c '$BWAMEMEXE $bwa_options -v 1 '$bwadb_file' ".join(' ',@$query_files)."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'\n"; exit;
	system("bash -c '$BWAMEMEXE $bwa_options -v 1 '$bwadb_file' ".join(' ',@$query_files)."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'");
	
	if ($bwadb_file =~ /\/tmp\/align_seqs_/) {
		`rm $bwadb_file.*.ebwt`;
	}
	
	return "$outfile.bam";

}

#################################################################################

# Runs BOWTIE (fast read aligner) against two files with FASTA sequences (Query: reads or short seqs, Sbjct: reference genome)
sub execute_bowtie {

	my ($query_files,$sbjct_file,$bowtie_options,$outfile) = @_;

	# $query_files should be a reference array
	if (ref($query_files) ne "ARRAY"){
		$query_files = [$query_files];
	}

	if (!defined($bowtie_options)){ $bowtie_options=''; }

	my $random_name = random_file_name();

	if (defined($outfile)){
		$outfile =~ s/\.bam$//; # Bowtie adds the extension automatically
	} else {
		$outfile="/tmp/$random_name.bowtie";
	}
	unlink($outfile);

	my $bowtiedb_file;
	if ($sbjct_file eq $NCBI_NT_DATABASE) {
		$bowtiedb_file = $NCBI_NT_DATABASE;
	} elsif (-e "$sbjct_file.1.ebwt" && -e "$sbjct_file.2.ebwt" && -e "$sbjct_file.3.ebwt"
		&& -e "$sbjct_file.4.ebwt" && -e "$sbjct_file.rev.1.ebwt"  && -e "$sbjct_file.rev.2.ebwt") {
		$bowtiedb_file = $sbjct_file;
	} elsif ($sbjct_file =~ /(.+)\.(fa|fasta|fq|fastq)/ && -e "$1.1.ebwt" && -e "$1.2.ebwt" && -e "$1.3.ebwt"
		&& -e "$1.4.ebwt" && -e "$1.rev.1.ebwt"  && -e "$1.rev.2.ebwt") {
		$bowtiedb_file = $1;
	} elsif ($sbjct_file =~ /(.+\.(fa|fasta|fq|fastq))/ && -e "$1.1.ebwt" && -e "$1.2.ebwt" && -e "$1.3.ebwt"
		&& -e "$1.4.ebwt" && -e "$1.rev.1.ebwt"  && -e "$1.rev.2.ebwt") {
		$bowtiedb_file = $1;
	} elsif (-e "$sbjct_file.1.ebwt" && -e "$sbjct_file.2.ebwt" && -e "$sbjct_file.3.ebwt"
		&& -e "$sbjct_file.4.ebwt" && -e "$sbjct_file.rev.1.ebwt"  && -e "$sbjct_file.rev.2.ebwt") {
		$bowtiedb_file = $1;
	} elsif (is_compressed($sbjct_file)) {
		$bowtiedb_file = "/tmp/align_seqs_$random_name.db";
		my $sbjct_file_ = decompress($sbjct_file);
		# print "\n$BOWTIEBUILDEXE '$sbjct_file_' '$bowtiedb_file'\n";exit;
		`$BOWTIEBUILDEXE '$sbjct_file_' '$bowtiedb_file' 2>&1`;
		`rm $sbjct_file_`;
	} else {
		$bowtiedb_file = "/tmp/align_seqs_$random_name.db";
		# print "\n$BOWTIEBUILDEXE '$sbjct_file' '$bowtiedb_file'\n"; exit;
		`$BOWTIEBUILDEXE '$sbjct_file' '$bowtiedb_file' 2>&1`;
	}
	# Bowtie '-S' option outputs standard SAM/BAM format with @SQ headers
	# If there are 2 files of paired-end reads
	if (scalar @$query_files == 2){
# 		# Bowtie does not deal with GZIP files
		my $file_type1 = is_compressed($query_files->[0]);
		if ($file_type1){
			$query_files->[0] = "<(zcat '".$query_files->[0]."')";
		} else {
			$query_files->[0] = "'".$query_files->[0]."'";
		}
		my $file_type2 = is_compressed($query_files->[1]);
		if ($file_type2){
			$query_files->[1] = "<(zcat '".$query_files->[1]."')";
		} else {
			$query_files->[1] = "'".$query_files->[1]."'";
		}
		# print "bash -c '$BOWTIEEXE $bowtie_options -S --quiet '$bowtiedb_file' -1 ".$query_files->[0]." -2 ".$query_files->[1]."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'\n"; exit;
		system("bash -c '$BOWTIEEXE $bowtie_options -S --quiet '$bowtiedb_file' -1 ".$query_files->[0]." -2 ".$query_files->[1]."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'");
	} else {
		my $is_fasta = is_fasta($query_files->[0]);
		my $file_type = is_compressed($query_files->[0]);
		if ($file_type){
			$query_files->[0] = "<(zcat '".$query_files->[0]."')";
		} else {
			$query_files->[0] = "'".$query_files->[0]."'";
		}
		if ($is_fasta){
			# print "bash -c '$BOWTIEEXE $bowtie_options -S --quiet '$bowtiedb_file' -f ".$query_files->[0]."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'\n"; exit;
			system("bash -c '$BOWTIEEXE $bowtie_options -S --quiet '$bowtiedb_file' -f ".$query_files->[0]."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'");
		} else {
			# print "bash -c '$BOWTIEEXE $bowtie_options -S --quiet '$bowtiedb_file' ".$query_files->[0]."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'\n"; exit;
			system("bash -c '$BOWTIEEXE $bowtie_options -S --quiet '$bowtiedb_file' ".$query_files->[0]."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'");
		}
	}
	if ($bowtiedb_file =~ /\/tmp\/align_seqs_/) {
		`rm $bowtiedb_file.*.ebwt`;
	}
	print "$outfile.bam\n";exit;
	return "$outfile.bam";

}

#################################################################################

# Runs BOWTIE2 (fast read aligner) against two files with FASTA sequences (Query: reads or short seqs, Sbjct: reference genome)
sub execute_bowtie2 {

	my ($query_files,$sbjct_file,$bowtie_options,$outfile) = @_;

	# $query_files should be a reference array
	if (ref($query_files) ne "ARRAY"){
		$query_files = [$query_files];
	}

	if (!defined($bowtie_options)){$bowtie_options='--no-unal --very-fast';}

	my $random_name = random_file_name();

	if (defined($outfile)){
		$outfile =~ s/\.bam$//; # Bowtie adds the extension automatically
	} else {
		$outfile="/tmp/$random_name.bowtie2";
	}
	unlink($outfile);

	my $bowtie2db_file;
	if ($sbjct_file eq $NCBI_NT_DATABASE) {
		$bowtie2db_file = $NCBI_NT_DATABASE;
	} elsif (-e "$sbjct_file.1.bt2" && -e "$sbjct_file.2.bt2" && -e "$sbjct_file.3.bt2"
		&& -e "$sbjct_file.4.bt2" && -e "$sbjct_file.rev.1.bt2"  && -e "$sbjct_file.rev.2.bt2") {
		$bowtie2db_file = $sbjct_file;
	} elsif ($sbjct_file =~ /(.+)\.(fa|fasta|fq|fastq)/ && -e "$1.1.bt2" && -e "$1.2.bt2" && -e "$1.3.bt2"
		&& -e "$1.4.bt2" && -e "$1.rev.1.bt2"  && -e "$1.rev.2.bt2") {
		$bowtie2db_file = $1;
	} elsif ($sbjct_file =~ /(.+\.(fa|fasta|fq|fastq))/ && -e "$1.1.bt2" && -e "$1.2.bt2" && -e "$1.3.bt2"
		&& -e "$1.4.bt2" && -e "$1.rev.1.bt2"  && -e "$1.rev.2.bt2") {
		$bowtie2db_file = $1;
	} elsif (-e "$sbjct_file.1.bt2" && -e "$sbjct_file.2.bt2" && -e "$sbjct_file.3.bt2"
		&& -e "$sbjct_file.4.bt2" && -e "$sbjct_file.rev.1.bt2"  && -e "$sbjct_file.rev.2.bt2") {
		$bowtie2db_file = $1;
	} elsif (is_compressed($sbjct_file)) {
		$bowtie2db_file = "/tmp/align_seqs_$random_name.db";
# 		print "\nzcat '$sbjct_file' | $BOWTIE2BUILDEXE '$bowtie2db_file'\n"; exit;
		`zcat '$sbjct_file' | $BOWTIE2BUILDEXE '$bowtie2db_file' 2>&1`;
	} else {
		$bowtie2db_file = "/tmp/align_seqs_$random_name.db";
		# print "\n$BOWTIE2BUILDEXE '$sbjct_file' '$bowtie2db_file'\n"; exit;
		`$BOWTIE2BUILDEXE '$sbjct_file' '$bowtie2db_file' 2>&1`;
	}
	# If there are 2 files of paired-end reads
	if (scalar @$query_files == 2){
# 		# New Bowtie2 version deals with GZIP files
# 		my $file_type1 = is_compressed($query_files->[0]);
# 		if ($file_type1 && $file_type1 ne 'gzip'){  # Bowtie2 accepts gzip files as input
# 			print "\nUncompressing '".$query_files->[0]."'\n";
# 			$query_files->[0] = decompress($query_files->[0],$file_type1);
# 		}
# 		my $file_type2 = is_compressed($query_files->[1]);
# 		if ($file_type2 && $file_type2 ne 'gzip'){
# 			print "\nUncompressing '".$query_files->[1]."'\n";
# 			$query_files->[1] = decompress($query_files->[1],$file_type2);
# 		}
		# print "\n$BOWTIE2EXE $bowtie_options -x '$bowtie2db_file' -1 '".$query_file->[0]."' -2 '".$query_file->[1]."' -S '$outfile'\n";exit;
		# system("$BOWTIE2EXE $bowtie_options -x '$bowtie2db_file' -1 '".$query_files->[0]."' -2 '".$query_files->[1]."' -S '$outfile'");
		# print "$BOWTIE2EXE $bowtie_options -x '$bowtie2db_file' -1 '".$query_files->[0]."' -2 '".$query_files->[1]."' | $SAMTOOLS view -Sb - | $SAMTOOLS sort - '$outfile'\n"; exit;
		system("$BOWTIE2EXE $bowtie_options --quiet -x '$bowtie2db_file' -1 '".$query_files->[0]."' -2 '".$query_files->[1]."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'");
# 		if ($file_type1 && $file_type1 ne 'gzip'){
# 			system('rm '.$query_files->[0]);
# 		}
# 		if ($file_type2 && $file_type2 ne 'gzip'){
# 			system('rm '.$query_files->[1]);
# 		}
	} else {
# 		# New Bowtie2 version deals with GZIP files
# 		my $file_type = is_compressed($query_files->[0]);
# 		if ($file_type && $file_type ne 'gzip'){  # Bowtie2 accepts gzip files as input
# 			print "\nUncompressing '".$query_files->[0]."'\n";
# 			$query_files->[0] = decompress($query_files->[0],$file_type);
# 		}
		# Eg. bowtie2 -p 20 --no-unal -x lambda_virus -U reads.fq.gz -S lambda_virus.sam
		if (is_fasta($query_files->[0])){
# 			print "$BOWTIE2EXE $bowtie_options -x '$bowtie2db_file' -f '".$query_files->[0]."' | $SAMTOOLS view -Sb - | $SAMTOOLS sort - '$outfile'\n"; exit;
			system("$BOWTIE2EXE $bowtie_options --quiet -x '$bowtie2db_file' -f '".$query_files->[0]."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'");
		} else {
			# print "\n$BOWTIE2EXE $bowtie_options -x '$bowtie2db_file' -U '".$query_files->[0]."' -S '$outfile'\n";exit;
			# system("$BOWTIE2EXE $bowtie_options -x '$bowtie2db_file' -U '".$query_files->[0]."' -S '$outfile'");
			# print "$BOWTIE2EXE $bowtie_options --quiet  -x '$bowtie2db_file' -U '".$query_files->[0]."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'\n"; exit;
			system("$BOWTIE2EXE $bowtie_options --quiet -x '$bowtie2db_file' -U '".$query_files->[0]."' | $SAMTOOLS view -Sb -F 4 - | $SAMTOOLS sort - '$outfile'");
		}
# 		if ($file_type && $file_type ne 'gzip'){
# 			system('rm '.$query_files->[0]);
# 		}
	}
	if ($bowtie2db_file =~ /\/tmp\/align_seqs_/) {
		`rm $bowtie2db_file.*.bt2`;
	}
	
	return "$outfile.bam";

}

#################################################################################

# Runs SAMTOOLS
sub execute_samtools {

	my ($infile,$options,$outfile) = @_;

	if (!defined($options) || $options !~ /view|sort|index|mpileup/){
		return undef;
	}

	my $random_name = random_file_name();
	if (!defined($outfile)){
		$outfile="/tmp/$random_name";
	}

	if ($options =~/^view/){
# 		print "$SAMTOOLS $options -o $outfile $infile\n"; exit;
		system("$SAMTOOLS $options -o $outfile $infile");
	} elsif ($options =~/^(sort|index)/){
# 		print "$SAMTOOLS $options $infile $outfile\n"; exit;
		system("$SAMTOOLS $options $infile $outfile");
	} elsif ($options =~/^mpileup/){
# 		print "$SAMTOOLS $options $infile > $outfile\n"; exit;
		system("$SAMTOOLS $options $infile > $outfile");
	}

	return $outfile;

}

#################################################################################

# Runs BCFTOOLS
sub execute_bcftools {

	my ($infile,$options,$outfile) = @_;

	if (!defined($options) || $options !~ /view|index/){
		return undef;
	}

	my $random_name = random_file_name();
	if (!defined($outfile)){
		$outfile="/tmp/$random_name";
	}
	unlink($outfile);

	if ($options =~/^(view|index)/){
# 		print "$BCFTOOLS $options $infile > $outfile\n"; exit;
		system("$BCFTOOLS $options $infile > $outfile");
	}

	return $outfile;

}

#################################################################################

# Runs AWK against two files with one sequence per line (Sbjct: short seqs or patterns to match in query seq file)
sub execute_regex_match {

	my ($query_file,$sbjct_file,$regex_options,$outfile) = @_;

	if (!defined($regex_options)){$regex_options='';}

	if (!defined($outfile)){
		my $random_name = random_file_name();
		$outfile="/tmp/$random_name.regex";
	}
	if (-e $outfile){
		`rm "$outfile"`;
	}
	
# 	print "perl -e 'open(PRIMERFILE,\$ARGV[0]);while(<PRIMERFILE>){chomp;if(/>/){\$primer_name=\$_;}elsif(\$_){\$primers{\$primer_name}=\$_;}}close(PRIMERFILE);open(SEQSFILE,\$ARGV[1]);while(<SEQSFILE>){chomp;if(/>/){\$query_name=\$_;}elsif(\$_){foreach \$primer_name(keys%primers){if(/\$primers{\$primer_name}/){printf(\"\%s\\t\%s\\t\%s\\t\%s\\t\%s\\n\",\$query_name,\$primer_name,\$-[0]+1,length(\$&),\$&);}}}}close(SEQSFILE);' $sbjct_file $query_file > $outfile\n";exit;
	system("perl -e 'open(PRIMERFILE,\$ARGV[0]);while(<PRIMERFILE>){chomp;if(/>/){\$primer_name=\$_;}elsif(\$_){\$primers{\$primer_name}=\$_;}}close(PRIMERFILE);open(SEQSFILE,\$ARGV[1]);while(<SEQSFILE>){chomp;if(/>/){\$query_name=\$_;}elsif(\$_){foreach \$primer_name(keys%primers){if(/\$primers{\$primer_name}/){printf(\"\%s\\t\%s\\t\%s\\t\%s\\t\%s\\n\",\$query_name,\$primer_name,\$-[0]+1,length(\$&),\$&);}}}}close(SEQSFILE);' $sbjct_file $query_file > $outfile");
# 	system("awk 'FNR==NR{a[\$0];next} {for (i in a) {if (p=match(\$0, i)) print FNR, RSTART, RLENGTH, substr(\$0, RSTART, RLENGTH) } }' $sbjct_file $query_file > $outfile");
# 	print "awk 'FNR==NR{ if (\$0~/>/){primer_name=\$0;} else if (\$0) {primers[primer_name]=\$0;} next;} { if (\$0~/>/) {query_name=\$0;} else if (\$0) { for (primer_name in primers) { if (p=match(\$0, primers[primer_name])) { print query_name\"\\t\"primer_name\"\\t\"RSTART\"\\t\"RLENGTH\"\\t\"substr(\$0, RSTART, RLENGTH); } } } }' $sbjct_file $query_file > $outfile\n";exit;
# 	system("awk 'FNR==NR{ if (\$0~/>/){primer_name=\$0;} else if (\$0) {primers[primer_name]=\$0;} next;} { if (\$0~/>/) {query_name=\$0;} else if (\$0) { for (primer_name in primers) { if (p=match(\$0, primers[primer_name])) { print query_name\"\\t\"primer_name\"\\t\"RSTART\"\\t\"RLENGTH\"\\t\"substr(\$0, RSTART, RLENGTH); } } } }' $sbjct_file $query_file > $outfile");
	
	# Output example: Query, Sbjct, Query_start, Length, Matched_string
	# >QK4PW:00031:00135       >Sbjct 1        3       21      GTTGTGTCTTTAACTCCACTG
	# >QK4PW:00039:00123       >Sbjct 1        7       21      GTTGTGTCTTTAACTCCACTG
	# >QK4PW:00022:00283       >Sbjct 1        27      21      GTTGTGTCTTTAACTCCACTG

	return $outfile;

}

#################################################################################

# Runs CD-HIT (clustering tool) against a file with FASTA sequences
sub execute_cdhit {

	my ($query_file,$identity_threshold,$cdhit_options,$outfile) = @_;


# 	if (!defined($cdhit_options) || $cdhit_options eq ''){ $cdhit_options='-n 2 -g 1'; }
	if (!defined($cdhit_options)) {
		$cdhit_options = '';
	}
	if (defined($identity_threshold)){
		if ($identity_threshold>1 && $identity_threshold<=100) {
			$identity_threshold = $identity_threshold/100;
		} 
		$cdhit_options .= " -c $identity_threshold";
	}


	if (!defined($outfile)){
		my $random_name = random_file_name();
		$outfile="/tmp/$random_name.cdhit.fasta";
	}
	unlink($outfile);

# 	print "\n$CDHITEXE -i '$query_file' -o '$outfile' -c $identity_threshold $cdhit_options\n";
	`$CDHITEXE -d 0 -i '$query_file' -o '$outfile'  $cdhit_options 2>&1`;

	return ($outfile, "$outfile.clstr");


}

#################################################################################

# Runs CD-HIT-EST (DNA clustering tool) against a file with FASTA sequences
sub execute_cdhit_est {

	my ($query_file,$identity_threshold,$cdhit_options,$outfile) = @_;

	if (!defined($cdhit_options)) {
		$cdhit_options = '';
	}
	if (defined($identity_threshold)){
		if ($identity_threshold>1 && $identity_threshold<=100) {
			$identity_threshold = $identity_threshold/100;
		} 
		$cdhit_options .= " -c $identity_threshold";
	}

	if (!defined($outfile)){
		my $random_name = random_file_name();
		$outfile="/tmp/$random_name.cdhit";
	}
	unlink($outfile);

# 	print "\n$CDHITESTEXE -d 0 -i '$query_file' -o '$outfile' $cdhit_options\n";
# 	exit;
	`$CDHITESTEXE -d 0 -i '$query_file' -o '$outfile' $cdhit_options 2>&1`;

	return ($outfile, "$outfile.clstr");

}

#################################################################################

# Runs CD-HIT-454 (clustering tool) against a file with FASTA sequences
sub execute_cdhit_454 {

	my ($query_file,$identity_threshold,$cdhit_options,$outfile) = @_;

	if (!defined($cdhit_options)) {
		$cdhit_options = '';
	}
	if (defined($identity_threshold)){
		if ($identity_threshold>1 && $identity_threshold<=100) {
			$identity_threshold = $identity_threshold/100;
		} 
		$cdhit_options .= " -c $identity_threshold";
	}

	if (!defined($outfile)){
		my $random_name = random_file_name();
		$outfile="/tmp/$random_name.cdhit";
	}
	unlink($outfile);

	`$CDHIT454EXE -d 0 -i '$query_file' -o '$outfile' $cdhit_options 2>&1`;

	return ($outfile, "$outfile.clstr");

}

#################################################################################

# Runs MAFFT multiple alignment against a file with FASTA sequences
sub execute_mafft {

	my ($query_file,$mafft_options,$outfile) = @_;

	if (!defined($mafft_options)){$mafft_options='';}

	if (!defined($outfile)){
		my $random_name = random_file_name();
		$outfile="/tmp/$random_name.mafft";
	}
	unlink($outfile);

# 	print "\n$MAFFTEXE --quiet $mafft_options '$query_file' > '$outfile'\n";exit;
	system("$MAFFTEXE --quiet $mafft_options '$query_file' > '$outfile' ");

	return $outfile;

}

#################################################################################

# Runs FLASH merging program against two files with FASTQsequences
sub execute_flash {

	my ($reads_file1,$reads_file2,$flash_options,$outfile) = @_;

	if (!defined($flash_options)){$flash_options='';}

	my $random_name = random_file_name();
	if (!defined($outfile)){
		my $random_name2 = random_file_name();
		$outfile="/tmp/$random_name2.flash";
	}
	unlink($outfile);

# 	print("$FLASHEXE -q -d /tmp $flash_options '$reads_file1' '$reads_file2' -o $random_name --suffix='' 1>&- 2>&-\n");exit;
	`$FLASHEXE -q -d /tmp $flash_options '$reads_file1' '$reads_file2' -o $random_name --suffix='' 2>&1`;
	system("mv /tmp/$random_name.extendedFrags.fastq $outfile");
# 	system("mv /tmp/$random_name.notCombined_1.fastq $outfile.nc_1");
# 	system("mv /tmp/$random_name.notCombined_2.fastq $outfile.nc_2");
	system("rm /tmp/$random_name.*");

# 	open(FLASH,"$FLASHEXE -q -c $flash_options '$reads_file1' '$reads_file2' |") || die "# $0 : cannot run $FLASHEXE\n";
# 	my ($sequences,$headers,$qualities) = read_fastq([<FLASH>]);
# 	close(FLASH);


	return $outfile ; # "$outfile.nc_1", "$outfile.nc_2";

}

#################################################################################

# Extracts data from BLASTN/BLASTP results
sub extract_blast_data {

	# parses blast alignments in order to extract alignment data,

	my ($infile,$limit,$format,$verbose) = @_;

	if (!defined($limit)){$limit=0;}
	if (!defined($format)){$format=0;} # 0 = Default Pairwise format  
	if (!defined($verbose)){$verbose=0;}

	my $blast_data;

	my ($read1,$read2) = (0,0);
	open(BLAST,$infile) || die "# $0 : #extract_blast_data cannot read $infile\n";

	my ($aux0,$aux1,$aux2,$aux3);
	my ($Q,$S,$count);
	my ($firstQ,$lastQ,$seqQ,$lenQ);
	my ($firstS,$lastS,$seqS,$lenS);
	my ($long,$name,$sbjct,$new_sbjct,$align,$ident,$strand,$gaps,$mismatches,$ev,$score);
	my ($ID,$IID,$sim,$Isim,$total,$totalI,$unalignedI);

	$Q = $S = 0;
	$count = 0;
	my $count_queries=0;
	my ($read_query_name, $read_sbjct_name) = (0,0);

	if ($verbose == 1){
		print "\nReading '$infile' Blast alignments ";
	}
	my $count_alignments = 0;

	while (my $row=<BLAST>) {
		$row =~ s/\012\015?|\015\012?//;
		
		if ($format == 10){ # 10 = Comma-separated values
			# moaC,gi|15800534|ref|NP_286546.1|,100.00,161,0,0,1,161,1,161,3e-114,330
			# moaC,gi|170768970|ref|ZP_02903423.1|,,99.38,161,1,,0,1,161,1,161,9e-114,329 
			($name,$sbjct,$ident,$long,$mismatches,$gaps,$firstQ,$lastQ,$firstS,$lastS,$ev,$score) = split(',',$row);
			if ($count >= 0){
				push(@{$blast_data->{$name}},{ 'sbjct' => $sbjct,
								'firstQ' => $firstQ,
								'lastQ' => $lastQ,
								'firstS' => $firstS,
								'lastS' => $lastS,
								'id' => sprintf("%.0f", $ident*$long/100),
								'long' => $long,
								'gaps' => $gaps,
								'ev' => $ev,
								'score' => $score,
								});
			}
		} elsif ($format == 6 || $format == 7){ # 6 = tabular ; 7 = tabular with comment lines
			# # 17 hits found
			# SHAA004TF	sp|P0A916|OMPW_SHIFL	50.43	115	49	2	389	733	1	107	1e-26	108
			# SHAA004TF	sp|P0A915|OMPW_ECOLI	50.43	115	49	2	389	733	1	107	1e-26	108
			if ($row =~ /^#/) {
				next;
			}
			($name,$sbjct,$ident,$long,$mismatches,$gaps,$firstQ,$lastQ,$firstS,$lastS,$ev,$score) = split("\t",$row);
			if ($count >= 0){
				push(@{$blast_data->{$name}},{ 'sbjct' => $sbjct,
								'firstQ' => $firstQ,
								'lastQ' => $lastQ,
								'firstS' => $firstS,
								'lastS' => $lastS,
								'id' => $ident,
								'long' => $long,
								'gaps' => $gaps,
								'ev' => $ev,
								'score' => $score,
								});
			}
		} else {  # 0 = Default Pairwise format  
			if (!defined($limit) || $count<$limit || $limit==0 || $row =~ /^Query=\s+([^\s]+)/ || $row =~ /^Query=\s+([^\|]+)/) {
# 				if($row =~ /^Query=/ || $row =~ /^>\s?(.+)/ || $row =~ /^  Database:/){
# 				print ''
# 				}
				if (defined($score) && $score ne '' && ($row =~ /^ Score =/ || $row =~ /^Query=/ || $row =~ /^>\s?(.+)/ || $row =~ /^  Database:/ || $row =~ /^Lambda/ )){
				# Can exist several alignments of the same sequence
				#
				# > HLA10490 HLA-B*46:01:10, Human MHC Class I sequence (partial)
				# Length=1012
				# 
				#  Score = 89.8 bits (48),  Expect = 1e-17
				#  Identities = 61/67 (91%), Gaps = 2/67 (3%)
				#  Strand=Plus/Minus
				# 
				# Query  1577  GGCTCCCATCTCCGGATCAGGGGC-TCAGGCAGCCCCTCATGCTGCACATGGCATGTGTA  1635
				#              |||||||||||| || | |||||| || ||||||||||||||||||||||||||||||||
				# Sbjct  899   GGCTCCCATCTCAGGGTGAGGGGCTTC-GGCAGCCCCTCATGCTGCACATGGCATGTGTA  841
				# 
				# Query  1636  TTTCTGC  1642
				#              | |||||
				# Sbjct  840   TCTCTGC  834
				# 
				# 
				# > HLA03558 HLA-G*01:01:19, Human MHC Class I sequence (partial)
				# Length=822
				# 
				#  Score = 86.1 bits (46),  Expect = 1e-16
				#  Identities = 56/61 (92%), Gaps = 0/61 (0%)
				#  Strand=Plus/Minus
				# 
				# Query  1582  CCATCTCCGGATCAGGGGCTCAGGCAGCCCCTCATGCTGCACATGGCATGTGTATTTCTG  1641
				#              ||||||| | || |||||||| ||||||||||||||||||||||||||||||||| ||||
				# Sbjct  821   CCATCTCAGCATGAGGGGCTCCGGCAGCCCCTCATGCTGCACATGGCATGTGTATCTCTG  762
				# 
				# Query  1642  C  1642
				#              |
				# Sbjct  761   C  761

					if ($count >= 0) {
						$blast_data->{$name}[$count]{'sbjct'}=$sbjct;
						$blast_data->{$name}[$count]{'firstQ'}=$firstQ;
						$blast_data->{$name}[$count]{'lastQ'}=$lastQ;
						$blast_data->{$name}[$count]{'firstS'}=$firstS;
						$blast_data->{$name}[$count]{'lastS'}=$lastS;
						$blast_data->{$name}[$count]{'id'}=$ident;
						$blast_data->{$name}[$count]{'long'}=$long;
						if ($strand) { # Only DNA alignment
							$blast_data->{$name}[$count]{'strand'}=$strand;
						}
						$blast_data->{$name}[$count]{'gaps'}=$gaps;
						$blast_data->{$name}[$count]{'ev'}=$ev;
						$blast_data->{$name}[$count]{'score'}=$score;
						$blast_data->{$name}[$count]{'seqQ'}=$seqQ;
						$blast_data->{$name}[$count]{'align'}=$align;
						$blast_data->{$name}[$count]{'seqS'}=$seqS;
						$blast_data->{$name}[$count]{'lenQ'}=$lenQ;
						$blast_data->{$name}[$count]{'lenS'}=$lenS;

						$count++;

						# Each 100 alignments processed, a notification message is printed
						$count_alignments++;
						if ($verbose == 1 && $count_alignments % 1000 == 0) {
							print "$count_alignments " ;
						}

					}

					$Q = $S = 0;  $firstQ=$lastQ=$seqQ=$firstS=$lastS=$seqS=$align=$long=$ident=$strand=$gaps=$ev=$score='';
					$ID=$IID=$sim=$Isim=$total=$totalI=$unalignedI='';

					if($count<0){$count++;}

				}
				if($row =~ /^(Query|Sbjct)\s+/ && $row !~ /^[\s\w\-\d=]+$/){
					$row =~ s/[^\s\w\-\d]/-/g;
				}
				if($row =~ /^Query\s+(\d+)\s+([a-zA-Z\-]+)\s+(\d+)/)
				{
					#($aux0, $aux1, $aux2, $aux3) = split(/\s+/,$row);
					if(!$Q){   $firstQ = $1;  }
	# 					if(!$Q || ($Q && $lastQ==$1-1))
	# 					{
						$lastQ = $3;
						$seqQ .= $2;
						$Q++;
	# 					}
				}
				elsif($row =~ /^\s+(.+)/ && $seqQ && !$seqS)
				{
					$align .= $1;
				}
				elsif($row =~ /^Sbjct\s+(\d+)\s+([a-zA-Z\-]+)\s+(\d+)/)
				{
					#($aux0, $aux1, $aux2, $aux3) = split(/\s+/,$row);
					if(!$S){   $firstS = $1;  }
	# 					if(!$S || ($S && $lastS==$1-1))
	# 					{
						$lastS = $3;
						$seqS .= $2;
						$S++;
	# 					}
				}
				elsif($row =~ /^Length=(\d+)/){
					if ($read_query_name) {
						$lenQ = $1;
						$name =~ s/\012\015?|\015\012?//;
						$read_query_name = 0;
						# Create an empty array of results
						$blast_data->{$name} = [];
					}
					if ($read_sbjct_name) {
						$lenS = $1;
						$sbjct =~ s/\012\015?|\015\012?//;
						$read_sbjct_name = 0;
					}
				}
				elsif ($read_sbjct_name && $row =~ /\w+/){
					$row =~ s/\012\015?|\015\012?//;
					$sbjct .= "".$row;
				}
				elsif ($read_query_name && $row =~ /\w+/){
					$row =~ s/\012\015?|\015\012?//;
					$name .= " ".$row;
				}
				# } elsif ($row =~ /^Query=\s+([^\s]+)/ || $row =~ /^Query=\s+([^\|]+)/){
				elsif ($row =~ /^Query=\s*(.+)/){
					#$row =~ s/\012\015?|\015\012?|^\s+|\s+$//;
					#$row =~ /^Query=\s*([^\t|^\|]+)/;
					#$1 =~ /^([^\|]+)/;
					#$name=$1;
	# 					if ($1 eq 'D6X84:00005:00010'){
	# 						print ''
	# 					}
					# If Query is new one, for ex. in Blast2 the query is always the same
					my $new_query = $1;
					if (!defined($name) || $name !~ /^\Q$new_query/) {
						$name = $new_query;
						$read_query_name = 1;
						$count = 0;

						# Each 100 alignments processed, a notification message is printed
						$count_queries++;
						if ($verbose && ($count_queries) % 100==0) {
							print "$count_queries " ;
						}
					}
					$lenQ='';

				}
				elsif($row =~ /^(>|Subject=)\s*(.+)/) 
				{
					$sbjct = $2;
					$read_sbjct_name = 1;
					$lenS = '';
					#$sbjct = read_sequence_name($1,'blast');
	# if ($sbjct =~ /comp5_c1_seq2/){
	# print '';
	# }
				}
				elsif($row =~ /^ Score =/) ## process previous alignment and last one $row =~ /^  Database:/ ||
				{
					$row =~ /^ Score =\s+([\d\.]+) bits \((\d+)\),  Expect = ([\d\.e\-]+)/; #/Expect = ([^,]+)/;
					$ev = $3;
					$score = $1;
				}
				#  Identities = 61/67 (91%), Gaps = 2/67 (3%)
				elsif($row =~ /Identities = (\d+)\/(\d+) \(\d+%\), Gaps = (\d+)\/\d+ \(\d+%\)/){
					$ident = $1;
					$long = $2;
					$gaps = $3;
					#$id =~ s/(\(|\)|\%|,)//g;
				}
				# PROTEIN ALIGNMENT: Identities = 20/20 (100%), Positives = 20/20 (100%), Gaps = 0/20 (0%)
				elsif($row =~ /Identities = (\d+)\/(\d+) \(\d+%\), Positives = \d+\/\d+ \(\d+%\), Gaps = (\d+)\/\d+ \(\d+%\)/){
					$ident = $1;
					$long = $2;
					$gaps = $3;
					#$id =~ s/(\(|\)|\%|,)//g;
				}
				elsif($row =~ /Strand=(\w+\/\w+)/){
					$strand = $1;
				}

			}
		}
	}
	close(BLAST);

	if ($verbose == 1){
		print "$count_alignments \n\n" ;
	}


	return $blast_data;
}

#################################################################################

# Extracts only alignment data from BLASTN/BLASTP results
sub extract_align_blast_data{

	my ($blast_data) = @_;

	my $align_blast_data;

	foreach my $query_name (keys %{$blast_data}){
# 		my %sbjct_names;
		my $count_sbjct = 0;
# 		if ($query_name eq 'D6X84:00005:00010'){
# 			print ''
# 		}
		foreach my $blast_result (@{$blast_data->{$query_name}}){
			# Store only one Sbjct result per query (the best alignment result)
# 			if (!defined($sbjct_names{$blast_result->{'sbjct'}})){
# 				$sbjct_names{$blast_result->{'sbjct'}} = $count_sbjct;
				$count_sbjct++;
				my $blast_data_;
# 				# If desired, preserve previous Blastp data values
# 				if (defined($previous_blast_data) && $previous_blast_data){
# 					$blast_data_ = $blast_result;
# 				}
				$blast_data_->{'NAME'} = $blast_result->{'sbjct'};
				if (defined($blast_result->{'seqQ'}) && defined($blast_result->{'seqS'})){
					$blast_data_->{'ALIGN'} = $blast_result->{'seqQ'}."\n".$blast_result->{'seqS'};
					my ($query_cols, $sbjct_cols) = obtain_cols_from_alignment($blast_result->{'seqQ'}."\n".$blast_result->{'seqS'}, $blast_result->{'firstQ'}, $blast_result->{'firstS'}, $blast_result->{'strand'});
					$blast_data_->{'COLS'} = join(",", @{$query_cols})."\n".join(",", @{$sbjct_cols});
				}
# 				$blast_data_->{'SIMIL'} = sprintf('%.0f', $blast_result->{'posit'}*$blast_result->{'long'}/100);
				$blast_data_->{'IDENT'} = $blast_result->{'id'};
				$blast_data_->{'EVALUE'} = $blast_result->{'ev'};
				$blast_data_->{'ALIGNED'} = $blast_result->{'long'};
				push(@{$align_blast_data->{$query_name}}, $blast_data_);
# 			} #else {
# 				my $sbjct_pos = $sbjct_names{$blast_result->{'sbjct'}};
# 				my %blast_data_ = %{$align_blast_data->{$query_name}[$sbjct_pos]};
# 				$blast_data_{'ALIGN'} = (split("\n",$blast_data_{'ALIGN'}))[0].'...'.$blast_result->{'seqQ'}."\n".(split("\n",$blast_data_{'ALIGN'}))[1].'...'.$blast_result->{'seqS'};
# 				my ($query_cols, $sbjct_cols) = obtain_cols_from_alignment($blast_result->{'seqQ'}."\n".$blast_result->{'seqS'}, $blast_result->{'firstQ'}, $blast_result->{'firstS'});
# 				$blast_data_{'COLS'} = (split("\n",$blast_data_{'COLS'}))[0].','.join(",", @{$query_cols})."\n".(split("\n",$blast_data_{'COLS'}))[1].','.join(",", @{$sbjct_cols});
# 				$blast_data_{'SIMIL'} = sprintf('%.0f', $blast_data_{'SIMIL'} + ($blast_result->{'id'}*$blast_result->{'long'}/100));
# 				push(@{$align_blast_data->{$query_name}}, \%blast_data_);
# 			}
		}
	}

	return $align_blast_data;

}

#################################################################################

# Extracts only alignment data from SAM/BAM format alignments
sub extract_align_bam_data {

	my ($infile,$ref_data,$limit) = @_;

	if (!defined($limit)){$limit=0;}

	my $align_data;

	my $ref_seqs;
	if (ref($ref_data) eq "HASH"){
		$ref_seqs = $ref_data;
	} elsif (is_fasta($ref_data)){
		$ref_seqs = read_fasta_file_hash($ref_data,1);
	} else {
		print "# extract_align_bam_data: Unrecognized FASTA format file '$ref_data'.\n";
	}

	# SAM data example:
	# 232d9c3aa0b15a91b46eaa4642c79193        16      ENST00000265620.11      537     1       15S73M  *       0       0       GACCTCAATTTTGTTTCAGGACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAAC        IIIIIIIIIIIIII
	# 8bd2c52cadac2ef8a3ff5154d013e204        16      ENST00000265620.11      538     1       72M     *       0       0       CAGGACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAAC        IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	# 4cb9965bd13db4572fc14c22e4266510        0       ENST00000371085.7       587     1       61M     *       0       0       ACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAA   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	# 39a367591c84a3f373c7df2ba38b0b41        16      ENST00000371095.7       541     1       77M4S   *       0       0       CAGGACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAACTTCCAGTAA       IIIIIIIIIIIIIIIIIIIIII
	# e775a2ee8a41d2b252d4c92d5e8038b6        0       ENST00000371095.7       545     1       72M     *       0       0       ACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAACTTCC        IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	# 0c463aad1591beb227f46282763d5bc0        0       ENST00000371102.8       2469    1       11S77M  *       0       0       TCAATTTTGTTTCAGGACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAACTTCC        IIIIIIIIIIIIII
	# 5e3f383ccd14705fdd80535369d4c558        0       ENST00000371102.8       2483    1       63M     *       0       0       GCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAACTTCC IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

	if (is_sam($infile)){
		open(ALIGN_FILE,$infile) || die "# $0 : # extract_align_bam_data cannot read '$infile'.\n";
	} elsif (is_bam($infile)){
		open(ALIGN_FILE,"$SAMTOOLS view $infile |") || die "# $0 : # extract_align_bam_data cannot read '$infile'.\n";
	} else {
		print "# extract_align_bam_data: Unrecognized BAM or SAM format file '$infile'.\n";
		exit;
	}
	my $count = 0;
	while(<ALIGN_FILE>){
		s/\012\015?|\015\012?//;
		my @line = split("\t");
		if (scalar @line < 10){
			next;
		} elsif ($line[2] eq '*'){ # Empty sbjct name - no alignment
			next;
		}

		my $query_name = $line[0];
		my $sbjct_name = $line[2];
# if ($query_name eq 'SRR953264.37100' && $sbjct_name  eq 'ENST00000412167.6'){
# print '';
# }
		my $query_seq = $line[9];
		my $sbjct_seq = $ref_seqs->{$sbjct_name};
		my $firstS = $line[3];
		my $cigar = $line[5];
 		my $firstQ = 1;
		my $lastQ = length($line[9]);
		my (@query_cols, @sbjct_cols);
		my ($seqQ, $seqS) = ('','');
		my ($sbjct_pos,$query_pos) = ($firstS,$firstQ);
		my ($ident,$alignedQ,$alignedS) = (0,0,0);
		while ($cigar =~ /(\d+)([MIDNSHP=X]+)/g){
			# Ex. CIGAR: 3M1I3M1D5M
			# CCATACT-GAACTGACTAAC
			#     ACTAGAA-TGGCT
			foreach (1..$1){
				if ($2 eq 'M' || $2 eq '=' || $2 eq 'X'){
					push(@sbjct_cols, $sbjct_pos);
					my $S = substr($sbjct_seq,$sbjct_pos-1,1);
					$seqS .= $S;
					$sbjct_pos++;
					push(@query_cols,$query_pos);
					my $Q = substr($query_seq,$query_pos-1,1);
					$seqQ .= $Q;
					$query_pos++;
					if ($S eq $Q) {
						$ident++;
					}
					$alignedS++;
					$alignedQ++;
				} elsif ($2 eq 'I') {
					push(@query_cols,$query_pos);
					$seqQ .= substr($query_seq,$query_pos-1,1);
					$query_pos++;
					$seqS .= '-';
					$alignedQ++;
				} elsif ($2 eq 'D') {
					push(@sbjct_cols,$sbjct_pos);
					$seqS .= substr($sbjct_seq,$sbjct_pos-1,1);
					$sbjct_pos++;
					$seqQ .= '-';
					$alignedS++;
				} elsif ($2 eq 'S') {
					# Soft-clipping: bases in 5' and 3' of the read are NOT part of the alignment.
					$query_pos++;
				} elsif ($2 eq 'H') {
					# Hard-clipping: bases in 5' and 3' of the read are NOT part of the alignment AND those bases have been removed from the read sequence in the BAM file. The 'real' sequence length would be length(SEQ)+ count-of-hard-clipped-bases
					$query_pos++;
					$firstQ++;
				}
			}
		}
		# Flag 16 means reverse
		if ($line[1] eq 16) {
			# Reverses query columns from the original query sequence (query sequence in the SAM file is already reversed)
			# Sbjct is the reference sequence that is always in direct sense
			my $query_len = length($query_seq);
			@query_cols = map $query_len-$_+1, @query_cols;
		}
# if ($seqS ne $seqQ){
# print '';
# }
		my %align_data_;
		$align_data_{'NAME'} = $sbjct_name;
		$align_data_{'IDENT'} = $ident;
		$align_data_{'COLS'} = join(",", @query_cols)."\n".join(",", @sbjct_cols);
		$align_data_{'ALIGN'} = $seqQ."\n".$seqS;
		if ($alignedS >= $alignedQ){
			$align_data_{'ALIGNED'} = $alignedS;
		} else {
			$align_data_{'ALIGNED'} = $alignedQ;
		}
		push(@{$align_data->{$query_name}}, \%align_data_);
		$count++;

		if ($limit && $limit<=$count){
			last;
		}
	}
	close(ALIGN_FILE);

	# Reorder results by identity, by default they are not ordered
	foreach my $seq_name (keys %{$align_data}) {
		my $aligns = $align_data->{$seq_name};
		my @sorted_aligns = sort { $aligns->[$b]{'IDENT'}/ $aligns->[$b]{'ALIGNED'}<=> $aligns->[$a]{'IDENT'}/$aligns->[$a]{'ALIGNED'} } 0 .. $#{$aligns};
# if (scalar @sorted_aligns > 1){
# print '';
# }
		my $sorted_aligns = [map $aligns->[$_] , @sorted_aligns];
		$align_data->{$seq_name} = $sorted_aligns;
	}

	return $align_data;

}

#################################################################################

# Extracts data from GASSST primer alignments
sub extract_align_gassst_data {

	my ($infile,$limit) = @_;

	if (!defined($limit)){$limit=0;}

	my $gassst_data;

	# GASST output example:
	# BANK           7  ATCCTCGCCC-GACCACC          23                4CB0D:00854:00654
	#                   |||||||||| |||||||                    # mismatche(s): 0    gap(s): 1    e-value: 7e-04
	# QUERY          1  ATCCTCGCCCTGACCACC          18                MHC2DQA_F001
	# 

	# Open __gassst_match_pairs.txt and save data in %alignments
	open(ALIGN,$infile) || die "# $0 : # extract_align_gassst_data cannot read '$infile'.\n";
	# open(ALIGN,$out_align_file) || die "# $0 : cannot read $out_align_file, gassst failed (seq)!\n";
	my $count_matches=0;
	my ($seq_name, $seq_align, $seq_start_align, $seq_strand, $match_substitutions, $match_gaps, $match_evalue);
	my $read=0;
	while(<ALIGN>){
		s/\012\015?|\015\012?//;
		if (/^BANK/){
			my @line = split("\t");
			$seq_name = $line[2];
			$read=1;
			$count_matches=0;
			# if ($include && !in_array($include,$seq_name)){
				# $read=0;
				# next;
			# }
			$line[0] =~ /^BANK\s+(\d+)\s+([a-zA-Z\-]+)\s+(\d+)/;
			$seq_start_align = $1;
			$seq_align = $2;
			if ($1 < $3) {
				$seq_strand = 'Plus';
			} else {
				$seq_strand = 'Minus';
			}
		} elsif ($read && /# mismatche\(s\): (\d+)\s+gap\(s\): (\d+)\s+e-value: ([e\-\d\.]+)/) {
			$match_substitutions = $1;
			$match_gaps = $2;
			$match_evalue = $3;
		} elsif ($read && /^QUERY/) {
			my @line = split("\t");
			my $primer_name = $line[2];
			$line[0] =~ /^QUERY\s+(\d+)\s+([a-zA-Z\-]+)\s+(\d+)/;
			my $primer_start_align = $1;
			my $primer_align = $2;
			my $primer_strand;
			if ($1 < $3) {
				$primer_strand = 'Plus';
			} else {
				$primer_strand = 'Minus';
			}
			my %gassst_data_;
			$gassst_data_{'NAME'} = $primer_name;
			$gassst_data_{'EVALUE'} = $match_evalue;
			$gassst_data_{'ALIGNED'} = length($primer_align);
			$gassst_data_{'IDENT'} = $gassst_data_{'ALIGNED'}-$match_gaps-$match_substitutions;
			my ($query_cols, $sbjct_cols) = obtain_cols_from_alignment($seq_align."\n".$primer_align, $seq_start_align, $primer_start_align, "$seq_strand/$primer_strand");
			if ($seq_strand eq 'Minus'){
				$query_cols = [ reverse @$query_cols ];
				$sbjct_cols = [ reverse @$sbjct_cols ];
				$seq_align = iupac_reverse_complementary($seq_align);
				$primer_align = iupac_reverse_complementary($primer_align);
			}
			$gassst_data_{'COLS'} = join(",", @{$query_cols})."\n".join(",", @{$sbjct_cols});
			if (defined($seq_align) && defined($primer_align)){$gassst_data_{'ALIGN'}  = "$seq_align\n$primer_align";}else{$gassst_data_{'ALIGN'}  = '';}
			push (@{$gassst_data->{$seq_name}},\%gassst_data_);
			$count_matches++;
			map { undef $_ } ($seq_name, $seq_align, $seq_start_align, $seq_strand, $match_substitutions, $match_gaps, $match_evalue);
		}
		if ($limit && $limit<=$count_matches){
				last;
		}
	}
	close(ALIGN);

	# Reorder results by evalue, by default they are not ordered
	# Results were ordered by Query (short sequence), but we annotated Query as Sbjct and the long sequence as Query
	foreach my $seq_name (keys %{$gassst_data}) {
		my $aligns = $gassst_data->{$seq_name};
		my @sorted_aligns = sort { $aligns->[$b]{'IDENT'}/ $aligns->[$b]{'ALIGNED'}<=> $aligns->[$a]{'IDENT'}/$aligns->[$a]{'ALIGNED'} } 0 .. $#{$aligns};
		my $sorted_aligns = [map $aligns->[$_] , @sorted_aligns];
		$gassst_data->{$seq_name} = $sorted_aligns;
	}

	return $gassst_data;

}

#################################################################################

# Extracts alignment style data from REGEX matches
sub extract_align_regex_data {

	my ($infile,$limit) = @_;

	if (!defined($limit)){$limit=0;}

	my $regex_data;

	# Infile example: Query, Sbjct, Query_start, Length, Matched_string
	# >QK4PW:00031:00135       >Sbjct 1        3       21      GTTGTGTCTTTAACTCCACTG
	# >QK4PW:00039:00123       >Sbjct 1        7       21      GTTGTGTCTTTAACTCCACTG
	# >QK4PW:00022:00283       >Sbjct 1        27      21      GTTGTGTCTTTAACTCCACTG

	open(ALIGN,$infile) || die "# $0 : # extract_align_regex_data cannot read '$infile'.\n";
	my $count_matches=0;
	my ($seq_name, $seq_align, $seq_start_align, $seq_strand, $match_substitutions, $match_gaps, $match_evalue);
	my $read=0;
	while(<ALIGN>){
		s/\012\015?|\015\012?//;
		my ($query_name, $sbjct_name, $start, $length, $seq_align) = split("\t");
		$query_name =~ s/^>\s*//;
		$sbjct_name =~ s/^>\s*//;
		my %regex_data_;
		$regex_data_{'NAME'} = $sbjct_name;
		$regex_data_{'EVALUE'} = 0;
		$regex_data_{'ALIGNED'} = $regex_data_{'IDENT'} = $length;
		$regex_data_{'COLS'} = join(",", $start .. $start+$length-1)."\n".join(",", 1 .. $length);
		$regex_data_{'ALIGN'}  = "$seq_align\n$seq_align";
		push (@{$regex_data->{$query_name}},\%regex_data_);
		$count_matches++;
		if ($limit && $limit<=$count_matches){
				last;
		}
	}
	close(ALIGN);

	# Reorder results by alignment length, by default they are not ordered
	# Results were ordered by Query (short sequence), but we annotated Query as Sbjct and the long sequence as Query
	foreach my $seq_name (keys %{$regex_data}) {
		my $aligns = $regex_data->{$seq_name};
		my @sorted_aligns = sort { $aligns->[$b]{'IDENT'}<=> $aligns->[$a]{'IDENT'} } 0 .. $#{$aligns};
		my $sorted_aligns = [map $aligns->[$_] , @sorted_aligns];
		$regex_data->{$seq_name} = $sorted_aligns;
	}

	return $regex_data;

}

#################################################################################

# Extracts data from CD-HIT clustering
sub extract_cdhit_data {

	my ($clusters_file,$seqs_file) = @_;

	my ($seqs,$headers) = read_fasta_file($seqs_file);

	my $clusters;

	# Open __gassst_match_pairs.txt and save data in %alignments
	open(FILE,$clusters_file) || die "# $0 : # extract_cdhit_data cannot read '$clusters_file'.\n";
	my $read_cluster = 0;
	my $cluster;
	while(<FILE>){
		s/\012\015?|\015\012?//;
		if (/^>Cluster (\d+)/){
			$cluster = $1;
		# The following pattern recognizes results from cd-hit, cd-hit-est and cd-hit-454
		# cd-hit-454: 0       331nt, >aa7cb7c7a9c40b7628535cf5c1f6ac23... at 1:331:1:338/+/99.70%
		# cd-hit-est: 0       333nt, >99a2a23f451af6cedbf5372f908b6ae0... at +/99.70%
		# cd-hit: 0       328aa, >27f617de2e938b3a94509c5321fc7067... at 99.70%
		} elsif (/^\d+\s+(\d+)(nt|aa)\, >?\s?(.+)\.\.\. ((\*)|at.+?\+?\/?([\d\.]+)%)/){
			my $cluster_data_;
			$cluster_data_->{'length'} = $1;
			my ($found, $pos) = in_array($headers,"/^$3/",1);
			if ($found) {
				$cluster_data_->{'name'} = $headers->[$pos->[0]];
			} else{
				$cluster_data_->{'name'} = $3;
				print "Cluster sequence '$3' couldn't be found in original sequences file '$seqs_file'.\n";
			}
			if (defined($5) && $5 eq '*') {
				$cluster_data_->{'identity'} = 100;
				unshift(@{$clusters->[$cluster]}, $cluster_data_);
			} else {
				$cluster_data_->{'identity'} = int($6);
				push(@{$clusters->[$cluster]}, $cluster_data_);
			}
		}
	}
	close(FILE);
	
	return $clusters;
	
}

#################################################################################

# Extracts data from Variant Call Format (VCF) file
sub extract_vcf_data {

	my ($file) = @_;

	my $variants;

	my @data_headers;
	open(FILE,$file) || die "# $0 : # extract_vcf_data cannot read '$file'.\n";
	while(my $line = <FILE>){
		$line =~ s/\012\015?|\015\012?//;
		##fileformat=VCFv4.1
		if ($line =~ /^##/){
			next;
		#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DRX015047.sorted.bam
		} elsif ($line =~ /^#(.+)/){
			@data_headers = split("\t",$1);
		#BRAF	1023	.	A	G	212	.	DP=93;VDB=4.902075e-01;RPB=-6.858533e-01;AF1=0.5;AC1=1;DP4=32,34,12,12;MQ=42;FQ=215;PV4=1,1,1,1	PL	242,0,255
		#braf	1084	.	A	G	215	.	DP=123;VDB=5.236580e-01;RPB=-6.874843e-01;AF1=0.5;AC1=1;DP4=40,47,14,18;MQ=42;FQ=218;PV4=1,0.39,1,1	GT:PL:DP:GQ	0/1:245,0,255:119:99
		} else {
			my @values = split("\t",$line,-1); # -1 forces to split also empty data_headers at the end
			# Reads data_headers
			my %variant_data;
			map $variant_data{$data_headers[$_]} = $values[$_] , 0..$#data_headers;
			# Reference/Chromosome name
			my $ref_name = $variant_data{'CHROM'};
# if ($ref_name eq 'NF1' && $variant_data{'POS'} == 4111){
# print '';
# }
			# E.g. pseudo-VCF files generated by 'common_variants.pl'
			if (!defined($variant_data{'FORMAT'})){
				$variant_data{'GENOTYPE'} = $variant_data{'REF'}.','.$variant_data{'ALT'};
				push(@{$variants->{$ref_name}}, \%variant_data);
				next;
			}
			# 'FORMAT' data column contains the abbreviations/formats of the information in the last columns (GT: genotpye; PL: genotype likelihoods; DP: depth of alleles; GQ: genotype quality)
			# Last data columns contains genotyping data separated by ':' (ex. PL: genotype likelihoods in Phred format - 0 or lowest values are the most probable)
			my @genotyping_headers = split(":",$variant_data{'FORMAT'});
			my $index_data = (grep { $data_headers[$_] eq 'FORMAT' } 0..$#data_headers)[0]+1;
			my (@samples, $genotyping_data);
			for (my $i=$index_data; $i<=$#data_headers; $i++){
				my $sample = $data_headers[$i];
				push(@samples,$sample);
				my @genotyping_values = split(":",$variant_data{$sample});
				for (my $j=0; $j<=$#genotyping_headers; $j++){
					$genotyping_data->{$sample}{$genotyping_headers[$j]} = $genotyping_values[$j];
				}
			}
# 			# Annotates the variant info if genotype is different to the reference one
			foreach my $sample (@samples){
				if (defined($genotyping_data->{$sample}{'GT'}) && $genotyping_data->{$sample}{'GT'} ne '1/0'){
					my %variant;
					if ($genotyping_data->{$sample}{'GT'} eq '1/1'){
						$variant{'GENOTYPE'} = $variant_data{'REF'}.','.$variant_data{'ALT'};
					} elsif ($genotyping_data->{$sample}{'GT'} eq '0/1'){
						$variant{'GENOTYPE'} = $variant_data{'ALT'};
					}
					if (defined($genotyping_data->{$sample}{'DP'})){
						$variant{'DP'} = $genotyping_data->{$sample}{'DP'};
					}
					if (defined($genotyping_data->{$sample}{'GQ'})){
						$variant{'GQ'} = $genotyping_data->{$sample}{'GQ'};
					}
					$variant{'POS'} = $variant_data{'POS'};
					$variant{'REF'} = $variant_data{'REF'};
					#$variant{'QUAL'} = $variant_data{'QUAL'};
					#$variant{'GENOTYPE'} = join(",",@$genotype);
					push(@{$variants->{$sample}{$ref_name}}, \%variant);
				}
	# 			my $min_phred_likelihood = min(@genotype_phred_likelihoods);
	# 			# Generates indexes to order genotypes as likelihoods
	# 			my $index_genotypes = generate_genotypes([split(",",$variant_data{'REF'})],[split(",",$variant_data{'ALT'})]),
	# 			my $genotype;
	# 			# Finds the most probable genotype (lowest Phred likelihood)
	# 			foreach my $index (sort {$a<=>$b} keys %$index_genotypes) {
	# 				if (shift @genotype_phred_likelihoods == $min_phred_likelihood){
	# 					$genotype = $index_genotypes->{$index};
	# 					last;
	# 				}
	# 			}
	# 			# Annotates the variant info if genotype is different to the reference one
	# 			my ($identical,$diff_pos) = compare_arrays([split(",",$variant_data{'REF'})],$genotype);
	# 			if (!$identical){
	# 				my %variant;
	# 				$variant{'POS'} = $variant_data{'POS'};
	# 				$variant{'REF'} = $variant_data{'REF'};
	# 				$variant{'QUAL'} = $variant_data{'QUAL'};
	# 				$variant{'GENOTYPE'} = join(",",@$genotype);
	# 				push(@{$variants->{$sample}{$ref_name}}, \%variant);
	# 			}
			}
		}
	}
	close(FILE);
	
	return $variants;
	
}

#################################################################################

# Extracts coverage data from BAM file
sub extract_coverage_data {

	my ($infile) = @_;

	my $coverages;

# 	print '$SAMTOOLS view -b \''.$infile.'\' | $SAMTOOLS mpileup - | awk \'{print $1"\t"$4}\''."\n";exit;
	open(COVERAGE_DATA, '$SAMTOOLS view -b \''.$infile.'\' | $SAMTOOLS mpileup - | awk \'{print $1"\t"$4}\' |') || die "# $0 : extract_coverage_data cannot read '$infile'\n";
	# Outputs name of chromosome/reference and coverage (each line is one position in the reference)
	# P53     156
	# P53     143
	# P53     121
	# ...
	my ($ref_name,$prev_ref,$pos_coverage,@pos_coverages);
	while(<COVERAGE_DATA>){
		if (/(.+)\t(\d+)/){
			($ref_name,$pos_coverage) = ($1,$2);
			if (defined($prev_ref) && $prev_ref ne $ref_name){
				$coverages->{$prev_ref}{'mean'} = sprintf("%.0f",mean(@pos_coverages));
				$coverages->{$prev_ref}{'min'} = min(@pos_coverages);
				$coverages->{$prev_ref}{'max'} = max(@pos_coverages);
				$coverages->{$prev_ref}{'p90'} = percentile(\@pos_coverages,90);
				$coverages->{$prev_ref}{'p95'} = percentile(\@pos_coverages,95);
				undef(@pos_coverages);
			}
			push(@pos_coverages,$pos_coverage);
			$prev_ref = $ref_name;
		}
	}
	$coverages->{$prev_ref}{'mean'} = mean(@pos_coverages);
	$coverages->{$prev_ref}{'min'} = min(@pos_coverages);
	$coverages->{$prev_ref}{'max'} = max(@pos_coverages);
	$coverages->{$prev_ref}{'p90'} = percentile(\@pos_coverages,90);
	$coverages->{$prev_ref}{'p95'} = percentile(\@pos_coverages,95);
	close(COVERAGE_DATA);

	return $coverages;
	
}

#################################################################################

# Based on the reference and the alternative nucleotide variants,
# genotype permutations are calculated in the same order that SAMtools calculates genotype likelihoods
sub generate_genotypes {

	my ($ref_variants,$alt_variants) = @_;

	my $index_genotypes;
	
	# Assigns index values to the variants in order
	my %index_variants;
	my $index_value = 0;
	foreach my $variant ((@$ref_variants, @$alt_variants)){
		$index_variants{$variant} = $index_value;
		$index_value++;
	}
	
	# Assigns index values to the variant permutations
	# Formula: F(j,k) = (k*(k+1)/2)+j
	foreach my $variant1 ((@$ref_variants, @$alt_variants)){
		foreach my $variant2 ((@$ref_variants, @$alt_variants)){
			$index_value = ($index_variants{$variant2}*($index_variants{$variant2}+1)/2)+$index_variants{$variant1};
			# Annotates the genotypes only once, in the correct order, 1st ref and 2nd alt
			if (!defined($index_genotypes->{$index_value})){
				if ($variant1 ne $variant2){
					$index_genotypes->{$index_value} = [$variant1,$variant2];
				} else {
					$index_genotypes->{$index_value} = [$variant1];
				}
			}
		}
	}
	
	return $index_genotypes;

}

#################################################################################

# Extracts column numbers from an alignment
sub obtain_cols_from_alignment{

	my ($alignment, $query_first, $sbjct_first, $strand) = @_;

	my ($query_seq, $sbjct_seq) = split("\n",$alignment);

	my @query_seq = split('', $query_seq);
	my @sbjct_seq = split('', $sbjct_seq);

	my ($query_gaps, $sbjct_gaps) = (0,0);

	my (@query_cols, @sbjct_cols);

	for (my $i=0; $i<=$#query_seq; $i++){
		if ($query_seq[$i] ne '-' && $sbjct_seq[$i] ne '-'){
			if (!defined($strand) || $strand eq '' || $strand eq 'Plus/Plus'){
				push(@query_cols, $query_first+$i-$query_gaps);
				push(@sbjct_cols, $sbjct_first+$i-$sbjct_gaps);
			} elsif ($strand eq 'Plus/Minus'){
				push(@query_cols, $query_first+$i-$query_gaps);
				push(@sbjct_cols, $sbjct_first-$i+$sbjct_gaps);
			} elsif ($strand eq 'Minus/Plus'){
				push(@query_cols, $query_first-$i+$query_gaps);
				push(@sbjct_cols, $sbjct_first+$i-$sbjct_gaps);
			} elsif ($strand eq 'Minus/Minus'){
				push(@query_cols, $query_first-$i+$query_gaps);
				push(@sbjct_cols, $sbjct_first-$i+$sbjct_gaps);
			}
		} elsif ($query_seq[$i] eq '-' && $sbjct_seq[$i] ne '-'){
			$query_gaps++;
		} elsif ($query_seq[$i] ne '-' && $sbjct_seq[$i] eq '-'){
			$sbjct_gaps++;
		} 
	}

	return(\@query_cols, \@sbjct_cols);

}

#################################################################################

# Creates a single line alignment from two lines input one
sub alignment_to_single_line {

	my ($seq1, $seq2) = split("\n",$_[0]);

	if (!defined($seq1) || !defined($seq2)){
		return;
	}
	
	my $seq3;

	if (length($seq1)==length($seq2)){

		my @seq1 = split(//,$seq1);
		my @seq2 = split(//,$seq2);

		for (my $i=0; $i<=$#seq1; $i++) {
			if (lc($seq1[$i]) ne lc($seq2[$i])){
				$seq3 .= '['.$seq1[$i].'|'.$seq2[$i].']';
			} else {
				$seq3 .= $seq1[$i];
			}
		}

	} else  {
		die "\nERROR 'alignment_to_single_line': Both sequences must have the same length to be compared.\n\n";
	}

	return $seq3;

}

#################################################################################

# Stores data from a variable into a file  
sub store_data_dump {  
	my ($data, $file, $method) = @_;
	if (!defined($method)) { $method = 'Dumper'; }
	if ($method eq 'Dumper') {
		open(FILE, ">$file");  
		print FILE Dumper($data);  
		close FILE;  
	} elsif ($method eq 'Storable') {
		store $data, $file;
	}
}

#################################################################################

# Recovers data from a file into a variable  
sub recover_data_dump {  
	my ($file, $method) = @_;
	if (!defined($method)) { $method = 'Dumper'; }
	if ($method eq 'Dumper') {
		return do $file;
	} elsif ($method eq 'Storable') {
		return retrieve($file);
	}
}


#################################################################################

# Generates a MD5 hash from a string
sub generate_md5 {

	my ($string, $length, $options) = @_;
# if (!defined($string)){
# print '';
# }
	# Generate an unique MD5 string from the sequence to identify unequivocally the TF
	my $md5 = new Digest::MD5;
	$md5->add(uc($string));
	my $md5_hex = $md5->hexdigest;

	if (defined($options) && in_array($options, 'base64')){
		$md5_hex = encode_base64($md5_hex);
	}

	if (defined($length)){
		$md5_hex = substr($md5_hex,0,$length);
	}
	
	return $md5_hex;
}

#################################################################################

# Checks if a folder doesn't contain any file
sub is_folder_empty {
	my $dirname = shift;
	opendir(my $dh, $dirname) or die "Not a directory";
	return scalar(grep { $_ ne "." && $_ ne ".." } readdir($dh)) == 0;
}

#################################################################################

# Lists folder files
sub list_dir {
        my $dir = shift @_;
        my @files;
        find({ wanted => sub { push @files, $_ } , no_chdir => 1 }, ($dir));
        return @files;
}

#################################################################################

# Converts all characters into a string into alphanumeric ones
sub convert_to_alphanum {

	my $string = shift;
	
	$string = trim($string);
	$string =~ s/[^a-zA-Z0-9_]/_/g;
	
	return $string;
}

#################################################################################

# Retrieves a common Taxonomy from a list of different (or similar) OTUs
sub group_taxo {

	my ($input_taxa,$taxo_level) = @_;
	
# 	my @input_taxa = unique([@$input_taxa]);
	if (defined($taxo_level)){
		$taxo_level = lc($taxo_level);
	}
	
	# Defines which taxo fields to group
	my @taxo_fields = ('kingdom','domain','phylum','class','order','family','genus','species');
	my @group_fields = ('order','family','genus','species'); # ('kingdom','domain','phylum','class','order','family','genus','species');
	my @nogroup_fields = ('kingdom','domain','phylum','class');

	my $taxa_grouped;
	for (my $i=0; $i<=$#{$input_taxa}; $i++) {
		my $taxo = $input_taxa->[$i];
		# OTU ex.: tax=d:Bacteria,p:Actinobacteria,c:Actinobacteria,o:Actinomycetales,f:Nocardiaceae,g:Smaragdicoccus;
		# k: kingdom, d:domain, p: phylum, c: class, o: order, f: family, g: genus, s: species
		my $taxo_data = read_taxo_data($taxo);
		foreach my $taxo_field (@taxo_fields){
			if (defined($taxo_data->{$taxo_field})){
				$taxa_grouped->{$taxo_field}{$taxo_data->{$taxo_field}}++;
			}
			if (defined($taxo_level) && $taxo_field eq $taxo_level){
				last;
			}
		}
	}
	my $taxa_number = scalar @$input_taxa;
	my %taxa_consensus;
	foreach my $taxo_field (@taxo_fields){
		if (defined($taxa_grouped->{$taxo_field})){
			my @taxa_assigns = keys %{$taxa_grouped->{$taxo_field}};
			if ($#taxa_assigns == 0 && $taxa_grouped->{$taxo_field}{$taxa_assigns[0]} == $taxa_number){
				$taxa_consensus{$taxo_field} = $taxa_assigns[0];
			} else {
				last;
			}
# 			if ($#taxa_assigns == 0){
# 				$taxa_data_consensus{$taxa_field} = $taxa_assigns[0];
# 			} else {
# 				last;
# 			}
		}
	}

	my $taxa_consensus = print_taxo_data(\%taxa_consensus);
	
# 	$grouped_taxa->{$taxa_consensus} = $input_taxa;
	
	return $taxa_consensus;

}

################################################################################

# # Retrieves a common Taxonomy from a list of different (or similar) OTUs
# sub group_taxo {
# 
# 	my @input_otus = @_;
# 
# 	my $grouped_otus;
# 
# 	# Defines which taxo fields to group
# 	my @taxo_fields = ('kingdom','domain','phylum','class','order','family','genus','species');
# 	my @group_fields = ('order','family','genus','species'); # ('kingdom','domain','phylum','class','order','family','genus','species');
# 	my @nogroup_fields = ('kingdom','domain','phylum','class');
# 
# 	for (my $i=0; $i<=$#input_otus; $i++) {
# 		my $otu1 = $input_otus[$i];
# 		# OTU ex.: tax=d:Bacteria,p:Actinobacteria,c:Actinobacteria,o:Actinomycetales,f:Nocardiaceae,g:Smaragdicoccus;
# 		# k: kingdom, d:domain, p: phylum, c: class, o: order, f: family, g: genus, s: species
# 		my $otu1_data = read_taxo_data($otu1);
# 		my $common_fields = 0;
# 		my @common_otus;
# 		# Stores the first OTU data into a hash to remove non-common fields later
# # 		my %common_otu_data;
# # 		map $common_otu_data{$_} = $otu1_data{$_}, @group_fields;
# 		for (my $j=0; $j<=$#input_otus; $j++) {
# 			my $otu2 = $input_otus[$j];
# 			if ($otu1 eq $otu2) {
# 				next;
# 			}
# 			my $otu2_data = read_taxo_data($otu2);
# 			my $field_matches = 0;
# 			for (my $n=0; $n<=$#group_fields; $n++) {
# 				if (!defined($otu1_data{$group_fields[$n]}) || !defined($otu2_data{$group_fields[$n]})){
# 					last;
# 				} elsif ($otu1_data{$group_fields[$n]} ne $otu2_data{$group_fields[$n]}){
# 					last;
# 				} else {
# 					$field_matches++;
# 				}
# 			}
# 			 # Removes matched and annotates the lowest field matched
# 			if ($field_matches){
# 				push(@common_otus, $otu2);
# 				splice(@input_otus,$j,1); $j--;
# 				if (!$common_fields || $field_matches<$common_fields){
# 					$common_fields = $field_matches;
# 				}
# 			}
# 		}
# 		if ($common_fields){ # If they have taxo data in common then simplify
# 			unshift(@common_otus, $otu1);
# 			my %otu_data_consensus;
# 			map $otu_data_consensus{$_} = $otu1_data{$_} , @nogroup_fields;
# 			map $otu_data_consensus{$group_fields[$_]} = $otu1_data{$group_fields[$_]} , 0..$common_fields-1;
# 			my $otu_consensus = print_taxo_data(\%otu_data_consensus);
# 			push(@{$grouped_otus->{$otu_consensus}}, @common_otus);
# 		} else {
# 			my $otu_consensus = print_taxo_data(\%otu1_data);
# 			push(@{$grouped_otus->{$otu_consensus}}, $otu1);
# 		}
# 	}
# 	
# 	return ($grouped_otus);
# 
# }

#################################################################################

# Reads taxonomic information from a FASTA header in UTAX format
sub read_taxo_data {

	my ($data, $format) = @_;
	
	if (!defined($format)){
		$format = 'utax';
	}
	
	my %taxo_data;
	
	if ($format eq 'utax'){
		# AmpliTAXO (UTAX) header format:
		# >AB243007_S000622964;tax=d:Bacteria,p:Actinobacteria,c:Actinobacteria,o:Actinomycetales,f:Nocardiaceae,g:Smaragdicoccus;
		# k: kingdom, d:domain, p: phylum, c: class, o: order, f: family, g: genus, s: species
		$data =~ s/["']//g;
		if ($data =~ /.+?;\s?taxo?=(.+);?/ || $data =~ /taxo?=(.+);?/ || $data =~ /(.+);?/) {
			foreach my $taxo_field (split(/,\s?/,$1)){
				if ($taxo_field =~ /(\w:)?unidentified/i) {
					last;
				} elsif ($taxo_field =~ /k:(.+)/) {
					$taxo_data{'kingdom'} = $1;
				} elsif ($taxo_field =~ /d:(.+)/) {
					$taxo_data{'domain'} = $1;
				} elsif ($taxo_field =~ /p:(.+)/) {
					$taxo_data{'phylum'} = $1;
				} elsif ($taxo_field =~ /c:(.+)/) {
					$taxo_data{'class'} = $1;
				} elsif ($taxo_field =~ /o:(.+)/) {
					$taxo_data{'order'} = $1;
				} elsif ($taxo_field =~ /f:(.+)/) {
					$taxo_data{'family'} = $1;
				} elsif (!defined($taxo_data{'genus'}) && $taxo_field =~ /g:(.+)/) {
					$taxo_data{'genus'} = $1;
				} elsif (!defined($taxo_data{'species'}) && defined($taxo_data{'genus'}) && $taxo_field =~ /s:${taxo_data{'genus'}}[\s_](.+)/) {
					$taxo_data{'species'} = $1;
				} elsif (!defined($taxo_data{'species'}) && $taxo_field =~ /s:(.+);|s:(.+)/) {
					$taxo_data{'species'} = $1;
				}
			}
		}
		
	}

	return \%taxo_data;

}

#################################################################################

# Prints taxonomic information in UTAX format
sub print_taxo_data {

	my ($taxo_data, $format) = @_;
	
	if (!defined($format)){
		$format = 'utax';
	}
	
	my ($data,@data_fields);
	
	if ($format eq 'utax'){
		# AmpliTAXO (UTAX) header format:
		# >AB243007_S000622964;tax=d:Bacteria,p:Actinobacteria,c:Actinobacteria,o:Actinomycetales,f:Nocardiaceae,g:Smaragdicoccus;
		# k: kingdom, d:domain, p: phylum, c: class, o: order, f: family, g: genus, s: species
		# UTAX format fields
		my @utax_fields = ('kingdom','domain','phylum','class','order','family','genus','species');
		foreach my $utax_field (@utax_fields) {
			if (defined($taxo_data->{$utax_field})){
				push(@data_fields,sprintf("%s:%s",substr($utax_field,0,1),$taxo_data->{$utax_field}));
			}
		}
		if (@data_fields) {
			$data = sprintf("tax=%s;",join(',',@data_fields));
		} else {
			$data = "tax=Unassigned;";
		}
	}
	
	return $data;

}


#################################################################################

# Converts numeric matrix to a TRANSFAC matrix
sub matrices_to_transfac {

	my ($matrices) = @_;

	my $pwms;

	foreach my $name (keys %$matrices) {

		my $matrix = $matrices->{$name};
		
		my $pwm = "DE  $name\n";

		my $count_row = 0;
		foreach my $row (@{$matrix}){
			$count_row++;
			$pwm .= sprintf('%02s',$count_row);
			$pwm .= "\t".join("\t",@{$row})."\n";
		}
		$pwm .= "XX";
		chomp $pwm;
		push(@{$pwms}, $pwm);
	}


	return $pwms;

}

#################################################################################

# Converts a TRANSFAC matrix to a numeric matrix
sub transfac_to_matrices {

	my ($pwms,$scale) = @_;

	my $matrices;

	foreach my $pwm (@{$pwms}){
		#$pwm = normalize_transfac($pwm);
		my @pwm = split("\n",$pwm);
		my ($matrix,$name);
		foreach my $row (@pwm){
			if ($row=~/^(\d+)\s+(.+)/){
				#$row_=$2;
				#$row=~s/\012\015?|\015\012?|^\s+|\s+$//;
				my @row = split(/\s+/,$2);
				push(@{$matrix},[$row[0],$row[1],$row[2],$row[3]]);
			}elsif ($row=~/^DE\s+([^\t|^\||^;]+)/){
				$name=$1;
			}
		}
		$matrices->{$name} = $matrix;
	}


	return $matrices;

}

#################################################################################

# Normalizes a Transfac PWM matrix
sub normalize_transfac {

	my ($pwm,$scale) = @_;

	my $npwm;

	foreach my $row (split(/\n/,$pwm)){
		if ($row=~/^(\d+)\s+(.+)/){
			my $data=$2;
			my $pos=$1;
			$data=~s/^\s+|\s+$//;
			my @data = split(/\s+/,$data);
			my $old_scale = $data[0]+$data[1]+$data[2]+$data[3];
			my $abs_old_scale = abs($data[0])+abs($data[1])+abs($data[2])+abs($data[3]);
			
			my $new_row = sprintf('%02s',$pos);
			
			my @nts;
			# If old scale is lower than 1, we make it positive before changing
			if (defined($scale)){
				# If old scale is lower than -1, we make it between -1 and 0
				if ($data[0] < -1 || $data[1] < -1 || $data[2] < -1 || $data[3] < -1) {
					for (my $i=0;$i<4;$i++){
						$data[$i] = $data[$i] / $abs_old_scale;
					}
					$old_scale = $data[0]+$data[1]+$data[2]+$data[3];
				}
				# If scale is between 0 and -1
				if (($data[0] <= 0 && $data[0] > -1) && ($data[1] <= 0 && $data[1] > -1) && ($data[2] <= 0 && $data[2] > -1) && ($data[3] <= 0 && $data[3] > -1) ) {
					my $min_data = min($data[0], $data[1], $data[2], $data[3]);
					for (my $i=0;$i<4;$i++){
# 						$data[$i] = $data[$i] - ( $min_data );
						$data[$i] = ($data[$i]-$min_data ) * ($data[$i]-$min_data);
					}
					$old_scale = $data[0]+$data[1]+$data[2]+$data[3];
				}
				# If scale is between 0 and 1
				if (eval($old_scale) > 0 && eval($old_scale) < 1) {
					for (my $i=0;$i<4;$i++){
						$data[$i] = $data[$i] * 10;
					}
					$old_scale = $data[0]+$data[1]+$data[2]+$data[3];
				} 
			}
			for (my $i=0;$i<4;$i++){
				if (!defined($scale) && length($data[$i]) < 7){
					$new_row .= sprintf('%7s',$data[$i]);

				} else {
					my $new_value;
					if (!defined($scale) && length($data[$i]) >= 7){
						$new_value = $data[$i];
					} else {
						$new_value = $data[$i]*($scale/$old_scale);
					}
					if ($new_value > 10 || $new_value == 0 || int($new_value) == $new_value || sprintf("%.2f",($new_value - int($new_value))) == 0 || sprintf("%.2f",($new_value - int($new_value))) == 1) {#(($data[$i]*$scale) % $old_scale == 0) ){
						$new_row .= sprintf(' %6.0f', $new_value);
					} else {
						$new_row .= sprintf(' %6.2f', $new_value);
					}
				}
# 				if ($data[$i]>0 && $data[$i]/$old_scale>0.1){
# 					if ($i==0){
# 						push(@nts,'A');
# 					} elsif ($i==1){
# 						push(@nts,'C');
# 					} elsif ($i==2){
# 						push(@nts,'G');
# 					} elsif ($i==3){
# 						push(@nts,'T');
# 					}
# 				}
			}
			my $major_nt = transfac_to_iupac($new_row);
			if (defined($major_nt)) {
				$npwm .= $new_row.sprintf('%7s',transfac_to_iupac($new_row))."\n";
			} else {
				$npwm .= $new_row."\n";
			}
# 			$npwm .= sprintf('%7s',convert_nt_to_iupac(@nts))."\n";
		} elsif ($row =~ /^DE\s+(.+)/){
			$npwm .= "DE  $1\n";
			$npwm .= sprintf('P0%7s%7s%7s%7s','A','C','G','T')."\n";
		} elsif ($row =~ /^XX/){
			$npwm .= "XX\n";
		}
	}

	$npwm =~ s/\n$//s;

	return  $npwm;

}

#################################################################################

# Converts a PWM matrix (array of arrays) into a list of sequences
sub matrix_to_sequences {

	my ($matrix,$maxseqs,$alphabet) = @_;

	my @seqs;

	my @alphabet;
	if (defined($alphabet)){
		@alphabet = split('',$alphabet);
	} else {
		@alphabet = ('A','C','G','T');
	}

	if (!defined($maxseqs)){
		$maxseqs = 100;
	}
	
	# Convert negative values in positives if necessary
	my $is_negative = 0;
	for(my $pos=0; $pos<=$#{$matrix}; $pos++){
		for(my $col=0; $col<=$#{$matrix->[$pos]}; $col++){
			if (defined($matrix->[$pos][$col]) && $matrix->[$pos][$col] < 0) {
				$is_negative = 1;
				last;
			}
		}
		if ($is_negative) {
			last;
		}
	}
	if ($is_negative) {
		for(my $pos=0; $pos<=$#{$matrix}; $pos++){
			for(my $col=0; $col<=$#{$matrix->[$pos]}; $col++){
				$matrix->[$pos][$col] = exp($matrix->[$pos][$col]);
			}
		}
	}
	
	# Normalizes each row to the maximum number of sequences
	for(my $pos=0; $pos<=$#{$matrix}; $pos++){
		my $row_total = 0;
		map $row_total += $_ , @{$matrix->[$pos]};
		my @new_row;
		my $new_total = 0;
		for(my $col=0; $col<=$#{$matrix->[$pos]}; $col++){
			my $new_value = sprintf("%.0f",$matrix->[$pos][$col]*$maxseqs/$row_total);
			if ($new_total+$new_value >= $maxseqs){
				$new_row[$col] = $maxseqs - $new_total;
				last;
			}
			$new_row[$col] = $new_value;
			$new_total += $new_value;
		}
# 		$row_total = 0;
# 		map $row_total += $_ , @new_row;
# 		if ($row_total != $maxseqs){
# 			for(my $col=0; $col<=$#{$matrix->[$pos]}; $col++){
# 				$new_row[$col] = sprintf("%.0f",$new_row[$col]*$maxseqs/$row_total);
# 			}
# 		}
		$matrix->[$pos] = [ @new_row ];
	}

	for(my $pos=0; $pos<=$#{$matrix}; $pos++){
		my $total = 0;
		my $count_seq = 0;
		for(my $col=0; $col<=$#{$matrix->[$pos]}; $col++){
			for(my $count = 0; $count < $matrix->[$pos][$col]; $count++){ 
				$seqs[$count_seq] .= $alphabet[$col];
				$total++;
				$count_seq++;
			} 
		}
		for(my $count = $total; $count < $maxseqs; $count++){
			$seqs[$count_seq] .= '-';
			$count_seq++;
		}
	} #foreach $total (@PWMseqs){ print "$total->[0]\n" }
	
	return @seqs;
}

#################################################################################

# Converts a PWM matrix (array of arrays) into a PNG image logo
sub execute_weblogo {

	my ($seqs,$outfile,$outformat,$options)=@_;

	if (!defined($outformat)){
		$outformat = 'PNG';
	} else {
		$outformat = uc($outformat);
	}
	if (!defined($outfile)){
		$outfile = "/tmp/".random_file_name();
	}
	if (!defined($options)){
		$options = "-c classic -X NO -Y NO --fineprint '' -s large --resolution 300 --aspect-ratio 3 --errorbars NO";
	}

	my $seqsfile = $outfile.".seqs";
	open(SEQFILE,">$seqsfile") || die "# 'matrix_to_logo' : cannot create $seqsfile\n";
	print SEQFILE join("\n",@$seqs);
	close(SEQFILE);

	if ($outfile !~ /\.${outformat}$/i){
		$outfile .= '.'.lc($outformat);
	}
	system("$WEBLOGOEXE -f $seqsfile -F $outformat -o $outfile $options");
# 	print "\n$WEBLOGOEXE  -f $seqsfile -F $outformat -o $outfile $options\n";exit;
	`rm $seqsfile`;

	return $outfile;

}

#################################################################################

1;
