#!/usr/bin/perl -w
#
################################################################
#
# Name: ampliMERGE.pl
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
#   Merges paired-end reads from next-generation sequencing experiments using FLASH tool
#   Generates a unique FASTQ format file with only merged reads
#
# Requires as input several FASTQ files with sequences/reads and several CVS files in the same order with primer/tag data information
#
# Example: perl ampliMERGE.pl -i reads_R1.fq.gz reads_R2.fq.gz -d amplicon_data.csv -o merged_reads
#
# If the reads have been already demultiplexed into separate files (one file per sample), they can be packed into a single .zip or .tar.gz file and use it as input
# The packed paired-end files should have identical prefix names followed by the suffixes '_R1' and '_R2' plus dot and the file extension (e.g. reads_R1.fastq and reads_R2.fastq)
# Example: perl ampliMERGE.pl -i reads.tar.gz -o merged_reads
#

my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Merges Illumina paired-end reads from amplicon sequencing experiments.";


# Modules are in folder 'lib' in the path of the script
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;
use Statistics::Descriptive;
use warnings;
no warnings ('uninitialized', 'substr');

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options
# Maximum number of allowed primer+tag sequences
# my $MAX_PRIMER_TAG_SEQS = 70000;
# Minimum overlapping between reads
my $MIN_READ_OVERLAP = 5;

my (@INP_reads_files, $INP_amplicons_file, $INP_outfile, $INP_threads, $INP_zip, $INP_gzip, $INP_direct, $INP_uselength, $INP_concatenate, $INP_concatenate_norevcomp, $INP_minoverlap, $INP_maxoverlap, $INP_misoverlap, $INP_maxmismatch);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s{,}' => \@INP_reads_files,
# 	'd|data=s' => \$INP_amplicons_file,
	'o|output=s' => \$INP_outfile,
# 	'di|direct' => \$INP_direct,
# 	'l|length' => \$INP_uselength,
	'c|concat' => \$INP_concatenate,
	'cn|concatnorev' => \$INP_concatenate_norevcomp,
	'min=i' => \$INP_minoverlap,
	'max=i' => \$INP_maxoverlap,
	'mis=i' => \$INP_misoverlap,
	'thr|threads=i' => \$INP_threads,
	'z|zip' => \$INP_zip,
	'gz|gzip' => \$INP_gzip,
	'<>' => \&usage,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i (<file1> <file2>|<path>) -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -i (<file1> <file2>|<path>)\n\t\tInput paired-end read files in FASTQ format (compressed or uncompressed) or path.\n";
# 	print "  -d <file>\tCSV file with primer/amplicon data.\n";
	print "  -o <file>\tOutput file name.\n";
	print "  -min <len>\tMinimum overlapping length (default=automatic).\n";
	print "  -max <len>\tMaximum overlapping length (default=automatic).\n";
	print "  -mis <ratio>\tMaximum allowed ratio between the number of mismatched base pairs and the overlap length (default=automatic).\n";
	print "  -c\t\tConcatenate reads (reverse complementary) instead of overlap.\n";
	print "  -cn\t\tConcatenate reads without changing paired-read orientation.\n";
# 	print "  -l\t\tUse marker length info to calculate optimum overlapping parameters (default=automatic).\n";
# 	print "  -di\t\tAnalyze reads only in direct sense.\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
	print "  -gz\t\tCompress results in GZIP format.\n";
	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}

# Prints usage help if no input file is specified
if (!defined($INP_reads_files[0]) || (!defined($INP_reads_files[1]) && !is_multifile($INP_reads_files[0]) && !-d $INP_reads_files[0]) ){
	print "\nERROR: You must specify two paired-end input files or multifile or path.\n\n";
	usage();
	exit;
}
# if (!defined($INP_amplicons_file)){
# 	print "\nERROR: You must specify primer/amplicon CSV file/s.\n\n";
# 	usage();
# 	exit;
# }
# if (defined($INP_concatenate) && ($INP_concatenate<0 || $INP_concatenate>100)){
# 	print "\nERROR: You must specify a minimum percentage of merged reads, below which reads will be concatenated.\n\n";
# 	usage();
# 	exit;
# }

# Creates name for the output file
if (!defined($INP_outfile)){
	if ($INP_reads_files[0] =~ /(.+\/)?(.+?)(_?R1.+)/ || $INP_reads_files[0] =~ /(.+\/)?(.+?)(_1\..+)/ || $INP_reads_files[0] =~ /(.+\/)?(.+?)\./){
		$INP_outfile = $2.".merged";
	} else {
		$INP_outfile = $INP_reads_files[0].".merged";
	}
}

print "\nRunning '$COMMAND_LINE'\n";


# Checks and reads amplicons file
# my ($markerdata,$primers,$sampledata,$tags,$paramsdata,$alleledata)
# = parse_amplicon_file($INP_amplicons_file,['verbose', 'skip samples']); # ,['skip errors']);
# if (scalar @{$primer_tag_headers} > $MAX_PRIMER_TAG_SEQS ){
# 	print "\nERROR: Too many primer and/or tag sequences to be processed, decrease the number of primers and/or samples in the analysis.\n\n";
# 	exit;
# }
# foreach my $paramname (keys %{$paramsdata}){
# 	if ($paramname eq 'min_overlap_length' && !defined($INP_minoverlap)){
# 		$INP_minoverlap = $paramsdata->{$paramname};
# 	} elsif ($paramname eq 'max_overlap_length' && !defined($INP_maxoverlap)){
# 		$INP_maxoverlap = $paramsdata->{$paramname};
# 	} elsif ($paramname eq 'mis_overlap_ratio' && !defined($INP_misoverlap)){
# 		$INP_misoverlap = $paramsdata->{$paramname};
# 	} elsif ($paramname eq 'concatenate_reads' && !defined($INP_concatenate)){
# 		$INP_concatenate = $paramsdata->{$paramname};
# # 	} elsif ($paramname eq 'use_marker_length' && !defined($INP_uselength)){
# # 		$INP_uselength = $paramsdata->{$paramname};
# 	}
# }

print "\nMerging paired-end files.\n";

my $paired_read_files = [ [@INP_reads_files] ];

# Checks if the input file is a set of packed files (.zip or .tar.gz/.tgz)
my $INP_multifile = 0;
my $multifile_tmpdir;
if (is_multifile($INP_reads_files[0]) || (defined($INP_reads_files[1]) && is_multifile($INP_reads_files[1]))) {
	$multifile_tmpdir = "/tmp/".random_file_name();
	$paired_read_files = extract_paired_read_files_from_multifiles(\@INP_reads_files,$multifile_tmpdir);
	$INP_multifile = 1;
} elsif (-d $INP_reads_files[0]) {
	$paired_read_files = extract_paired_read_files_from_path($INP_reads_files[0]);
}
# # If only one is multifile
# } else {
# 	print "\nERROR: '".$multifiles->[$i]."' doesn't contain multiple files.\n\n";
# 	usage();
# 	exit;

my $merged_tmp_folder = "/tmp/".random_file_name();
mkdir($merged_tmp_folder);

my @merged_files;
my $total_merged_seqs = 0;
foreach my $paired_reads_file_pair (@$paired_read_files) {
	my $paired_file1 = $paired_reads_file_pair->[0];
	my $paired_file2 = $paired_reads_file_pair->[1];
	my $file1_name = $paired_file1;
	my $file2_name = $paired_file2;
	my $number_total_seqs = count_seqs_from_fastq($paired_file1);
	my $common_name = $file1_name;
	if ($paired_file1 =~ /(.+\/)?(.+?)(_?R[12].+)/ || $paired_file1 =~ /(.+\/)?(.+?)(_[12]\..+)/){
		$file1_name = $2.$3;
		$common_name = $2;
	}
	if ($paired_file2 =~ /(.+\/)?(.+?)(_?R[12].+)/ || $paired_file2 =~ /(.+\/)?(.+?)(_[12]\..+)/){
		$file2_name = $2.$3;
	}
	# Defines a temporal file to store the merged results
	my $merged_tmp_file = $merged_tmp_folder."/$common_name.fastq";
	# Merges reads with FLASH2 default parameters if others are not specified
	if (!defined($INP_concatenate) && !defined($INP_concatenate_norevcomp)){
		# Set parameters and merge reads with FLASH program
		my $flash_options = '';
		if (defined($INP_minoverlap)) {
			$flash_options .= " -m $INP_minoverlap";
			printf("\tMinimum overlap: %d\n",$INP_minoverlap);
		}
		if (defined($INP_maxoverlap)) {
			$flash_options .= " -M $INP_maxoverlap";
			printf("\tMaximum overlap: %d\n",$INP_maxoverlap);
		}
		if (defined($INP_misoverlap)) {
			$flash_options .= " -x $INP_misoverlap";
			printf("\tMismatch ratio: %.2f\n",$INP_misoverlap);
		}
		if (defined($INP_threads) && $INP_threads>1){
			$flash_options .= " -t $INP_threads";
		}
		printf("\tMerging '%s' and '%s'\n",$file1_name,$file2_name);
		my $merged_file = execute_flash($paired_file1,$paired_file2,$flash_options,$merged_tmp_file);
		my $number_merged_seqs = count_seqs_from_fastq($merged_file);
		$total_merged_seqs += $number_merged_seqs;
		if ($number_merged_seqs>0){
			push(@merged_files,$merged_file);
			printf("\tMerged %d reads from %d\n", $number_merged_seqs, $number_total_seqs);
		} else {
			printf("\tWARNING: Some error occurred, no reads were merged.\n");
			next;
		}
	} else {
		printf("\tConcatenating '%s' and '%s'\n",$file1_name,$file2_name);

		# Opens reads files
		if (is_gzip($paired_file1) || is_zip($paired_file1)) {
			open(PAIREDFILE1, "zcat $paired_file1 |") || die "\nERROR: cannot open compressed '$paired_file1'\n\n";
		} else {
			open(PAIREDFILE1, $paired_file1) || die "\n ERROR: cannot open '$paired_file1'\n\n";
		}
		if (is_gzip($paired_file2) || is_zip($paired_file2)) {
			open(PAIREDFILE2, "zcat $paired_file2 |") || die "\nERROR: cannot open compressed '$paired_file2'\n\n";
		} else {
			open(PAIREDFILE2, $paired_file2) || die "\n ERROR: cannot open '$paired_file2'\n\n";
		}
		my $concatenated_file = $merged_tmp_file;
		open(OUTFILE,">$concatenated_file") || die "\n ERROR: cannot create '$merged_tmp_file'\n\n";
		# Concatenates reads
		my $number_concatenated_seqs = 0;
		my $field_count = 0;
		while (my $line1 = <PAIREDFILE1>) {
			my $line2 = <PAIREDFILE2>;
			if ($field_count) {
				$field_count++;
			}
			if ($field_count == 2){
				chomp $line1;
				if (defined($INP_concatenate)){
					chomp $line2;
					print OUTFILE $line1.iupac_reverse_complementary($line2)."\n";
				} elsif (defined($INP_concatenate_norevcomp)){
					print OUTFILE $line1.$line2;
				}
				$number_concatenated_seqs++;
			} elsif ($field_count == 4){ 
				chomp $line1;
				if (defined($INP_concatenate)){
					chomp $line2;
					print OUTFILE $line1.reverse_sequence($line2)."\n";
				} elsif (defined($INP_concatenate_norevcomp)){
					print OUTFILE $line1.$line2;
				}
				$field_count = 0;
			} elsif ($line1 =~ /^@/){
				print OUTFILE $line1;
				$field_count = 1;
			} else {
				print OUTFILE $line1;
			}
		}
		$total_merged_seqs += $number_concatenated_seqs;
		close PAIREDFILE1;
		close PAIREDFILE2;
		close OUTFILE;
		if ($number_concatenated_seqs == $number_total_seqs){
			push(@merged_files,$concatenated_file);
			print "\tConcatenated $number_concatenated_seqs reads\n";
		} elsif ($number_concatenated_seqs>0) {
			push(@merged_files,$concatenated_file);
			printf("\tWARNING: Some error occurred, only %d reads were concatenated from .\n", $number_concatenated_seqs, $number_total_seqs);
		} else {
			printf("\tWARNING: Some error occurred, no reads were concatenated.\n");
			next;
		}
	}
}
# Removes temporal folder
if (defined($multifile_tmpdir)) {
	system("rm -rf $multifile_tmpdir");
}

# # OBSOLETE:
# 		my (@concatenated_seqs, @concatenated_headers, @concatenated_qualities);
# 		my ($seqs1,$headers1,$qualities1) = read_fastq_file($paired_file1,1);
# 		my ($seqs2,$headers2,$qualities2) = read_fastq_file($paired_file2,1);
# 		if ($#{$headers1} != $#{$headers1}){
# 			printf("\tWARNING: Paired-end files have different number of sequences.\n", $file1_name);
# 			next;
# 		}
# 		for (my $i=0; $i<=$#{$headers1}; $i++) {
# 			push(@concatenated_headers,$headers1->[$i]);
# 			if (defined($INP_concatenate)){
# 				push(@concatenated_seqs,$seqs1->[$i].iupac_reverse_complementary($seqs2->[$i]));
# 			} elsif (defined($INP_concatenate_norevcomp)){
# 				push(@concatenated_seqs,$seqs1->[$i].$seqs2->[$i]);
# 			}
# 			push(@concatenated_qualities,$qualities1->[$i].$qualities2->[$i]);
# 		}
# 		my $number_total_seqs = scalar @$seqs1;
# 		my $number_concatenated_seqs = scalar @concatenated_headers;
# 		$total_merged_seqs += $number_concatenated_seqs;
# 		if ($number_concatenated_seqs == $number_total_seqs){
# 			my $concatenated_file = create_fastq_file(\@concatenated_seqs,\@concatenated_headers,\@concatenated_qualities,"$file1_prefix.fastq");
# 			push(@merged_files,$concatenated_file);
# 			print "\tConcatenated $number_concatenated_seqs reads\n";
# 		} elsif ($number_concatenated_seqs>0) {
# 			my $concatenated_file = create_fastq_file(\@concatenated_seqs,\@concatenated_headers,\@concatenated_qualities,"$file1_prefix.fastq");
# 			push(@merged_files,$concatenated_file);
# 			printf("\tWARNING: Some error occurred, only %d reads were concatenated from .\n", $number_concatenated_seqs, $number_total_seqs);
# 		} else {
# 			printf("\tWARNING: Some error occurred, no reads were concatenated.\n");
# 			next;
# 		}
# 	}


# Saves merged reads into output file
if (@merged_files) {
	my $type = 'merged';
	if (defined($INP_concatenate) || defined($INP_concatenate_norevcomp)){
		$type = 'concatenated';
	}
	if (defined($INP_gzip) && scalar @merged_files > 1){
		if (-f "$INP_outfile.tar.gz"){
			`rm $INP_outfile.tar.gz`;
		}
		system("tar -cvzf $INP_outfile.tar.gz ".join(" ",@merged_files)." --remove-files 1>&- 2>&-");
		printf("\nSaved %d files with %d %s sequences into '%s'.\n", scalar @merged_files, $total_merged_seqs, $type, "$INP_outfile.tar.gz");
	} elsif (defined($INP_gzip)){
		system("mv ".$merged_files[0]." $INP_outfile.fq");
		system("gzip -f $INP_outfile.fq");
		printf("\nSaved %d %s sequences into '%s'.\n", $total_merged_seqs, $type, "$INP_outfile.fq.gz");
	} elsif (scalar @merged_files > 1){
		if (-f "$INP_outfile.zip"){
			`rm $INP_outfile.zip`;
		}
		system("zip -jqm $INP_outfile.zip ".join(" ",@merged_files));
		printf("\nSaved %d files with %d %s sequences into '%s'.\n", scalar @merged_files, $total_merged_seqs, $type, "$INP_outfile.zip");
	} elsif (defined($INP_zip)){
		if (-f "$INP_outfile.fq.zip"){
			`rm $INP_outfile.fq.zip`;
		}
		system("zip -jqm $INP_outfile.fq.zip ".join(" ",@merged_files));
		printf("\nSaved %d %s sequences into '%s'.\n", $total_merged_seqs, $type, "$INP_outfile.fq.zip");
	} else {
		system("mv ".$merged_files[0]." $INP_outfile.fq");
		printf("\nSaved %d %s sequences into '%s'.\n", $total_merged_seqs, $type, "$INP_outfile.fq");
	}
} else {
	print "\nThere was some error in the merging process and no sequences were retrieved.\n";
}

# Removes temporal merged files
if (defined($merged_tmp_folder)) {
	system("rm -rf $merged_tmp_folder");
}

print "\n";

exit;



# ## OBSOLETE:
# my ($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads);
# # Proccess each paired-read file (R1 and R2 usually)
# for (my $i=0; $i<=1; $i++) {
# 	# Check and read some random sequences from read files
# 	($reads_file_format,$read_seqs->[$i],$read_headers->[$i],$read_qualities->[$i],$total_reads->[$i])
# 	= parse_sequence_file($INP_reads_files[$i],undef,['verbose','qualities']);
# 	if ($reads_file_format ne 'fastq'){
# 		print "\nERROR: '".$INP_reads_files[$i]."' file is not FASTQ format.\n\n";
# 		usage();
# 		exit;
# 	}
# }
# # Creates files with all sequences
# my @raw_seqs_files;
# $raw_seqs_files[0] = write_to_file("/tmp/".random_file_name(),join("\n",@{$read_seqs->[0]}));
# $raw_seqs_files[1] = write_to_file("/tmp/".random_file_name(),join("\n",@{$read_seqs->[1]}));
# 
# # Defines match amplicon options
# my @match_options;
# if (defined($INP_direct)) {
# 	push(@match_options, 'direct');
# }
# 
# 
# # Merges marker reads separately, to optimize the merging parameters to the marker length
# my (%pos_seqs_merged,@merged_files);
# foreach my $marker_name (@$primers){
# 
# 	print "\nProcessing '$marker_name' reads.\n";
# 
# 	# Extracts reads containing the amplicons
# 	my $pos_seqs_matched;
# 		
# 	# Changes $markerdata to look for only 1 primer in each file
# 	my $markerdata_;
# 	my $primers_ = [ "$marker_name file 1", "$marker_name file 2" ];
# 	foreach my $primer (@{$markerdata->{$marker_name}{'primer_f'}}, @{$markerdata->{$marker_name}{'primer_rc'}}) {
# 		push(@{$markerdata_->[0]{"$marker_name file 1"}{'primer_f'}}, $primer);
# 		push(@{$markerdata_->[1]{"$marker_name file 2"}{'primer_f'}}, iupac_reverse_complementary($primer));
# 	}
# # 	$markerdata_->[0]{$marker_name}{'primer_rc'} = [''];
# # 	$markerdata_->[1]{$marker_name}{'primer_f'} = [iupac_reverse_complementary($markerdata->{$marker_name}{'primer_f'}), iupac_reverse_complementary($markerdata->{$marker_name}{'primer_r'})];
# # 	$markerdata_->[1]{$marker_name}{'primer_f'} = [''];
# 	# Proccess each paired-read file (R1 and R2)
# 	for (my $i=0; $i<=1; $i++) {
# 		# Parses reads to find matching primer sequences (only primers, not tags)
# # 		# Threads doesn't increase speed in this step because usually there are few markers
# 		if (defined($INP_threads) && $INP_threads>1){
# 			$pos_seqs_matched->[$i]
# 			= find_amplicon_reads_with_threads($raw_seqs_files[$i],$markerdata_->[$i],undef,[$primers_->[$i]],undef,undef,\@match_options,$INP_threads);
# 		} else {
# 			$pos_seqs_matched->[$i]
# 			= find_amplicon_reads($raw_seqs_files[$i],$markerdata_->[$i],undef,[$primers_->[$i]],undef,undef,\@match_options);
# 		}
# 	}
# 
# 	# Recovers matched reads
# 	my ($read_headers_matched, $read_seqs_matched, $read_qualities_matched, @read_fastq_files_matched);
# 	for (my $pos_seq=1; $pos_seq<=$#{$read_headers->[0]}+1; $pos_seq++){
# 		# If the reads have not been merged before and both reads match the primers
# 		if (!defined($pos_seqs_merged{$pos_seq}) && defined($pos_seqs_matched->[0]{$pos_seq}) && defined($pos_seqs_matched->[1]{$pos_seq})){
# 			$pos_seqs_merged{$pos_seq} = 1;
# 			# Proccess each paired-read file (R1 and R2 usually)
# 			for (my $i=0; $i<=1; $i++) {
# 				push(@{$read_headers_matched->[$i]}, $read_headers->[$i][$pos_seq-1]);
# 				push(@{$read_seqs_matched->[$i]}, $read_seqs->[$i][$pos_seq-1]);
# 				if (defined($read_qualities) && @{$read_qualities->[$i]}){
# 					push(@{$read_qualities_matched->[$i]}, $read_qualities->[$i][$pos_seq-1]);
# 				}
# 			}
# 		}
# 	}
# 
# 	if (defined($read_headers_matched->[0])) {
# 	
# 		my $number_matched_seqs = scalar @{$read_headers_matched->[0]};
# 
# 		# Creates a FASTQ file with the matched reads
# 		for (my $i=0; $i<=1; $i++) {
# 			if (defined($read_seqs_matched->[$i]) && @{$read_seqs_matched->[$i]}){
# 				$read_fastq_files_matched[$i] = create_fastq_file($read_seqs_matched->[$i],$read_headers_matched->[$i],$read_qualities_matched->[$i],undef,1);
# 			}
# 		}
# 
# 		# Calculates optimum merging parameters and read lengths (percentile 90 of read lengths)
# 		$markerdata_ = { $marker_name => $markerdata->{$marker_name} }; # Only data from the processed marker
# 		my ($min_overlap,$max_overlap,$read_lengths) = calculate_merging_params($markerdata_,$sampledata,$read_seqs_matched);
# 		my $min_read_length = min(@$read_lengths);
# 
# 		# Sets default merging parameters in case there is any error or we don't use marker length to calculate them
# 		if (!defined($INP_uselength) || !defined($min_overlap) || !defined($max_overlap) || $max_overlap < $min_overlap || $min_overlap<$MIN_READ_OVERLAP){
# 			$min_overlap = $MIN_READ_OVERLAP;
# 			$max_overlap = $min_read_length;
# 		}
# 
# 		# If defined minimum or maximum overlapping lenghts, use them instead
# 		if (defined($INP_minoverlap)) {
# 			$min_overlap = $INP_minoverlap;
# 		}
# 		if (defined($INP_maxoverlap)) {
# 			$max_overlap = $INP_maxoverlap;
# 		}
# 		# Defines default 'Maximum allowed ratio between the number of mismatched base pairs and the overlap length'
# 		# FLASH default: 0.25.
# 		my $mis_overlap;
# 		if (defined($INP_misoverlap)) {
# 			$mis_overlap = $INP_misoverlap;
# 		} else {
# 			$mis_overlap = 0.25;
# 		}
# # 		my $mis_overlap;
# # 		if (defined($INP_misoverlap)) {
# # 			$mis_overlap = $INP_misoverlap;
# # 		} else {
# # 			$mis_overlap = sprintf("%.2f", 0.02*$max_overlap/100);
# # 		}
# 
# # 		my $flash_options = " -x 0";
# # 		if (defined($min_overlap) && defined($max_overlap) && $max_overlap >= $min_overlap && $min_overlap>=$MIN_READ_OVERLAP){
# # 			print "\tMerging $number_matched_seqs reads (min_overlap=$min_overlap, max_overlap=$max_overlap).\n";
# # 			$flash_options .= " -m $min_overlap -M $max_overlap";
# # 		} else {
# # 			# print "\tERROR: overlapping lengths couldn't be calculated.\n";
# # 			print "\tMerging $number_matched_seqs reads with default parameters (min_overlap=$MIN_READ_OVERLAP, max_overlap=$min_read_length).\n";
# # 			$flash_options .= " -m $MIN_READ_OVERLAP -M $min_read_length";
# # 		}
# # 		} elsif (!defined($min_overlap) || !defined($max_overlap)) {
# # 			print "\nERROR: marker '$marker_name' overlapping lengths couldn't be calculated, merging reads with default FLASH parameters.\n\n";
# # 		} elsif ($min_overlap<$MIN_READ_OVERLAP) {
# # 			print "\nERROR: marker '$marker_name' minimum overlapping length ($min_overlap) cannot be lower than $MIN_READ_OVERLAP, merging reads with default FLASH parameters.\n\n";
# # 		} elsif ($max_overlap < $min_overlap) {
# # 			print "\nERROR: marker '$marker_name' maximum overlapping length ($max_overlap) cannot be lower than minimum overlapping length ($min_overlap), merging reads with default FLASH parameters.\n\n";
# # 		}
# 		
# 		my ($merged_file,$number_merged_seqs);
# 		if (!defined($INP_concatenate) || $INP_concatenate<100){
# 			# Set parameters and merge reads with FLASH program
# 			printf("\tRead lengths: %d+%d\n\tMinimum overlap: %d\n\tMaximum overlap: %d\n\tMismatch ratio: %.2f\n",@$read_lengths,$min_overlap,$max_overlap,$mis_overlap);
# 			printf("\tMerging %d reads.\n",$number_matched_seqs);
# 			my $flash_options = " -m $min_overlap -M $max_overlap -x $mis_overlap";
# 			if (defined($INP_threads) && $INP_threads>1){
# 				$flash_options .= " -t $INP_threads";
# 			}
# 			$merged_file = execute_flash($read_fastq_files_matched[0],$read_fastq_files_matched[1],$flash_options);
# 			system("rm ".$read_fastq_files_matched[0]." ".$read_fastq_files_matched[1]);
# 			$number_merged_seqs = count_seqs_from_fastq($merged_file);
# 		}
# 		if (defined($INP_concatenate) && $number_merged_seqs < $INP_concatenate/100*$number_matched_seqs){
# 			if (defined($merged_file)) {
# 				system("rm $merged_file");
# 				if ($number_merged_seqs == 0){
# 					print "\tNo marker reads overlap.\n";
# 				} else {
# 					printf("\tOnly %d marker reads overlap (%.2f%%).\n", $number_merged_seqs, 100*$number_merged_seqs/$number_matched_seqs);
# 				}
# 			}
# 			# my $merged_reads = read_fastq_file_hash($merged_file);
# 			my (@concatenated_seqs, @concatenated_headers, @concatenated_qualities);
# 			for (my $i=0; $i<=$#{$read_headers_matched->[0]}; $i++) {
# 				$read_headers_matched->[0][$i] =~ s/\/1$|\/2$//;
# 				#if (!defined($merged_reads->{$read_headers_matched->[0][$i]})) {
# 				push(@concatenated_headers,$read_headers_matched->[0][$i]);
# 				push(@concatenated_seqs,$read_seqs_matched->[0][$i].$read_seqs_matched->[1][$i]);
# 				push(@concatenated_qualities,$read_qualities_matched->[0][$i].$read_qualities_matched->[1][$i]);
# 				#}
# 			}
# 			my $number_concatenated_seqs = scalar @concatenated_headers;
# 			#if ($number_merged_seqs + $number_concatenated_seqs == $number_matched_seqs){
# 			print "\tConcatenated $number_concatenated_seqs reads.\n";
# 			my $concatenated_file = create_fastq_file(\@concatenated_seqs,\@concatenated_headers,\@concatenated_qualities);
# 			push(@merged_files,$concatenated_file);
# 			#}
# 		} else {
# 			printf("\tMerged %d reads.\n", $number_merged_seqs);
# 			push(@merged_files,$merged_file);
# 		}
# 	} else {
# 		print "\tERROR: no reads contain the marker '$marker_name'.\n";
# 		next;
# 	}
# }
# 
# # # Merge remaining of the reads
# # if (scalar @{$read_headers->[0]} != scalar keys %pos_seqs_merged){
# # 	my ($read_headers_matched, $read_seqs_matched, $read_qualities_matched, @read_fastq_files_matched);
# # 	for (my $pos_seq=1; $pos_seq<=$#{$read_headers->[0]}+1; $pos_seq++){
# # 		if (!defined($pos_seqs_merged{$pos_seq})){
# # 			$pos_seqs_merged{$pos_seq} = 1;
# # 			# Proccess each paired-read file (R1 and R2 usually)
# # 			for (my $i=0; $i<=1; $i++) {
# # 				push(@{$read_headers_matched->[$i]}, $read_headers->[$i][$pos_seq-1]);
# # 				push(@{$read_seqs_matched->[$i]}, $read_seqs->[$i][$pos_seq-1]);
# # 				if (defined($read_qualities) && @{$read_qualities->[$i]}){
# # 					push(@{$read_qualities_matched->[$i]}, $read_qualities->[$i][$pos_seq-1]);
# # 				}
# # 			}
# # 		}
# # 	}
# # 	# Creates a FASTQ file with the matched reads
# # 	for (my $i=0; $i<=1; $i++) {
# # 		if (defined($read_seqs_matched->[$i]) && @{$read_seqs_matched->[$i]}){
# # 			$read_fastq_files_matched[$i] = create_fastq_file($read_seqs_matched->[$i],$read_headers_matched->[$i],$read_qualities_matched->[$i],undef,1);
# # 		}
# # 	}
# # 	if ($#read_fastq_files_matched == 1) {
# # 		print "\nMerging remaining reads.\n\n";
# # 		# Running FLASH to merge the reads
# # 		my $flash_options = "-x 0";
# # 		if (defined($INP_threads) && $INP_threads>1){
# # 			$flash_options .= " -t $INP_threads";
# # 		}
# # 		my $merged_file = execute_flash($read_fastq_files_matched[0],$read_fastq_files_matched[1],$flash_options);
# # 		push(@merged_files,$merged_file);
# # 		system("rm ".$read_fastq_files_matched[0]." ".$read_fastq_files_matched[1]);
# # 	}
# # }
# 
# # Removes temporal files with reads
# system("rm ".$raw_seqs_files[0]." ".$raw_seqs_files[1]);
# 
# # Writes merged reads into a unique file
# if (@merged_files) {
# 
# 	system("cat ".join(" ",@merged_files)." > $outfile");
# 	system("rm ".join(" ",@merged_files));
# 
# 	my $number_merged_seqs = count_seqs_from_fastq($outfile);
# 
# 	if ($number_merged_seqs > 0){
# 		if (defined($INP_zip) && defined($outfile)){
# 			`zip -jqm $outfile.zip $outfile` ;
# 			printf("\nSaved %d merged sequences into '%s'.\n", $number_merged_seqs, "$outfile.zip");
# 		} elsif (defined($INP_gzip) && defined($outfile)){
# 			`gzip -f $outfile` ;
# 			printf("\nSaved %d merged sequences into '%s'.\n", $number_merged_seqs, "$outfile.gz");
# 		} elsif (defined($outfile)){
# 			printf("\nSaved %d merged sequences into '%s'.\n", $number_merged_seqs, $outfile);
# 		}
# 	} else {
# 		`rm $outfile`;
# 		print "\nThere was some error in the merging process and no sequences were retrieved.\n\n";
# 	}
# 
# } else {
# 	print "\nThere was some error in the merging process and no sequences were retrieved.\n\n";
# }
# 
# print "\n";
# 
# exit;
# 
# ################################################################################
# 
# # Calculates the optimum parameters for merging reads of an amplicon sequencing experiment
# sub calculate_merging_params {
# 
# 	my ($markerdata,$sampledata,$reads) = @_;
# 	
# 	# Max. read lengths
# 	my (@read_lengths, @max_read_lengths);
# 	for (my $i=0; $i<=1; $i++) {
# 		if (scalar @{$reads->[$i]} > 1) {
# 			my $read_lengths_ = Statistics::Descriptive::Full->new();
# 			$read_lengths_->add_data( map length($_), @{$reads->[$i]} );
# 			$read_lengths[$i] = $read_lengths_->percentile('90');
# 			$max_read_lengths[$i] = $read_lengths_->max();
# 		} else {
# 			$read_lengths[$i] = $max_read_lengths[$i] = length($reads->[$i][0]);
# 		}
# 	}
# 
# 	my ($min_overlap,$max_overlap);
# 
# 	if (defined($markerdata) && defined($sampledata)) {
# 
# 		# Max. primer lengths
# 		my (@primer_fwd_lengths, @primer_rev_lengths, @marker_lengths);
# 		foreach my $marker (keys $markerdata){
# 			push(@primer_fwd_lengths, map length($_), @{$markerdata->{$marker}{'primer_f'}} );
# 			push(@primer_rev_lengths, map length($_), @{$markerdata->{$marker}{'primer_r'}} );
# 			if (defined($markerdata->{$marker}{'length'})){
# 				push(@marker_lengths, @{$markerdata->{$marker}{'length'}});
# 			}
# 		}
# 
# 		# Max. tag lengths
# 		my (@tag_fwd_lengths, @tag_rev_lengths);
# 		foreach my $sample (keys $sampledata){
# 			push(@tag_fwd_lengths, length($sampledata->{$sample}{'tag_f'}) );
# 			push(@tag_rev_lengths, length($sampledata->{$sample}{'tag_r'}) );
# 		}
# 		
# 		my $read_length = $read_lengths[0]+$read_lengths[1];
# 		
# 		# If marker lengths are not defined returns undef overlapping values
# 		if (@marker_lengths) {
# 			my $max_read_length = $max_read_lengths[0]+$max_read_lengths[1];
# 			my $max_primer_tag_length = max(@tag_fwd_lengths)+max(@tag_rev_lengths)+max(@primer_fwd_lengths)+max(@primer_rev_lengths);
# 			$min_overlap = $read_length - $max_primer_tag_length - max(@marker_lengths);
# 			$max_overlap = $max_read_length - $max_primer_tag_length - min(@marker_lengths);
# 		}
# 	}
# 	
# 	return ($min_overlap,$max_overlap,\@read_lengths);
# }
# 
# ################################################################################

