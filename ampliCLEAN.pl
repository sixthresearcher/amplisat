#!/usr/bin/perl -w
#
################################################################
#
# Name: ampliCLEAN.pl
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
#   Filters sequences/reads containing the primers and tags specified in a CSV with amplicon data
#   Generates a FASTA format file with the filtered sequences
#
# Requires as input a FASTA or FASTQ file with sequences/reads and a CSV format file with primer/tag data information
# Example: perl ampliCLEAN.pl -d amplicon_data.csv -i reads.fq.gz -o filtered_reads
#
# If the reads have been already demultiplexed into separate files (one file per sample), they can be packed into a single .zip or .tar.gz file and use it as input
# Example: perl ampliCLEAN.pl -i reads.tar.gz -o filtered_reads
#

my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Cleans FASTQ files, filtering reads not belonging to any amplicon and removing low quality and/or anomalous short/long ones.";

# Modules are in folder 'lib' in the path of the script
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use List::Util qw(shuffle);
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;


my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options
# Maximum number of allowed sequences per amplicon
my $INP_nreads_amplicon; # = 100000;
# Maximum number of allowed sequences
# my $MAX_READS_NUMBER = 7000000;
# Maximum number of allowed primer+tag sequences
# my $MAX_PRIMER_TAG_SEQS = 10000;


my ($INP_amplicons_file, $INP_reads_file, $INP_outfile, $INP_threads, $INP_nreads, $INP_direct, $INP_shuffle, $INP_trim, $INP_minlength, $INP_maxlength, $INP_minqual, $INP_zip, $INP_gzip);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s' => \$INP_reads_file,
	'd|data=s' => \$INP_amplicons_file,
	'o|output=s' => \$INP_outfile,
	'di|direct' => \$INP_direct,
	'n|number=i' => \$INP_nreads,
	'na|nampli=i' => \$INP_nreads_amplicon,
	's|shuffle' => \$INP_shuffle,
	'min=i' => \$INP_minlength, # Removes sequences shorter than threshold
	'max=i' => \$INP_maxlength, # Removes sequences longer than threshold
	'mqual=i' => \$INP_minqual, # Removes sequences with lower mean Phred quality
	'trim=i' => \$INP_trim, # Trims sequences to the desired length
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
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -i <file>\tInput FASTQ or FASTA file (compressed or uncompressed) or set of files packed into a unique .ZIP or .TAR.GZ file.\n";
	print "  -d <file>\tCSV file with primer/amplicon data.\n";
	print "  -o <file>\tOutput file name.\n";
	print "  -mqual <qual>\tMinimum mean Phred quality score (Sanger)\n";
	print "  -min <len>\tMinimum sequence length\n";
	print "  -max <len>\tMaximum sequence length\n";
	print "  -trim <len>\tTrims reads/sequences to the desired length.\n";
	print "  -n <number>\tTotal number of reads/sequences to analyze.\n";
	print "  -na <number>\tNumber of reads/sequences per amplicon to extract.\n";
	print "  -s\t\tShuffle/randomize reads/sequences.\n";
	print "  -di\t\tAnalyze reads only in direct sense.\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -gz\t\tCompress results in GZIP format.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}

my $INP_multifile;
# Checks if a set of demultiplexed files is given as input into a compressed file
if (defined($INP_reads_file) && is_multifile($INP_reads_file)){
	$INP_multifile = $INP_reads_file;
# Prints usage help if no input file is specified
} elsif (!defined($INP_reads_file) || !-f $INP_reads_file){
	print "\nERROR: You must specify a sequence input file.\n\n";
	usage();
	exit;
}
if (!defined($INP_amplicons_file) || !-f $INP_amplicons_file){
	print "\nWARNING: You didn't specify an amplicon data input file.\n\n";
# 	usage();
# 	exit;
}
if (defined($INP_trim) && (!is_numeric($INP_trim) || !$INP_trim)){
	print "\nERROR: A trimming valid length must be specified.\n\n";
	exit;
}
# Creates name for the output file
if (!defined($INP_outfile)){
	if ($INP_reads_file =~ /.+\/(.+?)\./){
		$INP_outfile = $1.".clean";
	} else {
		$INP_outfile = "$INP_reads_file.clean";
	}
}
# Defines compression method
my $compression;
if (defined($INP_zip)){
	$compression = 'zip';
} elsif (defined($INP_gzip)){
	$compression = 'gzip';
}

print "\nRunning '$COMMAND_LINE'\n";

# 2 PRIMERS => MARKER
# 1/2 TAGS => SAMPLE
# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

# Check and read amplicons file
my ($markerdata,$primers,$sampledata,$tags,$paramsdata,$alleledata);
if (defined($INP_amplicons_file)){
	($markerdata,$primers,$sampledata,$tags,$paramsdata,$alleledata)
	= parse_amplicon_file($INP_amplicons_file,['skip errors']);
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

my $total_reads_matched = 0;
my $total_reads_checked = 0;
my %paired_files;

while (my $file = shift @file_list){

	my $file_name = $file;

	if (!-f $file || (!is_fasta($file) && !is_fastq($file))){
		next;
	}

	# If there are paired-end read files in multifile
	my ($file2, $file2_name, $file1_prefix, $file2_prefix);
	if ($file =~ /(.+\/)?(.+?)(_?R1.+)/ || $file =~ /(.+\/)?(.+?)(_1\..+)/){
		$file1_prefix = $2;
		if (defined($INP_multifile)){
			foreach my $file2_ (@file_list){
				if (defined($paired_files{$file2_})){ # Already merged file
					next;
				} elsif ($file2_ =~ /(.+\/)?(.+?)(_?R2.+)/ || $file2_ =~ /(.+\/)?(.+?)(_2\..+)/){
					$file2_prefix = $2;
				}
				if ($file1_prefix eq $file2_prefix){
					$file2 = $file2_name = $file2_;
					$paired_files{$file2}=1;
					last;
				}
			}
		}
	}

# 	# If there are paired-end read files in multifile
# 	my ($file2, $file2_name);
# 	if (defined($INP_multifile)){
# 		my ($match,$pos) = in_array(\@file_list,'/(.+\/)?(.+?)(_?R[12].+)|(.+\/)?(.+?)(_?[12]\..+)/',1);
# 		if ($match){
# 			$file2 = $file2_name = $file_list[$pos->[0]];
# 			splice(@file_list,$pos->[0],1);
# 		}
# 	}

	# Checks and reads sequences file
	my ($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads);
	my ($reads_file_format2,$read_seqs2,$read_headers2,$read_qualities2,$total_reads2);
	if (!defined($INP_multifile)){
		($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads)
		= parse_sequence_file($file,$INP_nreads,['qualities','verbose']);
	} elsif (defined($file2)) {
		my ($read_seqs_,$read_headers_,$read_qualities_);
		my ($read_seqs2_,$read_headers2_,$read_qualities2_);
		# Do not take random reads if there are paired-end read files
		($reads_file_format,$read_seqs_,$read_headers_,$read_qualities_,$total_reads)
		= parse_sequence_file($file,undef,['qualities','verbose']);
		if (!defined($read_seqs_)){
			next;
		}
		if ($file =~ /\/tmp.+\/(.+)/){
			$file_name = $1;
		}
		($reads_file_format2,$read_seqs2_,$read_headers2_,$read_qualities2_,$total_reads2)
		= parse_sequence_file($file2,undef,['qualities']);
		if (!defined($read_seqs2_)){
			next;
		}
		if ($file2 =~ /\/tmp.+\/(.+)/){
			$file2_name = $1;
		}
		# Extracts the same random sequences from both paired-end read files
		if (defined($INP_nreads) || defined($INP_shuffle)){
			my @order = shuffle 0..$#{$read_seqs_};
			if (defined($INP_nreads)){
				@order = splice(@order, 0, $INP_nreads);
			}
			$read_seqs = [map $read_seqs_->[$_], @order];
			$read_headers = [map $read_headers_->[$_], @order];
			$read_seqs2 = [map $read_seqs2_->[$_], @order];
			$read_headers2 = [map $read_headers2_->[$_], @order];
			if (defined($read_qualities_)){
				$read_qualities = [map $read_qualities_->[$_], @order];
				$read_qualities2 = [map $read_qualities2_->[$_], @order];
			}
		} else {
			($read_seqs,$read_headers,$read_qualities) = ($read_seqs_,$read_headers_,$read_qualities_);
			($read_seqs2,$read_headers2,$read_qualities2) = ($read_seqs2_,$read_headers2_,$read_qualities2_);
		}
	} else {
		($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads)
		= parse_sequence_file($file,$INP_nreads,['qualities','verbose']);
		if (!defined($read_seqs)){
			next;
		}
		if ($file =~ /\/tmp.+\/(.+)/){
			$file_name = $1;
		}
	}
	if (defined($INP_minqual) && $reads_file_format ne 'fastq'){
		print "\nERROR: Input file must be in FASTQ format to read quality values.\n\n";
		exit;
	}

	my $pos_seqs_matched;
	if (defined($INP_amplicons_file)){
		# Defines match amplicon options
		my @match_options;
		if (defined($INP_direct)) {
			push(@match_options, 'direct');
		}
		# Extracts reads containing the amplicons
		print "\nDe-multiplexing amplicon sequences from reads.\n";
		# Randomizes/shuffles reads if $INP_nreads_amplicon limit exists
		if (defined($INP_shuffle)) {
			($read_headers,$read_seqs) = shuffle_seqs($read_headers,$read_seqs);
		}
		# Creates a file with one line per sequence
		my $raw_seqs_file = write_to_file("/tmp/".random_file_name(),join("\n",@{$read_seqs}));
		# Parses reads and primers+tags to find matching amplicons
		if (defined($INP_threads) && $INP_threads>1){
			$pos_seqs_matched->[0]
			= find_amplicon_reads_with_threads($raw_seqs_file,$markerdata,$sampledata,$primers,$tags,$INP_nreads_amplicon,\@match_options,$INP_threads);
		} else {
			$pos_seqs_matched->[0]
			= find_amplicon_reads($raw_seqs_file,$markerdata,$sampledata,$primers,$tags,$INP_nreads_amplicon,\@match_options);
		}
		`rm $raw_seqs_file`;
		if (defined($file2)){
			my $raw_seqs_file2 = write_to_file("/tmp/".random_file_name(),join("\n",@{$read_seqs2}));
			if (defined($INP_threads) && $INP_threads>1){
				$pos_seqs_matched->[1]
				= find_amplicon_reads_with_threads($raw_seqs_file2,$markerdata,$sampledata,$primers,$tags,$INP_nreads_amplicon,\@match_options,$INP_threads);
			} else {
				$pos_seqs_matched->[1]
				= find_amplicon_reads($raw_seqs_file2,$markerdata,$sampledata,$primers,$tags,$INP_nreads_amplicon,\@match_options);
			}
			`rm $raw_seqs_file2`;
		}
	} else {
		map $pos_seqs_matched->[0]{$_} = 1, 1..(scalar @$read_headers);
		if (defined($file2)){
			map $pos_seqs_matched->[1]{$_} = 1, 1..(scalar @$read_headers2);
		}
	}

	if (!defined($INP_multifile)){
		print "\nCleaning reads.\n";
	} else {
		if (!defined($file2)){
			printf("\nCleaning %d reads from file '%s'.\n", $total_reads, $file_name);
		} else {
			printf("\nCleaning %d reads from paired-end read files '%s' and '%s'.\n", $total_reads, $file_name, $file2_name);
		}
	}

	
	
	my ($read_headers_matched, $read_seqs_matched, $read_qualities_matched);
	my ($read_headers_matched2, $read_seqs_matched2, $read_qualities_matched2);
	# Recovers matched reads
	foreach my $pos_seq (sort {$a<=>$b} keys %{$pos_seqs_matched->[0]}){
		# Skip if the matched read is not matched in the respective paired-end file
		if (defined($pos_seqs_matched->[1]) && !defined($pos_seqs_matched->[1]{$pos_seq})){
			next;
		}
		my $seq = $read_seqs->[$pos_seq-1];
		my $qual;
		if (defined($read_qualities)){
			$qual = $read_qualities->[$pos_seq-1];
		}
		my $seq2;
		if (defined($file2)){
			$seq2 = $read_seqs2->[$pos_seq-1];
		}
		my $qual2;
		if (defined($file2) && defined($read_qualities)){
			$qual2 = $read_qualities->[$pos_seq-1];
		}
		if (defined($INP_trim)){
			$seq = substr($seq,0,$INP_trim);
			if (defined($qual)){
				$qual = substr($qual,0,$INP_trim);
			}
			if (defined($seq2)){
				$seq2 = substr($seq2,0,$INP_trim);
				if (defined($qual2)){
					$qual2 = substr($qual2,0,$INP_trim);
				}
			}
		}
		if (defined($INP_minlength) || defined($INP_maxlength)){
			my $len = length($seq);
			if (defined($INP_minlength) && $len<$INP_minlength){
				next;
			}
			if (defined($INP_maxlength) && $len>$INP_maxlength){
				next;
			}
			if (defined($seq2)){
				my $len2 = length($seq2);
				if (defined($INP_minlength) && $len2<$INP_minlength){
					next;
				}
				if (defined($INP_maxlength) && $len2>$INP_maxlength){
					next;
				}
			}
		}
		if (defined($INP_minqual)){
			my $mean_qual = mean_phred_ascii_quality($qual);
			if ($mean_qual<$INP_minqual){
				next;
			}
			if (defined($qual2)){
				my $mean_qual2 = mean_phred_ascii_quality($qual2);
				if ($mean_qual2<$INP_minqual){
					next;
				}
			}
		}
		push(@$read_headers_matched, $read_headers->[$pos_seq-1]);
		push(@$read_seqs_matched, $seq);
		push(@$read_qualities_matched, $qual);
		if (defined($file2)){
			push(@$read_headers_matched2, $read_headers2->[$pos_seq-1]);
			push(@$read_seqs_matched2, $read_seqs2->[$pos_seq-1]);
			push(@$read_qualities_matched2, $read_qualities2->[$pos_seq-1]);
		}
	}
	
	if (!defined($read_seqs_matched)){
		print "\tERROR: No sequences passed the filters, check if amplicon data is correct or relax the cleaning parameters.\n";
		next;
	}

	$total_reads_matched += scalar @$read_seqs_matched;
	$total_reads_checked += $total_reads;
	my $cleaned_reads = scalar @$read_seqs - scalar @$read_seqs_matched;

	printf("\t%d reads removed.\n", $cleaned_reads);
	printf("\t%d reads kept.\n", scalar @$read_seqs_matched);

	if (!defined($INP_multifile)){
		my $outfile;
		if (defined($read_seqs_matched) && @$read_seqs_matched && $reads_file_format eq 'fastq'){
			if ($INP_outfile !~ /\.(fastq|fq)$/){
				$INP_outfile = "$INP_outfile.fq";
			}
			$outfile = create_fastq_file($read_seqs_matched,$read_headers_matched,$read_qualities_matched,$INP_outfile,$compression);
		} elsif (defined($read_seqs_matched) && @$read_seqs_matched) {
			if ($INP_outfile !~ /\.(fasta|fa|fna)$/){
				$INP_outfile = "$INP_outfile.fa";
			}
			$outfile = create_fasta_file($read_seqs_matched,$read_headers_matched,$INP_outfile,$compression);
		}
		if (defined($outfile)){
			printf("\nSaved %d sequences from %d into '%s'.\n", scalar @$read_seqs_matched, $total_reads, $outfile);
		} else {
			print "\nThere was some error in the extraction and no sequences were retrieved.\n";
		}
		if (defined($outfile)){
			push(@outfiles,$outfile);
		}
	} else {
		my $outfile;
		`rm $file`;
		if (defined($read_seqs_matched) && @$read_seqs_matched && $reads_file_format eq 'fastq'){
			$outfile = create_fastq_file($read_seqs_matched,$read_headers_matched,$read_qualities_matched,$file);
		} elsif (defined($read_seqs_matched) && @$read_seqs_matched) {
			$outfile = create_fasta_file($read_seqs_matched,$read_headers_matched,$file);
		}
		if (defined($outfile)){
			push(@outfiles,$outfile);
		}
		if (defined($file2)){
			my $outfile2;
			`rm $file2`;
			if (defined($read_seqs_matched2) && @$read_seqs_matched2 && $reads_file_format eq 'fastq'){
				$outfile2 = create_fastq_file($read_seqs_matched2,$read_headers_matched2,$read_qualities_matched2,$file2);
			} elsif (defined($read_seqs_matched2) && @$read_seqs_matched2) {
				$outfile2 = create_fasta_file($read_seqs_matched2,$read_headers_matched2,$file2);
			}
			if (defined($outfile2)){
				push(@outfiles,$outfile2);
			}
		}
	}
}

if (defined($INP_multifile)){

	my $outfile;
	if (@outfiles){
		my $outmultifile;
		if (defined($compression) && $compression eq 'tar'){
			$outmultifile = "$INP_outfile.tar";
		} elsif (defined($compression) && ($compression eq 'gzip' || $compression eq 'tgz')){
			$outmultifile = "$INP_outfile.tar.gz";
		} else {
			$compression = 'zip';
			$outmultifile = "$INP_outfile.zip";
		}
		$outfile = compress(\@outfiles,$compression,$outmultifile);
	}
	if (defined($outfile)){
		printf("\nSaved %d sequences from %d into '%s'.\n", $total_reads_matched, $total_reads_checked, $outfile);
	} else {
		print "\nThere was some error in the extraction and no sequences were retrieved.\n";
	}

	`rm -rf $tmp_dir`;
}

print "\n";

exit;
