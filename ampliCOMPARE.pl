#!/usr/bin/perl -w
#
################################################################
#
# Name: ampliCOMPARE.pl
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
#   Compares two AmpliSAS/AmpliCHECK Excel result files
#
# Examples:
# perl ampliCOMPARE.pl -1 results_file1.xlsx -2 results_file2.xlsx -o comparison.xlsx

my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Performs a comparison between two genotyping Excel files.\nExcel files to combine must have been created by AmpliSAS, AmpliCHECK or AmpliLEGACY tools.";


# Modules are in folder 'lib' in the path of the script
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;

use warnings;
no warnings ('uninitialized', 'substr');

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my ($INP_file1,$INP_file2,$INP_demultiplexed,$INP_clustered,$INP_minreads,$INP_allele_file,$INP_output,$INP_verbose);

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

GetOptions(
	'h|help|?' =>  \&usage,
	'1=s' => \$INP_file1,
	'2=s' => \$INP_file2,
	'a|alleles=s' => \$INP_allele_file,
	'd|demul=s' => \$INP_demultiplexed,
	'c|clus=s' => \$INP_clustered,
	'o|output=s' => \$INP_output,
	'm=f' => \$INP_minreads,
	'v|verbose' => \$INP_verbose,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -1 <file1> -2 <file2> [options]\n";
	print "\nOptions:\n";
	print "  -1 <file>\tFirst genotyping Excel file\n";
	print "  -2 <file>\tSecond genotyping Excel file\n";
	print "  -o <file>\tOutput file name\n";
	print "  -a <file>\tFASTA file with allele names and sequences.\n";
	print "  -m <number>\tMinimum reads per sample\n";
	print "  -d\t\tOnly de-multiplexed Excel file\n";
	print "  -c\t\tOnly clustered Excel file\n";
	print "  -v\t\tVerbose output, prints additional information about comparison issues.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}

# Prints usage help if no input file is specified
if (!defined($INP_file1) || !defined($INP_file2)){
	print "\nERROR: You must specify 2 valid input Excel files to compare.\n\n";
	usage();
	exit;
}
if (!defined($INP_output)){
	$INP_output = 'comparison.xlsx';
} elsif ($INP_output !~ /\.xlsx$/) {
	$INP_output .= '.xlsx';
}

# Check and read alleles file (optional)
# Amplicon data file has preference over alleles in FASTA file
my $alleledata;
if (defined($INP_allele_file) && -e $INP_allele_file){
	print "\nReading allele sequences from '$INP_allele_file'.\n";
	$alleledata = read_allele_file($INP_allele_file);
}

my @files = ($INP_file1, $INP_file2);

my $verbose_output = '';
my $results;

print "\nRunning '$COMMAND_LINE'\n";
print "\n";
printf("Reading File '%s'.\n", $INP_file1);
$results->{$INP_file1} = read_amplisas_file_results($INP_file1);
print "\n";
printf("Reading File '%s'.\n", $INP_file2);
$results->{$INP_file2} = read_amplisas_file_results($INP_file2);
print "\n";
if (defined($INP_demultiplexed)){
	printf("Reading File '%s'.\n", $INP_demultiplexed);
	$results->{$INP_demultiplexed} = read_amplisas_file_results($INP_demultiplexed);
	print "\n";
}
if (defined($INP_clustered)){
	printf("Reading File '%s'.\n", $INP_clustered);
	$results->{$INP_clustered} = read_amplisas_file_results($INP_clustered);
	print "\n";
}

my $workbook  = Excel::Writer::XLSX->new($INP_output);
$workbook->set_properties(
	title    => "AmpliSAS results comparison",
	author   => "Alvaro Sebastian",
	comments => "AmpliSAS results comparison",
	company  => "SixthResearcher",
);
$workbook->compatibility_mode();
my $bold = $workbook->add_format( bold => 1 );
my $red = $workbook->add_format(bg_color => 'red');
my $green = $workbook->add_format(bg_color => 'green');
my $blue = $workbook->add_format(bg_color => 'blue');
my $yellow = $workbook->add_format(bg_color => 'yellow');
my $magenta = $workbook->add_format(bg_color => 'magenta');
my $cyan = $workbook->add_format(bg_color => 'cyan');

my @unique_markers = intersect(([keys %{$results->{$INP_file1}}], [keys %{$results->{$INP_file2}}]));

foreach my $marker_name (@unique_markers){
	if (!defined($results->{$INP_file1}{$marker_name})){
		$verbose_output .= sprintf("Marker '%s' is not present in file '%s'.\n", $marker_name, $INP_file1);
		next;
	} elsif (!defined($results->{$INP_file2}{$marker_name})){
		$verbose_output .= sprintf("Marker '%s' is not present in file '%s'.\n", $marker_name, $INP_file2);
		next;
	}

	my %stats;
	($stats{'total_samples'}, $stats{'file1_samples'}, $stats{'file2_samples'}) = (0,0,0);
	($stats{'total_seqs'}, $stats{'file1_seqs'}, $stats{'file2_seqs'}) = (0,0,0);
	($stats{'total_compared_samples'}, $stats{'file1_excluded_samples'}, $stats{'file2_excluded_samples'}) = (0,0,0);
	($stats{'total_compared_seqs'}, $stats{'file1_missing_seqs'}, $stats{'file2_missing_seqs'}) = (0,0,0);
	($stats{'common_assigns'},$stats{'file1_missing_assigns'},$stats{'file2_missing_assigns'}) = (0,0,0);

	my @samples1 = @{$results->{$INP_file1}{$marker_name}{'samples'}};
	my @samples2 = @{$results->{$INP_file2}{$marker_name}{'samples'}};
	my @unique_samples = unique((@samples1, @samples2));
# 	my @unique_samples = intersect((\@samples1, \@samples2));
	my @final_samples;

	my @seq_md5s1 = @{$results->{$INP_file1}{$marker_name}{'seq_md5s'}};
	my @seq_md5s2 = @{$results->{$INP_file2}{$marker_name}{'seq_md5s'}};
	my @unique_seq_md5s = unique((@seq_md5s1, @seq_md5s2));

	$stats{'file1_samples'} = scalar @samples1;
	$stats{'file2_samples'} = scalar @samples2;
	$stats{'total_samples'} = scalar @unique_samples;
	$stats{'file1_seqs'} = scalar @seq_md5s1;
	$stats{'file2_seqs'} = scalar @seq_md5s2;
	$stats{'total_seqs'} = scalar @unique_seq_md5s;

	my (@depth_amplicon, @depth_alleles, @count_alleles);
	for (my $i=0; $i<=$#unique_samples; $i++){
		my $sample = $unique_samples[$i];
		my $final_sample = 1;
		if (!in_array(\@samples1,$sample)){
			$verbose_output .= sprintf("Sample '%s' is not present for marker '%s' in file '%s'.\n", $sample, $marker_name, $INP_file1);
			$final_sample = 0;
			$stats{'file1_excluded_samples'}++;
		} elsif (defined($INP_minreads) && $results->{$INP_file1}{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'} < $INP_minreads){
			$verbose_output .= sprintf("Sample '%s' has less than %d sequences for marker '%s' in file '%s'.\n", $sample, $INP_minreads, $marker_name, $INP_file1);
			$final_sample = 0;
			$stats{'file1_excluded_samples'}++;
		}
		if (!in_array(\@samples2,$sample)){
			$verbose_output .= sprintf("Sample '%s' is not present for marker '%s' in file '%s'.\n", $sample, $marker_name, $INP_file2);
			$final_sample = 0;
			$stats{'file2_excluded_samples'}++;
		} elsif (defined($INP_minreads) && $results->{$INP_file2}{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'} < $INP_minreads){
			$verbose_output .= sprintf("Sample '%s' has less than %d sequences for marker '%s' in file '%s'.\n", $sample, $INP_minreads, $marker_name, $INP_file2);
			$final_sample = 0;
			$stats{'file2_excluded_samples'}++;
		}
		if ($final_sample) {
			push(@final_samples,$sample);
			if (defined($results->{$INP_file1}{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'}) && defined($results->{$INP_file2}{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'})){
				push(@depth_amplicon,sprintf("%d/%d",$results->{$INP_file1}{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'},$results->{$INP_file2}{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'}));
			}
			if (defined($results->{$INP_file1}{$marker_name}{'sample_data'}{$sample}{'depth_alleles'}) && defined($results->{$INP_file2}{$marker_name}{'sample_data'}{$sample}{'depth_alleles'})){
				push(@depth_alleles,sprintf("%d/%d",$results->{$INP_file1}{$marker_name}{'sample_data'}{$sample}{'depth_alleles'},$results->{$INP_file2}{$marker_name}{'sample_data'}{$sample}{'depth_alleles'}));
			}
			if (defined($results->{$INP_file1}{$marker_name}{'sample_data'}{$sample}{'count_alleles'}) && defined($results->{$INP_file2}{$marker_name}{'sample_data'}{$sample}{'count_alleles'})){
				push(@count_alleles,sprintf("%d/%d",$results->{$INP_file1}{$marker_name}{'sample_data'}{$sample}{'count_alleles'},$results->{$INP_file2}{$marker_name}{'sample_data'}{$sample}{'count_alleles'}));
			}
# 			push(@depth_amplicon,$results->{$INP_file1}{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'});
# 			push(@depth_alleles,$results->{$INP_file1}{$marker_name}{'sample_data'}{$sample}{'depth_alleles'});
# 			push(@count_alleles,$results->{$INP_file1}{$marker_name}{'sample_data'}{$sample}{'count_alleles'});
			$stats{'total_compared_samples'}++;
		} else {
			$stats{'total_missing_samples'}++;
		}

	}

	my $worksheet = $workbook->add_worksheet("$marker_name");
	$worksheet->set_column('F:F', undef, $bold);
	my $ws_row = 0;
	my @seq_data_headers = ('SEQUENCE', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'NAME');
	if (@depth_amplicon){
		unshift(@depth_amplicon,'DEPTH_AMPLICON');
		$worksheet->write_row($ws_row, $#seq_data_headers, \@depth_amplicon); $ws_row++;
	}
	if (@depth_alleles){
		unshift(@depth_alleles,'DEPTH_ALLELES');
		$worksheet->write_row($ws_row, $#seq_data_headers, \@depth_alleles); $ws_row++;
	}
	if (@count_alleles){
		unshift(@count_alleles,'COUNT_ALLELES');
		$worksheet->write_row($ws_row, $#seq_data_headers, \@count_alleles); $ws_row++;
	}
	$worksheet->write_row($ws_row, 0, ['SEQUENCE', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'NAME', @final_samples], $bold); $ws_row++;

	# Remove sequences that are not present in the final samples
	my %md5_to_sequence;
	for (my $i=0; $i<=$#unique_seq_md5s; $i++){
		my $seq_md5 = $unique_seq_md5s[$i];
		my @sample_assignments1 = keys %{$results->{$INP_file1}{$marker_name}{'assignments'}{$seq_md5}};
		my @sample_assignments2 = keys %{$results->{$INP_file2}{$marker_name}{'assignments'}{$seq_md5}};
		my $is_seq_in_samples = 0;
		foreach my $sample_assignment (unique((@sample_assignments1, @sample_assignments2))){
			if (in_array(\@final_samples,$sample_assignment)){
				$is_seq_in_samples = 1;
				last;
			}
		}
		if (!$is_seq_in_samples) {
			$verbose_output .= sprintf("Sequence '%s' is only present in removed samples.\n", $seq_md5);
			splice(@unique_seq_md5s,$i,1);
			$i--;
		} else {
			$stats{'total_compared_seqs'}++;
			if (defined($results->{$INP_file1}{$marker_name}{'seq_data'}{$seq_md5}{'sequence'})){
				$md5_to_sequence{$seq_md5} = $results->{$INP_file1}{$marker_name}{'seq_data'}{$seq_md5}{'sequence'};
			} else {
				$md5_to_sequence{$seq_md5} = $results->{$INP_file2}{$marker_name}{'seq_data'}{$seq_md5}{'sequence'};
			}
		}

	}

	$verbose_output .= "\n";

	# Align sequences to alleles and assign allele names to sequences
	my %md5_to_name;
	if (defined($alleledata) && %$alleledata){
		print "\nMatching allele sequences.\n";
		# %md5_to_name = %{match_alleles($alleledata,\%md5_to_sequence,\%md5_to_name,{'alignment' => 'dna match minlen 10 revcomp'})};
		%md5_to_name = %{match_alleles($alleledata,\%md5_to_sequence,\%md5_to_name,{'alignment' => 'dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 100', 'aligned' => 1, 'ident' => 1 })};
	}

	foreach my $seq_md5 (@unique_seq_md5s) {

		my %seq_data;
		my $format = undef;
		if (!in_array(\@seq_md5s1,$seq_md5)){
# 			$worksheet->set_row($ws_row, undef, $magenta);
			$format = $magenta;
			%seq_data = %{$results->{$INP_file2}{$marker_name}{'seq_data'}{$seq_md5}};
			$verbose_output .= sprintf("Sequence '%s' is not present in file '%s'.\n", $seq_md5, $INP_file1);
			$stats{'file1_missing_seqs'}++;
		} elsif (!in_array(\@seq_md5s2,$seq_md5)){
# 			$worksheet->set_row($ws_row, undef, $cyan);
			$format = $cyan;
			%seq_data = %{$results->{$INP_file1}{$marker_name}{'seq_data'}{$seq_md5}};
			$verbose_output .= sprintf("Sequence '%s' is not present in file '%s'.\n", $seq_md5, $INP_file2);
			$stats{'file2_missing_seqs'}++;
		} else {
			foreach my $key (keys %{$results->{$INP_file1}{$marker_name}{'seq_data'}{$seq_md5}}) {
				if ($key =~ /name/ && defined($md5_to_name{$seq_md5})){
					$seq_data{$key} = $md5_to_name{$seq_md5};
				} elsif ($key =~ /seq|md5|length|name/ && $results->{$INP_file1}{$marker_name}{'seq_data'}{$seq_md5}{$key} eq $results->{$INP_file2}{$marker_name}{'seq_data'}{$seq_md5}{$key}){
					$seq_data{$key} = $results->{$INP_file1}{$marker_name}{'seq_data'}{$seq_md5}{$key};
				} else {
					$seq_data{$key} = sprintf("%s/%s", $results->{$INP_file1}{$marker_name}{'seq_data'}{$seq_md5}{$key},$results->{$INP_file2}{$marker_name}{'seq_data'}{$seq_md5}{$key});
				}
			}
		}

		my (%seq_data1,%assignments1);
		if (defined($results->{$INP_file1}{$marker_name}{'seq_data'}{$seq_md5})){
			%seq_data1 = %{$results->{$INP_file1}{$marker_name}{'seq_data'}{$seq_md5}};
			%assignments1 = %{$results->{$INP_file1}{$marker_name}{'assignments'}{$seq_md5}};
		}
		my (%seq_data2,%assignments2);
		if (defined($results->{$INP_file2}{$marker_name}{'seq_data'}{$seq_md5})){
			%seq_data2 = %{$results->{$INP_file2}{$marker_name}{'seq_data'}{$seq_md5}};
			%assignments2 = %{$results->{$INP_file2}{$marker_name}{'assignments'}{$seq_md5}};
		}

		my @worksheet_values;
		foreach my $seq_data_header (@seq_data_headers) {
			if (defined($seq_data{lc($seq_data_header)})){
				push(@worksheet_values, $seq_data{lc($seq_data_header)});
			} else {
				push(@worksheet_values, '');
			}
		}
# 		push(@worksheet_values, $seq_md5);
		$worksheet->write_row($ws_row, 0, \@worksheet_values, $format);

		my $ws_col = scalar @worksheet_values;
		foreach my $sample (@final_samples) {

			if (!in_array(\@samples1,$sample)){
				next;
			} elsif (!in_array(\@samples2,$sample)){
				next;
			}

			if (defined($assignments1{$sample}) && defined($assignments2{$sample})){
# 				$worksheet->write($ws_row, $ws_col, $assignments1{$sample});
				$worksheet->write($ws_row, $ws_col, $assignments1{$sample}."/".$assignments2{$sample});
				$stats{'common_assigns'}++;
			} elsif (defined($assignments2{$sample}) && !defined($assignments1{$sample})){
				# Red: If cluster exists, but it has been filtered:
				if (defined($INP_clustered) && defined($results->{$INP_clustered}{$marker_name}{'assignments'}{$seq_md5}{$sample})) {
# 					$worksheet->write($ws_row, $ws_col, $results->{$INP_clustered}{$marker_name}{'assignments'}{$seq_md5}{$sample}, $red);
					$worksheet->write($ws_row, $ws_col, $results->{$INP_clustered}{$marker_name}{'assignments'}{$seq_md5}{$sample}."/".$assignments2{$sample}, $red);
				# Yellow: If sequence exists, but it has been included in another cluster:
				} elsif (defined($INP_demultiplexed) && defined($results->{$INP_demultiplexed}{$marker_name}{'assignments'}{$seq_md5}{$sample})) {
# 					$worksheet->write($ws_row, $ws_col, $results->{$INP_demultiplexed}{$marker_name}{'assignments'}{$seq_md5}{$sample}, $yellow);
					$worksheet->write($ws_row, $ws_col, $results->{$INP_demultiplexed}{$marker_name}{'assignments'}{$seq_md5}{$sample}."/".$assignments2{$sample}, $yellow);
				# Magenta: If sequence doesn't exist in file 1:
				} else {
# 					$worksheet->write($ws_row, $ws_col, 0, $magenta);
					$worksheet->write($ws_row, $ws_col, $assignments2{$sample}, $magenta);
				}
				$stats{'file1_missing_assigns'}++;
			} elsif (defined($assignments1{$sample}) && !defined($assignments2{$sample})){
				$worksheet->write($ws_row, $ws_col, $assignments1{$sample}, $cyan);
				$stats{'file2_missing_assigns'}++;
			} else {
				$worksheet->write($ws_row, $ws_col, '');
			}
			$ws_col++;
		}
		$ws_row++;
	}
	
	$stats{'total_assigns'} = $stats{'common_assigns'}+$stats{'file2_missing_assigns'}+$stats{'file1_missing_assigns'};
	printf("\nMARKER '%s':\n", $marker_name);
	printf("Total unique samples: %d (file1: %d, file2: %d)\n", $stats{'total_samples'}, $stats{'file1_samples'}, $stats{'file2_samples'});
	printf("Total seqs: %d (file1: %d, file2: %d)\n", $stats{'total_seqs'}, $stats{'file1_seqs'}, $stats{'file2_seqs'});
	printf("Compared samples: %d (excluded from file1: %d, from file2: %d)\n", $stats{'total_compared_samples'}, $stats{'file1_excluded_samples'}, $stats{'file2_excluded_samples'});
	printf("Compared seqs: %d (missing in file1: %d, in file2: %d)\n", $stats{'total_compared_seqs'}, $stats{'file1_missing_seqs'}, $stats{'file2_missing_seqs'});
	printf("Total assignments: %d (missing in file1: %d, missing in file2: %d)\n", $stats{'total_assigns'}, $stats{'file1_missing_assigns'}, $stats{'file2_missing_assigns'});

}

$workbook->close();

if (defined($INP_verbose)){
	print $verbose_output;
}

print "\nComparison results written into '$INP_output'.\n\n";

exit;



















