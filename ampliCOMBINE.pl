#!/usr/bin/perl -w
#
################################################################
#
# Name: ampliCOMBINE.pl
#
# Version: 1.0
#
# License: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
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
#   Combines several AmpliSAS/AmpliCHECK Excel result files
#
# Examples:
# perl ampliCOMBINE.pl -i results_file1.xlsx results_file2.xlsx -o combined.xlsx

my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Combines the results from two or more genotyping Excel files into a unique one.\nExcel files to combine must have been created by AmpliSAS, AmpliCHECK or AmpliLEGACY tools.";


# Modules are in folder 'lib' in the path of the script
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

use warnings;
no warnings ('uninitialized', 'substr');

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my (@INP_files,$INP_allele_file,$INP_minreads,$INP_concatenate,$INP_output,$INP_verbose,$INP_only_alleles);

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s{,}' => \@INP_files,
	'o|output=s' => \$INP_output,
	'a|alleles=s' => \$INP_allele_file,
	'm=f' => \$INP_minreads,
	'c|concat' => \$INP_concatenate,
	'oa|onlyalleles' => \$INP_only_alleles,
	'v|verbose' => \$INP_verbose,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file1> <file2> [ ... <fileN> ] [options]\n";
	print "\nOptions:\n";
	print "  -i <file1> <file2> [ ... <fileN> ]\n\t\tGenotyping Excel files to combine\n";
	print "  -o <file>\tOutput file name\n";
	print "  -a <file>\tFASTA file with allele names and sequences.\n";
	print "  -m <number>\tMinimum reads per sample\n";
	print "  -c\t\tJust concatenate results in Excel files, without combining samples with the same names.\n";
# 	print " -oa Combine only allele matching variants.\n";
	print "  -v\t\tVerbose output, prints additional information about combination issues.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}

# Prints usage help if no input file is specified
if (!@INP_files || !-e $INP_files[0] || !-e $INP_files[1]){
	print "\nERROR: You must specify at least 2 valid input Excel files to combine.\n\n";
	usage();
	exit;
}
if (!defined($INP_output)){
	$INP_output = 'combination.xlsx';
} elsif ($INP_output !~ /\.xlsx$/) {
	$INP_output .= '.xlsx';
}
# Check and read alleles file (optional)
# Amplicon data file has preference over alleles in FASTA file
my $alleledata;
if (defined($INP_allele_file) && -e $INP_allele_file){
	print "\nReading allele sequences from '$INP_allele_file'.\n";
	$alleledata = read_allele_file($INP_allele_file);
} elsif (defined($INP_only_alleles)){
	print "\nERROR: You must specify a FASTA file with allele sequences.\n\n";
	usage();
	exit;
}

my $verbose_output = '';
my $results;

print "\nRunning '$COMMAND_LINE'\n";
foreach my $file (@INP_files) {
	printf("\nReading File '%s'.\n", $file);
	$results->{$file} = read_amplisas_file_results($file);
}

my $workbook  = Excel::Writer::XLSX->new($INP_output);
$workbook->set_properties(
	title    => "AmpliSAS results combination",
	author   => "Alvaro Sebastian",
	comments => "AmpliSAS results combination",
	company  => "Evolutionary Biology Group, Adam Mickiewicz University",
);
$workbook->compatibility_mode();
my $bold = $workbook->add_format( bold => 1 );
my $red = $workbook->add_format(bg_color => 'red');
my $green = $workbook->add_format(bg_color => 'green');
my $blue = $workbook->add_format(bg_color => 'blue');
my $yellow = $workbook->add_format(bg_color => 'yellow');
my $magenta = $workbook->add_format(bg_color => 'magenta');
my $cyan = $workbook->add_format(bg_color => 'cyan');

my @unique_markers;
foreach my $file (@INP_files) {
	@unique_markers = unique(@unique_markers, keys %{$results->{$file}});
}

foreach my $marker_name (@unique_markers){

	print "\nCombining '$marker_name' data:\n";

	# Obtains all the variants in all the files
	my ($unique_samples,%seq_depths,%final_seq_depths,%final_seq_samples);
	my (%md5_to_sequence,$md5_to_names,%names_to_md5s);
	my $generate_new_names = 0;
	my (%common_samples,%common_seqs);
	foreach my $file (@INP_files) {
		printf("\t'%s': %d samples, %d variants.\n", $file, scalar @{$results->{$file}{$marker_name}{'samples'}}, scalar @{$results->{$file}{$marker_name}{'seq_md5s'}});
		if (defined($results->{$file}{$marker_name})){
			foreach my $md5 (keys %{$results->{$file}{$marker_name}{'seq_data'}}){
				$seq_depths{$md5} += $results->{$file}{$marker_name}{'seq_data'}{$md5}{'depth'};
				my $name = $results->{$file}{$marker_name}{'seq_data'}{$md5}{'name'};
				if (!defined($md5_to_sequence{$md5})){
					$md5_to_sequence{$md5} = $results->{$file}{$marker_name}{'seq_data'}{$md5}{'sequence'};
				} elsif (defined($md5_to_sequence{$md5}) && $md5_to_sequence{$md5} ne $results->{$file}{$marker_name}{'seq_data'}{$md5}{'sequence'}){
					print "\nERROR: Sequence ID '$md5' is not unique.\n";
					print "SEQ1: ".$md5_to_sequence{$md5}."\n";
					print "SEQ2: ".$results->{$file}{$marker_name}{'seq_data'}{$md5}{'sequence'}."\n\n";
					exit;
				} else {
					$common_seqs{$md5}++;
				}
				if (!$generate_new_names && defined($name) && defined($names_to_md5s{$name}) && $names_to_md5s{$name} ne $md5 ){
					$generate_new_names = 1;
				} else {
					$names_to_md5s{$name} = $md5
				}
				if (defined($name) && !in_array($md5_to_names->{$md5}, $name)){
					push(@{$md5_to_names->{$md5}}, $name);
				}
# if ($md5 eq '569965416b819382419dc9d3b3469e17'){
# print '';
# }


			}
			foreach my $sample (@{$results->{$file}{$marker_name}{'samples'}}){
				if (defined($unique_samples->{$sample})){
					$common_samples{$sample}++;
				}
				push(@{$unique_samples->{$sample}},$file);
			}
		}
	}
	if (!%seq_depths) { next; }
# 	printf("\tCommon: %d samples, %d variants.\n", scalar keys %common_samples, scalar keys %common_seqs);


	my @seq_md5s = sort { $seq_depths{$b} <=> $seq_depths{$a} } keys %seq_depths;
	
	# Aligns sequences to alleles and assign allele names to sequences
	my %md5_to_name;
	if (defined($alleledata) && %$alleledata){
		print "\nMatching allele sequences.\n";
		%md5_to_name = %{match_alleles($alleledata,\%md5_to_sequence)};
	}
	my $unique_seq_number = 0;
	foreach my $seq_md5 (@seq_md5s) {
# if ($seq_md5 eq '569965416b819382419dc9d3b3469e17'){
# print '';
# }
		$unique_seq_number++;
		# Skip if the sequence is in alleles file or only alleles mode is activated
		if (defined($md5_to_name{$seq_md5}) || defined($INP_only_alleles)){
			next;
		}
		# Generates new names if different variants have the same name in the files to combine
		# Or if the same variant has different names in the different files
		if ($generate_new_names || !defined($md5_to_names->{$seq_md5}) || $#{$md5_to_names->{$seq_md5}} > 0){
			$md5_to_name{$seq_md5} = sprintf("%s-%07d-C", $marker_name, $unique_seq_number);
		# Uses name from the input files if both coincide
# 		} elsif (defined($md5_to_names->{$seq_md5}) && $#{$md5_to_names->{$seq_md5}} == 0){
		} else {
			$md5_to_name{$seq_md5} = $md5_to_names->{$seq_md5}[0];
		}
	}
	
	# Creates worksheet and writes data headers
	my $worksheet = $workbook->add_worksheet("$marker_name");
	$worksheet->set_column('F:F', undef, $bold);
	my $ws_row = 0;
	my @seq_data_headers = ('SEQUENCE', 'MD5', 'LENGTH', 'DEPTH', 'SAMPLES', 'NAME');
	#$worksheet->write($ws_row, $#seq_data_headers, 'FILE'); $ws_row++;
	$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_AMPLICON'); $ws_row++;
	$worksheet->write($ws_row, $#seq_data_headers, 'DEPTH_ALLELES'); $ws_row++;
	$worksheet->write($ws_row, $#seq_data_headers, 'COUNT_ALLELES'); $ws_row++;
	if (defined($INP_concatenate)){
		$worksheet->write($ws_row, $#seq_data_headers, 'FILE'); $ws_row++;
	}
	$worksheet->write_row($ws_row, 0, \@seq_data_headers, $bold); $ws_row++;

	my $ws_row_first = $ws_row;
	my $ws_row_last = $ws_row_first+$#seq_md5s;
	my $ws_col = $#seq_data_headers;
	my $ws_col_first = $ws_col+1;
	my (%final_samples,%final_seqs);
	foreach my $file (@INP_files) {
		if (!defined($results->{$file}{$marker_name})){
			$verbose_output .= sprintf("Marker '%s' is not present in file '%s'.\n", $marker_name, $file);
			next;
		}
		my @samples = @{$results->{$file}{$marker_name}{'samples'}};
		my %stats;
# 		$stats{'file_samples'} = scalar @samples;
# 		$stats{'file_seqs'} = scalar @seq_md5s;
# 		$stats{'file_excluded_samples'} = 0;
# 		$stats{'file_combined_seqs'} = 0;
		foreach my $sample (@samples){
			if (!defined($unique_samples->{$sample})){
				next;
			}
			$final_samples{$sample}++;
# 			$stats{'total_combined_samples'}++;
			my ($depth_alleles, $count_alleles) = (0,0);
			my (@depth_amplicon, @depths, $seq_depths);
			# If results from samples with the same name from different files must be combined
			if (!defined($INP_concatenate)){
				foreach my $file_ (@{$unique_samples->{$sample}}) {
					if (defined($INP_minreads) && $results->{$file_}{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'} < $INP_minreads){
						$verbose_output .= sprintf("Sample '%s' has less than %d sequences for marker '%s' in file '%s'.\n", $sample, $INP_minreads, $marker_name, $file_);
		# 				$stats{'file_excluded_samples'}++;
						next;
					}
		# 			if ($results->{$file_}{$marker_name}{'sample_data'}{$sample}{'depth_alleles'} < 5000){
		# 				$verbose_output .= sprintf("Sample '%s' has less than %d allele depth for marker '%s' in file '%s'.\n", $sample, 5000, $marker_name, $file_);
		# 				next;
		# 			}
					push(@depth_amplicon, $results->{$file_}{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'});
					$depth_alleles += $results->{$file_}{$marker_name}{'sample_data'}{$sample}{'depth_alleles'};
					$count_alleles += $results->{$file_}{$marker_name}{'sample_data'}{$sample}{'count_alleles'};
					foreach my $seq_md5 (@seq_md5s) {
# if ($seq_md5 eq '569965416b819382419dc9d3b3469e17'){
# print ''; # $seq_depths->{'569965416b819382419dc9d3b3469e17'}
# }
						if (defined($results->{$file_}{$marker_name}{'assignments'}{$seq_md5}{$sample})){
							push(@{$seq_depths->{$seq_md5}}, $results->{$file_}{$marker_name}{'assignments'}{$seq_md5}{$sample});
						}
						$final_seqs{$seq_md5}++;
					}
				}
				delete($unique_samples->{$sample});
			# If samples from different files just must be concatenated
			} else {
				push(@depth_amplicon, $results->{$file}{$marker_name}{'sample_data'}{$sample}{'depth_amplicon'});
				$depth_alleles = $results->{$file}{$marker_name}{'sample_data'}{$sample}{'depth_alleles'};
				$count_alleles = $results->{$file}{$marker_name}{'sample_data'}{$sample}{'count_alleles'};
				foreach my $seq_md5 (@seq_md5s) {
					if (defined($results->{$file}{$marker_name}{'assignments'}{$seq_md5}{$sample})){
						push(@{$seq_depths->{$seq_md5}}, $results->{$file}{$marker_name}{'assignments'}{$seq_md5}{$sample});
					}
					$final_seqs{$seq_md5}++;
				}
			}
			$ws_row = 0; $ws_col++;
			#$worksheet->write($ws_row, $ws_col, $file); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=%s',join('+',@depth_amplicon)), undef, eval(join('+',@depth_amplicon))); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $depth_alleles); $ws_row++;
			$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row_first,$ws_col), xl_rowcol_to_cell($ws_row_last,$ws_col)), undef, $count_alleles); $ws_row++;
			if (defined($INP_concatenate)){
				$worksheet->write($ws_row, $ws_col, $file); $ws_row++;
			}
			$worksheet->write($ws_row, $ws_col, $sample, $bold); $ws_row++;
			foreach my $seq_md5 (@seq_md5s) {
# if ($seq_md5 eq '569965416b819382419dc9d3b3469e17'){
# print '';
# }
				# Skip variants not matching alleles if only alleles mode is activated
				if (defined($INP_only_alleles) && !defined($md5_to_name{$seq_md5})){
					next;
				}
				if (defined($seq_depths->{$seq_md5})){
					$worksheet->write($ws_row, $ws_col, sprintf("=%s",join("+",@{$seq_depths->{$seq_md5}})), undef, eval(join("+",@{$seq_depths->{$seq_md5}}))); $ws_row++;
					$final_seq_depths{$seq_md5} += eval(join("+",@{$seq_depths->{$seq_md5}}));
					$final_seq_samples{$seq_md5}++;
				} else {
					$worksheet->write($ws_row, $ws_col, ''); $ws_row++;
				}
			}			
		}
	}
	my $ws_col_last = $ws_col;
	
	printf("\nCombined %d samples (%d shared), %d variants (%d shared).\n", scalar keys %final_samples, scalar keys %common_samples, scalar keys %final_seqs, scalar keys %common_seqs);

	# Writes seq information
	$ws_row = $ws_row_first;
	foreach my $seq_md5 (@seq_md5s) {
		# Skip variants not matching alleles if only alleles mode is activated
		if (defined($INP_only_alleles) && !defined($md5_to_name{$seq_md5})){
			next;
		}
		$ws_col = 0;
# 		foreach my $seq_data_header (@seq_data_headers) {
# 			$worksheet->write_row($ws_row, $ws_col,$results->{$file}{$marker_name}{'seq_data'}{$seq_md5}{lc($seq_data_header)}); $ws_col++;
# 		}
# 		$ws_row++;
		$worksheet->write($ws_row, $ws_col,$md5_to_sequence{$seq_md5}); $ws_col++;
		$worksheet->write($ws_row, $ws_col,$seq_md5); $ws_col++;
		$worksheet->write($ws_row, $ws_col,length($md5_to_sequence{$seq_md5})); $ws_col++;
		$worksheet->write_formula($ws_row, $ws_col, sprintf('=SUM(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $final_seq_depths{$seq_md5}); $ws_col++;
		$worksheet->write_formula($ws_row, $ws_col, sprintf('=COUNT(%s:%s)', xl_rowcol_to_cell($ws_row,$ws_col_first), xl_rowcol_to_cell($ws_row,$ws_col_last)), undef, $final_seq_samples{$seq_md5}); $ws_col++;
		$worksheet->write($ws_row, $ws_col,$md5_to_name{$seq_md5}); $ws_col++;
		$ws_row++;
	}

	
# 	
# 	$stats{'total_assigns'} = $stats{'common_assigns'}+$stats{'file2_missing_assigns'}+$stats{'file1_missing_assigns'};
# 	printf("\nMARKER '%s':\n", $marker_name);
# 	printf("Total unique samples: %d (file1: %d, file2: %d)\n", $stats{'total_samples'}, $stats{'file1_samples'}, $stats{'file2_samples'});
# 	printf("Total seqs: %d (file1: %d, file2: %d)\n", $stats{'total_seqs'}, $stats{'file1_seqs'}, $stats{'file2_seqs'});
# 	printf("Compared samples: %d (excluded from file1: %d, from file2: %d)\n", $stats{'total_combined_samples'}, $stats{'file1_excluded_samples'}, $stats{'file2_excluded_samples'});
# 	printf("Compared seqs: %d (missing in file1: %d, in file2: %d)\n", $stats{'total_combined_seqs'}, $stats{'file1_missing_seqs'}, $stats{'file2_missing_seqs'});
# 	printf("Total assignments: %d (missing in file1: %d, missing in file2: %d)\n", $stats{'total_assigns'}, $stats{'file1_missing_assigns'}, $stats{'file2_missing_assigns'});
# 
}

$workbook->close();

if (defined($INP_verbose)){
	print $verbose_output;
}

print "\nCombined results written into '$INP_output'.\n\n";

exit;


















