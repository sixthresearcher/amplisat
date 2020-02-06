#!/usr/bin/perl -w
#
################################################################
#
# Name:  compare_variants.pl
#
# Version: 1.0
#
# Author: Alvaro Sebastian
#
# Support: Alvaro Sebastian (sixthresearcher@gmail.com)
#
# Sixth Researcher
# http://www.sixthresearcher.com
#
# Description:
#   Compares variants (each variant one line) between 2 files
#
# Examples:
# perl  compare_variants.pl -i seqs1.fq -u compare_variants.fa

my $VERSION = "v1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Compares variants (each variant one line) between 2 files";

# Perl modules necessaries for the correct working of the script
use File::Basename;
use lib dirname(__FILE__).'/lib';
use Getopt::Long;
use Bio::Sequences;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

# my $COMMAND_LINE = $0." ".join(" ",@ARGV);
my $COMMAND_LINE = $SCRIPT_NAME." ".join(" ",@ARGV);

my ($INP_file1,$INP_file2,$INP_noncommon,$INP_output);

my @argv = @ARGV;

GetOptions(
	'h|help|?' =>  \&usage,
	'1=s' => \$INP_file1,
	'2=s' => \$INP_file2,
	'n|noncommon' => \$INP_noncommon,
	'o|output=s' => \$INP_output,
	'<>' => \&usage,
);

# Usage help
sub usage
{
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <infile> -o <outfile> [options]\n";
	print "\nOptions:\n";
	print "  -1 <file>\tFirst variants file\n";
	print "  -2 <file>\tSecond variants file\n";
	print "  -n <file>\tPrint non-common variants\n";
	print "  -o <file>\tOutput file name\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}

# Prints usage help if no input file is specified
if (!defined($INP_file1) || !defined($INP_file2)){
	print "\nERROR: You must specify 2 valid text files with variants to compare.\n\n";
	usage();
	exit;
}

print "\nRunning '$COMMAND_LINE'\n";

printf("Reading File '%s'.\n", $INP_file1);
my ($variants1, $count_variants1) = extract_variants_from_file($INP_file1);
printf("\t%d variants found.\n", $count_variants1);

printf("Reading File '%s'.\n", $INP_file2);
my ($variants2, $count_variants2) = extract_variants_from_file($INP_file2);
printf("\t%d variants found.\n", $count_variants2);

printf("\nComparing variants.\n");

my ($common_variants, $noncommon_variants);
my ($count_common_variants,$count_noncommon_variants) = (0,0);
foreach my $gene1 (keys %$variants1){
if ($gene1 =~ /ENSG00000197172|ENSG00000108510|ENSG00000121892/){
print '';
}
	foreach my $variant1 (keys %{$variants1->{$gene1}}){
		if (defined($variants2->{$gene1}) && defined($variants2->{$gene1}{$variant1})){
			$common_variants->{$gene1}{$variant1}++;
			$count_common_variants++;
		} else {
			$noncommon_variants->{$gene1}{$variant1}++;
			$count_noncommon_variants++;
		}
	}
}

printf("\nCOMMON VARIANTS:\n\n");
foreach my $gene (keys %$common_variants){
	foreach my $variant (keys %{$common_variants->{$gene}}){
		printf("%s\t%s\n",$gene,$variant);
	}
}
if (defined($INP_noncommon)){
	printf("\n\nNON-COMMON VARIANTS:\n\n");
	foreach my $gene (keys %$noncommon_variants){
		foreach my $variant (keys %{$noncommon_variants->{$gene}}){
			printf("%s\t%s\n",$gene,$variant);
		}
	}
}

printf("\n%d common variants and %d non-common found.\n\n", $count_common_variants, $count_noncommon_variants);





#################################################################################

# Reads content from a file
sub extract_variants_from_file {
	my ($file) = @_;
	
	my $variants;
	my $count_variants = 0;

	my $lines = read_from_file($file);

	foreach my $line (@$lines) {
		my ($gene,$variant) = ('','');
		if ($line =~ /^#/){ # Skips comments
			next;
		}
# 		# OBSOLETE: Genome positions change in different genome versions
# 		if ($line =~ /(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)\:(\d{3,})/) {
# 			$gene = "$1:$2";
# 		}
		if ($line =~ /(ENSG\d{11})/){
			$gene = $1;
		}
		if ($line =~ /c\.(A|C|G|T)(\d+)(A|C|G|T|\*)/){
			$variant = "$1$2$3";
		} elsif ($line =~ /(\d+) (A|C|G|T)->(A|C|G|T|\*)/) {
			$variant = "$2$1$3";
		}
# if ($gene eq 'ENSG00000084774'){
# print '';
# }
		if ($gene && $variant){
			if (!defined($variants->{$gene}{$variant})){
				$count_variants++;
			}
			$variants->{$gene}{$variant}++;
		}
		print '';
	}

	return ($variants,$count_variants);
}

#################################################################################
