#!/usr/bin/perl -w
#
################################################################
#
# Name:  common_variants.pl
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
#   Calculates the number of common and non common variants between several VCF files
#   Creates a VCF file with the common and non common (optional) variants 
#
# Examples:
# perl common_variants.pl -min 2 variants1.vcf variants2.vcf variants3.vcf -o common_variants.vcf -non noncommon_variants.vcf

my $VERSION = "v1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Calculates the number of common variants between several VCF files.";

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

my (@INP_infiles,$INP_outfile,$INP_min_common,$INP_noncommon);

my @argv = @ARGV;

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s{,}' => \@INP_infiles,
	'o|output:s' => \$INP_outfile,
	'min|mincommon=i' => \$INP_min_common,
	'non|noncommon:s' => \$INP_noncommon,
);

# Usage help
sub usage
{
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -i <file1> <file2> [ ... <fileN> ]\n\t\tVCF files to combine\n";
	print "  -o <file>\tOutput VCF file name\n";
	print "  -min <number>\tMinimum number of files containing the variant (default=all)\n";
	print "  -non <file>\tPrints output file with non common variants\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}

# Prints usage help if no input file is specified
if (!@INP_infiles || !-e $INP_infiles[0] || !-e $INP_infiles[1]){
	print "\nERROR: You must specify at least 2 valid input VCF files to compare.\n\n";
	usage();
	exit;
}
if (!defined($INP_outfile)){
	$INP_outfile = 'common_variants.vcf';
}
if (defined($INP_noncommon) && is_numeric($INP_noncommon) && $INP_noncommon == 1){
	$INP_noncommon = 'noncommon_variants.vcf';
}
if (!defined($INP_min_common)){
	$INP_min_common = scalar @INP_infiles;
}

print "\nRunning '$COMMAND_LINE'\n";

# Put variants from every input file into a hash
my ($vcf_variants, $vcf_headers);
foreach my $file (@INP_infiles) {
	open(FILE,$file) || die "# $0 : # 'filter_vcf_file' cannot read '$file'.\n";
	my @data_headers;
	while(my $line = <FILE>){
		$line =~ s/\012\015?|\015\012?//;
		if ($line =~ /^##/){
			$vcf_headers->{$file} .= $line;
			next;
		} elsif ($line =~ /^#(.+)/ || $line =~ /^(chrom.+)/){ #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR	NORMAL
			@data_headers = split("\t",$1);
			$vcf_headers->{$file} .= $line;
		} elsif (split("\t",$line) > 4) {
			my @values = split("\t",$line,-1); # -1 forces to split also empty data_headers at the end
			# Reads data_headers
			my %variant_data;
			map $variant_data{$data_headers[$_]} = $values[$_] , 0..$#data_headers;
			# Stores the variants and their positions into the $all_variants HASH
			my $variant;
			if (defined($variant_data{'CHROM'})){
				$variant = sprintf("%s\t%s\t%s\t%s\t%s",$variant_data{'CHROM'},$variant_data{'POS'},$variant_data{'ID'},$variant_data{'REF'},$variant_data{'ALT'});
			} elsif  (defined($variant_data{'chrom'})){ # VarScan2 format
				$variant = sprintf("%s\t%s\t%s\t%s\t%s",$variant_data{'chrom'},$variant_data{'position'},'.',$variant_data{'ref'},$variant_data{'var'});
			}
			if ($variant) {
				$vcf_variants->{$file}{$variant} = $line;
			} else {
				print "\nERROR: It was an error reading the variant information at line:\n$line\n\n";
			}
		}
	}
	close(FILE);
}

# Annotates total variants
my $total_variants;
foreach my $file (@INP_infiles) {
	foreach my $variant (keys %{$vcf_variants->{$file}}){
		my @values = split("\t",$variant);
		my $chrom = shift @values;
		my $pos = shift @values;
		my $variant_ = join("\t",@values);
		if (!defined($total_variants->{$chrom}{$pos}{$variant_}{$file})){
			# Groups variants by chromosome to sort them when printing
			$total_variants->{$chrom}{$pos}{$variant_}{$file} = 1;
		}
	}
}

# Annotates common variants
my ($count_common_variants, $count_noncommon_variants) = (0,0);
my %common_variants_per_file;
map $common_variants_per_file{$_} = 0, @INP_infiles;
my ($vcf_common_data, $vcf_noncommon_data);
$vcf_common_data .= "#CHROM	POS	ID	REF	ALT\n";
foreach my $chrom (sort {$a cmp $b} keys %$total_variants){
	foreach my $pos (sort {$a<=>$b} keys %{$total_variants->{$chrom}}){
		foreach my $variant_ (keys %{$total_variants->{$chrom}{$pos}}){
			my $variant = "$chrom\t$pos\t$variant_";
			if (scalar keys %{$total_variants->{$chrom}{$pos}{$variant_}} >= $INP_min_common){
				$vcf_common_data .= "$variant\n";
				$count_common_variants++;
				# Annotates the number of common variants per file
				foreach my $file (keys %{$total_variants->{$chrom}{$pos}{$variant_}}){
					$common_variants_per_file{$file}++;
				}
			} else {
				$vcf_noncommon_data .= "$variant\n";
				$count_noncommon_variants++;
			}
		}
	}
}

printf("\n%d unique total variants.\n", $count_common_variants+$count_noncommon_variants);
# Prints the number of total variants per file
foreach my $file (@INP_infiles) {
	printf("\t%d total variants in '%s'.\n", scalar keys %{$vcf_variants->{$file}}, $file);
}


write_to_file($INP_outfile,$vcf_common_data);
if (scalar @INP_infiles == $INP_min_common){
	printf("\n%d common variants present in %d files written into '%s'.\n", $count_common_variants, $INP_min_common, $INP_outfile);
} else {
	printf("\n%d common variants present in %d files or more written into '%s'.\n", $count_common_variants, $INP_min_common, $INP_outfile);
}
# Prints the number of common variants per file
foreach my $file (@INP_infiles) {
	printf("\t%d common variants in '%s'.\n", $common_variants_per_file{$file}, $file);
}


if (defined($INP_noncommon)){
	write_to_file($INP_noncommon,$vcf_noncommon_data);
	printf("\n%d noncommon variants written into '%s'.\n", $count_noncommon_variants, $INP_noncommon);
}

print "\n";

exit;







