#!/usr/bin/perl -w
#
################################################################
#
# Name:  filter_somatic_variants.pl
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
#   Filters variants in VCF files
#
# Examples:
# perl  filter_somatic_variants.pl -i seqs1.fq -u filter_somatic_variants.fa

my $VERSION = "v1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Filters variants in VCF files.";

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

my ($INP_infile,$INP_outfile,$INP_only_snps,$INP_only_indels,$INP_only_pass,$INP_min_total_depth,$INP_min_var_depth,$INP_min_var_freq);

my @argv = @ARGV;

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s' => \$INP_infile,
	'o|outfile=s' => \$INP_outfile,
	'snps' => \$INP_only_snps,
	'indels' => \$INP_only_indels,
	'pass' => \$INP_only_pass,
	'totaldepth=i' => \$INP_min_total_depth,
	'vardepth=i' => \$INP_min_var_depth,
	'varfreq=f' => \$INP_min_var_freq,
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
	print "  -i <file>\tInput VCF file.\n";
	print "  -o <file>\tOutput VCF file.\n";
	print "  -snps\tFilters variants that are not SNPs.\n";
	print "  -indels\tFilters variants that are not indels.\n";
	print "  -pass\tFilters variants that do not 'PASS' all the filters specified in the VCF file.\n";
	print "  -totaldepth <number>\tFilters variants in positions with lower total read depth.\n";
	print "  -vardepth <number>\tFilters variants with lower depth.\n";
	print "  -varfreq <freq>\tFilters variants with lower frequency.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}

# Prints usage help if no input file is specified
if (!defined($INP_infile)){
	print "\nERROR: You must specify a valid input VCF file.\n\n";
	usage();
	exit;
}

print "\nRunning '$COMMAND_LINE'\n";

# Sets output file name if not defined
if (!defined($INP_outfile)){
	# Default filename for result files
	my ($output_folder,$output_name);
	if ($INP_infile =~ /(.+\/)?(.+)\./){
		$output_folder = $1;
		$output_name = $2;
	} else {
		$output_folder = './';
		$output_name = 'variants';
	}
	if (!defined($output_folder)){
		$output_folder = '';
	}
	$INP_outfile = $output_folder.$output_name.'.filtered.vcf';
}

# Extracts single nucleotide variants (SNVs) information
my $vcf_filtered_data = filter_vcf_file($INP_infile);

write_to_file($INP_outfile,$vcf_filtered_data);

print '';



#################################################################################

# Filters data from a Variant Call Format (VCF) file
sub filter_vcf_file {

	my ($file) = @_;

	my @vcf_filtered_data;
	
	# Example VCF4 file:
	# 
	# ##fileformat=VCFv4.0
	# ##fileDate=20090805
	# ##source=myImputationProgramV3.1
	# ##reference=1000GenomesPilot-NCBI36
	# ##phasing=partial
	# ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
	# ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
	# ##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
	# ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
	# ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
	# ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
	# ##FILTER=<ID=q10,Description="Quality below 10">
	# ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
	# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
	# ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
	# ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
	# ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
	# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO					FORMAT		NA00001		NA00002		NA00003
	# 20		14370	.	G	A	29	PASS	NS=3;DP=14;AF=0.5;DB;H2			GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
	# 20		17330	.	T	A	3	q10	NS=3;DP=11;AF=0.017			GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3
	# 20		1110696	.	A	G,T	67	PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4
	# 20		1230237	.	T	.	47	PASS	NS=3;DP=13;AA=T				GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2
	# 20		1234567	.	GTCT	G,GTACT	50	PASS	NS=3;DP=9;AA=G				GT:GQ:DP	0/1:35:4	0/2:17:2	1/1:40:3

	my $standard_genotyping_data = {
			'GT' => "Genotype",
			'GQ' => "Genotype quality",
			'DP' => "Total read depth",
			'AD' => "Allele depths",
			'AF' => "Allele frequency",
			'DB' => "dbSNP membership",
			'H2' => "HapMap2 membership",
			'PL' => "Phred-scaled normalized likelihoods for genotypes"

		
	};

	my $count_total_variants = 0;
	my $count_filtered_variants = 0;
	my $vcf_type;
	my @data_headers;
	open(FILE,$file) || die "# $0 : # 'filter_vcf_file' cannot read '$file'.\n";
	while(my $line = <FILE>){
		$line =~ s/\012\015?|\015\012?//;
		if ($line =~ /^##/){
			push(@vcf_filtered_data,$line);
			if ($line =~ /^##FILTER=<ID=alt_allele_in_normal/i){
				$vcf_type = 'mutect2';
			} elsif ($line =~ /^##FORMAT=<ID=DP4,Number=4/){
				$vcf_type = 'somaticsniper';
			} elsif ($line =~ /^##INFO=<ID=QSS,Number=1/){
				$vcf_type = 'strelka';
			} elsif ($line =~ /^##source=VarScan2/){
				$vcf_type = 'varscan2';
			}
			next;
		} elsif ($line =~ /^#(.+)/ || $line =~ /^(chrom.+)/){ #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR	NORMAL
			push(@vcf_filtered_data,$line);
			@data_headers = split("\t",$1);
# 			if ($line =~ /^chrom\tposition\tref\tvar\tnormal_reads1\tnormal_reads2/){
# 				$vcf_type = 'varscan2';
# 			}
		} elsif ($line =~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)/) {
			$count_total_variants++;
			my @values = split("\t",$line,-1); # -1 forces to split also empty data_headers at the end
			# Reads data_headers
			my %variant_data;
			map $variant_data{$data_headers[$_]} = $values[$_] , 0..$#data_headers;
			# Reference/Chromosome name
			my $ref_name = $variant_data{'CHROM'};
			# Reads variant information (only the TUMOR data, depth and frequency)
			my ($pass,$somatic) = (0,0);
			my ($total_depth,$variant_depth,$variant_frequency);
			my (%tumor_data,%normal_data);
			if (in_array(['mutect2','somaticsniper','strelka','varscan2'],$vcf_type)){
				# FORMAT data header contains the abbreviations/formats of the information in the TUMOR and NORMAL data_headers (GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1)
				# TUMOR and NORMAL data_headers contain genotyping data separated by ':' (ex. 0/1:10,2:0.182:0:2:0.00:353,68:3:7)
				my @headers = split(":",$variant_data{'FORMAT'});
				my @tumor_data = split(":",$variant_data{'TUMOR'});
				# my @normal_data = split(":",$variant_data{'NORMAL'});
				for (my $i=0; $i<=$#headers; $i++){
					$tumor_data{$headers[$i]} = $tumor_data[$i];
					# $normal_data{$headers[$i]} = $normal_data[$i];
				}
			}
			if ($vcf_type eq 'mutect2') {
				if ($variant_data{'FILTER'} eq 'PASS'){
					$pass = 1;
					$somatic = 1;
				} elsif ($variant_data{'FILTER'} !~ /germline_risk/) {
					$somatic = 1;
				}
				my ($ref_depth,$alt_depth) = split(',',$tumor_data{'AD'});
				$total_depth = $ref_depth+$alt_depth;
				$variant_depth = $alt_depth;
			} elsif ($vcf_type eq 'somaticsniper') {
				# All variants are already filtered in Somatic Sniper output
				$pass = 1;
				if ($tumor_data{'SS'} == 2){
					$somatic = 1;
				} # elsif ($tumor_data{'SS'} == 1){
					# $germline = 1;
				# }
				$total_depth = $tumor_data{'DP'};
				my %base_counts;
				($base_counts{'A'},$base_counts{'C'},$base_counts{'G'},$base_counts{'T'}) = split(',',$tumor_data{'BCOUNT'});
				$variant_depth = $base_counts{$variant_data{'ALT'}};
			} elsif ($vcf_type eq 'strelka') {
				if ($variant_data{'FILTER'} eq 'PASS'){
					$pass = 1;
				}
				if ($variant_data{'INFO'} =~ /SOMATIC/) {
					$somatic = 1;
				}
				$total_depth = $tumor_data{'DP'};
				my %base_counts;
				map $base_counts{$_} = (split(',',$tumor_data{$_.'U'}))[0] , ('A','C','G','T');
				$variant_depth = $base_counts{$variant_data{'ALT'}};
			} elsif ($vcf_type eq 'varscan2') {
				if ($variant_data{'FILTER'} eq 'PASS'){
					$pass = 1;
				}
				if ($variant_data{'INFO'} =~ /SOMATIC/) {
					$somatic = 1;
				}
				$total_depth = $tumor_data{'RD'}+$tumor_data{'AD'};
				$variant_depth = $tumor_data{'AD'};
			}

			# Few variants have strange annotations with zero depth
			if (!$total_depth || !$variant_depth){
				next;
			}
			$variant_frequency = sprintf("%.2f",$variant_depth/$total_depth);

			# Filters non somatic variants
			if (!$somatic) {
				next;
			}
			# Filters indels
			my @refs = split(',',$variant_data{'REF'});
			my @alts = split(',',$variant_data{'ALT'});
			if (defined($INP_only_snps)){
				my $indel = 0;
				for (my $i=0; $i<=$#refs; $i++){
					if (length($refs[$i]) != length($alts[$i])){
						$indel = 1;
						last;
					}
				}
				if ($indel){
					next;
				}
			}
			# Filters sSNPs
			if (defined($INP_only_indels)){
				my $snp = 0;
				for (my $i=0; $i<=$#refs; $i++){
					if (length($refs[$i]) == length($alts[$i])){
						$snp = 1;
						last;
					}
				}
				if ($snp){
					next;
				}
			}
			# Filters variants that do not 'PASS' all the filters specified in the VCF file
			if (defined($INP_only_pass) && !$pass){
				next;
			}
			# Filters variants in positions with lower depth.
			if (defined($INP_min_total_depth) && $total_depth < $INP_min_total_depth){
				next;
			}
			# Filters variants with lower depth
			if (defined($INP_min_var_depth) && $variant_depth < $INP_min_var_depth){
				next;
			}
			# Filters variants with lower frequency
			if (defined($INP_min_var_freq) && $variant_frequency < $INP_min_var_freq){
				next;
			}

			push(@vcf_filtered_data,$line);
			$count_filtered_variants++;
		}
	}
	close(FILE);
	
	printf("\n%d tumor SNP variants are kept from the initial %d.\n",$count_filtered_variants,$count_total_variants);
	print "\n";
	
	return join("\n",@vcf_filtered_data);
	
}


#################################################################################



