#!/usr/bin/perl -w
#
################################################################
#
# Name:  annotate_variants.pl
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
#   Annotates variants from a VCF file into protein mutations in reference sequences
#
# Examples:
# perl  annotate_variants.pl -i variants.vcf -r reference_file.fa

my $VERSION = "v1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Annotates variants from a VCF file into protein mutations in reference sequences.";

# Perl modules necessaries for the correct working of the script
use File::Basename;
use lib dirname(__FILE__).'/lib';
use Getopt::Long;
use Bio::Sequences;
use Bio::Onco;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

# my $COMMAND_LINE = $0." ".join(" ",@ARGV);
my $COMMAND_LINE = $SCRIPT_NAME." ".join(" ",@ARGV);

# Default options:
# Output format
my $INP_outformat = 'text short';
# Refences and sequence information
my %paramsdata = (
	'reference_sequences_file' => '/home/alvaro/db/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
	'reference_annotations_file' => '/home/alvaro/db/ensembl/Homo_sapiens.GRCh38.89.gff3.gz',
	# Uniprot ID mapping information file, with Ensembl IDs equivalences
	'idmapping_file' => '/home/alvaro/db/uniprot/HUMAN_9606_idmapping.dat.gz',
	# Folder with UniProt reference sequences and data from oncogenes, tumor suppressors...
	# It incorporates any new protein found in new data
	'oncodata_folder' => './oncodata',
);


my ($INP_infile,$INP_outfile);

my @argv = @ARGV;

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s' => \$INP_infile,
# 	'o|output:s' => \$INP_outfile,
	'f|format=s' => \$INP_outformat,
);

# Usage help
sub usage
{
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -i <file>\tInput VCF file.\n";
# 	print "  -o <file>\tOutput file name\n";
	print "  -f <format>\tOutput format (default='$INP_outformat', 'text long', 'html').\n";
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

# Comment program output if we desire HTML final data
if ($INP_outformat =~ /htm/i){
	$INP_outformat = 'html';
	print "<!--\n";
}

print "\nRunning '$COMMAND_LINE'\n";

# Extracts single nucleotide variants (SNVs) information
my $variants = extract_vcf_data($INP_infile);

# Converts $variants to the format required by 'retrieve_variant_mutations'
# $ref_variants->{$ref_name}{$ref_pos}{sprintf("%s->%s",$ref_nt,$seq_nt)} = $coverage;
my $ref_variants = { 'all' => {}, 'filtered' => {}, 'coding' => {}, 'annotations' => {} };
my %mutation_classes;
foreach my $ref_name (sort {$a cmp $b} keys %$variants){
	if (!defined($variants->{$ref_name}) || !@{$variants->{$ref_name}}){
		next;
	}
	foreach my $variant (@{$variants->{$ref_name}}){
		$ref_variants->{'all'}{$ref_name}{$variant->{'POS'}}{sprintf("%s->%s",$variant->{'REF'},$variant->{'ALT'})} = 1;
		$mutation_classes{$ref_name.':'.$variant->{'POS'}} = 'somatic';
	}
}

# Compares variants with reference sequences to retrieve full gene variant genotype description (mutations)
print "\nAnnotating variant information.\n\n";
$ref_variants = retrieve_variant_mutations($ref_variants,$paramsdata{'reference_sequences_file'},$paramsdata{'reference_annotations_file'});

# Formats and prints variant mutation information
# print "\nFormating and printing variant information.\n\n";
print "\nFormating and printing variant information.\n\n";
my $variants_output = print_variants($ref_variants->{'annotations'},\%mutation_classes,\%paramsdata,[$INP_outformat]);

# Comment program output if we desire HTML final data
if ($INP_outformat eq 'html'){
	print "\n-->\n\n";
}

# Prints variant mutation information
if (defined($variants_output)){
	if ($INP_outformat eq 'html'){
		print "<html>\n<head>\n<link rel=\"stylesheet\" type=\"text/css\" href=\"styles.css\">\n</head>\n<body>\n\n";
		print $variants_output;
		print "\n</body>\n</html>\n\n";
	} else {
		print $variants_output;
	}
} else {
	print "\nERROR: It was an error retrieving variants.\n\n";
	exit;
}



exit;


#################################################################################

# Compares variants with reference sequences to retrieve full gene variant genotype description:
sub retrieve_variant_mutations {

	my ($ref_variants,$ref_sequences_file,$ref_annotations_file,$params) = @_;

	if (!defined($params)){
		$params = {};
	}

	# Opens the Genome FASTA files
	print ("\tProcessing reference file '$ref_sequences_file'.\n");
	my $reference_data;
	if (is_gzip($ref_sequences_file)){
		$reference_data = Bio::SeqIO->new( -file => "zcat '$ref_sequences_file' |", -format => 'fasta') or die "# cannot read $ref_sequences_file\n";
	} elsif (is_zip($ref_sequences_file)){
		$reference_data = Bio::SeqIO->new( -file => "zcat '$ref_sequences_file' |", -format => 'fasta') or die "# cannot read $ref_sequences_file\n";
	} elsif (is_bzip2($ref_sequences_file)){
		$reference_data = Bio::SeqIO->new( -file => "bzcat '$ref_sequences_file' |", -format => 'fasta') or die "# cannot read $ref_sequences_file\n";
	} else {
		$reference_data = Bio::SeqIO->new( -file => "$ref_sequences_file", -format => 'fasta') or die "# cannot read $ref_sequences_file\n";
	}

	# Opens GFF file if the reference file used for mapping is an annotated genome
	my $gff_file_handle;
	if (defined($ref_annotations_file) && -e $ref_annotations_file) {
		# Parses the GFF annotations file
		print ("\tProcessing GFF file '$ref_annotations_file'.\n");
		if (is_gzip($ref_annotations_file)){
			open($gff_file_handle, "zcat '$ref_annotations_file' |") or die "# cannot read $ref_annotations_file\n";
		} elsif (is_zip($ref_annotations_file)){
			open($gff_file_handle, "zcat '$ref_annotations_file' |") or die "# cannot read $ref_annotations_file\n";
		} elsif (is_bzip2($ref_annotations_file)){
			open($gff_file_handle, "bzcat '$ref_annotations_file' |") or die "# cannot read $ref_annotations_file\n";
		} else {
			open($gff_file_handle,$ref_annotations_file) || die "# cannot read $ref_annotations_file\n";
		}
	}

	# Compares variants with reference sequences
	while (my $ref_obj = $reference_data->next_seq){ # Goes to next reference sequence (eg. next chromosome) if the alignments of the previous one have finished
		
		my $ref_id = $ref_obj->primary_id();

		my ($ref_variants_coding,$ref_variants_annotations,$ref_coding_mapped_length) = annotate_reference_coding_variants($ref_variants->{'all'}{$ref_id},$ref_obj,$gff_file_handle);
		$ref_variants->{'coding'}{$ref_id} = $ref_variants_coding;
		# $ref_mapped_lengths->{'coding'} += $ref_coding_mapped_length;
		if (%$ref_variants_annotations) {
			$ref_variants->{'annotations'} = { %{$ref_variants->{'annotations'}},  %$ref_variants_annotations };
		}
	}



	return $ref_variants;




}

#################################################################################


