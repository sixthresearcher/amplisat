#!/usr/bin/perl -w
#
################################################################
#
# Name: ampliON.pl
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
#   Fast analysis of NGS data (RNA-seq) to identify common mutations and alterations in therapeutic target cancer genes
#
# Requires as input one FASTQ file (compressed or uncompressed) with RNA-seq data from tumour cells
#
# Example:
# perl ampliON.pl -i reads.fq.gz -t lung -o results
#

my $VERSION = "1.0";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Fast analysis of NGS data (RNA-seq) to identify common mutations and alterations in therapeutic target cancer genes.";

# Modules are in folder 'lib' in the path of the script
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;
use Bio::Onco;
use Data::Dumper;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $SCRIPT_DIR = dirname(__FILE__);

# my $COMMAND_LINE = $0." ".join(" ",@ARGV);
my $COMMAND_LINE = $SCRIPT_NAME." ".join(" ",@ARGV);

# Default options
# Output format
my $INP_format = 'text short';

# Folder with UniProt reference sequences and data from oncogenes, tumor suppressors...
# It incorporates any new protein found in new data
my $INP_oncodata = $SCRIPT_DIR.'/oncodata';

# Ensembl human reference CDs/transcripts sequences file
my $INP_transcript_seqs_file = "/home/alvaro/ngs/human/Homo_sapiens.GRCh38.cds.all.fa.gz";
# Uniprot human reference protein sequences file
my $INP_protein_seqs_file = "/home/alvaro/evobio/scripts/uniprot/uniprot_human_refs.fa.gz";
# Uniprot ID mapping information file, with Ensembl IDs equivalences
my $INP_uniprot_idmapping_file = "/home/alvaro/ngs/human/HUMAN_9606_idmapping.dat.gz";

# Mapping options
my $INP_align = 'bowtie2 --very-fast';
# my $INP_align = 'bowtie2 --sensitive';
# my $INP_align = 'bowtie2 --very-fast-local';
# my $INP_align = 'bowtie2 --sensitive-local';

my (@INP_reads_files, $INP_outpath, $INP_threads, $INP_nreads, $INP_shuffle, $INP_synonymous, $INP_all_mutations, $INP_test, $INP_zip);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s{,}' => \@INP_reads_files,
	'o|output=s' => \$INP_outpath,
	'f|format=s' => \$INP_format,
	'a|align=s' => \$INP_align,
# 	's|shuffle' => \$INP_shuffle,
	'n|number=i' => \$INP_nreads,
	'syn' => \$INP_synonymous,
	'all' => \$INP_all_mutations,
	'thr|threads=i' => \$INP_threads,
	'test' => \$INP_test,
	'z|zip' => \$INP_zip,
	'<>' => \&usage,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -i <file1> [<file2>]\n\t\tInput single-end file or paired-end read files in FASTQ format (compressed or uncompressed).\n";
	print "  -o <path>\tOutput folder name.\n";
	print "  -f <format>\tOutput format (default='$INP_format', 'text short', 'html').\n";
	print "  -a <string>\tAlignment/mapping type and parameters (default='$INP_align').\n";
# 	print "  -s\t\tShuffle/randomize reads/sequences to analyze.\n";
	print "  -n <number>\tNumber of reads/sequences to analyze.\n";
	print "  -syn\t\tAnalyze also synonymous mutations.\n";
	print "  -all\t\tAnalyze all kind of mutations.\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	#  print " -auto\t\tParameters are automatically adjusted by a preliminar analysis of amplicons (DEFAULT).\n";
	#  print " -test\t\t Run in 'test' mode and create a dump file of the alignments.\n";
	exit;
}

# Checks if there are input files
if (!@INP_reads_files || !-e $INP_reads_files[0]){
	print "\nERROR: You must specify one single-end or two paired-end read input files.\n\n";
	usage();
	exit;
}
# Default output path
if (!defined($INP_outpath)){
	$INP_outpath  = lc((split('\.',$SCRIPT_NAME))[0]);
}
# Creates output folder
if (!-d $INP_outpath){
	mkdir($INP_outpath);
}
print "\nRunning '$COMMAND_LINE'\n";

# Defines the filename for output files
my ($filename,$filepath,$filesuffix) = fileparse($INP_reads_files[0],qr/\.[^.]*/); 

# Reads the alignment type and parameters
my ($align_type,$align_params) = ('','');
if ($INP_align =~ /(.+?)\s+(.+)/i){
	$align_type = $1;
	$align_params = $2;
}
$align_type =~ s/\s+//g;


# Performs the alignment/mapping of reads to references
my ($infile,$outfile);
if (is_fastq($INP_reads_files[0]) || is_fasta($INP_reads_files[0])){
	# Ex. $INP_align = 'bowtie2 --very-fast-local';
	if ($align_type eq 'bowtie2'){
		# Skips unaligned reads from results
		$align_params .= ' --no-unal';
		if (defined($INP_threads)){
			$align_params .= " -p $INP_threads";
		}
		$outfile = "$INP_outpath/$filename.sam";
		print "\nMapping ";
		execute_bowtie2(\@INP_reads_files,$INP_transcript_seqs_file,$align_params,$outfile);
		# Print mapping results into a file
		print "\nMapping results stored into '$outfile'.\n";
	}
} else {
	$outfile = $INP_reads_files[0];
}

# Converts SAM to BAM
$infile = $outfile;
if (is_sam($infile)){
	$outfile = "$INP_outpath/$filename.bam";
	print "\nCreating '$outfile' file.\n";
	system("samtools view -b -S -o $outfile $infile");
	# `rm $infile`;
} elsif ($infile ne $INP_reads_files[0]) {
	print "\nERROR: It was an error in the read mapping and SAM file could not be generated.\n\n";
	exit;
}

# Sorts and indexes the BAM file
$infile = $outfile;
if (is_bam($infile)){
	$outfile = "$INP_outpath/$filename";
	print "\nSorting and indexing '$infile' file.\n";
	system("samtools sort $infile $outfile");
	system("samtools index $infile");
	# Creates a BCF file with the probabilities of the different genotypes for each position in the reference sequences
	$outfile = "$INP_outpath/$filename.bcf";
	print "\nCreating '$outfile' file.\n";
	system("samtools mpileup -C50 -D -g -f $INP_transcript_seqs_file $infile > $outfile");
# 	system("samtools mpileup -g -f $INP_transcript_seqs_file $infile > $outfile");
} elsif ($infile ne $INP_reads_files[0]) {
	print "\nERROR: It was an error in the read mapping and BAM file could not be generated.\n\n";
	exit;
}

# Computes genotype likelihoods and creates the VCF file
$infile = $outfile;
if (is_bcf($infile)){
	$outfile = "$INP_outpath/$filename.vcf";
	print "\nCreating '$outfile' file.\n\n";
	system("bcftools view -vcg $infile > $outfile");
# 	system("bcftools view -c -v $infile > $outfile");
} elsif ($infile ne $INP_reads_files[0]) {
	print "\nERROR: It was an error in the identification of genotype likelihoods and BCF file could not be generated.\n\n";
	exit;
}
if (!is_vcf($outfile)){
	if ($outfile eq $INP_reads_files[0]) {
		print "\nERROR: '".$INP_reads_files[0]."' format is not supported.\n\n";
		usage();
		exit;
	} else {
		print "\nERROR: It was an error in the identification of genotype variants and VCF file could not be generated.\n\n";
		exit;
	}
}
# Extracts single nucleotide variants (SNVs) information
my $variants = extract_vcf_data($outfile);

# # Reads coverage data from .BAM file 
# my $coverages;
# if (-e "$INP_outpath/$filename.bam"){
# 	$coverages = extract_coverage_data("$INP_outpath/$filename.bam");
# } else {
# 	print "\nERROR: Coverages could not be read, '$INP_outpath/$filename.bam' file does not exist.\n\n";
# }

# Compares variants with reference sequences to retrieve full gene variant genotype description (mutations)
my $mutations = retrieve_variant_mutations(read_fasta_file_hash($INP_transcript_seqs_file,'ensembl'),$variants);

# Formats variant mutation information
my @print_options = ($INP_format);
if (defined($INP_all_mutations)){
	push(@print_options,'all');
} elsif (defined($INP_synonymous)){
	push(@print_options,'synonymous');
}

# Reads Uniprot ID mapping information file (cross-references to external databases)
my $ensembl_to_uniprot = read_uniprot_idmapping_file($INP_uniprot_idmapping_file,'Ensembl_TRS',-1);
my $uniprot_to_genename = read_uniprot_idmapping_file($INP_uniprot_idmapping_file,'Gene_Name');

# Convert Ensembl IDs into UniProt ones and retrieves protein information from UniProt if it has not been retrieved before
my $uniprot_data;
foreach my $ensembl_full_id (keys %$mutations){
	my $ensembl_id = $ensembl_full_id;
	if ($ensembl_id =~ /(ENST\d+)\.\d+/){
		$ensembl_id = $1;
	}
	# Converts the Ensembl ID into Uniprot one
	if (!defined($ensembl_to_uniprot->{$ensembl_id})){
		print "\nWARNING: '$ensembl_id' has no UniProt ID associated.\n";
		next;
	} elsif ($#{$ensembl_to_uniprot->{$ensembl_id}} > 0){
		print "\nWARNING: '$ensembl_id' has multiple UniProt ID associated.\n";
	}
	my $uniprot_id = $ensembl_to_uniprot->{$ensembl_id}[0];
	if ($uniprot_id =~ /(\w+?)\-\d+/){
		$uniprot_id = $1;
	}
	$mutations->{$uniprot_id} = $mutations->{$ensembl_full_id};
	delete($mutations->{$ensembl_full_id});
	# Retrieves protein data from UniProt (if not already present in the 'oncodata' folder)
	if (defined($uniprot_to_genename->{$uniprot_id})){
		my $gene_name = $uniprot_to_genename->{$uniprot_id}[0];
		my $uniprot_data_file = sprintf("%s/%s.txt", $INP_oncodata, lc($gene_name));
		if (!-e $uniprot_data_file){
			write_to_file($uniprot_data_file,retrieve_uniprot_data($uniprot_id,'txt'));
		}
		$uniprot_data->{$uniprot_id} = read_uniprot_single_file($uniprot_data_file)
	}
}

my $mutations_output = print_variant_mutations($mutations,$uniprot_data,\@print_options);
# Prints variant mutation information
if (defined($mutations_output)){
	if ($INP_format =~ /text/){
		print "\nVARIANT REPORT:\n";
		print "$mutations_output\n";
	}
} else {
	print "\nNo detected mutations in the analysis.\n\n";
	exit;
}


#################################################################################

# Compares variants with reference sequences to retrieve full gene variant genotype description:
sub retrieve_variant_mutations {

	my ($references, $variants, $params) = @_;
	
	my $mutations;
	
	if (!defined($params)){
		$params = {};
	}

	# Compares variants with reference sequences
# 	for (my $i=0; $i<=$#{$ref_names}; $i++) {
	foreach my $ref_name (keys %$variants){
		# Reference sequence data
		my $ref_seq = $references->{$ref_name};
		my $ref_prot = dna_to_prot($ref_seq);
# 		my $ref_name = $ref_names->[$i];
# 		my $ref_seq = $ref_seqs->[$i];
# 		my $ref_prot = dna_to_prot($ref_seq);
		# Counts the number of stop codons
		# my $ref_stops = () = $alt_prot =~ /\*/gi;

		# Read variants and stores variant sequence data
		my (@synonymous_info,@nonsynonymous_info);
		foreach my $variant (@{$variants->{$ref_name}}){
			my @genotypes = split(",",$variant->{'GENOTYPE'});
			# Reference nucleotides
			my $ref_variant_nuc = substr($ref_seq,$variant->{'POS'}-1,length($variant->{'REF'}));
			foreach my $genotype (@genotypes){
				my $type;
				if (scalar @genotypes == 1) {
					$type = 'Homozygous ';
				} else {
					$type = 'Heterozygous ';
				}
				# Skips low coverage genotypes
				if (defined($params->{'min_depth'}) && $variant->{'DP'} < $params->{'min_depth'}){
					next;
				}
				# Skips low quality genotypes
				if (defined($params->{'min_qual'}) && $variant->{'GQ'} < $params->{'min_qual'}){
					next;
				}
				# Skips extreme positions with alignment problems
				if (defined($params->{'min_pos'}) && ($variant->{'POS'} < $params->{'min_pos'} || $variant->{'POS'} > length($ref_seq)-$params->{'min_pos'})){
					next;
				}
				my $alt_seq = substr($ref_seq,0,$variant->{'POS'}-1).$genotype.substr($ref_seq,$variant->{'POS'}+length($variant->{'REF'})-1);
				my $alt_prot = dna_to_prot($alt_seq);
				# Variant nucleotides
				my $alt_variant_nuc = substr($alt_seq,$variant->{'POS'}-1,length($genotype));
				my $alt_variant_nuc_first = $variant->{'POS'};
				my $alt_variant_nuc_last = $variant->{'POS'}+length($genotype)-1;
				if ($ref_variant_nuc eq $alt_variant_nuc){
					next;
				}
				# First and last position of the variant amino acids
				my $alt_variant_prot_first = int(($alt_variant_nuc_first+2)/3);
				my $alt_variant_prot_last = int(($alt_variant_nuc_last+2)/3);
				my $ref_variant_prot = substr($ref_prot,$alt_variant_prot_first-1,$alt_variant_prot_last-$alt_variant_prot_first+1);
				my $alt_variant_prot = substr($alt_prot,$alt_variant_prot_first-1,$alt_variant_prot_last-$alt_variant_prot_first+1);
				# $alt_stops = () = $alt_prot =~ /\*/gi;
				# missense mutation: a single nucleotide change results in a codon that codes for a different amino acid.
				# nonsense mutation: a codon is changed to a premature stop codon that results in truncation of the resulting protein.
				# silent mutation: mutations in DNA that do not significantly alter the phenotype of the organism in which they occur:
				#                  1. Do not result in a change to the amino acid sequence of a protein (synonymous substitution)
				#                  2. Result in the insertion of an alternative amino acid with similar properties to that of the original amino acid
				my $mutation = {};
				if ($ref_variant_prot ne $alt_variant_prot){ # Non-synonymous substitution/indel
					if (length($ref_variant_nuc) == 1 && length($alt_variant_nuc) == 1){
						if ($alt_variant_prot =~ /\*/) {
							$type .= 'nonsense'; # Mutation introduces an stop codon
						} else {
							$type .= 'missense'; # Mutation changes protein sequence
						}
						# if ($ref_name eq 'EGFR' && $variant->{'POS'}==1620){
						# 	print '';
						# }
						$mutation->{'type'} = "$type substitution";
						$mutation->{'description'} = sprintf("%d %s->%s (p.%s%d%s)",$variant->{'POS'},$variant->{'REF'},$genotype,$ref_variant_prot,$alt_variant_prot_first,$alt_variant_prot);
						# $mutation->{'description'} = sprintf("%s%d%s (p.%s%d%s, cd.%s%d%s)",$ref_variant_prot,$alt_variant_prot_first,$alt_variant_prot,convert_aa_1_to_3($ref_variant_prot),$alt_variant_prot_first,convert_aa_1_to_3($alt_variant_prot),$ref_variant_nuc,$variant->{'POS'},$alt_variant_nuc);
					} elsif (length($ref_variant_nuc) < length($alt_variant_nuc)){
						my $insertion_nuc_length = length($alt_variant_nuc) - length($ref_variant_nuc);
						my $insertion_nuc_first = $variant->{'POS'}+length($ref_variant_nuc)-1;
						my $insertion_nuc_last = $insertion_nuc_first+1;
						my $insertion_nuc_seq = substr($alt_variant_nuc,length($ref_variant_nuc),$insertion_nuc_length);
						if ($alt_variant_prot =~ /\*/) {
							$type .= 'nonsense'; # Mutation introduces an stop codon
						} elsif ($insertion_nuc_first % 3 ==0 && $insertion_nuc_length % 3 == 0)  {
							$type .= 'inframe'; # Mutation changes the reading frame
						} else {
							$type .= 'frameshift'; # Mutation changes the reading frame
						}
						$mutation->{'type'} = "$type insertion";
						if ($type !~ /inframe/){
							$mutation->{'description'} = sprintf("%d %s->%s",$variant->{'POS'},$variant->{'REF'},$genotype);
							# $mutation->{'description'} = sprintf("cd.%d\_%dins%s",$insertion_nuc_first,$insertion_nuc_last,$insertion_nuc_seq);
						} else {
							my $ref_insertion_aa_first = int(($insertion_nuc_first+2)/3);
							my $ref_insertion_aa_last = $ref_insertion_aa_first+1;
							my $ref_insertion_aa_first_seq = substr($ref_prot,$ref_insertion_aa_first-1,1);
							my $ref_insertion_aa_last_seq = substr($ref_prot,$ref_insertion_aa_last-1,1);
							my $alt_insertion_aa_seq = substr($alt_prot,$ref_insertion_aa_first,$insertion_nuc_length/3);
							$mutation->{'description'} = sprintf("%d %s->%s (p.%s%d\_%s%dins%s)",$variant->{'POS'},$variant->{'REF'},$genotype,$ref_insertion_aa_first_seq,$ref_insertion_aa_first,$ref_insertion_aa_last_seq,$ref_insertion_aa_last,$alt_insertion_aa_seq);
							# $mutation->{'description'} = sprintf("p.%s%d\_%s%dins%s (cd.%d\_%dins%s)",$ref_insertion_aa_first_seq,$ref_insertion_aa_first,$ref_insertion_aa_last_seq,$ref_insertion_aa_last,$alt_insertion_aa_seq,$insertion_nuc_first,$insertion_nuc_last,$insertion_nuc_seq);
						}
					} else {
						my $deletion_nuc_length = length($ref_variant_nuc) - length($alt_variant_nuc);
						my $deletion_nuc_first = $variant->{'POS'}+length($alt_variant_nuc);
						my $deletion_nuc_last = $deletion_nuc_first+$deletion_nuc_length-1;
						my $deletion_nuc_seq = substr($ref_variant_nuc,length($alt_variant_nuc),$deletion_nuc_length);
						if ($alt_variant_prot =~ /\*/) {
							$type .= 'nonsense'; # Mutation introduces an stop codon
						} elsif (($deletion_nuc_first-1) % 3 ==0 && $deletion_nuc_length % 3 == 0)  {
							$type .= 'inframe'; # Mutation changes the reading frame
						} else {
							$type .= 'frameshift'; # Mutation changes the reading frame
						}
						$mutation->{'type'} = "$type deletion";
						if ($type !~ /inframe/){
							$mutation->{'description'} = sprintf("%d %s->%s",$variant->{'POS'},$variant->{'REF'},$genotype);
							# $mutation->{'description'} = sprintf("cd.%d\_%ddel",$deletion_nuc_first,$deletion_nuc_last);
						} else {
							my $ref_deletion_aa_first = int(($deletion_nuc_first+2)/3);
							my $ref_deletion_aa_last = int(($deletion_nuc_last+2)/3);
							my $ref_deletion_aa_first_seq = substr($ref_prot,$ref_deletion_aa_first-1,1);
							my $ref_deletion_aa_last_seq = substr($ref_prot,$ref_deletion_aa_last-1,1);
							$mutation->{'description'} = sprintf("%d %s->%s (p.%s%d\_%s%ddel)",$variant->{'POS'},$variant->{'REF'},$genotype,$ref_deletion_aa_first_seq,$ref_deletion_aa_first,$ref_deletion_aa_last_seq,$ref_deletion_aa_last);
							# $mutation->{'description'} = sprintf("p.%s%d\_%s%ddel (cd.%d\_%ddel)",$ref_deletion_aa_first_seq,$ref_deletion_aa_first,$ref_deletion_aa_last_seq,$ref_deletion_aa_last,$deletion_nuc_first,$deletion_nuc_last);
						}
	# 					push(@nonsynonymous_info, sprintf("\t%s deletion: %s%d%s (cd.%s%d%s)",$type,$ref_variant_prot,$alt_variant_prot_first,$alt_variant_prot,$ref_variant_nuc,$variant->{'POS'},$alt_variant_nuc));
					}
				} else { # Synonymous substitution/indel (protein sequence doesn't change)
					$type .= 'synonymous';
					if (length($ref_variant_prot) == 1 && length($alt_variant_prot) == 1){
						$mutation->{'type'} = "$type substitution";
						$mutation->{'description'} = sprintf("%d %s->%s (p.%s%d=)",$variant->{'POS'},$variant->{'REF'},$genotype,$ref_variant_prot,$alt_variant_prot_first,convert_aa_1_to_3($ref_variant_prot),$alt_variant_prot_first,$ref_variant_nuc,$variant->{'POS'},$alt_variant_nuc);
						# $mutation->{'description'} = sprintf("%s%d= (p.%s%d=, cd.%s%d%s)",$ref_variant_prot,$alt_variant_prot_first,convert_aa_1_to_3($ref_variant_prot),$alt_variant_prot_first,$ref_variant_nuc,$variant->{'POS'},$alt_variant_nuc);
					} else {
						$mutation->{'type'} = "$type mutation";
						$mutation->{'description'} = sprintf("%d %s->%s (p.%s%d=)",$variant->{'POS'},$variant->{'REF'},$genotype,$ref_variant_prot,$alt_variant_prot_first);
					}
				}
				push(@{$mutations->{$ref_name}},$mutation);
	# 			if ($ref_name eq 'KRAS') {
	# 				print "\n>REF\n$ref_seq\n>ALT\n$alt_seq\n";
	# 				print "\n>REF\n$ref_prot\n>ALT\n$alt_prot\n";
	# 				print "\n>REF\n$ref_variant_prot\n>ALT\n$alt_variant_prot\n";
	# 				print '';
	# 			}
			}
		}
	}

	return $mutations;

}

#################################################################################

# Prints variant mutation information
sub print_variant_mutations {

	my ($mutations,$ref_data,$options) = @_;

	my $output;

	# Sets options
	my ($format,$detailed,$all,$synonymous)=('text',0,0,0,0);
	foreach my $option (@$options){
		if ($option =~ /(text|html) (short|long)/){
			$format = $1;
			if ($2 eq 'long'){
				$detailed = 1;
			}
		}
		if ($option eq 'all'){
			$all = 1;
		}
		if ($option eq 'synonymous'){
			$synonymous = 1;
		}
	}

	# Loops all the variants
	foreach my $ref_name (sort {$a cmp $b} keys %$mutations){

		my $gene_output;
		
		my $gene_name = sprintf("%s - %s",$ref_data->{$ref_name}{'gene_names'}[0],$ref_data->{$ref_name}{'names'}[0]);
		my $gene_description = $ref_data->{$ref_name}{'gene_names'}[0];
		my $gene_function = $ref_data->{$ref_name}{'function'}[0];
		my (@gene_diseases, @gene_variants);
		if (defined($ref_data->{$ref_name}{'disease'})){
			@gene_diseases = @{$ref_data->{$ref_name}{'disease'}};
		}
		if (defined($ref_data->{$ref_name}{'variants'})){
			@gene_variants = @{$ref_data->{$ref_name}{'variants'}};
		}

		# Prints full gene variant genotype descriptions:
		foreach my $mutation (@{$mutations->{$ref_name}}){
			if (!$all && !$synonymous && $mutation->{'type'} =~ / synonymous/) {
				next;
			}
			if ($format eq 'text'){
				if ($mutation->{'type'} !~ / synonymous/){
					$gene_output .= sprintf("\t%-35s %s\n",$mutation->{'type'}.':',$mutation->{'description'});
				} else {
					$gene_output .= sprintf("\t# %-37s %s\n",$mutation->{'type'}.':',$mutation->{'description'});
				}
			}
		}
		if (($all || $synonymous) && !defined($gene_output)){
			if ($format eq 'text'){
				$gene_output .= sprintf("\t# No mutations detected\n");
			}
		}
		if (defined($gene_output)){
			if ($format eq 'text'){
				# $gene_output = sprintf("\nGene %s:\n%s", $ref_name, $gene_output);
# 				$gene_output = sprintf("\nGene %s:\n\n%s", $gene_name, $gene_output);
				$gene_output = sprintf("\nGene %s:\n\n%s", $gene_name, $gene_output);
			}
			if ($detailed){
				$gene_output .= "\n";
				$gene_output .= sprintf("\tProtein function:\n\t\t%s\n", $gene_function);
				$gene_output .= "\n";
				$gene_output .= sprintf("\tAssociated diseases: \n");
				foreach my $disease (@gene_diseases){
					$gene_output .= sprintf("\t\t%s\n", $disease);
				}
			}
		}
		$output .= $gene_output;
	}
	
	return $output;
	
	# 	if (defined($coverages)){
	# 		my %housekeeping_gene_coverages;
	# 		my @coverage_measures = ('mean', 'min', 'max', 'p90', 'p95');
	# 		foreach my $cov_measure (@coverage_measures){
	# 			$housekeeping_gene_coverages{$cov_measure} = mean(map $coverages->{$_}{$cov_measure}, @HOUSEKEEPING_GENES);
	# 		}
	# 		if ($coverages->{$ref_name}{'p95'} > 2*$housekeeping_gene_coverages{'p95'}){
	# 			printf("\tPotential overexpression. Coverage P95: %d (housekeeping transcripts: %d)\n",$coverages->{$ref_name}{'p95'},$housekeeping_gene_coverages{'p95'});
	# # 			printf("\tPotential overexpression. Coverage values: P95: %d, P90: %d, mean: %d, min: %d, max: %d\n",$coverages->{$ref_name}{'p95'},$coverages->{$ref_name}{'p90'},$coverages->{$ref_name}{'mean'},$coverages->{$ref_name}{'min'},$coverages->{$ref_name}{'max'});
	# 		}
	# 	}
}


#################################################################################












exit;
