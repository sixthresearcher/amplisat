#!/usr/bin/perl -w
#
################################################################
#
# Name:  reformat_amplitaxo.pl
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
#   Reformats FASTA/FASTQ files from taxonomy databases to AmpliTAXO database format (UTAX format)
#   AmpliTAXO (UTAX) header format:
#   >AB243007_S000622964;tax=d:Bacteria,p:Actinobacteria,c:Actinobacteria,o:Actinomycetales,f:Nocardiaceae,g:Smaragdicoccus;
#   k: kingdom, d:domain, p: phylum, c: class, o: order, f: family, g: genus, s: species
#
# Examples:
# perl  reformat_amplitaxo.pl -i SILVA.fa  -f silva -gz

# Modules are in folder 'lib' in the path of the script
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;
use warnings;
no warnings ('uninitialized', 'substr');

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

# Default options
# UTAX format fields
my @utax_fields = ('kingdom','domain','phylum','class','order','family','genus','species');
# Taxonomic universal format keys
my %taxo_keys = ('k'=>'kingdom','d'=>'domain','p'=>'phylum','c'=>'class','o'=>'order','f'=>'family','g'=>'genus','s'=>'species');
# Allele matching parameters
my $INP_allele_align = {'alignment' => 'dna blastn -evalue 1E-5 -ungapped -word_size 10 -perc_identity 90', 'aligned' => 0.9, 'ident' => 0.9 };
# rRNA and ITSs default databases
my $TAXO_DATABASES = {
	'oraldb' => { 'description' => 'A phylogenetically-curated 16S rDNA database of the core oral microbiome',
		'size' => 1262, 'sequence_file' => 'taxo/oralDB_012814.fasta.gz', 'version' => '20140128',
		'url' => 'http://microbiome.osu.edu', 'pubmed' => '21544197' },
	'homd' => { 'description' => 'The Human Oral Microbiome Database',
		'size' => 831, 'sequence_file' => 'taxo/HOMD_16S_rRNA_RefSeq_V13.2.fasta.gz', 'taxonomy_file' => 'taxo/homd_taxonomy_table.txt', 'version' => 'v13.2',
		'url' => 'http://www.homd.org', 'pubmed' => '20624719' },
	'greengenes' => { 'description' => 'Greengenes, a chimera-checked 16S rRNA gene database (only 90% identity core set)', # Greengenes clustered with cd-hit-est to 90% identity
		'size' => 4799, 'sequence_file' => 'taxo/greengenes_gg16S_90perc.fasta.gz', 'version' => '20110509',
		'url' => 'http://greengenes.lbl.gov', 'pubmed'=> '16820507' },
	'silva' => { 'description' => 'A comprehensive online resource for quality checked and aligned ribosomal RNA sequence data',
		'size' => 96642, 'sequence_file' => 'taxo/SILVA_123_LSURef_tax_silva.fasta.gz', 'version' => 'LSU Ref 123',
		'url' => 'http://www.arb-silva.de', 'pubmed'=> '23193283' },
# 	'silva' => { 'description' => '',
# 		'size' => 597607, 'sequence_file' => 'taxo/SILVA_123_SSURef_Nr99_tax_silva.fasta.gz', 'version' => 'SSU Ref NR99',
# 		'url' => 'http://www.arb-silva.de', 'pubmed'=> '23193283' },
	'unite' => { 'description' => 'Unified system for the DNA based fungal species linked to the classification',
		'size' => 22774, 'sequence_file' => 'taxo/sh_general_release_dynamic_01.08.2015.fasta.gz', 'version' => 'v7',
		'url' => 'https://unite.ut.ee', 'pubmed'=> '24112409' },
	'mothur' => { 'description' => 'MOTHUR RDP/PDS reference files',
		'size' => 10172, 'sequence_file' => 'trainset9_032012.pds.fasta', 'taxonomy_file' => 'trainset9_032012.pds.tax', 'version' => 'v1.38',
		'url' => 'http://www.mothur.org', 'pubmed'=> '19801464' },
	# QIIME separates in sequences and taxonomy files other formats
# 	'qiime' => { 'description' => 'QIIME OTUs files',
# 		'size' => 782906, 'sequence_file' => 'qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta', 'version' => 'v13.8',
# 		'url' => 'http://qiime.org', 'pubmed'=> '20383131' },
};

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

my ($INP_reads_file, $INP_taxonomy_file, $INP_format, $INP_outfile, $INP_nreads, $INP_shuffle, $INP_threads, $INP_gzip, $INP_zip);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s' => \$INP_reads_file,
	't|taxo=s' => \$INP_taxonomy_file,
	'f|format=s' => \$INP_format,
	'o|output=s' => \$INP_outfile,
	'n|number=i' => \$INP_nreads,
	's|shuffle' => \$INP_shuffle,
	'thr|threads=i' => \$INP_threads,
	'gz|gzip' => \$INP_gzip,
	'z|zip' => \$INP_zip,
	'<>' => \&usage,
);

# Usage help
sub usage
{
  print "\nusage: $0 -i READSFILE -d CVSFILE [-o OUTFILE] [options]\n\n";
  print " -h This message.\n";
  print " -i Input files in FASTA/FASTQ format (compressed or uncompressed).\n";
  print " -f Input format (".join(',',sort keys %$TAXO_DATABASES).").\n";
  print " -o Output filename.\n";
  print " -n Total number of reads to retrieve.\n";
  print " -s Shuffle and take random order reads.\n";
  print " -thr Number of threads to calculate the alignments.\n";
  print " -gz Compress results with GZIP.\n";
  print " -z Compress results with ZIP.\n";
  print "\n";
  #exit(-1);
  print "\nUnknown option: @_\n" if ( @_ );
#   print "\nusage: program [--url URL] [--size SIZE] [--help|-?]\n";
  #exit;
}

# Prints usage help if no input file is specified
if (!defined($INP_reads_file)){
	print "\nERROR: You must specify a FASTA/FASTQ input file.\n\n";
	usage();
	exit;
}
if (!defined($INP_format) || !in_array([keys %$TAXO_DATABASES],lc($INP_format))){
	print "\nERROR: You must specify a correct input format (".join(',',sort keys %$TAXO_DATABASES).").\n\n";
	usage();
	exit;
}
if ($INP_format =~ /homd/i && !defined($INP_taxonomy_file)) {
	print "\nERROR: You must specify a taxonomy file.\n\n";
	usage();
	exit;
}

# Creates name for the output file
if (!defined($INP_outfile)){
	if ($INP_reads_file =~ /(.+)\.(fa|fasta)/){
		$INP_outfile = $1.".utax";
	} elsif ($INP_reads_file =~ /(.+?)\./){
		$INP_outfile = $1.".utax";
	} else {
		$INP_outfile = $INP_reads_file.".utax";
	}
}

print "\nRunning '$COMMAND_LINE'\n";

# Check and read sequences file
my ($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads)
= parse_sequence_file($INP_reads_file,$INP_nreads,['verbose']);

my $id_to_taxo;
if (defined($INP_taxonomy_file) && -e $INP_taxonomy_file){
	# HOMD specific taxonomy info file
	if ($INP_format =~ /homd/i) {
		my $taxo_data = read_from_file($INP_taxonomy_file);
		my @col_names;
		foreach my $line (@$taxo_data){
			# Ex. 
			# HOT_ID  Domain  Phylum  Class   Order   Family  Genus   Species ...
			# 001     Bacteria        Proteobacteria  Alphaproteobacteria     Rhizobiales     Bartonellaceae  Bartonella      schoenbuchensis ...                                       
			if ($line =~ /^HOT_ID/){
				@col_names = split("\t",lc($line));
				next;
			} elsif ($line !~ /^\d+/) {
				next;
			}
			my @cols = split("\t",$line);
			for (my $j=1; $j<=7; $j++){
				if ($cols[$j] ne '') {
					$id_to_taxo->{$cols[0]}{$col_names[$j]}=$cols[$j];
				}
			}
		}
	# Standard taxonomy info file, with 2 tab separated columns: sequence ID and taxonomy data
	# Sometimes mothur format has sequences whose taxonomies are only annotated in an extra .tax file
	} else {
		my $taxo_data = read_from_file($INP_taxonomy_file);
		my @col_names;
		foreach my $line (@$taxo_data){
			# Ex. 
			# EU861894_S001148199	Bacteria;"Armatimonadetes";Armatimonadia;Armatimonadales;Armatimonadaceae;Armatimonas_Armatimonadetes_gp1;
			# EF115542_S001020530	Bacteria;Cyanobacteria_Chloroplast;Chloroplast;Chloroplast_order_incertae_sedis;Chloroplast;Streptophyta;
			my @cols = split("\t",$line);
			if (@cols) {
				$id_to_taxo->{$cols[0]} = $cols[1];
			}
		}
	}
}


my $count_seqs = 0;
my %previous_accessions;
my $acc_len = length($total_reads);
for (my $i=0; $i<=$#{$read_headers}; $i++) {

	$count_seqs++;

	# Hash to store taxo information, including UTAX fields
	my ($taxo_data_,%taxo_data);

	# If taxonomy data is given in an additional file
	if (defined($id_to_taxo)){
		my $seq_id = trim($read_headers->[$i]);
		if (defined($id_to_taxo->{$seq_id})) {
			$taxo_data{'accession'} = $seq_id;
			$taxo_data_ = $id_to_taxo->{$read_headers->[$i]};
		}
	}

	# Read the taxonomy information from the different formats
	if ($INP_format =~ /oraldb/i) {
		# >X80413; Actinomyces georgiae; georgiae; yes; Bacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Actinomycetaceae; Actinomyces;
		if (!defined($taxo_data_)){
			($taxo_data{'accession'}, $taxo_data{'otu_name'}, $taxo_data{'species'}, $taxo_data{'cultivated'}, $taxo_data{'kingdom'}, $taxo_data{'phylum'}, $taxo_data{'class'}, $taxo_data{'order'}, $taxo_data{'family'}, $taxo_data{'genus'}) = split(/\s?;\s?/,$taxo_data_);
		} else {
			($taxo_data{'otu_name'}, $taxo_data{'species'}, $taxo_data{'cultivated'}, $taxo_data{'kingdom'}, $taxo_data{'phylum'}, $taxo_data{'class'}, $taxo_data{'order'}, $taxo_data{'family'}, $taxo_data{'genus'}) = split(/\s?;\s?/,$taxo_data_);
		}
	} elsif ($INP_format =~ /silva/i) {
		# >JQ911630.2242.5087 Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Bradyrhizobiaceae;Bradyrhizobium;Bradyrhizobium elkanii
		if (!defined($taxo_data_) && $read_headers->[$i] =~ /(.+?)\.\d+\.\d+\s+(.+)/){
			$taxo_data{'accession'} = $1;
			$taxo_data_ = $2;
		}
		if (defined($taxo_data_)){
			($taxo_data{'kingdom'}, $taxo_data{'phylum'}, $taxo_data{'class'}, $taxo_data{'order'}, $taxo_data{'family'}, $taxo_data{'genus'}, $taxo_data{'species'}) = split(';',$taxo_data_);
			if (defined($taxo_data{'species'}) && $taxo_data{'species'} =~ /${taxo_data{'genus'}}\s+(.+)/ && $1 ne 'unidentified') {
				$taxo_data{'species'} = $1;
			} elsif (defined($taxo_data{'species'}) && $taxo_data{'species'} =~ /(.+?)\s+(.+)/) {
				$taxo_data{'genus'} = $1;
				$taxo_data{'species'} = $2;
			}
		}
		# Replace 'U' by 'T' in the sequences
		$read_seqs->[$i] =~ s/U/T/ig;
	} elsif ($INP_format =~ /homd/i) {
		# >001A28SC; Bartonella sp.; HOT_001; Strain_A28SC; GQ422708; Unnamed
		# HOMD format requires taxonomy data file and particular $id_to_taxo
		my ($id_seq, $organism, $taxo_id, $clone, $accession, $status) = split(/\s?;\s?/,$taxo_data_);
		if ($taxo_id =~/HOT_(\w+)/){
			%taxo_data = %{$id_to_taxo->{$1}};
			$taxo_data{'accession'} = $accession;
		}
	} elsif ($INP_format =~ /greengenes/i) {
		# >M32222.1; taxo=Methanothermus fervidus k__Archaea; p__Euryarchaeota; c__Methanobacteria; o__Methanobacteriales; Unclassified; otu_120
		if (!defined($taxo_data_) && $read_headers->[$i] =~ /(.+?);\s?(taxo=.+)/){
			$taxo_data{'accession'} = $1;
			$taxo_data_ = $2;
		}
		if (defined($taxo_data_) && $taxo_data_ =~ /taxo=(.+?) (k__.+)/){
			my @species = split(" ",$1);
			$taxo_data_ = $2;
			if ($#species == 1) {
				$taxo_data{'genus'} = $species[0];
				$taxo_data{'species'} = $species[1];
			}
		}
		if (defined($taxo_data_)){
			$taxo_data_ =~ s/["']//g;
			foreach my $taxo_field (split(/;\s?/,$taxo_data_)){
				if ($taxo_field =~ /(\w__)?unidentified/i) {
					last;
				} elsif ($taxo_field =~ /k__(.+)/) {
					$taxo_data{'kingdom'} = $1;
				} elsif ($taxo_field =~ /p__(.+)/) {
					$taxo_data{'phylum'} = $1;
				} elsif ($taxo_field =~ /c__(.+)/) {
					$taxo_data{'class'} = $1;
				} elsif ($taxo_field =~ /o__(.+)/) {
					$taxo_data{'order'} = $1;
				} elsif ($taxo_field =~ /f__(.+)/) {
					$taxo_data{'family'} = $1;
				} elsif (!defined($taxo_data{'genus'}) && $taxo_field =~ /g__(.+)/) {
					$taxo_data{'genus'} = $1;
				} elsif (!defined($taxo_data{'species'}) && $taxo_field =~ /s__${taxo_data{'genus'}}[\s_](.+)/) {
					$taxo_data{'species'} = $1;
				} elsif (!defined($taxo_data{'species'}) && $taxo_field =~ /s__(.+)/) {
					$taxo_data{'species'} = $1;
				}
			}
		}
	} elsif ($INP_format =~ /unite/i) {
		# >Sordariomycetes_sp|FJ613078|SH174117.07FU|reps|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__unidentified;f__unidentified;g__unidentified;s__Sordariomycetes_sp
		if (!defined($taxo_data_)){
			my @data = split('\|',$taxo_data_);
			$taxo_data{'accession'} = $data[1];
			$taxo_data_ = $data[4];
		}
		if (defined($taxo_data_)){
			$taxo_data_ =~ s/["']//g;
			foreach my $taxo_field (split(/;\s?/,$taxo_data_)){
				if ($taxo_field =~ /(\w__)?unidentified/i) {
					last;
				} elsif ($taxo_field =~ /k__(.+)/) {
					$taxo_data{'kingdom'} = $1;
				} elsif ($taxo_field =~ /p__(.+)/) {
					$taxo_data{'phylum'} = $1;
				} elsif ($taxo_field =~ /c__(.+)/) {
					$taxo_data{'class'} = $1;
				} elsif ($taxo_field =~ /o__(.+)/) {
					$taxo_data{'order'} = $1;
				} elsif ($taxo_field =~ /f__(.+)/) {
					$taxo_data{'family'} = $1;
				} elsif ($taxo_field =~ /g__(.+)/) {
					$taxo_data{'genus'} = $1;
				} elsif ($taxo_field =~ /s__${taxo_data{'genus'}}[\s_](.+)/) {
					$taxo_data{'species'} = $1;
				} elsif ($taxo_field =~ /s__(.+)/) {
					$taxo_data{'species'} = $1;
				}
			}
		}
	} elsif ($INP_format =~ /mothur/i) {
		# >EU861894_S001148199	Root;Bacteria;"Armatimonadetes";Armatimonadia;Armatimonadales;Armatimonadaceae;Armatimonas/Armatimonadetes_gp1
		if (!defined($taxo_data_) && $read_headers->[$i] =~ /(.+?)\s(Root;)?(.+)/){
			$taxo_data{'accession'} = $1;
			$taxo_data_ = $3;
		}
		if (defined($taxo_data_)){
			$taxo_data_ =~ s/["']//g;
			($taxo_data{'kingdom'}, $taxo_data{'phylum'}, $taxo_data{'class'}, $taxo_data{'order'}, $taxo_data{'family'}, $taxo_data{'genus'}, $taxo_data{'species'}) = split(/;/,$taxo_data_);
			if (%taxo_data){
				my $remove_uncertain = 0;
				foreach my $taxo ( ('kingdom','phylum','class','order','family','genus','species') ) {
					if ($remove_uncertain){
						delete($taxo_data{$taxo});
					} elsif ($taxo_data{$taxo} =~ /incertae/i) {
						$remove_uncertain = 1;
						delete($taxo_data{$taxo});
					}
				}
			}
		}
	}

	# Creates the header in UTAX format
	my @utax_data;
	foreach my $utax_field (@utax_fields){
		if (defined($taxo_data{$utax_field})){
			# Skip general unidentified fields
			if ($taxo_data{$utax_field} =~ /unidentified/i || $taxo_data{$utax_field} =~ /^sp\.?$/i){
				next;
			}
			push(@utax_data, sprintf("%s:%s",substr($utax_field,0,1),$taxo_data{$utax_field}));
		}
	}
	if (!@utax_data){
		if (defined($taxo_data_)){
			print "\tERROR: '$INP_format' format doesn't match the taxonomy data: '$taxo_data_'.\n";
		} else {
			print "\tERROR: '$INP_format' format doesn't match the taxonomy data: '".$read_headers->[$i]."'.\n";
		}
		splice(@{$read_headers},$i,1);
		splice(@{$read_seqs},$i,1);
		$i--;
	}
	
	# Redefines header
	if (defined($taxo_data{'accession'})){
		if (!defined($previous_accessions{$taxo_data{'accession'}})){
			$read_headers->[$i] = sprintf("%s;tax=%s;", $taxo_data{'accession'}, join(',',@utax_data));
		} else {
			$read_headers->[$i] = sprintf("%s.%d;tax=%s;", $taxo_data{'accession'}, $previous_accessions{$taxo_data{'accession'}}+1, join(',',@utax_data));
		}
		$previous_accessions{$taxo_data{'accession'}}++;
	} else {
		$read_headers->[$i] = sprintf("%0${acc_len}d;tax=%s;", $count_seqs, join(',',@utax_data));
	}
}

my $outfile;
if ($reads_file_format eq 'fastq'){
	$outfile = create_fastq_file($read_seqs,$read_headers,$read_qualities,"$INP_outfile.fq");
} else {
	$outfile = create_fasta_file($read_seqs,$read_headers,"$INP_outfile.fa");
}
if (defined($INP_zip) && defined($outfile)){
	`zip -jqm $outfile.zip $outfile` ;
	printf("\nReformatted %d sequences into '%s'.\n\n", scalar @$read_seqs, "$outfile.zip");
} elsif (defined($INP_gzip) && defined($outfile)){
	`gzip -f $outfile` ;
	printf("\nReformatted %d sequences into '%s'.\n\n", scalar @$read_seqs, "$outfile.gz");
} elsif (defined($outfile)){
	printf("\nReformatted %d sequences into '%s'.\n\n", scalar @$read_seqs, $outfile);
} else {
	print "\nThere was some error in the input file or format and no sequences were reformatted.\n\n";
}


exit;













