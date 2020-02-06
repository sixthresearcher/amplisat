#!/usr/bin/perl -w
#
# Name: ampliCANCER.pl
#
# Version: 1.0
#
# License: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
#          CC BY-NC-SA 4.0 (http://creativecommons.org/licenses/by-nc-sa/4.0/)
#
# Author: Alvaro Sebastian
#
# Support: Alvaro Sebastian (sixthresearcher@gmail.com)
#
# Sixth Researcher - www.sixthresearcher.com
#
# Description:
#   Analyzes raw data from Amplicon Sequencing Targeted Cancer Gene Panels and uses qualitative genotyping to detect clinically relevant variants in the human DNA/RNA from Amplicon Sequencing data for the purpose of reporting and interpreting main genetic alterations in the tumor cells.
#
# Example data:
# 1. Frequent mutations in human lung cancer tumors detected by Ion Torrent DNA sequencing. In order to identify genetic mutations in individual human lung cancer, we sequenced 13,500 loci from 45 cancer-related genes in 82 human lung cancer samples using the Ion Torrent Ampliseq Cancer Panel. The sequencing analysis revealed frequent missense mutations, unique combination mutations, and hotspot loci in the genes BRAF, EGFR, K-ras, PIK3CA, and TP53 in the human lung cancer samples of various histologic types.
#    http://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0095228
#    http://sra.dbcls.jp/search/view/SRP028756
#    http://www.ncbi.nlm.nih.gov/Traces/sra/?study = SRP028756
# 2. Exome capture was performed on the normal mucosa, adenoma, and adenocarcinoma tissues from the same patient by using NimbleGen 2.1 M Human Exome Array. The target regions of exome capture include 180,000 coding exon (28.4 Mb) and adjacent regions (5.7 M). Combining parallel sequencing on an Illumina GAII platform, we generated approximately 4.6, 4.4 and 4.3 billion bases of effective sequence data with an average read length of 69 bases for normal mucosa, adenoma and carcinoma, respectively.
#    http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0053310
#    https://trace.ddbj.nig.ac.jp/DRASearch/submission?acc=SRA052805
# 3. Mate-pair and targeted exome sequencing of primary and metastatic colorectal cancer of four patients
#    https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-10-r103
#    http://sra.dbcls.jp/search/view/ERP000875
#    
# 
# Requires as input a FASTA or FASTQ file with sequences/reads from one sample/individual
#
# Examples:
#    perl ampliCANCER.pl -i SRR953273.fastq.gz -f html -d panel -c lung
#    perl ampliCANCER.pl -i SRR953299.fastq.gz -d panel -c lung 


my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Analyzes raw data from Amplicon Sequencing Targeted Cancer Gene Panels by genotyping amplicon sequencing data with AmpliSAS clustering algorithm and reports aberrant variants vs. protein references from UniProt database.";

# Modules are in folder 'lib' in the path of the script
use lib "lib";
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;
use Bio::Onco;
use Spreadsheet::XLSX;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Sort::Naturally;
use Data::Dumper;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $SCRIPT_DIR = dirname(__FILE__);
my $COMMAND_LINE = $0." ".join(" ",@ARGV);
my $MAKEBLASTDBEXE = $SCRIPT_DIR.'/lib/Bio/tools/makeblastdb';

# Default options:
# Output format
my $INP_outformat = 'text short';
# Analysis parameters
my %DEFAULT_PARAMS_DNA = (
	'reference_data_folder' => "/home/alvaro/db/ensembl/Homo_sapiens.GRCh38",
	'reference_sequences_file' => "/home/alvaro/db/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
	'reference_sequences_url' => "ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
	'reference_annotations_file' => "/home/alvaro/db/ensembl/Homo_sapiens.GRCh38.89.gff3.gz",
	'reference_annotations_url' => "ftp://ftp.ensembl.org/pub/release-89/gff3/homo_sapiens/Homo_sapiens.GRCh38.89.gff3.gz",
# 	'reference_sequences_file' => "/home/alvaro/db/ensembl/Homo_sapiens.GRCh38.dna.chromosome.4.fa.gz",
# 	'reference_annotations_file' => "/home/alvaro/db/ensembl/Homo_sapiens.GRCh38.89.chromosome.4.gff3.gz",
	# Ensembl human reference gene sequences file (only introns and exons from the beginning to the end of the mRNA region)
	# This file is generated with the script 'extract_genome_seqs.pl' using as input the Primary Assembly genome and GFF3 file with annotations (extract_genome_seqs.pl -i Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -g Homo_sapiens.GRCh38.89.gff3.gz -t genecod -gz)
	#'reference_genes_file' => "/home/alvaro/db/ensembl/Homo_sapiens.premrna.fasta.gz",
	# Minimum frequency to annotate/discard low frequency variant mutations
	'min_variant_frequency' => 10, # 10%
	# Minimum coverage to annotate/discard low coverage variant mutations
	'min_variant_coverage' => 10,

);
my %DEFAULT_PARAMS_RNA = (
	# Ensembl human reference CDs/transcripts sequences file (http://www.ensembl.org/info/data/ftp/index.html)
	'reference_sequences_file' => "/home/alvaro/db/ensembl/Homo_sapiens.GRCh38.cds.all.fa",
	'reference_sequences_url' => "ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz",
# 	'reference_sequences_file' => $SCRIPT_DIR.'/human/Homo_sapiens.GRCh38.cds.all.fa.gz',
	# Minimum frequency to annotate/discard low frequency variant mutations
	'min_variant_frequency' => 1, # 1%
	# Minimum coverage to annotate/discard low coverage variant mutations
	'min_variant_coverage' => 3,
);
my %DEFAULT_PARAMS = (
	# Folder with UniProt reference sequences and data from oncogenes, tumor suppressors...
	# It incorporates any new protein found in new data
	'oncodata_folder' => "/home/alvaro/www/amplisat/bin/oncodata",
# 	'oncodata_folder' => $SCRIPT_DIR.'/oncodata',
	# COSMIC annotations file (coding mutations in genes listed in the Cancer Gene Census)
	'mutation_data_folder' => "/home/alvaro/db/cosmic",
	'mutation_annotations_file' => "/home/alvaro/db/cosmic/CosmicMutantExportCensus.tsv.gz",
	'mutation_annotations_url' => "/files/grch38/cosmic/v82/CosmicMutantExportCensus.tsv.gz",
	# Uniprot human reference protein sequences file
	# 'uniprot_file' => "/home/alvaro/db/uniprot/uniprot_human_refs.fa.gz",
	# Uniprot ID mapping information file, with Ensembl IDs equivalences
	'idmapping_file' => "/home/alvaro/db/uniprot/HUMAN_9606_idmapping.dat.gz",
	# Mapping and variant calling parameters
	'align_params' => 'bowtie2 --sensitive-local -k 2',
# 	'align_params' => 'bowtie2 --sensitive-local',
	# 'mpileup_params' => 'mpileup -Dgf /home/alvaro/db/ensembl/Homo_sapiens.GRCh38.cds.all.fa',
	# 'mpileup_params' => "mpileup -C50 -D -g -f /home/alvaro/db/ensembl/Homo_sapiens.GRCh38.cds.all.fa",
	# 'mpileup_params' => "samtools mpileup -g -f /home/alvaro/db/ensembl/Homo_sapiens.GRCh38.cds.all.fa",
	# 'bcftools_params' => 'view -vcg',
	# Highlight cancer type specific variants
	'cancer_type' => 'any',
	# Max. number of allowed errors in mapping
	'max_errors' => 2,
	# Annotate insertions and deletions
	'indels' => 1,
	# Keep or discard rare variant mutations
	'rare_variants' => 0,
	# Minimum distance in bps to the start/end of the read to annotate the mutations
	'min_variant_position' => 2,
# 	# Reference matching parameters
# 	'blastn_params' => 'dna blastn-short -qcov_hsp_perc 100',
#  	'blastn_params' => 'dna blastn -evalue 1E-5 -ungapped',
# 	'blastn_params' => 'dna blastn-short -evalue 1E-3',
# 	# Minimum % of sequence aligned to a reference
# 	'min_allele_align' => 100,
# 	# Minimum % of identity in the portion aligned to a reference
# 	'min_allele_ident' => 100,
# 	# Amplicons with lower total coverage will be discarded
# 	'min_amplicon_depth' => 100,
# 	# Minimum frequency of variants after clustering
# 	'min_amplicon_seq_frequency' => 0.002, # 0.005%
# 	# Variants with equal or higher identity will be clustered together
# 	'identity_threshold' => undef,
# 	# Keep singletons after clustering
# 	'keep_singletons' => undef,
);
# To print the options
my %yes_no = ( 0 => 'no', 1 => 'yes' );
# Important genes per tumor type (Sources: https://www.mycancergenome.org and literature)
my $CANCER_GENES = {
	'lung' => [ 'AKT1', 'ALK', 'BRAF', 'CD274', 'DDR2', 'EGFR', 'FGRF1', 'FGFR3', 'HER2', 'KRAS', 'MAP2K1', 'MEK1', 'MET', 'NRAS', 'NTRK1', 'P16', 'PIK3CA', 'PTEN', 'RET', 'RICTOR', 'ROS1', 'TP53', 'VEGF', 'VEGFR1', 'VEGFR2' ],
	'breast' => [ 'AKT1', 'AR', 'BRCA1', 'BRCA2', 'CCND1', 'CDK4', 'CDK6', 'ERBB2', 'ESR1','FGFR1', 'FGFR2', 'PARP', 'PGR', 'PIK3CA', 'PTEN', 'RB1', 'TP53' ], # HER2=ERBB2
	'colorectal' => [ 'AKT1', 'APC', 'ATRX', 'BRAF', 'DDR2', 'DNAH9', 'EPHA4', 'ERBB3', 'EXOC4', 'FGFR2', 'KDR', 'KR​AS', 'MLL3', 'NRAS', 'NUP98', 'PARP', 'PIK3CA', 'PRKCD', 'PTEN', 'PTPRF', 'RASA3', 'RFC1', 'SMAD2', 'SMAD4', 'TAOK1', 'TGFBR2', 'TP53', 'TTN' ],
};
# Important genes per panel type (several sources)
my $PANEL_GENES = {
	'foundationone' => ['ABL1','ABL2','ACVR1B','AKT1','AKT2','AKT3','ALK','AMER1','APC','AR','ARAF','ARFRP1','ARID1A','ARID1B','ARID2','ASXL1','AT M','AT R','ATRX','AURKA','AURKB','AXIN1','AXL','BAP1','BARD1','BCL2','BCL2L1','BCL2L2','BCL6','BCOR','BCORL1','BCR','BLM','BRAF','BRCA1','BRCA2','BRD4','BRIP1','BTG1','BTK','C17orf39','CARD11','CBFB','CBL','CCND1','CCND2','CCND3','CCNE1','CD274','CD79A','CD79B','CDC73','CDH1','CDK12','CDK4','CDK6','CDK8','CDKN1A','CDKN1B','CDKN2A','CDKN2B','CDKN2C',
				'CEBPA','CHD2','CHD4','CHEK1','CHEK2','CIC','CREBBP','CRKL','CRLF2','CSF1R','CTCF','CTNNA1','CTNNB1','CUL3','CYLD','DAXX','DDR2','DICER1','DNMT3A','DOT1L','EGFR','EMSY','EP300','EPHA3','EPHA5','EPHA7','EPHB1','ERBB2','ERBB3','ERBB4','ERG','ERRFl1','ESR1','ETV1','ETV4','ETV5','ETV6','EZH2','FAM123B','FAM46C','FANCA','FANCC','FANCD2','FANCE','FANCF','FANCG','FANCL','FAS','FAT 1','FBXW7','FGF10','FGF14','FGF19','FGF23','FGF3','FGF4','FGF6','FGFR1','FGFR2','FGFR3','FGFR4',
				'FH','FLCN','FLT1','FLT3','FLT4','FOXL2','FOXP1','FRS2','FUBP1','GABRA6','GATA1','GATA2','GATA3','GATA4','GATA6','GID4','GLI1','GNA11','GNA13','GNAQ','GNAS','GPR124','GRIN2A','GRM3','GSK3B','H3F3A','HGF','HNF1A','HRAS','HSD3B1','HSP90AA1','IDH1','IDH2','IGF1R','IGF2','IKBKE','IKZF1','IL7R','INHBA','INPP4B','IRF2','IRF4','IRS2','JAK1','JAK2','JAK3','JUN','KAT6A','KDM5A','KDM5C','KDM6A','KDR','KEAP1','KEL','KIT','KLHL6','KMT2A','KMT2C','KMT2d','KRAS','LMO1','LRP1B',
				'LYN','LZTR1','MAGI2','MAP2K1','MAP2K2','MAP2K4','MAP3K1','MCL1','MDM2','MDM4','MED12','MEF2B','MEN1','MET','MITF','MLH1','MLL','MLL2','MLL3','MPL','MRE11A','MSH2','MSH6','MTOR','MUTYH','MYB','MYC','MYCL','MYCL1','MYCN','MYD88','MYST3','NF1','NF2','NFE2L2','NFKBIA','NKX2-1','NOTCH1','NOTCH2','NOTCH3','NPM1','NRAS','NSD1','NTRK1','NTRK2','NTRK3','NUP93','PAK3','PALB2','PARK2','PAX5','PBRM1','PDCD1LG2','PDGFRA','PDGFRB','PDK1','PIK3C2B','PIK3CA','PIK3CB','PIK3CG','PIK3R1',
				'PIK3R2','PLCG2','PMS2','POLD1','POLE','PPP2R1A','PRDM1','PREX2','PRKAR1A','PRKCI','PRKDC','PRSS8','PTCH1','PTEN','PTPN11','QKI','RAC1','RAD50','RAD51','RAF1','RANBP2','RARA','RB1','RBM10','RET','RICTOR','RNF43','ROS1','RPTOR','RUNX1','RUNX1T1','SDHA','SDHB','SDHC','SDHD','SETD2','SF3B1','SLIT2','SMAD2','SMAD3','SMAD4','SMARCA4','SMARCB1','SMO','SNCAIP','SOCS1','SOX10','SOX2','SOX9','SPEN','SPOP','SPTA1','SRC','STAG2','STAT 4','STAT3','STK11','SUFU','SYK','TAF1','TBX3',
				'TERC','TERT','TET2','TGFBR2','TMPRSS2','TNFAIP3','TNFRSF14','TOP1','TOP2A','TP53','TSC1','TSC2','TSHR','U2AF1','VEGFA','VHL','WISP3','WT1','XPO1','ZBTB2','ZNF217','ZNF703'],
	'foundationact' => ['ABL1','AKT1','ALK','ARAF','BRAF','BRCA1','BRCA2','BTK','CCND1','CD274','CDH1','CDK4','CDK6','CDKN2A','CRKL','CTNNB1','DDR2','EGFR','ERBB2','ERRFI1','ESR1','EZH2','FGFR1','FGFR2','FGFR3','FLT3','FOXL2','GNA11','GNAQ','GNAS','HRAS','IDH1','IDH2','JAK2','JAK3','KIT','KRAS','MAP2K1','MAP2K2','MDM2','MET','MPL','MTOR','MYC','MYCN','MYD88','NF1','NPM1','NRAS','PDCD1LG2','PDGFRA','PDGFRB','PIK3CA','PTEN','PTPN11','RAF1','RET','ROS1','SMO','TERT','TP53','VEGFA'],
	'oncominefocusassay' => ['ABL1','AKT1','AKT3','ALK','AR','AXL','BRAF','CCND1','CDK4','CDK6','CTNNB1','DDR2','EGFR','ERBB2','ERBB3','ERBB4','ERG','ESR1','ETV1','ETV4','ETV5','FGFR1','FGFR2','FGFR3','FGFR4','GNA11','GNAQ','HRAS','IDH1','IDH2','JAK1','JAK2','JAK3','KIT','KRAS','MAP2K1','MAP2K2','MET','MTOR','MYC','MYCN','NRAS','NTRK1','NTRK2','NTRK3','PDGFRA','PIK3CA','PPARG','RAF1','RET','ROS1','SMO'],
	'oncominecomprehensiveassayv3' => ['AKT1','AKT2','AKT3','ALK','AR','ARAF','ARID1A','ATM','ATR','ATRX','AXL','BAP1','BRAF','BRCA1','BRCA2','BTK','CBL','CCND1','CCND2','CCND3','CCNE1','CDK12','CDK2','CDK4','CDK6','CDKN1B','CDKN2A','CDKN2B','CHEK1','CHEK2','CREBBP','CSF1R','CTNNB1','DDR2','EGFR','ERB83','ERB84','ERBB2','ERBB4','ERCC2','ERG','ESR1','ETV1','ETV4','ETV5','EZH2','FANCA','FANCD2','FANCI','FBXW7','FGF19','FGF3','FGFR1','FGFR2','FGFR3','FGFR4','FGR','FLT3','FOXL2','GATA2','GNA11',
						'GNAQ ','GNAS','H3F3A','HIST1H3B','HNF1A','HRAS','IDH1','IDH2','IGF1R','JAK1','JAK2','JAK3','KDR','KIT','KNSTRN','KRAS','MAGOH','MAP2K1','MAP2K2','MAP2K4','MAPK1','MAX','MDM2','MDM4','MED12','MET','MLH1','MRE11A','MSH2','MSH6','MTOR','MYB','MYBL1','MYC','MYCL','MYCN','MYD88','NBN','NF1','NF2','NFE2L2','NOTCH1','NOTCH2','NOTCH3','NOTCH4','NRAS','NRG1 ','NTRK1','NTRK2','NTRK3','NUTM1','PALB2','PDGFRA','PDGFRB','PIK3CA','PIK3CB','PIK3R1','PMS2 ','POLE',
						'PPARG','PPP2R1A ','PRKACA','PRKACB','PTCH1','PTEN','PTEN ','PTPN11','RAC1','RAD50','RAD51','RAD51B','RAD51C','RAD51D','RAF1','RB1','RELA','RET','RHEB','RHOA','RICTOR','RNF43 ','ROS1','ROS1 ','RSPO2','RSPO3','SETD2','SF3B1','SLX4','SMAD4','SMARCA4','SMARCB1','SMO','SPOP','SRC','STAT3','STK11','TERT','TERT ','TOP1','TP53','TSC1','TSC2','TSC2 ','U2AF1','XPO1'],
	'truseqamplicon' => ['ABL1','AKT1','ALK','APC','ATM','BRAF','CDH1','CDKN2A','CSF1R','CTNNB1','EGFR','ERBB2','ERBB4','FBXW7','FGFR1','FGFR2','FGFR3','FLT3','GNA11','GNAQ','GNAS','HNF1A','HRAS','IDH1','JAK2','JAK3','KDR','KIT','KRAS','MET','MLH1','MPL','NOTCH1','NPM1','NRAS','PDGFRA','PIK3CA','PTEN','PTPN11','RB1','RET','SMAD4','SMARCB1','SMO','SRC','STK11','TP53','VHL'],
	'trusightcancer' => ['AIP','ALK','APC','ATM','BAP1','BLM','BMPR1A','BRCA1','BRCA2','BRIP1','BUB1B','CDC73','CDH1','CDK4','CDKN1C','CDKN2A','CEBPA','CEP57','CHEK2','CYLD','DDB2','DICER1','DIS3L2','EGFR','EPCAM','ERCC2','ERCC3','ERCC4','ERCC5','EXT1','EXT2','EZH2','FANCA','FANCB','FANCC','FANCD2','FANCE','FANCF','FANCG','FANCI','FANCL','FANCM','FH','FLCN','GATA2','GPC3','HNF1A','HRAS','KIT','MAX','MEN1','MET','MLH1','MSH2','MSH6','MUTYH','NBN','NF1','NF2','NSD1','PALB2','PHOX2B','PMS1','PMS2',
				'PRF1','PRKAR1A','PTCH1','PTEN','RAD51C','RAD51D','RB1','RECQL4','RET','RHBDF2','RUNX1','SBDS','SDHAF2','SDHB','SDHC','SDHD','SLX4','SMAD4','SMARCB1','STK11','SUFU','TMEM127','TP53','TSC1','TSC2','VHL','WRN','WT1','XPA','XPC'],
	'trusightrnapancancer' => ['ABCC3','ABI1','ABL1','ABL2','ABLIM1','ACACA','ACE','ACER1','ACKR3','ACSBG1','ACSL3','ACSL6','ACVR1B','ACVR1C','ACVR2A','ADD3','ADM','AFF1','AFF3','AFF4','AGR3','AHCYL1','AHI1','AHR','AHRR','AIP','AK2','AK5','AKAP12','AKAP6','AKAP9','AKR1C3','AKT1','AKT2','AKT3','ALDH1A1','ALDH2','ALDOC','ALK','AMER1','AMH','ANGPT1','ANKRD28','ANLN','APC','APH1A','APLP2','APOD','AR','ARAF','ARFRP1','ARHGAP20','ARHGAP26','ARHGEF12','ARHGEF7','ARID1A','ARID2','ARIH2','ARNT','ARRDC4',
					'ASMTL','ASPH','ASPSCR1','ASTN2','ASXL1','ATF1','ATF3','ATG13','ATG5','ATIC','ATL1','ATM','ATP1B4','ATP8A2','ATR','ATRNL1','ATRX','AURKA','AURKB','AUTS2','AXIN1','AXL','BACH1','BACH2','BAG4','BAIAP2L1','BAP1','BARD1','BAX','BAZ2A','BCAS3','BCAS4','BCL10','BCL11A','BCL11B','BCL2','BCL2A1','BCL2L1','BCL2L2','BCL3','BCL6','BCL7A','BCL9','BCOR','BCORL1','BCR','BDNF','BHLHE22','BICC1','BIN1','BIRC3','BIRC6','BLM','BMP4','BMPR1A','BRAF','BRCA1','BRCA2','BRD1','BRD3',
					'BRD4','BRIP1','BRSK1','BRWD3','BTBD18','BTG1','BTG2','BTK','BTLA','BUB1B','C11orf1','C11orf30','C11orf54','C11orf95','C2CD2L','C2orf44','C3orf27','CACNA1F','CACNA1G','CACNA2D3','CAD','CALR','CAMK2A','CAMK2B','CAMK2G','CAMTA1','CANT1','CAPRIN1','CAPZB','CARD11','CARM1','CARS','CASC5','CASP3','CASP7','CASP8','CAV1','CBFA2T3','CBFB','CBL','CBLB','CBLC','CCAR2','CCDC28A','CCDC6','CCDC88C','CCK','CCL2','CCNA2','CCNB1IP1','CCNB3','CCND1','CCND2','CCND3','CCNE1',
					'CCNG1','CCT6B','CD19','CD22','CD274','CD28','CD36','CD44','CD58','CD70','CD74','CD79A','CD79B','CD8A','CDC14A','CDC14B','CDC25A','CDC25C','CDC42','CDC73','CDH1','CDH11','CDK1','CDK12','CDK2','CDK4','CDK5RAP2','CDK6','CDK7','CDK8','CDK9','CDKL5','CDKN1A','CDKN1B','CDKN1C','CDKN2A','CDKN2B','CDKN2C','CDKN2D','CDX1','CDX2','CEBPA','CEBPB','CEBPD','CEBPE','CENPF','CENPU','CEP170B','CEP57','CEP85L','CHCHD7','CHD2','CHD6','CHEK1','CHEK2','CHIC2','CHL1','CHMP2B',
					'CHN1','CHST11','CHUK','CIC','CIITA','CIRH1A','CIT','CKB','CKS1B','CLP1','CLTA','CLTC','CLTCL1','CMKLR1','CNBP','CNOT2','CNTN1','CNTRL','COG5','COL11A1','COL1A1','COL1A2','COL3A1','COL6A3','COL9A3','COMMD1','COX6C','CPNE1','CPS1','CPSF6','CRADD','CREB1','CREB3L1','CREB3L2','CREBBP','CRKL','CRLF2','CRTC1','CRTC3','CSF1','CSF1R','CSF3','CSF3R','CSNK1G2','CSNK2A1','CTCF','CTDSP2','CTLA4','CTNNA1','CTNNB1','CTNND2','CTRB1','CTSA','CUX1','CXCL8','CXCR4','CXXC4',
					'CYFIP2','CYLD','CYP1B1','CYP2C19','DAB2IP','DACH1','DACH2','DAXX','DCLK2','DCN','DDB2','DDIT3','DDR2','DDX10','DDX20','DDX39B','DDX3X','DDX5','DDX6','DEK','DGKB','DGKI','DGKZ','DICER1','DIRAS3','DIS3L2','DKK1','DKK2','DKK4','DLEC1','DLL1','DLL3','DLL4','DMRT1','DMRTA2','DNAJB1','DNM1','DNM2','DNM3','DNMT1','DNMT3A','DOCK1','DOT1L','DPM1','DPYD','DST','DTX1','DTX4','DUSP2','DUSP22','DUSP26','DUSP9','DUX4','E2F1','EBF1','ECT2L','EDIL3','EDNRB','EED','EEFSEC',
					'EGF','EGFR','EGR1','EGR2','EGR3','EGR4','EIF4A2','EIF4E','ELF4','ELK4','ELL','ELN','ELOVL2','ELP2','EML1','EML4','ENPP2','EP300','EP400','EPC1','EPCAM','EPHA10','EPHA2','EPHA3','EPHA5','EPHA7','EPHB1','EPHB6','EPO','EPOR','EPS15','ERBB2','ERBB3','ERBB4','ERC1','ERCC1','ERCC2','ERCC3','ERCC4','ERCC5','ERCC6','ERG','ERLIN2','ESR1','ETS1','ETS2','ETV1','ETV4','ETV5','ETV6','EWSR1','EXOSC6','EXT1','EXT2','EYA1','EYA2','EZH2','EZR','FAF1','FAM127C','FAM19A2',
					'FAM19A5','FAM46C','FAM64A','FANCA','FANCB','FANCC','FANCD2','FANCE','FANCF','FANCG','FANCI','FANCL','FANCM','FAS','FASLG','FBN2','FBXO11','FBXO31','FBXW7','FCGBP','FCGR2B','FCRL4','FEN1','FEV','FGF1','FGF10','FGF13','FGF14','FGF19','FGF2','FGF23','FGF3','FGF4','FGF6','FGF8','FGF9','FGFR1','FGFR1OP','FGFR1OP2','FGFR2','FGFR3','FGFR4','FH','FHIT','FHL2','FIGF','FIP1L1','FLCN','FLI1','FLNA','FLNC','FLT1','FLT3','FLT3LG','FLT4','FLYWCH1','FNBP1','FOS','FOSB',
					'FOSL1','FOXL2','FOXO1','FOXO3','FOXO4','FOXP1','FRK','FRMPD4','FRS2','FRYL','FSTL3','FUS','FUT1','FZD10','FZD2','FZD3','FZD6','FZD7','FZD8','GAB1','GABRG2','GADD45B','GANAB','GAS1','GAS5','GAS7','GATA1','GATA2','GATA3','GATA6','GBP2','GDF6','GFAP','GHR','GID4','GIT2','GLI1','GLI3','GMPS','GNA11','GNA12','GNA13','GNAI1','GNAQ','GNAS','GNG4','GOLGA5','GOPC','GOSR1','GOT1','GPC3','GPHN','GPR124','GPR128','GPR34','GRB10','GRB2','GRHPR','GRID1','GRIN2A','GRIN2B',
					'GRM1','GRM3','GSK3B','GSN','GSTT1','GTF2I','GTSE1','H2AFX','H3F3A','HAS2','HDAC1','HDAC2','HDAC3','HDAC4','HDAC5','HDAC6','HDAC7','HECW1','HEPH','HERPUD1','HES1','HES5','HEY1','HGF','HHEX','HIF1A','HIP1','HIPK1','HIPK2','HIST1H1C','HIST1H1D','HIST1H1E','HIST1H2AC','HIST1H2AG','HIST1H2AL','HIST1H2AM','HIST1H2BC','HIST1H2BJ','HIST1H2BK','HIST1H2BO','HIST1H3B','HIST1H4I','HLF','HMGA1','HMGA2','HMGB1','HMGN2P46','HNF1A','HNRNPA2B1','HOOK3','HOXA10','HOXA11',
					'HOXA13','HOXA3','HOXA9','HOXC11','HOXC13','HOXD11','HOXD13','HOXD9','HRAS','HSP90AA1','HSP90AB1','HSPA1A','HSPA2','HSPA4','HSPA5','HTRA1','HUWE1','IBSP','ICAM1','ICK','ID1','ID3','ID4','IDH1','IDH2','IFNG','IFRD1','IGF1','IGF1R','IGFBP2','IGFBP3','IKBKB','IKBKE','IKZF1','IKZF2','IKZF3','IL12RB2','IL13','IL13RA2','IL15','IL1B','IL1R1','IL1RAP','IL2','IL21R','IL2RA','IL3','IL6','IL7R','INHBA','INPP4A','INPP4B','INPP5A','INPP5D','IQCG','IRF1','IRF2BP2','IRF4',
					'IRF8','IRS1','IRS2','IRS4','ITGA5','ITGA7','ITGA8','ITGAV','ITGB3','ITK','ITPKA','JAG2','JAK1','JAK2','JAK3','JARID2','JAZF1','JUN','KALRN','KANK1','KAT2B','KAT6A','KAT6B','KCNB1','KDM1A','KDM2B','KDM4C','KDM5A','KDM5C','KDM6A','KDR','KDSR','KEAP1','KIAA0232','KIAA1524','KIAA1549','KIAA1598','KIF5B','KIT','KLF4','KLHL6','KLK2','KLK7','KMT2A','KMT2B','KMT2C','KMT2D','KPNB1','KRAS','KSR1','KTN1','LAMA1','LAMA5','LAMP2','LASP1','LCK','LCP1','LEF1','LEFTY2',
					'LFNG','LGALS3','LGR5','LHFP','LHX2','LHX4','LIFR','LINC00598','LINC00982','LINGO2','LMBRD1','LMO1','LMO2','LMO7','LNP1','LOX','LPAR1','LPP','LPXN','LRIG3','LRMP','LRP1B','LRP5','LRPPRC','LRRC37B','LRRC59','LRRC7','LRRK2','LTBP1','LYL1','LYN','MACROD1','MAD2L1','MADD','MAF','MAFB','MAGED1','MAGEE1','MALAT1','MALT1','MAML1','MAML2','MAP2','MAP2K1','MAP2K2','MAP2K3','MAP2K4','MAP2K5','MAP2K6','MAP2K7','MAP3K1','MAP3K14','MAP3K6','MAP3K7','MAPK1','MAPK3','MAPK8',
					'MAPK8IP2','MAPK9','MAPRE1','MATK','MAX','MB21D2','MBNL1','MBTD1','MCL1','MDC1','MDH1','MDM2','MDM4','MDS2','MEAF6','MECOM','MED12','MEF2B','MEF2C','MEF2D','MELK','MEN1','MET','METTL18','METTL7B','MFNG','MGEA5','MGMT','MIB1','MIPOL1','MITF','MKI67','MKL1','MKL2','MLF1','MLH1','MLLT1','MLLT10','MLLT11','MLLT3','MLLT4','MLLT6','MMP7','MMP9','MN1','MNAT1','MNX1','MPL','MRE11A','MSH2','MSH3','MSH6','MSI2','MSN','MTCP1','MTOR','MTUS2','MUC1','MUTYH','MYB','MYBL1',
					'MYC','MYCL','MYCN','MYD88','MYH11','MYH9','MYO18A','MYO1F','NAB2','NACA','NAPA','NAV3','NBEAP1','NBN','NBR1','NCAM1','NCKIPSD','NCOA1','NCOA2','NCOA3','NCOA4','NCOR2','NCSTN','NDC80','NDE1','NDRG1','NDUFAF1','NEDD4','NEURL1','NF1','NF2','NFATC1','NFATC2','NFE2L2','NFIB','NFKB1','NFKB2','NFKBIA','NGF','NGFR','NIN','NIPBL','NKX2-1','NKX2-5','NOD1','NODAL','NONO','NOS3','NOTCH1','NOTCH2','NOTCH3','NOTCH4','NPM1','NPM2','NR3C1','NR4A3','NR6A1','NRAS','NSD1',
					'NT5C2','NTF3','NTF4','NTRK1','NTRK2','NTRK3','NUMA1','NUP107','NUP214','NUP93','NUP98','NUTM1','NUTM2A','NUTM2B','OFD1','OLIG1','OLIG2','OLR1','OMD','P2RY8','PAFAH1B2','PAG1','PAK1','PAK3','PAK6','PAK7','PALB2','PAPPA','PASK','PATZ1','PAX3','PAX5','PAX7','PAX8','PBRM1','PBX1','PC','PCBP1','PCLO','PCM1','PCNA','PCSK7','PDCD1','PDCD11','PDCD1LG2','PDE4DIP','PDGFA','PDGFB','PDGFD','PDGFRA','PDGFRB','PDK1','PEG3','PER1','PFDN5','PHB','PHF1','PHF23','PHF6','PHOX2B',
					'PI4KA','PICALM','PIK3CA','PIK3CB','PIK3CD','PIK3CG','PIK3R1','PIK3R2','PIM1','PKM','PLA2G2A','PLA2G5','PLAG1','PLAT','PLAU','PLCB1','PLCB4','PLCG1','PLCG2','PLEKHM2','PML','PMS1','PMS2','POFUT1','POLD1','POLD4','POLR2H','POM121','POMGNT1','POSTN','POT1','POU2AF1','POU5F1','PPAP2B','PPARG','PPARGC1A','PPFIA2','PPFIBP1','PPM1D','PPP1CB','PPP1R13B','PPP1R13L','PPP2CB','PPP2R1A','PPP2R1B','PPP2R2B','PPP2R4','PPP3CA','PPP3CB','PPP3CC','PPP3R1','PPP3R2','PPP4C',
					'PQLC3','PRCC','PRDM1','PRDM16','PRDM7','PRF1','PRG2','PRICKLE1','PRKACA','PRKACG','PRKAR1A','PRKCA','PRKCB','PRKCD','PRKCG','PRKDC','PRKG2','PRMT1','PRMT8','PROM1','PRRX1','PRRX2','PRSS8','PSD3','PSEN1','PSIP1','PSMD2','PTBP1','PTCH1','PTCRA','PTEN','PTGS2','PTK2','PTK2B','PTK7','PTPN11','PTPN2','PTPN6','PTPRA','PTPRK','PTPRO','PTPRR','PTTG1','PVT1','RABEP1','RAC1','RAC2','RAC3','RAD21','RAD50','RAD51','RAD51B','RAD51C','RAD51D','RAD52','RAF1','RALGDS','RANBP17',
					'RANBP2','RAP1GDS1','RARA','RASAL1','RASGEF1A','RASGRF1','RASGRF2','RASGRP1','RB1','RBM15','RBM6','RCHY1','RCOR1','RCSD1','RECQL4','REEP3','RELA','RELN','RERG','RET','RGS7','RHBDF2','RHOA','RHOD','RHOH','RICTOR','RLTPR','RMI2','RNF213','RNF43','ROBO1','ROBO2','ROS1','RPA3','RPL22','RPN1','RPN2','RPS21','RPS6KA1','RPS6KA2','RPS6KA3','RPTOR','RREB1','RRM1','RRM2B','RTEL1','RTN3','RUNX1','RUNX1T1','RUNX2','RYR3','S1PR2','SARNP','SBDS','SCN8A','SDC4','SDHA','SDHAF2',
					'SDHB','SDHC','SDHD','SEC31A','SEPT2','SEPT5','SEPT6','SEPT9','SERP2','SERPINE1','SERPINF1','SET','SETBP1','SETD2','SETD7','SF3B1','SFPQ','SFRP2','SFRP4','SGK1','SGPP2','SH2D5','SH3BP1','SH3D19','SH3GL1','SH3GL2','SHC1','SHC2','SIK3','SIN3A','SIRT1','SKP2','SLC1A2','SLC34A2','SLC45A3','SLC7A5','SLCO1B3','SLX4','SMAD2','SMAD3','SMAD4','SMAD6','SMAP1','SMARCA1','SMARCA4','SMARCA5','SMARCB1','SMC1A','SMC3','SMO','SNAPC3','SNCG','SNHG5','SNW1','SNX29','SNX9',
					'SOCS1','SOCS2','SOCS3','SOD2','SORBS2','SORT1','SOS1','SOX10','SOX11','SOX2','SP1','SP3','SPECC1','SPEN','SPOP','SPP1','SPRY2','SPRY4','SPTAN1','SPTBN1','SQSTM1','SRC','SRF','SRGAP3','SRRM3','SRSF2','SRSF3','SS18','SS18L1','SSBP2','SSX1','SSX2','SSX4','ST6GAL1','STAG2','STAT1','STAT3','STAT4','STAT5A','STAT5B','STAT6','STIL','STK11','STL','STRN','STX5','STYK1','SUFU','SUGP2','SULF1','SUV39H2','SUZ12','SYK','SYP','TACC1','TACC2','TACC3','TAF1','TAF15','TAL1',
					'TAL2','TAOK1','TBL1XR1','TBX15','TCEA1','TCF12','TCF3','TCF7L2','TCL1A','TCL6','TCTA','TEAD1','TEAD2','TEAD3','TEAD4','TEC','TENM1','TERF1','TERF2','TERT','TET1','TET2','TFAP2A','TFDP1','TFE3','TFEB','TFG','TFPT','TFRC','TGFB2','TGFB3','TGFBI','TGFBR2','TGFBR3','THADA','THBS1','THRAP3','TIAM1','TIRAP','TLL2','TLR4','TLX1','TLX3','TMEM127','TMEM230','TMEM30A','TMPRSS2','TNC','TNF','TNFAIP3','TNFRSF10B','TNFRSF10D','TNFRSF11A','TNFRSF14','TNFRSF17','TNFRSF6B',
					'TOP1','TOP2A','TOP2B','TP53','TP53BP1','TP63','TP73','TPD52L2','TPM3','TPM4','TPO','TPR','TRAF2','TRAF3','TRAF5','TRHDE','TRIM24','TRIM27','TRIM33','TRIP11','TRPS1','TSC1','TSC2','TSHR','TTK','TTL','TUSC3','TYK2','TYMS','U2AF1','U2AF2','UBE2B','UBE2C','UFC1','UFM1','USP16','USP42','USP5','USP6','USP7','VCAM1','VEGFA','VEGFC','VGLL3','VHL','VTI1A','WASF2','WDFY3','WDR1','WDR18','WDR70','WDR90','WEE1','WHSC1','WHSC1L1','WIF1','WISP3','WNT10A','WNT10B','WNT11',
					'WNT16','WNT2B','WNT3','WNT4','WNT5B','WNT6','WNT7B','WNT8B','WRN','WSB1','WT1','WWOX','WWTR1','XBP1','XIAP','XKR3','XPA','XPC','XPO1','XRCC6','YAP1','YPEL5','YTHDF2','YWHAE','YY1AP1','ZBTB16','ZC3H7A','ZC3H7B','ZFP64','ZFPM2','ZFYVE19','ZIC2','ZMIZ1','ZMYM2','ZMYM3','ZMYND11','ZNF207','ZNF217','ZNF24','ZNF331','ZNF384','ZNF444','ZNF521','ZNF585B','ZNF687','ZNF703','ZRSR2'],
	'trusighttumor170' => ['ABL1','AKT1','AKT2','AKT3','ALK','APC','AR','ARID1A','ATM','ATR','AXL','BAP1','BARD1','BCL2','BCL6','BRAF','BRCA1','BRCA2','BRIP1','BTK','CARD11','CCND1','CCND2','CCND3','CCNE1','CD79A','CD79B','CDH1','CDK12','CDK4','CDK6','CDKN2A','CEBPA','CHEK1','CHEK2','CREBBP','CSF1R','CTNNB1','DDR2','DNMT3A','EGFR','EML4','EP300','ERBB2','ERBB3','ERBB4','ERCC1','ERCC2','ERG','ESR1','ETS1','ETV1','ETV4','ETV5','EWSR1','EZH2','FAM175A','FANCI','FANCL','FBXW7','FGF1','FGF10','FGF14',
				'FGF19','FGF2','FGF23','FGF3','FGF4','FGF5','FGF6','FGF7','FGF8','FGF9','FGFR1','FGFR2','FGFR3','FGFR4','FLI1','FLT1','FLT3','FOXL2','GEN1','GNA11','GNAQ','GNAS','HNF1A','HRAS','IDH1','IDH2','INPP4B','JAK2','JAK3','KDR','KIF5B','KIT','KMT2A (MLL)','KRAS','LAMP1','MAP2K1','MAP2K2','MCL1','MDM2','MDM4','MET','MLH1','MLLT3','MPL','MRE11A','MSH2','MSH3','MSH6','MTOR','MUTYH','MYC','MYCL1','MYCN','MYD88','NBN','NF1','NOTCH1','NOTCH2','NOTCH3','NPM1','NRAS','NRG1','NTRK1',
				'NTRK2','NTRK3','PALB2','PAX3','PAX7','PDGFRA','PDGFRB','PIK3CA','PIK3CB','PIK3CD','PIK3CG','PIK3R1','PMS2','PPARG','PPP2R2A','PTCH1','PTEN','PTPN11','RAD51','RAD51B','RAD51C','RAD51D','RAD54L','RAF1','RB1','RET','RICTOR','ROS1','RPS6KB1','SLX4','SMAD4','SMARCB1','SMO','SRC','STK11','TERT ','TET2','TFRC','TMPRSS2','TP53','TSC1','TSC2','VHL','XRCC2'],
	'ionampliseqcomprehensivecancerpanel' => ['ABL1','ABL2','ACVR2A','ADAMTS20','AFF1','AFF3','AKAP9','AKT1','AKT2','AKT3','ALK','APC','AR','ARID1A','ARID2','ARNT','ASXL1','ATF1','ATM','ATR','ATRX','AURKA','AURKB','AURKC','AXL','BAI3','BAP1','BCL10','BCL11A','BCL11B','BCL2','BCL2L1','BCL2L2','BCL3','BCL6','BCL9','BCR','BIRC2','BIRC3','BIRC5','BLM','BLNK','BMPR1A','BRAF','BRD3','BRIP1','BTK','BUB1B','CARD11','CASC5','CBL','CCND1','CCND2','CCNE1','CD79A','CD79B','CDC73','CDH1','CDH11','CDH2',
							'CDH20','CDH5','CDK12','CDK4','CDK6','CDK8','CDKN2A','CDKN2B','CDKN2C','CEBPA','CHEK1','CHEK2','CIC','CKS1B','CMPK1','COL1A1','CRBN','CREB1','CREBBP','CRKL','CRTC1','CSF1R','CSMD3','CTNNA1','CTNNB1','CYLD','CYP2C19','CYP2D6','DAXX','DCC','DDB2','DDIT3','DDR2','DEK','DICER1','DNMT3A','DPYD','DST','EGFR','EML4','EP300','EP400','EPHA3','EPHA7','EPHB1','EPHB4','EPHB6','ERBB2','ERBB3','ERBB4','ERCC1','ERCC2','ERCC3','ERCC4','ERCC5','ERG','ESR1',
							'ETS1','ETV1','ETV4','EXT1','EXT2','EZH2','FAM123B','FANCA','FANCC','FANCD2','FANCF','FANCG','FAS','FBXW7','FGFR1','FGFR2','FGFR3','FGFR4','FH','FLCN','FLI1','FLT1','FLT3','FLT4','FN1','FOXL2','FOXO1','FOXO3','FOXP1','FOXP4','FZR1','G6PD','GATA1','GATA2','GATA3','GDNF','GNA11','GNAQ','GNAS','GPR124','GRM8','GUCY1A2','HCAR1','HIF1A','HLF','HNF1A','HOOK3','HRAS','HSP90AA1','HSP90AB1','ICK','IDH1','IDH2','IGF1R','IGF2','IGF2R','IKBKB','IKBKE',
							'IKZF1','IL2','IL21R','IL6ST','IL7R','ING4','IRF4','IRS2','ITGA10','ITGA9','ITGB2','ITGB3','JAK1','JAK2','JAK3','JUN','KAT6A','KAT6B','KDM5C','KDM6A','KDR','KEAP1','KIT','KLF6','KRAS','LAMP1','LCK','LIFR','LPHN3','LPP','LRP1B','LTF','LTK','MAF','MAFB','MAGEA1','MAGI1','MALT1','MAML2','MAP2K1','MAP2K2','MAP2K4','MAP3K7','MAPK1','MAPK8','MARK1','MARK4','MBD1','MCL1','MDM2','MDM4','MEN1','MET','MITF','MLH1','MLL','MLL2','MLL3','MLLT10','MMP2',
							'MN1','MPL','MRE11A','MSH2','MSH6','MTOR','MTR','MTRR','MUC1','MUTYH','MYB','MYC','MYCL1','MYCN','MYD88','MYH11','MYH9','NBN','NCOA1','NCOA2','NCOA4','NF1','NF2','NFE2L2','NFKB1','NFKB2','NIN','NKX2-1','NLRP1','NOTCH1','NOTCH2','NOTCH4','NPM1','NRAS','NSD1','NTRK1','NTRK3','NUMA1','NUP214','NUP98','PAK3','PALB2','PARP1','PAX3','PAX5','PAX7','PAX8','PBRM1','PBX1','PDE4DIP','PDGFB','PDGFRA','PDGFRB','PER1','PGAP3','PHOX2B','PIK3C2B','PIK3CA',
							'PIK3CB','PIK3CD','PIK3CG','PIK3R1','PIK3R2','PIM1','PKHD1','PLAG1','PLCG1','PLEKHG5','PML','PMS1','PMS2','POT1','POU5F1','PPARG','PPP2R1A','PRDM1','PRKAR1A','PRKDC','PSIP1','PTCH1','PTEN','PTGS2','PTPN11','PTPRD','PTPRT','RAD50','RAF1','RALGDS','RARA','RB1','RECQL4','REL','RET','RHOH','RNASEL','RNF2','RNF213','ROS1','RPS6KA2','RRM1','RUNX1','RUNX1T1','SAMD9','SBDS','SDHA','SDHB','SDHC','SDHD','SEPT9','SETD2','SF3B1','SGK1','SH2D1A','SMAD2',
							'SMAD4','SMARCA4','SMARCB1','SMO','SMUG1','SOCS1','SOX11','SOX2','SRC','SSX1','STK11','STK36','SUFU','SYK','SYNE1','TAF1','TAF1L','TAL1','TBX22','TCF12','TCF3','TCF7L1','TCF7L2','TCL1A','TET1','TET2','TFE3','TGFBR2','TGM7','THBS1','TIMP3','TLR4','TLX1','TNFAIP3','TNFRSF14','TNK2','TOP1','TP53','TPR','TRIM24','TRIM33','TRIP11','TRRAP','TSC1','TSC2','TSHR','UBR5','UGT1A1','USP9X','VHL','WAS','WHSC1','WRN','WT1','XPA','XPC','XPO1','XRCC2','ZNF384',
							'ZNF521'],
	'ionampliseqcancerhotspotpanelv2' => ['ABL1','AKT1','ALK','APC','ATM','BRAF','CDH1','CDKN2A','CSF1R','CTNNB1','EGFR','ERBB2','ERBB4','EZH2','FBXW7','FGFR1','FGFR2','FGFR3','FLT3','GNA11','GNAQ','GNAS','HNF1A','HRAS','IDH1','IDH2','JAK2','JAK3','KDR','KIT','KRAS','MET','MLH1','MPL','NOTCH1','NPM1','NRAS','PDGFRA','PIK3CA','PTEN','PTPN11','RB1','RET','SMAD4','SMARCB1','SMO','SRC','STK11','TP53','VHL'],
	'carismolecularintelligence' => ['ABI1','ABL1','ABL2','ACKR3','ACSL3','ACSL6','AFF1','AFF3','AFF4','AKAP9','AKT1','AKT2','AKT3','ALDH2','ALK','AMER1','APC','AR','ARAF','ARFRP1','ARHGAP26','ARHGEF12','ARID1A','ARID2','ARNT','ASPSCR1','ASXL1','ATF1','ATIC','ATM','ATP1A1','ATP2B3','ATR','ATRX','AURKA','AURKB','AXIN1','AXL','BAP1','BARD1','BCL10','BCL11A','BCL11B','BCL2','BCL2L11','BCL2L2','BCL3','BCL6','BCL7A','BCL9','BCOR','BCORL1','BCR','BIRC3','BLM','BMPR1A','BRAF','BRCA1','BRCA2','BRD3',
						'BRD4','BRIP1','BTG1','BTK','BUB1B','C11orf30','C15orf65','C17orf39','C2orf44','CACNA1D','CALR','CAMTA1','CANT1','CARD11','CARS','CASC5','CASP8','CBFA2T3','CBFB','CBL','CBLB','CBLC','CCDC6','CCNB1IP1','CCND1','CCND2','CCND3','CCNE1','CD274','CD74','CD79A','CD79B','CDC73','CDH1','CDH11','CDK12','CDK4','CDK6','CDK8','CDKN1B','CDKN2A','CDKN2B','CDKN2C','CDX2','CEBPA','CHCHD7','CHEK1','CHEK2','CHIC2','CHN1','CIC','CIITA','CLP1','CLTC','CLTCL1','cMET',
						'CNBP','CNOT3','CNTRL','COL1A1','COPB1','COX6C','CREB1','CREB3L1','CREB3L2','CREBBP','CRKL','CRLF2','CRTC1','CRTC3','CSF1R','CSF3R','CTCF','CTLA4','CTNNA1','CTNNB1','CYLD','CYP2D6','DAXX','DDB2','DDIT3','DDR2','DDX10','DDX5','DDX6','DEK','DICER1','DNM2','DNMT3A','DOT1L','EBF1','ECT2L','EGFR','EIF4A2','ELF4','ELK4','ELL','ELN','EML4','EMSY','EP300','EPHA3','EPHA5','EPHB1','EPS15','ERBB2','ERBB3','ERBB4','ERC1','ERCC1','ERCC2','ERCC3','ERCC4','ERCC5',
						'ERG','ESR1','ETV1','ETV4','ETV5','ETV6','EWSR1','EXT1','EXT2','EZH2','EZR','FAM123B','FAM46C','FANCA','FANCC','FANCD2','FANCE','FANCF','FANCG','FANCL','FAS','FBXO11','FBXW7','FCRL4','FEV','FGF10','FGF14','FGF19','FGF23','FGF3','FGF4','FGF6','FGFR1','FGFR1OP','FGFR2','FGFR3','FGFR4','FH','FHIT','FIP1L1','FLCN','FLI1','FLT1','FLT3','FLT4','FNBP1','FOXA1','FOXL2','FOXO1','FOXO3','FOXO4','FOXP1','FSTL3','FUBP1','FUS','GAS7','GATA1','GATA2','GATA3','GID4',
						'GMPS','GNA11','GNA13','GNAQ','GNAS','GOLGA5','GOPC','GPC3','GPHN','GPR124','GRIN2A','GSK3B','H3F3A','H3F3B','HER2','HER3','HER4','HERPUD1','HEY1','HGF','HIP1','HIST1H3B','HIST1H4I','HLF','HMGA1','HMGA2','HMGN2P46','HNF1A','HNRNPA2B1','HOOK3','HOXA11','HOXA13','HOXA9','HOXC11','HOXC13','HOXD11','HOXD13','HRAS','HSP90AA1','HSP90AB1','IDH1','IDH2','IGF1R','IKBKE','IKZF1','IL2','IL21R','IL6ST','IL7R','INHBA','IRF4','IRS2','ITK','JAK1','JAK2','JAK3','JAZF1',
						'JUN','KAT6A','KAT6B','KCNJ5','KDM5A','KDM5C','KDM6A','KDR','KDSR','KEAP1','KIAA1549','KIF5B','KIT','KLF4','KLHL6','KLK2','KMT2A','KMT2C','KMT2D','KRAS','KTN1','LASP1','LCK','LCP1','LGR5','LHFP','LIFR','LK','LMO1','LMO2','LPP','LRIG3','LRP1B','LYL1','MAF','MAFB','MALT1','MAML2','MAP2K1','MAP2K2','MAP2K4','MAP3K1','MAX','MCL1','MDM2','MDM4','MDS2','MECOM','MED12','MEF2B','MEN1','MET','MITF','MKL1','MLF1','MLH1','MLL','MLL2','MLL3','MLLT1','MLLT10',
						'MLLT11','MLLT3','MLLT4','MLLT6','MN1','MNX1','MPL','MRE11A','MSH2','MSH6','MSI2','MSN','MTCP1','MTOR','MUC1','MUTYH','MYB','MYC','MYCL','MYCL1','MYCN','MYD88','MYH11','MYH9','MYST3','NACA','NBN','NCKIPSD','NCOA1','NCOA2','NCOA4','NDRG1','NF1','NF2','NFE2L2','NFIB','NFKB2','NFKBIA','NIN','NKX2-1','NONO','NOTCH1','NOTCH2','NPM1','NR4A3','NRAS','NSD1','NT5C2','NTRK1','NTRK2','NTRK3','NUMA1','NUP214','NUP93','NUP98','NUTM1','NUTM2B','OLIG2','OMD','P2RY8',
						'PAFAH1B2','PAK3','PALB2','PATZ1','PAX3','PAX5','PAX7','PAX8','PBRM1','PBX1','PCM1','PCSK7','PD1','PDCD1','PDCD1LG2','PDE4DIP','PDGFB','PDGFRA','PDGFRB','PDK1','PDL1','PDL2','PER1','PHF6','PHOX2B','PICALM','PIK3CA','PIK3CG','PIK3R1','PIK3R2','PIM1','PLAG1','PML','PMS1','PMS2','POLE','POT1','POU2AF1','POU5F1','PPARG','PPP2R1A','PRCC','PRDM1','PRDM16','PRF1','PRKAR1A','PRKDC','PRRX1','PSIP1','PTCH1','PTEN','PTPN11','PTPRC','RABEP1','RAC1','RAD21','RAD50',
						'RAD51','RAD51B','RAF1','RALGDS','RANBP17','RAP1GDS1','RARA','RB1','RBM15','RECQL4','REL','RET','RHOH','RICTOR','RMI2','RNF213','RNF43','ROS1','RPL10','RPL22','RPL5','RPN1','RPTOR','RSPO3','RUNX1','RUNX1T1','SBDS','SDC4','SDHAF2','SDHB','SDHC','SDHD','SEPT5','SEPT6','SEPT9','SET','SETBP1','SETD2','SF3B1','SFPQ','SH2B3','SH3GL1','SLC34A2','SLC45A3','SMAD2','SMAD4','SMARCA4','SMARCB1','SMARCE1','SMO','SNX29','SOCS1','SOX10','SOX2','SPECC1','SPEN','SPOP',
						'SRC','SRGAP3','SRSF2','SRSF3','SS18','SS18L1','SSX1','STAG2','STAT3','STAT4','STAT5B','STIL','STK11','SUFU','SUZ12','SYK','TAF15','TAL1','TAL2','TBL1XR1','TCEA1','TCF12','TCF3','TCF7L2','TCL1A','TERT','TET1','TET2','TFE3','TFEB','TFG','TFPT','TFRC','TGFBR2','THRAP3','TLX1','TLX3','TMPRSS2','TNFAIP3','TNFRSF14','TNFRSF17','TOP1','TP53','TPM3','TPM4','TPR','TRAF7','TRIM26','TRIM27','TRIM33','TRIP11','TRRAP','TSC1','TSC2','TSHR','TTL','U2AF1','UBR5',
						'USP6','VEGFA','VEGFB','VEGFR2','VHL','VTI1A','WAS','WHSC1','WHSC1L1','WIF1','WISP3','WRN','WT1','WWTR1','XPA','XPC','XPO1','YWHAE','ZBTB16','ZMYM2','ZNF217','ZNF331','ZNF384','ZNF521','ZNF703','ZRSR2'],
	'oncodeepandtrace' => ['ABL1','ACVRL1','AKT1','AKT3','ALK','APC','APEX1','AR','ARAF','ASXL1','ATM','ATP11B','AURKA','AXL','BAP1','BCL2L1','BCL9','BIRC2','BIRC3','BRAF','BRCA1','BRCA2','BRD3','BTK','CBL','CCND1','CCNE1','CD274','CD44','CDH1','CDK4','CDK6','CDKN2A','CHEK2','CSF1R','CSNK2A1','CTNNB1','CYP2C19','CYP2D6','DCUN1D1','DDR2','DNMT3A','DPYD','EGFR','EIF3E','EPHA3','EPHA5','EPHB1','ERBB2','ERBB3','ERBB4','ERG','ERRFI1','ESR1','ETV1','ETV4','ETV5','ETV6','EWSR1','EZH2','FANCA','FANCC',
				'FANCD2','FANCE','FANCF','FANCL','FAS','FBXW7','FGFR1','FGFR2','FGFR3','FGFR4','FLCN','FLT1','FLT3','FLT4','FOXL2','GAS6','GATA1','GATA2','GATA3','GATA6','GNA11','GNA13','GNAQ','GNAS','HGF','HNF1A','HRAS','IDH1','IDH2','IFITM1','IFITM3','IGF1R','IKBKE','IL6','INHBA','IRF2','JAK1','JAK2','JAK3','KDR','KEAP1','KEL','KIT','KNSTRN','KRAS','LYN','MAGOH','MAML2','MAP2K1','MAP2K2','MAPK1','MAX','MCL1','MDM2','MDM4','MED12','MET','MLH1','MPL','MRE11A','MSH2','MTOR','MYB',
				'MYC','MYCL','MYCN','MYD88','MYO18A','NCOA2','NF1','NF2','NFE2L2','NKX2-1','NKX2-8','NOTCH1','NPM1','NRAS','NTRK1','NTRK3','PALB2','PAX5','PD-1','PD-L1','PDCD1LG2','PDGFRA','PDGFRB','PIK3CA','PIK3CB','PIK3CG','PIK3R1','PIK3R2','PLAG1','PNP','POLD1','POLE','PPARG','PPP2R1A','PRKCI','PRKDC','PTCH1','PTEN','PTPN11','PTPRD','RAC1','RAD51','RAD51C','RAF1','RARA','RB1','RET','RHEB','RHOA','RNF43','ROS1','RPS6KB1','RPTOR','RUNX1','RUNX1T1','SETD2','SF3B1','SMAD4','SMARCA4',
				'SMARCB1','SMO','SOX2','SPOP','SRC','STAT3','STK11','TERT','TET2','TFE3','TGFBR2','TIAF1','TOP1','TOP2A','TP53','TPMT','TSC1','TSC2','TSHR','U2AF1','UGT1A1','VHL','WT1','XPO1','ZNF217'],
	'oncodeepclinical' => ['ABL1','ABL2','ACVR2A','ADAMTS20','AFF1','AFF3','AKAP9','AKT1','AKT2','AKT3','ALK','APC','AR','ARID1A','ARID2','ARNT','ASXL1','ATF1','ATM','ATR','ATRX','AURKA','AURKB','AURKC','AXL','BAI3','BAP1','BCL10','BCL11A','BCL11B','BCL2','BCL2L1','BCL2L2','BCL3','BCL6','BCL9','BCR','BIRC2','BIRC3','BIRC5','BLM','BLNK','BMPR1A','BRAF','BRD3','BRIP1','BTK','BUB1B','CARD11','CASC5','CBL','CCND1','CCND2','CCNE1','CD79A','CD79B','CDC73','CDH1','CDH11','CDH2','CMPK1','CDH20','CDH5',
				'CDK12','CDK4','CDK6','CDK8','CDKN2A','CDKN2B','CDKN2C','CEBPA','CHEK1','CHEK2','CIC','CKS1B','COL1A1','CRBN','CREB1','CREBBP','CRKL','CRTC1','CSF1R','CSMD3','CTNNA1','CTNNB1','CYLD','CYP2C19','CYP2D6','DAXX','DCC','DDB2','DDIT3','DDR2','DEK','DICER1','DNMT3A','DPYD','DST','EGFR','EML4','EP300','EP400','EPHA3','EPHA7','EPHB1','EPHB4','EPHB6','ERBB2','ERBB3','ERBB4','ERCC1','ERCC2','ERCC3','ERCC4','ERCC5','ERG','ESR1','ETS1','ETV1','ETV4','EXT1','EXT2','EZH2','FAM123B',
				'FANCA','FANCC','FANCD2','FANCF','FANCG','FAS','FBXW7','FGFR1','FGFR2','FGFR3','FGFR4','FH','FLCN','FLI1','FLT1','FLT3','FLT4','FN1','FOXL2','FOXO1','FOXO3','FOXP1','FOXP4','FZR1','G6PD','GATA1','HRAS','GATA2','GATA3','GDNF','GNA11','GNAQ','GNAS','GPR124','GRM8','GUCY1A2','HCAR1','HIF1A','HLF','HNF1A','HOOK3','HSP90AA1','HSP90AB1','ICK','IDH1','IDH2','IGF1R','IGF2','IGF2R','IKBKB','IKBKE','IKZF1','IL2','IL21R','IL6ST','IL7R','ING4','IRF4','IRS2','ITGA10','ITGA9',
				'ITGB2','ITGB3','JAK1','JAK2','JAK3','JUN','KAT6A','KAT6B','KDM5C','KDM6A','KDR','KEAP1','KIT','KLF6','KMT2A','KMT2C','KMT2D','KRAS','LAMP1','LCK','LIFR','LPHN3','LPP','LRP1B','LTF','LTK','MAF','MAFB','MAGEA1','MAGI1','MALT1','MAML2','MAP2K1','MAP2K2','MAP2K4','MAP3K7','MAPK1','MAPK8','MARK1','MARK4','MBD1','MCL1','MDM2','MDM4','MEN1','MET','MITF','MLH1','MLLT10','MMP2','MN1','MPL','MRE11A','MSH2','MSH6','NCOA2','MTOR','MTR','MTRR','MUC1','MUTYH','MYB','MYC','MYCL1',
				'MYCN','MYD88','MYH11','MYH9','NBN','NCOA1','NCOA4','NF1','NF2','NFE2L2','NFKB1','NFKB2','NIN','NKX2-1','NLRP1','NOTCH1','NOTCH2','NOTCH4','NPM1','NRAS','NSD1','NTRK1','NTRK3','NUMA1','NUP214','NUP98','PAK3','PALB2','PARP1','PAX3','PAX5','PAX7','PAX8','PBRM1','PBX1','PDE4DIP','PDGFB','PDGFRA','PDGFRB','PER1','PGAP3','PHOX2B','PIK3C2B','PIK3CA','PIK3CB','PIK3CD','PIK3CG','PIK3R1','PIK3R2','PIM1','PKHD1','PLAG1','PLCG1','PLEKHG5','PML','PMS1','PMS2','POT1','POU5F1',
				'PPARG','PPP2R1A','PRDM1','PRKAR1A','PRKDC','PSIP1','PTCH1','PTEN','PTGS2','PTPN11','PTPRD','PTPRT','RAD50','RAF1','RALGDS','RARA','RB1','RECQL4','REL','RET','RHOH','RNASEL','RNF2','RNF213','ROS1','RPS6KA2','RRM1','RUNX1','RUNX1T1','SAMD9','SBDS','SDHA','SDHB','SDHC','SDHD','SEPT9','SETD2','SF3B1','SGK1','SH2D1A','SMAD2','SMAD4','SMARCA4','SMARCB1','SMO','SMUG1','SOCS1','SOX11','SOX2','SRC','SSX1','STK11','STK36','SUFU','SYK','SYNE1','TAF1','TAF1L','TAL1','TBX22',
				'TCF12','TCF3','TCF7L1','TCF7L2','TCL1A','TET1','TET2','TFE3','TGFBR2','TGM7','THBS1','TIMP3','TLR4','TLX1','TNFAIP3','TNFRSF14','TNK2','TOP1','TP53','TPR','TRIM24','TRIM33','TRIP11','TSHR','UBR5','UGT1A1','USP9X','VHL','WAS','WHSC1','WRN','WT1','X PA','XPC','XPO1','TRRAP','TSC1','TSC2','XRCC2','ZNF384','ZNF521'],
	'oncodeepdx' => ['ABL1','AKT1','ALK','APC','ATM','AURKA','BRAF','CCDN1','CCNE1','CDH1','CDK4','CDKN2A','cMYC','CSF1R','CTNNB1','DDR2','DPYD','EGFR','ERBB2','ERBB4','ESR1','EZH2','FBXW7','FGFR1','FGFR2','FGFR3','FLT3','GNA11','GNAQ','GNAS','HNF1A','HRAS','IDH1','IDH2','JAK2','JAK3','KDR','KIT','KRAS','MAP2K1','MCL1','MDM2','MET','MLH1','MPL','NF1','NOTCH1','NPM1','NRAS','PDGFRA','PDGFRB','PIK3CA','PTEN','PTPN11','RB1','RET','SMAD4','SMARCB1','SMO','SRC','STK11','TP53','TPMT','UGT1A','VHL'],
	'guardant360' => ['AKT1','ALK','APC','AR','ARAF','ARID1A','ATM','BRAF','BRCA1','BRCA2','CCND1','CCND2','CCNE1','CDH1','CDK4','CDK6','CDKN2A','CTNNB1','DDR2','EGFR','ERBB2','ESR1','EZH2','FBXW7','FGFR1','FGFR2','FGFR3','GATA3','GNA11','GNAQ','GNAS','HNF1A','HRAS','IDH1','IDH2','JAK2','JAK3','KIT','KRAS','MAP2K1','MAP2K2','MAPK1','MAPK3','MET','MLH1','MPL','MTOR','MYC','NF1','NFE2L2','NOTCH1','NPM1','NRAS','NTRK1','NTRK3','PDGFRA','PIK3CA','PTEN','PTPN11','RAF1','RB1','RET','RHEB','RHOA',
				'RIT1','ROS1','SMAD4','SMO','STK11','TERT','TP53','TSC1','VHL'],
};

map push(@{$CANCER_GENES->{'any'}},@{$CANCER_GENES->{$_}}) , keys %$CANCER_GENES;
# Data types supported
my @DATA_TYPES = ('exome', 'rna', 'panel');

my (@INP_reads_files,%INP_params,$INP_params_file,$INP_data_type,$INP_outpath,$INP_shuffle,$INP_tech,$INP_nreads,$INP_nreads_amplicon,$INP_threads,$INP_zip,$INP_verbose,$INP_update,$INP_test);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s{,}' => \@INP_reads_files,
	'd|data=s' => \$INP_data_type,
	'o|output=s' => \$INP_outpath,
	'f|format=s' => \$INP_outformat,
	'p|params=s' => \$INP_params_file,
	'c|cancer=s' => \$INP_params{'cancer_type'},
	'al|align=s' => \$INP_params{'align_params'},
	'in|indels' => \$INP_params{'indels'},
	'ra|rare' => \$INP_params{'rare_variants'},
	'fr|freq=i' => \$INP_params{'min_variant_frequency'},
	'cv|cover=i' => \$INP_params{'min_variant_coverage'},
	'thr|threads=i' => \$INP_threads,
# 	'z|zip' => \$INP_zip,
# 	'v|verbose' => \$INP_verbose,
	'u|update' => \$INP_update,
	'test' => \$INP_test,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <type> [options]\n";
	print "\nOptions:\n";
	print "  -i <file1> <file2>\n\t\tInput one single-end or two paired-end FASTQ or FASTA files (compressed or uncompressed).\n";
	print "  -d <type>\tData type to analyze ('".join("', '",@DATA_TYPES)."').\n";
	print "  -o <path>\tOutput folder name.\n";
	print "  -f <format>\tOutput format (default='$INP_outformat', 'text long', 'html').\n";
	print "  -p <file>\tCSV file with params data.\n";
	print "  -c <type>\tCancer type (default='$DEFAULT_PARAMS{'cancer_type'}', '".join("', '", sort {$a cmp $b} keys %$CANCER_GENES)."').\n";
	print "  -al <type>\tRead alignment type (default='".$DEFAULT_PARAMS{'align_params'}."', 'bowtie2 --sensitive', 'bowtie2 --very-sensitive-local', 'bowtie2 --very-sensitive', 'bowtie2 --fast-local', 'bowtie2 --fast', ...)\n";
	print "  -in\t\tAnnotates insertions and deletions (default=".$DEFAULT_PARAMS{'indels'}.").\n";
	print "  -ra\t\tAnnotates rare variants (default=".$DEFAULT_PARAMS{'rare_variants'}.").\n";
	print "  -fr <number>\tMinimum frequency required to annotate variants (default=".$DEFAULT_PARAMS{'min_variant_frequency'}.").\n";
	print "  -cv <number>\tMinimum coverage required to annotate variants (default=".$DEFAULT_PARAMS{'min_variant_coverage'}.").\n";
# 	print "  -n <number>\tNumber of reads/sequences to analyze.\n";
# 	print "  -na <number>\tNumber of reads/sequences per amplicon to analyze.\n";
# 	print "  -al <perc>\tMinimum sequence percentage required to align to a reference to be annotated (default=".$DEFAULT_PARAMS{'min_allele_align'}."%).\n";
# 	print "  -id <perc>\tMinimum percentage of the aligned fragment required to be identical to a reference to be annotated (default=".$DEFAULT_PARAMS{'min_allele_ident'}."%).\n";
# 	print "  -fr <freq>\tFilter sequences with lower frequency after clustering.\n";
# 	print "  -ci <number>\tCluster together sequences with higher or equal identity.\n";
# 	print "  -min <number>\tAmplicons with lower total depth/coverage will be discarded (default=".$DEFAULT_PARAMS{'min_amplicon_depth'}.").\n";
# 	print "  -ks\t\tKeep singletons after clustering.\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
# 	print "  -z\t\tCompress results in ZIP format.\n";
# 	print "  -v\t\tPrint AmpliSAS output and verbose FASTA files.\n";
# 	print "  -u\t\tUpdates the human reference files to the latest versions.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}



if (defined($INP_update)){

	print "\nRunning '$COMMAND_LINE'\n";

	# Updates reference sequences files

	print "\nUpdating human genome reference files to the latest version from the Ensembl database.\n";
	if (!-d $DEFAULT_PARAMS{'reference_data_folder'}){
		mkdir($DEFAULT_PARAMS{'reference_data_folder'});
	} else {
		printf("\nWARNING: Destination folder '%s' already exists.\n",$DEFAULT_PARAMS{'reference_data_folder'});
		print "Do you want to overwrite its contents? (yes,no): ";
		my $option = <STDIN>;
		chomp $option;
		if ($option ne 'y' && $option ne 'yes') {
			print "\nExiting without changes\n\n";
			exit;
		}
	}
	# CORRECT BECAUSE URL IS GZIP (.fasta.gz) AND FILE IS NOT (.fa):
	my $reference_genome_sequences_file = download_file($DEFAULT_PARAMS_DNA{'reference_sequences_url'},$DEFAULT_PARAMS_DNA{'reference_sequences_file'});
	my $reference_cds_sequences_file = download_file($DEFAULT_PARAMS_RNA{'reference_sequences_url'},$DEFAULT_PARAMS_RNA{'reference_sequences_file'});
	my $reference_genome_annotations_file = download_file($DEFAULT_PARAMS_DNA{'reference_annotations_url'},$DEFAULT_PARAMS_DNA{'reference_annotations_file'});
	if (!is_fasta($reference_genome_sequences_file) || !is_fasta($reference_cds_sequences_file)) {
		print "\nERROR: Reference sequence files cannot be downloaded from Ensebl database server.\n\n";
		exit;
	}
	printf("\nUpdated '%s' Ensembl genome reference file.\n",$reference_genome_sequences_file);
	printf("\nUpdated '%s' Ensembl CDS reference file.\n",$reference_cds_sequences_file);
	printf("\nUpdated '%s' Ensembl annotations file.\n",$reference_genome_annotations_file);
	print "\n";

	# Updates COSMIC annotations file (coding mutations in genes listed in the Cancer Gene Census)
	
	print "\nUpdating coding mutations in genes listed in the Cancer Gene Census to the latest version from the COSMIC database.\n";
	if (!-d $DEFAULT_PARAMS{'mutation_data_folder'}){
		mkdir($DEFAULT_PARAMS{'mutation_data_folder'});
	} else {
		printf("\nWARNING: Destination folder '%s' already exists.\n",$DEFAULT_PARAMS{'reference_data_folder'});
		print "Do you want to overwrite its contents? (yes,no): ";
		my $option = <STDIN>;
		chomp $option;
		if ($option ne 'y' && $option ne 'yes') {
			print "\nExiting without changes\n\n";
			exit;
		}
	}
	use Net::SFTP::Foreign;
	my $host = 'sftp-cancer.sanger.ac.uk';
	my $port = 22;
	my $userid = 'sebalv@amu.edu.pl';
	my $pword = 'amplicancer';
	my $sftp = Net::SFTP::Foreign->new($host, timeout => 240,user => $userid, password => $pword) or die "\nERROR: Could NOT create SFTP connection\n\n";
	$sftp->get($DEFAULT_PARAMS{'mutation_annotations_url'},$DEFAULT_PARAMS{'mutation_annotations_file'}) or die "\nERROR: Mutation annotation file cannot be downloaded from COSMIC server.\n\n";
	$sftp->disconnect;
	printf("\nUpdated '%s' COSMIC mutation annotation file.\n",$DEFAULT_PARAMS{'mutation_annotations_file'});
	print "\n";

	exit;

}

# Downloads files if links are provided
my @INP_file_urls;
if ($INP_reads_files[0] =~ /(ftp|http|https)\:/ ){
	for (my $i=0; $i<=$#INP_reads_files; $i++) {
		$INP_file_urls[$i] = $INP_reads_files[$i];
		printf("\nDownloading file from '%s'\n",$INP_file_urls[$i]);
		$INP_reads_files[$i] = download_file($INP_file_urls[$i]);
		if (!defined($INP_reads_files[$i])){
			printf("\nERROR: broken link: '%s'.\n\n",$INP_file_urls[$i]);
			exit;
		}
	}
}

# Checks if an excel file or multifle are given as input with sequences previously processed
if (!@INP_reads_files){
	print "\nERROR: You must specify at least one input file.\n\n";
	usage();
	exit;
} elsif (!defined($INP_data_type) || !in_array(\@DATA_TYPES,"/$INP_data_type/")){
	print "\nERROR: Please specify a valid data type.\n\n";
	usage();
	exit;
}

# Default output path
if (!defined($INP_outpath)){
	$INP_outpath  = lc((split('\.',$SCRIPT_NAME))[0]);
}
# Default filename for result files
my $infile_name;
if ($INP_reads_files[0] =~ /(.+\/)?(.+)\.(fastq|fastq|fa|fq|fna)/ || $INP_reads_files[0] =~ /(.+\/)?(.+?)\./){
	$infile_name = $2;
}

# Comment program output if we desire HTML final data
if ($INP_outformat =~ /htm/i){
	$INP_outformat = 'html';
	print "<!--\n";
}

print "\nRunning '$COMMAND_LINE'\n";

# Creates output folder
if (!-d $INP_outpath){
	mkdir($INP_outpath);
}

# Sets default params for the type of data to analyze
if ($INP_data_type =~ /rna/i){
	%DEFAULT_PARAMS = ( %DEFAULT_PARAMS, %DEFAULT_PARAMS_RNA );
} else {
	%DEFAULT_PARAMS = ( %DEFAULT_PARAMS, %DEFAULT_PARAMS_DNA );
}


# Reads parameters data from CVS input file if defined
my %paramsdata;
if (defined($INP_params_file)){
	my $paramsdata_ = read_amplicon_data($INP_params_file,'params',['only params']);
	# Stores the params in a simple hash
	map $paramsdata{$_} = $paramsdata_->{$_}[0] , keys %$paramsdata_;
}
# Command line params have priority over CSV file (if there is no command line neither CSV param, then use default value)
foreach my $param (keys %DEFAULT_PARAMS) {
	# If specified cancer type is not supported
	if (defined($INP_params{$param})){
		$paramsdata{$param} = lc($INP_params{$param});
	} elsif (!defined($paramsdata{$param})){
		$paramsdata{$param} = $DEFAULT_PARAMS{$param};
	}
	if ($param eq 'cancer_type' && !defined($CANCER_GENES->{$paramsdata{'cancer_type'}})){
		print "\nWARNING: '".$paramsdata{'cancer_type'}."' cancer type not recognized or not implemented in the program.\n\n";
		$paramsdata{'cancer_type'} = 'any'
	}

}
write_to_file("$INP_outpath/amplicancer_params.csv", print_amplicon_data(\%paramsdata,'params',[keys %DEFAULT_PARAMS]));

# Sets genomic as default data type if nothing specified
my ($reference_sequences_type, $reference_sequences_file, $reference_annotations_file);
if ($INP_data_type =~ /rna/i){
	$reference_sequences_type = 'cds';
	$reference_sequences_file = $paramsdata{'reference_sequences_file'};
} else {
	$reference_sequences_type = 'genome';
	$reference_sequences_file = $paramsdata{'reference_sequences_file'};
	$reference_annotations_file = $paramsdata{'reference_annotations_file'};
}

# Checks if it exists reference sequences file
if (!defined($reference_sequences_file) || !-e $reference_sequences_file){
	print "\nERROR: '$reference_sequences_file' FASTA file with reference sequences cannot be found.\n\n";
	usage();
	exit;
}

my $variants;
my $ref_mapped_lengths;
if ( $paramsdata{'align_params'} =~ /bowtie2 (.+)/i ){

	my $align_params = $1;
	# No output messages:
	$align_params .= " --quiet";
	# Maximum 2 alignments per read:
	#$align_params .= " -k 2";

	if (defined($INP_threads) && $INP_threads>1){
		$align_params .= " --threads $INP_threads";
	}
	my $bamfile;
	if (is_bam($INP_reads_files[0])){
		$bamfile = $INP_reads_files[0];
	} elsif (defined($infile_name)){
		$bamfile = "$INP_outpath/$infile_name.$reference_sequences_type.bam";
	} else {
		$bamfile = "$INP_outpath/variants.$reference_sequences_type.bam";
	}
	unless (defined($INP_test) && is_bam($bamfile)){
		print "\nMapping variant sequences to references with Bowtie2.\n";
		execute_bowtie2(\@INP_reads_files,$reference_sequences_file,$align_params,$bamfile);
		if (!is_bam($bamfile)){
			print "\nERROR: It was an error in the read mapping and BAM file could not be generated.\n\n";
			exit;
		}
		# Indexes BAM file for variant retrieval from references and visualization
		execute_samtools($bamfile, 'index', "$bamfile.bai");
	}

	# Retrieves information about variant mutations from the alignment of reads/variants to references
	print ("\nRetrieving variant information.\n\n");
	($variants,$ref_mapped_lengths) = retrieve_variants_from_bam_data($bamfile,$reference_sequences_file,$reference_annotations_file,\%paramsdata);

	# Remove downloaded files:
	if (@INP_file_urls) {
		foreach my $file (@INP_reads_files) {
			# print "rm $file\n";
			`rm $file`;
		}
	}

} else {
	print "\nERROR: Alignment method not supported: '".$paramsdata{'align_params'}."'.\n\n";
	exit;
}

# if (!defined($variants->{'annotations'})) {
# 	print "\nERROR: It was an error retrieving variants.\n\n";
# 	exit;
# }
# printf("\n%d variants found.\n", scalar keys %{$variants->{'annotations'}});
# write_to_file('/tmp/kk',Dumper($variants));
# $variants = do '/tmp/kk';
# $ref_mapped_lengths->{'coding'} = 91630;

# Reads the genome positions of the variants/mutations and their frequencies
my %variant_genome_pos_to_freq;
my $count_total_variants = 0;
# For small panels it takes the coding and non-coding regions of variants to increase the mapped length (higher accuracy)
if ($INP_data_type =~ /panel/i && $ref_mapped_lengths->{'coding'} < 500000){
	foreach my $ref_id (keys %{$variants->{'filtered'}}){
		foreach my $ref_pos (keys %{$variants->{'filtered'}{$ref_id}}){
			foreach my $ref_mut (keys %{$variants->{'filtered'}{$ref_id}{$ref_pos}}){
				# TO IMPROVE: annotates only one mutation per position
				if (!defined($variant_genome_pos_to_freq{"$ref_id:$ref_pos"})){
					$variant_genome_pos_to_freq{"$ref_id:$ref_pos"} = sprintf("%.2f",eval($variants->{'filtered'}{$ref_id}{$ref_pos}{$ref_mut}));
					$count_total_variants++;
				}
			}
		}
	}
# Only coding regions:
} else {
	foreach my $ensembl_id (keys %{$variants->{'annotations'}}){
		foreach my $variant (@{$variants->{'annotations'}{$ensembl_id}}){
			# TO IMPROVE: annotates only one mutation per position
			if (!defined($variant_genome_pos_to_freq{$variant->{'genome_position'}})){
				$variant_genome_pos_to_freq{$variant->{'genome_position'}} = sprintf("%.2f",eval($variant->{'coverage'}));
				$count_total_variants++;
			}
		}
	}
}

my ($tmb, %tmb_messages, %mutation_classes);
# RNA-seq data does not allow to retrieve somatic/germline mutation information
if ($INP_data_type =~ /rna/i){

	printf("\tTotal variants: %d\n\n", $count_total_variants);
	while(my ($genome_pos,$allele_freq) = each %variant_genome_pos_to_freq){
		$mutation_classes{$genome_pos} .= 'unknown';
	}

# Discovers and annotaes somatic/germline mutations and COSMIC annotations
} else {

	print "\nDiscovering and filtering germline variants.\n\n";

	# Retrieves information from the COSMIC database of coding mutations in genes listed in the Cancer Gene Census
	my $cosmic_data = read_cosmic_data_file($paramsdata{'mutation_annotations_file'},[keys %{$variants->{'annotations'}}],['Gene name','Primary site','Primary histology','Mutation ID','Mutation CDS','Mutation AA','Mutation Description','Mutation zygosity','Mutation genome position','Pubmed_PMID']);
	# write_to_file('/tmp/kk2',Dumper($cosmic_data)); exit;
	# my $cosmic_data = do '/tmp/kk2';
	# Annotates COSMIC IDs and their genome positions from COSMIC annotations
	my %cosmic_mutation_ids;
	foreach my $ensembl_id (keys %$cosmic_data){
		foreach my $cosmic_mutation (@{$cosmic_data->{$ensembl_id}}){
			if ($cosmic_mutation->{'Mutation genome position'} =~ /(\d+\:\d+)\-\d+/){
				$cosmic_mutation_ids{$1} = $cosmic_mutation->{'Mutation ID'};
			}
		}
	}
	# Finds variant mutations in COSMIC database and their frequencies
	my (@cosmic_variants_allele_frequencies,@somatic_variants_allele_frequencies,@germline_variants_allele_frequencies,@unknown_variants_allele_frequencies);
	my ($count_cosmic_variants,$count_somatic_variants,$count_germline_variants,$count_unknown_variants) = (0,0,0,0);
	while(my ($genome_pos,$allele_freq) = each %variant_genome_pos_to_freq){
		if (defined($cosmic_mutation_ids{$genome_pos})){
			$mutation_classes{$genome_pos} = $cosmic_mutation_ids{$genome_pos}."-"; # COSMIC ID, ex. COSM6224
			if ($variant_genome_pos_to_freq{$genome_pos} < 0.4){
				push(@cosmic_variants_allele_frequencies,$variant_genome_pos_to_freq{$genome_pos});
			}
			$count_cosmic_variants++;
		}
	}
	# Establishes the allele frequency thresholds to consider the variant a somatic mutation as COSMIC ones
	# The median excludes extreme frequencies from COSMIC mutations that could be germline
	my $median_somatic_allele_frequency = median(@cosmic_variants_allele_frequencies);
	while(my ($genome_pos,$allele_freq) = each %variant_genome_pos_to_freq){
		if ($variant_genome_pos_to_freq{$genome_pos} > $median_somatic_allele_frequency*.8
		&& $variant_genome_pos_to_freq{$genome_pos} < $median_somatic_allele_frequency*1.2) {
			$mutation_classes{$genome_pos} .= 'somatic';
			push(@somatic_variants_allele_frequencies,$variant_genome_pos_to_freq{$genome_pos});
			$count_somatic_variants++;
		} elsif ($variant_genome_pos_to_freq{$genome_pos} > 0.4){
			$mutation_classes{$genome_pos} .= 'germline';
			push(@germline_variants_allele_frequencies,$variant_genome_pos_to_freq{$genome_pos});
			$count_germline_variants++;
		} else {
			$mutation_classes{$genome_pos} .= 'unknown';
			push(@unknown_variants_allele_frequencies,$variant_genome_pos_to_freq{$genome_pos});
			$count_unknown_variants++;
		}
	}

	printf("\tTotal variants: %d (%d in COSMIC database), freqs: %s\n", $count_total_variants, $count_cosmic_variants, join(', ',@cosmic_variants_allele_frequencies));
	printf("\tSomatic mutations: %d, freqs: %s\n", $count_somatic_variants,  join(', ',@somatic_variants_allele_frequencies));
	printf("\tGermline mutations: %d, freqs: %s\n", $count_germline_variants,  join(', ',@germline_variants_allele_frequencies));
	printf("\tOther mutations: %d, freqs: %s\n", $count_unknown_variants,  join(', ',@unknown_variants_allele_frequencies));

	# Calculates Tumor Mutation Burden (mapped length is multiplied by 2 becasuse both DNA strands are consider for mutations)
	# my $tmb = sprintf("%.0f", 1000000*($count_total_variants-$count_cosmic_variants-$count_germline_variants)/$ref_mapped_lengths->{'coding'});
	my $tmb_mapped_length;
	if ($INP_data_type =~ /panel/i && $ref_mapped_lengths->{'coding'} < 500000){
		$tmb_mapped_length = $ref_mapped_lengths->{'filtered'};
	} else {
		$tmb_mapped_length = $ref_mapped_lengths->{'coding'};
	}
	$tmb = sprintf("%.0f", 1000000*$count_somatic_variants/$tmb_mapped_length);
	if ($tmb) {
		printf("\nTumor Mutation Burden (*TMB): %d mutations per megabase\n\n", $tmb);
		if ($tmb_mapped_length < 500000) {
			$tmb_messages{'warning'} = sprintf("TMB calculation may be unaccurate because there is only %.2f megabase mapped (recommended: >0.5Mb).",$tmb_mapped_length/1000000);
			printf("\tWARNING: %s\n\n", $tmb_messages{'warning'});
		}
		if ($INP_data_type =~ /panel/i && $ref_mapped_lengths->{'coding'} < 500000){
			$tmb_messages{'calculation'} = sprintf("TMB is calculated by dividing the number of somatic mutations (%d, selected by an automatic algorithm) by the mapped length (%d, after filtering low freq. mutations, including coding and non-coding regions).", $count_somatic_variants, $tmb_mapped_length);
		} else {
			$tmb_messages{'calculation'} = sprintf("TMB is calculated by dividing the number of somatic mutations (%d, selected by an automatic algorithm) by the mapped length (%d, after filtering low freq. mutations, including only coding regions).", $count_somatic_variants, $tmb_mapped_length);
		}
		printf("\t%s\n", $tmb_messages{'calculation'});
	} else {
		$tmb_messages{'error'} = sprintf("\nThere was an error and TMB could not be calculated.\n\n");
		printf("\t%s\n", $tmb_messages{'error'});
	}
	# print "\t*TMB is calculated by dividing the number of somatic mutations (selected by an in-house algorithm and excluding COSMIC ones) by the total mapped length in coding regions of the genome.\n";
	# print "\t*TMB formula: (GENOMESUBS-COSMIC-GERMLINE)/GENOMEMAP
	# \t'GENOMESUBS' is the total number of genome positions in which references have nucleotide substitutions,
	# \t'COSMIC' is the number of these positions present in genes from the Cancer Gene Census list and annotated by the COSMIC database,
	# \t'GERMLINE' is the number of positions assigned to germline mutations by an in-house algorithm and
	# \t'GENOMEMAP' is the total mapped length in the coding regions of the genome.\n";
	# exit;
}


# Formats and prints variant mutation information
print "\nFormating and printing variant information.\n\n";
my $variants_output = print_variants($variants->{'annotations'},\%mutation_classes,\%paramsdata,[$INP_outformat]);

# Comment program output if we desire HTML final data
if ($INP_outformat eq 'html'){
	print "\n-->\n\n";
}
# Prints additional information about mapping, variants and tumor mutation burden
if ($INP_outformat eq 'html' && $INP_data_type !~ /rna/i){
# 	print "<b>\n";
# 	printf("Total variants: %d<br>\n", $count_total_variants);
# 	printf("In COSMIC database: %d<br>\n", $count_cosmic_variants);
# 	printf("Somatic mutations: %d<br>\n", $count_somatic_variants);
# 	printf("Germline mutations: %d<br>\n", $count_germline_variants);
# 	printf("Other mutations: %d<br>\n", $count_unknown_variants);
# 	print "</b><br>\n";
	if ($tmb) {
		printf("<font size=\"+2\"><b>Tumor Mutation Burden (TMB<sup>*</sup>): <font color=\"red\">%d mutations per megabase</font></b></font><br>\n",$tmb);
		printf("<b><sup>*</sup></b>%s<br>\n", $tmb_messages{'calculation'});
		if (defined($tmb_messages{'warning'})) {
			printf("<font color=\"red\"><b>%s</b></font><br>\n",$tmb_messages{'warning'});
		}
	} else {
		printf("<font color=\"red\"><b>%s</b></font><br>\n",$tmb_messages{'error'});
	}
	print "<br><br>";
}

# Prints variant mutation information
if (defined($variants_output)){
	if ($INP_outformat eq 'html'){
# 		print "<html>\n<head>\n<link rel=\"stylesheet\" type=\"text/css\" href=\"styles.css\">\n</head>\n<body>\n\n";
		print $variants_output;
# 		print "\n</body>\n</html>\n\n";
	} else {
		print $variants_output;
	}
} else {
	print "\nERROR: It was an error retrieving variants.\n\n";
	exit;
}



exit;


#########################################################################################

