package Bio::Onco;

######################################################################################################
#
# Bio::Onco - Package of subroutines for tumour genotyping
#
# Author: Alvaro Sebastian
#
# Support: Alvaro Sebastian (sixthresearcher@gmail.com)
#
# Sixth Researcher
# http://www.sixthresearcher.com
#
# Description:
# Package of subroutines for tumour genotyping
#
######################################################################################################

# # Returns full path and file name
# sub dirname {
# 	my $path = shift;
# 	if ($path =~ m{^(.*/)?.*}s){
# 		if (defined($1)) {
# 			return $1;
# 		} else {
# 			return '';
# 		}
# 	} else {
# 		print "\nERROR: The path of the script couldn't be found, run it in a Linux system.\n\n";
# 		exit;
# 	}
# }
# # Libraries are in folder '../' in the path of the script
# use lib dirname(__FILE__).'../';

# Modules are in folder '../' in the path of the script
use File::FindLib '../';
use Sort::Naturally;
use threads;
use threads::shared;
use Bio::Sequences;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;
use 5.004;
use strict;
use warnings;
no warnings ('uninitialized', 'substr');
use Time::HiRes qw(gettimeofday);
use File::Basename;
use HTTP::Tiny;
use Data::Dumper;

# Turn autoflush on
local $| = 1;

use Exporter;
 
our @ISA       = qw(Exporter);
our $VERSION   = '1.0';
our @EXPORT = qw(
read_oncodata_file
retrieve_variants
retrieve_variants_from_bam_data
annotate_reference_coding_variants
filter_variants
retrieve_transcript_variants
retrieve_alleles_from_bam_data
print_variants
retrieve_genbank_data
retrieve_uniprot_data
read_uniprot_single_file
read_genbank_single_file
read_uniprot_idmapping_file
read_cosmic_data_file

);

######################################################################################################

# Routes to binary files and databases
my $TOOLSDIR = dirname (__FILE__).'/tools/';
my $SAMTOOLS = $TOOLSDIR.'samtools';
# my $TOOLSDIR = dirname (__FILE__).'/tools/';
# my $NEEDLEALLEXE = $TOOLSDIR.'needleall';
# my $NEEDLEMANWUNSCHEXE = $TOOLSDIR.'needleman_wunsch';
# my $SMITHWATERMANEXE = $TOOLSDIR.'smith_waterman';
# my $MAFFTEXE = $TOOLSDIR.'mafft';

#################################################################################


# Reads information from .xlsx file with cancer related genes, alterations, drugs...:
sub read_oncodata_file {

	my ($file,$options) = @_;

	if (!-e $file){
		print "\nERROR 'read_oncodata_file': file '$file' doesn't exist.\n\n";
		exit;
	}
	
	my ($verbose) = (0);
	if (in_array($options, 'verbose')){
		$verbose = 1;
	}
	
	my $oncodata;
	
	my $converter = Text::Iconv -> new ("utf-8", "windows-1251");
	my $excel = Spreadsheet::XLSX -> new ($file, $converter);

	if (!defined($verbose)){
		$verbose = 0;
	}

	my $headers_row = 0; # Row that contains the header names
	my $index_col = 0; # Column that contains indexes (values that are different for every data row)
	my %headers;
	foreach my $sheet (@{$excel -> {Worksheet}}) {
		my $data_type = lc($sheet->{Name});
		if ($verbose){
			printf("\tReading Sheet '%s'\n", $sheet->{Name});
		}
		$sheet -> {MaxRow} ||= $sheet -> {MinRow};
		foreach my $row ($sheet -> {MinRow} .. $sheet -> {MaxRow}) {
			my $index;
			my %row_data;
			$sheet -> {MaxCol} ||= $sheet -> {MinCol};
			foreach my $col ($sheet -> {MinCol} ..  $sheet -> {MaxCol}) {
				my $cell_data = $sheet -> {Cells} [$row] [$col];
				my $cell = $cell_data->{Val};
				if ($row == $headers_row) {
					if (defined($cell)){
						$headers{$col}=lc($cell);
					} else {
						last; # Doesn't read data after last empty column (ex. Abreviations, Definitions...)
					}
				} elsif ($col == $index_col && defined($cell)){
					$index = lc($cell);
				} elsif (defined($headers{$col}) && defined($cell)){
					$row_data{$headers{$col}}=$cell;
				}
			}
			if (defined($index) && %row_data){
				$oncodata->{$data_type}{$index} = \%row_data;
			}
		}
	}

	return $oncodata;

}

#################################################################################

# Compares variants with reference sequences to retrieve full gene variant genotype description:
sub retrieve_variants {

	my ($align_data, $ref_seqs, $md5_to_depth, $params) = @_;
	
	my $variants;
	
	if (!defined($params)){
		$params = {};
	}

	# Changes the order of $align_data, Ensembl references will be the keys of the hash
	my %references_depth;
	foreach my $md5 (keys %$align_data){
		foreach my $result (@{$align_data->{$md5}}) {
			my $ref_name = $result->{'NAME'};
			$result->{'NAME'} = $md5;
			$result->{'COLS'} = (split("\n",$result->{'COLS'}))[1]."\n".(split("\n",$result->{'COLS'}))[0];
			$result->{'ALIGN'} = (split("\n",$result->{'ALIGN'}))[1]."\n".(split("\n",$result->{'ALIGN'}))[0];
			push(@{$align_data->{$ref_name}}, $result);
			if (defined($md5_to_depth->{$md5})){
				$references_depth{$ref_name} += $md5_to_depth->{$md5};
			} else {
				$references_depth{$ref_name}++;
			}
		}
		delete($align_data->{$md5});
	}
	my @sorted_references = sort { $references_depth{$b} <=> $references_depth{$a} } keys %references_depth;

	# Loops alignments to annotate mutations
	foreach my $ref_name (@sorted_references){

		my $ref_seq = $ref_seqs->{$ref_name};
		my @ref_codons = $ref_seq =~ /\w{3}/g;

		my $ref_variants;
		my %ref_coverage;
		my %var_terminal;
		my %var_insertions;
		foreach my $result (@{$align_data->{$ref_name}}) {
			if ($result->{'ALIGNED'} == $result->{'IDENT'}){
				my $var_md5 = $result->{'NAME'};
				my ($ref_aligned_cols, $var_aligned_cols) = split("\n",$result->{'COLS'});
				my @ref_aligned_cols = split(",",$ref_aligned_cols);
				# Annotates the depth in every position of the alignment, aligned columns are the real positions in the reference sequence
				if (defined($md5_to_depth->{$var_md5})){
					map $ref_coverage{$_} += $md5_to_depth->{$var_md5}, @ref_aligned_cols;
				} else {
					map $ref_coverage{$_}++, @ref_aligned_cols;
				}
				next;
			} else {
				my $var_md5 = $result->{'NAME'};
# if ($ref_name eq 'ENST00000580074.1' && $var_md5 eq 'SRR953264.20270'){
# print '';
# }
				# Annotates sequence errors respect to the reference (2nd function argument)
				my ($ref_aligned_seq, $var_aligned_seq) = split("\n",$result->{'ALIGN'});
				my ($substitutions, $insertions, $deletions) = detect_sequence_errors($var_aligned_seq,$ref_aligned_seq);
				# Discards wrong alignments
				if (scalar @{[@$substitutions, @$insertions, @$deletions]} > 2){
					next;
				}
				my @ref_aligned_seq = split("",$ref_aligned_seq);
				my @var_aligned_seq = split("",$var_aligned_seq);
				my ($ref_aligned_cols, $var_aligned_cols) = split("\n",$result->{'COLS'});
				my @ref_aligned_cols = split(",",$ref_aligned_cols);
				my @var_aligned_cols = split(",",$var_aligned_cols);
				my $first_ref_aligned_col = $ref_aligned_cols[0] - 1;
				# Annotates the depth of every error among all variants
				if (defined($md5_to_depth->{$var_md5})){
					map $ref_variants->{$_+$first_ref_aligned_col}{sprintf("%s->%s",$ref_aligned_seq[$_-1],$var_aligned_seq[$_-1])} += $md5_to_depth->{$var_md5}, @$substitutions;
					map $ref_variants->{$_+$first_ref_aligned_col}{'insertion'} += $md5_to_depth->{$var_md5}, @$insertions;
					map $ref_variants->{$_+$first_ref_aligned_col}{'deletion'} += $md5_to_depth->{$var_md5}, @$deletions;
					# Annotates the depth in every position of the alignment, aligned columns are the real positions in the reference sequence
					map $ref_coverage{$_} += $md5_to_depth->{$var_md5}, @ref_aligned_cols;
				} else {
					map $ref_variants->{$_+$first_ref_aligned_col}{sprintf("%s->%s",$ref_aligned_seq[$_-1],$var_aligned_seq[$_-1])}++, @$substitutions;
					map $ref_variants->{$_+$first_ref_aligned_col}{'insertion'}++, @$insertions;
					map $ref_variants->{$_+$first_ref_aligned_col}{'deletion'}++, @$deletions;
					map $ref_coverage{$_}++, @ref_aligned_cols;
				}
# if ($ref_name eq 'ENST00000412167.6' && defined($ref_variants->{1865})){
# print '';
# }
				# Annotates nts in insertion positions
				map $var_insertions{$_+$first_ref_aligned_col} = $var_aligned_seq[$_], @$insertions;
				# Annotates errors in read terminal positions
				if (defined($params->{'min_variant_position'})){
					foreach my $pos (unique(@$substitutions,@$insertions,@$deletions)){
						if ($pos < $params->{'min_variant_position'} || $pos > (scalar @var_aligned_seq - $params->{'min_variant_position'})){
							$var_terminal{$pos+$first_ref_aligned_col}++;
						}
					}
				}
			}
		}
# 		my $var_length = sprintf("%.0f", mean(@var_lengths));
		my (@deletion, @insertion, $deletion_comments, $insertion_comments);
		foreach my $pos (sort {$a<=>$b} keys %$ref_variants){
			foreach my $error_type (keys %{$ref_variants->{$pos}}) {
				# Annotates previous deletions/insertions
				if (@deletion && ($error_type ne 'deletion' || $pos != $deletion[-1]+1)){
					my $del_first_pos = $deletion[0];
					my $del_len = scalar @deletion;
					my $del_seq =  substr($ref_seq,$del_first_pos-1,scalar @deletion);
					my $del_fist_codon = $ref_codons[int(($del_first_pos-1)/3)];
					my $del_first_aa = dna_to_prot($del_fist_codon);
					my $del_first_pos_aa = int(($del_first_pos+2)/3);
					my $variant_;
					$variant_->{'mutation'} = sprintf("del_%d-%d(%s) (p.%s%d?del)",$del_first_pos,$del_first_pos+$del_len,$del_seq,$del_first_aa,$del_first_pos_aa);
					$variant_->{'position'} = $del_first_pos;
					$variant_->{'coverage'} = sprintf("%d/%d",$ref_variants->{$del_first_pos}{'deletion'},$ref_coverage{$del_first_pos});
					if (defined($deletion_comments)) {
						$variant_->{'comments'} = $deletion_comments;
					}
					push(@{$variants->{$ref_name}}, $variant_);
					undef(@deletion);
					undef($deletion_comments);
				}
				if (@insertion && ($error_type ne 'insertion' || $pos != $insertion[-1]+1)){
					my $ins_first_pos = $insertion[0];
					my $ins_len = scalar @insertion;
					my $ins_seq;
					map $ins_seq .= $var_insertions{$_}, @insertion;
					my $ins_fist_codon = $ref_codons[int(($ins_first_pos-1)/3)];
					my $ins_first_aa = dna_to_prot($ins_fist_codon);
					my $ins_first_pos_aa = int(($ins_first_pos+2)/3);
					my $variant_;
					$variant_->{'mutation'} = sprintf("ins_%d-%d(%s) (p.%s%d?ins)",$ins_first_pos,$ins_first_pos+$ins_len,$ins_seq,$ins_first_aa,$ins_first_pos_aa);
					$variant_->{'position'} = $ins_first_pos;
					$variant_->{'coverage'} = sprintf("%d/%d",$ref_variants->{$ins_first_pos}{'insertion'},$ref_coverage{$ins_first_pos});
					if (defined($insertion_comments)) {
						$variant_->{'comments'} = $insertion_comments;
					}
					push(@{$variants->{$ref_name}}, $variant_);
					undef(@insertion);
					undef($insertion_comments);
				}
				my @comments;
				# Excludes rare variant mutations if desired:
				# Variants with lower frequency than threshold
				if (defined($params->{'min_variant_frequency'}) && $ref_variants->{$pos}{$error_type} < $params->{'min_variant_frequency'}/100*$ref_coverage{$pos}){
					if (defined($params->{'rare_variants'}) && !$params->{'rare_variants'}){
						next;
					}
					push(@comments, sprintf('Low frequency variant (%.2f%%).',100*$ref_variants->{$pos}{$error_type}/$ref_coverage{$pos}));
				}
				# Variants with lower coverage than threshold
				if (defined($params->{'min_variant_coverage'}) && $ref_variants->{$pos}{$error_type} < $params->{'min_variant_coverage'}){
					if (defined($params->{'rare_variants'}) && !$params->{'rare_variants'}){
						next;
					}
					push(@comments, sprintf('Low coverage variant (%d reads).',$ref_variants->{$pos}{$error_type}));
				}
				# Variants mapped by mostly terminal positions of reads
				if (defined($params->{'min_variant_position'}) && $var_terminal{$pos} > 0.9*$ref_variants->{$pos}{$error_type} && $var_terminal{$pos} < 0.2*$ref_coverage{$pos} ){
					if (defined($params->{'rare_variants'}) && !$params->{'rare_variants'}){
						next;
					}
					push(@comments, sprintf('Terminal mapping variant (%d/%d).',$var_terminal{$pos},$ref_variants->{$pos}{$error_type}));
				}
				if ($error_type =~ /(\w)->(\w)/){
					my $ref_nt = $1;
					my $var_nt = $2;
					my $ref_codon = $ref_codons[int(($pos-1)/3)];
					my $ref_aa = dna_to_prot($ref_codon);
# 					my $mut_seq = substr($ref_seq,0,$pos-1).$var_nt.substr($ref_seq,$pos+1);
# 					my @mut_codons = $mut_seq =~ /\w{3}/g;
# 					my $mut_codon = $mut_codons[int(($pos-1)/3)];
# 					my $mut_aa = dna_to_prot($mut_codon);
					my $mut_codon = $ref_codon;
					substr($mut_codon,($pos-1) % 3,1) = $2;
					my $mut_aa = dna_to_prot($mut_codon);
					my $pos_aa = int(($pos+2)/3);
					my $variant_;
					if ($mut_aa ne $ref_aa){
						$variant_->{'mutation'} = sprintf("%d %s->%s (p.%s%d%s)", $pos, $ref_nt, $var_nt, convert_aa_1_to_3($ref_aa), $pos_aa, convert_aa_1_to_3($mut_aa));
						$variant_->{'variant'} = sprintf("%s%d%s", $ref_aa, $pos_aa, $mut_aa);
						$variant_->{'type'} = "non-synonymous substitution";
					} else {
						# Skips synonymous substitutions
						next;
						$variant_->{'mutation'} = sprintf("%d %s->%s (p.%s%d=)", $pos, $ref_nt, $var_nt, convert_aa_1_to_3($ref_aa), $pos_aa);
						$variant_->{'type'} = "synonymous substitution";
					}
					$variant_->{'position'} = $pos;
					$variant_->{'coverage'} = sprintf("%d/%d",$ref_variants->{$pos}{$error_type},$ref_coverage{$pos});
					if (@comments) {
						$variant_->{'comments'} = join("; ",@comments);
					}
					push(@{$variants->{$ref_name}}, $variant_);
				} elsif ($error_type eq 'deletion'){
					# next;
					push(@deletion,$pos);
					if (@comments) {
						$deletion_comments = join("; ",@comments);
					}
				} elsif ($error_type eq 'insertion'){
					# next;
					push(@insertion,$pos);
					if (@comments) {
						$insertion_comments = join("; ",@comments);
					}
				} else {
					next;
				}
			}
		}
	}

	# print Dumper($variants);
	# exit;

	return $variants;

}


#################################################################################

# Retrieves variant mutations directly from a BAM/SAM file, faster than running extract_align_bam_data()+retrieve_variants()
sub retrieve_variants_from_bam_data {

	my ($infile,$ref_sequences_file,$ref_annotations_file,$params) = @_;
	
	my $ref_variants = { 'all' => {}, 'filtered' => {}, 'coding' => {}, 'annotations' => {} };
	my $ref_mapped_lengths = { 'all' => 0, 'filtered' => 0, 'coding' => 0 };

	if (!defined($params)){
		$params = {};
	}

	# SAM data example (samtools view alignment.bam | less -S):
	# 232d9c3aa0b15a91b46eaa4642c79193        16      ENST00000265620.11      537     1       15S73M  *       0       0       GACCTCAATTTTGTTTCAGGACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAAC        IIIIIIIIIIIIII
	# 8bd2c52cadac2ef8a3ff5154d013e204        16      ENST00000265620.11      538     1       72M     *       0       0       CAGGACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAAC        IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	# 4cb9965bd13db4572fc14c22e4266510        0       ENST00000371085.7       587     1       61M     *       0       0       ACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAA   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	# 39a367591c84a3f373c7df2ba38b0b41        16      ENST00000371095.7       541     1       77M4S   *       0       0       CAGGACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAACTTCCAGTAA       IIIIIIIIIIIIIIIIIIIIII
	# e775a2ee8a41d2b252d4c92d5e8038b6        0       ENST00000371095.7       545     1       72M     *       0       0       ACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAACTTCC        IIIIIIIIIIIIIIIIIIIIIIIIIIIIII
	# 0c463aad1591beb227f46282763d5bc0        0       ENST00000371102.8       2469    1       11S77M  *       0       0       TCAATTTTGTTTCAGGACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAACTTCC        IIIIIIIIIIIIII
	# 5e3f383ccd14705fdd80535369d4c558        0       ENST00000371102.8       2483    1       63M     *       0       0       GCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAACTTCC IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

	# Opens the SAM/BAM file, it should be in the same order than the reference genome/sequences
	print ("\tProcessing alignment file '$infile'.\n");
	if (is_sam($infile)){
		open(ALIGN_FILE,$infile) || die "# $0 : # 'retrieve_variants_from_bam_data' cannot read '$infile'.\n";
	} elsif (is_bam($infile)){
		open(ALIGN_FILE,"$SAMTOOLS view $infile |") || die "# $0 : # 'retrieve_variants_from_bam_data' cannot read '$infile'.\n";
	} else {
		print "# 'retrieve_variants_from_bam_data': Unrecognized BAM or SAM format file '$infile'.\n";
		exit;
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
	# Writes each entry in the input to the output file
	# Ex:
	# >ENST00000632684.1 cds chromosome:GRCh38:7:142786213:142786224:1 gene:ENSG00000282431.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene gene_symbol:TRBD1 description:T-cell receptor beta diversity 1 [Source:HGNC Symbol;Acc:HGNC:12158]
	# GGGACAGGGGGC....
	#
	# $line[18] is the MD tag: String for mismatching positions.The MD field aims to achieve SNP/indel calling without looking at the reference.  For example, a string ‘10A5^AC6’ means  from the  leftmost reference  base in the  alignment,  there  are 10 matches followed by an A on the reference which is different from the aligned read base; the next 5 reference bases are matches followed by a 2bp deletion from the reference; the deleted sequence is AC; the last 6 bases are matches. The MD field ought to match the CIGAR string
	# But MD tag doesn't give information about insertions
	# 15S73M  ....  AS:i:140        XS:i:140        XN:i:0  XM:i:1  XO:i:0  XG:i:0  NM:i:1  MD:Z:61G11      YT:Z
	# GACCTCAATTTTGTTTCAGGACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGAACAAAGTCAAC 
	#                                                                             -
	#                TCAGGACCTGCTTCGCTGCCGTGTCCTGACTTCTGGAATCTTTGAGACCAAGTTCCAGGTGGACAAAGTCAAC

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

	my (%ref_coverage,%var_terminal,%var_insertions);
	my $ref_obj;
	my $ref_id = '';
	my $read_count = 0;
	while(<ALIGN_FILE>){
		# s/\012\015?|\015\012?//;
		my @line = split("\t");
		if (scalar @line < 10){
			next;
		} elsif ($line[2] eq '*'){ # Empty sbjct name - no alignment
			next;
		}
		my $seq_name = $line[0];
		my $ref_name = $line[2];
		my $cigar = $line[5];
		my $ref_pos = $line[3];
		# Annotates results from previous reference
		# It's not neccessary to annotate again at the end of the while loop, because after chromosomes there are more sequences in the genome file (ex. KI270392.1)
		if ($ref_id ne $ref_name && defined($ref_variants->{'all'}{$ref_id})){
			# Annotates total reference mapped length
			$ref_mapped_lengths->{'all'} += scalar keys %ref_coverage;
			# Filters low coverage/frequency mutations found in the previous reference sequence
			$ref_variants->{'filtered'}{$ref_id} = filter_variants($ref_variants->{'all'}{$ref_id},\%ref_coverage,\%var_terminal,$params);
			# Annotates variants
			if (%{$ref_variants->{'filtered'}{$ref_id}}){
				# Annotates filtered variants reference mapped length
				$ref_mapped_lengths->{'filtered'} += scalar keys %ref_coverage;
				# Detects and annotates coding region variants
				my ($ref_variants_coding,$ref_variants_annotations,$ref_coding_mapped_length) = annotate_reference_coding_variants($ref_variants->{'filtered'}{$ref_id},$ref_obj,$gff_file_handle,\%ref_coverage,\%var_insertions,$params);
				$ref_variants->{'coding'}{$ref_id} = $ref_variants_coding;
				$ref_mapped_lengths->{'coding'} += $ref_coding_mapped_length;
				if (%$ref_variants_annotations) {
					$ref_variants->{'annotations'} = { %{$ref_variants->{'annotations'}},  %$ref_variants_annotations };
				}
			}
			# Annotates variant frequencies instead of coverage
			foreach my $type (('all','filtered','coding')) {
# if (%{$ref_variants->{'filtered'}{$ref_id}}){
# print '';
# }
				foreach my $pos (keys %{$ref_variants->{$type}{$ref_id}}){
					foreach my $mut (keys %{$ref_variants->{$type}{$ref_id}{$pos}}){
						$ref_variants->{$type}{$ref_id}{$pos}{$mut} = sprintf("%d/%d", $ref_variants->{$type}{$ref_id}{$pos}{$mut},$ref_coverage{$pos});
					}
				}
			}
			undef(%ref_coverage);
			undef(%var_insertions);
			undef(%var_terminal);
		}
		# Goes to next reference sequence (eg. next chromosome) if the alignments of the previous one have finished
		while ($ref_id ne $ref_name){
			$ref_obj = $reference_data->next_seq;
			$ref_id = $ref_obj->primary_id();
		}
# if ($seq_name eq 'SRR953264.37100' && $ref_name  eq 'ENST00000412167.6'){
# print '';
# }
		my ($mismatches,$indels) = (0,0);
		# Detects mismatches in SAM format additional fields
		if (/XM:i:(\d+)/){
			$mismatches = $1;
		}
		# Detects indels (gap openings) in SAM format additional fields
		if (/XO:i:(\d+)/){
			$indels = $1;
		}
		# Continues annotation only if there are not mismatches or gap openings in the alignment
		if (defined($params->{'max_errors'}) && $mismatches+$indels > $params->{'max_errors'}){
			next;
		}
		# Annotates the depth in the aligned positions
		if (!$mismatches && !$indels){
			# Annotates the depth in every position of the alignment, aligned columns are the real positions in the reference sequence
			my $seq_pos = 1;
			while ($cigar =~ /(\d+)([MSH=X]+)/g){
				foreach (1..$1){
					if ($2 eq 'M' || $2 eq '=' || $2 eq 'X'){
						$ref_coverage{$ref_pos}++;
						$ref_pos++;
					} elsif ($2 eq 'S') {
						# Soft-clipping: bases in 5' and 3' of the read are NOT part of the alignment.
						$seq_pos++;
					} elsif ($2 eq 'H') {
						# Hard-clipping: bases in 5' and 3' of the read are NOT part of the alignment AND those bases have been removed from the read sequence in the BAM file. The 'real' sequence length would be length(SEQ)+ count-of-hard-clipped-bases
						$seq_pos++;
					}
				}
			}
			# next;
		} else {
			# Annotates the depth in every position of the alignment, aligned columns are the real positions in the reference sequence
			my $seq_seq = $line[9];
			my $seq_pos = 1;
			my $new_ins = 0;
			my $seq_len = 0;
			$seq_len += $1 while $cigar =~ /(\d+)([MIDNSHP=X]+)/g;
			while ($cigar =~ /(\d+)([MIDNSHP=X]+)/g){
				# Ex. CIGAR: 3M1I3M1D5M
				# CCATACT-GAACTGACTAAC
				#     ACTAGAA-TGGCT
				foreach (1..$1){
					my $var_pos = 0;
					if ($2 eq 'M' || $2 eq '=' || $2 eq 'X'){
						my $ref_nt = $ref_obj->subseq($ref_pos,$ref_pos);
						my $seq_nt = substr($seq_seq,$seq_pos-1,1);
						if ($ref_nt ne $seq_nt){
							# Annotates the mutation coverage in the reference position
							$ref_variants->{'all'}{$ref_name}{$ref_pos}{sprintf("%s->%s",$ref_nt,$seq_nt)}++;
							$var_pos = $seq_pos;
						}
# # 						# Annotates the total reference mapped length (to calculate TMB)
# 						if (!defined($ref_coverage{$ref_pos})){
# 							$ref_mapped_lengths->{'all'}++;
# 						}
						# Annotates the total reference coverage
						$ref_coverage{$ref_pos}++;
						$new_ins = 0;
						$seq_pos++;
						$ref_pos++;
					} elsif ($2 eq 'I') {
# if ($ref_name eq '9' && $ref_pos == 130872112){
# print Dumper(@line)."\n";
# print '';
# }
						if (defined($params->{'indels'}) && $params->{'indels'}){
							# Skips homopolymer sequencing errors
							if (substr($seq_seq,$seq_pos-4,6) !~ /[ACTG]{3}/){
								# Annotates the insertion coverage only once in the reference insertion position
								if (!$new_ins){
									$ref_variants->{'all'}{$ref_name}{$ref_pos}{'insertion'}++;
									$var_pos = $seq_pos;
									$new_ins = $ref_pos;
								}
								# Annotates nts in insertion positions
								$var_insertions{$new_ins} .= substr($seq_seq,$seq_pos-1,1);
							}
						}
						$seq_pos++;
					} elsif ($2 eq 'D') {
# if ($ref_name eq '9' && $ref_pos == 130872112){
# print Dumper(@line)."\n";
# print '';
# }
						if (defined($params->{'indels'}) && $params->{'indels'}){
							# Skips homopolymer sequencing errors
							if (substr($seq_seq,$seq_pos-4,6) !~ /[ACTG]{3}/){
								# Annotates the deletion coverage in the reference deleted positions
								$ref_variants->{'all'}{$ref_name}{$ref_pos}{'deletion'}++;
								$var_pos = $seq_pos;
								$new_ins = 0;
							}
						}
						$ref_pos++;
					} elsif ($2 eq 'S' || $2 eq 'H') {
						# Soft-clipping: bases in 5' and 3' of the read are NOT part of the alignment.
						$seq_pos++;
					}
					# Annotates errors in read terminal positions
					if (defined($params->{'min_variant_position'}) && $var_pos) {
						if ($var_pos < $params->{'min_variant_position'} || $var_pos > ($seq_len - $params->{'min_variant_position'})){
							$var_terminal{$ref_pos}++;
						}
					}
				}
			}
		}
		$read_count++;
		if (defined($params->{'limit'}) && $params->{'limit'}<=$read_count){
			last;
		}
	}
	close(ALIGN_FILE);
	if (defined($gff_file_handle)){
		close $gff_file_handle;
	}

	# Stores genome positions from variants
	my $variant_genome_pos = {};

	# Counts total unique genome positions of variants
	foreach my $ref_id (keys %{$ref_variants->{'all'}}){
		foreach my $ref_pos (keys %{$ref_variants->{'all'}{$ref_id}}){
			$variant_genome_pos->{'all'}{"$ref_id:$ref_pos"} = 1;
		}
	}

	# Counts unique genome positions of filtered variants
	foreach my $ref_id (keys %{$ref_variants->{'filtered'}}){
		foreach my $ref_pos (keys %{$ref_variants->{'filtered'}{$ref_id}}){
			$variant_genome_pos->{'filtered'}{"$ref_id:$ref_pos"} = 1;
		}
	}

	# Counts unique genome positions of filtered variants coding regions
	foreach my $ref_id (keys %{$ref_variants->{'coding'}}){
		foreach my $ref_pos (keys %{$ref_variants->{'coding'}{$ref_id}}){
			$variant_genome_pos->{'coding'}{"$ref_id:$ref_pos"} = 1;
		}
	}
	
	# print Dumper($ref_variants_all);
	# print Dumper($ref_variants->{'filtered'});
	# print Dumper($ref_variants->{'coding'});

# 	print "\n";
# 	printf("\t%d reads mapped to references\n",$read_count);
# 	printf("\tTotal length mapped: %d bp.\n",$ref_mapped_lengths->{'all'});
# 	printf("\t%d mutated positions found.\n", scalar keys %variant_total_genome_pos);
# 	printf("\t%d mutated positions after filtering low frequency variants.\n", scalar keys %variant_filtered_genome_pos);
# 	printf("\t%d mutated positions found in coding regions.\n", scalar keys %variant_coding_genome_pos);
# 	printf("\tTotal coding length mapped: %d bp.\n", $ref_mapped_lengths->{'coding'});
	# exit;
	print "\n";
	printf("\t%d reads mapped to references\n",$read_count);
	printf("\tTotal length mapped: %d bp.\n",$ref_mapped_lengths->{'all'});
	printf("\tTotal mutated positions: %d.\n", scalar keys %{$variant_genome_pos->{'all'}});
	printf("\tLength mapped after filtering low frequency variants: %d bp.\n",$ref_mapped_lengths->{'filtered'});
	printf("\tMutated positions after filtering low frequency variants: %d.\n", scalar keys %{$variant_genome_pos->{'filtered'}});
	printf("\tLength mapped by coding regions: %d bp.\n",$ref_mapped_lengths->{'coding'});
	printf("\tMutated positions found in coding regions: %d.\n", scalar keys %{$variant_genome_pos->{'coding'}});



	return ($ref_variants, $ref_mapped_lengths);

}




#################################################################################

# Annotates variant mutations in coding regions
sub annotate_reference_coding_variants {

	my ($ref_variants,$ref_obj,$gff_file_handle,$ref_coverage,$var_insertions,$params) = @_;

	# Variable to store annotated coding variants
	my $ref_coding_variants = {};
	my $ref_coding_variants_annotations = {};

	# Variable to store mapped length in coding regions
	my $ref_coding_mapped_length = 0;
	
	my $ref_id = $ref_obj->primary_id();

	# If the reference file used for mapping is an annotated genome
	if (defined($gff_file_handle)) {

		my %gene_mutations;
		my $ref_to_transcript_variants;
		# Variables to annotate coding regions
		my ($gene_id,$gene_coords,$gene_name,$gene_description,$transcript_id);
		my $ref_found = 0;
		my $read_transcript = 0;
		my $transcript_seq = '';
		my $transcript_first_pos = 0;
		my $transcript_last_pos = 0;
		my $transcript_intron_count = 0;
		while (<$gff_file_handle>) {
			if (/^#/){
				next;
			}
			# s/\012\015?|\015\012?//; # trim line break
			my @gff_cols = split("\t");
			# if ($#gff_cols < 8){
				# next;
			# }
			# Reads the GFF file till it finds the current reference
			if ($ref_id eq $gff_cols[0]){
				$ref_found = 1;
			# Finish annotation
			} elsif ($ref_found && $ref_id ne $gff_cols[0]){
				# Goes back one GFF line to the first line of the next reference
				# seek($gff_file_handle, -length($_), 1); # Not working code, but is not neccessary
				last;
			} else {
				next;
			}

			# Skips genes without variants
			# if (!defined($ref_variants)){
				# next;
			# }
			# If we desire only sequences with manual annotation (HAVANA) to avoid automatic transcript annotations
			if ($gff_cols[1] !~ /havana/i){
				next;
			}
			my $ref_name = $gff_cols[0];
			my $type = $gff_cols[2];
			# Gene example:
			# 7       ensembl_havana  gene    55019021        55211628        .       +       .       ID=gene:ENSG00000146648;Name=EGFR;biotype=protein_coding;description=epidermal growth factor receptor [Source:HGNC Symbol%3BAcc:HGNC:3236];gene_id=ENSG00000146648;logic_name=ensembl_havana_gene;version=17
			if ($type eq 'gene' || $type eq 'mt_gene'){
				# Reset variables with gene info
				$read_transcript = 0;
				undef($gene_id);
				undef($gene_name);
				undef($gene_description);
				undef($gene_coords);
				undef(%gene_mutations);
				my @attrs = split( ";", $gff_cols[8] );
				# Only protein coding genes, exclude pseudogenes and highly polimorphyc ones
				if ($attrs[2] !~ /protein_coding/ || $gene_description =~ /olfactory receptor|T-cell receptor|B-cell receptor|immunoglobulin|major histocompatibility complex/){
					next;
				}
		

				if ($gff_cols[8] =~ /^ID=(gene:)?(\w+)/){
					$gene_id = $2;
				}
				if ($gff_cols[8] =~ /Name=(.+?);/){
					$gene_name = $1;
				}
				if ($gff_cols[8] =~ /description=(.+?);/){
					$gene_description = $1;
				}
				$gene_coords = sprintf("chr%s:%d-%d",$ref_name,$gff_cols[3],$gff_cols[4]);
# if ($gene_id eq 'ENSG00000156920'){
# print '';
# }
				foreach my $mut_pos (keys %$ref_variants){
					if ($mut_pos >= $gff_cols[3] && $mut_pos <= $gff_cols[4]) {
						$gene_mutations{$mut_pos} = 1;
						$read_transcript = 1;
					}
				}
			# CDS example:
			# 7       havana  CDS     55019278        55019365        .       +       0       ID=CDS:ENSP00000415559;Parent=transcript:ENST00000455089;protein_id=ENSP00000415559
			} elsif ($read_transcript && $type eq 'CDS' ) {
				if (!defined($transcript_id) && $gff_cols[8] =~ /Parent=(transcript:)?(\w+)/){
					$transcript_id = $2;
				}
if ($transcript_id eq 'ENST00000394141'){
print '';
}

				$transcript_seq .= $ref_obj->subseq($gff_cols[3],$gff_cols[4]);
				if ($transcript_last_pos){
					$transcript_intron_count += $gff_cols[3] - $transcript_last_pos -1;
				} else {
					$transcript_first_pos = $gff_cols[3];
				}
				# Assigns transcript relative positions to the mutations
				foreach my $mut_pos (keys %gene_mutations){
					if ($mut_pos < $gff_cols[3] || $mut_pos > $gff_cols[4]) {
						next;
					}
					# Annnotates mutation position, type an coverage
					$ref_coding_variants->{$mut_pos} = { %{$ref_variants->{$mut_pos}} };
					# If it is reverse strand, annotates negative position in $ref_to_transcript_variants
					if ($gff_cols[6] eq '-'){
						$ref_to_transcript_variants->{-$mut_pos} = $mut_pos - $transcript_first_pos - $transcript_intron_count;
					} else {
						$ref_to_transcript_variants->{$mut_pos} = $mut_pos - $transcript_first_pos - $transcript_intron_count;
					}
				}
				$transcript_last_pos = $gff_cols[4];
				# Annotates the coding region mapped
				foreach my $pos ($gff_cols[3]..$gff_cols[4]){
					if (defined($ref_coverage->{$pos})){
						$ref_coding_mapped_length++;
					}
				}
				
			# mRNA example:
			# 7       havana  mRNA    55019021        55203076        .       +       .       ID=transcript:ENST00000455089;Parent=gene:ENSG00000146648;Name=EGFR-207;biotype=protein_coding;tag=basic;transcript_id=ENST00000455089;transcript_support_level=1;version=5
			# five_prime_UTR example:
			# 7       havana  five_prime_UTR  55019021        55019277        .       +       .       Parent=transcript:ENST00000455089
			} elsif ($read_transcript && ($type eq 'five_prime_UTR' || $type eq 'mRNA' || $type eq 'NMD_transcript_variant')) {
if ($transcript_id eq 'ENST00000394141'){
print '';
}
				my $transcript_len = length($transcript_seq);

 				# Reverse transcript sequence if it is reverse strand
				if ($gff_cols[6] eq '-'){
					$transcript_seq = iupac_reverse_complementary($transcript_seq);
					foreach my $mut_pos (keys %$ref_to_transcript_variants){
						$ref_to_transcript_variants->{$mut_pos} = $transcript_len - $ref_to_transcript_variants->{$mut_pos};
					}
				}

				# Annotates mutations
				if (defined($ref_to_transcript_variants)){
					my $ref_coding_variants_annotations_ = retrieve_transcript_variants($gene_id,$transcript_id,$transcript_seq,$ref_to_transcript_variants,$ref_name,$ref_variants,$ref_coverage,$var_insertions,$gene_description,$params);
					$ref_coding_variants_annotations = { %$ref_coding_variants_annotations, %$ref_coding_variants_annotations_ };
				}

				# Reset variables with transcript info
				undef($ref_to_transcript_variants);
				undef($transcript_id);
				$transcript_seq = '';
				$transcript_first_pos = 0;
				$transcript_last_pos = 0;
				$transcript_intron_count = 0;
			}
		}

	# If the reference file used for mapping is not annotated (Ex. CDS sequences)
	} else {

		# Sequence example:
		# >ENST00000421351.3 cds chromosome:GRCh38:6:138088505:138107511:-1 gene:ENSG00000112378.11 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:PERP description:PERP, TP53 apoptosis effector [Source:HGNC Symbol;Acc:HGNC:17637]
		# $ref_id = 'ENST00000421351.3'
		# $ref_description = 'cds chromosome:GRCh38:6:138088505:138107511:-1 gene:ENSG00000112378.11 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:PERP description:PERP, TP53 apoptosis effector [Source:HGNC Symbol;Acc:HGNC:17637]'

		# Skips genes without variants
		# if (!defined($ref_variants)){
			# next;
		# }

		# Reference name
		my $ref_name = $ref_id;

		# Reference descriptions
		my $ref_description = $ref_obj->desc();
		# Only protein coding genes, exclude pseudogenes and highly polimorphyc ones
		if ($ref_description =~ /gene_biotype:protein_coding/ && $ref_description =~ /gene:(\w+)/
		&& $ref_description !~/olfactory receptor|T-cell receptor|B-cell receptor|immunoglobulin|major histocompatibility complex/){
		# $ref_description !~ /gene_biotype:(IG_\w_gene|TR_\w_gene|polymorphic_pseudogene|)/
			my $gene_id = $1;
			my $transcript_id;
			# Annotates Ensembl ID without version number: ENST00000390588.1 -> ENST00000390588
			if ($ref_name =~ /^(ENST\d+)\.\d+/){
				$transcript_id = $1;
			} else {
				$transcript_id = $ref_name;
			}
			my $transcript_seq = $ref_obj->seq();;

			# All the transcript positions will be the same than references (the same transcripts)
			my $ref_to_transcript_variants;
			# Assigns transcript relative positions to the mutations (the same than mapped)
			foreach my $mut_pos (keys %$ref_variants){
				$ref_to_transcript_variants->{$mut_pos} = $mut_pos;
				# Annnotates mutation position, type an coverage (copies hash as a value, dereferencing)
				$ref_coding_variants->{$mut_pos} = { %{$ref_variants->{$mut_pos}} };
			}

			# Annotates mutation detailed information
			my $ref_coding_variants_annotations_ = retrieve_transcript_variants($gene_id,$transcript_id,$transcript_seq,$ref_to_transcript_variants,$ref_name,$ref_variants,$ref_coverage,$var_insertions,$ref_description,$params);
			$ref_coding_variants_annotations = { %$ref_coding_variants_annotations, %$ref_coding_variants_annotations_ };
			
			# Annotates the coding region mapped
			$ref_coding_mapped_length = scalar keys %$ref_coverage;
		}

	}
	
	return ($ref_coding_variants,$ref_coding_variants_annotations,$ref_coding_mapped_length);
}

#################################################################################
	
# Filters low coverage/frequency variants
sub filter_variants {

	my ($variants,$coverage,$terminal,$params) = @_;

	my $filtered_variants = {};

	#my $comments;
	foreach my $mut_pos (keys %$variants){
		foreach my $error_type (keys %{$variants->{$mut_pos}}){
			my $total_coverage = $coverage->{$mut_pos}; #+$variants->{$mut_pos}{$error_type};
			my $var_coverage = $variants->{$mut_pos}{$error_type};
			# Excludes rare variant mutations if desired:
			my $keep_mut = 1;
			# Variants with lower frequency than threshold
			if (defined($params->{'min_variant_frequency'}) && $var_coverage < $params->{'min_variant_frequency'}/100*$total_coverage){
				unless (defined($params->{'rare_variants'}) && $params->{'rare_variants'}){
					$keep_mut = 0;
				} # else {
					# push(@{$comments->{$mut_pos}{$error_type}}, sprintf('Low frequency variant (%.2f%%).',100*$var_coverage/$total_coverage));
				# }
			}
			# Variants with lower coverage than threshold
			if (defined($params->{'min_variant_coverage'}) && $var_coverage < $params->{'min_variant_coverage'}){
				unless (defined($params->{'rare_variants'}) && $params->{'rare_variants'}){
					$keep_mut = 0;
				} # else {
					# push(@{$comments->{$mut_pos}{$error_type}}, sprintf('Low coverage variant (%d reads).',$var_coverage));
				# }
			}
			# Variants mapped by mostly terminal positions of reads
			if (defined($params->{'min_variant_position'}) && $terminal->{$mut_pos} > 0.9*$var_coverage && $terminal->{$mut_pos} < 0.2*$total_coverage ){
				unless (defined($params->{'rare_variants'}) && $params->{'rare_variants'}){
					$keep_mut = 0;
				} # else {
					# push(@{$comments->{$mut_pos}{$error_type}}, sprintf('Terminal mapping variant (%d/%d).',$terminal->{$mut_pos},$var_coverage));
				# }
			}
			if ($keep_mut){
				$filtered_variants->{$mut_pos}{$error_type} = $variants->{$mut_pos}{$error_type};
				#delete($variants->{$mut_pos}{$error_type})
			}
# if ($keep_mut && $error_type eq 'insertion'){
# print '';
# }
# if ($keep_mut && $error_type eq 'deletion'){
# print '';
# }

		}
		#if (!%{$variants->{$mut_pos}}){
			#delete($variants->{$mut_pos});
		#}
	}
	
	return $filtered_variants;
	
}

#################################################################################

# Annotates transcript sequence mutations
sub retrieve_transcript_variants {

	my ($gene_id,$transcript_id,$transcript_seq,$ref_to_transcript_mutations,$ref_name,$ref_mutations,$ref_coverage,$var_insertions,$gene_description,$params) = @_;
	
	my $variants = {};

if ($transcript_id eq 'ENST00000368948'){
print '';
}

	# Annotates mutations in each position of the CDS
	my ($deletion, $insertion, $deletion_comments, $insertion_comments);
	my %annotated_variants;
	my @transcript_codons = $transcript_seq =~ /\w{3}/g;
	my $reverse = 0;
	foreach my $mut_pos (keys %$ref_to_transcript_mutations){
		# $mut_pos is negative if transcript is in reverse strand
		# $transcript_mut_pos is always positive
		my $transcript_mut_pos = $ref_to_transcript_mutations->{$mut_pos};
		# Checks if sequence is reverse
		if ($mut_pos < 0){
			$reverse = 1;
		}
		# But positions are absolute in $ref_mutations, they do not indicate the sense of the transcript
		foreach my $error_type (keys %{$ref_mutations->{abs($mut_pos)}}){
			my $total_coverage = $ref_coverage->{abs($mut_pos)}; #+$ref_mutations->{abs($mut_pos)}{$error_type};
			my $var_coverage = $ref_mutations->{abs($mut_pos)}{$error_type};
			# First, annotates previous deletions/insertions
			if (defined($deletion) && @{$deletion->{'reference'}} && ($error_type ne 'deletion' || $transcript_mut_pos != $deletion->{'transcript'}[-1]+1)){
				my $del_first_pos_ref = $deletion->{'reference'}[0];
				my $del_first_pos = $deletion->{'transcript'}[0];
				my $del_len = scalar @{$deletion->{'transcript'}};
				my $del_seq =  substr($transcript_seq,$del_first_pos-1,$del_len);
				my $del_fist_codon = $transcript_codons[int(($del_first_pos-1)/3)];
				my $del_first_aa = dna_to_prot($del_fist_codon);
				my $del_first_pos_aa = int(($del_first_pos+2)/3);
				my %variant;
				$variant{'gene'} = $gene_id;
				$variant{'gene_description'} = $gene_description;
				$variant{'genome_position'} = sprintf("%s:%d",$ref_name,abs($mut_pos));
				$variant{'mutation'} = sprintf("del_%d-%d(%s) (p.%s%d?del)",$del_first_pos,$del_first_pos+$del_len,$del_seq,$del_first_aa,$del_first_pos_aa);
				$variant{'position'} = $del_first_pos;
				$variant{'coverage'} = sprintf("%d/%d",$ref_mutations->{$del_first_pos_ref}{'deletion'},$ref_coverage->{$del_first_pos_ref}); #+$ref_mutations->{$del_first_pos_ref}{'deletion'});
				if (defined($deletion_comments)) {
					$variant{'comments'} = $deletion_comments;
				}
				if (!defined($annotated_variants{$variant{'mutation'}})){
					push(@{$variants->{$transcript_id}}, \%variant);
					$annotated_variants{$variant{'mutation'}} = 1;
				}
				undef($deletion);
				undef($deletion_comments);
			} elsif (defined($insertion) && @{$insertion->{'reference'}} && ($error_type ne 'insertion' || $transcript_mut_pos != $insertion->{'transcript'}[-1]+1)){
				my $ins_first_pos_ref = $insertion->{'reference'}[0];
				my $ins_first_pos = $insertion->{'transcript'}[0];
				my $ins_len = scalar @{$insertion->{'transcript'}};
				my $ins_seq = $var_insertions->{$ins_first_pos_ref};
				my $ins_fist_codon = $transcript_codons[int(($ins_first_pos-1)/3)];
				my $ins_first_aa = dna_to_prot($ins_fist_codon);
				my $ins_first_pos_aa = int(($ins_first_pos+2)/3);
				my %variant;
				$variant{'gene'} = $gene_id;
				$variant{'gene_description'} = $gene_description;
				$variant{'genome_position'} = sprintf("%s:%d",$ref_name,abs($mut_pos));
				$variant{'mutation'} = sprintf("ins_%d-%d(%s) (p.%s%d?ins)",$ins_first_pos,$ins_first_pos+$ins_len,$ins_seq,$ins_first_aa,$ins_first_pos_aa);
				$variant{'position'} = $ins_first_pos;
				$variant{'coverage'} = sprintf("%d/%d",$ref_mutations->{$ins_first_pos_ref}{'insertion'},$ref_coverage->{$ins_first_pos_ref}); #+$ref_mutations->{$ins_first_pos_ref}{'insertion'});
				if (defined($insertion_comments)) {
					$variant{'comments'} = $insertion_comments;
				}
				if (!defined($annotated_variants{$variant{'mutation'}})){
					push(@{$variants->{$transcript_id}}, \%variant);
					$annotated_variants{$variant{'mutation'}} = 1;
				}
				undef($insertion);
				undef($insertion_comments);
			}
			if ($error_type =~ /(\w)->(\w)/){
				my $ref_nt = $1;
				my $var_nt = $2;
				# If it is reverse strand
				if ($reverse){
					$ref_nt = complementary_sequence($ref_nt);
					$var_nt = complementary_sequence($var_nt);
				}
				my $ref_codon = $transcript_codons[int(($transcript_mut_pos-1)/3)];
# if (!defined($ref_codon)){
# print "";
# # next;
# }
				my $ref_aa = dna_to_prot($ref_codon);
				my $mut_codon = $ref_codon;
				substr($mut_codon,($transcript_mut_pos-1) % 3,1) = $var_nt;
				my $mut_aa = dna_to_prot($mut_codon);
				my $transcript_mut_pos_aa = int(($transcript_mut_pos+2)/3);
				my %variant;
				$variant{'gene'} = $gene_id;
				$variant{'gene_description'} = $gene_description;
				$variant{'genome_position'} = sprintf("%s:%d",$ref_name,abs($mut_pos));
				if ($mut_aa ne $ref_aa){
					$variant{'mutation'} = sprintf("%d %s->%s p.%s%d%s (%s%d%s)", $transcript_mut_pos, $ref_nt, $var_nt, convert_aa_1_to_3($ref_aa), $transcript_mut_pos_aa, convert_aa_1_to_3($mut_aa), $ref_aa, $transcript_mut_pos_aa, $mut_aa);
					$variant{'variant'} = sprintf("%s%d%s", $ref_aa, $transcript_mut_pos_aa, $mut_aa);
					$variant{'type'} = "non-synonymous substitution";
				} else {
					next;
					$variant{'mutation'} = sprintf("%d %s->%s p.%s%d= (%s%d=)", $transcript_mut_pos, $ref_nt, $var_nt, convert_aa_1_to_3($ref_aa), $transcript_mut_pos_aa, $ref_aa, $transcript_mut_pos_aa);
					$variant{'type'} = "synonymous substitution";
				}
				$variant{'position'} = $transcript_mut_pos;
				$variant{'coverage'} = sprintf("%d/%d",$var_coverage,$total_coverage);
				# if (defined($var_comments) && defined($var_comments->{abs($mut_pos)}{$error_type})) {
					# $variant{'comments'} = join("; ",@{$var_comments->{abs($mut_pos)}{$error_type}});
				# }
				if (!defined($annotated_variants{$variant{'mutation'}})){
					push(@{$variants->{$transcript_id}}, \%variant);
					$annotated_variants{$variant{'mutation'}} = 1;
				}
# UNFINISHED:
			} elsif (defined($params->{'indels'}) && $params->{'indels'} && $error_type eq 'deletion'){
				# next;
				push(@{$deletion->{'reference'}},abs($mut_pos));
				push(@{$deletion->{'transcript'}},$transcript_mut_pos);
				# if (defined($var_comments) && defined($var_comments->{abs($mut_pos)}{$error_type})) {
					# $deletion_comments = join("; ",@{$var_comments->{abs($mut_pos)}{$error_type}});
				# }
			} elsif (defined($params->{'indels'}) && $params->{'indels'} && $error_type eq 'insertion'){
				# next;
				push(@{$insertion->{'reference'}},abs($mut_pos));
				push(@{$insertion->{'transcript'}},$transcript_mut_pos);
				# if (defined($var_comments) && defined($var_comments->{abs($mut_pos)}{$error_type})) {
					# $insertion_comments = join("; ",@{$var_comments->{abs($mut_pos)}{$error_type}});
				# }
			}
		}
	}
	return $variants;
}

#################################################################################

# Retrieves alleles directly from a BAM/SAM file
sub retrieve_alleles_from_bam_data {

	my ($infile,$ref_sequences_file,$params) = @_;
	
	my $alleles = {};

	if (!defined($params)){
		$params = {};
	}

	# SAM data example (samtools view alignment.bam | less -S):
	# ERR351515.1879	256	HLA:HLA04427	270	9	35S210M249S	*	0	0	CTCGGGGGACCGGGCTGACCGCGGGGTCCGGGCCAGGTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCAGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGCCCAATTGTCTCCCCTCCTTGTGGGAGGCCAGCCCGGGAGATCTACAGGCGATCAGGGAGGCGCCCCGTGGCCCCTGGTACCCGTGCGCTGAGCGTCTCCTTCCCGTTCTCCAGGTATCTGCGGAGCCACTCCACGCACGTGCCCTCCAGGTAGGCTCTCAACTGCTCCGCCACCTGGGCCGCCTCCCACTTGTGCTTGGTGGTCTGAGCTGCCATGTCCCCCGCGGTCCAAAGCCCCCGGGCCCCTTT	?BBADDDEEEDEGGGGGGIIHHHHHHHHHHHHHHHGGGGGGGGGGGEGGGGGECEEGGGCEGGGGGGGGGGG>GGGGDEDGGGGGGGGGGGGGEGGGGGGGAAEECCCEGGGGCGGGGGCEGADGGGEGGCECCEGGCEEC>CEECCC?*9?CECC<G88CCCECC:8A>>2<DDDD??:C?CGEEGEE?8:?CEGC?C:??E:::8.)'.'28888*::***.'4'8EE*?C::8:C?C?C*:1?A??ABBBDDDDDDDDEEFFFFFHHHHHHHHHHHHHHHHEHHIHHIIHFHHCCHHHFFFDEH<CHHHHEHBEHHFFDF4BFFFD*3@EEEEEEEB@AABCBCEEFAE**8EFEEFCEEE?EE'8;8?EFFFFF?EEE88>**?EE*?E::A:AEEEEE*A:*::*:*:A)?A88*8:C?CA44<)8CA:CEEE**:*::8*88)1*:**1*11A****10*.'.''').40*****.)'.4'.4''*0*	AS:i:313	XS:i:313	XN:i:0	XM:i:16	XO:i:0	XG:i:0	NM:i:16	MD:Z:12A6T0A21C5G21G0G3G34C43T8G14T5C10C0G5T7	YT:Z:UU
	# ERR351515.31106	256	HLA:HLA05726	270	9	35S210M249S	*	0	0	CTTGGGGGACTGGGCTGACCGCGGGGTCCGGGCCAGGTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCGGCAGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGCCCAATTGTCTCCCCTCCTTGTGGGAGGCCAGCCCGGGAGATCTACAGGCGATCAGGGAGGCGCCCCGTGGCCCCTGGTACCCGTGCGCTGAGCGTCTCCTTCCCGTTCTCCAGGTATCTGCGGAGCCACTCCACGCACGTGCCCTCCAGGTAGGCTCTCAACTGCTCCGCCACATGGGCCGCCTCCCACTTGTGCTTGGTGGTCTGAGCTGCCATGTCTGCCGCGGTCCAAGGGCGCAGGTCCTCTTT	ABBBDDDDDDEEGGGGGGIIHHHHHHHHIHHHHHHIHGHHHHHHHHHHHHHHHHHGGHGGGGGGGHGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCEGGGGGGGGGGGGGGGDGGGGGGCGCEEEGGEEEGDGCCCCCCCCEGGGED?ACCECGEEC<2<?AD?GEGGEGEGGGGGGGGEEC??EECEECCCGC??888<8'22.<?C?*::*482<?CG111?*8?CC??:9:1?????BBBDDDDDDEEGGGGGGIHHIIHHIIHHHHHHHHHHIIHGHIGHIHHHHDFHHIHHHGEHHHEEECDFHGFFFFFFGDDDDD>@@ECEA@>C;CCECGGGGGEE8C?CGG:CGGE:?DGDDGGEGCGGCGCE2>DGEGGGCEGEG?:::8C*?*CCCC:CE*CE8CE8EC?CEEGC2<A)?EG*0*:0*1:??CG?CE)1*1::C:?C?*:*11*0:?.''4)48*:?)))'''420:???:C:	AS:i:317	XS:i:317	XN:i:0	XM:i:16	XO:i:0	XG:i:0	NM:i:16	MD:Z:12A6T0A21C5G21G0G3G34C43T8G14T5C10C0G5T7	YT:Z:UU

	# Opens the SAM/BAM file, it should be in the same order than the reference genome/sequences
	if (is_sam($infile)){
		open(ALIGN_FILE,$infile) || die "# $0 : # 'retrieve_alleles_from_bam_data' cannot read '$infile'.\n";
	} elsif (is_bam($infile)){
		open(ALIGN_FILE,"$SAMTOOLS view $infile |") || die "# $0 : # 'retrieve_alleles_from_bam_data' cannot read '$infile'.\n";
	} else {
		print "# 'retrieve_alleles_from_bam_data': Unrecognized BAM or SAM format file '$infile'.\n";
		exit;
	}

	# Opens the Reference FASTA files
	print ("\nParsing reference file '$ref_sequences_file'.\n");
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

	print ("\nParsing alignment file '$infile'.\n");
	my $ref_coverage;
	my $ref_obj;
	my $ref_id = '';
	my $ref_seq = '';
	my %unique_reads_mapped;
	my $total_reads_mapped = 0;
	while(<ALIGN_FILE>){
		# s/\012\015?|\015\012?//;
		my @line = split("\t");
		if (scalar @line < 10){
			next;
		} elsif ($line[2] eq '*'){ # Empty sbjct name - no alignment
			next;
		}
		my $seq_name = $line[0];
		my $ref_name = $line[2];
		my $cigar = $line[5];
		my $ref_pos = $line[3];
		# Go to next reference sequence if the alignments of the previous one have finished
		while ($ref_id !~ /^\Q$ref_name/){
			$ref_obj = $reference_data->next_seq;
			$ref_id = $ref_obj->primary_id();
		}
# if ($seq_name eq 'ERR351515.1879' && $ref_name  eq 'HLA:HLA04427'){
# print '';
# }
		# Detects mismatches (substitutions) and indels (gap openings) in SAM format additional fields
		my ($mismatches,$indels) = (0,0);
		if (/XM:i:(\d+).+XO:i:(\d+)/){
			$mismatches = $1;
			$indels = $2;
		}
# 		} else {
# 			# Annotates mismatches (substitutions) and indels (gap openings) when SAM format has not additional fields
# 			my $seq_seq = $line[9];
# 			my $seq_pos = 1;
# 			my $new_ins = 0;
# 			my $seq_len = 0;
# 			while ($cigar =~ /(\d+)([MIDNSHP=X]+)/g){
# 				# Ex. CIGAR: 3M1I3M1D5M
# 				# CCATACT-GAACTGACTAAC
# 				#     ACTAGAA-TGGCT
# 				foreach (1..$1){
# 					if ($2 eq 'M' || $2 eq '=' || $2 eq 'X'){
# 						my $ref_nt = $ref_obj->subseq($ref_pos,$ref_pos);
# 						my $seq_nt = substr($seq_seq,$seq_pos-1,1);
# 						if ($ref_nt ne $seq_nt){
# 							$mismatches++;
# 						}
# 						$new_ins = 0;
# 						$seq_pos++;
# 						$ref_pos++;
# 					} elsif ($2 eq 'I') {
# 						# Skips homopolymer sequencing errors
# 						if (substr($seq_seq,$seq_pos-4,6) !~ /[ACTG]{3}/){
# 							# Annotates the insertion coverage only once in the reference insertion position
# 							if (!$new_ins){
# 								$indels++;
# 								$new_ins = $ref_pos;
# 							}
# 						}
# 						$seq_pos++;
# 					} elsif ($2 eq 'D') {
# 						# Skips homopolymer sequencing errors
# 						if (substr($seq_seq,$seq_pos-4,6) !~ /[ACTG]{3}/){
# 							$indels++;
# 							$new_ins = 0;
# 						}
# 						$ref_pos++;
# 					} elsif ($2 eq 'S' || $2 eq 'H') {
# 						# Soft-clipping: bases in 5' and 3' of the read are NOT part of the alignment.
# 						$seq_pos++;
# 					}
# 				}
# 			}
# 		}
		# Continues annotation only if there are not mismatches or gap openings in the alignment
		if (defined($params->{'max_errors'}) && $mismatches+$indels > $params->{'max_errors'}){
			next;
		}
		# Detects mismatches in SAM format additional fields and discards wrong alignments
		if (defined($params->{'max_substitutions'}) && $mismatches > $params->{'max_substitutions'}){
			next;
		}
		# Detects indels in SAM format additional fields and discards wrong alignments
		if (defined($params->{'max_indels'}) && $indels > $params->{'max_indels'}){
			next;
		}
		# Annotates allele sequence and coverage
		$alleles->{$ref_id}{'coverage'}++;
		if (!defined($alleles->{$ref_id}{'sequence'})){
			$alleles->{$ref_id}{'sequence'} = $ref_obj->seq();
		}
		# Annotates mapped reads ids 
		push(@{$alleles->{$ref_id}{'mapped_seqs'}}, $seq_name);
		# Annotates mapped reads with errors
		$alleles->{$ref_id}{'errors'} += $mismatches+$indels;
		# Annotates the depth in every position of the alignment
		my $seq_pos = 1;
		while ($cigar =~ /(\d+)([MSH=X]+)/g){
			foreach (1..$1){
				if ($2 eq 'M' || $2 eq '=' || $2 eq 'X'){
					$ref_coverage->{$ref_id}{$ref_pos}++;
					$ref_pos++;
				} elsif ($2 eq 'S') {
					# Soft-clipping: bases in 5' and 3' of the read are NOT part of the alignment.
					$seq_pos++;
				} elsif ($2 eq 'H') {
					# Hard-clipping: bases in 5' and 3' of the read are NOT part of the alignment AND those bases have been removed from the read sequence in the BAM file. The 'real' sequence length would be length(SEQ)+ count-of-hard-clipped-bases
					$seq_pos++;
				}
			}
		}
# 		# Detects longest mapped length
# 		while ($cigar =~ /(\d+)[M=X]/g){
# 			if ($alleles->{$ref_id}{'mapped_length'} < $1) {
# 				$alleles->{$ref_id}{'mapped_length'} = $1;
# 			}
# 		}
		$total_reads_mapped++;
		$unique_reads_mapped{$seq_name} = 1;
		if (defined($params->{'limit'}) && $params->{'limit'}<=$total_reads_mapped){
			last;
		}
	}
	close(ALIGN_FILE);

	my $unique_reads_mapped = scalar keys %unique_reads_mapped;
	printf("\n%d reads mapped to %d reference alleles.\n",$unique_reads_mapped,scalar keys %$alleles);
	
	# Annotates mapped lengths and filters low coverage/frequency alleles
# 	foreach my $ref_id (sort {$alleles->{$b}{'coverage'}<=>$alleles->{$a}{'coverage'}} keys %$alleles){
	foreach my $ref_id (keys %$ref_coverage){
		# Annotates mapped lengths
		$alleles->{$ref_id}{'mapped_length'} = scalar keys %{$ref_coverage->{$ref_id}};
		# Annotates total number of mapped nts
		map $alleles->{$ref_id}{'mapped_area'} += $ref_coverage->{$ref_id}{$_} , keys %{$ref_coverage->{$ref_id}};
		# Excludes low frequency alleles (references)
		my $keep_allele = 1;
		# Alleles with lower frequency than threshold
		if (defined($params->{'min_allele_frequency'}) && $alleles->{$ref_id}{'coverage'} < $params->{'min_allele_frequency'}/100*$unique_reads_mapped){
			unless (defined($params->{'rare_variants'}) && $params->{'rare_variants'}){
				$keep_allele = 0;
			} else {
				push(@{$alleles->{$ref_id}{'comments'}}, sprintf('Low frequency allele (%.2f%%).',100*$alleles->{$ref_id}{'coverage'}/$unique_reads_mapped));
			}
		}
		# Alleles with lower coverage than threshold
		if (defined($params->{'min_allele_coverage'}) && $alleles->{$ref_id}{'coverage'} < $params->{'min_allele_coverage'}){
			unless (defined($params->{'rare_variants'}) && $params->{'rare_variants'}){
				$keep_allele = 0;
			} else {
				push(@{$alleles->{$ref_id}{'comments'}}, sprintf('Low coverage allele (%d reads).',$alleles->{$ref_id}{'coverage'}));
			}
		}
		if (!$keep_allele){
			delete($alleles->{$ref_id})
		}
	}

	printf("\n%d alleles are kept after filtering artifacts.\n", scalar keys %$alleles);

	return $alleles;

}


# #################################################################################
# 
# OBSOLETE:
# # Retrieves HLA alleles directly from a BAM/SAM file
# sub retrieve_hla_alleles_from_bam_data {
# 
# 	my ($infile,$ref_sequences_file,$params,$md5_to_depth) = @_;
# 	
# 	my $alleles = {};
# 
# 	if (!defined($params)){
# 		$params = {};
# 	}
# 
# 	# Reads references file
# 	my $references = read_fasta_file_hash($ref_sequences_file);
# 	# Parses reference names
# 	my %ref_name_to_allele;
# 	foreach my $ref_header (keys %$references){
# 		# Ex. >HLA:HLA00005 A*02:01:01:01 1098 bp
# 		my @header_fields = split('\s',$ref_header);
# 		# Annotates the classical HLA loci (and more)
# 		# if ($header_fields[1] =~ /(\w+\*\d+\:\d+)
# 		if ($header_fields[1] =~ /(A|B|C|DQA1|DQB1|DRB1)(\*\d+\:\d+)/){
# 			$ref_name_to_allele{$header_fields[0]} = $1.$2;
# 		} else {
# 			print '';
# 		}
# 	}
# 
# 	# SAM data example (samtools view alignment.bam | less -S):
# 	# ERR351515.1879	256	HLA:HLA04427	270	9	35S210M249S	*	0	0	CTCGGGGGACCGGGCTGACCGCGGGGTCCGGGCCAGGTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCAGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGCCCAATTGTCTCCCCTCCTTGTGGGAGGCCAGCCCGGGAGATCTACAGGCGATCAGGGAGGCGCCCCGTGGCCCCTGGTACCCGTGCGCTGAGCGTCTCCTTCCCGTTCTCCAGGTATCTGCGGAGCCACTCCACGCACGTGCCCTCCAGGTAGGCTCTCAACTGCTCCGCCACCTGGGCCGCCTCCCACTTGTGCTTGGTGGTCTGAGCTGCCATGTCCCCCGCGGTCCAAAGCCCCCGGGCCCCTTT	?BBADDDEEEDEGGGGGGIIHHHHHHHHHHHHHHHGGGGGGGGGGGEGGGGGECEEGGGCEGGGGGGGGGGG>GGGGDEDGGGGGGGGGGGGGEGGGGGGGAAEECCCEGGGGCGGGGGCEGADGGGEGGCECCEGGCEEC>CEECCC?*9?CECC<G88CCCECC:8A>>2<DDDD??:C?CGEEGEE?8:?CEGC?C:??E:::8.)'.'28888*::***.'4'8EE*?C::8:C?C?C*:1?A??ABBBDDDDDDDDEEFFFFFHHHHHHHHHHHHHHHHEHHIHHIIHFHHCCHHHFFFDEH<CHHHHEHBEHHFFDF4BFFFD*3@EEEEEEEB@AABCBCEEFAE**8EFEEFCEEE?EE'8;8?EFFFFF?EEE88>**?EE*?E::A:AEEEEE*A:*::*:*:A)?A88*8:C?CA44<)8CA:CEEE**:*::8*88)1*:**1*11A****10*.'.''').40*****.)'.4'.4''*0*	AS:i:313	XS:i:313	XN:i:0	XM:i:16	XO:i:0	XG:i:0	NM:i:16	MD:Z:12A6T0A21C5G21G0G3G34C43T8G14T5C10C0G5T7	YT:Z:UU
# 	# ERR351515.4697	256	HLA:HLA04427	270	9	35S210M249S	*	0	0	CTTGGGGGACTGGGCTGACCGCGGGGTCCGGGCCAGGTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCAGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGCCCAATTGTCTCCCCTCCTTGTGGGAGGCCAGCCCGGGAGATCTACAGGCGATCAGGGAGGCGCCCCGTGGCCCCTGGTACCCGTGCGCTGAGCGTCTCCTTCCCGTTCTCCAGGTATCTGCGGAGCCACTCCACGCACGTGCCCTCCAGGTAGGCTCTCAACTGCTCCGCCACATGGGCCGCCTCCCCCTTGTGCTTGGTGGTCTGAGCTGCCATGTCCGCTGCGGTCCAAGAGCGCAGGTCCTCTTT	ABBBDDDDDDDDGGGGGGIIHHHHHHHHIHHHHHHIHHGHHHHHFHHHHHHHHHHGGGGGGGFGGGGGGGGGBEGGGGGGGGG2CEGGGGGGGGGGGGGGGG>GGEEGGGCEEEDDEEEGEGGGGGGGEGECCCEGGCEGGDCC:C*C*8?CEGGE>AG?EEEECC:8AAA8>DD8>C?CCEGEEGEEGCC:88CE98:::?C*::**8'82222<CCC*::*542'4*:*00:?*??E:?C??:?????BBBDDDDDDDDGGGGGFFEHHIHHIIIIIHGHHHHHHHIIIIHHHHH>FHHHIHH7CCEHHHE+CCHHHGFGGFGGGGGGEDBGGGCEAB8EEGGEGGGGGG???EG:CE1?CGGEG28;>GGGGCEGGGCGGGGE8?CE*08:CC:?CEE?C*:CEGEE?:CE8CC'00?CEEE82<<8CCC)..0*1*1:?C80.?C?*:0*:**19*0*:?088)'..'048::?***'.'.2*1?*:*1*	AS:i:321	XS:i:321	XN:i:0	XM:i:16	XO:i:0	XG:i:0	NM:i:16	MD:Z:12A6T0A21C5G21G0G3G34C43T8G14T5C10C0G5T7	YT:Z:UU
# 	# ERR351515.7514	256	HLA:HLA04427	270	2	35S210M249S	*	0	0	CTTGGGGGACCGGGCTGACCGCGGGGTCCGGGCCAGGTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCAGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGCTGAGAGCCTACCTGCCCAATTGTCTCCCCTCCTTGTGGGAGGCCAGCCCGGGAGATCTACAGGCGATCAGGGAGGCGCCCCGTGGCCCCTGGTACCCGTGCGCTGAGCGTCTCCTTCCCGTTCTCCAGGTATCTGCGGAGCCACTCCACGCACGTGCCCTCCAGGTAGGCTCTCAACTGCTCCGCCACATGGGCCGCCTCCCACTTGTGCTTGGTGGTCTGAGCTGCCATGTCCGCTGCGGTCCCAGAGCGCAGGTCCTCTTT	ABBBDDEEDDDDGGGGGGIIHHHHHHHHIHHHHHHHHHHHHHHHHHHHGGGGGGGGGGGEGGGEGGGGGGGHBGGGGGGGGGGGGGGGGGGGGGGGGGG;AD>GGEEGECCCECGDGCEEGEGGGGAEGEEGECGGEEEEGGGECC?CC09?CGEEGGACEECCCCC2>D><<A>>DGECEEEEGGGEECCCC8?CCEE:?EE:?:8.82.2'45<CE::*?*'.4'4CC*1*:*?:CGE:?CE:AAAAABBBDDDDDDDBGGGGGFHHHIIHHIGHHIHDHHHHHHGGGHIEEFHHHIEFGFHH*C<CHHHAEHHEHHGFFGFDFGGGDGE@GBEGGCEEECEEEEEGEECEGCEGCEE*::CCEGGG>GGGGGEEG:CCEGD29CCEECEGEEEEEG8C:C:CCEGGCE?EC)CE'80:CC*?C2>A)CEG?CCC*?*:C:?E:CE):*?C*?**1CC**1:C8<)'8<8)54****0:4>28>9?C*00**	AS:i:320	XS:i:320	XN:i:0	XM:i:15	XO:i:0	XG:i:0	NM:i:15	MD:Z:12A6T0A21C5G21G0G3G34C43T8G14T5C11G5T7	YT:Z:UU
# 	# ERR351515.13525	272	HLA:HLA04427	381	2	3S155M1D10M326S	*	0	0	AAAGAGGACCTGCGCTCTTGGACCGCAGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTCAGCGCACGGGTACCAGGGGCCACGGGGCGCCTCCCTGATCGCCTGTAGATCTCCCGGGCTGGCCTCCCACAAGGAGGGGAGACAATTGGGCAGGTAGGCTCTCAACTGCTCCGCCACATGGGCCGCCTCCCACTTGTGCTTGGTGGTCTGAGCTGCCATGTCCGCTGCGGTCCAAGAGCGCAGGTCCTCTTTCAGGGCGATGTAATCCTTGCCGTCGTAGGCGTACTGGTGGTACCCGCGGAGGAAGCGCCAGTCCGACCCCACGTCGCAGCCATACATCCTCTGGACGGTGTGAGAACCTGGCCCGGACCCCGCGGTCAGCCCGGTCCCCCAAG	ECEECC:98<>A<:0*10CGGCADD<C>8??0*C:9*?::*CCGECC?C?CEEGCEC:0?*GGE8G><CCCCEC?C2C8):C?GEGGGGGEC:C?CEECGEEGEGGGGCA>8CCC*EGGGGGGGGGGGGEECCCGGEGECECGGGGGGGGGEEGDGGGGGE@GEDEGGFFGGGGGHHEEEEHHHHHEEHEGHHFDFHHHIHHHFFHDFHCHHHHHIHFIHHCFHHHHGFGGGGDEDDDDDB?BB?????1:*:?C98:.**?C:EC82'4.*1**CCC4'''.)').0*::*CC?CEGECC?CGGGEGGEGEECGGGGAD<2DDGG>EEGECC?CD<>?C?GEC8CC?:ECGAEEEEGGGGGGGGGGGGGDGGGGGEGGGGGGGGGGECGGAGG>GGEGGGGGGGGGGEGGGGGGE;GGGGEGGGGGGGGGGGEEEGGGGGGGGGGHHHFHHHHHHHHHHHHHHHIHHHHHHHHIIGGGGGGDDDDDDDDBBB?	AS:i:247	XS:i:247	XN:i:0	XM:i:11	XO:i:1	XG:i:1	NM:i:12	MD:Z:43T8G14T5C10C0G5T13C0G9C0G37^G10	YT:Z:UU
# 	# ERR351515.34281	256	HLA:HLA05157	1	2	27S218M249S	*	0	0	CGGGTCTCAGCCACTGCTCGTCCCCAGGCTCTCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCAAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGTCCGGAGTATTGGGACGGGGAGACACGGAAAGTGAAGGCCCACTCACAGACTCTCGGACCCGGAGACTGTGGGCGACCTGGCCCGTCCGTGGGGGATGAGAGGTCGTGACCTGCGCCCCGGGCCGGGGTCACTCACCGGCCTCGCTCTGGTTGTAGTAGCCGCGCAGGGTCCCCAGGTCCACTCGGTGAGTCTGTGAGTGGGCCTTCACTTTCCGTGTCTCCCCGTCCCAATACTCCGGACCCTCCTGCTCTATCCACGGCGCCCGCGGCTCCATCCCCTGGCTTGCGGCGTCGCTGTCGC	AAA=DDDDDDDDGGGGGGIHHHHHIIHIIIIIIIIIIHIIHHHHAEHIHHIIIIIIIIHHFHIIIIHHHHHHGHHBEGGEGGGGGGGGGGGGGGGGGGG?CCEEEEGGCGGGGAGGGGGGGGC?EGCEEGGC>EGGCGDGGGA>DG<AGGGEGGE9C?CCE?CCEG>A'8>88ADD<<DECEE:ECEGC80?>CEE>D8>C?:?::?CC4442'4???C882)'9*?*:::80))09?:*:?CE:?????@@@DDDDDDDDFFFFEFHHHHHHIIIIHHHHDEFHHHHCEHFDFGD,DEHEFHHHHHHHHFFEEEEEE@BE'?CEEFFEEE4?DEFEDFFEFFEEAA*:?ACEA?D?DDDDDCCEEFFEFEEEFF?C?EAE0AEFFFAAEE?EEEEFFFFFFEAAEFE:?EF:AAEEDAADEEAE:E?CEED8<?EAAACEFEEFFEEE:CEDD??D>.882228)*0:?*0*?8?*04.'.2.4888)))*0*	AS:i:387	XS:i:387	XN:i:0	XM:i:8	XO:i:0	XG:i:0	NM:i:8	MD:Z:4C65C58A36G15C0A12T0A20	YT:Z:UU
# 	# ERR351515.42072	0	HLA:HLA05157	1	9	27S218M249S	*	0	0	CGGGTCTCAGCCACTGCTCGTCCCCAGGCTCTCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCAAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGTCCGGAGTATTGGGACGGCGAGACACGGAAAGTGAAGGCCCACTCACAGACTCTCGGACCCGGAGACTGTGGGCGACCTGGCCCGTCCGTGGGGGATGAGAGGTCGTGACCTGCGCCCCGGGCCGGGGTCACTCACCGGCCTCGCTCTGGTTGTAGTAGCCGCGCAGGGTCCCCAGGTCCACTCGGTGAGTCTGTGAGTGGGCCTTCACTTTCCGTGTCTCCCCGTCCCAATACTCCGGCCCCTCCTGCCCTATCCACGGCGCCCGCGGCTCCATCCCCTGGCTTGCGGCGTCGCTGTCGA	??@@DDDDDDDDFFFFFFIHHHHHHHHHHDHHIIIIIIIIIIIIIIIIIFGHHHIIIIIHHFFHHGAHHHHHHHFHDEEEB@@BEEEDDDEEFFEEEEEEE:AEEFEECA?AEDD2;DEDD??CEEAE?><?>EDE8ADEDD>D888?><8*8A?8?8CC?8?AC882.2<<><<8'48E?EAC?AEAA))0'8:*'<?20:*0?*0?A2<'48'4ACAA'40'0::?*1**8?)080?:?:ACC?????@@@DDBDDDDDFCFCCFCHHHHIIIIHHHHHHDFDBEC<CD=FFDH,C=)4CDHHFHHHHFF==@EEB@@@.8?C:CEEE?;?DAA?>CA*:C?CEAC??C:AC;D;D;>;2:AAAEA?C:CACE:ACE.8AAAEFE?ACEECEFEEEEEFF*?*:C?:CAE::AAE>??;AE::AE0?:.4'.8>A88:*1*01*1:C?CA;>288242'.4'88**10*00)*0*8'*.''''.'')))00*	AS:i:385	XS:i:268	XN:i:0	XM:i:9	XO:i:0	XG:i:0	NM:i:9	MD:Z:4C65C58A36G15C0A0G11T0A20	YT:Z:UU
# 	# ERR351515.25819	256	HLA:HLA05332	1	9	27S218M249S	*	0	0	CAGGTCTCAGCCACTCCTCGTCCCCAGGCTCTCACTCCATGAGGTATTTCTTCACATCCGTGTCCAGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGTCCGGAGTATTGGGACGGGGAGACACGGAAAGTGAAGGCCCACTCACAGACTCTCGGACCCGGAGACTGTGGGCGACCTGGCCCGTCCGTGGGGGATGAGAGGTCGTGACCTGCGCCCCGGGCCGGGGTCACTCACCGGCCTCGCTCTGGTTGTAGTAGCCGCGCAGGGTCCCCAGGTCCACTCGGTGAGTCTGTGAGTGGGCCTTCACTTTCCGTGTCTCCCCGTCCCAATACTCCGGACCCTCCTGCTCTATCCACGGCGCCCGCGGCTCCATCCTCTGGCTCGCGGCGTCGCTGTCGA	?BBBBBDDDDDDFGGGGGIHIHHHIIHIIIIIHHIHIIIIIIIIEFHIIIIIIIIIIHIHHFHHHHHHHIHHHHHHGHHHGGGEGGGGGGGGGGGGGGGEGCGGGEGGEEEGEGGDGGGGGGCEGGGG?ADCDGGGEGGGGGG>GD>DDDD8<CE9CCEGECCEEC8824<''8>82<DECCECEGGGGC8?>CECA<E>E*::::CCCD'''42?8:?E28).?:C:C::C*).89?*::?CEC?????BAAD@DBBDDDGFFFEGHHHHHIIGIIHHHHHFFHHHHCEHGFFGHCFHCHHHHHHHHHHGGGGGGGGGGG4CCEGGGGEEDGGGEEGEGEGGEEGEGEEC?EED>;;;>GG?CEEGGCEGEGGGGECE8CCCCCCGEGEGEGEEGCCEGGGEEEGEGGGCE?EGEGGCGGGGGGCEGCEGG?8>ECEEC0?*:?C:CECCE222><8488>G'8C:CEE*:*?:::*'.8''.'..8).8*)*	AS:i:384	XS:i:384	XN:i:0	XM:i:8	XO:i:0	XG:i:0	NM:i:8	MD:Z:4C65C58A36G15C0A12T0A20	YT:Z:UU
# 	# ERR351515.16126	16	HLA:HLA05546	381	2	3S155M1D10M326S	*	0	0	AAAGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGGTACCTGGAGAACGGAAAGGAGACGCTCAGCGCACGGGTACCAGGGGCCACGGGGCGCCTCCCTGATCGCCTGTAGATCTCCCGGGCTGGCCTCCCACAAGGAGGGGAGACAATTGGGCAGGTAGGCTCTCAACTGCTCCGCCACATGGGCCGCCTCCCACTTGTGCTTGGTGGTCTGAGCTGCCATGTCCGCCGCGGTCCAAGAGCGCAGGTCCTCTTTCAGGGCGATGTAATCCTTGCCGTCGTAGGCGTACTGGTGGTACCCGCGGAGGAAGCGCCAGTCCGACCCCACGTCGCAGCCATACATCCTCTGGACGGTGTGAGAACCTGGCCCGGACCCCGCGGTCAGCCCGGTCCCCCAAG	EECC??:*A82><:09CCGD>8>4.<GA?C?:?CGEC:*CEGGCCEC?:CEC?EC:*C?9*GGC8><GEECGGE??2C8.?C:ECEE?EC:CEGGE?:EC*GGGGGECCGGGGGECEECGGGDGGGGGGGGECCGGEGGGGGGEGGGGGGCCGEGEGEGGGGGGGGGGFFFGGGGHHHHHEHHHHDHHHHFIIHHHHHHIIIIHHIIIHHBHHHHIHIIIHHIHHHIGGGGGGEDEEEEDDBBBAA???:ECCC:C?:8:*E?:EE84'44**::?CC..''4)80*C:::CCC:EGEEC88:CGGEGECGGGCEGGG?A?<42?88EEEEECE?<DDGEGGGECECCEGCGDGCEEGGGEGGGEGEGGDDGGGGGEGGGEGGGGGGGGGGDGGDGGGGGGGGGGGGGGEGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGHHHHHHHHHHHHHHHHHHIHHHHHHHHIIGGGGGGDDEEEEDDBABA	AS:i:241	XS:i:239	XN:i:0	XM:i:12	XO:i:1	XG:i:1	NM:i:13	MD:Z:43T8G14T5C10C0G5T13C0G9C0G10A26^G10	YT:Z:UU
# 	# ERR351515.16904	272	HLA:HLA05546	381	2	3S155M1D10M326S	*	0	0	AAAGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGTGGAGTGGCTCCGCAGATACCTGGAGAACGGAAAGGAGACGCTCAGCGCACGGGTACCAGGGGCCACGGGGCGCCTCCCTGATCGCCTGTAGATCTCCCGGGCTGGCCTCCCACAAGGAGGGGAGACAATTGGGCAGGTAGGCTCTCAACTGCTCCGCCACATGGGCCGCCTCCCACTTGTGCTTGGTGGTCTGAGCTGCCATGTCCGCCGCGGTCCAAGAGCGCAGGTCCTCTTTCAGGGCGATGTAATCCTTGCCGTCGTAGGCGTACTGGTGGTACCCGCGGAGGAAGCGCCAGTCCGACCCCACGTCGCAGCCATACATCCTCTGGACGGTGTGAGAACCTGGCCCGGACCCCGCGGTCAGCCCGGTCCCCCAAG	?C:0*CC?<8.''*?0*E?84)44.2C>8EE?::GEE:*:**0*?.EC*GC:E?C0?:CC?GGGC>D>EECEGGEGDGGECE?GECEEC:**:E?C:1EGGGECEGGGED><C?EGGGGGGGDAAGDCCCGEC?EEGGGEGGGGGGEGGEE?C8ECC;EEGED:DGGGGFFGFFGHHHHHHHHHHHEHDEHHFFHIHHHHHHHHHIIIHHHHHEHIIIIIHHHHHHHGGGGGGDDDDDDDDBBBA????:CC:??C?C?::?*1:?<828.**:*CEC>84'5)88*?::CC:C:?C?C?88?CCEEGGCGGEECCEDAD><2>><>CECCCCCC8DA?ECEC?C?:CC??GDCCEEGGECEEGEGGGGDAGEGGEGEGDGGGEGGECEEGGGGGGGGGGGGGGGGGGGEGGGGGGGGGGGGHGHGGGGGGGGEDEGGGGGGGGGFHHHHHHHHHHHHHHHHHHHIHHHHHHHHIIGGGGGGDDDDDDDDBBBA	AS:i:247	XS:i:247	XN:i:0	XM:i:11	XO:i:1	XG:i:1	NM:i:12	MD:Z:43T8G14T5C10C0G5T13C0G9C0G37^G10	YT:Z:UU
# 	# ERR351515.32843	16	HLA:HLA05546	381	2	3S155M1D10M326S	*	0	0	AAAGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGACCACCAAGCGCAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGGAGGGCACGTGCGGGGAGTGGCTCCGCAGATACCTGGAGAACGGAAAGGAGACGCTCAGCGCACGGGTACCAGGGGCCACGGGGTGCCTCCCTGATCGCCTGTAGATCACCCGGGCTGGCCTCCCACAAGGAGGGGAGACAATTGGGCAGGTAGGCTCTCAACTGCTCCGCCACATGGGCCGCCTCCCACTTGCGCTTGGTGGTCTGAGCTGCCATGTCCGCCGCGGTCCAAGAGCGCAGGTCCTCTTTCAGGGCGATGTAATCCTTGCCGTCGTAGGCGTACTGGTGGTACCCGCGGAGGAAGCGCCAGTCCGACCCCACGTCGCAGCCATACATCCTCTGGACGGTGTGAGAACCTGGCCCGGACCCCGCGGTCAGCCCGGTCCCCCAAG	66?<E;/(48'4'69(<E=;8/'8';<;6;;;<?6EE;EEE?<<6/66(??9;<;EEEEEEE9E?42--<E<<(-=;<9EE;(?<;66((.6EEEE<6E<;<6<<6<E<44262-(=E<;6<;44;;E=E<;(/*;EEE:892*@EEE@@@8:@=91358))9+DCCEEDD@DC)@DCDEEC<<A5++*++,--+CCACEBA--.EA8+5++++7A7++6++++++CC>CC@A@<9@@@<@<><5,55,110::C80***0A?0EE;8'4.**:*:??.''2.)'8**::0?8>88??8.0**.*AEEEA?ECCEAEDD;D;82;2>CEFCAAAA>;88ECA?C??0**0:E;A8EEEEFEFEFEEEDDDD?A8EEECDDAEEFEEAC?:*;DD;DEECFEEEEDEDEFFDD?DDEDD>EEEEEEBFFFFEEEEFEDDFEDFEFEFHHHHHHHHFFFHHHHHEHDHHHHHHHHHHHFEFFDDDDDDDDDDBBB?	AS:i:259	XS:i:240	XN:i:0	XM:i:11	XO:i:1	XG:i:1	NM:i:12	MD:Z:43T23T5C10C0G5T13C0G5T3C0G37^G10	YT:Z:UU
# 	# ERR351515.11293	272	HLA:HLA05725	65	9	62M1D5M1I64M1I2M1D73M286S	*	0	0	CATCGCCGTGGGCTACGTGGACGACCCGCAGTTCGTGCGTTTCGACAGCGACGCCGCGAGTCCAAGAGGGGAGCCGCGGGAGCCGTGGGTGGAGCAGGAGGGGCCGGAGTATTGGGACCGGGAGACACAGAAGTACAAGCGCCAGGCACAGGCTGACCGAGTGAACCTGCGGAAACTGCGCGGCTACTACAACCAGAGCGAGGACGGTGAGTGACCCCGGCCCGGGGCGCAGGTCACGACTCCCCATCCTTGGACTCGCGGCGTCGCTGTCGAACCGCACGAACTGCGTGTCGTCCACGTAGCCCACTGCGATGAAGCGGGGCTCCCCGCGGCCGGGCCAGGACACGGATGTGGAGAAATACCTCATGGAGTGGGAGCCTGGGGACGAGGAGTGGCTGAGACCCGCCCGACCCTCCTCCCTGCGCGGCTCCCCGGGTCCTGCGCCCTCGCCGGGCGGGCCCCTCGCTCCTCTCCCCAGAGGCCATTTCCCTGCC	0''2.''0*0:8*0.)*88)..'.''EEA.;8)28;88.).*.1)?A8'>848';;81*:*C?*:E:EFEA2DD?>8EA?'D;DACEFFFFFFEAAAA*DD>>>ECFFFFFFEE>>8?DFEEACE?EEFFFFEA8EEDE>EEAA?A*C?CCCCEEEBCEEECCAFEEEBEEEDFEEEEFHHHHFF?HHHHF?HHHHFCHIHHHHHFFHIHHHHHHHBHHHHHHHHHHD@FFFFDDDDDDDD????????*:EE>8'4''.'2<<2CAC?E>?..28888C??8GD<C<8..':?8:9GECCCC::8CEEEGECC?'<.00?C882''<4'<<8ACCECEEC*EC?C:CEC??C?CEC:ECCGEEGGGGGEECECECECGEGEECGGGGGEEECGGGGGGDGDD;GGGDECGGGEDC8GGGGGGDGGGGGGGGGGGGGGGGGGGHHHHHHHHHEHHHHHHHHHHHIHIHHHHIHIIHGGGGGGEDDDDDDDBB@?	AS:i:257	XS:i:257	XN:i:0	XM:i:19	XO:i:4	XG:i:4	NM:i:23	MD:Z:25A34C1^A0G4T11C7A1A28A8G5^T0G3G0C3C0T5A9C7G3C0C33	YT:Z:UU
# 	# ERR351515.31106	256	HLA:HLA05726	270	9	35S210M249S	*	0	0	CTTGGGGGACTGGGCTGACCGCGGGGTCCGGGCCAGGTTCTCACACCGTCCAGAGGATGTATGGCTGCGACGTGGGGTCGGACTGGCGCTTCCTCCGCGGGTACCACCAGTACGCCTACGACGGCAAGGATTACATCGCCCTGAAAGAGGACCTGCGCTCTTGGACCGCGGCAGACATGGCAGCTCAGACCACCAAGCACAAGTGGGAGGCGGCCCATGTGGCGGAGCAGTTGAGAGCCTACCTGCCCAATTGTCTCCCCTCCTTGTGGGAGGCCAGCCCGGGAGATCTACAGGCGATCAGGGAGGCGCCCCGTGGCCCCTGGTACCCGTGCGCTGAGCGTCTCCTTCCCGTTCTCCAGGTATCTGCGGAGCCACTCCACGCACGTGCCCTCCAGGTAGGCTCTCAACTGCTCCGCCACATGGGCCGCCTCCCACTTGTGCTTGGTGGTCTGAGCTGCCATGTCTGCCGCGGTCCAAGGGCGCAGGTCCTCTTT	ABBBDDDDDDEEGGGGGGIIHHHHHHHHIHHHHHHIHGHHHHHHHHHHHHHHHHHGGHGGGGGGGHGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGCEGGGGGGGGGGGGGGGDGGGGGGCGCEEEGGEEEGDGCCCCCCCCEGGGED?ACCECGEEC<2<?AD?GEGGEGEGGGGGGGGEEC??EECEECCCGC??888<8'22.<?C?*::*482<?CG111?*8?CC??:9:1?????BBBDDDDDDEEGGGGGGIHHIIHHIIHHHHHHHHHHIIHGHIGHIHHHHDFHHIHHHGEHHHEEECDFHGFFFFFFGDDDDD>@@ECEA@>C;CCECGGGGGEE8C?CGG:CGGE:?DGDDGGEGCGGCGCE2>DGEGGGCEGEG?:::8C*?*CCCC:CE*CE8CE8EC?CEEGC2<A)?EG*0*:0*1:??CG?CE)1*1::C:?C?*:*11*0:?.''4)48*:?)))'''420:???:C:	AS:i:317	XS:i:317	XN:i:0	XM:i:16	XO:i:0	XG:i:0	NM:i:16	MD:Z:12A6T0A21C5G21G0G3G34C43T8G14T5C10C0G5T7	YT:Z:UU
# 
# 	# Opens the SAM/BAM file, it should be in the same order than the reference genome/sequences
# 	if (is_sam($infile)){
# 		open(ALIGN_FILE,$infile) || die "# $0 : # 'retrieve_alleles_from_bam_data' cannot read '$infile'.\n";
# 	} elsif (is_bam($infile)){
# 		open(ALIGN_FILE,"$SAMTOOLS view $infile |") || die "# $0 : # 'retrieve_alleles_from_bam_data' cannot read '$infile'.\n";
# 	} else {
# 		print "# 'retrieve_alleles_from_bam_data': Unrecognized BAM or SAM format file '$infile'.\n";
# 		exit;
# 	}
# 
# 	print ("\nParsing alignment file '$infile'.\n");
# 	my $total_coverage = 0;
# 	while(<ALIGN_FILE>){
# 		# s/\012\015?|\015\012?//;
# 		my @line = split("\t");
# 		if (scalar @line < 10){
# 			next;
# 		} elsif ($line[2] eq '*'){ # Empty sbjct name - no alignment
# 			next;
# 		}
# 		my $seq_name = $line[0];
# 		my $ref_name = $line[2];
# 		my $cigar = $line[5];
# 
# 		# Skips loci without standard annotated alleles
# 		if (!defined($ref_name_to_allele{$ref_name})){
# 			next;
# 		}
# 		my $reference_allele = $ref_name_to_allele{$ref_name};
# 
# 		# Detects maximum number of matches
# 		while ($cigar =~ /(\d+)[M=X]/g){
# 			if ($alleles->{$reference_allele}{'mapped_length'} < $1) {
# 				$alleles->{$reference_allele}{'mapped_length'} = $1;
# 			}
# 		}
# 		my $max_mapped_length = 0;
# 		# Detects maximum number of matches
# 		while ($cigar =~ /(\d+)[M=X]/g){
# 			if ($max_mapped_length < $1) {
# 				$max_mapped_length = $1;
# 			}
# 		}
# # 		# Skips reads previously annotated
# # 		if (defined($allele_mapped_reads->{$reference_allele}{$seq_name})){
# # 			next;
# # 		}
# # 		$allele_mapped_reads->{$reference_allele}{$seq_name} = 1;
# 
# 		# Detects mismatches in SAM format additional fields and discards wrong alignments
# 		if (defined($params->{'max_substitutions'}) && /XM:i:(\d+)/ && $1 > $params->{'max_substitutions'}){
# 			next;
# 		}
# 		# Detects indels in SAM format additional fields and discards wrong alignments
# 		if (defined($params->{'max_indels'}) && /XO:i:(\d+)/ && $1 > $params->{'max_indels'}){
# 			next;
# 		}
# 		$alleles->{$reference_allele}{'sequence'} = $references->{$reference_allele};
# 		$alleles->{$reference_allele}{'coverage'}++;
# 		$total_coverage++;
# 		if (defined($params->{'limit'}) && $params->{'limit'}<=$total_coverage){
# 			last;
# 		}
# 	}
# 	close(ALIGN_FILE);
# 
# 	printf("\n%d reads mapped to %d reference alleles.\n",$total_coverage,scalar keys %$alleles);
# 	
# 	# Filters low coverage/frequency alleles
# 	my $var_comments;
# 	foreach my $reference_allele (sort {$alleles->{$b}{'coverage'}<=>$alleles->{$a}{'coverage'}} keys %$alleles){
# printf("%s\t%d\t%d\n", $reference_allele,  $alleles->{$reference_allele}{'mapped_length'}, $alleles->{$reference_allele}{'coverage'});
# 		# Excludes low frequency alleles (references)
# 		my $keep_allele = 1;
# 		# Alleles with lower frequency than threshold
# 		if (defined($params->{'min_allele_frequency'}) && $alleles->{$reference_allele}{'coverage'} < $params->{'min_allele_frequency'}/100*$total_coverage){
# 			unless (defined($params->{'rare_variants'}) && $params->{'rare_variants'}){
# 				$keep_allele = 0;
# 			} else {
# 				push(@{$var_comments->{$reference_allele}}, sprintf('Low frequency allele (%.2f%%).',100*$alleles->{$reference_allele}{'coverage'}/$total_coverage));
# 			}
# 		}
# 		# Alleles with lower coverage than threshold
# 		if (defined($params->{'min_allele_coverage'}) && $alleles->{$reference_allele}{'coverage'} < $params->{'min_allele_coverage'}){
# 			unless (defined($params->{'rare_variants'}) && $params->{'rare_variants'}){
# 				$keep_allele = 0;
# 			} else {
# 				push(@{$var_comments->{$reference_allele}}, sprintf('Low coverage allele (%d reads).',$alleles->{$reference_allele}{'coverage'}));
# 			}
# 		}
# 		if (!$keep_allele){
# 			delete($alleles->{$reference_allele}{'coverage'})
# 		} else {
# 			print '';
# 		}
# 	}
# 
# 	printf("\n%d alleles are kept after filtering artifacts.\n", scalar keys %$alleles);
# 	
# exit;
# 	return $alleles;
# 
# }

#################################################################################

# Prints variant mutation information
sub print_variants {

	my ($variants,$mutation_classes,$params,$options) = @_;

	my $output;

	# Sets options
	my ($format,$detailed)=('text',0);
	foreach my $option (@$options){
		if ($option =~ /(text|html)\s?(short|long)?/){
			$format = $1;
			if ($2 eq 'long'){
				$detailed = 1;
			}
		}
	}

	# Reads Uniprot ID mapping information file (cross-references to external databases)
	my $ensembl_to_uniprot = read_uniprot_idmapping_file($params->{'idmapping_file'},'Ensembl_TRS',-1);
	# my $uniprot_to_ensembl_gene = read_uniprot_idmapping_file($params->{'idmapping_file'},'Ensembl');
	# my $uniprot_to_genename = read_uniprot_idmapping_file($params->{'idmapping_file'},'Gene_Name');

	# Retrieves protein information from UniProt
	my ($genename_to_uniprot, $genename_to_names, %ensembl_to_uniprot, $uniprot_to_ensembl);
	my (%cosmic_mutation_genes, %somatic_mutation_genes, %germline_mutation_genes, %unknown_mutation_genes);
	my $uniprot_data;
	foreach my $ensembl_id (keys %$variants){
		if ($ensembl_id !~ /^ENST\d+$/){ # Ensembl ID without version
			printf("\tWARNING: Ensembl header format not recognized: %s.\n", $ensembl_id);
			next;
		}
# 		# Annotates Ensembl reference data
# 		# ENST00000390588.1 cds chromosome:GRCh38:14:105912257:105912276:-1 gene:ENSG00000211928.1 gene_biotype:IG_D_gene transcript_biotype:IG_D_gene gene_symbol:IGHD5-5 description:immunoglobulin heavy diversity 5-5 [Source:HGNC Symbol;Acc:HGNC:5511]
# 		my $ensembl_id;
# 		if ($ensembl_full_id =~ /^(ENST\d+)\.\d+/){
# 			$ensembl_id = $1;
# 			# Renames the $variants hash key with the $ensembl_id without version number
# 			push(@{$variants->{$ensembl_id}}, @{$variants->{$ensembl_full_id}});
# 			delete($variants->{$ensembl_full_id});
# 		} elsif ($ensembl_full_id =~ /^(ENST\d+)/){ # Ensembl ID without version
# 			$ensembl_id = $1;
# 		} else {
# 			printf("\nWARNING: Ensembl header format not recognized: %s.\n", $ensembl_full_id);
# 			next;
# 		}
		# Converts the Ensembl ID into Uniprot one
		if (!defined($ensembl_to_uniprot->{$ensembl_id})){
			print "\tWARNING: '$ensembl_id' has not UniProt ID associated.\n";
			# $ensembl_to_genename{$ensembl_id} = "ZZZ$ensembl_id";
			next;
		} elsif ($#{$ensembl_to_uniprot->{$ensembl_id}} > 0){
			print "\tWARNING: '$ensembl_id' has multiple UniProt ID associated.\n";
		}
		# Each Ensembl transcript ID will be linked to a unique UniProt protein ID without isoform number
		my $uniprot_id = $ensembl_to_uniprot->{$ensembl_id}[0];
		if ($uniprot_id =~ /(\w+?)\-\d+/){
			$uniprot_id = $1;
		}
		$ensembl_to_uniprot{$ensembl_id} = $uniprot_id;
		push(@{$uniprot_to_ensembl->{$uniprot_id}},$ensembl_id);
		# Retrieves protein data from UniProt (if not already present in the 'oncodata' folder)
		if (!defined($uniprot_data->{$uniprot_id})){
			if (!-d $params->{'oncodata_folder'}){
				mkdir($params->{'oncodata_folder'});
			}
			my $uniprot_data_file = sprintf("%s/%s.txt", $params->{'oncodata_folder'}, $uniprot_id);
			if (!-e $uniprot_data_file){
				my $uniprot_txt_data = retrieve_uniprot_data($uniprot_id,'txt');
				if (defined($uniprot_txt_data)){
					write_to_file($uniprot_data_file,$uniprot_txt_data);
				}
			}
			$uniprot_data->{$uniprot_id} = read_uniprot_single_file($uniprot_data_file);
			# If file is erroneous (empty for ex.) tries to retrieve the data once more:
			if (!defined($uniprot_data->{$uniprot_id})){
				my $uniprot_txt_data = retrieve_uniprot_data($uniprot_id,'txt');
				if (defined($uniprot_txt_data)){
					write_to_file($uniprot_data_file,$uniprot_txt_data);
					$uniprot_data->{$uniprot_id} = read_uniprot_single_file($uniprot_data_file);
				}
			}
		}
		# Annotates the gene names, only the standard part of the name (ex: TP53 {ECO:0000313|Ensembl:ENSP00000425104})
		my @gene_names = map { (split(/\s/,$_))[0] } @{$uniprot_data->{$uniprot_id}{'gene_names'}};
		my $gene_name = $gene_names[0];
		if (!defined($gene_name)){
			$gene_name = "ZZ$ensembl_id"; # To show in last positions
			@gene_names = ( $ensembl_id );
		}
		
		# Checks if the gene is in the list of important genes for the cancer type
# 		my $important = 0;
# 		foreach my $cancer_gene (@$cancer_genes){
# 			if (in_array($uniprot_data->{$uniprot_id}{'gene_names'},"/$cancer_gene/")){
# 				$important = 1;
# 				last;
# 			}
# 		}
		my ($cosmic_mutation, $somatic_mutation, $germline_mutation, $unknown_mutation) = (0,0,0,0);
		foreach my $variant (@{$variants->{$ensembl_id}}){
			my $mut_pos = $variant->{'genome_position'};
			if (defined($mutation_classes->{$mut_pos})){
				# Ex. 'COSM6224-somatic'
				if ($mutation_classes->{$mut_pos} =~ /^COSM/i){
					$cosmic_mutation = 1;
				}
				# Ex. 'germline'
				if ($mutation_classes->{$mut_pos} =~ /somatic/){
					$somatic_mutation = 1;
				} elsif ($mutation_classes->{$mut_pos} =~ /germline/){
					$germline_mutation = 1;
				} elsif ($mutation_classes->{$mut_pos} =~ /unknown/){
					$unknown_mutation = 1;
				}
			}
		}
		# COSMIC genes will be annotated in 2 groups
		if ($cosmic_mutation){
			$cosmic_mutation_genes{$gene_name} = 1;
		}
		if ($somatic_mutation) {
			$somatic_mutation_genes{$gene_name} = 1;
		} elsif ($germline_mutation) {
			$germline_mutation_genes{$gene_name} = 1;
		} else {
			$unknown_mutation_genes{$gene_name} = 1;
		}
		
		# Stores the transcripts and proteins associated with the same gene name
		# push(@{$genename_to_ensembl->{$gene_name}}, $ensembl_id); # Can be several transcripts for the same gene
		push(@{$genename_to_uniprot->{$gene_name}}, $uniprot_id); # Can be several proteins for the same gene
		push(@{$genename_to_names->{$gene_name}}, @gene_names); # Can be several names for the same gene

	}

	# Stores information to print
	my @data_headers = ('Name', 'Gene', 'Transcripts', 'Proteins', 'Position', 'Variant', 'Class', 'Coverage', 'Cross-references', 'Drugs', 'Gene function');
	if ($format eq 'html') {
		$output .= "<table width=\"80%\" class=\"bordered\">\n\t<tbody>\n";
		# my @table_headers = ('Gene name', 'Ensembl ID', 'Mutation info', 'Gene function', 'Associated diseases');
		$output .= "\t\t<tr>\n\t\t\t<th>".join("</th>\n\t\t\t<th>", @data_headers)."</th>\n\t\t</tr>\n";
	} else {
		$output .= "\nVARIANT REPORT:\n\n";
		if (!$detailed) {
			$output .= join("\t", @data_headers)."\n"
		}
	}

	# Sorts mutations by gene names
	my @cosmic_mutation_genes = sort {$a cmp $b} keys %cosmic_mutation_genes;
	my @somatic_mutation_genes = sort {$a cmp $b} keys %somatic_mutation_genes;
	my @germline_mutation_genes = sort {$a cmp $b} keys %germline_mutation_genes;
	my @unknown_mutation_genes = sort {$a cmp $b} keys %unknown_mutation_genes;

	# Loops all the variant mutations and prints the information
	my %annotated_genes;
	foreach my $gene_name (@cosmic_mutation_genes,@somatic_mutation_genes,@unknown_mutation_genes,@germline_mutation_genes){
		# Skips already annotated COSMIC genes that will be also in the other classes
		if (defined($annotated_genes{$gene_name})){
			next;
		}
		$annotated_genes{$gene_name} = 1;
		my @gene_names = unique(@{$genename_to_names->{$gene_name}});
		my @uniprot_ids = unique(@{$genename_to_uniprot->{$gene_name}});
		my @ensembl_gene_ids = ();
		my $ensembl_ids = [];
		my @position_info = ();
		my @variant_info = ();
		my @class_info = ();
		my @coverage_info = ();
		my @crossref_info = ();
		my @drugs_info = ();
		my $function_info = '';
		my @gene_diseases = ();
		foreach my $uniprot_id (@uniprot_ids){

			# Annotates the Ensembl transcript IDs, Ensembl Gene IDs and Uniprot IDs
			
			# Can be several transcripts for the same protein
			push(@$ensembl_ids, [@{$uniprot_to_ensembl->{$uniprot_id}}]);
			# All the transcripts of the same protein will have the same annotations
			my $ensembl_id = $uniprot_to_ensembl->{$uniprot_id}[0];
# if ($ensembl_id eq "ENST00000257430"){
# print '';
# }
			
			# Annotates variant mutations details
			my (@variant_info_, @position_info_, @class_info_, @coverage_info_);
			foreach my $variant (@{$variants->{$ensembl_id}}){
				if (defined($uniprot_data->{$uniprot_id}{'variants'}) && $variant->{'variant'} =~ /^([A-Z])(\d+)([A-Z\*])?/){
					my ($ref_aa, $ref_pos, $var_aa) = ($1,$2,$3);
					# Retrieves UniProt identical variant mutation annotations
					foreach my $variant_ (@{$uniprot_data->{$uniprot_id}{'variants'}}){
						if ($variant_ =~ /^($ref_aa)($ref_pos)(\Q$var_aa\E)(.+)/){ # $var_aa can be *
							if ($format eq 'html') {
								push(@crossref_info,"<b>".convert_aa_1_to_3($1)."$2".convert_aa_1_to_3($3)."</b>$4");
							} else {
								push(@crossref_info,convert_aa_1_to_3($1)."$2".convert_aa_1_to_3($3)."$4");
							}
							#@crossref_info = split (/;\s?/,$variant_);
							last;
						}
					}
					# Try to find other mutations in the same position
					if (!@crossref_info) {
						foreach my $variant_ (@{$uniprot_data->{$uniprot_id}{'variants'}}){
							if ($variant_ =~ /^([A-Z])($ref_pos)([A-Z\*])(.+)/){
								if ($format eq 'html') {
									push(@crossref_info,"<b>".convert_aa_1_to_3($1)."$2".convert_aa_1_to_3($3)."</b>$4");
								} else {
									push(@crossref_info,convert_aa_1_to_3($1)."$2".convert_aa_1_to_3($3)."$4");
								}
							}
						}
					}
					# Creates links to cross-reference databases
					foreach my $crossref_info (@crossref_info) {
						$crossref_info =~ s/;\s+?;/;/g;
						if ($format eq 'html') {
							# Creates link to dbSNP (ex. dbSNP:rs587780076;)
							$crossref_info =~ s/dbSNP:rs(\d+)/<a target="_blank" href="https:\/\/www.ncbi.nlm.nih.gov\/SNP\/snp_ref.cgi?type=rs&rs=$1">$1<\/a>/g;
							# Creates link to PubMed (ex. PubMed:15489334,17344846)
							$crossref_info =~ s/PubMed:([\d,]+)/<a target="_blank" href="https:\/\/www.ncbi.nlm.nih.gov\/pubmed\/?term=$1\[uid\]">$1<\/a>/g;
						}
					}
				}
				if (defined($variant->{'comments'})){
					push(@variant_info_, sprintf("%s (%s)",$variant->{'mutation'},$variant->{'comments'}));
				} else {
					push(@variant_info_, sprintf("%s",$variant->{'mutation'}));
				}
				if (defined($variant->{'genome_position'})){
					if ($format eq 'html'){
						push(@position_info_, sprintf('<a target="_blank" href="http://www.ensembl.org/Homo_sapiens/Location/View?r=%s">%s</a>', $variant->{'genome_position'}, $variant->{'genome_position'}));
					} else {
						push(@position_info_, $variant->{'genome_position'});
					}
					my $mut_class = $mutation_classes->{$variant->{'genome_position'}};
					if ($format eq 'html' && $mut_class =~ /^COSM(\d+)\-(.+)/i){
						push(@class_info_, sprintf('%s (<a target="_blank" href="http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=%d">%s</a>)',ucfirst($2),$1,'COSM'.$1)); # uc($mut_class)
					} elsif ($format eq 'html' && $mut_class eq 'somatic'){
						push(@class_info_, sprintf("<font color=\"green\">%s</font>",ucfirst($mut_class)));
					}elsif ($mut_class =~ /^COSM(\d+)\-(.+)/i){
						push(@class_info_, sprintf("%s (%s)",ucfirst($2),'COSM'.$1));
					} else {
						push(@class_info_, sprintf("%s",ucfirst($mut_class)));
					}
				} else {
					push(@position_info_, "");
					push(@class_info_, "");
				}
# if (!defined($variant->{'coverage'})){
# print '';
# }
				if ($format eq 'html') {
					push(@ensembl_gene_ids, sprintf('<a target="_blank" href="http://www.ensembl.org/homo_sapiens/Gene/Summary?g=%s&db=core">%s</a>',$variant->{'gene'},$variant->{'gene'}));
					if (eval($variant->{'coverage'})*100 >= 20){
						push(@coverage_info_, sprintf("<font color=\"green\">%s</font>", $variant->{'coverage'}));
					} elsif (eval($variant->{'coverage'})*100 >= 5){
						push(@coverage_info_, sprintf("<font color=\"orange\">%s</font>", $variant->{'coverage'}));
					} else {
						push(@coverage_info_, sprintf("<font color=\"red\">%s</font>", $variant->{'coverage'}));
					}
				} else {
					push(@ensembl_gene_ids, $variant->{'gene'});
					push(@coverage_info_, $variant->{'coverage'});
				}
			}
			if ($format eq 'html') {
				push(@variant_info, join('<br><br>', @variant_info_));
				push(@position_info, join('<br><br>', @position_info_));
				push(@class_info, join('<br><br>', @class_info_));
				push(@coverage_info, join('<br><br>', @coverage_info_));
			} else {
				push(@variant_info, join('; ', @variant_info_));
				push(@position_info, join('; ', @position_info_));
				push(@class_info, join('; ', @class_info_));
				push(@coverage_info, join('; ', @coverage_info_));
			}
			# Creates link to KEGG Pathway entry
			if ($format eq 'html' && defined($uniprot_data->{$uniprot_id}{'KEGG'})){
				push(@crossref_info,sprintf('<a target="_blank" href="http://www.genome.jp/dbget-bin/www_bget?%s">KEGG Pathway</a>', $uniprot_data->{$uniprot_id}{'KEGG'}));
			}
			# Annotates the available drugs
			if (defined($uniprot_data->{$uniprot_id}{'DRUGBANK'})){
				foreach my $drug (@{$uniprot_data->{$uniprot_id}{'DRUGBANK'}}){
					if ($drug =~ /DB\d+; (.+)/ && $1 !~ /^(\d|[NOS]-)/){
						push(@drugs_info, $drug);
					}
				}
			}
			if ($format eq 'html') {
				map $drugs_info[$_] =~ s/(DB\d+); (.+)/<a target="_blank" href="https:\/\/www.drugbank.ca\/drugs\/$1">$2<\/a>/ , 0..$#drugs_info;
			}
			
			# Annotates detailed gene function
			if (defined($uniprot_data->{$uniprot_id}{'function'})){
				my $function_info_ = $uniprot_data->{$uniprot_id}{'function'}[0];
				if ($format eq 'html') {
					$function_info_ =~ s/PubMed:([\d,]+)/<a target="_blank" href="https:\/\/www.ncbi.nlm.nih.gov\/pubmed\/?term=$1\[uid\]">$&<\/a>/g;
					$function_info_ =~ s/\.\s+/.<br>/g;
				}
				if (length($function_info_) > length ($function_info)){
					$function_info = $function_info_;
				}
			}
			
# 			# Annotates cancer related diseases
# 			if (defined($uniprot_data->{$uniprot_id}{'disease'})){
# 				foreach my $disease (@{$uniprot_data->{$uniprot_id}{'disease'}}){
# 					if ($disease =~ /tumor|tumour|cancer|carcinoma|melanoma|leukemia|neoplasm|angiogenesis/i){
# 						push(@gene_diseases, $disease);
# 					}
# 				}
# 			}

		}
		my $rowspan = scalar @uniprot_ids;
		my $highlight = '';
		if (defined($cosmic_mutation_genes{$gene_name})){
			$highlight = ' style="background: #fbe9ec;"';
		}
		if ($format eq 'html') {
			$output .= "\t\t<tr$highlight>\n";
			$output .= sprintf("\t\t\t<td rowspan=\"%d\"><b>%s</b></td>\n", $rowspan, join("<br>",unique(@gene_names)));
			$output .= sprintf("\t\t\t<td rowspan=\"%d\">%s</td>\n", $rowspan, join("<br>",unique(@ensembl_gene_ids)));
			$output .= sprintf("\t\t\t<td>%s</td>\n", join("<br>",map sprintf('<a target="_blank" href="http://www.ensembl.org/homo_sapiens/Transcript/Summary?t=%s&db=core">%s</a>',$_,$_), @{$ensembl_ids->[0]}));
			$output .= sprintf("\t\t\t<td><a target=\"_blank\"href=\"http://www.uniprot.org/uniprot/%s\">%s</a></td>\n", $uniprot_ids[0], $uniprot_ids[0]);
			$output .= sprintf("\t\t\t<td><b>%s</b></td>\n", $position_info[0]);
			$output .= sprintf("\t\t\t<td><b>%s</b></td>\n", $variant_info[0]);
			$output .= sprintf("\t\t\t<td><b>%s</b></td>\n", $class_info[0]);
			$output .= sprintf("\t\t\t<td><b>%s</b></td>\n", $coverage_info[0]);
			$output .= sprintf("\t\t\t<td width=\"100\" rowspan=\"%d\">%s<br></td>\n", $rowspan, join("<br><br>",unique(@crossref_info)));
			$output .= sprintf("\t\t\t<td rowspan=\"%d\">%s<br></td>\n", $rowspan, join("<br>",unique(@drugs_info)));
			$output .= sprintf("\t\t\t<td rowspan=\"%d\"><font size=-1>%s</font></td>\n", $rowspan, $function_info);
			$output .= "\t\t</tr>\n";
			for (my $i=1; $i<=$#uniprot_ids; $i++){
				$output .= "\t\t<tr$highlight>\n";
				$output .= sprintf("\t\t\t<td>%s</td>\n", join("<br>",map sprintf('<a target="_blank" href="http://www.ensembl.org/homo_sapiens/Transcript/Summary?t=%s&db=core">%s</a>',$_,$_), @{$ensembl_ids->[$i]}));
				$output .= sprintf("\t\t\t<td><a target=\"_blank\" href=\"http://www.uniprot.org/uniprot/%s\">%s</a></td>\n", $uniprot_ids[$i], $uniprot_ids[$i]);
				$output .= sprintf("\t\t\t<td><b>%s</b></td>\n", join("<br><br>",$position_info[$i]));
				$output .= sprintf("\t\t\t<td><b>%s</b></td>\n", join("<br><br>",$variant_info[$i]));
				$output .= sprintf("\t\t\t<td><b>%s</b></td>\n", join("<br><br>",$class_info[$i]));
				$output .= sprintf("\t\t\t<td><b>%s</b></td>\n", join("<br><br>",$coverage_info[$i]));
				$output .= "\t\t</tr>\n";
			}
		} else {
			#$output .= "\n";
			if ($detailed) {
				$output .= sprintf("%s - Gene: %s:\n", join(",",unique(@gene_names)), join(",",unique(@ensembl_gene_ids)));
			}
			for (my $i=0; $i<=$#uniprot_ids; $i++){
				if ($detailed) {
					$output .= sprintf("\n\tPosition: %s\n\n", $position_info[$i]);
					$output .= sprintf("\n\tVariants: %s\n\n", $variant_info[$i]);
					$output .= sprintf("\tProtein: %s, Transcripts: %s\n", $uniprot_ids[$i], join(',',@{$ensembl_ids->[$i]}));
					$output .= sprintf("\tClass: %s\n", $class_info[$i]);
					$output .= sprintf("\tCoverage: %s\n", $coverage_info[$i]);
					$output .= sprintf("\tCross-references: %s\n", join("; ",unique(@crossref_info)));
					$output .= sprintf("\tProtein function: %s\n", (split(/\. /,$function_info))[0]);
					# $output .= sprintf("\t\tAssociated diseases: %s\n", $diseases_info);
				} else {
					$output .= sprintf("%s\t", join(",",unique(@gene_names)));
					$output .= sprintf("%s\t", join(",",unique(@ensembl_gene_ids)));
					$output .= sprintf("%s\t", join(',',@{$ensembl_ids->[$i]}));
					$output .= sprintf("%s\t", $uniprot_ids[$i]);
					$output .= sprintf("%s\t", $position_info[$i]);
					$output .= sprintf("%s\t", $variant_info[$i]);
					$output .= sprintf("%s\t", $class_info[$i]);
					$output .= sprintf("%s\t", $coverage_info[$i]);
					$output .= sprintf("%s\t", join(";",unique(@crossref_info)));
					$output .= sprintf("%s\t", join(";",unique(@drugs_info)));
					if ($function_info) {
						$output .= sprintf("%s\t", (split('\. ',$function_info))[0]);
					} else {
						$output .= "\t";
					}
					$output .= "\n";
				}
			}
		}

	}
	
	if ($format eq 'html') {
		$output .= "\t</tbody>\n</table>";
	}
	
	return "$output\n";

}

#################################################################################
 
# OBSOLETE:
# # Retrieves protein reference sequences from UniProt database
# sub download_amplicancer_refs {
# 
# 	my $outfile = shift @_;
# 	
# 	if (!defined($outfile)){
# 		$outfile = sprintf("/tmp/%s.fa.gz",random_file_name());
# 	}
# 	
# 	my $tmpdir = sprintf("/tmp/%s",random_file_name());
# 	rmdir($tmpdir);
# 
# 	`wget -q -P $tmpdir ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`;
# 	`wget -q -P $tmpdir ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz`;
# 	# Non-curated proteins:
# 	# `wget -q -P $tmpdir ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz`;
# 	
# 	# Keep only sequences from human
# 	system("zcat $tmpdir/*.fasta.gz | perl -ne '".'if(/^>.+OS=Homo sapiens/){$header=$_;print $_}elsif(/^>/){$header=0}elsif($header){print $_}'."' | gzip > $tmpdir/uniprot.fa.gz");
# # 	`zcat $tmpdir/*.fasta.gz | grep --no-group-separator -B1 -A2 'OS=Homo sapiens' | gzip > $tmpdir/uniprot.fa.gz`;
# 
# 	if (!is_fasta("$tmpdir/uniprot.fa.gz")) {
# 		`rm -rf $tmpdir`;
# 		return 0;
# 	} else {
# 		`rm $outfile*`; # Removes Blast indexes
# 		`mv $tmpdir/uniprot.fa.gz $outfile`;
# 		`rm -rf $tmpdir`;
# 		return $outfile;
# 	}
# 
# } 

#################################################################################

# Retrieves TXT format Uniprot data
sub retrieve_genbank_data {

	my ($seq_id,$format) = @_;
	
	if (!defined($format) || $format =~ /genbank/i){
		$format = 'gb'; # Genbank
	}

	# Ex. https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NM_005228&rettype=gb&retmode=text
	my $url;
	if ($format =~ /gb|fasta|xml/i){
		$url = sprintf("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=%s&retmode=text",uc($seq_id),lc($format));
	}
	my $response = HTTP::Tiny->new->get($url);
	if ($response->{success}) {
		return $response->{content};
	} else {
		print "\nERROR 'retrieve_genbank_data': '$seq_id' Genbank entry could not be retrieved.\n\n";
		return undef;
	}


}

#################################################################################

# Retrieves TXT format Uniprot data
sub retrieve_uniprot_data {

	my ($uniprot_id,$format) = @_;

	if (!defined($format)){
		$format = 'txt';
	}

	my $url;
	if ($format =~ /txt|fasta|xml/i){
		$url = sprintf("http://www.uniprot.org/uniprot/%s.%s",uc($uniprot_id),lc($format));
	}
	my $response = HTTP::Tiny->new->get($url);
	if ($response->{success}) {
		return $response->{content};
	} else {
		print "\nERROR 'retrieve_uniprot_data': '$uniprot_id' UniProt entry could not be retrieved.\n\n";
		return undef;
	}

}

#################################################################################

# Reads data from a TXT format Uniprot file
sub read_uniprot_single_file {

	my ($file) = @_;

	my $prot_data;

	# Reads data
	my ($prot_id, $count_biblio, $seq, @pubmed, @biblio, @cc_lines, $cc_type, @ft_lines);
	my $count_prot = 0;
	open(FILE,$file) || die "\n# $0 : cannot open '$file' file\n\n";
	while (my $row = <FILE>){
		$row =~ s/\012\015?|\015\012?//;
# 		if ($row =~ /^ID   (\w+)\s/) {
# PROVISIONAL:
		if ($row =~ /^ID   (\w+?)_/ || $row =~ /^ID   (\w+)\s/) {
			$prot_id = $1;
			$count_prot++;
		} elsif ($row =~ /^AC   (.+)/){
			push(@{$prot_data->{'accessions'}}, split(/\s*;\s*/, $1));
		} elsif ($row =~ /^DE   \w+Name: Full=(.+);/){
			push(@{$prot_data->{'names'}},$1);
		} elsif ($row =~ /^GN   Name=(.+?);( Synonyms=(.+);)?/){
			push(@{$prot_data->{'gene_names'}},$1);
			if (defined($3)){
				push(@{$prot_data->{'gene_names'}}, split(/\s*,\s*/, $3));
			}
		} elsif ($row =~ /^CC   -!- (.+?): (.+)/){
			if (@cc_lines){ # Annotates previous entry
				my $cc_info = join(' ',@cc_lines);
				my $pubmed_ids = read_uniprot_pubmed_ids($cc_info);
				my $omim_ids = read_uniprot_omim_ids($cc_info);
				if ($cc_info =~ /(.+) \{/){
					$cc_info = $1;
				}
				$cc_info =~ s/Note=//g;
				$cc_info =~ s/\s?\[.+?\]//g;
				if (defined($pubmed_ids)){
					$cc_info .= ' '.$pubmed_ids.';';
				}
				if (defined($omim_ids)){
					$cc_info .= ' '.$omim_ids.';';
				}
				push(@{$prot_data->{$cc_type}},sprintf("%s",$cc_info));
				undef($cc_type);
				undef(@cc_lines);
			}
			# Annotates only function and disease information
			if ($row =~ /^CC   -!- (FUNCTION|DISEASE): (.+)/){
				$cc_type = lc($1);
				push(@cc_lines, $2);
			}
		} elsif (@cc_lines && $row =~ /^CC\s{4,}(.+)/){
			push(@cc_lines, $1);
		} elsif ($row =~ /^DR   (.+)./){
			my @cols = split('; ',$1);
			if ($cols[0] eq 'EMBL' && $cols[4] eq 'mRNA'){
				push(@{$prot_data->{'RNA'}},$cols[1]);
			} elsif ($cols[0] eq 'EMBL' && $cols[4] eq 'Genomic_DNA'){
				push(@{$prot_data->{'DNA'}},$cols[1]);
			} elsif ($cols[0] eq 'PDB' && $cols[2] eq 'X-ray'){
				push(@{$prot_data->{'PDB'}},$cols[1]);
			} elsif ($cols[0] eq 'BindingDB'){
				push(@{$prot_data->{'BINDINGDB'}},$cols[1]);
			} elsif ($cols[0] eq 'ChEMBL'){
				$prot_data->{'CHEMBL'} = $cols[1];
			} elsif ($cols[0] eq 'DrugBank'){
				push(@{$prot_data->{'DRUGBANK'}},$cols[1]."; ".$cols[2]);
			} elsif ($cols[0] eq 'KEGG'){
				$prot_data->{'KEGG'} = $cols[1];
			} elsif ($cols[0] eq 'GuidetoPHARMACOLOGY'){
				$prot_data->{'IUPHAR'} = $cols[1];
			} elsif ($cols[0] eq 'BindingDB'){
				$prot_data->{'BINDINGDB'} = $cols[1];
			} elsif ($cols[0] eq 'BioMuta'){
				$prot_data->{'BIOMUTA'} = $cols[1];
			}
		} elsif ($row =~ /^FT   VARIANT\s+(\d+)\s+(\d+)\s+(.+)/){
			my $first_aa = $1;
			my $last_aa = $2;
			my $ft_line = $3;
			if (@ft_lines){ # Annotates previous entry
				my $variant_info = join(' ',@ft_lines);
				my $pubmed_ids = read_uniprot_pubmed_ids($variant_info);
				$variant_info =~ s/\(in dbSNP/( ; dbSNP/gi;
				if ($variant_info =~ /(.+?); \((.+)(; dbSNP:|\))/ ){
					$variant_info = sprintf("%s; %s;", $1, ucfirst($2));
				}
				if (defined($pubmed_ids)){
					$variant_info .= ' '.$pubmed_ids.';';
				}
				push(@{$prot_data->{'variants'}},sprintf("%s",$variant_info));
				undef(@ft_lines);
			}
			if ($ft_line =~ /^(.+?) (\(.+)/){
				push(@ft_lines, sprintf("%s; %s",read_uniprot_variant($1,$first_aa,$last_aa), $2));
			} else {
				push(@ft_lines, $ft_line);
			}
		} elsif (@ft_lines && $row =~ /^FT\s{4,}(.+)/){
			push(@ft_lines, $1);
		} elsif ($row =~ /^SQ   SEQUENCE/){
			$seq = '';
		} elsif (defined($seq) && $row =~ /^\s{5}(.+)/){
			$seq .= $1;
		} elsif ($row =~ /^\/\//){ # End of the file
			# Annotates sequence
			if (defined($seq)){
				$seq =~ s/\s//g;
				$prot_data->{'sequence'} = $seq;
			}
			# Annotates last entries
			if (@cc_lines){ # Annotates previous entry
				my $cc_info = join(' ',@cc_lines);
				my $pubmed_ids = read_uniprot_pubmed_ids($cc_info);
				my $omim_ids = read_uniprot_omim_ids($cc_info);
				if ($cc_info =~ /(.+) \{/){
					$cc_info = $1;
				}
				$cc_info =~ s/Note=//g;
				$cc_info =~ s/\s?\[.+?\]//g;
				if (defined($pubmed_ids)){
					$cc_info .= ' '.$pubmed_ids.';';
				}
				if (defined($omim_ids)){
					$cc_info .= ' '.$omim_ids.';';
				}
				push(@{$prot_data->{$cc_type}},sprintf("%s",$cc_info));
				undef($cc_type);
				undef(@cc_lines);
			}
			if (@ft_lines){ # Annotates previous entry
				my $variant_info = join(' ',@ft_lines);
				my $pubmed_ids = read_uniprot_pubmed_ids($variant_info);
				$variant_info =~ s/In dbSNP/dbSNP/gi;
				if ($variant_info =~ /(.+?); \((.+)(;dbSNP:|\))/ ){
					$variant_info = sprintf("%s; %s;", $1, ucfirst($2));
				}
				if (defined($pubmed_ids)){
					$variant_info .= ' '.$pubmed_ids.';';
				}
				push(@{$prot_data->{'variants'}},sprintf("%s",$variant_info));
				undef(@ft_lines);
			}
			undef($prot_id);
		}
	}
	close(FILE);

	return $prot_data;

}

# Reads PubMed IDs from Uniprot data
sub read_uniprot_variant {
	my ($data,$first_aa,$last_aa) = @_;
	if ($first_aa == $last_aa && $data =~ /^([A-Z]) -> ([A-Z])/){ # Amino acid substitution
		return sprintf("%s%d%s",$1,$first_aa,$2);
	} elsif ($data =~ /^Missing/){ # Amino acid deletion
		return sprintf("%d\_%ddel",$first_aa,$1);
	} elsif ($data =~ /^([A-Z]) -> ([A-Z])/){ # Amino acid insertion
		return sprintf("%d\_%dins%s",$first_aa,$last_aa,$2);
	} elsif ($data =~ /^([A-Z]+) -> ([A-Z]+)/){ # Amino acid deletion+insertion
		return sprintf("%d\_%ddelins%s",$first_aa,$last_aa,$2);
	} else {
		return $data;
	}
}

# Reads PubMed IDs from Uniprot data
sub read_uniprot_pubmed_ids {
	my ($data) = @_;
	my @pubmed_ids;
	while ($data =~ /PubMed:(\d+)/gi){
		push(@pubmed_ids,$1);
	}
	if (@pubmed_ids){
		return "PubMed:".join(',',@pubmed_ids);
	} else {
		return undef;
	}
}

# Reads OMIM IDs from Uniprot data
sub read_uniprot_omim_ids {
	my ($data) = @_;
	my @omim_ids;
	while ($data =~ /MIM:(\d+)/gi){
		push(@omim_ids,$1);
	}
	if (@omim_ids){
		return "OMIM:".join(',',@omim_ids);
	} else {
		return undef;
	}
}

#################################################################################

# Reads data from a Genbank format file
sub read_genbank_single_file {

	my ($file) = @_;

	my $seq_data;

	# Processes data
	my $seqIO = Bio::SeqIO->new( -file => $file, -format => 'genbank');
# 	my $seqIO = Bio::SeqIO->new( -string => $data, -format => 'genbank');
	while (my $seq = $seqIO->next_seq){
		# my $accession = $seq->accession_number;
		# my $version = $seq->version;

		foreach my $feat_object ($seq->get_SeqFeatures) {

			my $primary_tag = $feat_object->primary_tag;

			if ($feat_object->strand == -1) {
				$seq_data->{$primary_tag}{'location_gene'} = "c".$feat_object->location->end."-".$feat_object->location->start;
			} else {
				$seq_data->{$primary_tag}{'location_gene'} = $feat_object->location->start."-".$feat_object->location->end;
			}
			my @location_;
			if ($feat_object->location->isa('Bio::Location::SplitLocationI')){
				for my $location_object ( $feat_object->location->sub_Location ) {
					push(@location_,$location_object->start."-".$location_object->end);
				}
				if ($feat_object->strand == -1) {
					@location_ = reverse(@location_);
				}
			} else {
				push(@location_,$feat_object->location->start."-".$feat_object->location->end);
			}
			$seq_data->{$primary_tag}{'location_tra'} = join(',',@location_);

			if ($primary_tag eq 'source'){
				$seq_data->{$primary_tag}{'organism'} = join('',$feat_object->get_tag_values('organism'));
			}

			if ($primary_tag eq 'gene'){
				if (grep { $_ eq 'db_xref' } $feat_object->get_all_tags){
					foreach my $value ($feat_object->get_tag_values('db_xref')){
						if ($value =~ /^GI:(\d+)/){
							$seq_data->{$primary_tag}{'gene_id'} = $1;
						}
					}
				}
				if (grep { $_ eq 'note' } $feat_object->get_all_tags){
					$seq_data->{$primary_tag}{'product'} = join('',$feat_object->get_tag_values('note'));
				}
				if (grep { $_ eq 'gene' } $feat_object->get_all_tags){
					$seq_data->{$primary_tag}{'gene_name'} = join('',$feat_object->get_tag_values('gene'));
				}
				my $gene_seq = $feat_object->seq->seq;
				if(length($gene_seq) > 0) {
					# $gene_seq =~ s/(.{1,70})/$1\n/g;
					$seq_data->{$primary_tag}{'sequence'} = $gene_seq;
				}
			} elsif ($primary_tag =~ /RNA$/ && $primary_tag ne "mRNA"){ #($primary_tag =~ /RNA$/){
				if ($primary_tag eq "tRNA" && grep { $_ eq 'product' } $feat_object->get_all_tags){
					my $product = join('',$feat_object->get_tag_values('product'));
					if ($product =~ /^(tRNA)-(.+)$/){
						$seq_data->{$primary_tag}{'product'} = "$2 $1";
					}
				}
				if (grep { $_ eq 'gene' } $feat_object->get_all_tags){
					$seq_data->{$primary_tag}{'gene_name'} = join('',$feat_object->get_tag_values('gene'));
				}
				if (grep { $_ eq 'locus_tag' } $feat_object->get_all_tags){
					$seq_data->{$primary_tag}{'locus_tag'} = join('',$feat_object->get_tag_values('locus_tag'));
				}
				my $transcription = $feat_object->spliced_seq->seq;
				if(length($transcription) > 0) {
					# $transcription =~ s/(.{1,70})/$1\n/g;
					$seq_data->{$primary_tag}{'sequence'} = $transcription;
				}
			} elsif ($primary_tag eq 'CDS'){
# Stops the script if more than one CD is found
if (defined($seq_data->{$primary_tag}{'product'})){
	print "\nERROR 'read_genbank_single_file': Multiple CDS in file '$file'.\n\n";
	exit;
}
				if (grep { $_ eq 'db_xref' } $feat_object->get_all_tags){
					foreach my $value ($feat_object->get_tag_values('db_xref')){
						if ($value =~ /^GI:(\d+)/){
							$seq_data->{$primary_tag}{'gene_id'} = $1;
						}
					}
				}
				if (grep { $_ eq 'product' } $feat_object->get_all_tags){
					$seq_data->{$primary_tag}{'product'} = join('',$feat_object->get_tag_values('product'));
				}
				if (grep { $_ eq 'gene' } $feat_object->get_all_tags){
					$seq_data->{$primary_tag}{'gene_name'} = join('',$feat_object->get_tag_values('gene'));
				}
				if (grep { $_ eq 'locus_tag' } $feat_object->get_all_tags){
					$seq_data->{$primary_tag}{'locus_tag'} = join('',$feat_object->get_tag_values('locus_tag'));
				}
				if (grep { $_ eq 'protein_id' } $feat_object->get_all_tags){
					$seq_data->{$primary_tag}{'protein_id'} = join('',$feat_object->get_tag_values('protein_id'));
				}
				# Codificant Domain (CD) cDNA
				my $CD = $feat_object->spliced_seq->seq;
				# Protein sequence
				my $translation;
				if(length($CD) > 0) {
					my $seq_obj = Bio::Seq->new(-seq => $CD,
								-alphabet => 'dna' );
					$translation = $seq_obj->translate->seq;
					$translation =~ s/\*$//;
					# $translation =~ s/(.{1,70})/$1\n/g;
					# $CD =~ s/(.{1,70})/$1\n/g;

					$seq_data->{$primary_tag}{'sequence'} = $CD;
					$seq_data->{$primary_tag}{'translation'} = $translation;
				}
			}
		}
	}

	return $seq_data;

}


#################################################################################

# Reads Uniprot ID mapping information file (cross-references to external databases)
sub read_uniprot_idmapping_file {

	my ($file,$database,$sense) = @_;

	my $id_data;
	
	# Sense 1 means Uniprot => Ext_database
	# Sense -1 means Ext_database => Uniprot
	if (!defined($sense)){
		$sense = 1;
	}

	if (is_gzip($file)){
		open(FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_zip($file)){
		open(FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_bzip2($file)){
		open(FILE, "bzcat '$file' |") or die "# cannot read $file\n";
	} else {
		open(FILE,$file) || die "# cannot read $file\n";
	}
	while (<FILE>){
		s/\012\015?|\015\012?//;
		my @cols = split("\t");
		if (!defined($database) || $cols[1] eq $database) {
			if ($sense>0){
				push(@{$id_data->{$cols[0]}}, $cols[2]);
			} else {
				push(@{$id_data->{$cols[2]}}, $cols[0]);
			}
		}
	}
	close(FILE);

	return $id_data;

	# Example data:
	# P00519  UniProtKB-ID    ABL1_HUMAN
	# P00519  Gene_Name       ABL1
	# P00519  Gene_Synonym    ABL
	# P00519  Gene_Synonym    JTK7
	# P00519  GI      576865144
	# P00519  GI      185178034
	# ...
	# P00519  UniRef100       UniRef100_P00519
	# P00519  UniRef90        UniRef90_P00519
	# P00519  UniRef50        UniRef50_P00519
	# P00519-1        UniParc UPI000013D6D4
	# P00519-2        UniParc UPI000013E4DE
	# P00519  UniParc UPI000013D6D4
	# P00519  EMBL    M14752
	# P00519  EMBL    X16416
	# ...
	# P00519  EMBL-CDS        AAA51561.1
	# P00519  EMBL-CDS        CAA34438.1
	# ...
	# P00519  EMBL-CDS        AAA51561.1
	# P00519  EMBL-CDS        CAA34438.1
	# ...
	# P00519  DIP     DIP-1042N
	# P00519  MINT    MINT-7236141
	# P00519  STRING  9606.ENSP00000361423
	# P00519  ChEMBL  CHEMBL1862
	# P00519  DrugBank        DB08043
	# P00519  DrugBank        DB00171
	# ...
	# P00519  GuidetoPHARMACOLOGY     1923
	# P00519  BioMuta ABL1
	# P00519  DMDM    85681908
	# P00519  DNASU   25
	# P00519  Ensembl ENSG00000097007
	# P00519-1        Ensembl_TRS     ENST00000318560
	# P00519-1        Ensembl_PRO     ENSP00000323315
	# P00519-2        Ensembl_TRS     ENST00000372348
	# P00519-2        Ensembl_PRO     ENSP00000361423
	# P00519  GeneID  25
	# P00519  KEGG    hsa:25
	# P00519-1        UCSC    uc004bzv.4
	# P00519  GeneCards       ABL1
	# P00519  HGNC    HGNC:76
	# P00519  HPA     CAB002686
	# P00519  HPA     HPA027251
	# P00519  HPA     HPA027280
	# P00519  HPA     HPA028409
	# P00519  MIM     189980
	# P00519  MIM     608232
	# P00519  neXtProt        NX_P00519
	# P00519  Orphanet        521
	# P00519  Orphanet        99860
	# P00519  Orphanet        99861
	# P00519  PharmGKB        PA24413
	# ...

}

#################################################################################

# Reads .tsv file with COSMIC database information about gene mutations involved in cancer
sub read_cosmic_data_file {

	my ($file,$ids,$fields) = @_;

	my $data;
	
	# Extracts only specific Ensembl transcript IDs data
	if (defined($ids) && @$ids){
		$ids = { array_to_hash(@$ids) };
	} else {
		$ids = {};
	}
	# Extracts only specific fields data
	if (defined($fields) && @$fields){
		$fields = { array_to_hash(@$fields) };
	} else {
		$fields = {};
	}

	if (is_gzip($file)){
		open(FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_zip($file)){
		open(FILE, "zcat '$file' |") or die "# cannot read $file\n";
	} elsif (is_bzip2($file)){
		open(FILE, "bzcat '$file' |") or die "# cannot read $file\n";
	} else {
		open(FILE,$file) || die "# cannot read $file\n";
	}
	my %headers;
	while (<FILE>) {
		if (/^#/){
			next;
		}
		s/\012\015?|\015\012?//; # trim line break
		my @cols = split("\t");
		# Retrieves column headers from the first line
		if (!%headers && $cols[1] =~ /Accession/){
			for (my $i=0; $i<=$#cols; $i++){
				# If we desire only specific fields data
				if (!defined($fields->{$cols[$i]})){
					next;
				}
				$headers{$i} = $cols[$i];
			}
			next;
		}
		my $ensembl_id = $cols[1];
		# If we desire only specific Ensembl transcript IDs data
		if (!defined($ids->{$ensembl_id})){
			next;
		}
		my %data_;
		foreach my $i (keys %headers){
			$data_{$headers{$i}} = $cols[$i];
		}
		push(@{$data->{$ensembl_id}}, \%data_);
	}
	close FILE;
	
	return $data;

	# Example data:
	#
	# Gene name	Accession Number	Gene CDS length	HGNC ID	Sample name	ID_sample	ID_tumour	Primary site	Site subtype 1	Site subtype 2	Site subtype 3	Primary histology	Histology subtype 1	Histology subtype 2	Histology subtype 3	Genome-wide screen	Mutation ID	Mutation CDS	Mutation AA	Mutation Description	Mutation zygosity	LOH	GRCh	Mutation genome position	Mutation strand	SNP	Resistance Mutation	FATHMM prediction	FATHMM score	Mutation somatic status	Pubmed_PMID	ID_STUDY	Sample source	Tumour origin	Age
	# APBB1IP	ENST00000376236	2001	17379	NCI-H1184	688000	616160	lung	NS	NS	NS	carcinoma	small_cell_carcinoma	NS	NS	n	COSM18232	c.1089C>T	p.N363N	Substitution - coding silent	het	n	38	10:26541626-26541626	+	n	-	PATHOGENIC	0.83092	Confirmed somatic variant		24	cell-line	metastasis	42
	# KRAS	ENST00000311936	567	6407	S32338	711700	637403	large_intestine	colon	sigmoid	NS	adenoma	NS	NS	NS	n	COSM521	c.35G>A	p.G12D	Substitution - Missense		y	38	12:25245350-25245350	-	n	-	PATHOGENIC	0.97875	Reported in another cancer sample as somatic	11215737		surgery-fixed	primary	54
	# KRAS	ENST00000311936	567	6407	S47542	727500	651688	large_intestine	NS	NS	NS	adenoma	NS	NS	NS	n	COSM520	c.35G>T	p.G12V	Substitution - Missense		y	38	12:25245350-25245350	-	n	-	PATHOGENIC	0.98367	Reported in another cancer sample as somatic	12438234		surgery fresh/frozen	primary	
	# KRAS	ENST00000311936	567	6407	S50521	730700	654758	pancreas	NS	NS	NS	carcinoma	ductal_carcinoma	NS	NS	n	COSM521	c.35G>A	p.G12D	Substitution - Missense		y	38	12:25245350-25245350	-	n	-	PATHOGENIC	0.97875	Reported in another cancer sample as somatic	8121190		surgery-fixed	primary	55
	# KRAS	ENST00000311936	567	6407	S54724	735100	658746	large_intestine	NS	NS	NS	carcinoma	adenocarcinoma	NS	NS	n	COSM516	c.34G>T	p.G12C	Substitution - Missense		y	38	12:25245351-25245351	-	n	-	PATHOGENIC	0.98367	Confirmed somatic variant	2666051		surgery fresh/frozen	primary	
	# KRAS	ENST00000311936	567	6407	S71858	749600	666189	thyroid	NS	NS	NS	carcinoma	follicular_carcinoma	NS	NS	n	COSM532	c.38G>A	p.G13D	Substitution - Missense		y	38	12:25245347-25245347	-	n	-	PATHOGENIC	0.97875	Reported in another cancer sample as somatic	12947056		surgery-fixed	primary	
	# KIT	ENST00000288135	2931	6342	E36297	819800	739677	soft_tissue	fibrous_tissue_and_uncertain_origin	stomach	NS	gastrointestinal_stromal_tumour	spindle	NS	NS	n	COSM1349	c.?_?ins?	p.R588_L589insDHKWEFPRN	Insertion - In frame		y					-			Variant of unknown origin	12960119		surgery fresh/frozen	NS	
	# CTNNB1	ENST00000349496	2346	2514	E5297	869700	788943	large_intestine	NS	NS	NS	carcinoma	adenocarcinoma	NS	NS	n	COSM6019	c.14_241del228	p.A5_A80del	Deletion - In frame		n	38	3:41224526-41224753	+		-			Variant of unknown origin	9500465		NS	primary	
	# KRAS	ENST00000311936	567	6407	S80668	875600	794842	large_intestine	colon	sigmoid	NS	carcinoma	adenocarcinoma	NS	NS	n	COSM516	c.34G>T	p.G12C	Substitution - Missense	hom	u	38	12:25245351-25245351	-	n	-	PATHOGENIC	0.98367	Confirmed somatic variant	14633673		surgery - NOS	primary	71


}

#################################################################################



1;
