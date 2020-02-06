#!/usr/bin/perl -w
#
################################################################
#
# Name: ampliMIX.pl
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
#   Mixes several FASTA/FASTQ files with same samples and different tags specified in CSV files
#   Generates a unique FASTA/FASTQ format file with unique tags per sample
#
# Requires as input several FASTA/FASTQ files and several CVS files in the same order with tag/sample data
#
# Example:
# perl ampliMIX.pl -d primers_tags1.cvs primers_tags2.cvs primers_tags3.cvs -i exp1.fq.gz exp2.fq.gz exp3.fq.gz -o mixed
#

my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Combines amplicon sequencing data from two different experiments/runs that share common tag sequences and/or samples.";

# Modules are in folder 'lib' in the path of the script
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;


my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options
# Maximum number of allowed primer+tag sequences
# my $MAX_PRIMER_TAG_SEQS = 10000;

my (@INP_amplicons_files, @INP_reads_files, $INP_outfile, $INP_direct, $INP_zip, $INP_gzip);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s{,}' => \@INP_reads_files,
	'd|data=s{,}' => \@INP_amplicons_files,
	'o|output=s' => \$INP_outfile,
	'di|direct' => \$INP_direct,
	'z|zip' => \$INP_zip,
	'gz|gzip' => \$INP_gzip,
	'<>' => \&usage,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file1> <file2> [...<fileN>] -d <file1> <file2> [...<fileN>] [options]\n";
	print "\nOptions:\n";
	print "  -i <file1> <file2> [...<fileN>]\n\t\tInput FASTQ or FASTA file (compressed or uncompressed).\n";
	print "  -d <file1> <file2> [...<fileN>]\n\t\tCSV file with primer/amplicon data.\n";
	print "  -o <name>\tOutput file name.\n";
	print "  -di\t\tAnalyze reads only in direct sense.\n";
	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -gz\t\tCompress results in GZIP format.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}


# Prints usage help if no input file is specified
if (!@INP_reads_files || !@INP_amplicons_files){
	print "\nERROR: You must specify input files.\n\n";
	usage();
	exit;
}
my ($INP_reads_outfile,$INP_amplicons_outfile);
if (!defined($INP_outfile)){
	if ($INP_reads_files[0] =~ /(.+)\./){
		$INP_reads_outfile = $1.".comb";
	} else {
		$INP_reads_outfile = $INP_reads_files[0].".comb";
	}
	if ($INP_amplicons_files[0] =~ /(.+)\./){
		$INP_amplicons_outfile = $1.".comb";
	} else {
		$INP_amplicons_outfile = $INP_amplicons_files[0].".comb";
	}
} else {
	$INP_reads_outfile = $INP_outfile;
	$INP_amplicons_outfile = $INP_outfile;
}

print "\nRunning '$COMMAND_LINE'\n";


my ($mixed_read_seqs,$mixed_read_headers,$mixed_read_qualities,%previous_samples,%previous_tags);
my ($mixed_markerdata,$mixed_sampledata) = ({},{});

my $mixed_reads_file_format = 'fastq';

for (my $i=0; $i<=$#INP_reads_files; $i++) {

	my $INP_reads_file = $INP_reads_files[$i];
	my $INP_amplicons_file;
	if (defined($INP_amplicons_files[$i])) {
		$INP_amplicons_file = $INP_amplicons_files[$i];
	} else { # If all the amplicon data is in one file
		$INP_amplicons_file = $INP_amplicons_files[0];
	}

	# my ($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads);
	# Check and read sequences file
	my ($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads)
	= parse_sequence_file($INP_reads_file,undef,['qualities','verbose']);
	# If any of the files is not in FASTQ format
	if ($reads_file_format ne $mixed_reads_file_format){
		$mixed_reads_file_format = 'fasta';
	}

	# 2 PRIMERS => MARKER
	# 1/2 TAGS => SAMPLE
	# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

	# Check and read amplicons file
	my ($markerdata,$primers,$sampledata,$tags,$paramsdata,$alleledata)
	= parse_amplicon_file($INP_amplicons_file,['verbose']);
# 	if (scalar @{$primer_tag_headers} > $MAX_PRIMER_TAG_SEQS ){
# 		print "\nERROR: Too many primer and/or tag sequences to be processed, decrease the number of primers and/or samples in the analysis.\n\n";
# 		exit;
# 	}

	print "\nMixing data from '$INP_reads_file' and '$INP_amplicons_file'.\n";
	my $count_mixed = 0;

	# Include marker data
	#$mixed_markerdata = { %$mixed_markerdata, %$markerdata };
	foreach my $marker (keys %$markerdata) {
		foreach my $field (keys %{$markerdata->{$marker}}) {
			if (in_array(['primer_f', 'primer_r', 'primer_rc', 'length'], $field)){
				foreach my $value (@{$markerdata->{$marker}{$field}}){
					if (!defined($mixed_markerdata->{$marker}) ||  !defined($mixed_markerdata->{$marker}{$field}) || !in_array($mixed_markerdata->{$marker}{$field},$value) ){
						if ($i>0 && defined($mixed_markerdata->{$marker}) && defined($mixed_markerdata->{$marker}{$field})){
							printf("#WARNING: '%s' '%s' value '%s' is added to previous values.\n",$marker,$field,$value);
						}# elsif (!defined($mixed_markerdata->{$marker}) && !defined($mixed_markerdata->{$marker}{$field})){
						
						push(@{$mixed_markerdata->{$marker}{$field}},$value);
					}
				}
			} else {
				my $value = $markerdata->{$marker}{$field};
				if (!defined($mixed_markerdata->{$marker}) ||  !defined($mixed_markerdata->{$marker}{$field})){
					$mixed_markerdata->{$marker}{$field} = $value;
				} elsif ($i>0 && defined($mixed_markerdata->{$marker}) && defined($mixed_markerdata->{$marker}{$field})){
					if (ref($mixed_markerdata->{$marker}{$field}) eq 'ARRAY' && ref($value) eq 'ARRAY') {
						if (join(',',@{$mixed_markerdata->{$marker}{$field}}) ne join(',',@$value) ){
							printf("#WARNING: '%s' '%s' values '%s' are different to previous values '%s'.\n",$marker,$field,join(',',@$value),join(',',@{$mixed_markerdata->{$marker}{$field}}));
						}
					} elsif ($mixed_markerdata->{$marker}{$field} ne $value){
						printf("#WARNING: '%s' '%s' value '%s' is different to previous value '%s'.\n",$marker,$field,$value,$mixed_markerdata->{$marker}{$field});
					}
				}
			}
		}
	}

	# Include tag data and reorganize tags
	my $tag_replacements;
	foreach my $sample (@$tags) {
		my $modified_tags = 0;
		# If the sample is in a previous file with different tags, change tags to the previous one
		if (defined($previous_samples{$sample})){
			if ($sampledata->{$sample}{'tag_f'} ne $mixed_sampledata->{$sample}{'tag_f'}
			|| $sampledata->{$sample}{'tag_r'} ne $mixed_sampledata->{$sample}{'tag_r'} ){
				$modified_tags = 1;
			}
		} else {
			# If the sample tags are already used in a previous file, change tags to new ones
			$mixed_sampledata->{$sample} = { %{$sampledata->{$sample}} };
			while (defined($previous_tags{sprintf("%s-%s",$mixed_sampledata->{$sample}{'tag_f'},$mixed_sampledata->{$sample}{'tag_r'})})){
				$mixed_sampledata->{$sample}{'tag_f'} = generate_tag(length($sampledata->{$sample}{'tag_f'}));
				$mixed_sampledata->{$sample}{'tag_r'} = generate_tag(length($sampledata->{$sample}{'tag_r'}));
				$mixed_sampledata->{$sample}{'tag_rc'} = iupac_reverse_complementary($mixed_sampledata->{$sample}{'tag_r'});
				$modified_tags = 1;
			}
		}
		# @$tag_replacements stores in indexes 0 and 1 the old tags (fwd and revcomp) to be replaced and in 2 and 3 the new ones
		# First replacement corresponds to direct DNA sense and second to reverse sense.
		if ($modified_tags){
			$tag_replacements->{$sample}{'tag_f'} = $sampledata->{$sample}{'tag_f'};
			$tag_replacements->{$sample}{'tag_r'} = $sampledata->{$sample}{'tag_r'};
			$tag_replacements->{$sample}{'tag_rc'} = $sampledata->{$sample}{'tag_rc'};
		}
		$previous_samples{$sample} = 1;
		$previous_tags{sprintf("%s-%s",$mixed_sampledata->{$sample}{'tag_f'},$mixed_sampledata->{$sample}{'tag_r'})} = $sample;
	}
	
	# Replaces old tags
	foreach my $marker (@$primers){
		my @primers_f = @{$markerdata->{$marker}{'primer_f'}};
		my @primers_rc = @{$markerdata->{$marker}{'primer_rc'}};
		foreach my $sample (@$tags) {

			print "\t$marker-$sample merging\n";

			my ($tag_f, $tag_r, $tag_rc, $tag_fc) = ('', '', '', '');
			
			if (defined($tag_replacements->{$sample})){
				$tag_f = $tag_replacements->{$sample}{'tag_f'};
				$tag_r = $tag_replacements->{$sample}{'tag_r'};
				$tag_rc = $tag_replacements->{$sample}{'tag_rc'};
				$tag_fc = iupac_reverse_complementary($tag_replacements->{$sample}{'tag_f'});
			} else {
				$tag_f = $sampledata->{$sample}{'tag_f'};
				$tag_r = $sampledata->{$sample}{'tag_r'};
				$tag_rc = $sampledata->{$sample}{'tag_rc'};
				$tag_fc = iupac_reverse_complementary($sampledata->{$sample}{'tag_f'});
			}
			my (@patterns_dir, @patterns_rev);
			foreach my $primer_f (@primers_f) {
				foreach my $primer_rc (@primers_rc) {
					# Very important the parenthesis in the regex for later locate the start position and the length of the sequence between primers and tags
					push(@patterns_dir, regex("$tag_f($primer_f.+$primer_rc)$tag_rc"));
					unless (defined($INP_direct)) {
						push(@patterns_rev, regex("$tag_r(".iupac_reverse_complementary($primer_rc).".+".iupac_reverse_complementary($primer_f).")$tag_fc"));
					}
				}
			}
			# $i is already in use to loop the input files
			for (my $j=0; $j<=$#{$read_seqs}; $j++){
				my $mixed_seqs = 0;
				for (my $p=0; $p<=$#patterns_dir; $p++) {
					if ($read_seqs->[$j] =~ /$patterns_dir[$p]/){
						$read_seqs->[$j] = $`.$mixed_sampledata->{$sample}{'tag_f'}.$1.$mixed_sampledata->{$sample}{'tag_rc'}.$';
						$mixed_seqs = 1;
						# $read_headers->[$j] = "MIXED $sample ".$read_headers->[$j];
					} elsif (!defined($INP_direct) && $read_seqs->[$j] =~ /$patterns_rev[$p]/){
						$read_seqs->[$j] = $`.$mixed_sampledata->{$sample}{'tag_r'}.$1.iupac_reverse_complementary($mixed_sampledata->{$sample}{'tag_f'}).$';
						$mixed_seqs = 1;
						# $read_headers->[$j] = "MIXED $sample ".$read_headers->[$j];
					}
					if ($mixed_seqs){
						$count_mixed++;
						push(@$mixed_read_seqs,$read_seqs->[$j]);
						push(@$mixed_read_headers,$read_headers->[$j]);
						if ($mixed_reads_file_format eq 'fastq'){
							push(@$mixed_read_qualities,$read_qualities->[$j]);
						}
						last;
					}
				}
# 				if (!$mixed_seqs) {
# 					splice(@$read_seqs,$j,1);
# 					splice(@$read_headers,$j,1);
# 					splice(@$read_qualities,$j,1);
# 					$j--;
# 				}
			}
			print "\t$marker-$sample mixed\n";
		}
	}
	
	print "\nMixed $count_mixed sequences.\n";


# 	# To save all the remaining reads, but they may contain the same tags than the mixed ones and fuck up the resulting file
# 	push(@$mixed_read_seqs,@$read_seqs);
# 	push(@$mixed_read_headers,@$read_headers);
# 	if ($mixed_reads_file_format eq 'fastq'){
# 		push(@$mixed_read_qualities,@$read_qualities);
# 	}

}

# Prints amplicon parameters into a file
my $mixed_amplicon_data;
$mixed_amplicon_data .= print_amplicon_data($mixed_markerdata,'markers')."\n";
$mixed_amplicon_data .= print_amplicon_data($mixed_sampledata,'samples')."\n";
write_to_file("$INP_amplicons_outfile.csv",$mixed_amplicon_data);
print "\nSaved mixed amplicon data into '$INP_amplicons_outfile.csv'.\n";

# Defines compression method
my $compression;
if (defined($INP_zip)){
	$compression = 'zip';
} elsif (defined($INP_gzip)){
	$compression = 'gzip';
}

# Prints mixed sequences into a file
my $outfile;
if (@$mixed_read_seqs && $mixed_reads_file_format eq 'fastq'){
	$outfile = create_fastq_file($mixed_read_seqs,$mixed_read_headers,$mixed_read_qualities,"$INP_reads_outfile.fq",$compression);
} elsif (@$mixed_read_seqs) {
	$outfile = create_fasta_file($mixed_read_seqs,$mixed_read_headers,"$INP_reads_outfile.fa",$compression);
}

if (defined($outfile)){
	printf("\nSaved %d mixed sequences into '%s'.\n", scalar @$mixed_read_seqs, $outfile);
} else {
	print "\nThere was some error in the merging process and no sequences were retrieved.\n\n";
}

# if (defined($INP_zip) && defined($outfile)){
# 	`zip -jqm $outfile.zip $outfile` ;
# 	printf("\nSaved %d mixed sequences into '%s'.\n", scalar @$mixed_read_seqs, "$outfile.zip");
# } elsif (defined($INP_gzip) && defined($outfile)){
# 	`gzip $outfile` ;
# 	printf("\nSaved %d mixed sequences into '%s'.\n", scalar @$mixed_read_seqs, "$outfile.gz");
# } elsif (defined($outfile)){
# 	printf("\nSaved %d mixed sequences into '%s'.\n", scalar @$mixed_read_seqs, $outfile);
# } else {
# 	print "\nThere was some error in the merging process and no sequences were retrieved.\n\n";
# }

print "\n";



























