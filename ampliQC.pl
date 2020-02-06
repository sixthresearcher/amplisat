#!/usr/bin/perl -w
#
# Name: ampliQC.pl
#
# Version: 1.0
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
# Analyzes quality of data in amplicon sequencing (AS) experiments.
# Quality assesment is performed comparing real primer sequences to experimental sequences or reads.
# Is independent of the platform and applicable to other several NGS methodologies when we previously know the sequences of used primers and tags.
# Gives as output an Excel file with quality report, read/amplicon length frequencies and read/primer alignments
# and a FASTA format file with primer sequences used in the analysis
#
# Requires as input a CVS format file with primer and amplicon data and a FASTQ file with sequences to analyze
#
# Example:
# perl ampliQC.pl -d amplicons.cvs -i reads.fq -o report
#

my $VERSION = "1.2";
my $SCRIPT_NAME = "ampliQC.pl";
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Analyzes the quality of the data in amplicon sequencing experiments.";


# Modules are in folder 'lib' in the path of the script
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Getopt::Long;
use Excel::Writer::XLSX;
use Excel::Writer::XLSX::Utility;
use Bio::Sequences;
use Bio::Ampli;
# use Data::Dumper;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

my ($INP_amplicons_file, $INP_reads_file, $INP_threads, $INP_nreads, $INP_direct, $INP_bootstrap, $INP_full, $INP_verbose, $INP_test);

# Default options
# Default output filename
my $INP_outfile = $SCRIPT_NAME;
# Default alignment algorithm
my $INP_align = 'gassst';
# Maximum number of allowed primer+tag sequences
my $MAX_PRIMER_TAG_SEQS = 10000;
# Maximum number of allowed sequences
my $MAX_READS_NUMBER = 200000;

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s' => \$INP_reads_file,
	'd|data=s' => \$INP_amplicons_file,
	'o|output:s' => \$INP_outfile,
	'a|align:s' => \$INP_align,
	'di|direct' => \$INP_direct,
	'n|number=i' => \$INP_nreads,
	'r|rounds=i' => \$INP_bootstrap,
	'f|full' => \$INP_full,
	'thr|threads=i' => \$INP_threads,
	'v|verbose' => \$INP_verbose,
# 	'test' => \$INP_test,
	'<>' => \&usage,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -i <file>\tInput FASTQ or FASTA file (compressed or uncompressed).\n";
	print "  -d <file>\tCSV file with primer/amplicon data.\n";
	print "  -o <path>\tOutput folder name.\n";
	print "  -a <type>\tPrimer alignment strategy: match, gassst, blast (default=$INP_align).\n";
	print "  -n <number>\tNumber of reads/sequences to analyze per round.\n";
	print "  -r <number>\tNumber of resampling rounds.\n";
	print "  -f\t\tAnalyze the full set of sequences.\n";
	print "  -v\t\tVerbose output, prints additional details about sequence alignments.\n";
	print "  -di\t\tAnalyze reads only in direct sense.\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
# 	print "  -z\t\tCompress results in ZIP format.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}

# Prints usage help if no input file is specified
if (!defined($INP_reads_file) || !defined($INP_amplicons_file)){
	print "\nERROR: You must specify input files.\n";
	usage();
	exit;
}


print "\nRunning '$COMMAND_LINE'\n";

# Assign dump file, only used in case of -test parameter chosen
my $dump_file = "$INP_outfile.ampli.dump";

# 2 PRIMERS => MARKER
# 1/2 TAGS => SAMPLE
# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

# my ($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads);
# Checks sequences file
my ($reads_file_format,$total_reads)
= parse_sequence_file($INP_reads_file,undef,['stats']);

# Sets default number of reads to 100000:
if (!defined($INP_full) && !defined($INP_nreads) && $total_reads >= 100000){
	print "\nSetting number of sequences to analyze per round to 10000 (default).\n";
	$INP_nreads = 10000;
} elsif (!defined($INP_full) && !defined($INP_nreads)){
	print "\nSetting number of sequences to analyze per round to $total_reads (total).\n";
	$INP_nreads = $total_reads;
	$INP_full = 1;
} elsif (!defined($INP_full) && defined($INP_nreads) && $total_reads < $INP_nreads){
	print "\nERROR: Number of sequences to analyze must be lower or equal to the number of reads in the file '$INP_reads_file' ($total_reads).\n\n";
	exit;
}

# Sets default alignment type to GASSST
if (!defined($INP_align)){
	$INP_align = 'gassst';
}

# Sets default aligment of DNA in both senses
my $INP_revcomp = 1;
if (!defined($INP_direct)){
	$INP_revcomp = 0;
}

# Sets default resampling to 10:
if (!defined($INP_full)){
	if (!defined($INP_bootstrap)){
		print "\nSetting number of resampling rounds to 10 (default).\n";
		$INP_bootstrap = 10;
	}
} else {
	$INP_nreads = $total_reads;
	$INP_bootstrap = 1;
}

if ($INP_nreads >= $MAX_READS_NUMBER){
	print "\tERROR: Analyzing '$INP_nreads' sequences can be slow and use too much memory, it's recommended to analyze a lower number of sequences with resampling option.\n";
	exit;
}


# 2 PRIMERS => MARKER
# 1/2 TAGS => SAMPLE
# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

# Check and read sequences/amplicons file
my ($find_seqs, $find_amplicons);
my ($markerdata,$primers,$sampledata,$tags,$primer_seqs,$primer_headers,$primer_tag_seqs,$primer_tag_headers,$paramsdata,$alleledata);
if (is_fasta($INP_amplicons_file)){
	my ($seqs,$headers) = read_fasta_file($INP_amplicons_file);
	for (my $i=0; $i<=$#{$headers}; $i++) {
		my $count_seqs = 0;
		push(@$primers,$headers->[$i]);
		my @seqs = unambiguous_dna_sequences($seqs->[$i]);
		if ($#seqs==0) {
			push(@$primer_headers,$headers->[$i]);
			push(@$primer_seqs,$seqs->[$i]);
		} else {
			foreach my $seq (@seqs){
				$count_seqs++;
				push(@$primer_headers,sprintf('%s_%03d', $headers->[$i], $count_seqs));
				push(@$primer_seqs,$seq);
				push(@{$markerdata->{$headers->[$i]}{'primer_f'}}, $seq);
			}
		}
	}
	print "\tNumber of unique primer+tag sequences: ".scalar @{$primer_headers}.".\n";
	if (scalar @{$primer_headers} > $MAX_PRIMER_TAG_SEQS ){
		print "\nERROR: Too many reference sequences to be processed, decrease the number of reference sequences in the analysis.\n\n";
		exit;
	}
	$find_seqs = 1;
	$find_amplicons = 0;
} else {
	($markerdata,$primers,$sampledata,$tags,$paramsdata,$alleledata) 
	= parse_amplicon_file($INP_amplicons_file,['verbose']);
# ERROR: should extract all primer+tag combinations
	($primer_tag_seqs, $primer_tag_headers) = extract_primer_tag_seqs($markerdata, $primers);
	print "\tNumber of unique primer+tag sequences: ".scalar @{$primer_tag_headers}.".\n";
	if (scalar @{$primer_tag_headers} > $MAX_PRIMER_TAG_SEQS ){
		print "\nERROR: Too many primer and/or tag sequences to be processed, decrease the number of primers and/or samples in the analysis.\n\n";
		exit;
	}
	$find_seqs = 0;
	$find_amplicons = 1;
	# Checks if $markerdata contains markers or normal sequences
	# To look for amplicons or single sequence matches
	foreach my $markername (keys %$markerdata) {
		# If one marker has only one primer sequence, no amplicons will be searched
		if ($#{$markerdata->{$markername}{'primer_f'}} < 0 || $#{$markerdata->{$markername}{'primer_r'}} < 0){
			$find_amplicons = 0;
		}
	}
}

# Stores primers into a hash
my (%primer_seqs_hash, %tag_seqs_hash);
map $primer_seqs_hash{$primer_headers->[$_]} = $primer_seqs->[$_], 0 .. $#{$primer_headers};
if (defined($tags)) {
	# Stores tags into a hash
	map $tag_seqs_hash{$_} = $sampledata->{$_}{'tag_f'}.'...'.$sampledata->{$_}{'tag_rc'}, @{$tags};
}

# Stores alignments to write to file at the end
my $output_alignment_data;
push(@{$output_alignment_data}, ['READ_NAME', 'AMPLICON_NAME', 'FW_NAME', 'FW_RANGE', 'READ-FW_PRIMER', 'FW_ERROR', 'RV_NAME', 'RV_RANGE', 'READ-RV_PRIMER', 'RV_ERROR']);
# Initialize counters
my @counter_headers = (	'aligned_reads', 
			'amplicons', 'amplicon_errors_fwd', 'amplicon_errors_rev', 'amplicon_errors_both', 'amplicon_errors_free',
			'amplicon_fwd_qual', 'amplicon_rev_qual', 'amplicon_errors_fwd_qual', 'amplicon_errors_rev_qual', 'amplicon_errors_free_qual',
			'amplicon_aligned_nts', 'amplicon_aligned_nts_errors', 'amplicon_aligned_nts_errors_fwd', 'amplicon_aligned_nts_errors_rev',
			'amplicon_aligned_nts_substitutions_fwd', 'amplicon_aligned_nts_substitutions_rev',
			'amplicon_aligned_nts_insertions_fwd', 'amplicon_aligned_nts_insertions_rev',
			'amplicon_aligned_nts_deletions_fwd', 'amplicon_aligned_nts_deletions_rev',
			'amplicon_aligned_nts_homoindels_fwd', 'amplicon_aligned_nts_homoindels_rev',
			'markers', 'marker_errors', 'marker_errors_free', 'marker_qual', 'marker_errors_qual', 'marker_errors_free_qual',
			'marker_aligned_nts', 'marker_aligned_nts_errors', 'marker_aligned_nts_substitutions', 'marker_aligned_nts_insertions', 'marker_aligned_nts_deletions', 'marker_aligned_nts_homoindels'
			);
my %counters;
map $counters{$_} = 0 , @counter_headers;
# Initialize counter for each primer sequence
my %count_primer_matches;
map $count_primer_matches{$_} = 0 , @{$primer_headers};
my %count_primer_errors;
map $count_primer_errors{$_} = 0 , @{$primer_headers};
# Initialize counter for each tag sequence
my %count_tag_matches;
my %count_tag_errors;
if (defined($tags)) {
	map $count_tag_matches{$_} = 0 , @{$tags};
	map $count_tag_errors{$_} = 0 , @{$tags};
}
# Initialize counters for read and amplicon lengths:
my %total_read_lengths;
my %total_amplicon_lengths;
# Initialize hash to store Sample/Amplicon assignments
my $amplicon_assignments;

for (my $n=1; $n<=$INP_bootstrap; $n++){


	# Select random reads to analyze
	if (!defined($INP_full)){
		print "\nRESAMPLING ROUND $n. Processing $INP_nreads sequences.\n";
	} else {
		print "\tProcessing $INP_nreads sequences.\n";
	}

	my ($read_seqs,$read_headers,$read_qualities);
	($reads_file_format,$read_seqs,$read_headers,$read_qualities,$total_reads)
	= parse_sequence_file($INP_reads_file,$INP_nreads,['qualities']);
# 	my ($read_seqs,$read_headers,$read_qualities);
# 	if (!defined($INP_full)){
# 		if ($reads_file_format eq 'fastq'){
# 			($read_seqs,$read_headers,$read_qualities) = extract_random_seqs_from_fastq($INP_reads_file, $INP_nreads, 1);
# 		} elsif ($reads_file_format eq 'fasta'){
# 			($read_seqs,$read_headers) = read_seqs_from_fasta($INP_reads_file, $INP_nreads);
# 		}
# 	} else {
# 		($read_seqs,$read_headers) = read_fasta_file($INP_reads_file,1);
# 	}

	# Aligns reads against primer/tag sequences
	my ($align_data, $align_primer_data, $align_primer_tag_data);
	if (defined($tags)) {
		print "\tAligning primer+tag sequences.\n";
		$align_primer_tag_data 
		= align_amplicons($read_headers,$read_seqs,$primer_tag_headers,$primer_tag_seqs,$INP_align,$INP_revcomp,$INP_threads);
		$align_data = $align_primer_tag_data;
	} else {
	print "\tAligning primer sequences.\n";
		$align_primer_data 
		= align_amplicons($read_headers,$read_seqs,$primer_headers,$primer_seqs,$INP_align,$INP_revcomp,$INP_threads);
		$align_data = $align_primer_data;
	}

	# Annotates stored alignments to avoid duplications
	my %output_alignment_reads;
	# Loops reads with alignment results
	print "\tExtracting and analyzing amplicon sequences from reads.\n";
	for (my $i=0; $i<=$#{$read_headers}; $i++){
		my $read_header = $read_headers->[$i];
		my $read_seq = $read_seqs->[$i];
		my $read_quality;
		if ($reads_file_format eq 'fastq'){
			$read_quality = $read_qualities->[$i];
		}
		# Counts read length frequency
		$total_read_lengths{length($read_seq)}++;
		if (!defined($align_data->{$read_header})){
			next;
		}
		# Counts reads aligned to any primer sequence
		$counters{'aligned_reads'}++;

		if ($find_amplicons) {
			my ($forward_primer_matches, $reverse_primer_matches);
			my $amplicon_found;
			my ($first_forward_match, $first_reverse_match);
			# Find between primer results matches with common forward and reverse sequences
			foreach my $result (@{$align_data->{$read_header}}) {
	# 			if ($result->{'ALIGNED'} < 14) { next; }
				my $primer_header = $result->{'NAME'};
				$primer_header =~ /(.+)_([F|R])\d+/;
				# If matched sequence is forward primer/tag
				if ($2 eq 'F'){
					# Annotates first match
					if (!defined($first_forward_match)){
						$first_forward_match = $1;
					}
					# Annotate forward sequences matched
					if (!defined($forward_primer_matches->{$1})){
						$forward_primer_matches->{$1} = $result;
					}
					# Stop checking results if the same amplicon has been detected in reverse sequences
					if (defined($reverse_primer_matches) && defined($reverse_primer_matches->{$1})){
						$amplicon_found = $1;
						last;
					}
				# If matched sequence is reverse primer/tag
				} elsif ($2 eq 'R'){
					# Annotates first match
					if (!defined($first_reverse_match)){
						$first_reverse_match = $1;
					}
					# Annotate reverse sequences matched
					if (!defined($reverse_primer_matches->{$1})){
						$reverse_primer_matches->{$1} = $result;
					}
					# Stop checking results if the same amplicon has been detected in forward sequences
					if (defined($forward_primer_matches) && defined($forward_primer_matches->{$1})){
						$amplicon_found = $1;
						last;
					}
				}
			}
			# Calculate statistics and print found primer/tag sequences aligned with the reads
			# If the two primers forward and reverse of the same amplicon are found
			if (defined($amplicon_found)) {
				# Count reads that match a pair of primer/tag forward and reverse sequences
				$counters{'amplicons'}++;
				# Extract names and sequences of the primers, it should be the same name for forward and reverse
				my ($sample_name,$marker_name,$forward_primer_name,$reverse_primer_name);
				if (defined($tags)) {
					$forward_primer_matches->{$amplicon_found}{'NAME'} =~ /(.+)-(.+)_([F|R]\d+)/;
					$sample_name = $1;
					$marker_name = $2;
					$forward_primer_name = "$2\_$3";
					$reverse_primer_matches->{$amplicon_found}{'NAME'} =~ /(.+)-(.+)_([F|R]\d+)/;
					$reverse_primer_name = "$2\_$3";
				} else {
					$forward_primer_matches->{$amplicon_found}{'NAME'} =~ /(.+)_([F|R]\d+)/;
					$marker_name = $1;
					$forward_primer_name = "$1\_$2";
					$reverse_primer_matches->{$amplicon_found}{'NAME'} =~ /(.+)_([F|R]\d+)/;
					$reverse_primer_name = "$1\_$2";
				}
				# Count how many times a primer amplicon is found
				$count_primer_matches{$forward_primer_name}++;
				$count_primer_matches{$reverse_primer_name}++;
				# Annotate aligned sequences of reads and primers
				my ($forward_read_aligned_seq, $forward_primer_aligned_seq) = split("\n",$forward_primer_matches->{$amplicon_found}{'ALIGN'});
				my ($reverse_read_aligned_seq, $reverse_primer_aligned_seq) = split("\n",$reverse_primer_matches->{$amplicon_found}{'ALIGN'});
				# Annotate aligned positions of reads and primers
				my ($forward_read_aligned_cols, $forward_primer_aligned_cols) = split("\n",$forward_primer_matches->{$amplicon_found}{'COLS'});
				my @forward_read_aligned_cols = split(",",$forward_read_aligned_cols);
				my $forward_read_aligned_range = $forward_read_aligned_cols[0]."-".$forward_read_aligned_cols[-1];
				my ($reverse_read_aligned_cols, $reverse_primer_aligned_cols) = split("\n",$reverse_primer_matches->{$amplicon_found}{'COLS'});
				my @reverse_read_aligned_cols = split(",",$reverse_read_aligned_cols);
				my $reverse_read_aligned_range = $reverse_read_aligned_cols[0]."-".$reverse_read_aligned_cols[-1];
				# Annotate mean quality of aligned positions
				my ($forward_read_aligned_qual,$reverse_read_aligned_qual);
				if ($reads_file_format eq 'fastq'){
					$forward_read_aligned_qual = mean(@{convert_phred_ascii_to_numeric(substr($read_quality,$forward_read_aligned_cols[0]-1,$forward_read_aligned_cols[-1]-$forward_read_aligned_cols[0]+1))});
					$reverse_read_aligned_qual = mean(@{convert_phred_ascii_to_numeric(substr($read_quality,$reverse_read_aligned_cols[0]-1,$reverse_read_aligned_cols[-1]-$reverse_read_aligned_cols[0]+1))});
					$counters{'amplicon_fwd_qual'} += $forward_read_aligned_qual;
					$counters{'amplicon_rev_qual'} += $reverse_read_aligned_qual;
				}
				# Count amplicon length frequency
				my $first_forward_primer_pos = $forward_read_aligned_cols[0];
				my $last_reverse_primer_pos = $reverse_read_aligned_cols[-1];
				$total_amplicon_lengths{($last_reverse_primer_pos-$first_forward_primer_pos+1)}++;
				# Counts the number of forward errors
				my $errors_fwd = $forward_primer_matches->{$amplicon_found}{'ALIGNED'} - $forward_primer_matches->{$amplicon_found}{'IDENT'};
				# If forward aligned sequences doesn't match perfectly, annotate errors
				my $forward_primer_errors = '';
				if ($errors_fwd){
					$counters{'amplicon_errors_fwd'}++;
					$count_primer_errors{$forward_primer_name}++;	
					if ($reads_file_format eq 'fastq'){
						$counters{'amplicon_errors_fwd_qual'} += $forward_read_aligned_qual;
					}
					# my ($identity,$total) = binary_score_nts($forward_read_aligned_seq,$forward_primer_aligned_seq)
					# my $errors_fwd = $total - $identity;
					$counters{'amplicon_aligned_nts_errors'} += $errors_fwd;
					$counters{'amplicon_aligned_nts_errors_fwd'} += $errors_fwd;
					# Counts the type of errors
					my ($substitutiones, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($forward_read_aligned_seq,$forward_primer_aligned_seq);
					$counters{'amplicon_aligned_nts_substitutions_fwd'} += scalar @$substitutiones;
					$counters{'amplicon_aligned_nts_insertions_fwd'} += scalar @$insertions;
					$counters{'amplicon_aligned_nts_deletions_fwd'} += scalar @$deletions;
					$counters{'amplicon_aligned_nts_homoindels_fwd'} += scalar @$homopolymer_indels;
					$forward_primer_errors = print_sequence_errors($substitutiones, $insertions, $deletions, $homopolymer_indels);
				}
				# Counts the number of reverse errors
				my $errors_rev = $reverse_primer_matches->{$amplicon_found}{'ALIGNED'} - $reverse_primer_matches->{$amplicon_found}{'IDENT'};
				# If reverse aligned sequences doesn't match perfectly, annotate errors
				my $reverse_primer_errors = '';
				if ($errors_rev){
					$counters{'amplicon_errors_rev'}++;
					$count_primer_errors{$reverse_primer_name}++;
					if ($reads_file_format eq 'fastq'){
						$counters{'amplicon_errors_rev_qual'} += $reverse_read_aligned_qual;
					}
					$counters{'amplicon_aligned_nts_errors'} += $errors_rev;
					$counters{'amplicon_aligned_nts_errors_rev'} += $errors_rev;
					# Counts the type of errors
					my ($substitutiones, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($reverse_read_aligned_seq,$reverse_primer_aligned_seq);
					$counters{'amplicon_aligned_nts_substitutions_rev'} += scalar @$substitutiones;
					$counters{'amplicon_aligned_nts_insertions_rev'} += scalar @$insertions;
					$counters{'amplicon_aligned_nts_deletions_rev'} += scalar @$deletions;
					$counters{'amplicon_aligned_nts_homoindels_rev'} += scalar @$homopolymer_indels;
					$reverse_primer_errors = print_sequence_errors($substitutiones, $insertions, $deletions, $homopolymer_indels);
				}
				# Count when both forward and reverse sequences have or are free of errors
				if ($errors_fwd && $errors_rev){
					$counters{'amplicon_errors_both'}++;
					if ($reads_file_format eq 'fastq'){
						$counters{'both_errors_qual'} += ($forward_read_aligned_qual+$reverse_read_aligned_qual)/2;
					}
				} elsif (!$errors_fwd && !$errors_rev){
					$counters{'amplicon_errors_free'}++;
	# 				$exact_total_amplicon_lengths{($last_reverse_primer_pos-$first_forward_primer_pos+1)}++;
					if ($reads_file_format eq 'fastq'){
						$counters{'amplicon_errors_free_qual'} += ($forward_read_aligned_qual+$reverse_read_aligned_qual)/2;
					}
				}
	# 	 		if ($read_header eq '8GJMT:00443:01178'){
	# 				print '';
	# 			}
				# Counts the number of total nucleotides
				$counters{'amplicon_aligned_nts'} += $forward_primer_matches->{$amplicon_found}{'ALIGNED'};
				$counters{'amplicon_aligned_nts'} += $reverse_primer_matches->{$amplicon_found}{'ALIGNED'};

				# Annotates reads alignment data (name of reads, matching amplicons and primer sequences aligned with the reads
				if (defined($INP_verbose)){
					my $forward_output_seq = alignment_to_single_line($forward_primer_matches->{$amplicon_found}{'ALIGN'});
					my $reverse_output_seq = alignment_to_single_line($reverse_primer_matches->{$amplicon_found}{'ALIGN'});
					push(@{$output_alignment_data}, [$read_header, $amplicon_found, $forward_primer_name, $forward_read_aligned_range, $forward_output_seq, $forward_primer_errors, $reverse_primer_name, $reverse_read_aligned_range, $reverse_output_seq, $reverse_primer_errors]);
				}
			} else {
				$counters{'no_amplicon'}++;
				if (!defined($forward_primer_matches) && !defined($reverse_primer_matches)){
					$counters{'no_primers'}++;
				} elsif (defined($forward_primer_matches) || defined($reverse_primer_matches)){
					$counters{'one_primer'}++;
				}
				# Annotates empty alignment data of non aligned reads
				if (defined($INP_verbose)){
					if (!defined($forward_primer_matches) && !defined($reverse_primer_matches)){
						push(@{$output_alignment_data}, [$read_header]);
					} elsif (defined($forward_primer_matches) && defined($reverse_primer_matches)){
						my $forward_primer_name = $forward_primer_matches->{$first_forward_match}{'NAME'};
						my $reverse_primer_name = $reverse_primer_matches->{$first_reverse_match}{'NAME'};
						my ($forward_read_aligned_seq, $forward_primer_aligned_seq) = split("\n",$forward_primer_matches->{$first_forward_match}{'ALIGN'});
						my ($reverse_read_aligned_seq, $reverse_primer_aligned_seq) = split("\n",$reverse_primer_matches->{$first_reverse_match}{'ALIGN'});
						my ($forward_read_aligned_cols, $forward_primer_aligned_cols) = split("\n",$forward_primer_matches->{$first_forward_match}{'COLS'});
						my @forward_read_aligned_cols = split(",",$forward_read_aligned_cols);
						my $forward_read_aligned_range = $forward_read_aligned_cols[0]."-".$forward_read_aligned_cols[-1];
						my ($reverse_read_aligned_cols, $reverse_primer_aligned_cols) = split("\n",$reverse_primer_matches->{$first_reverse_match}{'COLS'});
						my @reverse_read_aligned_cols = split(",",$reverse_read_aligned_cols);
						my $reverse_read_aligned_range = $reverse_read_aligned_cols[0]."-".$reverse_read_aligned_cols[-1];
						my $forward_primer_errors = '';
						if ($forward_read_aligned_seq ne $forward_primer_aligned_seq){
							my ($substitutiones, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($forward_read_aligned_seq,$forward_primer_aligned_seq);
							$forward_primer_errors = print_sequence_errors($substitutiones, $insertions, $deletions, $homopolymer_indels);
						}
						my $reverse_primer_errors = '';
						if ($reverse_read_aligned_seq ne $reverse_primer_aligned_seq){
							my ($substitutiones, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($reverse_read_aligned_seq,$reverse_primer_aligned_seq);
							$reverse_primer_errors = print_sequence_errors($substitutiones, $insertions, $deletions, $homopolymer_indels);
						}
						my $forward_output_seq = alignment_to_single_line($forward_primer_matches->{$first_forward_match}{'ALIGN'});
						my $reverse_output_seq = alignment_to_single_line($reverse_primer_matches->{$first_reverse_match}{'ALIGN'});
						push(@{$output_alignment_data}, [$read_header, '', $forward_primer_name, $forward_read_aligned_range, $forward_output_seq, $forward_primer_errors, $reverse_primer_name, $reverse_read_aligned_range, $reverse_output_seq, $reverse_primer_errors]);
					} elsif (defined($forward_primer_matches)){
						my $forward_primer_name = $forward_primer_matches->{$first_forward_match}{'NAME'};
						my ($forward_read_aligned_seq, $forward_primer_aligned_seq) = split("\n",$forward_primer_matches->{$first_forward_match}{'ALIGN'});
						my ($forward_read_aligned_cols, $forward_primer_aligned_cols) = split("\n",$forward_primer_matches->{$first_forward_match}{'COLS'});
						my @forward_read_aligned_cols = split(",",$forward_read_aligned_cols);
						my $forward_read_aligned_range = $forward_read_aligned_cols[0]."-".$forward_read_aligned_cols[-1];
						my $forward_primer_errors = '';
						if ($forward_read_aligned_seq ne $forward_primer_aligned_seq){
							my ($substitutiones, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($forward_read_aligned_seq,$forward_primer_aligned_seq);
							$forward_primer_errors = print_sequence_errors($substitutiones, $insertions, $deletions, $homopolymer_indels);
						}
						my $forward_output_seq = alignment_to_single_line($forward_primer_matches->{$first_forward_match}{'ALIGN'});
						push(@{$output_alignment_data}, [$read_header, '', $forward_primer_name, $forward_read_aligned_range, $forward_output_seq, $forward_primer_errors]);
					} elsif (defined($reverse_primer_matches)){
						my $reverse_primer_name = $reverse_primer_matches->{$first_reverse_match}{'NAME'};
						my ($reverse_read_aligned_seq, $reverse_primer_aligned_seq) = split("\n",$reverse_primer_matches->{$first_reverse_match}{'ALIGN'});
						my ($reverse_read_aligned_cols, $reverse_primer_aligned_cols) = split("\n",$reverse_primer_matches->{$first_reverse_match}{'COLS'});
						my @reverse_read_aligned_cols = split(",",$reverse_read_aligned_cols);
						my $reverse_read_aligned_range = $reverse_read_aligned_cols[0]."-".$reverse_read_aligned_cols[-1];
						my $reverse_primer_errors = '';
						if ($reverse_read_aligned_seq ne $reverse_primer_aligned_seq){
							my ($substitutiones, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($reverse_read_aligned_seq,$reverse_primer_aligned_seq);
							$reverse_primer_errors = print_sequence_errors($substitutiones, $insertions, $deletions, $homopolymer_indels);
						}
						my $reverse_output_seq = alignment_to_single_line($reverse_primer_matches->{$first_reverse_match}{'ALIGN'});
						push(@{$output_alignment_data}, [$read_header, '', '', '', '', '', $reverse_primer_name, $reverse_read_aligned_range, $reverse_output_seq, $reverse_primer_errors]);
					}
				}
			}
		# If there are not primers/tags forward and reverse (only one side)
		} elsif (!$find_amplicons && defined($align_data->{$read_header})) {
			# Annotates the statistics for the (tag)+marker alignments (if they exist)
			my $result = $align_data->{$read_header}[0];
			# Count reads that match a primer/tag sequence
			$counters{'markers'}++;
			# Extract names and sequences of the primers
			my ($sample_name,$marker_name,$primer_name,$match_name);
			if (defined($tags)) {
				$result->{'NAME'} =~ /(.+)-(.+)_([F|R]\d+)/;
				$sample_name = $1;
				$marker_name = $2;
				$primer_name = "$2\_$3";
				$match_name = "$1-$2";
			} elsif ($result->{'NAME'} =~ /(.+)_([F|R]\d+)/) {
				$marker_name = $1;
				$primer_name = "$1\_$2";
				$match_name = $1;
			} else {
				$primer_name = $result->{'NAME'};
				$match_name = $result->{'NAME'};
			}
			# Count how many times a primer amplicon is found
			$count_primer_matches{$primer_name}++;
			# Annotate aligned sequences of reads and primers
			my ($read_aligned_seq, $primer_aligned_seq) = split("\n",$result->{'ALIGN'});
			# Annotate aligned positions of reads and primers
			my ($read_aligned_cols, $primer_aligned_cols) = split("\n",$result->{'COLS'});
			my @read_aligned_cols = split(",",$read_aligned_cols);
			my $read_aligned_range = $read_aligned_cols[0]."-".$read_aligned_cols[-1];
			# Annotate mean quality of aligned positions
			my $read_aligned_qual;
			if ($reads_file_format eq 'fastq'){
				$read_aligned_qual = mean(convert_phred_ascii_to_numeric(substr($read_quality,$read_aligned_cols[0]-1,$read_aligned_cols[-1]-$read_aligned_cols[0]+1)));
				$counters{'marker_qual'} += $read_aligned_qual;
			}
			# Counts the number of errors
			my $errors = $result->{'ALIGNED'} - $result->{'IDENT'};
			# If aligned sequences doesn't match perfectly, annotate errors
			my $primer_errors = '';
			if ($errors){
				$counters{'marker_errors'}++;
				$count_primer_errors{$primer_name}++;	
				if ($reads_file_format eq 'fastq'){
					$counters{'marker_errors_qual'} += $read_aligned_qual;
				}
				# my ($identity,$total) = binary_score_nts($read_aligned_seq,$primer_aligned_seq)
				# my $errors_fwd = $total - $identity;
				$counters{'marker_aligned_nts_errors'} += $errors;
				# Counts the type of errors
				my ($substitutiones, $insertions, $deletions, $homopolymer_indels) = detect_sequence_errors($read_aligned_seq,$primer_aligned_seq);
				$counters{'marker_aligned_nts_substitutions'} += scalar @$substitutiones;
				$counters{'marker_aligned_nts_insertions'} += scalar @$insertions;
				$counters{'marker_aligned_nts_deletions'} += scalar @$deletions;
				$counters{'marker_aligned_nts_homoindels'} += scalar @$homopolymer_indels;
				$primer_errors = print_sequence_errors($substitutiones, $insertions, $deletions, $homopolymer_indels);
			} else {
				$counters{'marker_errors_free'}++;
# 				$exact_total_amplicon_lengths{($last_reverse_primer_pos-$first_forward_primer_pos+1)}++;
				if ($reads_file_format eq 'fastq'){
					$counters{'marker_errors_free_qual'} += $read_aligned_qual;
				}
			}
			# Counts the number of total nucleotides
			$counters{'marker_aligned_nts'} += $result->{'ALIGNED'};
			# Annotates reads alignment data (name of reads, matching amplicons and primer sequences aligned with the reads
			if (defined($INP_verbose)){
				my $output_seq = alignment_to_single_line($result->{'ALIGN'});
				push(@{$output_alignment_data}, [$read_header, $match_name, $primer_name, $read_aligned_range, $output_seq, $primer_errors]);
			}
		}
		if (defined($output_alignment_reads{$read_header})){
			print "\tERROR: sequence name '$read_header' is duplicated.\n";
		} else {
			$output_alignment_reads{$read_header} = 1;
		}
		if (defined($tags)) {
			# Find between primer+tag results matches with common forward and reverse sequences
			my ($forward_tag_matches, $reverse_tag_matches) = (undef,undef);
			my $sample_found;
			foreach my $result (@{$align_primer_tag_data->{$read_header}}) {
				my $primer_tag_header = $result->{'NAME'};
				$primer_tag_header =~ /(.+)_([F|R])\d+/;
				# If matched sequence is forward primer/tag
				if ($2 eq 'F'){
					# Annotate forward sequences matched
					if (!defined($forward_tag_matches->{$1})){
						$forward_tag_matches->{$1} = $result;
					}
					# Stop checking results if the same primer_tag has been detected in reverse sequences
					if (defined($reverse_tag_matches) && defined($reverse_tag_matches->{$1})){
						$sample_found = $1;
						last;
					}
				# If matched sequence is reverse primer/tag
				} elsif ($2 eq 'R'){
					# Annotate reverse sequences matched
					if (!defined($reverse_tag_matches->{$1})){
						$reverse_tag_matches->{$1} = $result;
					}
					# Stop checking results if the same primer_tag has been detected in forward sequences
					if (defined($forward_tag_matches) && defined($forward_tag_matches->{$1})){
						$sample_found = $1;
						last;
					}
				}
			}
			# Annotates the statistics for double tags
			if (defined($sample_found)) {
				# Extract names and sequences of the primers+tags, it should be the same name for forward and reverse
				$forward_tag_matches->{$sample_found}{'NAME'} =~ /(.+)-(.+)_([F|R])\d+/;
				my $sample_name = $1;
				my $marker_name = $2;
				# Count how many times a sample is found
				$count_tag_matches{$sample_name}++;
				# Count how many times a sample/amplicon assignment is found
				$amplicon_assignments->{$sample_name}{$marker_name}++;
				# Annotate aligned sequences of reads and primers/tags
				my ($forward_read_aligned_seq, $forward_sample_aligned_seq) = split("\n",$forward_tag_matches->{$sample_found}{'ALIGN'});
				my ($reverse_read_aligned_seq, $reverse_sample_aligned_seq) = split("\n",$reverse_tag_matches->{$sample_found}{'ALIGN'});
				if ($forward_read_aligned_seq ne $forward_sample_aligned_seq || $reverse_read_aligned_seq ne $reverse_sample_aligned_seq){
					$count_tag_errors{$sample_name}++;
				}
			} else {
				# Annotates the statistics for single tags
				my $tag_matches = {};
				if (defined($reverse_tag_matches) && defined($forward_tag_matches)){
					$tag_matches = { %$forward_tag_matches, %$reverse_tag_matches };
				} elsif (defined($forward_tag_matches)){
					$tag_matches = $forward_tag_matches;
				} elsif (defined($reverse_tag_matches)){
					$tag_matches = $reverse_tag_matches;
				}
				foreach my $match (keys %$tag_matches){
					# Extract names and sequences of the primers+tags, it should be the same name for forward and reverse
					$tag_matches->{$match}{'NAME'} =~ /(.+)-(.+)_([F|R])\d+/;
					my $sample_name = $1;
					my $marker_name = $2;
					# Count how many times a sample is found
					$count_tag_matches{$sample_name}++;
					# Count how many times a sample/amplicon assignment is found
					$amplicon_assignments->{$sample_name}{$marker_name}++;
					# Annotate aligned sequences of reads and primers/tags
					my ($read_aligned_seq, $sample_aligned_seq) = split("\n",$tag_matches->{$match}{'ALIGN'});
					if ($read_aligned_seq ne $sample_aligned_seq){
						$count_tag_errors{$sample_name}++;
					}

				}
			}
		}
	}
}

# Prints statistics
print "\nQuality report printed in file '$INP_outfile.xlsx'.\n";
my $workbook  = Excel::Writer::XLSX->new("$INP_outfile.xlsx");
$workbook->set_properties(
	title    => 'AmpliQC',
	author   => 'Alvaro Sebastian',
	comments => 'AmpliQC software report',
);
my %font_normal = (
	font  => 'Arial',
	size  => 9,
	color => 'black',
	bold  => 0,
);
my %font_bold = (
	font  => 'Arial',
	size  => 10,
	color => 'black',
	bold  => 1,
);
my %font_header = (
	font  => 'Arial',
	size  => 12,
	color => 'black',
	bold  => 1,
);
my $format_normal = $workbook->add_format( %font_normal );
my $format_bold = $workbook->add_format( %font_bold );
my $format_header = $workbook->add_format( %font_header );
my $format_percentage = $workbook->add_format( %font_normal );
$format_percentage->set_num_format( '0.00%' );
my $format_decimals = $workbook->add_format( %font_normal );
$format_decimals->set_num_format( '[=0]0;.##' );
my $format_nodecimals = $workbook->add_format( %font_normal );
$format_nodecimals->set_num_format( '#0' );
my $format_right = $workbook->add_format( %font_normal );
$format_right->set_align( 'right' );
my $format_right_bold = $workbook->add_format( %font_bold );
$format_right_bold->set_align( 'right' );
my $format_center = $workbook->add_format( %font_normal );
$format_center->set_align( 'center' );
my $format_center_bold = $workbook->add_format( %font_bold );
$format_center_bold->set_align( 'center' );

# Write sheet with quality statistics
my $worksheet = $workbook->add_worksheet('quality');
# Apply column formats
$worksheet->set_column( 'A:A', 40, $format_normal);
$worksheet->set_column( 'B:B', undef, $format_nodecimals);
$worksheet->set_column( 'C:C', undef, $format_percentage);
$worksheet->set_column( 'D:D', undef, $format_normal );
# Multiply results by scaling factor
my $scaling_factor = $total_reads / ($INP_nreads*$INP_bootstrap);
my $row = 1;
$worksheet->write( "A$row", 'QUALITY ANALYSIS STATISTICS (AmpliQC v1.0)', $format_header); $row++;
$row++;
$worksheet->write( "A$row", 'INPUT PARAMETERS:', $format_header); $row++;
$row++;
$worksheet->write( "A$row", 'Input file'); $worksheet->write( "B$row", $INP_reads_file); $row++;
$worksheet->write( "A$row", 'CSV data file'); $worksheet->write( "B$row", $INP_amplicons_file); $row++;
$worksheet->write( "A$row", 'Aligment type'); $worksheet->write( "B$row", $INP_align); $row++;
my $totalreads_pos = "B$row";
$worksheet->write( "A$row", 'Total reads in input file', $format_bold); $worksheet->write( "B$row", $total_reads, $format_bold); $row++;
my $nreads_pos = "B$row";
$worksheet->write( "A$row", 'Number of reads analyzed per round'); $worksheet->write( "B$row", $INP_nreads); $row++;
my $bootstrap_pos = "B$row";
$worksheet->write( "A$row", 'Number of resampling rounds'); $worksheet->write( "B$row", $INP_bootstrap); $row++;
$row++;
if (!defined($INP_full)){
	$worksheet->write( "A$row", 'EXTRAPOLATED RESULTS AFTER ANALYSIS:', $format_header); $row++;
} else {
	$worksheet->write( "A$row", 'RESULTS AFTER ANALYSIS:', $format_header); $row++;
}
$row++;
if ($counters{'amplicons'}>0){
	$worksheet->write( "A$row", 'Alignment of reads with primers:', $format_bold); $row++;
	$worksheet->write( "A$row", 'Reads aligned with at least one primer'); $worksheet->write( "B$row", sprintf("%.0f", $counters{'aligned_reads'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalreads_pos"); $row++;
	my $matchampli_pos = "B$row";
	$worksheet->write( "A$row", 'Reads matching both primers'); $worksheet->write( "B$row", sprintf("%.0f", $counters{'amplicons'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalreads_pos"); $row++;
	$worksheet->write( "A$row", 'With errors in forward primer', $format_right); $worksheet->write( "B$row", sprintf("%.0f", $counters{'amplicon_errors_fwd'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$matchampli_pos"); $row++;
	$worksheet->write( "A$row", 'With errors in reverse primer', $format_right); $worksheet->write( "B$row", sprintf("%.0f", $counters{'amplicon_errors_rev'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$matchampli_pos"); $row++;
	$worksheet->write( "A$row", 'With errors in both primers', $format_right); $worksheet->write( "B$row", sprintf("%.0f", $counters{'amplicon_errors_both'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$matchampli_pos"); $row++;
	$worksheet->write( "A$row", 'Without errors', $format_right); $worksheet->write( "B$row", sprintf("%.0f", $counters{'amplicon_errors_free'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$matchampli_pos"); $row++;
	$row++;
	$worksheet->write( "A$row", 'Errors in aligned nucleotides:', $format_bold); $row++;
	my $totalnts_pos = "B$row";
	$worksheet->write( "A$row", 'Nucleotides aligned with primers'); $worksheet->write( "B$row", sprintf("%.0f", $counters{'amplicon_aligned_nts'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Nucleotides with errors' ); $worksheet->write( "B$row", sprintf("%.0f", $counters{'amplicon_aligned_nts_errors'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Nucleotides with errors in forward primer'); $worksheet->write( "B$row", sprintf("%.0f", $counters{'amplicon_aligned_nts_errors_fwd'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Nucleotides with errors in reverse primer'); $worksheet->write( "B$row", sprintf("%.0f", $counters{'amplicon_aligned_nts_errors_rev'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Substitutions', $format_right); $worksheet->write( "B$row", sprintf("%.0f", ($counters{'amplicon_aligned_nts_substitutions_fwd'}+$counters{'amplicon_aligned_nts_substitutions_rev'})*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Insertions', $format_right); $worksheet->write( "B$row", sprintf("%.0f", ($counters{'amplicon_aligned_nts_insertions_fwd'}+$counters{'amplicon_aligned_nts_insertions_rev'})*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Deletions', $format_right); $worksheet->write( "B$row", sprintf("%.0f", ($counters{'amplicon_aligned_nts_deletions_fwd'}+$counters{'amplicon_aligned_nts_deletions_rev'})*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Homopolymer indels', $format_right); $worksheet->write( "B$row", sprintf("%.0f", ($counters{'amplicon_aligned_nts_homoindels_fwd'}+$counters{'amplicon_aligned_nts_homoindels_rev'})*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$row++;
	if ($reads_file_format eq 'fastq'){
		$worksheet->write( "A$row", 'Average PHRED quality of aligned reads:', $format_bold); $row++;
		if ($counters{'amplicons'}>0){
			$worksheet->write( "A$row", 'Quality of aligned reads in forward region'); $worksheet->write( "B$row", ($counters{'amplicon_fwd_qual'}/$counters{'amplicons'}), $format_decimals); $row++;
			$worksheet->write( "A$row", 'Quality of aligned reads in reverse region'); $worksheet->write( "B$row", ($counters{'amplicon_rev_qual'}/$counters{'amplicons'}), $format_decimals); $row++;
		}
		if ($counters{'amplicon_errors_fwd'}>0){
			$worksheet->write( "A$row", 'Quality of aligned reads with errors in forward region'); $worksheet->write( "B$row", ($counters{'amplicon_errors_fwd_qual'}/$counters{'amplicon_errors_fwd'}), $format_decimals); $row++;
		}
		if ($counters{'amplicon_errors_rev'}>0){
			$worksheet->write( "A$row", 'Quality of aligned reads with errors in reverse region'); $worksheet->write( "B$row", ($counters{'amplicon_errors_rev_qual'}/$counters{'amplicon_errors_rev'}), $format_decimals); $row++;
		}
		if ($counters{'amplicon_errors_free'}>0){
			$worksheet->write( "A$row", 'Quality of aligned reads without errors'); $worksheet->write( "B$row", ($counters{'amplicon_errors_free_qual'}/$counters{'amplicon_errors_free'}), $format_decimals); $row++;
		}
		$row++;
	}
} else {
	$worksheet->write( "A$row", 'Alignment of reads with markers:', $format_bold); $row++;
	$worksheet->write( "A$row", 'Reads aligned with at least one marker'); $worksheet->write( "B$row", sprintf("%.0f", $counters{'aligned_reads'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalreads_pos"); $row++;
	my $matchampli_pos = "B$row";
	$worksheet->write( "A$row", 'Reads matching a marker'); $worksheet->write( "B$row", sprintf("%.0f", $counters{'markers'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalreads_pos"); $row++;
	$worksheet->write( "A$row", 'With errors', $format_right); $worksheet->write( "B$row", sprintf("%.0f", $counters{'marker_errors'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$matchampli_pos"); $row++;
	$worksheet->write( "A$row", 'Without errors', $format_right); $worksheet->write( "B$row", sprintf("%.0f", $counters{'marker_errors_free'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$matchampli_pos"); $row++;
	$row++;
	$worksheet->write( "A$row", 'Errors in aligned nucleotides:', $format_bold); $row++;
	my $totalnts_pos = "B$row";
	$worksheet->write( "A$row", 'Nucleotides aligned with markers'); $worksheet->write( "B$row", sprintf("%.0f", $counters{'marker_aligned_nts'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Nucleotides with errors' ); $worksheet->write( "B$row", sprintf("%.0f", $counters{'marker_aligned_nts_errors'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Substitutions', $format_right); $worksheet->write( "B$row", sprintf("%.0f", $counters{'marker_aligned_nts_substitutions'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Insertions', $format_right); $worksheet->write( "B$row", sprintf("%.0f", $counters{'marker_aligned_nts_insertions'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Deletions', $format_right); $worksheet->write( "B$row", sprintf("%.0f", $counters{'marker_aligned_nts_deletions'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$worksheet->write( "A$row", 'Homopolymer indels', $format_right); $worksheet->write( "B$row", sprintf("%.0f", $counters{'marker_aligned_nts_homoindels'}*$scaling_factor)); $worksheet->write( "C$row", "=B$row/$totalnts_pos"); $row++;
	$row++;
	if ($reads_file_format eq 'fastq'){
		$worksheet->write( "A$row", 'Average PHRED quality of aligned reads:', $format_bold); $row++;
		if ($counters{'markers'}>0){
			$worksheet->write( "A$row", 'Quality of aligned reads'); $worksheet->write( "B$row", ($counters{'marker_qual'}/$counters{'markers'}), $format_decimals); $row++;
		}
		if ($counters{'marker_errors'}>0){
			$worksheet->write( "A$row", 'Quality of aligned reads with errors'); $worksheet->write( "B$row", ($counters{'marker_errors_qual'}/$counters{'marker_errors'}), $format_decimals); $row++;
		}
		if ($counters{'marker_errors_free'}>0){
			$worksheet->write( "A$row", 'Quality of aligned reads without errors'); $worksheet->write( "B$row", ($counters{'marker_errors_free_qual'}/$counters{'marker_errors_free'}), $format_decimals); $row++;
		}
		$row++;
	}

}
# $worksheet->write( "A$row", 'Results primer by primer:', $format_bold);
$row++;
$worksheet->write( "A$row", 'Primer', $format_bold); $worksheet->write( "B$row", 'Reads', $format_right_bold); $worksheet->write( "C$row", 'Errors', $format_right_bold); $worksheet->write( "D$row", '% Errors', $format_right_bold); $row++;
foreach my $primer_header (@{$primer_headers}){
	$worksheet->write( "A$row", sprintf("%s (%s)", $primer_header, $primer_seqs_hash{$primer_header}));
	$worksheet->write( "B$row", sprintf("%.0f", $count_primer_matches{$primer_header}*$scaling_factor), $format_nodecimals);
	$worksheet->write( "C$row", sprintf("%.0f", $count_primer_errors{$primer_header}*$scaling_factor), $format_nodecimals);
	if ($count_primer_matches{$primer_header}){
		$worksheet->write( "D$row", "=C$row/B$row", $format_percentage);
	}
	$row++;
}
if (defined($tags)) {
	$row++;
	$worksheet->write( "A$row", 'Sample', $format_bold); $worksheet->write( "B$row", 'Reads', $format_right_bold); $worksheet->write( "C$row", 'Errors', $format_right_bold); $worksheet->write( "D$row", '% Errors', $format_right_bold); $row++;
	foreach my $sample_name (@{$tags}){
		$worksheet->write( "A$row", sprintf("%s (%s)", $sample_name, $tag_seqs_hash{$sample_name}));
		$worksheet->write( "B$row", sprintf("%.0f", $count_tag_matches{$sample_name}*$scaling_factor), $format_nodecimals);
		$worksheet->write( "C$row", sprintf("%.0f", $count_tag_errors{$sample_name}*$scaling_factor), $format_nodecimals);
		if ($count_tag_matches{$sample_name}){
			$worksheet->write( "D$row", "=C$row/B$row", $format_percentage);
		}
		$row++;
	}
}

# Write sheet with read length statistics
$worksheet = $workbook->add_worksheet('lengths');
$worksheet->set_column( 'A:A', undef, $format_nodecimals);
$worksheet->set_column( 'B:B', undef, $format_nodecimals);
$worksheet->set_column( 'C:C', undef, $format_nodecimals);
$row = 1;
if (!defined($INP_full)){
	$worksheet->write( "A$row", 'EXTRAPOLATED READ/AMPLICON LENGTH STATISTICS (AmpliQC v1.0)', $format_header); $row++;
} else {
	$worksheet->write( "A$row", 'READ/AMPLICON LENGTH STATISTICS (AmpliQC v1.0)', $format_header); $row++;
}
$row++;
$worksheet->write( "A$row", 'Number of sequences:', $format_bold); $row++;
my $values_row = $row;
if ($counters{'amplicons'}>0){
	$worksheet->write( "A$row", 'Length (bps)', $format_right_bold); $worksheet->write( "B$row", 'Reads', $format_right_bold); $worksheet->write( "C$row", 'Amplicons', $format_right_bold); $row++;
} else {
	$worksheet->write( "A$row", 'Length (bps)', $format_right_bold); $worksheet->write( "B$row", 'Reads', $format_right_bold); $row++;
}
my @total_read_lengths = sort {$a<=>$b} keys %total_read_lengths;
my $min_length = $total_read_lengths[0];
my $max_length = $total_read_lengths[-1];
my $first_row = 5;
foreach my $length ($min_length .. $max_length){
	$worksheet->write( "A$row", $length);
	if (defined($total_read_lengths{$length})){
		$worksheet->write( "B$row", sprintf("%.0f", $total_read_lengths{$length}*$scaling_factor));
	} else {
		$worksheet->write( "B$row", 0);
	}
	if (defined($total_amplicon_lengths{$length})){
		$worksheet->write( "C$row", sprintf("%.0f", $total_amplicon_lengths{$length}*$scaling_factor));
	} else {
		$worksheet->write( "C$row", 0);
	}
	$row++;
}
my $last_row = $row-1;
# Create a new chart object. In this case an embedded chart.
my $chart = $workbook->add_chart( type => 'column', embedded => 1 );
# Add a chart title and some axis labels.
$chart->set_title ( name => 'Sequence length frequency' );
$chart->set_x_axis( name => 'Length (bps)' );
# $chart->set_y_axis( name => 'Number of sequences' );
# Set an Excel chart style. Colors with white outline and shadow.
$chart->set_style( 2 );
# Configure the series.
$chart->add_series(
	name       => "=lengths!\$B\$$values_row",
	categories => "=lengths!A$first_row:A$last_row",
	values => "=lengths!B$first_row:B$last_row",
);
$chart->add_series(
	name       => "=lengths!\$C\$$values_row",
	categories => "=lengths!A$first_row:A$last_row",
	values => "=lengths!C$first_row:C$last_row",
);
# Insert the chart into the worksheet (with an offset).
$worksheet->insert_chart( 'E4', $chart, 0, 0, 1.5, 1.5 );

# Write sheet with sample statistics
if (defined($tags)) {
	$worksheet = $workbook->add_worksheet('amplicons');
	$row = 0;
	if (!defined($INP_full)){
		$worksheet->write($row, 0 , 'EXTRAPOLATED AMPLICON ASSIGNMENTS (AmpliQC v1.0)', $format_header); $row++;
	} else {
		$worksheet->write($row, 0 , 'AMPLICON ASSIGNMENTS (AmpliQC v1.0)', $format_header); $row++;
	}
	$row++;
	for (my $i=0; $i<=$#{$primers}; $i++){
		$worksheet->write($row, $i+1, $primers->[$i], $format_right_bold);
	}
	$worksheet->write($row, $#{$primers}+2, "TOTAL", $format_right_bold);
	$row++;
	for (my $i=0; $i<=$#{$tags}; $i++){
		$worksheet->write($row, 0, $tags->[$i], $format_right_bold);
		for (my $j=0; $j<=$#{$primers}; $j++){
			if (defined($amplicon_assignments->{$tags->[$i]}{$primers->[$j]})){
				$worksheet->write($row, $j+1, sprintf("%.0f", $amplicon_assignments->{$tags->[$i]}{$primers->[$j]}*$scaling_factor), $format_normal);
			} else {
				$worksheet->write($row, $j+1, sprintf("%.0f", 0), $format_normal);
			}
		}
		$worksheet->write($row, $#{$primers}+2, "=SUM(B".($row+1).":".xl_col_to_name($#{$primers}+1).($row+1).")", $format_right_bold);
		$row++;
	}
	$worksheet->write($row, 0, "TOTAL", $format_right_bold);
	for (my $j=0; $j<=$#{$primers}+1; $j++){
		$worksheet->write($row, $j+1, "=SUM(".xl_col_to_name($j+1)."4:".xl_col_to_name($j+1).$row.")", $format_right_bold);
	}

}

if (defined($INP_verbose)){
# DON'T WRITE ALIGNMENTS IN THE EXCEL FILE BECAUSE WILL BE VERY BIG AND DIFFICULT TO OPEN
# 	# Write sheet with sequence alignments and errors
# 	$worksheet = $workbook->add_worksheet('alignments');
# 	$worksheet->set_column( 'A:A', 20, $format_center);
# 	$worksheet->set_column( 'B:B', 20, $format_center);
# 	$worksheet->set_column( 'C:C', 10, $format_center);
# 	$worksheet->set_column( 'D:D', 10, $format_center);
# 	$worksheet->set_column( 'E:E', 30, $format_center);
# 	$worksheet->set_column( 'F:F', 10, $format_center);
# 	$worksheet->set_column( 'G:G', 10, $format_center);
# 	$worksheet->set_column( 'H:H', 10, $format_center);
# 	$worksheet->set_column( 'I:I', 30, $format_center);
# 	$worksheet->set_column( 'J:J', 10, $format_center);
# 	# $worksheet->set_column( 'B:J', undef, $format_center);
# 	$row = 1;
# 	$worksheet->write( "A$row", 'SEQUENCE/PRIMER ALIGNMENTS (AmpliQC v1.0)', $format_header); $row++;
# 	for (my $col_=0; $col_<=$#{$output_alignment_data->[0]}; $col_++) {
# 		for (my $row_=0; $row_<=$#{$output_alignment_data}; $row_++){
# 			if ($row_ != 0){
# 				$worksheet->write($row_+$row, $col_, $output_alignment_data->[$row_][$col_] );
# 			} else {
# 				$worksheet->write($row_+$row, $col_, $output_alignment_data->[$row_][$col_], $format_center_bold );
# 			}
# 		}
# 	}

	print "\nPrinting sequence/primer alignments into file '$INP_outfile.seqs.txt'.\n";
	open(OUTFILE, ">$INP_outfile.seqs.txt");
	print OUTFILE "SEQUENCE/PRIMER ALIGNMENTS (AmpliQC v1.0)\n\n";
	for (my $row_=0; $row_<=$#{$output_alignment_data}; $row_++){
		print OUTFILE join("\t", @{$output_alignment_data->[$row_]})."\n";
	}
	close OUTFILE;
}

$workbook->close();

print "\n";




