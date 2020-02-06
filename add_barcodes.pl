#!/usr/bin/perl -w
#
################################################################
#
# Name: add_tags.pl
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
#   Add tags to a list of files taking as input an amplicon data file with primers
#
# Examples:
# perl add_tags.pl -i seqs_file.fa -d primers.csv -fl 6 -rl 6
# perl add_tags.pl -i fastq_folder -d primers.csv -fl 6

# Modules are in folder 'lib' in the path of the script
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use File::Basename;
use File::Find qw(find);
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options
# Default alignment algorithm
my $INP_align = 'match';
# By default align direct and reverse complementary reads
my $INP_revcomp = 1;

my (@INP_input,@INP_amplicons_files,$INP_output,$INP_length_fwd,$INP_length_rev,$INP_notrim,$INP_concatenate,$INP_threads);

my @argv = @ARGV;

GetOptions(
	'h|help|?' =>  \&usage,
	'i|input=s{1,}' => \@INP_input,
	'd|data=s{,}' => \@INP_amplicons_files,
	'o|output=s' => \$INP_output,
	'fl|flength=i' => \$INP_length_fwd,
	'rl|rlength=i' => \$INP_length_rev,
	'rc|revcomp' => \$INP_revcomp,
	'nt|notrim' => \$INP_notrim,
	'c|concat' => \$INP_concatenate,
	't|threads=i' => \$INP_threads,
);

# Usage help
sub usage
{
  print "usage: $0 -i FILE [options]\n";
  print " -h this message\n";
  print " -d CVS file with primer/amplicon data.\n";
  print " -i input sequence or FASTA file/s or folder\n";
  print " -o output FASTA file\n";
  print " -fl forward tag length\n";
  print " -rl reverse tag length\n";
  print " -rc Analyze direct and reverse complement (default=$INP_revcomp).\n";
  print " -nt Do not trim the sequences.\n";
  print " -c concatenate primers at the beginning and at the end of all the sequences".
        "    (use this option when primers are not present in the reads).\n";
  print " -t Number of threads to calculate the alignments.\n";
  exit(-1);
}

# Prints usage help if no input file is specified
if (!@INP_input || !@INP_amplicons_files){
	print "\nERROR: You must specify input files and primer data.\n\n";
	usage();
	exit;
}
if (!defined($INP_output)){
	print "\nERROR: You must specify output file.\n\n";
	usage();
	exit;
}

print "\nRunning '$COMMAND_LINE'\n";

# Reads sample/tag data from CVS input files if it is specified here
my $sampledata = {};
my $samples = [];
print "\n";
foreach my $INP_amplicons_file (@INP_amplicons_files) {
	print "Reading sample data from file '$INP_amplicons_file'.\n";
	my ($sampledata_, $samples_) = read_amplicon_data($INP_amplicons_file, 'samples');
	if (defined($sampledata_)) {
		$sampledata = { %$sampledata, %$sampledata_ };
		$samples = [ @$samples, @$samples_ ];
	}
}
if (@{$samples}){
	print "\tNumber of samples: ".scalar @{$samples}.".\n";
} elsif (!defined($INP_length_fwd) && !defined($INP_length_rev)){
	print "\nERROR: You must specify at least one tag length or provide tag data in primer/amplicon data file.\n\n";
	usage();
	exit;
}

# Reads markers/primers data from CVS input files
my $markerdata = {};
my $markers = [];
print "\n";
foreach my $INP_amplicons_file (@INP_amplicons_files) {
	print "Reading marker data from file '$INP_amplicons_file'.\n";
	my ($markerdata_, $markers_) = read_amplicon_data($INP_amplicons_file, 'markers');
	$markerdata = { %$markerdata, %$markerdata_ };
	$markers = [ @$markers, @$markers_ ];
}
my ($primer_seqs, $primer_headers);
if (@{$markers}){
	print "\tNumber of markers: ".scalar @$markers.".\n";
	($primer_seqs, $primer_headers) = extract_primer_tag_seqs($markerdata, $markers);
	print "\tNumber of unique primer sequences: ".scalar @{$primer_headers}.".\n";
} else {
	print "\nERROR: You must specify primer sequences in primer/amplicon data file.\n\n";
	usage();
	exit;
}
if ($INP_concatenate && scalar @$markers > 1){
	print "\nERROR: Only one marker is allowed if 'concatenate' option is activated.\n\n";
}

my @options;
if (!$INP_revcomp){
	push(@options,'direct');
}
if (defined($INP_notrim)){
	push(@options,'notrim');
}
if (defined($INP_concatenate)){
	push(@options,'concat');
}

my $added_tags;
my $tags_fwd = [];
my $tags_rev = [];
my ($read_headers_with_tags,$read_seqs_with_tags,$read_quals_with_tags);
foreach my $INP_input (@INP_input) {
	if (-d $INP_input){
		my @files = list_dir($INP_input);
		foreach my $file (@files) {
			if ($file =~ /([^\/]+)(.fa|.fas|.fasta|.fa.gz|.fasta.gz|.fq|.fastq|.fq.gz)$/){
				my $filename = $1;
				print "\nProcessing file '$file'\n";
				# Finds the sample name corresponding to the file
				my $sample_name = '';
				if (my $pos = in_array($samples,$filename,1)){
					$sample_name = $samples->[$pos->[0]];
				}
				# Generates new tags
				my ($tag_fwd, $tag_rev) = ('','');
				if (defined($sampledata->{$sample_name}) && defined($sampledata->{$sample_name}{'tag_f'})){
					$tag_fwd = $sampledata->{$sample_name}{'tag_f'};
				} elsif (defined($INP_length_fwd)){
					$tag_fwd = generate_tag($INP_length_fwd,$tags_fwd);
					push(@$tags_fwd, $tag_fwd);
				}
				if (defined($sampledata->{$sample_name}) && defined($sampledata->{$sample_name}{'tag_r'})){
					$tag_rev = $sampledata->{$sample_name}{'tag_r'};
				} elsif (defined($INP_length_rev)){
					$tag_rev = generate_tag($INP_length_rev,$tags_rev);
					push(@$tags_rev, $tag_rev);
				}
				if (!$tag_fwd && !$tag_rev){
					print "\nERROR: Incorrect tags, sample names must match file names in primer/amplicon data file.\n\n";
					exit;
				}
				# Check and read sequences file
				my ($reads_file_format,$read_seqs,$read_headers,$read_quals,$total_reads)
				= parse_sequence_file($file,undef,['qualities','verbose']);
				my ($added_tags_,$read_seqs_with_tags_,$read_headers_with_tags_,$read_quals_with_tags_) = add_tags($tag_fwd,$tag_rev,$read_headers,$read_seqs,$read_quals,$primer_headers,$primer_seqs,$INP_align,\@options,$INP_threads);
				if (defined($added_tags_)  && %$added_tags_){
					foreach my $marker_name (keys %$added_tags_){
						$added_tags->{$filename}{$marker_name} = $added_tags_->{$marker_name};
					}
					push(@$read_headers_with_tags, @$read_headers_with_tags_);
					push(@$read_seqs_with_tags, @$read_seqs_with_tags_);
					if (defined($read_quals_with_tags_) && @$read_quals_with_tags_) {
						push(@$read_quals_with_tags, @$read_quals_with_tags_);
					}
				} else {
					print "\nERROR: No tags were added to '$file', check if the file contains sequences matching Fwd and Rev primers.\n\n";
				}
			}
		}
	} elsif (-f $INP_input){
		my $file = $INP_input;
		if ($file =~ /([^\/]+)(.fa|.fas|.fasta|.fa.gz|.fasta.gz|.fq|.fastq|.fq.gz)$/){
			my $filename = $1;
			print "\nProcessing file '$file'\n";
			# Finds the sample name corresponding to the file
			my $sample_name = '';
			if (my $pos = in_array($samples,$filename,1)){
				$sample_name = $samples->[$pos->[0]];
			}
			# Generates new tags
			my ($tag_fwd, $tag_rev) = ('','');
			if (defined($sampledata->{$sample_name}) && defined($sampledata->{$sample_name}{'tag_f'})){
				$tag_fwd = $sampledata->{$sample_name}{'tag_f'};
			} elsif (defined($INP_length_fwd)){
				$tag_fwd = generate_tag($INP_length_fwd,$tags_fwd);
				push(@$tags_fwd, $tag_fwd);
			}
			if (defined($sampledata->{$sample_name}) && defined($sampledata->{$sample_name}{'tag_r'})){
				$tag_rev = $sampledata->{$sample_name}{'tag_r'};
			} elsif (defined($INP_length_rev)){
				$tag_rev = generate_tag($INP_length_rev,$tags_rev);
				push(@$tags_rev, $tag_rev);
			}
			if (!$tag_fwd && !$tag_rev){
				print "\nERROR: Incorrect tags, sample names must match file names in primer/amplicon data file.\n\n";
				exit;
			}
			# Check and read sequences file
			my ($reads_file_format,$read_seqs,$read_headers,$read_quals,$total_reads)
			= parse_sequence_file($file,undef,['qualities','verbose']);
			my ($added_tags_,$read_seqs_with_tags_,$read_headers_with_tags_,$read_quals_with_tags_) = add_tags($tag_fwd,$tag_rev,$read_headers,$read_seqs,$read_quals,$primer_headers,$primer_seqs,$INP_align,\@options,$INP_threads);
			if (defined($added_tags_)  && %$added_tags_){
				foreach my $marker_name (keys %$added_tags_){
					$added_tags->{$filename}{$marker_name} = $added_tags_->{$marker_name};
				}
				push(@$read_headers_with_tags, @$read_headers_with_tags_);
				push(@$read_seqs_with_tags, @$read_seqs_with_tags_);
				if (defined($read_quals_with_tags_) && @$read_quals_with_tags_) {
					push(@$read_quals_with_tags, @$read_quals_with_tags_);
				}
			} else {
				print "\nERROR: No tags were added to '$file', check if the file contains sequences matching Fwd and Rev primers.\n\n";
			}
		}
	}
}

print "\nFILE\t\tMARKER\ttag_F\ttag_R\tREADS\n";
my %duplicated_tags;
my $tags_csv;
if (@$tags_fwd && @$tags_rev) {
	$tags_csv = ">sample,tag_f,tag_r\n";
} elsif (@$tags_fwd) {
	$tags_csv = ">sample,tag_f\n";
} elsif (@$tags_rev) {
	$tags_csv = ">sample,tag_r\n";
}
foreach my $filename (sort keys %$added_tags){
	foreach my $marker_name (sort keys %{$added_tags->{$filename}}){
		foreach my $tag (sort { $added_tags->{$filename}{$marker_name}{$b} <=> $added_tags->{$filename}{$marker_name}{$a} } keys %{$added_tags->{$filename}{$marker_name}}){
# 			if (defined($duplicated_tags{$tag})){ next; }
# 			$duplicated_tags{$tag} = 1;
			# if ($added_tags->{$filename}{$marker_name}{$tag}<100) { next; }
			#if ($tag =~ /X/) { next; }
			if ($tag =~ /(.+)-(.+)/) {
				printf("%s\t\t%s\t%s\t%s\t%d\n", $filename, $marker_name, $1, $2, $added_tags->{$filename}{$marker_name}{$tag});
				$tags_csv .= sprintf("%s,%s,%s\n", $filename, $1, $2);
			} elsif ($tag =~ /(.+)-/) {
				printf("%s\t\t%s\t%s\t%s\t%d\n", $filename, $marker_name, $1, '', $added_tags->{$filename}{$marker_name}{$tag});
				$tags_csv .= sprintf("%s,%s\n", $filename, $1);
			} elsif ($tag =~ /-(.+)/) {
				printf("%s\t\t%s\t%s\t%s\t%d\n", $filename, $marker_name, '', $1, $added_tags->{$filename}{$marker_name}{$tag});
				$tags_csv .= sprintf("%s,%s\n", $filename, $1);
			}
		}
	}
}
print "\n";

write_to_file("$INP_output.csv", $tags_csv);
if (defined($read_quals_with_tags) && @$read_quals_with_tags){
	create_fastq_file($read_seqs_with_tags,$read_headers_with_tags,$read_quals_with_tags,"$INP_output.fq");
	`gzip -qf $INP_output.fq`;
	printf("\ntag info stored in '%s' and tagd sequences in '%s'\n\n","$INP_output.csv","$INP_output.fq.gz");
} else {
	create_fasta_file($read_seqs_with_tags,$read_headers_with_tags,"$INP_output.fa");
	`gzip -qf $INP_output.fa`;
	printf("\ntag info stored in '%s' and tagd sequences in '%s'\n\n","$INP_output.csv","$INP_output.fa.gz");
}
exit;


sub add_tags {

	my ($tag_fwd,$tag_rev,$read_headers,$read_seqs,$read_quals,$primer_headers,$primer_seqs,$aligntype,$options,$INP_threads) = @_;

	my ($added_tags,$read_headers_with_tags,$read_seqs_with_tags,$read_quals_with_tags);

	my ($revcomp,$trim,$concat)=(1,1,0);
	if (in_array($options, 'direct')){
		$revcomp = 0;
	}
	if (in_array($options, 'notrim')){
		$trim = 0;
	}
	if (in_array($options, 'concat')){
		$concat = 1;
	}

	# If only concatenation of primers and tags is desired
	if ($concat) {
	
		# Annotates the name of the marker
		$primer_headers->[0] =~ /(.+)_([F|R])\d+/;
		my $marker_name = $1;

		# Divides the primers into Fwd and Rev
		my (@primers_fwd, @primers_rev);
		for (my $i=0; $i<=$#{$primer_headers}; $i++){
			if ($primer_headers->[$i] =~ /F\d+$/){
				push(@primers_fwd, $primer_seqs->[$i]);
			} elsif  ($primer_headers->[$i] =~ /R\d+$/){
				push(@primers_rev, $primer_seqs->[$i]);
			}
		}

		# Loops reads
		for (my $i=0; $i<=$#{$read_headers}; $i++){
			my $read_seq_with_tags = $read_seqs->[$i];
			my ($primer_fwd, $primer_rev) = ('','');
			if (@primers_fwd){
				$primer_fwd = $primers_fwd[rand @primers_fwd];
				$read_seq_with_tags = $primer_fwd.$read_seq_with_tags;
			}
			# @primers_rev contains reverse complementary sequences of reverse primers
			if (@primers_rev){
				$primer_rev = $primers_rev[rand @primers_rev];
				$read_seq_with_tags = $read_seq_with_tags.$primer_rev;
			}
			$read_seq_with_tags = $tag_fwd.$read_seq_with_tags;
			if ($tag_rev) {
				$read_seq_with_tags = $read_seq_with_tags.iupac_reverse_complementary($tag_rev);
			}
			push(@$read_headers_with_tags, $read_headers->[$i]);
			push(@$read_seqs_with_tags, $read_seq_with_tags);
			if (defined($read_quals)){
				my $read_qual_with_tags = $read_quals->[$i];
				$read_qual_with_tags = '?'x(length($tag_fwd.$primer_fwd)).$read_qual_with_tags.'?'x(length($primer_rev.$tag_rev));
				push(@$read_quals_with_tags, $read_qual_with_tags);
			}
		}

		if ($tag_fwd && $tag_rev){
			$added_tags->{$marker_name}{"$tag_fwd-$tag_rev"} += scalar @$read_headers_with_tags;
		} elsif ($tag_fwd && !$tag_rev){
			$added_tags->{$marker_name}{"$tag_fwd-X"} += scalar @$read_headers_with_tags;
		} elsif (!$tag_fwd && $tag_rev){
			$added_tags->{$marker_name}{"X-$tag_rev"} += scalar @$read_headers_with_tags;
		}

	# If search for primers and tag insertion is desired
	} else {
		# Aligns reads against primer sequences
	# 	print "\nAligning primer sequences.\n";
		my $align_amplicon_data
		= align_amplicons($read_headers,$read_seqs,$primer_headers,$primer_seqs,$aligntype,$revcomp,$INP_threads);
		
		# Loops reads with alignment results
		for (my $i=0; $i<=$#{$read_headers}; $i++){
			my $read_header = $read_headers->[$i];

			if (!defined($align_amplicon_data->{$read_header})){
				next;
			}
			my $read_seq = $read_seqs->[$i];
			my $read_length = length($read_seq);

			# Find between results matches with common forward and reverse primer sequences
			my ($forward_seqs, $reverse_seqs);
			my $amplicon_found;
			foreach my $result (@{$align_amplicon_data->{$read_header}}) {
				my $primer_header = $result->{'NAME'};
				$primer_header =~ /(.+)_([F|R])\d+/;
				# If matched sequence is forward primer
				if ($2 eq 'F'){
					# Annotate forward sequences matched
					if (!defined($forward_seqs->{$1})){
						$forward_seqs->{$1} = $result;
					}
					# Stop checking results if the same amplicon has been detected in reverse sequences
					if (defined($reverse_seqs) && defined($reverse_seqs->{$1})){
						$amplicon_found = $1;
						last;
					}
				# If matched sequence is reverse primer
				} elsif ($2 eq 'R'){
					# Annotate reverse sequences matched
					if (!defined($reverse_seqs->{$1})){
						$reverse_seqs->{$1} = $result;
					}
					# Stop checking results if the same primer has been detected in forward sequences
					if (defined($forward_seqs) && defined($forward_seqs->{$1})){
						$amplicon_found = $1;
						last;
					}
				}
			}
			if (defined($amplicon_found)) {
				push(@$read_headers_with_tags, $read_header);
				# Extracts names and sequences of the primers, it should be the same name for forward and reverse
				$forward_seqs->{$amplicon_found}{'NAME'} =~ /(.+)_([F|R])\d+/;
				my $marker_name = $1;
				# Annotates the tag sequence positions
				my ($forward_read_aligned_cols, $forward_amplicon_aligned_cols) = split("\n",$forward_seqs->{$amplicon_found}{'COLS'});
				my @forward_read_aligned_cols = split(",",$forward_read_aligned_cols);
				my ($reverse_read_aligned_cols, $reverse_amplicon_aligned_cols) = split("\n",$reverse_seqs->{$amplicon_found}{'COLS'});
				my @reverse_read_aligned_cols = split(",",$reverse_read_aligned_cols);
				# Check if alignment is direct or reverse complementary
				my ($first_fwd_tag_pos,$first_rev_tag_pos);
				if ($forward_read_aligned_cols[0]<$reverse_read_aligned_cols[0]) {
					$first_fwd_tag_pos = $forward_read_aligned_cols[0];
					$first_rev_tag_pos = $reverse_read_aligned_cols[-1]+1;
					# Adds tag sequences
					if ($trim){
						my $read_seq_with_tags = $tag_fwd.substr($read_seq,$first_fwd_tag_pos-1, $first_rev_tag_pos-$first_fwd_tag_pos);
						if ($tag_rev) {
							$read_seq_with_tags .= iupac_reverse_complementary($tag_rev);
						}
						push(@$read_seqs_with_tags, $read_seq_with_tags);
					} else {
						my $read_seq_with_tags = substr($read_seq,0,$first_fwd_tag_pos-1).$tag_fwd;
						$read_seq_with_tags .= substr($read_seq,$first_fwd_tag_pos-1, $first_rev_tag_pos-$first_fwd_tag_pos);
						if ($tag_rev) {
							$read_seq_with_tags .= iupac_reverse_complementary($tag_rev).substr($read_seq,$first_rev_tag_pos-1);
						}
						push(@$read_seqs_with_tags, $read_seq_with_tags);
					}
					if (defined($read_quals)){
						my $read_qual = $read_quals->[$i];
						if ($trim){
							my $read_qual_with_tags = '?'x(length($tag_fwd)).substr($read_qual,$first_fwd_tag_pos-1, $first_rev_tag_pos-$first_fwd_tag_pos);
							$read_qual_with_tags .= '?'x(length($tag_rev));
							push(@$read_quals_with_tags, $read_qual_with_tags);
						} else {
							my $read_qual_with_tags = substr($read_qual,0,$first_fwd_tag_pos-1).'?'x(length($tag_fwd));
							$read_qual_with_tags .= substr($read_qual,$first_fwd_tag_pos-1, $first_rev_tag_pos-$first_fwd_tag_pos);
							$read_qual_with_tags .= '?'x(length($tag_rev)).substr($read_qual,$first_rev_tag_pos-1);
							push(@$read_quals_with_tags, $read_qual_with_tags);
						}
					}
				} else {
					$first_rev_tag_pos = $reverse_read_aligned_cols[0];
					$first_fwd_tag_pos = $forward_read_aligned_cols[-1]+1;
					# Adds tag sequences
					if ($trim){
						my $read_seq_with_tags = $tag_rev.substr($read_seq,$first_rev_tag_pos-1, $first_fwd_tag_pos-$first_rev_tag_pos);
						$read_seq_with_tags .= iupac_reverse_complementary($tag_fwd);
						push(@$read_seqs_with_tags, $read_seq_with_tags);
					} else {
						my $read_seq_with_tags = substr($read_seq,0,$first_rev_tag_pos-1).$tag_rev;
						$read_seq_with_tags .= substr($read_seq,$first_rev_tag_pos-1, $first_fwd_tag_pos-$first_rev_tag_pos);
						$read_seq_with_tags .= iupac_reverse_complementary($tag_fwd).substr($read_seq,$first_fwd_tag_pos-1);
						push(@$read_seqs_with_tags, $read_seq_with_tags);
					}
					if (defined($read_quals)){
						my $read_qual = $read_quals->[$i];
						if ($trim){
							my $read_qual_with_tags = '?'x(length($tag_rev)).substr($read_qual,$first_rev_tag_pos-1, $first_fwd_tag_pos-$first_rev_tag_pos);
							$read_qual_with_tags .= '?'x(length($tag_fwd));
							push(@$read_quals_with_tags, $read_qual_with_tags);
						} else {
							my $read_qual_with_tags = substr($read_qual,0,$first_rev_tag_pos-1).'?'x(length($tag_rev));
							$read_qual_with_tags .= substr($read_qual,$first_rev_tag_pos-1, $first_fwd_tag_pos-$first_rev_tag_pos);
							$read_qual_with_tags .= '?'x(length($tag_fwd)).substr($read_qual,$first_fwd_tag_pos-1);
							push(@$read_quals_with_tags, $read_qual_with_tags);
						}
					}
				}
				if ($tag_fwd && $tag_rev){
					$added_tags->{$marker_name}{"$tag_fwd-$tag_rev"}++;
				} elsif ($tag_fwd && !$tag_rev){
					$added_tags->{$marker_name}{"$tag_fwd-X"}++;
				} elsif (!$tag_fwd && $tag_rev){
					$added_tags->{$marker_name}{"X-$tag_rev"}++;
				}
			} elsif (!$trim){
				push(@$read_headers_with_tags, $read_header);
				push(@$read_seqs_with_tags, $read_seq);
				if (defined($read_quals)){
					my $read_qual = $read_quals->[$i];
					push(@$read_quals_with_tags, $read_qual);
				}
			}
		}
	}

	return ($added_tags,$read_seqs_with_tags,$read_headers_with_tags,$read_quals_with_tags);

}

sub list_dir {
        my $dir = shift @_;
        my @files;
        find({ wanted => sub { push @files, $_ } , no_chdir => 1 }, ($dir));
        return @files;
}





