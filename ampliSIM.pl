#!/usr/bin/perl -w
#
################################################################
#
# Name: ampliSIM.pl
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
#   Generates amplicon sequencing simulated sequences/reads containing the primers and tags specified in a CSV with amplicon data
#   Also desired alleles can be specified in an additional FASTA file
#   Generates a FASTA format file with the generated sequences
#
# Requires as input a CVS format file with primer/amplicon and tag/sample data.
#
# Example:
# perl ampliSIM.pl -d primers_tags.cvs -a alleles.fa -o simulated.reads
#

my $VERSION = "1.2";
my $SCRIPT_NAME = fileparse($0);
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Generates amplicon sequencing simulated sequences/reads containing the primers and tags specified in a CSV with amplicon data.";

# Modules are in folder 'lib' in the path of the script
use File::FindLib 'lib';
# Perl modules necessaries for the correct working of the script
use Cwd;
use File::Basename;
use Getopt::Long;
use Bio::Sequences;
use Bio::Ampli;
use List::Util qw(shuffle);
# use Data::Dumper;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;


my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options
# By default creates only direct reads
my $INP_revcomp = 0;
# Clustering parameters
my @ALL_SEQUENCING_PARAMS = ('substitution_threshold', 'indel_threshold', 'chimera_frequency', 'error_models');
# Sequence filtering parameters
my @ALL_AMPLICON_PARAMS = ('min_allele_number', 'max_allele_number', 'min_amplicon_depth', 'max_amplicon_depth');
# Default amplicon parameters
my $DEFAULT_COMMON_PARAMS = {
	'min_allele_number' => { 'all' => [ 1 ] },
	'max_allele_number' => { 'all' => [ 2 ] },
	'min_amplicon_depth' => { 'all' => [ 1000 ] },
	'max_amplicon_depth' => { 'all' => [ 1000 ] },
# 	'chimera_frequency' => { 'all' => [ 15 ] },
};
# Default recommended technology parameters
my $TECH_PARAMS = {
	'pyroseq' => { # 454 and IonTorrent
		'substitution_threshold' => { 'all' => [ '0.01' ] },
		'indel_threshold' => { 'all' => [ '1' ] },
		%$DEFAULT_COMMON_PARAMS,
	},
	'illumina' => { # Illumina
		'substitution_threshold' => { 'all' => [ '0.1' ] },
		'indel_threshold' => { 'all' => [ '0.001' ] },
		%$DEFAULT_COMMON_PARAMS,
	},
};
$TECH_PARAMS->{'454/iontorrent'} = $TECH_PARAMS->{'454'} = $TECH_PARAMS->{'iontorrent'} = $TECH_PARAMS->{'pyroseq'};
$TECH_PARAMS->{'solexa'} = $TECH_PARAMS->{'illumina'};

# Chimeras default frequency as in Biedriczka eat al 2016
my $CHIMERA_DEFAULT_FREQ = 15; # 15% of the reads (0.003% into an amplicon)

# Error prone sequence motifs by technology
my $ERROR_MODELS = {
	'pyroseq' => { # 454 and IonTorrent
		'AAACA' => { 'error' => 'AAAAA', 'freq' => 1.07 }, # Frequencies are %
		'CCCAC' => { 'error' => 'CCCCC', 'freq' => 1.02 },
		'CCCCG' => { 'error' => 'CCCAG', 'freq' => 0.75 },
		'AAAGG' => { 'error' => 'AAAAG', 'freq' => 0.70 },
		'AGGAA' => { 'error' => 'AGGGA', 'freq' => 0.52 },
	},
	'illumina' => { # Illumina
		'GGGTC' => { 'error' => 'GGGGC', 'freq' => 3.94 }, # Frequencies are %
		'CTCGG' => { 'error' => 'CTCCG', 'freq' => 2.92 },
		'GGGCG' => { 'error' => 'GGGGG', 'freq' => 2.03 },
		'CGGTG' => { 'error' => 'CGGGG', 'freq' => 3.18 },
		'GGGTA' => { 'error' => 'GGGGA', 'freq' => 1.60 },
		'AGGTG' => { 'error' => 'AGGGG', 'freq' => 1.84 },
		'GGGTG' => { 'error' => 'GGGGG', 'freq' => 1.22 },
		'CGGTC' => { 'error' => 'CGGGC', 'freq' => 0.99 },
	},
};
$ERROR_MODELS->{'454/iontorrent'} = $ERROR_MODELS->{'454'} = $ERROR_MODELS->{'iontorrent'} = $ERROR_MODELS->{'pyroseq'};
$ERROR_MODELS->{'solexa'} = $ERROR_MODELS->{'illumina'};



my ($INP_amplicons_file, $INP_outfile, $INP_allele_file, $INP_input_file, $INP_tech, $INP_chimeras, $INP_error_models, $INP_nreads, $INP_nreads_amplicon, $INP_shuffle, $INP_zip, $INP_gzip);

GetOptions(
	'h|help|?' =>  \&usage,
	'd|data=s{,}' => \$INP_amplicons_file,
	'a|alleles=s' => \$INP_allele_file,
	'i|input=s' => \$INP_input_file,
	't|tech=s' => \$INP_tech,
	'chi:i' => \$INP_chimeras,
	'mod|models' => \$INP_error_models,
	'o|output=s' => \$INP_outfile,
	'n|number=i' => \$INP_nreads,
	'na|nampli=i' => \$INP_nreads_amplicon,
	'r|revcomp=i' => \$INP_revcomp,
	's|shuffle' => \$INP_shuffle,
	'z|zip' => \$INP_zip,
	'gz|gzip' => \$INP_gzip,
	'<>' => \&usage,
);

# Usage help
sub usage {
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -d <file>\tCSV file with primer/amplicon data.\n";
	print "  -o <file>\tOutput file name.\n";
	print "  -a <file>\tFASTA file with allele names and sequences.\n";
	print "  -i <file>\tExcel file with alleles and frequencies to be used as simulation templates.\n";
	print "  -t <tech>\tSequencing technology ('454', 'IonTorrent', 'Illumina', 'Unknown').\n";
	print "  -chi <number>\tSimulates chimeras in the specified percentage of reads (default=$CHIMERA_DEFAULT_FREQ%).\n";
	print "  -mod\t\tIntroduces non-random substitutions using error model motifs (McElroy et al. 2012).\n";
	print "  -n <number>\tTotal number of reads/sequences to generate.\n";
	print "  -na <number>\tNumber of reads/sequences per amplicon to generate.\n";
	print "  -r\t\tCreate reads in both: direct and reverse complement senses (default=$INP_revcomp)\n";
	print "  -s\t\tShuffle/randomize generated reads.\n";
	print "  -z\t\tCompress output file in ZIP format.\n";
	print "  -gz\t\tCompress output file in GZIP format.\n";
	print "  -h\t\tHelp.\n";
	print "\n";
	exit;
}


# Prints usage help if no input file is specified
if (!defined($INP_amplicons_file)){
	print "\nERROR: You must specify amplicon data file.\n\n";
	usage();
	exit;
}
# Checks technology default parameters
if (defined($INP_tech)){
	$INP_tech =~ s/\s//g;
	$INP_tech = lc($INP_tech);
	if (!defined($TECH_PARAMS->{$INP_tech})){ 
		print "\nERROR: You must specify a valid sequencing technology ('Illumina', '454' or 'IonTorrent').\n\n";
		usage();
		exit;
	}
}

# Allele and Excel input files are not compatible
if (defined($INP_allele_file) && defined($INP_input_file)){
	print "\nERROR: Alleles and Excel input files are not compatible options, specify only one as input.\n\n";
}
# Gives default name to output file
if (!defined($INP_outfile)){
	if ($INP_amplicons_file =~ /^(.+?)\./) {
		$INP_outfile = $1.".reads";
	} else {
		$INP_outfile = "$INP_amplicons_file.reads";
	}
}

print "\nRunning '$COMMAND_LINE'\n";


# 2 PRIMERS => MARKER
# 1/2 TAGS => SAMPLE
# 2 PRIMERS + 1/2 TAGS => AMPLICON (Single PCR product)

# Check and read amplicons file
my ($markerdata,$primers,$sampledata,$tags,$paramsdata,$alleledata)
= parse_amplicon_file($INP_amplicons_file,['skip errors']);

# If a file with variants and frequencies is provided
my $template_allele_data;
if (defined($INP_input_file) && -e $INP_input_file){
	printf("\nReading File '%s'.\n", $INP_input_file);
	$template_allele_data = read_amplisas_file_results($INP_input_file);
}

if (!defined($INP_nreads) && !defined($INP_nreads_amplicon) && !defined($paramsdata->{'min_amplicon_depth'}) && !defined($template_allele_data)){
	print "\nERROR: You must specify the number of reads/sequences total or per amplicon to generate.\n\n";
	usage();
	exit;
}

# Check and read alleles file (optional)
# Amplicon data file has preference over alleles in FASTA file
my ($allele_seqs,$allele_headers);
if (defined($INP_allele_file) && -e $INP_allele_file){
	print "\nReading allele sequences from '$INP_allele_file'.\n";
	($allele_seqs,$allele_headers) = read_fasta_file($INP_allele_file);
}

# Retrieves sequencing and amplicon parameters from the CSV file
# If not defined and technology is specified, technogy defaults are used
my @SEQUENCING_PARAMS;
my @AMPLICON_PARAMS;
# Sets minimum per amplicon frequency if it's specified in command line
# Command line params have priority over CSV file and auto mode
if (defined($INP_chimeras) && $INP_chimeras==0){
	$paramsdata->{'chimera_frequency'}{'all'} = [ $CHIMERA_DEFAULT_FREQ ];
} elsif (defined($INP_chimeras)){
	$paramsdata->{'chimera_frequency'}{'all'} = [ $INP_chimeras ];
}
if (defined($INP_error_models)){
	$paramsdata->{'error_models'}{'all'} = [ ];
}
# Read always the params from CSV file, they will have priority over auto ones
# Sequencing parameters are mandatory
foreach my $paramname (@ALL_SEQUENCING_PARAMS){
	if (defined($paramsdata->{$paramname})) {
		push(@SEQUENCING_PARAMS,$paramname);
	} elsif (defined($INP_tech) && defined($TECH_PARAMS->{$INP_tech}{$paramname})){
		$paramsdata->{$paramname} = $TECH_PARAMS->{$INP_tech}{$paramname};
		push(@SEQUENCING_PARAMS,$paramname);
	}
}
if (defined($paramsdata->{'error_models'}) && (!defined($INP_tech) || !defined($ERROR_MODELS->{$INP_tech}))){
	print "\nERROR: You must specify a valid sequencing technology ('Illumina', '454' or 'IonTorrent').\n\n";
	usage();
	exit;
}
if (!@SEQUENCING_PARAMS) {
	print "\nERROR: No sequencing error parameters are specified. You can use default parameters by specifying the sequencing technology or include them into CSV data file.\n\n";
	usage();
	exit;
}
if (!defined($template_allele_data)){
	foreach my $paramname (@ALL_AMPLICON_PARAMS){
		if (defined($paramsdata->{$paramname})) {
			push(@AMPLICON_PARAMS,$paramname);
		} elsif (defined($INP_tech) && defined($TECH_PARAMS->{$INP_tech}{$paramname})){
			$paramsdata->{$paramname} = $TECH_PARAMS->{$INP_tech}{$paramname};
			push(@AMPLICON_PARAMS,$paramname);
		}
	}
	if (!@AMPLICON_PARAMS) {
		print "\nERROR: No amplicon parameters are specified. You can use default parameters by specifying the sequencing technology used or include them into CSV data file.\n\n";
		usage();
		exit;
	}
}



# Generates the reads, amplicon by amplicon
print "\nGenerating amplicon sequences with the following parameters ('threshold' 'marker' 'min value' 'max value'):\n";
# In $paramsdata hash there are also sequence filters
foreach my $paramname (@AMPLICON_PARAMS){
	if (defined($paramsdata->{$paramname})) {
		foreach my $marker_name (sort keys %{$paramsdata->{$paramname}}){
			printf("\t%s\t%s\t%s\n",$paramname,$marker_name,join("\t",@{$paramsdata->{$paramname}{$marker_name}}));
		}
	}
}

print "\n";
my $verbose = 1;
my ($read_seqs_generated,$read_headers_generated);
my %allele_lengths;
foreach my $marker_name (@$primers){

	if (!defined($markerdata->{$marker_name}{'primer_f'}) && !defined($markerdata->{$marker_name}{'primer_r'}) && !defined($markerdata->{$marker_name}{'primer_rc'})){
		printf("\nERROR: Marker '%s' has no primers defined.\n\n", $marker_name);
		exit;
	}

	my @primers_f;
	if (defined($markerdata->{$marker_name}{'primer_f'})){
		foreach my $primer_f (@{$markerdata->{$marker_name}{'primer_f'}}) {
			push(@primers_f, unambiguous_dna_sequences($primer_f));
		}
	} else {
		@primers_f = ('');
	}
	
	my @primers_rc;
	if (defined($markerdata->{$marker_name}{'primer_rc'})){
		foreach my $primer_rc (@{$markerdata->{$marker_name}{'primer_rc'}}) {
			push(@primers_rc, unambiguous_dna_sequences($primer_rc));
		}
	} elsif (defined($markerdata->{$marker_name}{'primer_r'})){
		foreach my $primer_r (@{$markerdata->{$marker_name}{'primer_r'}}){
			push(@primers_rc, unambiguous_dna_sequences(uc(iupac_reverse_complementary($primer_r))));
		}
	} else {
		@primers_rc = ('');
	}

	# If no tags are provided, it will include only primers
	if (!defined($tags) || !@$tags) {
		$sampledata->{''}{'tag_f'} = '';
		$sampledata->{''}{'tag_rc'} = '';
		@$tags = ('');
	}

	# Defines max and min amplicon depths
	my ($min_depth,$max_depth);
	# Command line params have priority over CSV file and auto mode
	if (defined($INP_nreads_amplicon)){
		$min_depth = $max_depth = $INP_nreads_amplicon;
	} elsif (defined($INP_nreads)){
		$min_depth = $max_depth = $INP_nreads / (scalar @$primers * scalar @$tags);
	} elsif (!defined($template_allele_data)){
		$min_depth = $paramsdata->{'min_amplicon_depth'}{'all'}[0];
		if (defined($paramsdata->{'max_amplicon_depth'})){
			$max_depth = $paramsdata->{'max_amplicon_depth'}{'all'}[0];
		} else {
			$max_depth = $paramsdata->{'min_amplicon_depth'}{'all'}[0];
		}
	}
	
	# Generates a finite number of alleles to use among the amplicons
	my $marker_alleles;
	if (!defined($template_allele_data)){
		if (defined($allele_seqs)){
			# if more than 1 marker, alleles will be selected if they contain in the header the marker name
			for (my $i=0; $i<$#{$allele_headers}; $i++){
				if (scalar @$primers == 1 || $allele_headers->[$i] =~ /$marker_name/){
					push(@$marker_alleles,[$allele_seqs->[$i],$allele_headers->[$i]]);
				}
			}
		} else {
			my $total_alleles = scalar @$tags * mean($paramsdata->{'min_allele_number'}{'all'}[0],$paramsdata->{'max_allele_number'}{'all'}[0]);
			my $allele_id_length = length($total_alleles);
			for (my $i=0; $i<$total_alleles; $i++){
				if (!defined($markerdata->{$marker_name}{'length'})){
					print "\nERROR: Allele length/s must be specified together with primers and tags into the amplicon data.\n\n";
					exit;
				}
				my $length = $markerdata->{$marker_name}{'length'}->[rand @{$markerdata->{$marker_name}{'length'}}];
				push(@$marker_alleles,[random_sequence($length), sprintf("%0".$allele_id_length."u",$i+1)]);
			}
		}
	}

	foreach my $sample_name (@$tags) {

		if ($verbose && $sample_name ne ''){
			print "\t$marker_name-$sample_name processing\n";
		} elsif ($verbose) {
			print "\t$marker_name processing\n";
		}

		# Retrieves amplicon primer and tag sequences
		my $primer_f = $primers_f[rand @primers_f];
		my $primer_rc = $primers_rc[rand @primers_rc];
		my ($tag_f, $tag_rc) = ('', '');
		if (defined($sampledata->{$sample_name}{'tag_f'})){
			$tag_f = $sampledata->{$sample_name}{'tag_f'};
		}
		if (defined($sampledata->{$sample_name}{'tag_rc'})){
			$tag_rc = $sampledata->{$sample_name}{'tag_rc'};
		}

		# Assigns a depth to the amplicon
		my $depth;
		if (!defined($INP_nreads_amplicon) && defined($template_allele_data) && defined($template_allele_data->{$marker_name}) && defined($template_allele_data->{$marker_name}{'sample_data'}{$sample_name}{'depth_amplicon'})){
			$depth = $template_allele_data->{$marker_name}{'sample_data'}{$sample_name}{'depth_amplicon'};
		} else {
			$depth = $min_depth + int(rand($max_depth-$min_depth+1));
		}

		# Retrieves alleles specified in the Excel input file
		my (@amplicon_allele_seqs,@amplicon_allele_names,@amplicon_allele_efficiencies);
		if (defined($template_allele_data) && defined($template_allele_data->{$marker_name}{'seq_md5s'})){
			foreach my $md5 (@{$template_allele_data->{$marker_name}{'seq_md5s'}}){
				if (defined($template_allele_data->{$marker_name}{'assignments'}{$md5}{$sample_name})){
					push(@amplicon_allele_seqs, $template_allele_data->{$marker_name}{'seq_data'}{$md5}{'sequence'});
					push(@amplicon_allele_names, $template_allele_data->{$marker_name}{'seq_data'}{$md5}{'name'});
					if (defined($template_allele_data->{$marker_name}{'seq_data'}{$md5}{'efficiency'})){
						push(@amplicon_allele_efficiencies,$template_allele_data->{$marker_name}{'seq_data'}{$md5}{'efficiency'});
					} elsif (defined($template_allele_data->{$marker_name}{'seq_data'}{$md5}{'mean_freq'})){
						push(@amplicon_allele_efficiencies,$template_allele_data->{$marker_name}{'seq_data'}{$md5}{'mean_freq'});
					}
				}
			}
		# Generates de novo a finite number of alleles for the amplicon
		} else {
			my $min_alleles = $paramsdata->{'min_allele_number'}{'all'}[0];
			my $max_alleles = $paramsdata->{'max_allele_number'}{'all'}[0];
			my $n_alleles = $min_alleles + int(rand($max_alleles-$min_alleles+1));
			
			foreach my $i (splice([shuffle(0..$#{$marker_alleles})],0,$n_alleles)){
				push(@amplicon_allele_seqs, $marker_alleles->[$i][0]);
				push(@amplicon_allele_names, $marker_alleles->[$i][1]);
			}
		}

		# Sets the length of the sequence identifiers
		my $seq_id_length = length($depth);
		my $seq_id = 0;

		# If there is only 1 allele or not efficiencies are specified
		# Generates random number of reads for each allele
		if ($#amplicon_allele_seqs == 0 || !@amplicon_allele_efficiencies || $#amplicon_allele_efficiencies != $#amplicon_allele_seqs){
			for (my $i=0; $i<$depth; $i++){
				$seq_id++;
				my $rand = int(rand @amplicon_allele_seqs);
				if ($INP_revcomp && int(rand(2)) == 1) { # int(rand(2)) generates 0 or 1
					push(@$read_seqs_generated, iupac_reverse_complementary($tag_f.$primer_f.$amplicon_allele_seqs[$rand].$primer_rc.$tag_rc));
				} else {
					push(@$read_seqs_generated, $tag_f.$primer_f.$amplicon_allele_seqs[$rand].$primer_rc.$tag_rc);
				}
				push(@$read_headers_generated, sprintf("SEQ%0".$seq_id_length."u | allele=%s | amplicon=%s-%s", $seq_id, $amplicon_allele_names[$rand], $marker_name, $sample_name));
				$allele_lengths{length($read_seqs_generated->[-1])} = 1;
			}
		# Generates reads based on allele efficiencies
		} else {
			# Calculates acumulated efficiencies
# 			my $sum = sum(@amplicon_allele_efficiencies);
			my (@cum_efficiencies,$cum_efficiency);
			foreach my $efficiency (@amplicon_allele_efficiencies){
				$cum_efficiency += $efficiency;
				push(@cum_efficiencies, $cum_efficiency);
			}
			for (my $i=0; $i<$depth; $i++){
				$seq_id++;
				my $rand_effi = rand($cum_efficiency);
				my $prev_effi = 0;
				for (my $j=0; $j<=$#cum_efficiencies; $j++){
					if ($rand_effi>=$prev_effi && $rand_effi<$cum_efficiencies[$j]) {
						if ($INP_revcomp && int(rand(2)) == 1) { # int(rand(2)) generates 0 or 1
							push(@$read_seqs_generated, iupac_reverse_complementary($tag_f.$primer_f.$amplicon_allele_seqs[$j].$primer_rc.$tag_rc));
						} else {
							push(@$read_seqs_generated, $tag_f.$primer_f.$amplicon_allele_seqs[$j].$primer_rc.$tag_rc);
						}
						push(@$read_headers_generated, sprintf("SEQ%0".$seq_id_length."u | allele=%s | amplicon=%s-%s", $seq_id, $amplicon_allele_names[$j], $marker_name, $sample_name));
						$allele_lengths{length($read_seqs_generated->[-1])} = 1;
						last;
					}
					$prev_effi = $cum_efficiencies[$j];
				}
			}
		}

		if ($verbose && $sample_name ne ''){
			printf("\t%s-%s generated with %d sequences (%d alleles)\n", $marker_name, $sample_name, $depth, scalar @amplicon_allele_seqs);
		} elsif ($verbose) {
			printf("\t%s generated with %d sequences (%d alleles)\n", $marker_name, $depth, scalar @amplicon_allele_seqs);
		}
	}
}

# if (defined($INP_nreads)){
# 	$read_seqs_generated = splice_seqs($INP_nreads,shuffle_seqs($read_seqs_generated,$read_headers_generated));
# }
if (defined($INP_shuffle)){
	($read_seqs_generated,$read_headers_generated) = shuffle_seqs($read_seqs_generated,$read_headers_generated);
}


# Insert errors into sequences
print "\nInserting sequencing errors with the following parameters ('threshold' 'marker' 'value'):\n";
# In $paramsdata hash there are also sequence filters
foreach my $paramname (@SEQUENCING_PARAMS){
	if (defined($paramsdata->{$paramname})) {
		foreach my $marker_name (sort keys %{$paramsdata->{$paramname}}){
			printf("\t%s\t%s\t%s\n",$paramname,$marker_name,join(',',@{$paramsdata->{$paramname}{$marker_name}}));
		}
	}
}
print "\n";

# Calculates accumulated probabilities for chimeras


# Calculates accumulated binomial probabilities for errors for each allele length
my $error_cum_probabilities;
my @allele_lengths = sort {$a<=>$b} keys %allele_lengths;
# Chimeras from sequences with insertions or deletions can have lower or longer lengths (+-5)
foreach my $len ($allele_lengths[0]-5..$allele_lengths[-1]+5){
	if (defined($error_cum_probabilities->{$len})){
		next;
	}
	foreach my $error_threshold ( ('substitution_threshold', 'indel_threshold') ){
		if (!defined($paramsdata->{$error_threshold})) {
			next;
		}
		my $error_prob = $paramsdata->{$error_threshold}{'all'}[0] / 100;
		my $errors = 0;
		my $cum_bin_prob = 0;
		while (1-$cum_bin_prob > 0.00001){
			$cum_bin_prob += binomial_probability($len,$errors,$error_prob);
			push(@{$error_cum_probabilities->{$len}{$error_threshold}}, $cum_bin_prob);
			$errors++;
		}
	}
}
# exit;

my $substitution_threshold = 100/$paramsdata->{'substitution_threshold'}{'all'}[0];
my $indel_threshold = 100/$paramsdata->{'indel_threshold'}{'all'}[0];
my $count_errors = {	'noerror' => 0,
			'chimeras' =>0,
			'modelled-substitutions' =>0,
			'substitutions' => {},
			'indels' => {},
};
my $chimera_frequency;
if (defined($paramsdata->{'chimera_frequency'})){
	$chimera_frequency = $paramsdata->{'chimera_frequency'}{'all'}[0];
}
my $error_models;
if (defined($paramsdata->{'error_models'})){
	$error_models = 1;
}
for (my $i=0; $i<=$#{$read_seqs_generated}; $i++) {
	
	my $len = length($read_seqs_generated->[$i]);

	# If chimeras are desired
	my $is_chimera = 0;
	if (defined($chimera_frequency)){
		my $rand_freq = rand(100);
		if ($rand_freq < $chimera_frequency) {
			my $rand_seq = $read_seqs_generated->[int(rand($i-1))];
			my $len_rand_seq = length($rand_seq);
			my $rand_pos = 0;
			while ($rand_pos<10){
				if ($len <= $len_rand_seq){
					$rand_pos = int(rand($len));
				} else {
					$rand_pos = int(rand($len_rand_seq));
				}
			}
			if (int(rand(2)) == 0){ # int(rand(2)) generates 0 or 1
				$read_seqs_generated->[$i] = substr($read_seqs_generated->[$i],0,$rand_pos).substr($rand_seq,$rand_pos);
			} else {
				$read_seqs_generated->[$i] = substr($rand_seq,0,$rand_pos).substr($read_seqs_generated->[$i],$rand_pos);
			}
			$count_errors->{'chimeras'}++;
			$is_chimera = 1;
# 			printf("\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",$len_seq1,$len_seq2,$rand_pos,$rand_seq1,$rand_seq2,substr($rand_seq1,0,$rand_pos),substr($rand_seq2,$rand_pos));
# 			exit;
		}
	}
	# Reads again the length of the read (can be different if a chimera has been introduced)
	$len = length($read_seqs_generated->[$i]);

# 	# Chimeras from sequences with insertions or deletions can have lower or longer lengths
# 	# Skip these sequences from introducing more errors
# 	if (!defined($error_cum_probabilities->{$len})){
# 		$count_errors->{'wrong-length-chimeras'}++;
# 	}

	my ($substitutions, $insertions, $deletions, $model_substitutions) = ([],[],[],[]);

	my ($subst_quot, $indel_quot);
	foreach my $error_threshold ( ('substitution_threshold', 'indel_threshold') ){
		my @cum_probabilities = @{$error_cum_probabilities->{$len}{$error_threshold}};
		my $rand_prob = rand($cum_probabilities[-1]);
		my $prev_prob = 0;
		for (my $j=0; $j<=$#cum_probabilities; $j++){
			if ($rand_prob>=$prev_prob && $rand_prob<$cum_probabilities[$j]) {
				if ($error_threshold eq 'substitution_threshold'){
					$subst_quot = $j;
				} elsif ($error_threshold eq 'indel_threshold'){
					$indel_quot = $j;
				}
				last;
			}
			$prev_prob = $cum_probabilities[$j];
		}
	}

	# If non-random error models are desired
	if (defined($error_models)){
		foreach my $model (shuffle keys %{$ERROR_MODELS->{$INP_tech}}){
			if ($read_seqs_generated->[$i] =~ /$model/){
				my $rand_freq = rand(100);
				if ($rand_freq < $ERROR_MODELS->{$INP_tech}{$model}{'freq'}) {
					my $error = $ERROR_MODELS->{$INP_tech}{$model}{'error'};
					$read_seqs_generated->[$i] =~ s/$model/$error/;
# 						printf("\n%s\n%s\n%s\n%s\n%s\n",$model,$error,$-[0],$kk,$read_seqs_generated->[$i]);
# 						exit;
					$count_errors->{'modelled-substitutions'}++;
# 					push(@$model_substitutions,$-[0]+1);
					push(@$model_substitutions,"$model->$error");
# 					$subst_quot--;
# 					if ($subst_quot <= 0) {
# 						last;
# 					}
				}
			}
		}
	}

	# Generates substitutions
	if ($subst_quot>=1) { # Inserts substitutions
		$count_errors->{'substitutions'}{$subst_quot}++;
		for (my $i=0; $i<$subst_quot; $i++){
			push(@$substitutions,int(rand($len+1)));
		}
	}

	# Generates indels
	if ($indel_quot>=1) { # Inserts indels
		$count_errors->{'indels'}{$indel_quot}++;
		for (my $i=0; $i<$indel_quot; $i++){
			if (int(rand(2)) == 0){ # int(rand(2)) generates 0 or 1
				push(@$insertions,int(rand($len+1)));
			} else {
				push(@$deletions,int(rand($len+1)));
			}
		}
	}
	
	$read_headers_generated->[$i] .= sprintf(" | md5=%s", generate_md5($read_seqs_generated->[$i]));
	if ($is_chimera){
		$read_headers_generated->[$i] .= sprintf(" | chimera");
	}
	if (@$model_substitutions){
		$read_headers_generated->[$i] .= sprintf(" | %s",join(',',@$model_substitutions));
	}
	if (@$substitutions || @$insertions || @$deletions) { # || @$model_substitutions) {
		$read_seqs_generated->[$i] = insert_sequence_errors($read_seqs_generated->[$i],$substitutions, $insertions, $deletions);
		$read_headers_generated->[$i] .= sprintf(" | errors=%s", print_sequence_errors($substitutions,$insertions,$deletions));
	} elsif (!$is_chimera) {
		$count_errors->{'noerror'}++;
	}

}
if ($verbose) {
	printf("\tTotal reads:\t\t%d\n", scalar @$read_seqs_generated);
	printf("\tWithout errors:\t\t%d\n", $count_errors->{'noerror'});
	if (defined($chimera_frequency)){
		printf("\tChimeras:\t\t%d\n", $count_errors->{'chimeras'});
	}
	if (defined($error_models)){
		printf("\tModelled substitutions:\t%d\n", $count_errors->{'modelled-substitutions'});
	}
	foreach my $n (sort {$a<=>$b} keys %{$count_errors->{'substitutions'}}) {
		printf("\tWith %d substitution/s:\t%d\n", $n, $count_errors->{'substitutions'}{$n});
	}
	foreach my $n (sort {$a<=>$b} keys %{$count_errors->{'indels'}}) {
		printf("\tWith %d indel/s:\t\t%d\n", $n, $count_errors->{'indels'}{$n});
	}

}



# Defines compression method
my $compression;
if (defined($INP_zip)){
	$compression = 'zip';
} elsif (defined($INP_gzip)){
	$compression = 'gzip';
}

my $outfile;
if (defined($read_seqs_generated)){
	$outfile = create_fasta_file($read_seqs_generated,$read_headers_generated,"$INP_outfile.fa",$compression);
}

if (defined($outfile)){
	printf("\nSaved %d simulated sequences into '%s'.\n\n", scalar @$read_seqs_generated, $outfile);
} else {
	print "\nThere was some error in the proccess and no sequences were generated.\n\n";
}




exit;
