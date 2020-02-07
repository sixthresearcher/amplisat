#!/usr/bin/perl -w
#
# Name: install_scripts.pl
#
# Version: 1.0
#
# License: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International
#          CC BY-NC-SA 4.0 (http://creativecommons.org/licenses/by-nc-sa/4.0/)
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
#   Checks and installs automatically the correct working of AmpliSAT scripts
#
# Requires test files in ./example folder
#
# Examples:
#
# Installing AmpliSAT scripts in the desired folder after checking their correct working:
#   perl install_scripts.pl -i ~/destination_folder -thr 20
#
# Installing after checking only a specific script:
#   perl install_scripts.pl -t amplitaxo -i ~/www/amplisat/bin/ -thr 20
#
# Checking only the correct working of a specific script indicating the folder with test files:
#   perl install_scripts.pl -t amplihla amplitaxo -e examples -thr 20
#
# Installing AmpliSAT scripts without testing:
#   perl install_scripts.pl -nt -i ~/www/amplisat/bin/
#
# Installing AmpliSAT scripts without testing (only advanced users: fast update of server scripts):
#   perl install_scripts.pl -nt -os -i ~/www/amplisat/bin/
#

my $VERSION = "1.0";
my $SCRIPT_NAME = 'install_scripts.pl';
my $AUTHOR = "Alvaro Sebastian";
my $DESCRIPTION = "Test automatically the correct working of AmpliSAT scripts.";

use Getopt::Long;
use File::Basename;

# All variables must be declared before their use
use strict;
# Turn autoflush on
local $| = 1;

my $COMMAND_LINE = $0." ".join(" ",@ARGV);

# Default options
# Tools to test
my @AMPLISAT_TOOLS = ('amplisas','amplicheck','amplilegacy','amplitcr','amplicdr3','amplimerge','ampliclean','amplimix','amplicompare','amplicombine', 'amplisim', 'install_scripts', 'reformat_amplitaxo'); # 'add_barcodes','ampliqc');
my @OPTIONAL_TOOLS = ('amplicancer','amplitaxo','amplihla' );
# Data examples location
my $EXAMPLE_DIR = dirname(__FILE__).'/examples';
# Routes to binary files
my $EXTERNAL_TOOLS_DIR = './lib/Bio/tools';
# Routes to required additional files
my $ADDITIONAL_FILES = {
	'amplihla'  => [ 'imgt/*.fa' ],
	'amplitaxo' => [ 'taxo/*.utax.fa.gz' ],
	'amplicdr3' => [ 'tcrefs/*.fna' ],
	# 'amplicancer' => [ 'human/*' ],
};

my ($INP_install,$INP_test,@INP_tools,$INP_example_folder,$INP_threads,$INP_onlyscripts,$INP_notest);

GetOptions(
	'h|help|?' =>  \&usage,
	'i|install=s' => \$INP_install,
	'test' => \$INP_test,
	't|tool=s{,}' => \@INP_tools,
	'e|examples=s' => \$INP_example_folder,
	'thr|threads=i' => \$INP_threads,
	'os|onlyscripts' => \$INP_onlyscripts,
	'nt|notest' => \$INP_notest,
);

# Usage help
sub usage
{
	print "\n$SCRIPT_NAME version $VERSION by $AUTHOR\n";
	print "\n$DESCRIPTION\n";
	print "\nUsage: ";
	print "$SCRIPT_NAME -i <file> -d <file> [options]\n";
	print "\nOptions:\n";
	print "  -i <folder>\tDestination folder to install the scripts if they pass the tests.\n";
	print "  -test \tOnly test the correct working of the scripts without installing.\n";
	print "  -t <tool1> [ <tool2> ... <toolN> ]\n\t\tTools to install/test ('".join("','",@AMPLISAT_TOOLS)."').\n\t\tDefault: ALL\n";
	print "  -e <folder>\tFolder containing the example data (Default: '$EXAMPLE_DIR')\n";
	print "  -thr <number>\tNumber of threads to calculate the alignments.\n";
	print "  -h\t\tHelp.\n";
# 	print "  -os\t\tInstall only perl scripts.\n";
	print "  -nt\t\tSkip tests.\n";
	print "\n";
	exit;
}

if (!defined($INP_install) && !defined($INP_test)){
	print "\nERROR: Please specify the destination folder to install the scripts.\n\n";
	usage();
	exit;
} elsif (defined($INP_install)){
	if (-d $INP_install){
		print "\nWARNING: Destination folder '$INP_install' already exists.\n";
		print "Do you want to overwrite its contents? (yes,no): ";
		my $option = <STDIN>;
		chomp $option;
		if ($option ne 'y' && $option ne 'yes') {
			print "\nERROR: Scripts couldn't be installed.\n\n";
			exit;
		}
	}
}
if (@INP_tools) {
	@INP_tools = map lc($_), @INP_tools;
} else {
	@INP_tools = @AMPLISAT_TOOLS;
}
if (defined($INP_example_folder) && -d $INP_example_folder){
	$EXAMPLE_DIR = $INP_example_folder;
}
if (!defined($INP_threads)){
	$INP_threads = 1;
}

# Data examples files
my $EXAMPLE_DATA = {
	'amplisas' => {
		'reads_file' => $EXAMPLE_DIR."/amplisas_example.fq.gz",
		'amplicons_file' => $EXAMPLE_DIR."/amplisas_example.csv",
		'alleles_file' => $EXAMPLE_DIR."/amplisas_example.alleles.fa",
		'results_file' => $EXAMPLE_DIR."/amplisas_example_amplisas.xlsx",
	},
	'amplimerge' => {
		'reads1_file' => $EXAMPLE_DIR."/amplisas_example_R1.fq.gz",
		'reads2_file' => $EXAMPLE_DIR."/amplisas_example_R2.fq.gz",
	},
	'amplihla' => {
		#'reads_file' => $EXAMPLE_DIR."/amplihla_example.zip",
		'reads_file' => $EXAMPLE_DIR."/amplihla_example.fq.gz",
		'amplicons_file' => $EXAMPLE_DIR."/amplihla_example.csv",
# 		'alleles_file' => $EXAMPLE_DIR."/amplihla_example.alleles.fa",
# 		'results_file' => $EXAMPLE_DIR."/amplihla_example.xlsx",
	},
	'amplitaxo' => {
		'reads_file' => $EXAMPLE_DIR."/amplitaxo_example.zip",
		'db_file' => $EXAMPLE_DIR."/amplitaxo_example.db.fa.gz",
# 		'reads_file' => $EXAMPLE_DIR."/amplitaxo_example.fq.gz",
# 		'amplicons_file' => $EXAMPLE_DIR."/amplitaxo_example.csv",
# 		'db_file' => $EXAMPLE_DIR."/amplitaxo_example.db.fa.gz",
# 		'results_file' => $EXAMPLE_DIR."/amplitaxo_example.xlsx",
	},
	'amplicancer' => {
		'reads_file' => $EXAMPLE_DIR."/amplicancer_example.fq.gz",
	},
	'amplitcr' => {
		'reads_file' => $EXAMPLE_DIR."/amplitcr_example.fq.gz",
	},
	'amplicdr3' => {
		'reads_file' => $EXAMPLE_DIR."/amplicdr3_example.zip",
		'amplicons_file' => $EXAMPLE_DIR."/amplicdr3_example.csv",
	},
	'amplilegacy' => {
		'reads2_file' => $EXAMPLE_DIR."/amplisas_example2.fq.gz",
		'amplicons2_file' => $EXAMPLE_DIR."/amplisas_example2.csv",
	},
	'amplimix' => {
		'reads2_file' => $EXAMPLE_DIR."/amplisas_example2.fq.gz",
		'amplicons2_file' => $EXAMPLE_DIR."/amplisas_example2.csv",
	},
	'amplicompare' => {
		'results2_file' => $EXAMPLE_DIR."/amplisas_example2_amplisas.xlsx",
	},
	'amplicombine' => {
		'results2_file' => $EXAMPLE_DIR."/amplisas_example2_amplisas.xlsx",
	},
	'amplisim' => {},,
};

if (!-d $EXAMPLE_DIR){
	print "\nERROR: '$EXAMPLE_DIR' doesn't contain example data.\n\n";
	usage();
	exit;
}

foreach my $tool (@INP_tools) {
	if (!in_array([@AMPLISAT_TOOLS,@OPTIONAL_TOOLS], $tool)){
		next;
	}
	if (defined($EXAMPLE_DATA->{$tool})){
		foreach my $example_file (keys %{$EXAMPLE_DATA->{$tool}}){
			if (!-e $EXAMPLE_DATA->{$tool}{$example_file} && !in_array(\@OPTIONAL_TOOLS,$tool)){
				printf("\nERROR: Example file '%s' doesn't exist.\n\n", $EXAMPLE_DATA->{$tool}{$example_file});
				exit;
			}
		}
	}
}

if (!defined($INP_notest)){
	if (!-d $EXTERNAL_TOOLS_DIR){
		print "\nERROR: '$EXTERNAL_TOOLS_DIR' folder couldn't be found. Please check that all AmpliSAT ZIP package contents are correctly extracted.\n\n";
		exit;
	}
	my $EXTERNAL_TOOLS = {
		"BLASTN" => { "command" => "blastn -version", "output" => "blastn: 2\." },
		"MAKEBLASTDB" => { "command" => "makeblastdb -version", "output" => "makeblastdb: 2\." },
		"MAFFT" => { "command" => "mafft --version", "output" => "v7\." },
		"FLASH" => { "command" => "flash2 -version", "output" => "FLASH v2\." },
		"NEEDLE" => { "command" => "needle -version", "output" => "EMBOSS:6.6\." },
		"NEEDLEALL" => { "command" => "needleall -version", "output" => "EMBOSS:6.6\." },
	};
	printf("\nChecking external dependencies.\n");
	foreach my $tool (sort {$a cmp $b} keys %$EXTERNAL_TOOLS) {
		if ($tool !~ /needle/i) {
			my $command = $EXTERNAL_TOOLS_DIR.'/'.$EXTERNAL_TOOLS->{$tool}{'command'};
			my $output = `$command 2>&1`;
			my $expected_output = $EXTERNAL_TOOLS->{$tool}{'output'};
			if (!defined($output) || $output !~ /$expected_output/){
				print "\nERROR: '$tool' couldn't be found in '$EXTERNAL_TOOLS_DIR' folder. Please check that all AmpliSAT ZIP package contents are correctly extracted and that your computer is a 64 bits Linux system.\n\n";
				exit;
			}
		} else {
			my $command = $EXTERNAL_TOOLS->{$tool}{'command'};
			my $output = `$command 2>&1`;
			my $expected_output = $EXTERNAL_TOOLS->{$tool}{'output'};
			if (!defined($output) || $output !~ /$expected_output/){
				print "\nERROR: EMBOSS Suite 6.6 or newer is required. Install in Ubuntu with the command: 'sudo apt-get install emboss'.\n\n";
				exit;
			}
		}
		printf("\t'%s' found and working.\n", $tool);
	}
	
	my @PERL_MODULES = ( 'File::FindLib' );
	printf("\nChecking Perl modules.\n");
	foreach my $module (sort {$a cmp $b} @PERL_MODULES) {
		my $output = `perldoc -l $module`;
		if (!defined($output) || $output =~ /No documentation found/){
			print "\nERROR: '$module' Perl module couldn't be found. Install in Ubuntu with the command: 'sudo cpan -i $module'.\n\n";
		}
		printf("\t'%s' found and working.\n", $module);
	}

}

my (%tool_names, %script_names);
foreach my $tool (@INP_tools) {

	if (!in_array([@AMPLISAT_TOOLS,@OPTIONAL_TOOLS], $tool)){
		next;
	}

	my $tool_name = $tool;
	my $script_name = $tool.".pl";
	if ($tool =~ /^(ampli)(.+)/i){
		$tool_name = ucfirst($1).uc($2);
		$script_name = $1.uc($2).".pl";
	}
	$tool_names{$tool} = $tool_name;
	$script_names{$tool} = $script_name;

	if (defined($INP_notest)){
		next;
	}

	my $output_path = $tool."_out";
	# my $output_path = "/tmp/".random_file_name();

	my $command;
	if ($tool eq 'amplicheck' || $tool eq 'amplisas'){
		$command = sprintf("./%s -i %s -d %s -o %s -thr %d -a %s", $script_name, $EXAMPLE_DATA->{'amplisas'}{'reads_file'}, $EXAMPLE_DATA->{'amplisas'}{'amplicons_file'}, $output_path, $INP_threads, $EXAMPLE_DATA->{'amplisas'}{'alleles_file'});
	} elsif ($tool eq 'ampliqc'){
		$command = sprintf("./%s -i %s -d %s -o %s -thr %d", $script_name, $EXAMPLE_DATA->{'amplisas'}{'reads_file'}, $EXAMPLE_DATA->{'amplisas'}{'amplicons_file'}, $output_path, $INP_threads);
	} elsif ($tool eq 'ampliclean'){
		$command = sprintf("./%s -i %s -d %s -o %s -thr %d -gz", $script_name, $EXAMPLE_DATA->{'amplisas'}{'reads_file'}, $EXAMPLE_DATA->{'amplisas'}{'amplicons_file'}, $output_path, $INP_threads);
	} elsif ($tool eq 'amplihla'){
		# $command = sprintf("./%s -i %s -o %s -thr %d -na 500", $script_name, $EXAMPLE_DATA->{'amplihla'}{'reads_file'}, $output_path, $INP_threads);
		$command = sprintf("./%s -i %s -d %s -o %s -thr %d -na 500", $script_name, $EXAMPLE_DATA->{'amplihla'}{'reads_file'}, $EXAMPLE_DATA->{'amplisas'}{'amplicons_file'}, $output_path, $INP_threads);
	} elsif ($tool eq 'amplitaxo'){
		$command = sprintf("./%s -i %s -o %s -thr %d -na 500 -a %s", $script_name, $EXAMPLE_DATA->{'amplitaxo'}{'reads_file'}, $output_path, $INP_threads, $EXAMPLE_DATA->{'amplitaxo'}{'db_file'});
		# Following command is not working because AmpliTAXO database UTAX files are not downloaded (bif files):
		# $command = sprintf("./%s -i %s -d %s -o %s -thr %d -na 500 -db greengenes", $script_name, $EXAMPLE_DATA->{'amplitaxo'}{'reads_file'}, $EXAMPLE_DATA->{'amplitaxo'}{'amplicons_file'}, $output_path, $INP_threads);
	} elsif ($tool eq 'amplicancer'){
		$command = sprintf("./%s -i %s -o %s -thr %d -d panel", $script_name, $EXAMPLE_DATA->{'amplicancer'}{'reads_file'}, $output_path, $INP_threads);
	} elsif ($tool eq 'amplitcr'){
		$command = sprintf("./%s -i %s -o %s -vpat '%s' -jpat '%s' -cpat '%s' -thr %d", $script_name, $EXAMPLE_DATA->{'amplitcr'}{'reads_file'}, $output_path, 'Qx[PS]x(14)Cx(10,11)WYx(39,42)[LM]x(14)C', 'GxGx(2)Lx[VI]', 'EDL*', $INP_threads);
	} elsif ($tool eq 'amplicdr3'){
		$command = sprintf("./%s -s vole -i %s -d %s -o %s -thr %d", $script_name, $EXAMPLE_DATA->{'amplicdr3'}{'reads_file'}, $EXAMPLE_DATA->{'amplicdr3'}{'amplicons_file'}, $output_path, $INP_threads);
	} elsif ($tool eq 'amplilegacy'){
		$command = sprintf("./%s -i %s -d %s -o %s -thr %d -a %s -m lighten", $script_name, $EXAMPLE_DATA->{'amplisas'}{'reads_file'}, $EXAMPLE_DATA->{'amplisas'}{'amplicons_file'}, $output_path, $INP_threads, $EXAMPLE_DATA->{'amplisas'}{'alleles_file'});
	} elsif ($tool eq 'amplimerge'){
		$command = sprintf("./%s -i %s %s -o %s -thr %d -gz", $script_name, $EXAMPLE_DATA->{'amplimerge'}{'reads1_file'}, $EXAMPLE_DATA->{'amplimerge'}{'reads2_file'}, $output_path, $INP_threads);
	} elsif ($tool eq 'amplimix'){
		$command = sprintf("./%s -i %s %s -d %s %s -o %s -gz", $script_name, $EXAMPLE_DATA->{'amplisas'}{'reads_file'}, $EXAMPLE_DATA->{'amplimix'}{'reads2_file'}, $EXAMPLE_DATA->{'amplisas'}{'amplicons_file'}, $EXAMPLE_DATA->{'amplimix'}{'amplicons2_file'}, $output_path);
	} elsif ($tool eq 'amplicompare'){
		$command = sprintf("./%s -1 %s -2 %s -o %s", $script_name, $EXAMPLE_DATA->{'amplisas'}{'results_file'}, $EXAMPLE_DATA->{'amplicompare'}{'results2_file'}, $output_path);
	} elsif ($tool eq 'amplicombine'){
		$command = sprintf("./%s -i %s %s -o %s", $script_name, $EXAMPLE_DATA->{'amplisas'}{'results_file'}, $EXAMPLE_DATA->{'amplicombine'}{'results2_file'}, $output_path);
	} elsif ($tool eq 'amplisim'){
		$command = sprintf("./%s -d %s -o %s -na 500 -t illumina -gz", $script_name, $EXAMPLE_DATA->{'amplisas'}{'amplicons_file'}, $output_path);
	} else {
		printf("\n'%s' test skipped.\n", $tool_name);
		next
	}

	printf("\nTesting '%s'.\n", $tool_name);
	printf("\tRunning '%s'.\n", $command);
	# exit;
	my $tool_output = `$command`;

	my $errors = 0;
	if ($tool_output =~ /(ERROR:?\s*.+)/) {
		print "\t$1\n\n";
		$errors = 1;
	} elsif (in_array(['amplisas','amplicheck','amplilegacy','amplihla','amplitaxo'],$tool)){
		if (!-d $output_path){
			printf("\tERROR: Not detected output folder.\n\n");
			$errors = 1;
		} elsif (!-e "$output_path/results.xlsx"){
			printf("\tERROR: Not detected results file.\n\n");
			$errors = 1;
		} elsif ($tool_output =~ /results (stored|written) into '.+'/){
			printf("\t%s test finished CORRECTLY.\n", $tool_name);
		} else {
			printf("\n%s test finished with ERRORS, check the ouput:\n", $tool_name);
			print $tool_output;
			$errors = 1;
		}
	} elsif (in_array(['amplicompare','amplicombine'],$tool)){
		if (!-e "$output_path.xlsx"){
			printf("\tERROR: Not detected results file.\n\n");
			$errors = 1;
		} elsif ($tool_output =~ /results (stored|written) into '.+'/){
			printf("\t%s test finished CORRECTLY.\n", $tool_name);
		} else {
			printf("\n%s test finished with ERRORS, check the ouput:\n", $tool_name);
			print $tool_output;
			$errors = 1;
		}
	} elsif (in_array(['amplimerge','ampliclean','amplimix','amplisim'],$tool)){
		if (!-e "$output_path.fq.gz" && !-e "$output_path.fa.gz"){
			printf("\tERROR: Not detected output file.\n\n");
			$errors = 1;
		} elsif (in_array(['amplimix'],$tool) && !-e "$output_path.csv"){
			printf("\tERROR: Not detected output file.\n\n");
			$errors = 1;
		} elsif ($tool_output =~ /\d+ (merged|concatenated|cleaned|extracted|mixed|simulated)? sequences into '.+'|Saved \d+ sequences from \d+ into '.+'/){
			printf("\t%s test finished CORRECTLY.\n", $tool_name);
		} else {
			printf("\n%s test finished with ERRORS, check the ouput:\n", $tool_name);
			print $tool_output;
			$errors = 1;
		}
	} elsif (in_array(['amplitcr'],$tool)){
		if (!-e "$output_path/amplitcr_example.tcr.fa.gz"){
			printf("\tERROR: Not detected output file.\n\n");
			$errors = 1;
		} elsif ($tool_output =~ /results (stored|written) into '.+'/){
			printf("\t%s test finished CORRECTLY.\n", $tool_name);
		} else {
			printf("\n%s test finished with ERRORS, check the ouput:\n", $tool_name);
			print $tool_output;
			$errors = 1;
		}
	} elsif (in_array(['amplicdr3'],$tool)){
		if (!-d $output_path){
			printf("\tERROR: Not detected output folder.\n\n");
			$errors = 1;
		} elsif (!-e "$output_path/amplicdr3_example.stats.xlsx"){
			printf("\tERROR: Not detected results file.\n\n");
			$errors = 1;
		} elsif ($tool_output =~ /results (stored|written) into '.+'/){
			printf("\t%s test finished CORRECTLY.\n", $tool_name);
		} else {
			printf("\n%s test finished with ERRORS, check the ouput:\n", $tool_name);
			print $tool_output;
			$errors = 1;
		}
	} elsif (in_array(['ampliqc'],$tool)){
		if (!-e "$output_path.xlsx"){
			printf("\tERROR: Not detected results file.\n\n");
			$errors = 1;
		} elsif ($tool_output =~ /printed in file '.+'/){
			printf("\t%s test finished CORRECTLY.\n", $tool_name);
		} else {
			printf("\n%s test finished with ERRORS, check the ouput:\n", $tool_name);
			print $tool_output;
			$errors = 1;
		}
	} elsif (in_array(['amplicancer'],$tool)){
		if ($tool_output =~ /Formating and printing mutation information/){
			printf("\t%s test finished CORRECTLY.\n", $tool_name);
		} else {
			printf("\n%s test finished with ERRORS, check the ouput:\n", $tool_name);
			print $tool_output;
			$errors = 1;
		}
	}
	`rm -rf $output_path*`;
	if ($errors) {
		exit;
	}

}

if (!defined($INP_notest)){
	printf("\nAll tested tools work CORRECTLY.\n\n");
}

if (defined($INP_install)){

	printf("\nInstalling the scripts into '$INP_install'.\n\n");
	if (!-d $INP_install){
		mkdir($INP_install);
	}

	if (!defined($INP_onlyscripts)){
		`cp -frp lib $INP_install/`;
		`cp -frp $EXAMPLE_DIR $INP_install/`;
	} else {
		if (!-d "$INP_install/lib"){
			mkdir("$INP_install/lib");
		}
		if (!-d "$INP_install/lib/Bio"){
			mkdir("$INP_install/lib/Bio");
		}
		`cp -fp lib/Bio/Ampli.pm lib/Bio/Sequences.pm $INP_install/lib/Bio/`;
		`cp -fp $SCRIPT_NAME $INP_install/`;
		if (in_array(\@INP_tools,'amplicancer') || in_array(\@INP_tools,'amplihla')){
			`cp -fp lib/Bio/Onco.pm $INP_install/lib/Bio/`;
		}
	}

	foreach my $tool (@INP_tools) {
		my $install_command = sprintf("cp -f %s %s/", $script_names{$tool}, $INP_install);
		`$install_command`;
		printf("\t'%s' installed.\n", $tool_names{$tool});
		# Install specific requirements
		if (!defined($INP_onlyscripts) && $tool eq 'amplihla'){
			if (!-d "$INP_install/imgt"){
				mkdir("$INP_install/imgt");
			}
# 			print "cp -frp ".join(" ",@{$ADDITIONAL_FILES->{'amplihla'}})." $INP_install/imgt/\n";
			system("cp -frp ".join(" ",@{$ADDITIONAL_FILES->{'amplihla'}})." $INP_install/imgt/");
# 			`cp -frp imgt $INP_install/`;
		} elsif (!defined($INP_onlyscripts) && $tool eq 'amplitaxo'){
			if (!-d "$INP_install/taxo"){
				mkdir("$INP_install/taxo");
			}
			system("cp -frp ".join(" ",@{$ADDITIONAL_FILES->{'amplitaxo'}})." $INP_install/taxo/");
		} elsif (!defined($INP_onlyscripts) && $tool eq 'amplicdr3'){
			if (!-d "$INP_install/tcrefs"){
				mkdir("$INP_install/tcrefs");
			}
			system("cp -frp ".join(" ",@{$ADDITIONAL_FILES->{'amplicdr3'}})." $INP_install/tcrefs/");
		}
# 		} elsif (!defined($INP_onlyscripts) && $tool eq 'amplicancer'){
# 			if (!-d "$INP_install/human"){
# 				mkdir("$INP_install/human");
# 			}
# 			system("cp -frp ".join(" ",@{$ADDITIONAL_FILES->{'amplicancer'}})." $INP_install/human/");
# 		}
	}
	
	printf("\nThe scripts have been CORRECTLY installed into '$INP_install'.\n\n");

}


exit;


#################################################################################

# Searches a string inside an array
# Example: in_array($sentence, "/$word/", 1);
sub in_array {

	my ($array,$search_for,$search_posic) = @_;

	my @posic;

	# If search a pattern
	if ($search_for =~ /^\/(.+)\/$/){
		for (my $i=0; $i<=$#{$array}; $i++) {
			if (defined($array->[$i]) && $array->[$i] =~ /$1/i){
				push(@posic, $i);
			}
		}
	} else {
		for (my $i=0; $i<=$#{$array}; $i++) {
			if (defined($array->[$i]) && $array->[$i] eq $search_for){
				push(@posic, $i);
			}
		}
	}

	if (defined($search_posic) && $search_posic){
		(scalar @posic)?return (1,\@posic):return (0);
	} else {
		(scalar @posic)?return 1:return 0;
	}
}


#################################################################################

# Creates a random file name
sub random_file_name {
	my ($length) = @_;

	if (!defined($length)){
		$length = 15;
	}

	my $max_number = "9"x$length;

	return sprintf("%0".$length."u",int(rand($max_number)));
}

#################################################################################
