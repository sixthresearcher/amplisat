
# AmpliSAT - Amplicon Sequence ASsignment Tools

## License

## Running online

Remember that it is not necessary to intall AmpliSAT, you can run it from the online server:
- http://evobiolab.biol.amu.edu.pl/amplisat/

## Running as a Docker container

If you need to automatize some pipelines or running some more customize commands,
then you can install [Docker](https://docs.docker.com/install/) in your server or computer
and run out-of-the-box the AmpliSAT Docker image from Docker Hub:
- `docker run -it sixthresearcher/amplisat`

## Installation in your local computer
If you cannot run AmpliSAT with the Docker image, here are the instructions to install it from scratch.

Clone from the Github repository and install later all the Perl dependencies:
1. Install git and few other utilities required:<br>
 `sudo apt-get install git wget curl unzip cpanminus expat libexpat1-dev`
2. Clone from the Github repo:<br>
 `git clone https://github.com/sixthresearcher/amplisat.git`
3. Install EMBOSS package:<br>
 `sudo apt-get install emboss`
4. Install required Perl modules:<br>
  ```
  cpanm File IO OLE Sort Sort::Naturally Archive::Zip Excel Excel::Writer::XLSX Spreadsheet Spreadsheet::XLSX
  cpanm Text Text::Iconv Net::SFTP::Foreign Statistics::Descriptive
  cpanm File::FindLib Inline::C threads
  cpanm Bio::Seq
  ```
5. Add to the PATH the folder where it was installed:<br>
 `export PATH=/amplisat-folder:$PATH`
6. Test that it is working properly:
 ``

Another way to install is to follow the instructions from here:
- http://evobiolab.biol.amu.edu.pl/amplisat/index.php?downloads


## Examples

1. Run AmpliCHECK to annotate the correct lenghts of your sequences and check which kind of sequencing errors do you have and in which proportion.<br>
   `perl ampliCHECK.pl -i guppy_example.fa.gz -d example.cvs -o amplicheck_guppy_example`
2. Execute AmpliSAS clustering (choose Illumina or 454/IonTorrent default parameters) specifying the correct length/s of your sequences (VERY IMPORTANT) and use in Filtering Parameters a 'Minimum per amplicon frequency (%)' of 3-5% (for a max of 10 expected alleles should work). Check results and compare then with AmpliCHECK ones and decide if you want a more strict filtering or if you have excluded real alleles that are very similar (you can tune the parameter 'Minimum frequency respect to the dominant').<br>
   `perl ampliSAS.pl -i guppy_example.fa.gz -d example.cvs -o amplisas_guppy_example -debug `
3. Compare AmpliSAS results obtained with different clustering/filtering parameters:<br>
   `perl ampliCOMPARE.pl -1 amplicheck_guppy_example/results.xlsx -2 amplisas_guppy_example/results.xlsx`
