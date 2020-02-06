
# AmpliSAT - Amplicon Sequence ASsignment Tools

## License

[See below paragraph](#License details)

## Documentation

[AmpliSAT manual](docs/amplisas_manual.pdf)

## Running online

Remember that it is not necessary to intall AmpliSAT, you can run it from the online server:
- http://evobiolab.biol.amu.edu.pl/amplisat/

## Running as a Docker container

If you need to automatize some pipelines or running some more customize commands,
then you can install [Docker](https://docs.docker.com/install/) in your server or computer
and run out-of-the-box the AmpliSAT Docker image from Docker Hub:
 `docker run -it sixthresearcher/amplisat`

## Installation in your local computer
If you cannot run AmpliSAT with the Docker image, here are the instructions to install it from scratch.

Clone from the Github repository and install all the Perl and system dependencies:
1. Install git and few other utilities required:<br>
   `sudo apt-get install git wget curl unzip cpanminus expat libexpat1-dev`
2. Go to the folder where you want to install AmpliSAT (recommended '/opt'):<br>
   `cd /opt`
3. Clone from the Github repo:<br>
   `git clone https://github.com/sixthresearcher/amplisat.git`
4. Install EMBOSS package:<br>
   `sudo apt-get install emboss`
5. Install required Perl modules:<br>
   ```
   cpanm File IO OLE Sort Sort::Naturally Archive::Zip Excel Excel::Writer::XLSX Spreadsheet Spreadsheet::XLSX
   cpanm Text Text::Iconv Net::SFTP::Foreign Statistics::Descriptive
   cpanm File::FindLib Inline::C threads
   cpanm Bio::Seq
   ```
6. Add to the PATH the folder where AmpliSAT was installed:<br>
   `export PATH=/amplisat-folder:$PATH`
7. Test that it is working properly:
   ``

Another way to install is to follow the instructions from here: http://evobiolab.biol.amu.edu.pl/amplisat/index.php?downloads

## Examples

1. Run AmpliCHECK to annotate the correct lenghts of your sequences and check which kind of sequencing errors do you have and in which proportion.<br>
   `perl ampliCHECK.pl -i guppy_example.fa.gz -d example.cvs -o amplicheck_guppy_example`
2. Execute AmpliSAS clustering (choose Illumina or 454/IonTorrent default parameters) specifying the correct length/s of your sequences (VERY IMPORTANT) and use in Filtering Parameters a 'Minimum per amplicon frequency (%)' of 3-5% (for a max of 10 expected alleles should work). Check results and compare then with AmpliCHECK ones and decide if you want a more strict filtering or if you have excluded real alleles that are very similar (you can tune the parameter 'Minimum frequency respect to the dominant').<br>
   `perl ampliSAS.pl -i guppy_example.fa.gz -d example.cvs -o amplisas_guppy_example -debug `
3. Compare AmpliSAS results obtained with different clustering/filtering parameters:<br>
   `perl ampliCOMPARE.pl -1 amplicheck_guppy_example/results.xlsx -2 amplisas_guppy_example/results.xlsx`

## License details

The AmpliSAT software has been developed partially as an academic project and partially as a private development.

- The scripts ampliSAS.pl, ampliCHECK.pl, ampliCLEAN.pl, ampliCOMBINE.pl, ampliCOMPARE.pl, ampliQC.pl, ampliLEGACY.pl,
ampliMERGE.pl, ampliMIX.pl, ampliSIM.pl, add_barcodes.pl and reformat_amplitaxo.pl have been developed
by Alvaro Sebastian, Magdalena Migalska and Jacek Radwan (Copyright 2015-2020) and they have the following license:
  Permission to use, copy, modify and distribute any part of this program for any educational, research, profit and non-profit purposes,
  without fee, and without a written agreement is hereby granted, provided that this paragraph and the following license paragraphs appear in all copies.
- The scripts ampliTCR.pl and ampliCDR3.pl have been developed
by Alvaro Sebastian, Magdalena Migalska and Jacek Radwan (Copyright 2015-2020) and they have the following license when used for the analysis of any data excluding human and mouse data:
  Permission to use, copy, modify and distribute any part of this program for any educational, research, profit and non-profit purposes,
  without fee, and without a written agreement is hereby granted, provided that this paragraph and the following license paragraphs appear in all copies.
- The scripts ampliTCR.pl and ampliCDR3.pl have been developed
by Alvaro Sebastian (Copyright 2015-2020) and they have the following license when used for the analysis of human or mouse data:
  Permission to use, copy, modify and distribute any part of this program for any educational, research and non-profit purposes, by non-profit institutions only,
  without fee, and without a written agreement is hereby granted, provided that this paragraph and the following license paragraphs appear in all copies.
- The scripts ampliCANCER.pl, ampliHLA.pl, ampliON.pl, ampliTAXO.pl, annotate_variants.pl, common_variants.pl, compare_variants.pl and filter_somatic_variants.pl have been developed
by Alvaro Sebastian (Copyright 2015-2020) and they have the following license when used for any kind of data:
  Permission to use, copy, modify and distribute any part of this program for any educational, research and non-profit purposes, by non-profit institutions only,
  without fee, and without a written agreement is hereby granted, provided that this paragraph and the following license paragraphs appear in all copies.

Those desiring to incorporate this work into commercial products or use for commercial purposes should contact Alvaro Sebastian using the following email address: sixthresearcher@gmail.com

IN NO EVENT SHALL THE INVENTORS (ALVARO SEBASTIAN, MAGDALENA MIGALSKA OR JACEK RADWAN) BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE INVENTORS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE SOFTWARE PROVIDED HEREIN IS ON AN “AS IS” BASIS, AND THE INVENTORS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. THE INVENTORS MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
