# amplisat
AmpliSAT - Amplicon Sequence ASsignment Tools

Examples of use:

1. Run AmpliCHECK to annotate the correct lenghts of your sequences and check which kind of sequencing errors do you have and in which proportion.

perl ampliCHECK.pl -i guppy_example.fa.gz -d example.cvs -o amplicheck_guppy_example

2. Execute AmpliSAS clustering (choose Illumina or 454/IonTorrent default parameters) specifying the correct length/s of your sequences (VERY IMPORTANT) and use in Filtering Parameters a 'Minimum per amplicon frequency (%)' of 3-5% (for a max of 10 expected alleles should work). Check results and compare then with AmpliCHECK ones and decide if you want a more strict filtering or if you have excluded real alleles that are very similar (you can tune the parameter 'Minimum frequency respect to the dominant').

perl ampliSAS.pl -i guppy_example.fa.gz -d example.cvs -o amplisas_guppy_example -debug 

3. Compare AmpliSAS results obtained with different clustering/filtering parameters:

perl ampliCOMPARE.pl -1 amplicheck_guppy_example/results.xlsx -2 amplisas_guppy_example/results.xlsx