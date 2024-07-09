# tff: Transcription Factor Finder
Pipeline developed to identify Transcription Factors (TFs) from protein fasta file using different strategies.

The scripts were created in order to identify TFs from protein sequence fasta files and hopefully can aid in the identification of a more reliable set of TFs specially for non-model organisms. First, an AnimalTFDB4 database is constructed for DIAMOND blastp search. Diamond blastp search results are (default) filtered by an e-value of 1e-3 and one maximum target sequence for each query-sequence. Then, a new fasta file is created filtering sequences that presented a blast hit. This new filtered fasta is used as query for a interproscan search against PFAM database using InterProScan to seek if the sequence with a hit also presents the conserved domain representative of the TF family from which the sequence it had a match is part of. Finally, sequences are filtered based on DeepTFactor that extracts protein features to classify protein sequences as probable TFs or not.

The diamond database is too heavy to add to github but it was constructed concatenating all fastas from AnimalTFDB4 and creating a diamond database using 'diamond make db'


Usage:
```
./tf_pipeline.sh <fasta_file> <path to atfdb4 diamond database>
```
