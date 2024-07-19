# TFA: Transcription Factor Annotator #

Bash and Python pipeline developed to annotate putative Transcription Factors (TFs) from protein sequence fasta files using complementary approaches. The pipeline uses popular softwares as Diamond, InterProScan and more recent methods as DeepTFactor.

The pipeline searches for putative TFs using two strategies. 
+ Search based on Homology: an AnimalTFDB4 database is constructed for DIAMOND blastp search. Diamond blastp search results are (default) filtered by an e-value of 1e-3 and one maximum target sequence for each query-sequence. Then, a new fasta file is created filtering sequences that presented a blast hit. This new filtered fasta is used as query for a interproscan search against PFAM database using InterProScan to seek if the sequence with a hit also presents the conserved domain representative of the TF family from which the sequence it had a match is part of. Lastly, sequences are filtered based on DeepTFactor that extracts protein features to classify protein sequences as probable TFs or not. This step results in a 
dataframe that contains query proteins identified by similarity and domain conservation and an extra column that's the result of DeepTFactor showing the probability of the protein being a TF based on protein sequence features (output_homologytfs.csv)

+ Search based only on protein features: The whole fasta file is used as input for DeepTFactor. This strategy can reveal putative TFs that would not be captured by homology in the other step. It returns a table with the query sequence and the plausability of the protein being a TF based on protein features (output_no_homolytfs.csv).

  
The diamond database is too heavy to add to github but it was constructed concatenating all fastas from AnimalTFDB4 and creating a diamond database using 'diamond make db'

Usage:
```
./tf_pipeline.sh <fasta_file> <path to atfdb4 diamond database>
```
