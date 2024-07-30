# tfa: transcription factor annotator #

TFa is  a pipeline developed for the annotation of putative TF sequences from protein sequence files.The pipeline searches for putative TFs using two strategies:

+ **homology-based search**: An AnimalTFDB4 database is constructed for DIAMOND blastp search. The results are filtered by an e-value of 1e-3 and limited to one maximum target sequence per query sequence. A new FASTA file is created with sequences that had a blast hit, which is then used for an InterProScan search against the PFAM database. This step checks if the hit sequences contain the conserved domain representative of the matched TF family. Finally, sequences are filtered using DeepTFactor, which classifies protein sequences based on features to determine their likelihood of being TFs. The resulting dataframe (output_homologytfs.csv) contains query proteins identified by similarity and domain conservation, along with a column showing the probability of the protein being a TF according to DeepTFactor.

+ **protein features-based search**:The entire FASTA file is used as input for DeepTFactor. This approach can identify potential TFs that might be missed by homology-based searches. The output is a table (output_no_homolytfs.csv) containing the query sequences and the probability of each protein being a TF based on its features.

  
Requirements:
+ [diamond](https://github.com/bbuchfink/diamond)  
+ [interproscan](https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html)
+ [deeptfactor](https://bitbucket.org/kaistsystemsbiology/deeptfactor/src/master/)

Usage:
```
./tf_pipeline.sh <fasta_file> <path to atfdb4 diamond database>
```
