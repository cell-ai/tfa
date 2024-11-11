# tfa: transcription factor annotator #

tfa is  a pipeline developed for the annotation of putative TF sequences from protein sequence files. The pipeline searches for putative TFs using two strategies:

+ **homology-based search**: An AnimalTFDB4 [1] database is constructed for DIAMOND blastp search [2]. The results are filtered by an e-value of 1e-3 and limited to one maximum target sequence per query sequence. A new FASTA file is created with sequences that had a blast hit, which is then used for an InterProScan search [3] against the PFAM database [4], using a reference list to consider the presence of TF-DNA binding domains (TFsdomains.csv). This step checks if the hit sequences contain the conserved domain representative of the matched TF family. Finally, sequences are filtered using DeepTFactor [5], which classifies protein sequences based on features to determine their likelihood of being TFs. The resulting dataframe (output_homologytfs.csv) contains query proteins identified by similarity and domain conservation, along with a column showing the probability of the protein being a TF according to DeepTFactor.

+ **protein features-based search**: The entire FASTA file is used as input for DeepTFactor [5]. This approach can identify potential TFs that might be missed by homology-based searches. The output is a table (output_no_homolytfs.csv) containing the query sequences and the probability of each protein being a TF based on its features.

The output is a dataframe (saved to file with the extension "tfa_results.csv"), which contains the information of what strategy was used to identify the protein ("homology" vs. "no homology"), the DeepTFactor score, and in the case of TFs identified by homology, additional information regarding protein, TF family and domains identified is also present.

**References:**

[1] Shen, W. K., Chen, S. Y., Gan, Z. Q., Zhang, Y. Z., Yue, T., Chen, M. M., ... & Guo, A. Y. (2023). AnimalTFDB 4.0: a comprehensive animal transcription factor database updated with variation and expression annotations. Nucleic acids research, 51(D1), D39-D45.

[2] Buchfink, B., Xie, C., & Huson, D. H. (2015). Fast and sensitive protein alignment using DIAMOND. Nature methods, 12(1), 59-60.

[3] Jones, P., Binns, D., Chang, H. Y., Fraser, M., Li, W., McAnulla, C., ... & Hunter, S. (2014). InterProScan 5: genome-scale protein function classification. Bioinformatics, 30(9), 1236-1240.

[4] Bateman, A., Coin, L., Durbin, R., Finn, R. D., Hollich, V., Griffiths‐Jones, S., ... & Eddy, S. R. (2004). The Pfam protein families database. Nucleic acids research, 32(suppl_1), D138-D141.

[5] Kim, G. B., Gao, Y., Palsson, B. O., & Lee, S. Y. (2021). DeepTFactor: A deep learning-based tool for the prediction of transcription factors. Proceedings of the National Academy of Sciences, 118(2), e2021171118.
  
**Requirements:**
+ [DIAMOND](https://github.com/bbuchfink/diamond)  
+ [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html)
+ [DeepTFactor](https://bitbucket.org/kaistsystemsbiology/deeptfactor/src/master/)

Usage:
```
#create database for similairty search
diamond makedb --in <fasta_file_animaltfdb4> --db <animaltfdb4_tf_db>
#run tfa
./tf_pipeline.sh <fasta_file> <animaltfdb4_tf_db>
```
