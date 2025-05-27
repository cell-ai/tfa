# tfa: transcription factor annotator #

Understanding how phenotypes emerge from regulatory interactions is one of the central challenges in Systems Biology. This challenge becomes more daunting when studying non-model organisms, which make up the vast majority of the existing biodiversity. To address this problem, a crucial step forward is the comprehensive annotation of Transcription Factors (TFs). Such annotations can significantly refine the search space of variables and of potential gene regulators analyzed. However, most existing TF databases are heavily biased toward model organisms, particularly vertebrates. To bridge this gap, we developed a transcription factor annotator (tfa) designed to broaden the scope of TF annotation beyond traditional model systems.

**tfa is  a pipeline developed for the annotation of putative TF sequences from protein sequence files**.

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
```
```
#run tfa
tfa_v2.0.sh --input_fasta <protein_fasta_file> --db_path <path_to_atfdb_db> --deeptfactor_folder <deeptfactor_folder> --tfs_domains <dataframe with tfs domains> --output_prefix <prefix used in of output files>"
```
