# 🧬 tfa: Transcription Factor Annotator

tfa is a pipeline for high-confidence transcription factor (TF) annotation designed to work beyond traditional model organisms. It integrates sequence similarity, domain annotation, and deep learning predictions to identify and classify TFs from protein sequences.

📖 **Overview**

Understanding how phenotypes emerge from regulatory interactions is a central challenge in Systems Biology. This challenge becomes more pronounced when studying non-model organisms, which make up the vast majority of biodiversity but lack comprehensive functional annotation.

A key step in addressing this gap is the accurate identification of transcription factors (TFs). However, most existing TF databases are heavily biased toward model organisms, particularly vertebrates.

tfa was developed to:

Expand TF annotation to diverse organisms
Reduce the search space of candidate regulators
Provide reproducible and scalable TF predictions

⚙️ **Pipeline Overview**

tfa combines multiple complementary approaches to generate robust TF annotations:
- 🔍 **Sequence Similarity Search**
  - Uses DIAMOND to align input proteins
  - Retains candidate TF sequences

- 🧩 **Domain Annotation**
  - Uses InterProScan
  - Detects conserved DNA-binding domains

- 🤖 **Machine Learning Prediction**
  - Uses DeepTFactor
  - Predicts TF likelihood

- 🧹 **Filtering & Integration**
  - Merges results and applies filters

**Results are merged and filtered based on:**

- Prediction scores
- Removal of transposases
- Presence of known TF-associated domains

**Produces a consolidated results table including:**

- TF predictions
- Similarity metrics
- Domain annotations
- TF family classification

📦 **Requirements**

Make sure the following tools are installed:

- DIAMOND
- InterProScan
- DeepTFactor

🚀 **Installation**

Clone the repository:

```
git clone https://github.com/yourusername/tfa.git

cd tfa
```

Install dependencies (example using pip):

pip install -r requirements.txt
🧪 Usage
1. Create DIAMOND Database
```diamond makedb --in <uniprot_fasta> --db <uniprot_db>```
2. Run tfa
```
tfa.py \
  --input_fasta <protein_fasta_file> \
  --db_path <path_to_tf_database> \
  --deeptfactor_folder <deeptfactor_folder> \
  --tfs_domains <tf_domains_dataframe> \
  --output_prefix <output_prefix>
```
📂 Output

tfa generates:

Annotated TF prediction table
Domain annotation results
Similarity search metrics
TF family classifications

All outputs are prefixed using --output_prefix.

📁 Example Results and Benchmark Data

This repository includes a results/ folder containing example outputs generated during the development of the tfa pipeline. These datasets were used to guide parameter selection, validate filtering criteria, and evaluate the integration of similarity, domain, and machine learning predictions.

Providing these results serves two purposes:

Transparency and reproducibility
Users can inspect intermediate and final outputs to understand how predictions are generated.
Practical reference
Users can compare their own results against these examples to verify correct pipeline execution and expected output structure.

These files represent real use cases that motivated the development of tfa and illustrate its performance on biological datasets.

📚 References
Shen et al. (2023) – AnimalTFDB 4.0
Buchfink et al. (2015) – DIAMOND
Jones et al. (2014) – InterProScan
Bateman et al. (2004) – Pfam
Kim et al. (2021) – DeepTFactor

🤝 Contributing

Contributions are welcome! Please open an issue or submit a pull request.

📄 License

This project is licensed under the MIT License.

✨ Citation

If you use tfa in your research, please cite the relevant tools and databases listed above.
