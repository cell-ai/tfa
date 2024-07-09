#!/bin/bash

# Check for input fasta file and database path
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <protein_fasta_file> <path_to_atfdb_db>"
    exit 1
fi

# Assign input file and database path to variables
input_fasta=$1
db_path=$2

# Define output files
diamond_output="$(basename ${input_fasta%.faa}).atfdb.1e3.txt"
filtered_fasta="$(basename ${input_fasta%.faa}).filtered_sequences.fasta"
interproscan_output="$(basename ${input_fasta%.faa}).ips_results.tsv"

# Step 1: Run diamond blastp
echo "Running diamond blastp..."
diamond blastp --query "$input_fasta" --max-target-seqs 1 --evalue 1e-3 --db "$db_path" --outfmt 6 qseqid qlen sseqid slen salltitles pident length mismatch gapopen qstart qend sstart send qcovhsp evalue bitscore --out "$diamond_output" --threads 50 --sensitive

# Check if diamond ran successfully
if [ $? -ne 0 ]; then
    echo "Error running diamond blastp"
    exit 1
fi

#Extract putative tf headers
cut -f 1 "$diamond_output" > "$(basename ${input_fasta%.faa}).tfheaders.txt"
#Filter original fasta file using seqtk and the list of tf headers
seqtk subseq "$input_fasta" "$(basename ${input_fasta%.faa}).tfheaders.txt" > "$filtered_fasta"

# Check if the filtering was successful
if [ $? -ne 0 ]; then
    echo "Error filtering fasta file"
    exit 1
fi

# Step 3: Run interproscan on the filtered fasta file
echo "Running interproscan..."
interproscan -appl Pfam -cpu 50 -f TSV -dra -iprlookup -i "$filtered_fasta" -o "$interproscan_output"

# Check if interproscan ran successfully
if [ $? -ne 0 ]; then
    echo "Error running interproscan"
    exit 1
fi
original_dir=$(pwd)
# Step 4: Run DeepTFactor
echo "Running DeepTFactor..."
original_dir=$(pwd)
og_fasta="$original_dir/$input_fasta"

cd /project/yutaka/spiderNetwork/local/deeptfactor || exit

python tf_running.py -cpu 50 -i "$og_fasta" -o "$original_dir/$(basename ${input_fasta%.faa}).deeptf.res"

# Check if deeptfactor ran successfully
if [ $? -ne 0 ]; then
    echo "Error running DeepTFactor"
    exit 1
fi

cd "$original_dir" || exit

# Step 5: Run Python analysis script
echo "Running Python analysis script..."
python analyze_ips_diam.py --diamond "$(basename ${input_fasta%.faa}).atfdb.1e3.txt" --ips "$(basename ${input_fasta%.faa}).ips_results.tsv" --output $(basename ${input_fasta%.faa})

python merge_all_res.py --diamond_ips_res $(basename ${input_fasta%.faa}).diamips.tfs.csv --deeptf_res "$original_dir/$(basename ${input_fasta%.faa}).deeptf.res/prediction_result.txt" --output $(basename ${input_fasta%.faa})

# Check if the Python script ran successfully
if [ $? -ne 0 ]; then
    echo "Error running Python analysis script"
    exit 1
fi

echo "Pipeline completed successfully"
