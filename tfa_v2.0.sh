#!/bin/bash

#
# ============================================
#  Program Name : tfa_v2.0.sh
#  Description  : Transcription Factor annotation based on amino acid sequences.
#  Version      : 2.0.0
#  Author       : Marcela Akemi Ishihara
#  Language     : Bash (GNU/Linux)
#  Created On   : 2025-05-27
#  Last Updated : 2025-11-20
#
#  Usage        : ./tfa_v2.0.sh [options]
# ============================================
#

##########-------------##########-------------##########-------------##########-------------##########-------------
#Step 0: Initialize and parameter check
# Initialize variables
input_fasta=""
db_path=""
deeptfactor_folder=""
tfs_domains=""

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_fasta) input_fasta="$2"; shift ;;
        --db_path) db_path="$2"; shift ;;
        --deeptfactor_folder) deeptfactor_folder="$2"; shift ;;
	--tfs_domains) tfs_domains="$2"; shift ;;
	--output_prefix) output_prefix="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check if required arguments are provided
if [[ -z "$input_fasta" || -z "$db_path" || -z "$deeptfactor_folder" || -z "$tfs_domains" || -z "$output_prefix" ]]; then
    echo "Usage: $0 --input_fasta <protein_fasta_file> --db_path <path_to_atfdb_db> --deeptfactor_folder <deeptfactor_folder> --tfs_domains <dataframe with tfs domains> --output_prefix <prefix used in of output files>"
    exit 1
fi

# Print assigned variables
echo "Input FASTA file: $input_fasta"
echo "TFS domains csv: $tfs_domains"
echo "DeepTFactor folder: $deeptfactor_folder"
echo "Database path: $db_path"
echo "Saved files with prefix: $output_prefix"

# Define output files
deeptfactor_output="${output_prefix}.deeptf.res/prediction_result.txt"
diamond_output="${output_prefix}.diamond.1e3.txt"
filtered_fasta="${output_prefix}.filtered_sequences.fasta"
interproscan_output="${output_prefix}.ips_results.tsv"

##########-------------##########-------------##########-------------##########-------------##########-------------
# Step 1: Run DeepTFactor (evaluate the probability of a sequence representing a putative TF)
echo "Running DeepTFactor..."
original_dir=$(pwd)
og_fasta="$original_dir/$input_fasta"

#%todo: fix deeptfactor installation so we can skip this step
cd "$deeptfactor_folder" || exit

python tf_running.py -cpu 50 -i "$og_fasta" -o "$original_dir/${output_prefix}.deeptf.res"

# Step 1.2: Check if deeptfactor ran successfully
if [ $? -ne 0 ]; then
    echo "Error running DeepTFactor"
    exit 1
fi

cd "$original_dir" || exit

# Step 1.3: Extract putative tf headers with score > 0.5
awk -v column_number=3 -v threshold=0.5 '$column_number > threshold' ${deeptfactor_output} | cut -f 1 | tail -n +2 > "${output_prefix}.tfheaders.txt"

# Step 1.4: Filter original fasta file using seqtk and the list of tf headers
seqtk subseq "$input_fasta" "${output_prefix}.tfheaders.txt" > "$filtered_fasta"

# Step 1.5: Check if the filtering was successful
if [ $? -ne 0 ]; then
    echo "Error filtering fasta file"
    exit 1
fi
# Step 1.6: Check if the filtered fasta file is empty
if [ ! -s "$filtered_fasta" ]; then
    echo "Filtered fasta file is empty. No putative TFs found."
    exit 1
fi

echo "Filtered fasta of putative TFs created: $filtered_fasta"

##########-------------##########-------------##########-------------##########-------------##########-------------
# Step 2: Run diamond blastp
echo "Running diamond blastp..."

diamond blastp --query "$filtered_fasta" --max-target-seqs 1 --evalue 1e-3 --db "$db_path" --outfmt 6 qseqid qlen sseqid slen salltitles pident length mismatch gapopen qstart qend sstart send qcovhsp evalue bitscore --out "$diamond_output" --threads 50 --sensitive

# Step 2.1: Check if diamond ran successfully
if [ $? -ne 0 ]; then
    echo "Error running diamond blastp"
    exit 1
fi

##########-------------##########-------------##########-------------##########-------------##########-------------
# Step 3: Run interproscan on the filtered fasta file
echo "Running interproscan..."
interproscan -appl Pfam -cpu 50 -f TSV -dra -iprlookup -i "$filtered_fasta" -o "$interproscan_output"

# Step 3.1: Check if interproscan ran successfully
if [ $? -ne 0 ]; then
    echo "Error running interproscan"
    exit 1
fi
original_dir=$(pwd)

# Step 3.2: Check if interproscan output is empty
if [ ! -s "$interproscan_output" ]; then
    echo "InterProScan output is empty. No domains found."
    exit 1
fi

##########-------------##########-------------##########-------------##########-------------##########-------------
# Step 4: Process results into unique table
echo "Running python script tf_running.py to process results..."

python3 tfa_process.py --diamond ${diamond_output} --tfs_domains "${tfs_domains}" --ips ${interproscan_output} --output ${output_prefix}  --deeptf_res ${deeptfactor_output}

# Check if the Python script ran successfully
if [ $? -ne 0 ]; then
    echo "Error running Python analysis script"
    exit 1
fi

echo "tfa pipeline finished without errors."
echo "final tfa result saved with extension *tfa_results.csv"
