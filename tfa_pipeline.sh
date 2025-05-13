#!/bin/bash

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
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check if required arguments are provided
if [[ -z "$input_fasta" || -z "$db_path" || -z "$deeptfactor_folder" || -z "$tfs_domains" ]]; then
    echo "Usage: $0 --input_fasta <protein_fasta_file> --db_path <path_to_atfdb_db> --deeptfactor_folder <deeptfactor_folder> --tfs_domains <dataframe with tfs domains>"
    exit 1
fi

# Print assigned variables
echo "Input FASTA file: $input_fasta"
echo "Database path: $db_path"
echo "DeepTFactor folder: $deeptfactor_folder"
echo "TFS domains csv: $tfs_domains"

# Define output files
diamond_output="$(basename ${input_fasta%.faa}).atfdb.1e3.txt"
filtered_fasta="$(basename ${input_fasta%.faa}).filtered_sequences.fasta"
interproscan_output="$(basename ${input_fasta%.faa}).ips_results.tsv"

# Step 1: Run diamond blastp
echo "Running diamond blastp..."
diamond blastp --query "$input_fasta" --max-target-seqs 1 --evalue 1e-3 --db "$db_path" --outfmt 6 qseqid qlen sseqid slen salltitles pident length mismatch gapopen qstart qend sstart send qcovhsp evalue bitscore --out "$diamond_output" --threads 50 --sensitive

# Step 1.1: Check if diamond ran successfully
if [ $? -ne 0 ]; then
    echo "Error running diamond blastp"
    exit 1
fi

# Step 1.2: Extract putative tf headers
cut -f 1 "$diamond_output" > "$(basename ${input_fasta%.faa}).tfheaders.txt"

# Step 1.3: Filter original fasta file using seqtk and the list of tf headers
seqtk subseq "$input_fasta" "$(basename ${input_fasta%.faa}).tfheaders.txt" > "$filtered_fasta"

# Step 1.3: Check if the filtering was successful
if [ $? -ne 0 ]; then
    echo "Error filtering fasta file"
    exit 1
fi

# Step 2: Run interproscan on the filtered fasta file
echo "Running interproscan..."
interproscan -appl Pfam -cpu 50 -f TSV -dra -iprlookup -i "$filtered_fasta" -o "$interproscan_output"

# Step 2.1: Check if interproscan ran successfully
if [ $? -ne 0 ]; then
    echo "Error running interproscan"
    exit 1
fi
original_dir=$(pwd)

# Step 3: Run DeepTFactor
echo "Running DeepTFactor..."
original_dir=$(pwd)
og_fasta="$original_dir/$input_fasta"

# Step 3.1: Move to deeptfactor folder %todo: fix deeptfactor installation so it doesnt have to do this
cd "$deeptfactor_folder" || exit

python tf_running.py -cpu 50 -i "$og_fasta" -o "$original_dir/$(basename ${input_fasta%.faa}).deeptf.res"

# Step 3.2: Check if deeptfactor ran successfully
if [ $? -ne 0 ]; then
    echo "Error running DeepTFactor"
    exit 1
fi

cd "$original_dir" || exit

# Step 4: Process results into unique table
echo "Running python script tf_running.py to merge results..."

python tfa_process.py  --diamond "$(basename ${input_fasta%.faa}).atfdb.1e3.txt" --tfs_domains "${tfs_domains}" --ips "$(basename ${input_fasta%.faa}).ips_results.tsv" --output $(basename ${input_fasta%.faa})  --deeptf_res "$original_dir/$(basename ${input_fasta%.faa}).deeptf.res/prediction_result.txt"

# Check if the Python script ran successfully
if [ $? -ne 0 ]; then
    echo "Error running Python analysis script"
    exit 1
fi

echo "tfa pipeline finished without errors."
echo "final tfa result saved with extension *tfa_results.csv"
