#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import pandas as pd

# --------------------------------------------------
# Utility
# --------------------------------------------------

def run(cmd, step_name):
    print(f"\n[{step_name}] Running:\n{' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        sys.exit(f"ERROR during step: {step_name}")

def basename_no_ext(path):
    return os.path.basename(path).rsplit(".", 1)[0]

# --------------------------------------------------
# TF family processing logic (from tfa_process.py)
# --------------------------------------------------

def process_results(diamond, ips, deeptf_res, tfs_domains, output_prefix):
    print("\n[RESULT MERGING] Processing TF results")

    known_DomTFS = pd.read_csv(
        tfs_domains, sep=",", header=None, names=["family", "domain"]
    )
    known = known_DomTFS.explode("domain")

    family_dict = (
        known_DomTFS.groupby("family")["domain"]
        .apply(lambda x: x.tolist() if len(x) > 1 else x.iloc[0])
        .to_dict()
    )

    ref_sets_list = [
        set(v) if isinstance(v, list) else {v}
        for v in family_dict.values()
    ]
    family_keys_list = list(family_dict.keys())

    # ---- InterProScan ----
    ips_df = pd.read_csv(
        ips, sep="\t", header=None,
        names=[
            "qseqid", "md5", "len", "db", "domain", "description",
            "start", "end", "evalue", "type", "date",
            "ipr", "ipr_description", "x1", "x2"
        ]
    )

    query_dict = (
        ips_df.groupby("qseqid")["domain"]
        .apply(lambda x: x.tolist() if len(x) > 1 else x.iloc[0])
        .to_dict()
    )

    def filter_tffam(value):
        value_set = set(value) if isinstance(value, list) else {value}
        matches = []

        for idx, ref_set in enumerate(ref_sets_list):
            if ref_set.issubset(value_set):
                fam = family_keys_list[idx]
                if fam == "Homeobox":
                    if "PF03529" in value_set:
                        matches.append("TF_Otx")
                    elif "PF00157" in value_set:
                        matches.append("Pou")
                    elif "PF02376" in value_set:
                        matches.append("CUT")
                    else:
                        matches.append("Homeobox")
                else:
                    matches.append(fam)

        return matches[0] if matches else None

    family_map = {k: filter_tffam(v) for k, v in query_dict.items()}

    # ---- DeepTFactor ----
    deeptf_df = pd.read_csv(deeptf_res, sep="\t")
    deeptf_df = deeptf_df.rename(columns={"sequence_ID": "qseqid"})
    deeptf_df = deeptf_df[deeptf_df["prediction"] == True]

    # ---- Diamond ----
    diamond_df = pd.read_csv(
        diamond, sep="\t", header=None,
        names=[
            "qseqid", "qlen", "sseqid", "slen", "salltitles",
            "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send",
            "qcovhsp", "evalue", "bitscore"
        ]
    )

    merged = deeptf_df.merge(
        diamond_df[["qseqid", "salltitles", "pident", "qcovhsp", "evalue"]],
        on="qseqid",
        how="left"
    )

    merged = merged.merge(
        ips_df[["qseqid", "domain", "description", "start", "end"]],
        on="qseqid",
        how="left"
    )

    merged["domain"] = (
        merged.groupby("qseqid")["domain"]
        .transform(lambda x: "; ".join(sorted(set(x.dropna()))))
    )
    merged["description"] = (
        merged.groupby("qseqid")["description"]
        .transform(lambda x: "; ".join(sorted(set(x.dropna()))))
    )

    merged["family_keys"] = merged["qseqid"].map(family_map)

    merged["salltitles"] = (
        merged["salltitles"]
        .str.split(" OS=").str[0]
        .str.split(" ").str[1:].str.join(" ")
    )

    merged = merged.drop_duplicates("qseqid")

    # Filter by DeepTFactor score
    merged = merged[merged["score"] >= 0.5]
    print(f"After DeepTFactor score filter: {merged.shape[0]}")

    mask_transposase = merged["description"].str.contains("transposase", case=False, na=False)

    merged = merged[~mask_transposase]
    print(f"After transposase removal: {merged.shape[0]}")
    
    # Remove sequences without TF family classification (no known TF DBDs)
    merged = merged[merged["family_keys"].notna()]
    print(f"After family classification filter: {merged.shape[0]}")

    # -----------------------------------------
    # Output formatting
    # -----------------------------------------

    final_cols = [
        "qseqid", "prediction", "score", "salltitles",
        "pident", "qcovhsp", "evalue",
        "domain", "description", "start", "end", "family_keys"
    ]

    out_file = f"{output_prefix}.tfa_results.csv"
    merged[final_cols].to_csv(out_file, index=False)

    
    print(f"[DONE] Results written to {out_file}")

# --------------------------------------------------
# Main pipeline
# --------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Unified TF annotation pipeline"
    )
    parser.add_argument("--input_fasta", required=True)
    parser.add_argument("--db_path", required=True)
    parser.add_argument("--deeptfactor_folder", required=True)
    parser.add_argument("--tfs_domains", required=True)
    parser.add_argument("--threads", type=int, default=50)
    parser.add_argument("--output_prefix", required=True,help="Prefix for all output files")

    args = parser.parse_args()

    prefix = args.output_prefix

    diamond_out = f"{prefix}.atfdb.1e3.txt"
    tf_headers = f"{prefix}.tfheaders.txt"
    filtered_fasta = f"{prefix}.filtered_sequences.fasta"
    ips_out = f"{prefix}.ips_results.tsv"
    deeptf_outdir = f"{prefix}.deeptf.res"
    deeptf_pred = f"{deeptf_outdir}/prediction_result.txt"
    
    # ---- Diamond ----
    run([
        "diamond", "blastp",
        "--query", args.input_fasta,
        "--db", args.db_path,
        "--max-target-seqs", "1",
        "--evalue", "1e-3",
        "--outfmt", "6",
        "qseqid", "qlen", "sseqid", "slen", "salltitles",
        "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send",
        "qcovhsp", "evalue", "bitscore",
        "--out", diamond_out,
        "--threads", str(args.threads),
        "--sensitive"
    ], "DIAMOND")

    # ---- Filter FASTA ----
    run(["cut", "-f1", diamond_out], "EXTRACT HEADERS")
    run(["seqtk", "subseq", args.input_fasta, tf_headers], "SEQTK")

    os.rename("stdout", filtered_fasta) if os.path.exists("stdout") else None

    # ---- InterProScan ----
    run([
        "interproscan",
        "-appl", "Pfam",
        "-cpu", str(args.threads),
        "-f", "TSV",
        "-dra",
        "-iprlookup",
        "-i", filtered_fasta,
        "-o", ips_out
    ], "INTERPROSCAN")

    # ---- DeepTFactor ----
    cwd = os.getcwd()
    os.chdir(args.deeptfactor_folder)

    run([
        sys.executable, "tf_running.py",
        "-cpu", str(args.threads),
        "-i", os.path.join(cwd, args.input_fasta),
        "-o", os.path.join(cwd, deeptf_outdir)
    ], "DEEPTFACTOR")

    os.chdir(cwd)

    # ---- Merge results ----
    process_results(
        diamond=diamond_out,
        ips=ips_out,
        deeptf_res=deeptf_pred,
        tfs_domains=args.tfs_domains,
        output_prefix=prefix
    )

    print("\nPipeline finished successfully.")

if __name__ == "__main__":
    main()
