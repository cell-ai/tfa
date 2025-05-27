import pandas as pd
import os, glob
import argparse

parser = argparse.ArgumentParser(description="Combine Diamond, InterProScan and DeepTFactor results")

# Add arguments
parser.add_argument("--diamond", type=str, help="diamond result")
parser.add_argument("--ips", type=str, help="interproscan result")
parser.add_argument("--deeptf_res", type=str, help="deeptfactor result")
parser.add_argument("--output", type=str, help="name to be used in output")
parser.add_argument("--tfs_domains", type = str, help = "dataframe with known TF domains (PFAM)")

# Parse arguments
args = parser.parse_args()

#Reference csv of tffam and pfam correspondence
known_DomTFS=pd.read_csv(args.tfs_domains, sep=",", header=None, names=["family", "domain"])
known=known_DomTFS.explode("domain")

family_dict = known_DomTFS.groupby("family")["domain"].apply(lambda x: x.tolist() if len(x) > 1 else x.iloc[0]).to_dict()

def extract_accession_id(file_path):
    return file_path.split("/")[-1].split(".")[0]


all_ips_TFres=[]

tf_files_ips = args.ips
#print(tf_files_ips)

query_dicts=[]
df=pd.read_csv(tf_files_ips, sep="\t", header=None,names=["qseqid", "md5", "len", "db","domain", "description", "start", "end", "evalue", "type", "date", "ipr", "ipr_description","nosei", "nosei2"])
dictio=df.groupby("qseqid")["domain"].apply(lambda x: x.tolist() if len(x) > 1 else x.iloc[0]).to_dict()
query_dicts.append(dictio)
all_ips_TFres.append(df)

ref_sets_list = [set(val) if isinstance(val, list) else {val} for val in family_dict.values()]

def filter_tffam_based_on_pfam(value):
    value_set = set(value) if isinstance(value, list) else {value}
    matching_keys = []
    family_keys_list = list(family_dict.keys())  # Convert to list for indexing

    for idx, ref_set in enumerate(ref_sets_list):
        if ref_set.issubset(value_set):
            if family_keys_list[idx] == "Homeobox":
                if "PF03529" in value_set:
                    matching_keys.append("TF_Otx")
                elif "PF00157" in value_set:
                    matching_keys.append("Pou")
                elif "PF02376" in value_set:
                    matching_keys.append("CUT")
                elif not any(pfam in value_set for pfam in ["PF03529", "PF00157", "PF02376"]):
                    matching_keys.append("Homeobox")
            else:
                matching_keys.append(family_keys_list[idx])
    
    return len(matching_keys) > 0, matching_keys

filtered_query = []

for search_dict in query_dicts:
    new_dict = {}
    for key, value in search_dict.items():
        matched, family_keys = filter_tffam_based_on_pfam(value)
        if matched:
            if key not in new_dict:
                new_dict[key] = {"query_value": value, "family_keys": family_keys}
            else:
                new_dict[key]["family_keys"].extend(family_keys)
    filtered_query.append(new_dict)

# Keep only the first family_key for each
for item in filtered_query:
    for key, value in item.items():
        value["family_keys"] = value["family_keys"][0] if value["family_keys"] else None

# Read deeptfactor results
deeptf_files = args.deeptf_res
deeptf_df = pd.read_csv(deeptf_files, sep="\t")
deeptf_df = deeptf_df.rename(columns = {"sequence_ID": "qseqid"})

# Read diamond-blastp results
diamond_files = args.diamond
diamond_df = pd.read_csv(diamond_files, sep="\t", header=None, names=["qseqid", "qlen", "sseqid", "slen", "salltitles", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qcovhsp", "evalue", "bitscore"])

merge_deeptf_diamond = deeptf_df.merge(diamond_df[["qseqid", "salltitles", "pident", "qcovhsp", "evalue"]], on="qseqid", how="outer")

merge_deeptf_diamond = merge_deeptf_diamond[merge_deeptf_diamond["prediction"] == True]

# Merge with InterProScan results
merge_deeptf_diam_ips = merge_deeptf_diamond.merge(all_ips_TFres[0][["qseqid", "domain", "description", "start", "end"]], on="qseqid", how="outer")

# Collapse the results of domain to the same row, based on qseqid value
merge_deeptf_diam_ips["domain"] = merge_deeptf_diam_ips.groupby("qseqid")["domain"].transform(lambda x: "; ".join(x.dropna().unique()))
merge_deeptf_diam_ips["description"] = merge_deeptf_diam_ips.groupby("qseqid")["description"].transform(lambda x: "; ".join(x.dropna().unique()))

# Add family keys to the merged DataFrame
merge_deeptf_diam_ips["family_keys"] = merge_deeptf_diam_ips["qseqid"].map(lambda x: filtered_query[0].get(x, {}).get("family_keys", None))

merged_deeptf_diam_ips = merge_deeptf_diam_ips[["qseqid", "prediction", "score", "salltitles", "pident", "qcovhsp", "evalue", "domain", "description", "start", "end", "family_keys"]]
merged_deeptf_diam_ips["salltitles"] = merged_deeptf_diam_ips["salltitles"].str.split(" OS=").str[0].str.split(" ").str[1:].str.join(" ")

merged_deeptf_diam_ips.to_csv(args.output + ".tfa_results.csv", index=False)

print(f"Output saved to {args.output}" + ".tfa_results.csv")
print("tfa completed successfully.")