import pandas as pd
import os, glob
import argparse

parser = argparse.ArgumentParser(description='Combine Diamond and InterProScan results')

# Add arguments
parser.add_argument('--diamond', type=str, help='diamond result')
parser.add_argument('--ips', type=str, help='interproscan result')
parser.add_argument('--output', type=str, help='name to be used in output')

# Parse arguments
args = parser.parse_args()

#####

#tfs of interest
#tfs_of_interest=pd.read_csv("tfs_of_interest_fam_path_2.csv")
#tfs_to_search=list(tfs_of_interest['Symbol'].unique())


known_DomTFS=pd.read_csv("/project/yutaka/spiderNetwork/data/AnimalTfdb_concat/TFsdomains.csv", sep=',', header=None, names=['family', 'domain'])
known=known_DomTFS.explode('domain')

family_dict = known_DomTFS.groupby('family')['domain'].apply(lambda x: x.tolist() if len(x) > 1 else x.iloc[0]).to_dict()

def extract_accession_id(file_path):
    return file_path.split('/')[-1].split('.')[0]


#folder_path = '../new_atfdb/tfs/*tsv'
#tf_files_ips = glob.glob(os.path.join(folder_path))
#tf_files_ips=sorted(tf_files_ips)
#lista_query_domains=[]
all_ips_TFres=[]

tf_files_ips = args.ips
print(tf_files_ips)

query_dicts=[]
df=pd.read_csv(tf_files_ips, sep='\t', header=None,names=['query', 'md5', 'len', 'db','domain', 'description', 'start', 'end', 'evalue', 'type', 'date', 'ipr', 'ipr_description','nosei', 'nosei2'])
dictio=df.groupby('query')['domain'].apply(lambda x: x.tolist() if len(x) > 1 else x.iloc[0]).to_dict()
query_dicts.append(dictio)
all_ips_TFres.append(df)


ref_sets_list = [set(val) if isinstance(val, list) else {val} for val in family_dict.values()]

# Define the filter_test function to return the family keys
def filter_test(value):
    value_set = set(value) if isinstance(value, list) else {value}
    matching_keys = []
    for idx, ref_set in enumerate(ref_sets_list):
        if value_set == ref_set:  # Check for exact match
            matching_keys.append(list(family_dict.keys())[idx])
    return len(matching_keys) > 0, matching_keys

# Initialize the filtered_query list
filtered_query = []

# Iterate through each search dictionary in query_dicts
for search_dict in query_dicts:
    new_dict = {}
    # Iterate through key-value pairs in each search dictionary
    for key, value in search_dict.items():
        matched, family_keys = filter_test(value)
        if matched:
            # Check if the key already exists in new_dict
            if key not in new_dict:
                new_dict[key] = {'query_value': value, 'family_keys': family_keys}
            # If the key exists, append the family keys to the existing list
            else:
                new_dict[key]['family_keys'].extend(family_keys)
    # Append the new_dict to the filtered_query list
    filtered_query.append(new_dict)

for item in filtered_query:
    for key, value in item.items():
        value['family_keys'] = value['family_keys'][0]
#print(filtered_query[-1])

#diam_path= '../new_atfdb/tfs/*1e3.txt'
#diam_path = glob.glob(os.path.join(diam_path))
#diam_path = sorted(diam_path)

diam_path=args.diamond

diam_dfs=[]
diam_df=pd.read_csv(diam_path, sep='\t', header=None, names=['qseqid', 'qlen', 'sseqid', 'slen', 'ensembl_id', 'ensembl_protein', 'protein', 'tffam', 'salltitles', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'qcovhsp', 'evalue'] )
accession_id = extract_accession_id(diam_path)
diam_df['source'] = accession_id
diam_dfs.append(diam_df)
print(diam_dfs[-1])


filtered_dfs=[]


for dictio in filtered_query:
    records=[]
    for query, details in dictio.items():
        query_values = details['query_value']
        if isinstance(query_values, list):
            for value in query_values:
                records.append((query, value, details['family_keys']))
        else:
            records.append((query, query_values, details['family_keys']))

        df = pd.DataFrame(records, columns=['qseqid', 'architecture', 'tffam'])
    filtered_dfs.append(df)

final_dfs=[]
for i, df in enumerate(filtered_dfs):
    df=df.merge(diam_dfs[i][['qseqid', 'protein', 'source']], on='qseqid')
    df['tffam'] = df['tffam'].fillna(df['protein'])
    final_dfs.append(df)

output_name = args.output
final_dfs[0].to_csv(output_name + ".diamips.tfs.csv", header=True, index=False)