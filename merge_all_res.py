import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Merging diamond, interproscan and deeptfactor results. Produces a table with the prediction of deeptfactor for TFs predicted by homology + domain conservation and one table of putative TFs based only on deeptfactor')

parser.add_argument('--diamond_ips_res', type=str, help='diamond and interproscan results table')
parser.add_argument('--deeptf_res', type=str, help='interproscan result')
parser.add_argument('--output', type=str, help='name to be used in output')

args = parser.parse_args()

diam_path = args.diamond_ips_res
deeptf_path = args.deeptf_res
output_name = args.output

diam_ips_df = pd.read_csv(diam_path)

deeptf_df = pd.read_csv(deeptf_path,sep='\t')
deeptf_df = deeptf_df.rename(columns={'sequence_ID': 'qseqid'})
deeptf_df = deeptf_df[deeptf_df['prediction'] == True]

def add_existence_column(df1, df2, column_name, new_column_name='deeptf_res'):
    df1[new_column_name] = df1[column_name].isin(df2[column_name])
    return df1

homology_df = add_existence_column(diam_ips_df, deeptf_df, 'qseqid')
homology_df.to_csv(output_name + ".homologytfs.csv", header=True, index=False)

def get_exclusive_values(df1, df2, column_name):
    exclusive_df = df2[~df2[column_name].isin(df1[column_name])]
    return exclusive_df

exclusive_deeptf_df = get_exclusive_values(diam_ips_df, deeptf_df, 'qseqid')
exclusive_deeptf_df.to_csv(output_name + ".no_homologytfs.csv", header=True, index=False)
