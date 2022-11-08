import pandas as pd

file = "Structure_Analysis_Descriptors_TSSIG_PIP_Strict_Complete_Cleaned.csv"

df = pd.read_csv(file, index_col='Ligand')
#print(df)
# Create bin range based on the range of activation energies
bins = list(range(-25, 125, 1))
#print(bins)

# Bin data based on activation energy
df['binned'] = pd.cut(df['Activation Energy (kcal/mol)'], bins)
#print(df)

# df_sample = df.sample(frac=0.5)
# print(df_sample)

list_of_binned_dfs = list(zip(df.groupby('binned')))

sampled_df = pd.DataFrame(columns=df.columns)

for data_bin in list_of_binned_dfs:
    if len(data_bin[0][1].index) < 250:
        sampled_df = sampled_df.append(data_bin[0][1])
    else:
        df_sample = data_bin[0][1].sample(250)
        sampled_df = sampled_df.append(df_sample, ignore_index=False)

print(sampled_df)

outfile = "Structure_Analysis_Descriptors_TSSIG_PIP_Strict_Complete_Cleaned_Trimmed.csv"

sampled_df.to_csv(outfile)
