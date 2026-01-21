import pandas as pd


silva_result = pd.read_csv('taxonomy/taxonomy_silva.tsv', sep='\t')    # [486 rows x 3 columns]
ncbi_result = pd.read_csv('taxonomy/taxonomy_blast.tsv', sep='\t') # [377 rows x 3 columns]

ncbi_result.Confidence = ncbi_result.Confidence / 100

silva = silva_result.loc[silva_result.Confidence >= 0.97, :]    # [329 rows x 3 columns]
ncbi = ncbi_result.loc[ncbi_result.Confidence >= 0.97, :] # [272 rows x 3 columns]


def extract_tax_levels(df, col='Taxon'):
    df['kingdom'] = df[col].str.extract(r'k__([^;]+)')
    df['phylum']  = df[col].str.extract(r'p__([^;]+)')
    df['class']   = df[col].str.extract(r'c__([^;]+)')
    df['order']   = df[col].str.extract(r'o__([^;]+)')
    df['family']  = df[col].str.extract(r'f__([^;]+)')
    df['genus']   = df[col].str.extract(r'g__([^;]+)')
    df['species'] = df[col].str.extract(r's__([^;]+)')
    return df


df_silva = extract_tax_levels(silva)    # [329 rows x 10 columns]
df_ncbi = extract_tax_levels(ncbi)  # [272 rows x 10 columns]

merged = pd.merge(df_silva, df_ncbi, on='Feature ID', how='inner', suffixes=('_silva', '_ncbi'))    # [222 rows x 19 columns]

merged['same_genus'] = merged['genus_silva'] == merged['genus_ncbi']
merged['same_family'] = merged['family_silva'] == merged['family_ncbi']
merged['same_order'] = merged['order_silva'] == merged['order_ncbi']


for level in ['genus', 'family', 'order']:
    merged[f'{level}_silva'] = merged[f'{level}_silva'].fillna('Unknown')
    merged[f'{level}_ncbi'] = merged[f'{level}_ncbi'].fillna('Unknown')


matched = merged[(merged['same_genus'] == True) & (merged['same_family'] == True) & (merged['same_order'] == True)] # [113 rows x 22 columns]
mismatch = merged[(merged['same_genus'] == False) | (merged['same_family'] == False) | (merged['same_order'] == False)] # [109 rows x 22 columns]


final = matched[['Feature ID', 'Taxon_ncbi', 'Confidence_ncbi']]    # [113 rows x 3 columns]
final.columns = ['Feature ID', 'Taxon', 'Confidence']

# final.to_csv('taxonomy/taxonomy.tsv', sep='\t', index=False)

ncbi_result.to_csv('taxonomy/taxonomy.tsv', sep='\t', index=False)
