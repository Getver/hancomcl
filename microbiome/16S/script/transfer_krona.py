import pandas as pd

# TSV 불러오기
df = pd.read_csv("taxonomy/taxonomy.tsv", sep='\t')

# Taxon 컬럼 ; 기준으로 분리
df_tax = df['Taxon'].str.split(';', expand=True)

# Count 컬럼 생성 (예: 1 read per feature)
df_out = pd.concat([pd.Series([1]*len(df), name='Count'), df_tax], axis=1)

# Krona용 TSV 저장
df_out.to_csv("taxonomy/krona_input.tsv", sep='\t', index=False, header=False)