import pandas as pd

# 데이터 불러오기
abun = pd.read_csv("abundance/abundance.tsv", sep="\t", index_col=0, header=1)
meta = pd.read_csv("/work/meta_data.tsv", sep="\t", index_col=0)[['day']]

abun_T = abun.T
# SampleID 기준으로 합치기
df = meta.join(abun_T).T # 풍부도 테이블을 transpose해서 SampleID가 index가 되도록

df = df.reset_index()  # index -> 컬럼 '#SampleID'로 이동
df.to_csv("LEfSe/lefse_input_combined.tsv", sep="\t", index=False)
# 합친 테이블 저장
# df.to_csv("LEfSe/lefse_input_combined.tsv", sep="\t", index=True)
