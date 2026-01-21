import pandas as pd

# 파일 경로
tax_path = "taxonomy/taxonomy.tsv"
abund_path = "abundance/abundance.tsv"

# taxonomy 파일 불러오기
tax = pd.read_csv(tax_path, sep="\t", index_col=0)
tax = tax.rename(columns={"Taxon": "clade_name"})

# abundance 파일 불러오기
abund = pd.read_csv(abund_path, sep="\t", skiprows=1, index_col=0)

# 평균 풍부도 계산 (샘플 여러 개일 경우)
abund["relative_abundance"] = abund.mean(axis=1) * 100  # 백분율

# 병합
merged = pd.merge(tax[["clade_name"]], abund[["relative_abundance"]], left_index=True, right_index=True)

# NCBI_tax_id는 없는 경우가 많으므로 placeholder 사용
merged["NCBI_tax_id"] = "NA"

# 열 순서 조정
merged = merged[["clade_name", "NCBI_tax_id", "relative_abundance"]]

# 헤더 추가
header_lines = [
    "#SampleID\tMetaphlan_Analysis",
    "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species"
]

# 출력
out_path = "metaphlan_style_profile.txt"
with open(out_path, "w") as f:
    f.write("\n".join(header_lines) + "\n")
    merged.to_csv(f, sep="\t", index=False)

