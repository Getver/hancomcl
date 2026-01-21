# import pandas as pd
# import subprocess
# import os

# ## blast 결과 중 ncbi_id 중복 없이 저장 후 tax_id mapping
# os.system("awk -F '|' '{print $4}' taxonomy/blast_results.tsv | sort -u > accessions.txt")
# a = pd.read_csv(subprocess.Popen("awk 'NR==FNR{keep[$1]; next} ($2 in keep)' accessions.txt /work/ncbi/nucl_gb.accession2taxid", stdout=subprocess.PIPE, shell=True).stdout, sep='\t', header=None, names=['name', 'name_v', 'tax_id', 'otherDBid'])

# ## blast 결과 가져오기
# blast = pd.read_csv('taxonomy/blast_results.tsv', sep='\t', header=None)
# blast.columns = ['ASV id', 'sseqid', 'Identity', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'staxids', 'sscinames']

# ## sseqid 저장 후 tax_id 붙임
# blast['name_v'] = blast['sseqid'].str.split('|').str[3]
# blast_anno = pd.merge(blast, a[['name_v', 'tax_id']], on='name_v', how='left')

# ## 풍부도 테이블 불러옴(ASV_id, 풍부도)
# abundance = pd.read_csv('abundance/abundance.tsv', sep='\t', header=1)

# ## fastq 파일에서 ASV_id와 sequence 가져옴
# records = []
# with open('filtered/seq/denoise_noMT_seq.fasta') as f:
#     name = None
#     seq_lines = []
#     for line in f:
#         line = line.strip()
#         if line.startswith(">"):
#             # 이전 시퀀스 저장
#             if name and seq_lines:
#                 records.append({"ASV id": name, "Sequence": "".join(seq_lines)})
#             name = line[1:]  # '>' 제거
#             seq_lines = []
#         else:
#             seq_lines.append(line)
#     # 마지막 시퀀스 저장
#     if name and seq_lines:
#         records.append({"ASV id": name, "Sequence": "".join(seq_lines)})

# seq = pd.DataFrame(records)

# ## 미생물 DB 불러옴
# taxonomy = pd.read_csv('/work/ncbi/TaxonomyDB.csv')[['tax_id', 'name_txt', 'rank', 'lineage_path', 'division_id']]
# division = pd.read_csv('/work/ncbi/DivisionDB.csv')

# ## 미생물에 division_id 붙임
# taxDB = pd.merge(taxonomy, division[['division_id', 'division_name']], on='division_id', how='left')

# ## blast 결과(ASV에 미생물 id 있는 것 + tax_id)에 sequence 붙임
# df = pd.merge(blast_anno, seq, on='ASV id', how='left')

# ## tax_id에 division_name 붙임
# df2 = pd.merge(df, taxDB, on='tax_id', how='left').drop('division_id', axis=1)

# ## 미생물 annotation 정확도 역순으로 정렬 후 겹치는 ASC_id 없앰
# df3 = df2.sort_values('Identity', ascending=False).drop_duplicates(subset='ASV id')

# ## 풍부도 붙임
# df4 = pd.merge(df3, abundance, left_on='ASV id', right_on='#OTU ID', how='left').drop('#OTU ID', axis=1)

# df4.to_csv('RESULT_ASV_Table.txt', sep='\t', index=False)

# os.system("rm accessions.txt")



####################################
import pandas as pd

# --------------------------------------------------------
# rank map
rank_map = {
    "domain": "k__",
    "superkingdom": "k__",
    "kingdom": "k__",
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
    "strain": "t__"
}

# --------------------------------------------------------
# lineage_path → dict
def lineage_to_dict(lineage_path):
    if lineage_path is None or pd.isna(lineage_path):
        return {}
    ids = [x.strip() for x in str(lineage_path).split(";") if x.strip().isdigit()]
    result = {}
    for taxid in ids:
        name = tax_lookup_name.get(str(taxid))
        rank = tax_lookup_rank.get(str(taxid))
        if not name or not rank:
            continue
        if rank.lower() == "no rank":
            continue
        result[rank.lower()] = name
    return result

# dict → QIIME2-style taxonomy
def lineage_dict_to_qiime2(lineage_dict):
    if not lineage_dict:
        return None
    order = ["domain","kingdom","phylum","class","order","family","genus","species","strain"]
    parts = []
    for r in order:
        if r in lineage_dict:
            prefix = rank_map.get(r, "")
            parts.append(f"{prefix}{lineage_dict[r]}")
    return "; ".join(parts) if parts else None

# dict → 개별 컬럼
def extract_tax_columns_from_dict(lineage_dict):
    return pd.Series({
        "Kingdom": lineage_dict.get("domain") or lineage_dict.get("kingdom"),
        "Phylum": lineage_dict.get("phylum"),
        "Class": lineage_dict.get("class"),
        "Order": lineage_dict.get("order"),
        "Family": lineage_dict.get("family"),
        "Genus": lineage_dict.get("genus"),
        "Species": lineage_dict.get("species"),
        "Strain": lineage_dict.get("strain")
    })

# --------------------------------------------------------
# 데이터 불러오기
blast = pd.read_csv('taxonomy/blast_results.tsv', sep='\t', header=None)
blast.columns = ['ASV id', 'sseqid', 'Identity', 'length', 'mismatch', 'gapopen', 
                 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'staxids', 'sscinames']

taxonomy = pd.read_csv('/work/ncbi/TaxonomyDB.csv')[['tax_id', 'name_txt', 'rank', 'lineage_path', 'division_id']]
taxonomy = taxonomy.astype({"tax_id": str})

abundance = pd.read_csv('abundance/abundance.tsv', sep='\t', header=1)

# sseqid → GI, NR ID
blast['GI ID'] = blast['sseqid'].str.split('|').str[1]
blast['NR ID'] = blast['sseqid'].str.split('|').str[3]

# --------------------------------------------------------
# tax_lookup
tax_lookup_name = dict(zip(taxonomy['tax_id'], taxonomy['name_txt']))
tax_lookup_rank = dict(zip(taxonomy['tax_id'], taxonomy['rank']))

# fastq에서 ASV id, sequence
records = []
with open('filtered/seq/denoise_noMT_seq.fasta') as f:
    name = None
    seq_lines = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if name and seq_lines:
                records.append({"ASV id": name, "Sequence": "".join(seq_lines)})
            name = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if name and seq_lines:
        records.append({"ASV id": name, "Sequence": "".join(seq_lines)})


seq = pd.DataFrame(records)

# --------------------------------------------------------
# merge
df1 = pd.merge(blast, seq, on='ASV id', how='inner')
df1['staxids'] = df1['staxids'].astype(str)

df2 = pd.merge(df1, taxonomy, left_on='staxids', right_on='tax_id', how='left')

# lineage dict, QIIME2 taxonomy, 개별 컬럼
lineage_dicts = df2["lineage_path"].apply(lineage_to_dict)
df2["Taxonomy"] = lineage_dicts.apply(lineage_dict_to_qiime2)
tax_cols = lineage_dicts.apply(extract_tax_columns_from_dict)

df_tax_final = pd.concat([df2, tax_cols], axis=1)


final_cols = ['ASV id','Taxonomy','Kingdom','Phylum','Class','Order','Family','Genus','Species','Strain',
              'Identity','Sequence','GI ID','NR ID','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore']

# abundance 붙이기
final_table = pd.merge(df_tax_final[final_cols], abundance, left_on='ASV id', right_on='#OTU ID', how='inner').drop('#OTU ID', axis=1)


final_table.to_csv('RESULT_ASV_TABLE.txt', index=False, sep='\t')


# # header=0으로 읽어서 첫 행을 컬럼 이름으로 사용
# count_table = pd.read_csv('filtered/txt/denoise_noMT_table.tsv', sep='\t', header=1, index_col=0, skiprows=1)

# # 모든 sample 열 합산해서 'Count' 컬럼 생성
# count_table['Count'] = count_table.sum(axis=1)

# print(count_table['Count'].sort_values(ascending=False).head())

qiime2_taxonomy = final_table[['ASV id', 'Taxonomy', 'Identity']].copy()
qiime2_taxonomy.columns = ['Feature ID', 'Taxon', 'Confidence']   # QIIME2 header 형식 맞추기

qiime2_taxonomy.to_csv('taxonomy/ncbi_qiime2.tsv', sep='\t', index=False)
# --------------------------------------------------------

result_asv_table = (final_table.sort_values(['bitscore', 'Identity'], ascending=False).drop_duplicates(subset='ASV id', keep='first'))

# 저장
result_asv_table.to_csv('RESULT_ASV_TABLE.txt', sep='\t', index=False)