import pandas as pd
# import sys
# from io import StringIO
import os, subprocess
import modin.pandas as mpd
# from functools import reduce



## 여기는 테이블 만든다(minfi 활용..)

dmp = pd.read_csv('/disk0/sm/methyl/03_pollution/total/DMP_Minfi_dmpFinder.txt', sep='\t', low_memory=False)
dmr = pd.read_csv('/disk0/sm/methyl/03_pollution/total/DMR_Minfi_bumphunter.txt', sep='\t', low_memory=False)
gtf = pd.read_csv(subprocess.Popen('cat /disk0/references/hg38/UCSC/RNA/all_chromosomes/gencode/gencode.v44.basic.annotation.gtf | grep -v "^#"', shell = True, stdout = subprocess.PIPE).stdout, sep='\t', names=['chr', 'data', 'type', 'start', 'end', '1', 'strand', '2', 'info'], low_memory=False)
probe = pd.read_csv('/disk0/sm/methyl/probe_v2_distance.txt', sep='\t', low_memory=False)


dmrbed = dmr[['chr', 'start', 'end', 'value']]
gtfbed = gtf[['chr', 'start', 'end', 'info', 'type']]

gtfbed['tmp'] = gtfbed['info'].str.split(';')

for key, value in gtfbed.iterrows():
    for v in range(len(value['tmp'])):
        if 'gene_name' in value['tmp'][v]:
            gtfbed.at[key, 'gene_name'] = value['tmp'][v].split(' ')[2]
        elif 'gene_type' in value['tmp'][v]:
            gtfbed.at[key, 'gene_type'] = value['tmp'][v].split(' ')[2]

gtfbed = gtfbed[['chr', 'start', 'end', 'type', 'gene_name', 'gene_type']]

## 임시파일
dmrbed.to_csv('/disk0/sm/methyl/bed_dmr.bed', sep='\t', header=None, index=False)
gtfbed.to_csv('/disk0/sm/methyl/bed_base.bed', sep='\t', header=None, index=False)


## 나온 dmr 정보에 gtf의 gene 정보를 붙여줌
df = pd.read_csv(subprocess.Popen('bedtools intersect -a /disk0/sm/methyl/bed_dmr.bed -b /disk0/sm/methyl/bed_base.bed -wa -wb', shell = True, stdout = subprocess.PIPE).stdout, sep='\t', names=['chr', 'start', 'end', 'value', 'chrr', 'genestart', 'geneend', 'gene', 'gene_name', 'gene_type'])
df2 = df.drop(['chrr', 'genestart', 'geneend'], axis = 1)
df3 = df2[df2.gene == 'gene']
df4 = pd.merge(dmrbed, df3, how='left', on=['chr', 'start', 'end', 'value']).drop_duplicates(subset=['chr', 'start', 'end', 'gene']).replace('"', '', regex=True)

## gwas annotation 한거를 가져옴
# probe = pd.read_csv('/disk0/sm/methyl/probe_v2.txt', sep='\t', low_memory=False)
probe['end'] = probe['pos']+1
probebed = probe[['chr', 'pos', 'end', 'strand', 'Name', 'Probe_rs', 'CpG_rs', 'Islands_Name', 'Relation_to_Island', 'UCSC_RefGene_Name', 'Methyl450_Loci', 'Methyl27_Loci', 'EPICv1_Loci', 'gwas_rs', 'distance']]

## 임시파일
df4.to_csv('/disk0/sm/methyl/dmrs_anno.txt', sep='\t', index=False, header=None)


## gwas annotation 한거를 가져옴
probe = pd.read_csv('/disk0/sm/methyl/probe_v2_distance.txt', sep='\t', low_memory=False)
probe['end'] = probe['pos']+1
probebed = probe[['chr', 'pos', 'end', 'strand', 'Name', 'Probe_rs', 'CpG_rs', 'Islands_Name', 'Relation_to_Island', 'UCSC_RefGene_Name', 'Methyl450_Loci', 'Methyl27_Loci', 'EPICv1_Loci', 'gwas_rs', 'distance']]

# 임시파일
probebed.to_csv('/disk0/sm/methyl/bed_probe.txt', sep='\t', index=False, header=None)


## bedtools로 DMR 안에 있는 cpg와 그 annotation을 구함
result = pd.read_csv(subprocess.Popen('bedtools intersect -a /disk0/sm/methyl/dmrs_anno.txt -b /disk0/sm/methyl/bed_probe.txt -wa -wb', shell = True, stdout = subprocess.PIPE).stdout, sep='\t', 
names=['DMRchr', 'DMRstart', 'DMRend', 'DMRvalue', 'type', 'Gencode', 'gene_type', 'chr', 'pos', 'end', 'strand', 'Name', 'Probe_rs', 'CpG_rs', 'Islands_Name', 'Relation_to_Island', 'UCSC_RefGene_Name', 'Methyl450_Loci', 'Methyl27_Loci', 'EPICv1_Loci', 'gwas_rs', 'distance'])

r2 = result.drop(['chr', 'pos', 'end'], axis = 1)


# ## 거기에 dmp분석 결과로 나온 F값과 qvalue를 가져옴(dmp 열이름 확인하고 가져와야함)

dmp['num'] = range(len(dmp))
dmp['DMP'] = 'DMP_' +  dmp['num'].astype(str)
dmp = dmp.drop('num', axis = 1)

# dmp = dmp[dmp.pval < 0.00000001]

dmp.to_csv('/disk0/sm/methyl/03_pollution/total/data/DMP_numbering.txt', sep='\t')

r3 = pd.merge(r2, dmp, how='left', left_on='Name', right_index=True)

## 필요없어보이는 열 뺌
r3 = r3.drop(['Probe_rs', 'CpG_rs', 'Islands_Name', 'UCSC_RefGene_Name', 'Methyl450_Loci', 'Methyl27_Loci', 'EPICv1_Loci', 'strand', 'gwas_rs', 'distance'], axis = 1).fillna('-')


# 결과파일
r3.to_csv('/disk0/sm/methyl/03_pollution/total/data/DMRS_table.txt', sep='\t', index=False)
r4 = r3.drop_duplicates(subset=['DMRstart', 'DMRend']) # 이거의 gene을 세는 게 더 맞을 것 같음

geneC = pd.DataFrame(r3['Gencode'].value_counts())
# r3 = r3[abs(r3.f) > 1465.95]

## 결과파일
geneC.to_csv('/disk0/sm/methyl/03_pollution/total/data/gene_count.txt', sep='\t')

os.system('rm /disk0/sm/methyl/bed_data.bed /disk0/sm/methyl/bed_base.bed /disk0/sm/methyl/dmrs_anno.txt /disk0/sm/methyl/bed_probe.txt /disk0/sm/methyl/bed_dmr.bed')









import pandas as pd
# import sys
# from io import StringIO
import os, subprocess
import modin.pandas as mpd
# from functools import reduce



## 여기는 테이블 만든다(minfi 활용..)

dmp = pd.read_csv('/disk0/sm/methyl/03_pollution/total/DMP_Minfi_dmpFinder.txt', sep='\t', low_memory=False)
gtf = pd.read_csv(subprocess.Popen('cat /disk0/references/hg38/UCSC/RNA/all_chromosomes/gencode/gencode.v44.basic.annotation.gtf | grep -v "^#"', shell = True, stdout = subprocess.PIPE).stdout, sep='\t', names=['chr', 'data', 'type', 'start', 'end', '1', 'strand', '2', 'info'], low_memory=False)
probe = pd.read_csv('/disk0/sm/methyl/probe_v2_distance.txt', sep='\t', low_memory=False)

gtfbed = gtf[['chr', 'start', 'end', 'info', 'type']]

gtfbed['tmp'] = gtfbed['info'].str.split(';')

for key, value in gtfbed.iterrows():
    for v in range(len(value['tmp'])):
        if 'gene_name' in value['tmp'][v]:
            gtfbed.at[key, 'gene_name'] = value['tmp'][v].split(' ')[2]
        elif 'gene_type' in value['tmp'][v]:
            gtfbed.at[key, 'gene_type'] = value['tmp'][v].split(' ')[2]

gtfbed = gtfbed[['chr', 'start', 'end', 'type', 'gene_name', 'gene_type']]

gtfbed.to_csv('/disk0/sm/gtfbed.txt', sep='\t', header=False, index=False)

epicv2 = pd.read_csv('/disk0/sm/methyl/mchip/EPICv2.csv', sep=',', low_memory=False, header=7)[['CHR', 'MAPINFO', 'IlmnID', 'Name']]
epicv2['end'] = epicv2['MAPINFO'].astype('Int64') + 1
epicv2['start'] = epicv2['MAPINFO'].astype('Int64')
epicv2_bed = epicv2[['CHR', 'start', 'end', 'IlmnID', 'Name']]


epicv2_bed.to_csv('/disk0/sm/epicv2_bed.txt', sep='\t', header=False, index=False)


## 나온 dmr 정보에 gtf의 gene 정보를 붙여줌
df = pd.read_csv(subprocess.Popen('bedtools intersect -a /disk0/sm/epicv2_bed.txt -b /disk0/sm/gtfbed.txt -wa -wb', shell = True, stdout = subprocess.PIPE).stdout, sep='\t', names=['CHR', 'START', 'END', 'IlmnID', 'Name', 'chr', 'start', 'end', 'type', 'gene_name', 'gene_type'])


cpgs = pd.read_csv('/disk0/sm/methyl/03_pollution/total/age_glm_Age_678.txt', sep='\t')


probe['cpg'] = probe['Name'].str.split('_').str[0]

table = pd.merge(cpgs, df, left_on='cpg', right_on='Name',  how='left')

a = table.drop_duplicates(['cpg', 'IlmnID'])

b = pd.merge(a, probe, on='cpg', how='left').drop_duplicates(['cpg', 'IlmnID'])

c = b[['cpg', 'IlmnID', 'Relation_to_Island', 'gene_name', 'gene_type']].replace('"', '', regex=True)

c.to_csv('/disk0/sm/S678.txt', sep='\t', index=False)