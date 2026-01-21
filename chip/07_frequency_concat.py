import json, pandas as pd
# import modin.pandas as pd
import sys
from io import StringIO
import os, subprocess

## 아뉘 .frq 파일 \t으로 구분 안되어있어서 따로 처리 해줘야하고 sample ID 는 그때그때 상황에 맞게 고쳐가면서 활용해야함
KEY='576_PangenomiX'
FOLDER='MACROGEN_1056/'+KEY


## table1은 /disk0/sm/PMDA/minjin/pharmCAToutput에서 phenotype.json파일에서 컬럼 가져와서 사용(알아봐야함)
# json1 = pd.read_json('/disk0/sm/PMDA/minjin/test.json', orient='value')

def make_df(json1, gene, key):
    df = pd.DataFrame()
    for g in gene:
        tmp = pd.DataFrame.from_dict(data=json1['geneReports'][key][g]['variants'], orient='columns')
        tmp_df = pd.concat([df, tmp], axis=0)
        df = tmp_df
    return df

def read_json(path, sampleID):
    with open('%s' %(path)) as file:
        json1 = json.load(file)
    #
    dpwg_gene = list(json1['geneReports']['DPWG'].keys())
    cpic_gene = list(json1['geneReports']['CPIC'].keys())
    df_dpwg = make_df(json1, dpwg_gene, 'DPWG')
    df_cpic = make_df(json1, cpic_gene, 'CPIC')
    df_dpwg['sampleID'] = sampleID
    df_cpic['sampleID'] = sampleID
    return df_dpwg, df_cpic

pheno_json = pd.read_csv(subprocess.Popen("ls /disk0/sm/PMDA/"+FOLDER+"/minjin/pharmCAToutput/*.phenotype.json", shell = True, stdout = subprocess.PIPE).stdout, sep='\t', low_memory=False, names=['file'])
pheno_json['sampleID'] = pheno_json['file'].str.split('/').str[-1].str.split('.').str[1].str.split('_').str[0]  #.str.split('.').str[0]

dpwg = pd.DataFrame()
cpic = pd.DataFrame()

for _, value in pheno_json.iterrows():
    tmp_dpwg, tmp_cpic= read_json(value[0], value[1])
    tmp2_dpwg = pd.concat([dpwg, tmp_dpwg], axis=0)
    tmp2_cpic = pd.concat([cpic, tmp_cpic], axis=0)
    dpwg = tmp2_dpwg.reset_index(drop=True)
    cpic = tmp2_cpic.reset_index(drop=True)

dpwg = dpwg.loc[dpwg.call != dpwg['wildtypeAllele']+'/'+dpwg['wildtypeAllele']].dropna(subset='call')
cpic = cpic.loc[cpic.call != cpic['wildtypeAllele']+'/'+cpic['wildtypeAllele']].dropna(subset='call')

## table2 제작
## frequency concat
name = pd.read_csv(subprocess.Popen("ls /disk0/sm/PMDA/frequency", shell = True, stdout = subprocess.PIPE).stdout, sep='\t', low_memory=False, names=['gene_file'])
name['gene'] = name['gene_file'].str.split('.').str[0]
gene = name['gene'].to_list()

tmp_df = pd.DataFrame()
for g in gene:
    tmp_frequency = pd.read_csv('/disk0/sm/PMDA/frequency/%s.txt' %(g), sep='\t', low_memory=False)
    tmp_frequency['Gene'] = g
    tmp_frequency.columns = ['star allele'] + list(tmp_frequency)[1:]
    tmp_df2 = pd.concat([tmp_df, tmp_frequency], axis=0)
    tmp_df = tmp_df2

# tmp_df.to_csv('/disk0/sm/PMDA/frequency.txt', sep='\t', index=False)

drug = pd.read_csv('/disk0/sm/PMDA/drug_merge.txt', sep='\t', low_memory=False)
frequency = tmp_df

df2 = pd.merge(drug, frequency, on=['Gene', 'star allele'], how='inner')
table2 = df2[['Gene', 'star allele', 'pharmvar', 'protein', 'position(GRCh38)', 'position(RefSeqGene)', 'rsID', 'ref', 'alt',
            'Activity Value (Optional)', 'Allele Biochemical Functional Status (Optional)', 'Allele Clinical Functional Status (Required)', 
            'Allele Clinical Function Substrate Specificity (Optional)', 'References (Required)', 'Strength of Evidence (Required)',
            'Summary of Findings (Required)', 'Comments', 'African American/Afro-Caribbean', 'Central/South Asian', 'East Asian', 'European',
            'Latino', 'Near Eastern', 'Oceanian', 'Sub-Saharan African']]


## table3 제작

frq = pd.read_csv('/disk0/sm/PMDA/'+FOLDER+'/minjin/plink/S2.frq', sep='\t', low_memory=False)
frq['AF'] = 1 - frq['MAF']
frq['CHR'] = 'chr'+frq['CHR'].astype(str)

pan_db = pd.read_csv(subprocess.Popen('zcat /disk0/sm/PMDA/'+FOLDER+'/03_vcf/'+KEY+'.vcf.gz | grep -v "^##"', shell = True, stdout = subprocess.PIPE, text=True).stdout, sep='\t', low_memory=False)
pan_db['rsID'] = pan_db['INFO'].replace('RSID=', '', regex=True)
pan = pan_db[['ID', 'rsID']]
frq_rs = pd.merge(frq, pan, left_on='SNP', right_on='ID', how='left')

table3 = frq_rs


## table merge

def merge_table(table1, table2, table3):
    COL = table1.columns
    for key, value in table1.iterrows():
        for i in range(len(value[5])):
            table1.at[key, str(i)] = value[5][i]
    mable1 = table1.melt(id_vars=COL)
    mable1['alleles'] = mable1['value']
    table1 = mable1.drop(['variable', 'value'], axis=1)
    tmp1 = pd.merge(table1, table2, left_on=['gene', 'dbSnpId', 'alleles'], right_on=['Gene', 'rsID', 'star allele'], how='left').drop(['Gene', 'star allele', 'rsID'], axis=1)
    tmp2 = pd.merge(tmp1, table3, left_on=['chromosome', 'dbSnpId'], right_on=['CHR', 'rsID'], how='left').drop(['ID', 'CHR', 'rsID', 'SNP'], axis=1)
    tmp3 = tmp2[['gene', 'chromosome', 'position', 'dbSnpId', 'alleles', 'call'] + list(tmp2.columns[6:])]
    new_name = list(['Gene', 'chromosome', 'position', 'rsID', 'starallele'] + list(tmp3.columns[5:11]+'_pharmcat') + list(tmp3.columns[11:33]+ '_pharmgkb') + list(tmp3.columns[33:]+ '_freq'))
    tmp3.columns = new_name
    return tmp3

D = merge_table(dpwg, table2, table3)
C = merge_table(cpic, table2, table3)

D.to_csv('/disk0/sm/PMDA/'+FOLDER+'/minjin/TABLE_PharmGKB_DPWG.txt', sep='\t', index=False)
C.to_csv('/disk0/sm/PMDA/'+FOLDER+'/minjin/TABLE_PharmGKB_CPIC.txt', sep='\t', index=False)

