import pandas as pd
import sys
from io import StringIO
import os, subprocess
from functools import reduce


KEY = sys.argv[1]
LIBRARY = sys.argv[2]


def make_recommend(KEY):
    data = pd.read_csv(subprocess.Popen("grep -v '^#' /workdir/%s/01_Genotype/AxiomGT1.calls.txt" %(KEY), shell = True, stdout = subprocess.PIPE).stdout, sep='\t', index_col=0, low_memory=False)
    recommend = open('/workdir/%s/02_SNPolisher/Recommended.ps' %(KEY), 'r').read().split('\n')
    snp = list(set(data.index) - set(recommend))
    df = data.drop(snp, axis=0)
    tdata = df.transpose()
    for key,value in tdata.iterrows():
        tdata.at[key, 'no'] = list(value).count(-1)
    tdata['recommend_callrate'] = (len(tdata.columns) - tdata['no'])/len(tdata.columns) * 100
    callrate = tdata[['recommend_callrate']]
    return callrate


def KORV2_table(KEY):
    callrate = pd.read_csv('/workdir/%s/CallRate.txt' %(KEY), sep='\t', index_col=0, low_memory=False)
    recallrate = pd.read_csv('/workdir/%s/CallRate_re.txt' %(KEY), sep='\t', index_col=0, low_memory=False)
    normalcallrate = pd.read_csv('/workdir/%s/CallRate_normal.txt' %(KEY), sep='\t', index_col=0, low_memory=False)
    gender = pd.read_csv('/workdir/%s/Gender.txt' %(KEY), sep='\t', index_col=0, low_memory=False)
    dqc = pd.read_csv('/workdir/%s/DQC_Rate.txt' %(KEY), sep='\t', index_col=0, low_memory=False)
    recommendcallrate = make_recommend(KEY)
    df = reduce(lambda x,y: pd.merge(x,y, left_index=True, right_index=True, how='inner'), [gender, callrate, recallrate, normalcallrate, recommendcallrate, dqc])
    df['cel_files'] = df.index
    df.columns = ['computed_gender', 'callrate', 're_callrate', 'normal_callrate', 'recommend_callrate', 'DQC', 'cel_files']
    sample_gender = pd.read_csv('/workdir/SAMPLE_GENDER.txt', sep='\t', low_memory=False)
    df['ID'] = df.index.str.split('_').str[6].str.split('.').str[0]
    df['plate'] = df.index.str.split('_').str[4]
    df2 = pd.merge(df, sample_gender, on='ID', how='inner')
    df2['gender_check'] = df2['gender'] == df2['computed_gender']
    table = df2[['cel_files', 'plate', 'ID', 'gender', 'computed_gender', 'gender_check', 'callrate', 're_callrate', 'normal_callrate', 'recommend_callrate', 'DQC']]
    return table


def TWO_table(KEY):
    qc_table = pd.read_csv('/workdir/%s/01_Genotype/QC_table.txt' %(KEY), sep=' ', index_col=0, low_memory=False).drop('Unnamed: 5', axis=1)
    recommendcallrate = make_recommend(KEY)
    df = reduce(lambda x,y: pd.merge(x,y, left_index=True, right_index=True, how='inner'), [qc_table, recommendcallrate])
    df['cel_files'] = df.index
    df.columns = ['computed_gender', 'callrate', 'het_rate', 'gender_ratio', 'DQC', 'recommend_callrate', 'cel_files']
    table = df[['cel_files', 'computed_gender', 'gender_ratio', 'callrate', 'recommend_callrate', 'DQC', 'het_rate']]
    return table



def KBA_table(KEY):
    callrate = pd.read_csv('/workdir/%s/CallRate.txt' %(KEY), sep='\t', index_col=0, low_memory=False)
    recallrate = pd.read_csv('/workdir/%s/CallRate_re.txt' %(KEY), sep='\t', index_col=0, low_memory=False)
    normalcallrate = pd.read_csv('/workdir/%s/CallRate_normal.txt' %(KEY), sep='\t', index_col=0, low_memory=False)
    gender = pd.read_csv('/workdir/%s/Gender.txt' %(KEY), sep='\t', index_col=0, low_memory=False)
    dqc = pd.read_csv('/workdir/%s/DQC_Rate.txt' %(KEY), sep='\t', index_col=0, low_memory=False)
    recommendcallrate = make_recommend(KEY)
    df = reduce(lambda x,y: pd.merge(x,y, left_index=True, right_index=True, how='inner'), [gender, callrate, recallrate, normalcallrate, recommendcallrate, dqc])
    df['cel_files'] = df.index
    df.columns = ['computed_gender', 'callrate', 're_callrate', 'normal_callrate', 'recommend_callrate', 'DQC', 'cel_files']
    sample_gender = pd.read_csv('/workdir/SAMPLE_GENDER.txt', sep='\t', low_memory=False)
    df['ID'] = df.index.str.split('_').str[6].str.split('.').str[0]
    df['plate'] = df.index.str.split('_').str[4]
    df2 = pd.merge(df, sample_gender, on='ID', how='inner')
    df2['gender_check'] = df2['gender'] == df2['computed_gender']
    table = df2[['cel_files', 'plate', 'ID', 'gender', 'computed_gender', 'gender_check', 'callrate', 're_callrate', 'normal_callrate', 'recommend_callrate', 'DQC']]
    return table


def run_pipeline(KEY, LIBRARY):
    if LIBRARY == 'KORV1' or LIBRARY == 'HCLV':
        table = TWO_table(KEY)
    elif LIBRARY == 'KORV2':
        table = KORV2_table(KEY)
    elif LIBRARY == 'PMDA' or LIBRARY == 'PangenomiX':
        table = TWO_table(KEY)
    elif LIBRARY == 'pharmacoFocus':
        table = TWO_table(KEY)
    elif LIBRARY == 'KBA_B':
        table = KORV2_table(KEY)
    elif LIBRARY == 'KBA_A':
        table = TWO_table(KEY)
    return table



table = run_pipeline(KEY, LIBRARY)
table.to_csv('/workdir/%s/QC_TABLE.txt' %(KEY), sep='\t', index=False)
