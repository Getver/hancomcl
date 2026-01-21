# import sys
# import re
import pandas as pd
from functions import *
# from io import StringIO
import os

# https://covariants.org/   covid virus reference data

# HOME_DIR = '/home/covid/'
HOME_DIR = '/home/root/'
REF_DIR = HOME_DIR + 'reference/'
SEQ_DIR = HOME_DIR + 'sequence/'
WORK_DIR = HOME_DIR + 'work/'

df_position, df_name, df_vari, df_pheno, df_tmp = manipulate_bundle_files()

reference_fasta = process_reference()

## fasta폴더의 spike protein부분 아미노산 서열 만들기
references = os.popen("ls " + SEQ_DIR + "fasta/* | awk -F '/' '{ print $NF }' | awk -F '.' '{print $1}' ").read().strip().split('\n')
query_dir = SEQ_DIR + 'fasta'
for seq_id in references:
    query_fasta, sample, file_name = run_alignment(seq_id, query_dir, reference_fasta)
    target = query_fasta[df_position.at['S:', 'Start']-1 : df_position.at['S:', 'End']]
    spike_amino = run_translate(target)
    amino = ['>'+sample, spike_amino]
    amino.append('>'+sample)
    multi_alignment_amino = open(SEQ_DIR + 'amino/' + sample + '_amino.aa', 'w')
    for i in amino:
        multi_alignment_amino.write(i+'\n')
    multi_alignment_amino.close()
    name = seq_id.split('/')[-1]

