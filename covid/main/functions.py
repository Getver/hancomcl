#!/usr/bin/env python3

import os
import re
import subprocess
import shutil
import pandas as pd

from Bio import SeqIO

HOME_DIR = '/home/root/'
REF_DIR = os.path.join(HOME_DIR, 'reference/')
SEQ_DIR = os.path.join(HOME_DIR, 'sequence/')
WORK_DIR = os.path.join(HOME_DIR, 'work/')
TMP_DIR = os.path.join(HOME_DIR, 'tmp/')

BUNDLE_DIR = '/files/'

def create_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def process_reference():
    with open(os.path.join(REF_DIR, 'Reference_MN908947.fa'), 'r') as reference_file:
        return ''.join(line.rstrip() for line in reference_file if not line.startswith('>'))

def manipulate_bundle_files():
    position_file = 'corona_position.txt'
    name_file = 'corona_name.txt'
    vari_file = 'corona_vari.txt'
    pheno_file = 'corona_pheno.txt'

    df_position = pd.read_csv(os.path.join(BUNDLE_DIR, position_file), sep='\t', low_memory=True, index_col=0)
    df_name = pd.read_csv(os.path.join(BUNDLE_DIR, name_file), sep='\t', low_memory=True).fillna(' ')
    df_vari = pd.read_csv(os.path.join(BUNDLE_DIR, vari_file), sep='\t', low_memory=True)
    df_pheno = pd.read_csv(os.path.join(BUNDLE_DIR, pheno_file), sep='\t', low_memory=True).to_dict('list')

    df_name.insert(0, 'Sample', '')
    df_name['Sample'] = df_name[df_name.columns[1:]].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
    df_vari.columns = df_name['Sample'].to_list()
    df_vari = df_vari.to_dict('list')
    df_tmp = pd.DataFrame(columns = list(df_name.keys()) + ['NTD', 'RBD', 'FCS'] + list(df_pheno.keys()) + ['No phenotypes', 'Total epitopes'])

    return df_position, df_name, df_vari, df_pheno, df_tmp

def run_translate(seq):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        '---':'-'
    }

    ambiguous_nucleotides = {'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N', '.'}
    protein = ""

    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]

        if '-' in codon and codon != '---':
            codon = codon.replace('-', 'N')

        if set(codon) & ambiguous_nucleotides:
            protein += '.'
        else:
            protein += codon_table.get(codon, '.')

    return protein

def retrieve_target(sequence, df_position, region):
    return sequence[int(df_position.loc[region, 'Start'])-1 : int(df_position.loc[region, 'End'])]

def parse_blast_output(output, reference_length):
    query_amino = ['.'] * (reference_length + 100)
    sample, file_name = None, None

    for line in output.split('\n'):
        if line.startswith('Query= '):
            sample = re.split(r'[ =]', line)[-1].split(' ')[0]
            file_name = sample.split(' Severe')[0]
        elif line.startswith('Sbjct '):
            position = int(re.split(r'\s+', line)[1])
            sequence = re.split(r'\s+', line)[2].rstrip()
            for i, amino in enumerate(sequence):
                query_amino[position + i] = amino

    query_fasta = ''.join(query_amino[1:-99])
    return query_fasta, sample, file_name

def run_blast(seq_id, query_dir):
    fasta_file = os.path.join(HOME_DIR, f'{query_dir}/{seq_id}.fasta')

    blast_cmd = f'blastn -db Wuhan_Hu_1 -query {fasta_file}'
    result = subprocess.run(blast_cmd, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

    return result.stdout

def parse_blast_output(blast_output, reference_fasta):
    query_amino = ['.' for _ in range(len(reference_fasta) + 100)]
    sample = ""
    file_name = ""
    lines = blast_output.strip().split('\n')

    for line_num, line in enumerate(lines):
        if line.startswith('Query= '):
            sample, file_name = extract_sample_info(line)
        if line.startswith('Sbjct '):
            query_amino = update_query_amino(query_amino, lines, line, line_num)

    query_amino = query_amino[1:-99]
    query_fasta = ''.join(query_amino)
    return query_fasta, sample, file_name

def extract_sample_info(line):
    sample = line.strip().split('= ')[1].split(' ')[0].replace('(', '').replace(')', '')
    file_name = line.strip().split('= ')[1].split(' Severe')[0].replace('(', '').replace(')', '')
    return sample, file_name

def update_query_amino(query_amino, lines, line, line_num):
    a = 0
    subject_start = int(re.split(r'\s+', line)[1].rstrip())
    for k in re.split(r'\s+', lines[line_num - 2])[2].rstrip():
        query_amino[subject_start + a] = k
        a += 1
    return query_amino

def run_alignment(seq_id, query_dir, reference_fasta):
    blast_output = run_blast(seq_id, query_dir)
    query_fasta, sample, file_name = parse_blast_output(blast_output, reference_fasta)
    return query_fasta, sample, file_name

def retrieve_target_and_translate(sequence, df_position, region):
    target = sequence[int(df_position.loc[region, 'Start'])-1 : int(df_position.loc[region, 'End'])]
    return run_translate(target.upper())

def find_mutations(query_amino, reference_amino):
    mutations = [
        f'{ref}{i + 1}{query}'
        for i, (ref, query) in enumerate(zip(reference_amino, query_amino))
        if query != '.' and ref != query
    ]
    n_index = [i for i, query in enumerate(query_amino) if query == '.']

    return mutations, n_index

def merge_deletions(mutations):
    merged = []
    tmp = []

    for i in range(len(mutations)):
        if i > 0 and int(mutations[i][1:-1]) - int(mutations[i-1][1:-1]) == 1:
            tmp.extend([i, i-1])
            tmp = list(set(tmp))
        
        if tmp and (i == len(mutations) - 1 or int(mutations[i][1:-1]) - int(mutations[i-1][1:-1]) != 1):
            if '-' in mutations[tmp[0]] or '-' in mutations[tmp[-1]]:
                mutations[tmp[-1]] = mutations[tmp[-1]][:-1] + mutations[tmp[0]][-1]
                mutations[tmp[0]] = mutations[tmp[0]][:-1] + '-'
            tmp = []

    merged = [mut for mut in mutations if mut[0] != mut[-1]]
    return merged

def run_mutation(query_amino, reference_amino, region, sample):
    mutations, _ = find_mutations(query_amino, reference_amino)
    merged_mutations = merge_deletions(mutations)
    return [f'{region}{mutation}' for mutation in merged_mutations]

def run_corona(df_vari, df_pheno, df_name, mutations, sample, df):
    calculation = [
        sum(1 for mutation in mutations if mutation in df_vari[key]) / len(df_vari[key])
        for key in df_vari
    ]

    corona = max(df_vari.keys(), key=lambda key: calculation[list(df_vari.keys()).index(key)])

    df.at[sample, 'Sample'] = sample
    df.at[sample, 'Total epitopes'] = ' '.join(mutations)

    for i, key in enumerate(list(df_name.keys())[1:], start=1):
        df.at[sample, key] = corona.split('|')[i - 1]

    smuts = [mutation[2:] for mutation in mutations if 'S:' in mutation]

    for smut in smuts:
        matched = False
        for key, value in df_pheno.items():
            if smut in value:
                df.at[sample, key] = f"{df.at[sample, key]}{smut} ".replace('nan', '')
                matched = True
        if not matched:
            df.at[sample, 'No phenotypes'] = f"{df.at[sample, 'No phenotypes']}{smut} ".replace('nan', '')

    for smut in smuts:
        smut_pos = int(smut[1:-1])
        if 14 <= smut_pos <= 306:
            df.at[sample, 'NTD'] = f"{df.at[sample, 'NTD']}{smut} ".replace('nan', '')
        elif 331 <= smut_pos <= 528:
            df.at[sample, 'RBD'] = f"{df.at[sample, 'RBD']}{smut} ".replace('nan', '')
        elif 677 <= smut_pos <= 690:
            df.at[sample, 'FCS'] = f"{df.at[sample, 'FCS']}{smut} ".replace('nan', '')

    return df

def run_additional(variant_label, seq_dir, reference_fasta, df_position, query_amino, reference_amino, df_total, sample):
    classified_fasta, _, _ = run_alignment(variant_label, seq_dir, reference_fasta)
    classified_amino = retrieve_target_and_translate(classified_fasta, df_position, "S:")

    df_total = update_rates(query_amino, reference_amino, df_total, sample, 'Wuhan_Hu')
    df_total = update_rates(query_amino, classified_amino, df_total, sample, 'Reference')

    return df_total

def update_rates(query, subject, df_total, sample, label):
    regions = {'S protein': (0, 1274), 'NTD': (13, 306), 'RBD': (330, 528), 'FCS': (676, 690)}

    for region, (start, end) in regions.items():
        query_segment = query[start:end]
        subject_segment = subject[start:end]

        coverage = calculate_coverage(query_segment, subject_segment)
        identity = calculate_identity(query_segment, subject_segment)

        df_total.at[sample, f'{label}:Identities({region})'] = identity
        df_total.at[sample, f'{label}:Coverage({region})'] = coverage

    return df_total

def calculate_coverage(query_segment, subject_segment):
    query_cleaned = query_segment.translate({ord(c): None for c in '-.'})
    coverage_percentage = (len(query_cleaned) / len(subject_segment)) * 100
    return f'{len(query_cleaned)}/{len(subject_segment)}({coverage_percentage:.2f}%)'

def calculate_identity(query_segment, subject_segment):
    matches = sum(1 for q, s in zip(query_segment, subject_segment) if q == s and q not in '-.')
    identity_percentage = (matches / len(subject_segment)) * 100
    return f'{matches}/{len(subject_segment)} ({identity_percentage:.2f}%)'

def create_alignment_file(sample, df, query_amino, file_name, seq_dir, ref_dir):
    who_label = df.at[sample, 'WHO Label']
    file_list = pd.read_csv(
        subprocess.Popen(f'ls {os.path.join(seq_dir, "amino/")}', shell=True, stdout=subprocess.PIPE).stdout, 
        sep='\t', names=['file']
    )['file'].str.replace('_amino.aa', '', regex=True).tolist()

    if who_label not in file_list:
        who_label = 'Alpha'

    file_path = subprocess.run(
        f'ls {os.path.join(seq_dir, f"amino/*{who_label}*")}', shell=True, check=True, text=True, stdout=subprocess.PIPE
    ).stdout.split('\n')[0]

    os.system(f'cp -pr {os.path.join(ref_dir, "Reference_MN908947_amino.aa")} {os.path.join(TMP_DIR, "alignment/")}')
    os.system(f'cp -pr {file_path} {os.path.join(TMP_DIR, "alignment/")}')

    amino_sequence = ['>' + file_name] + [
        query_amino[i:i+60] for i in range(0, len(query_amino), 60)
    ]
    
    alignment_file_path = os.path.join(TMP_DIR, 'alignment/', f'{sample}.aa')
    with open(alignment_file_path, 'w') as f:
        f.write('\n'.join(amino_sequence) + '\n')

    return file_path

def run_jalview(ref_dir, file_path, sample, work_dir):
    create_directory(work_dir)

    try:
        subprocess.run([
            'jalview', '--open', os.path.join(ref_dir, 'region.aa'), '--colour', 'gecos-blossom', 
            '--append', os.path.join(ref_dir, 'Reference_MN908947_amino.aa'), '--wrap', 
            '--append', file_path, os.path.join(TMP_DIR, 'alignment/', f'{sample}.aa'), 
            '--image', os.path.join(work_dir, f'{sample}.png'), '--scale', '4', '--overwrite'
        ], timeout=60)
        return True
    except subprocess.TimeoutExpired:
        return False
    
def run_aminofile(sample, df, query_amino, file_name, seq_dir, ref_dir, work_dir):
    file_path = create_alignment_file(sample, df, query_amino, file_name, seq_dir, ref_dir)
    
    if not run_jalview(ref_dir, file_path, sample, work_dir):
        if not run_jalview(ref_dir, file_path, sample, work_dir):
            return False
    return True

def split_fasta(input_fasta, query_dir=os.path.join(HOME_DIR, "input_dir")):
    create_directory(query_dir)

    reference_fasta = process_reference()
    run_name = os.path.basename(input_fasta).split('.')[0]

    df_position, df_name, df_vari, df_pheno, df_tmp = manipulate_bundle_files()

    for _, record in enumerate(SeqIO.parse(input_fasta, "fasta")):
        seq_id = record.id

        output_fasta = os.path.join(query_dir, f"{seq_id}.fasta")
        with open(output_fasta, "w") as f:
            SeqIO.write(record, f, "fasta")
        
        query_fasta, sample, file_name = run_alignment(seq_id, query_dir, reference_fasta)
        
        mutation = []

        refer_S, query_S = [], []

        for idx in df_position.index.to_list():
            reference_amino = retrieve_target_and_translate(reference_fasta, df_position, idx)
            query_amino = retrieve_target_and_translate(query_fasta, df_position, idx)

            if idx == "S:":
                refer_S, query_S = reference_amino, query_amino

            mutation.extend(run_mutation(query_amino, reference_amino, idx, sample))
        
        df_total = run_corona(df_vari, df_pheno, df_name, mutation, sample, df_tmp)

        variant_label = df_total.at[sample, 'WHO Label'] # seq_id
        seq_dir = 'sequence/fasta/' # query_dir

        df = run_additional(variant_label, seq_dir, reference_fasta, df_position, query_S, refer_S, df_total, sample)

        pass_status = run_aminofile(sample, df, query_S, file_name, SEQ_DIR, REF_DIR, WORK_DIR)

        if not pass_status:
            shutil.move(f"input_dir/{seq_id}.fasta", os.path.join(TMP_DIR, f'{seq_id}.fasta'))
            print()
            print('######################################################')
            print(f'{seq_id} Failed !!')
            print('######################################################')
            print()

    df.to_csv(os.path.join(WORK_DIR, f'{run_name}.txt'), sep='\t', index=False)

    command = 'ls ' + TMP_DIR + '*.fasta'
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


    if result.returncode == 0:
        os.system('cat ' + TMP_DIR + '*.fasta > ' + WORK_DIR + 'failed.fasta')
        print("Failed multiple alignment for some sequences.")
        print('Check "failed.fasta" file')
    else:
        print("All done !!")

    os.system('rm -rf ' + TMP_DIR + '*.fasta')
    os.system('rm -rf ' + query_dir)
