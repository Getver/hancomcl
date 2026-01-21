# https://link.springer.com/article/10.1186/s40168-020-00934-6
# https://github.com/ravel-lab/VALENCIA
# speciateIT


#VALENCIA
### The purpose of this tool is to classify vaginal microbial communities into community state types (CSTs)
### in a standardized and repeatable fashion, clusters are based on the 13,000+ women data
### A total of 13 CSTs are considered: CST I-A, I-B, II, III-A, III-B, IV-A, IV-B, IV-C0, IV-C1, IV-C2, IV-C3, IV-C4, V 
#### This script tests new samples based on similarity to already defined cluster centroids

import pandas as pd
import numpy as np
import re
import sys

taxonomy_path='taxonomy/taxonomy.tsv'
feature_table_path='filtered/txt/denoise_noMT_table.tsv'
reference_path='/disk0/sm/microbiome/16S/docker/DB/VALENCIA/CST_centroids_012920.csv'

# load
tax = pd.read_csv(taxonomy_path, sep="\t")

taxon_key = tax['Taxon'].str.split(';', expand=True).drop(columns=[1, 8])

# rank 컬럼 이름
taxon_key.columns = ['k','p','c','o','f','g','s']

# prefix 제거 (k__, p__ 등)
for col in taxon_key.columns:
    taxon_key[col] = taxon_key[col].str.replace(r'^[a-z]__', '', regex=True)

# index = Feature ID
taxon_key.index = tax['Feature ID']

#replacing taxon_key with 
#reading in the table of counts
counts_table = pd.read_csv(feature_table_path, sep="\t", header=1, index_col=0)

#function that determines the highest level of taxonimc specifity and then formats the condensed name
def strip_prefix(x):
    if pd.isna(x):
        return None
    x = str(x).strip()                 # ← 공백 제거 (중요)
    x = re.sub(r'^[a-z]__', '', x)     # g__, s__, f__ 등 제거
    return x

def taxon_condense(row):
    g = strip_prefix(row['g'])
    s = strip_prefix(row['s'])
    f = strip_prefix(row['f'])
    o = strip_prefix(row['o'])
    c = strip_prefix(row['c'])
    p = strip_prefix(row['p'])
    k = strip_prefix(row['k'])
    # species 우선 (VALENCIA focal taxa)
    if s is not None and g is not None:
        s_only = s.replace(f"{g} ", "")
        if g in ['Lactobacillus','Prevotella','Gardnerella','Atopobium','Sneathia','Fannyhessea']:
            return f"{g}_{s_only}"
        else:
            return f"g_{g}"
    # genus
    if g is not None:
        return f"g_{g}"
    if f is not None:
        return f"f_{f}"
    if o is not None:
        return f"o_{o}"
    if c is not None:
        return f"c_{c}"
    if p is not None:
        return f"p_{p}"
    if k is not None:
        return f"k_{k}"
    return "None"


#applying function to each row of the taxa key file
taxon_key['taxa'] = taxon_key.apply(lambda x : taxon_condense(x), axis=1)

#manual correction of names, these should be checked by looking at the ASV sequences and see how they match to the new name 
taxon_key['taxa'] = taxon_key['taxa'].replace({'g_Gardnerella':'Gardnerella_vaginalis','Lactobacillus_acidophilus/casei/crispatus/gallinarum':'Lactobacillus_crispatus'
                                                ,'Lactobacillus_fornicalis/jensenii':'Lactobacillus_jensenii','g_Escherichia/Shigella':'g_Escherichia.Shigella'
                                                ,'Lactobacillus_gasseri/johnsonii':'Lactobacillus_gasseri'})

#creating a dataframe for merging with just the information from the new condense column
taxon_merge = taxon_key[['taxa']]
#merging the counts table with the taxa table
counts_table_named = pd.merge(left=taxon_merge,right=counts_table,right_index=True,left_index=True,how="inner")
#grouping asvs with the same name and summing
counts_table_named = counts_table_named.groupby('taxa').sum()
#transposing to a table with samples are rows and counts as columns
counts_table_named = counts_table_named.T

#sorting the table by the study wide read count for each taxa
counts_table_named = counts_table_named.reindex(counts_table_named.sum().sort_values(ascending=False).index, axis=1)
#summing the read counts for each sample to be used by valencia in calculation of relative abundance
counts_table_named['read_count'] = counts_table_named.sum(axis=1)
#moving read count to first column
read_count_column = counts_table_named.pop('read_count')
counts_table_named.insert(0,'read_count',read_count_column)

#####################################################################
                    ######  ###### #######
                   #       #          #
                  #         #####     #
                   #             #    #
                    ######  #####     #
#####################################################################

#defining function to determine yue-clayton theta
def yue_distance(row, median):
    #creating a counting variable to index the median list
    taxon_count = 0
    #creating lists to iterativly store output
    median_times_obs = []
    median_minus_obs_sq = []    
    #looping through the row and calculating product and difference squared between row data and median data
    for taxon_abund in row:
        #calculate p * q
        median_times_obs.append(median[taxon_count]*taxon_abund)
        #calculate p-q squared
        median_minus_obs_sq.append((median[taxon_count]-taxon_abund)**2)
        taxon_count += 1    
    #calculate sum p* q
    product = np.nansum(median_times_obs)   
    #calculate sum p-q squared
    diff_sq = np.nansum(median_minus_obs_sq)
    #calculate yue_med_dist
    yue_med_dist = product / (diff_sq + product)
    #return the value of yue distance
    return yue_med_dist

#### defining fuction to determine the penalized similarity score, not currently in use
def penalized_simil_score(row):
    if row['subCST'] == 'I-A':
        row = row.drop('I-B_sim')
    elif row['subCST'] == 'I-B':
        row = row.drop('I-A_sim')
    elif row['subCST'] == 'III-A':
        row = row.drop('III-B_sim')
    elif row['subCST'] == 'III-B':
        row = row.drop('III-A_sim')
    row = row.drop('subCST')
    similarity_scores = list(row)
    similarity_scores.sort()
    score_len = len(similarity_scores)
    penalized_score = similarity_scores[score_len-1] * (similarity_scores[score_len-1]-similarity_scores[score_len-2]) ** (1./2)
    return penalized_score

#list of subCSTs 
CSTs = ['I-A','I-B','II','III-A','III-B','IV-A','IV-B','IV-C0','IV-C1','IV-C2','IV-C3','IV-C4','V']

reference_centroids=pd.read_csv(reference_path, sep=',')

rename_map = {'Atopobium_vaginae': 'Fannyhessea_vaginae'}
reference_centroids = reference_centroids.rename(columns={k: v for k, v in rename_map.items() if k in reference_centroids.columns})

sample_data_OG=counts_table_named
sample_data_OG.insert(0,'sampleID',sample_data_OG.index)


#checking if first two columns have appropriate headers
# if list(sample_data_OG.columns[0:2]) != ['sampleID','read_count']:
#     print('Input file expected to be a CSV with first two column headers: sampleID,read_count')

ref_cols = reference_centroids.columns
sample_data = sample_data_OG.reindex(columns=ref_cols, fill_value=0)

# set(sample_data_OG.columns).intersection(set(reference_centroids.columns))

#forcing the sample data to have the same number of columns as the reference and in the same order 
combined_data = pd.concat([sample_data_OG,reference_centroids], ignore_index=True,sort=False)
sample_data = combined_data[:-13].fillna(0).drop(['sub_CST'],axis=1)
# sample_data = sample_data
reference_centroids = combined_data.tail(13).fillna(0).drop(["sampleID","read_count"],axis=1).set_index('sub_CST')
# reference_centroids = reference_centroids

#converting all of the read counts to relative abundance data and adding back the first two columns
sample_data_rel = sample_data[sample_data.columns[2:]].div(sample_data['read_count'],axis=0)
sample_data_rel = pd.concat([sample_data[sample_data.columns[0:2]],sample_data_rel],axis=1)

# #loop measuring the similarity of each sample to each subCST centroid using yue + clayon theta
for CST in CSTs:
    target_col = CST + '_sim'
    sample_data_OG[target_col] = (sample_data_rel.apply(lambda x: yue_distance(x.iloc[2:], reference_centroids.loc[CST].values), axis=1).values)

# for CST in CSTs:
#     target_col = CST + '_sim'
    # sample_data_OG[target_col] = sample_data_rel.apply(lambda x: yue_distance(x.iloc[2:], reference_centroids.loc[CST].values),axis=1)
    
    
#outputting the acquired data with the new variability measure
#identify for each sample, which subCST was most similar, then correcting name of subCST to remove _sim
sample_data_OG['subCST'] = sample_data_OG.iloc[:,-13:].idxmax(axis=1)
sample_data_OG['subCST'] = sample_data_OG['subCST'].str.replace('_sim',"")

#applying function to calculate a penalized score, not currently in use
#sample_data_OG['penalized_score'] = sample_data_OG.apply(lambda x : penalized_simil_score(x[-14:]), axis=1)

#store the similarity between each sample and its as=ssigned 
sample_data_OG['score'] = sample_data_OG.iloc[:,-14:-1].max(axis=1)

#determine higher order CST assignment based on subCST assignment
sample_data_OG['CST'] = sample_data_OG['subCST'].replace({'I-A':'I','I-B':'I','III-A':'III','III-B':'III', 'IV-A':'IV', 'IV-B':'IV','IV-C0':'IV','IV-C1':'IV','IV-C2':'IV','IV-C3':'IV','IV-C4':'IV'})

old = ['I', 'II', 'III', 'IV', 'V']
new = ['CST1', 'CST2', 'CST3', 'CST4', 'CST5']


for i in range(0, len(new)):
    sample_data_OG['CST'].replace(old[i], new[i], inplace=True)


#output the assignments in new CSV file
sample_data_OG.to_csv('CST_VALENCIA.txt', sep='\t',index=None)





# #plotting distributions of similarity scores for each subCST agaist that for the reference dataset
import matplotlib.pyplot as plt

# defining figure
similarity_fig, similarity_axs = plt.subplots(
    1, 1, figsize=(10, 6), facecolor='w', edgecolor='k'
)
similarity_fig.subplots_adjust(
    left=0.15, right=0.85, bottom=0.15, top=0.9
)

# x-axis location
loc = 0
CST_labels = []

# reference average and std (VALENCIA paper)
ref_ave = {
    'I-A':0.995678,'I-B':0.858443,'II':0.811734,
    'III-A':0.972983,'III-B':0.810466,
    'IV-A':0.718100,'IV-B':0.659592,
    'IV-C0':0.321432,'IV-C1':0.745991,
    'IV-C2':0.695203,'IV-C3':0.758735,
    'IV-C4':0.681502,'V':0.734973
}

ref_std = {
    'I-A':0.007615,'I-B':0.151048,'II':0.161987,
    'III-A':0.058698,'III-B':0.133990,
    'IV-A':0.150829,'IV-B':0.162577,
    'IV-C0':0.137558,'IV-C1':0.231891,
    'IV-C2':0.246322,'IV-C3':0.213653,
    'IV-C4':0.186462,'V':0.140276
}

# CST color scheme
CST_color_scheme = {
    'I-A':'#ff6868','I-B':'#ffd4da','II':'#b4ff68',
    'III-A':'#ffbc6b','III-B':'#e4a67b',
    'IV-A':'#c1adec','IV-B':'#91a8ed',
    'IV-C0':'#989898','IV-C1':'#ffc0cb',
    'IV-C2':'#a8e5e5','IV-C3':'#9acc9a',
    'IV-C4':'#800080','V':'#ffff71'
}

# plotting
for CST in CSTs:
    subset = sample_data_OG[sample_data_OG['subCST'] == CST]
    if subset.shape[0] == 0:
        continue
    boxprops = dict(linewidth=1, color="k")
    medianprops = dict(linewidth=1, color="k")
    # boxplot of sample similarity scores
    similarity_axs.boxplot(
        subset['score'],
        positions=[loc - 0.25],
        notch=True,
        widths=[0.25],
        patch_artist=True,
        boxprops=boxprops,
        medianprops=medianprops
    )
    # reference average ± std
    similarity_axs.scatter(
        loc,
        ref_ave[CST],
        s=20,
        marker='d',
        color=CST_color_scheme[CST],
        zorder=3
    )
    similarity_axs.plot(
        [loc, loc],
        [ref_ave[CST] - ref_std[CST], ref_ave[CST] + ref_std[CST]],
        color=CST_color_scheme[CST],
        linewidth=2
    )
    CST_labels.append(CST)
    loc += 1

# labels and limits
similarity_axs.set_xticks(range(len(CST_labels)))
similarity_axs.set_xticklabels(CST_labels, rotation=45)
similarity_axs.set_ylim(0, 1)
similarity_axs.set_xlabel("subCST")
similarity_axs.set_ylabel("Similarity to assigned subCST")
similarity_axs.set_title("VALENCIA similarity score distributions")

plt.savefig(
    "VALENCIA_similarity_distribution.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()



## 결과 비교


valencia_result = sample_data_OG[['CST', 'score', 'subCST']].replace('')
sm_result = pd.read_json('Final_result.json').T

sm_result['sm_CST'] = sm_result['info'].apply(lambda x: x.get('CT') if isinstance(x, dict) else None)
sm_result['sm_score'] = sm_result['info'].apply(lambda x: x.get('DS', {}).get('CST', {}).get('score') if isinstance(x, dict) else None)

sm_result = sm_result[['sm_CST', 'sm_score']]

sm_result['sm_CST'].replace('CST4-BV', 'CST4', inplace=True)
sm_result['sm_CST'].replace('CST4-AV', 'CST4', inplace=True)


result = pd.merge(sm_result, valencia_result, left_index=True, right_index=True, how='inner')


result['check'] = result['sm_CST'] == result['CST']

print(result.loc[result['check'] == False, ])


result.to_csv('tmp_check.txt', sep='\t')



