import pandas as pd 
import modin.pandas as mpd
import os

# os.environ["MODIN_ENGINE"] = "dask"  # Modin will use Dask

# /disk0/sm/methyl/03_pollution/total/mid_report/two-way_normalized_beta_value_all_CpGs.txt
beta = mpd.read_csv('/disk0/sm/methyl/03_pollution/total/two-way_normalized_beta_value_all_CpGs.txt', sep='\t', low_memory=False, index_col=0)

## minfi
dmp = pd.read_csv('/disk0/sm/methyl/03_pollution/total/DMP_Minfi_dmpFinder.txt', sep='\t', low_memory=False)
dmp_good_100 = dmp.iloc[0:100, ]
beta_good_100 = beta.loc[dmp_good_100.index,]

## champ
dmp = pd.read_csv('/disk0/sm/methyl/03_pollution/total/DMP_ChAMP.txt', sep='\t', low_memory=False)
beta_tmp = beta
beta_tmp.index = beta.index.str.split('_').str[0]
beta_good_100 = beta_tmp.loc[dmp['G11_to_G1.Name'][:100],]
beta_good_100 = beta_good_100.loc[~beta_good_100.index.duplicated(keep='first')]


beta_good_100.to_csv('/disk0/sm/methyl/03_pollution/total/data/headmap_beta_100.txt', sep='\t')
