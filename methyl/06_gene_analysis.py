import pandas as pd
# import sys
# from io import StringIO
import os, subprocess
import modin.pandas as mpd
# from functools import reduce


## ChAmp data 활용해서 그림 만들기..
dmp = pd.read_csv('/disk0/sm/methyl/03_pollution/total/DMP_ChAMP.txt', sep='\t', low_memory=False)

dmp_pval = dmp[dmp['G11_to_G1.adj.P.Val'] <= 0.00000001]
dmp_island = dmp_pval[dmp_pval['G11_to_G1.cgi'] == 'Island']

# ## gene 세기

dmp_gene = dmp_island['G11_to_G1.gene'].value_counts()
dmp_gene.to_csv('/disk0/sm/methyl/03_pollution/total/data/dmp_gene.txt', sep='\t')

dmp_count = pd.DataFrame(dmp_pval[['G11_to_G1.cgi', 'G11_to_G1.feature']].value_counts()).unstack()
dmp_count.to_csv('/disk0/sm/methyl/03_pollution/total/data/dmp_count.txt', sep='\t')



