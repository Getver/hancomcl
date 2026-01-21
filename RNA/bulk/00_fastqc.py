# import modin.pandas as mpd
import pandas as pd
import sys
from io import StringIO
import os, subprocess
from functools import reduce

Folder = sys.argv[1]
Folder='/disk0/sm/bulk/2024/1114'

readQ = pd.read_csv(Folder+'/total_sequence.txt', sep=' ', names=['ID', 'Total Bases', 'Read Length', 'Total Reads', 'GC Ratio(%)']).drop_duplicates('ID').reset_index(drop=True)
atgQ = pd.read_csv(Folder+'/total_qulity.txt', sep=' ', names=['ID', 'Q20_R1', 'Q30_R1', 'Q20_R2', 'Q30_R2'])

table = pd.merge(readQ, atgQ, on='ID', how='inner')

table['Total Bases(Gb)'] = round(table['Total Bases'] / 1000000000, 2)
table['>Q20 Reads(%)'] = round((table['Q20_R1'] + table['Q20_R2'])/2, 2)
table['>Q30 Reads(%)'] = round((table['Q30_R1'] + table['Q30_R2'])/2, 2)
df = table[['ID', 'Read Length', 'Total Reads', 'Total Bases', 'Total Bases(Gb)', 'GC Ratio(%)', '>Q20 Reads(%)', '>Q30 Reads(%)']]

df.to_csv(Folder+'/total_QCtable.txt', sep='\t', index=False)
