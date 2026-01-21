import pandas as pd

data = pd.read_csv('/disk0/sm/methyl/00_input.txt', sep='\t')
data.to_csv('/disk0/sm/methyl/00_output.txt', sep=',', index=False)
data.info()
