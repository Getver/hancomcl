import pandas as pd
import sys
from io import StringIO
import os, subprocess
from functools import reduce

## 예상 나이
age_total = pd.read_csv('/disk0/sm/methyl/03_pollution/total/DNAmAge_predict.txt', sep='\t', low_memory=False)
# age_total = pd.read_csv('/disk0/sm/methyl/03_pollution/total/data/DNAmAge_predict_678_2024.txt.txt', sep='\t', low_memory=False)

age_total['Name'] = age_total['id'].str.split('_').str[0]
age_total['phase'] = age_total['id'].str.split('_').str[1]
age_total['phase'].replace("1", '1기', inplace=True)
age_total['phase'].replace("11", '11기', inplace=True)
age = age_total[['id', 'Name', 'phase', 'Horvath', 'Hannum', 'Levine', 'skinHorvath', 'PedBE', 'Wu', 'BLUP', 'EN', 'SM']]

age.to_csv('/disk0/sm/methyl/03_pollution/total/data/age.txt', sep='\t', index=False)

# folder = '01'

## 11기와 1기 나이차이 구하기
# age = pd.read_csv('/disk0/sm/methyl/03_pollution/total/data/age.txt', sep='\t')
# age.set_index('Name', inplace=True)

# age_11
age_11 = age[age.phase == '11기']
age_11 = age_11.drop(['phase', 'id'], axis=1)

# age_1
age_1 = age[age.phase == '1기']
age_1 = age_1.drop(['phase', 'id'], axis=1)

# age 1_11
age_1 = age_1.sort_values('Name').set_index('Name')
age_11 = age_11.sort_values('Name').set_index('Name')
age_1_11 = age_11 - age_1
age_1_11.to_csv('/disk0/sm/methyl/03_pollution/total/data/age_1_11.txt', sep='\t')



## 나이차이비교

actual = pd.read_csv('/disk0/sm/methyl/2023/01_1760/sample_sheet.csv', sep=',', low_memory=False, header=6)[['Sample_Name', 'Age']]
age = pd.read_csv('/disk0/sm/methyl/2023/01_1760/data/age.txt', sep='\t', low_memory=False)

tmp = pd.merge(actual, age, left_on='Sample_Name', right_on='id', how='inner')
tmp = tmp.set_index('id')

tmp.to_csv('/disk0/sm/methyl/2023/01_1760/data/age_table.txt', sep='\t', index=False)

age_coef = pd.DataFrame(index=tmp.index)

age_coef['Horvath'] = tmp['Horvath'] - tmp['Age']
age_coef['Hannum'] = tmp['Hannum'] - tmp['Age']
age_coef['Levine'] = tmp['Levine'] - tmp['Age']
age_coef['skinHorvath'] = tmp['skinHorvath'] - tmp['Age']
age_coef['PedBE'] = tmp['PedBE'] - tmp['Age']
age_coef['Wu'] = tmp['Wu'] - tmp['Age']
age_coef['BLUP'] = tmp['BLUP'] - tmp['Age']
age_coef['EN'] = tmp['EN'] - tmp['Age']
age_coef['SM'] = tmp['SM'] - tmp['Age']

age_coef.to_csv('/disk0/sm/methyl/03_pollution/total/data/age_algorithm.txt', sep='\t')





# # 예상
# age_1720 = pd.concat([age_1, age_11], axis=0)
# age_1720.to_csv('/disk0/sm/methyl/total/data/age_292.txt', sep='\t')





# # age_5.index = age_5.index.str.split('_').str[0]
# # age_11.index = age_11.index.str.split('_').str[0]


# # strange
# age['strange'] = age['EN'] - age['Actual age']

# age_dif = pd.DataFrame(index=age_1.index)

# age_dif['Actual age'] = age_1['Actual age']
# age_dif['1_sub'] = age_1['strange']
# age_dif['11_sub'] = list(age_11['strange'])
# age_dif['strange'] = age_dif['11_sub'] - age_dif['1_sub']

# age_dif.to_csv('/disk0/sm/methyl/total/data/age_EN.txt', sep='\t')
# age



