import pandas as pd
import sys
from io import StringIO
import os, subprocess
from functools import reduce
from sklearn.metrics import mean_squared_error


## 예상 나이
age_total = pd.read_csv('/disk0/sm/methyl/clock/mydata/DNAmAge_predict_678_2024.txt', sep='\t', low_memory=False)
age = age_total[['id', 'Horvath', 'Hannum', 'Levine', 'skinHorvath', 'PedBE', 'Wu', 'BLUP', 'EN', 'SM']] #, 'BNN', 'TL', 

## 실제나이

actual = pd.read_csv('/disk0/sm/methyl/2024/03_pollution/total/total_samplesheet_KoBBID.txt', sep='\t', low_memory=False)
actual = actual[actual.Pool_ID == 2024]
actual['id'] = actual['ID']

tmp = pd.merge(actual, age, on='id', how='inner')
tmp = tmp.set_index('id')

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


age_coef.SM.mean()
age_coef.SM.std()
mean_squared_error(tmp.Age, tmp.SM)




## SM만 6개 그리기

age1 = pd.read_csv('./clock/mydata/DNAmAge_predict_813.txt', sep='\t', low_memory=False)[['id', 'Horvath', 'Levine', 'SM']]
age2 = pd.read_csv('./clock/mydata/DNAmAge_predict_678.txt', sep='\t', low_memory=False)[['id', 'BLUP', 'Hannum', 'SM']]
age3 = pd.read_csv('./clock/TEST_glmnet/DNAmAge_predict_S2106.txt', sep='\t', low_memory=False)[['id', 'SM']]
age4 = pd.read_csv('./clock/TEST_glmnet/DNAmAge_predict_S2106_n015_2000.txt', sep='\t', low_memory=False)[['id', 'SM']]
age5 = pd.read_csv('./clock/TEST_glmnet/DNAmAge_predict_dmp_qval_5.txt', sep='\t', low_memory=False)[['id', 'SM']]

age1.columns = ['id', 'Horvath', 'Levine', 'S813']
age2.columns = ['id', 'BLUP', 'Hannum', 'S678'] 
age3.columns = ['id', 'S859']
age4.columns = ['id', 'S709']
age5.columns = ['id', 'S669']


# actual = pd.read_csv('/home/data/test/test2_age.txt', sep='\t', low_memory=False)
# actual['id'] = age1['id']

actual = pd.read_csv('/home/data/work/sample_sheet_new.csv', sep=',', low_memory=False, header=6)
actual = actual[actual.Pool_ID == 2023]
actual['id'] = actual['Sample_Name']

tmp = reduce(lambda x,y: pd.merge(x,y, on='id', how='inner'), [actual, age1, age2, age3, age4, age5])


age_coef = pd.DataFrame(index=tmp.index)

age_coef['Horvath'] = tmp['Horvath'] - tmp['Age']
age_coef['Hannum'] = tmp['Hannum'] - tmp['Age']
age_coef['Levine'] = tmp['Levine'] - tmp['Age']
age_coef['BLUP'] = tmp['BLUP'] - tmp['Age']
age_coef['S813'] = tmp['S813'] - tmp['Age']
age_coef['S678'] = tmp['S678'] - tmp['Age']
age_coef['S859'] = tmp['S859'] - tmp['Age']
age_coef['S709'] = tmp['S709'] - tmp['Age']
age_coef['S669'] = tmp['S669'] - tmp['Age']


age_coef.to_csv('/home/data/age_algorithm_check1.txt', sep='\t', index=False)




