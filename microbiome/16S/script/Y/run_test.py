import pandas as pd
import json
# import transfer_rank as TR
# import preprocessing_data as PD
import make_dictionary as MD
# import return_result as RR
# import calculation_score as CS
import pipeline as P
import sys
# /snap/bin/mmdc -i pipeline_flow.mmd -o pipeline.png

## 내가 만든 테이블
## 미생물 기본 설명 테이블
basicT = pd.read_csv('/disk0/sm/microbiome/16S/docker/info/basicTable_Y.txt', sep='\t')
# basicT = pd.read_csv('/work/info/basicTable_Y.txt', sep='\t')
## 미생물별 min, max 풍부도와 방향성 
abundanceT = pd.read_csv('/disk0/sm/microbiome/16S/docker/info/abundanceTable_Y.txt', sep='\t')
# abundanceT = pd.read_csv('/work/info/abundanceTable_Y.txt', sep='\t')
## 점수 별 percentage 부여
standardT = pd.read_csv('/disk0/sm/microbiome/16S/docker/info/standardSample_Y.txt', sep='\t')
# standardT = pd.read_csv('/work/info/standardSample_Y.txt', sep='\t')
standardT['info'] = pd.to_numeric(standardT['info'], errors='coerce')

# path = 'SM/TEST4'
path = sys.argv[1]

## 풍부도
# input_abundance = pd.read_csv('RESULT_ASV_TABLE.txt', sep='\t')
input_abundance = pd.read_csv('/disk0/sm/microbiome/16S/'+path+'/RESULT_ASV_TABLE.txt', sep='\t')
sampleL = input_abundance.columns[23:]

## 샘플별로 코드 돌려서 결과 얻기
results = {}

for sample in sampleL:
    results[sample] = P.pipeline_vagina(sample, input_abundance, basicT, abundanceT, standardT)

## mapping Table 불러와서 Key 값 API용으로 변환하기
mappingT = pd.read_csv('/disk0/sm/microbiome/16S/docker/info/mappingTable_Y.txt', sep='\t')[['key', 'value']]
# mappingT = pd.read_csv('/work/info/mappingTable_Y.txt', sep='\t')[['key', 'value']]
mapping_dict = mappingT.set_index('value')['key'].to_dict()

mapped_result = {sample: {'info': MD.map_keys(info, mapping_dict)} 
                 for sample, info in results.items()}

## 파일 저장하기
with open("Final_result.json", "w", encoding="utf-8") as f:
    json.dump(mapped_result, f, ensure_ascii=False, indent=2)


##########################
# import pandas as pd
# from pandas import json_normalize

# # 딕셔너리를 평탄화(flatten)
# df = json_normalize(mapped_result, sep='_').T.reset_index()
# df.columns = ['sample', 'info']

# df['col'] = df['sample'].str.split('_').str[1:].str.join('_')
# df['name'] = df['sample'].str.split('_').str[0]
# df.to_csv('/disk0/sm/microbiome/service/vagina/data/test_'+path+'.txt',sep='\t', index=False)

