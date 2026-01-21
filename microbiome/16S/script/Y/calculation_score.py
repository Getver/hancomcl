import pandas as pd
import numpy as np
import math

def score_minmax(sample, filterG, basicT, abundanceT) :
    T = pd.merge(basicT, abundanceT, on=['name_txt', 'rank'], how='left')
    G = T.loc[T['rank'] == 'genus',['name_txt', 'min', 'max', 'rank']].dropna()
    score = pd.merge(G, filterG, left_on='name_txt', right_on='name', how='left').dropna()
    ## min:10, max:90
    score['ten'] = (score['max']-score['min'])/10
    score['SRR_scaled'] = ((score[sample]*100 - (score['min']-score['ten'])) / ((score['max']+score['ten']) - (score['min']-score['ten']))).clip(0, 1)
    return score

def disease_score(score, CST_type, abundanceT) :
    disease_score_dict = {}
    for i in abundanceT.columns[5:-1]:  # 질환 컬럼 범위
        tmpT = abundanceT.loc[abundanceT[i] == '+', ['name_txt', i]]
        score_direct = pd.merge(tmpT, score, on='name_txt', how='left')
        score_direct['score'] = score_direct.apply(lambda row: row['SRR_scaled'] if row[i] == '-' else 1 - row['SRR_scaled'], axis=1)
        mean_score = score_direct['score'].mean()
        # score['score'] = score_direct.apply(lambda row: row['SRR_scaled'] if row[i] == '-' else 1 - row['SRR_scaled'], axis=1)
        # mean_score = score['score'].mean()
        weight = {'CST1' : 1, 
                'CST2' : 0.9, 
                'CST3' : 0.8, 
                'CST4-BF' : 0.8,
                'CST4-AV' : 0.7, 
                'CST4-BV' : 0.7, 
                'CST5' : 0.9}
        # 딕셔너리에 저장
        disease_score_dict[i] = mean_score * weight.get(CST_type, 1) * 100    ## min-max에 가중치 부여 하는 법
        # disease_score_dict[i] = mean_score * 100
        # disease_score_dict[i] = np.sqrt(mean_score) * 100   # 제곱근
        # disease_score_dict[i] = mean_score ** 2 * 100    # 제곱
        disease_score_dict = {k: round(v) for k, v in disease_score_dict.items()}
    return disease_score_dict


def antibiotic_type(sample, abundanceT, filterG) : 
    anti = abundanceT[['name_txt', '항생제 내성']].dropna()
    antiAB = pd.merge(anti, filterG, left_on='name_txt', right_on='name', how='inner')
    antiAB_expanded = (
    antiAB.assign(antibiotic=antiAB['항생제 내성'].str.split(';'))
            .explode('antibiotic')
            .reset_index(drop=True)
    )
    antibiotic_sum = (
        antiAB_expanded
        .groupby('antibiotic', as_index=False)[sample]
        .sum()
        .sort_values(sample, ascending=False)
    )
    antibiotic_dict = antibiotic_sum.set_index('antibiotic')[sample].to_dict()
    antibiotic_dict = {k: round(v, 4) for k, v in antibiotic_dict.items()}
    return antibiotic_dict


# def calc_percentile_from_reference(standardT, disease, new_score):
#     ref = (standardT.loc[standardT['disease_value'] == disease, 'info'].astype(float).dropna().values)
#     if len(ref) == 0:
#         return None
#     percentile = (1-(ref <= float(new_score)).mean()) * 100
#     return round(percentile, 1)

def calc_percentile_from_reference(standardT, disease, new_score):
    ref = (standardT.loc[standardT['disease_value'] == disease, 'info'].astype(float).dropna())
    if len(ref) == 0:
        return None
    rank = (ref < float(new_score)).sum() + 0.5 * (ref == float(new_score)).sum()
    percentile = 100 * (1 - rank / len(ref))
    return round(percentile, 1)


def add_disease_percentile(disease_score_dict, standardT):
    result = {}
    for disease, score in disease_score_dict.items():
        percentile = calc_percentile_from_reference(
            standardT,
            disease,
            score
        )
        result[disease] = {
            'score': round(float(score), 2),
            'percentile': percentile
        }
    return result

