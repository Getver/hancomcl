# import pandas as pd
# import json

def make_dict(sample, basicT, filterG, filterS):
    Lac_list = [
        'Lactobacillus',
        'Lactobacillus crispatus',
        'Lactobacillus gasseri',
        'Lactobacillus iners',
        'Lactobacillus jensenii',
        'OLac'
    ]
    # dict로 반복 생성
    Lac_dict = {}
    for sp in Lac_list:
        if sp == 'Lactobacillus':
            Lac_dict[sp] = next(iter(filterG.loc[filterG['name'] == sp, sample].values), 0.0)
        elif sp == 'OLac':
            # 주요 4종 제외한 나머지 Lactobacillus species 합계
            major_species = [
                'Lactobacillus crispatus',
                'Lactobacillus gasseri',
                'Lactobacillus iners',
                'Lactobacillus jensenii'
            ]
            other_sum = filterS.loc[filterS['Species'].isin(major_species), sample].sum()
            Lac_dict[sp] = Lac_dict['Lactobacillus'] - other_sum
        else:
            Lac_dict[sp] = next(iter(filterS.loc[filterS['Species'] == sp, sample].values), 0.0)
    ##
    BVM = basicT.loc[basicT['구분'] == '혐기성', 'name_txt'].to_list()
    BV = filterG.loc[filterG['name'].isin(BVM), ['name', sample]]
    BV_dict = dict(zip(BV['name'], BV[sample]))
    BVsum = BV[sample].sum()
    ##
    AVM = basicT.loc[basicT['구분'] == '호기성', 'name_txt'].to_list()
    AV = filterG.loc[filterG['name'].isin(AVM), ['name', sample]]
    AV_dict = dict(zip(AV['name'], AV[sample]))
    AVsum = AV[sample].sum()
    ##
    HVM = basicT.loc[basicT['구분'] == '유해균', 'name_txt'].to_list()
    HV = filterG.loc[filterG['name'].isin(HVM), ['name', sample]]
    HV_dict = dict(zip(HV['name'], HV[sample]))
    ##
    BifidoV = filterG.loc[filterG['name'] == 'Bifidobacterium', ['name', sample]]
    Bifido_dict = dict(zip(BifidoV['name'], BifidoV[sample]))
    M_dict = Lac_dict | AV_dict | BV_dict | HV_dict | Bifido_dict
    M_dict['AVsum'] = AVsum
    M_dict['BVsum'] = BVsum
    return M_dict



def map_keys(d, mapping):
    if isinstance(d, dict):
        return {mapping.get(k, k): map_keys(v, mapping) for k, v in d.items()}
    elif isinstance(d, list):
        return [map_keys(i, mapping) for i in d]
    else:
        return d
