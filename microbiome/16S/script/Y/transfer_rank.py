import pandas as pd

# input_abundance = pd.read_csv('/disk0/sm/microbiome/16S/1013_park/Run_0/abundance/abundance.tsv', header=1, sep='\t')

def genus_abundance(input_abundance, abundance_value) :
    # genus 이름만 따기
    input_abundance['name'] = input_abundance['Genus'].str.split(' ').str[0]
    # genus 이름 + 풍부도 테이블
    df = pd.concat([input_abundance[['name']], pd.DataFrame(abundance_value)], axis=1)
    # genus 이름 같으면 풍부도 합치기
    collapsed = df.groupby('name', as_index=False)[df.columns[1]].sum()
    return collapsed


def species_abundance(input_abundance, abundance_value) :
    # species 이름 + 풍부도 테이블
    df = pd.concat([input_abundance[['Species']], pd.DataFrame(abundance_value)], axis=1)
    # species 이름 같으면 풍부도 합치기
    collapsed = df.groupby('Species', as_index=False)[df.columns[1]].sum()
    return collapsed


def get_value(x):
    """NumPy array, list, or float 입력 처리하여 순수 Python float 반환"""
    if isinstance(x, (float, int)):  # 이미 숫자면 그대로 반환
        return float(x)
    try:
        # NumPy array 또는 list 처리
        if len(x) == 0:
            return 0.0
        return float(x[0])
    except TypeError:
        # 그 외 처리 불가 타입
        return 0.0


def U_abundance(sample, tmp_filterG, basicT):
    # tmp_filterG와 basicT 병합하여 '구분' 컬럼 추가
    filterG = pd.merge(
        tmp_filterG,
        basicT[['name_txt', '구분']],
        left_on='name',
        right_on='name_txt',
        how='left'
    )
    # sample이 문자열이면 단일 컬럼, 리스트면 다중 컬럼 처리
    abundance_value = filterG[[sample]].copy()
    # '구분' 컬럼 추가
    abundance_value['구분'] = filterG['구분'].values
    abundance_value['구분'] = abundance_value['구분'].replace({'상재균': '기타', '혐기성': '유해균', '호기성': '유해균', '병원균': '유해균'})
    collapsed = abundance_value.groupby('구분', as_index=False).sum()
    UU_dict = collapsed.set_index('구분')[sample].to_dict()
    total_sum = sum(UU_dict.values())
    if abs(total_sum - 1) > 1e-6:
        if '기타' in UU_dict:
            UU_dict['기타'] = max(0, 1 - (UU_dict.get('유익균', 0) + UU_dict.get('유해균', 0)))
        else:
            # 기타가 없으면 새로 추가
            UU_dict['기타'] = max(0, 1 - total_sum)
    UU_dict = {k: round(v, 2) for k, v in UU_dict.items()}
    return UU_dict
