# import pandas as pd
import transfer_rank as TR
import preprocessing_data as PD
import make_dictionary as MD
import return_result as RR
import calculation_score as CS
# import antibiotic_resistant as AR

def pipeline_vagina(sample, input_abundance, basicT, abundanceT, standardT) :
    ## 데이터 파싱 : genus/species 로 풍부도 계산하기
    genus_ab, species_ab = PD.prepare_data(sample, input_abundance) # Bifidobacterium
    ## target 가져오기
    filterG, filterS = PD.target_select(genus_ab, species_ab, abundanceT)   # Bifidobacterium
    ## @@ 미생물 별 풍부도 딕셔너리, 카테고리별 풍부도(혐기성, 호기성)
    result_dict = MD.make_dict(sample, basicT, filterG, filterS)
    ## @@ CST 분류 결과
    CST_type = RR.return_CST(result_dict)
    ## 유익균, 유해균, 기타 풍부도 계산
    UU_dict = TR.U_abundance(sample, genus_ab, basicT)
    ## @@ rate를 기준으로 점수 매김
    # scoreT, total_score = CS.score_minmax(sample, filterG, basicT, abundanceT)
    scoreT = CS.score_minmax(sample, filterG, basicT, abundanceT)
    ## @@ disease 별 점수 매겼음
    # disease_score_dict = CS.disease_score(scoreT, CST_type, diseaseT)
    raw_disease_score = CS.disease_score(scoreT, CST_type, abundanceT)
    disease_score_dict = CS.add_disease_percentile(raw_disease_score, standardT)
    ## @@ 항생제 풍부도 반환했음
    antibiotic_dict = CS.antibiotic_type(sample, abundanceT, filterG)
    ## 결과 딕셔너리를 깔끔하게 정리
    result_dict_clean = {k: round(float(v), 4) for k, v in result_dict.items()}
    ## 내가 원하는 형식으로 변환
    final_result = {
        'CST_type': CST_type,
        'total_rate': UU_dict,
        'disease_score': disease_score_dict,
        'antibiotic_abundance':antibiotic_dict, 
        'abundance': result_dict_clean
    }
    return final_result
