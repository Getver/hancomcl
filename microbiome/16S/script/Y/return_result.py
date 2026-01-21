# import pandas as pd
import transfer_rank as TR

# def return_CST(Lac_dict, AVsum, BVsum) :
#     Lac = TR.get_value(Lac_dict['Lactobacillus'])
#     L1 = TR.get_value(Lac_dict['Lactobacillus crispatus'])
#     L2 = TR.get_value(Lac_dict['Lactobacillus gasseri'])
#     L3 = TR.get_value(Lac_dict['Lactobacillus iners'])
#     L5 = TR.get_value(Lac_dict['Lactobacillus jensenii'])
#     ## dominant = 
#     CST = ''
#     if L1 >= 0.5: CST = 'CST1'
#     elif L2 >= 0.5: CST = 'CST2'
#     elif L3 >= 0.5: CST = 'CST3'
#     elif L5 >= 0.5: CST = 'CST5'
#     else:
#         if Lac >= 0.5: CST = 'CST3'
#         elif AVsum >= BVsum:
#             CST = 'CST4-AV'
#         else:
#             CST = 'CST4-BV'
#     ##
#     return CST

# def return_CST(Lac_dict) :
#     values = {
#         "CST1": TR.get_value(Lac_dict['Lactobacillus crispatus']),
#         "CST2": TR.get_value(Lac_dict['Lactobacillus gasseri']),
#         "CST3": TR.get_value(Lac_dict['Lactobacillus iners']),
#         "CST4-BF": TR.get_value(Lac_dict['Bifidobacterium']),
#         "CST5": TR.get_value(Lac_dict['Lactobacillus jensenii']),
#         "CST4-AV": Lac_dict['AVsum'],
#         "CST4-BV": Lac_dict['BVsum'],
#     }
#     CST = max(values, key=values.get)
#     return CST



def return_CST(Lac_dict):
    values = {
        "CST1": TR.get_value(Lac_dict.get('Lactobacillus crispatus', 0)),
        "CST2": TR.get_value(Lac_dict.get('Lactobacillus gasseri', 0)),
        "CST3": TR.get_value(Lac_dict.get('Lactobacillus iners', 0)),
        "CST5": TR.get_value(Lac_dict.get('Lactobacillus jensenii', 0)),

        # CST4 subtypes
        "CST4-BF": TR.get_value(Lac_dict.get('Bifidobacterium', 0)),
        "CST4-AV": TR.get_value(Lac_dict.get('AVsum', 0)),
        "CST4-BV": TR.get_value(Lac_dict.get('BVsum', 0)),
    }

    # 전부 0이면 fallback
    if max(values.values()) == 0:
        return "CST-Unknown"

    CST = max(values, key=values.get)
    return CST