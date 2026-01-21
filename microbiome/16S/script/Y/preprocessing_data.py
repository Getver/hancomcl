# import pandas as pd
import transfer_rank as TR

def prepare_data(sample,input_abundance) :
    # input_abundance['name_txt'].replace('uncultured ', '', regex=True, inplace=True)
    # input_abundance['name_txt'].replace('Candidatus ', '', regex=True, inplace=True)
    abundance_value = input_abundance.loc[:, sample]
    genus_ab = TR.genus_abundance(input_abundance, abundance_value)
    species_ab = TR.species_abundance(input_abundance, abundance_value)
    return genus_ab, species_ab

def target_select(genus_ab, species_ab, abundanceT) :
    targetM = abundanceT['name_txt'].to_list()
    filterG = genus_ab[genus_ab['name'].isin(targetM)]
    filterS = species_ab[species_ab['Species'].isin(targetM)]
    return filterG, filterS

