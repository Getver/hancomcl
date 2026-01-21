from Bio import Phylo
import pandas as pd

tree = Phylo.read("tree/tree.nwk", "newick")
tax = pd.read_csv("taxonomy/taxonomy.tsv", sep="\t")

# g__ 부분만 추출
tax["Genus"] = tax["Taxon"].str.extract(r"g__([^;]*)")

# 매핑
id_to_genus = dict(zip(tax["Feature ID"], tax["Genus"]))

for clade in tree.get_terminals():
    if clade.name in id_to_genus and pd.notna(id_to_genus[clade.name]):
        clade.name = id_to_genus[clade.name]

Phylo.write(tree, "tree/tree_genus_labeled.nwk", "newick")







