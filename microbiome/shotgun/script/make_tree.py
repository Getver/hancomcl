import sys

if len(sys.argv) != 2:
    print("Usage: python3 make_tree.py <metaphlan_profile_file>")
    sys.exit(1)

input_file = '/home/data/metaphlan/' + sys.argv[1] + '_profile.txt'
output_file = input_file.rsplit(".", 1)[0] + ".tree"

with open(input_file) as f, open(output_file, "w") as out:
    taxa_set = set()
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        taxon = line.split()[0]  # e.g. k__Bacteria|p__Firmicutes|c__Clostridia
        ranks = taxon.split("|")
        path_str = ".".join([r.split("__")[1] for r in ranks if "__" in r])
        if path_str not in taxa_set:
            out.write(path_str + "\n")
            taxa_set.add(path_str)

print(f"Tree file generated: {output_file}")
