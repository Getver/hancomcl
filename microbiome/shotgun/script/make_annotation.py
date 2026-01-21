import sys
import pandas as pd
from io import StringIO

if len(sys.argv) != 2:
    print("Usage: python3 make_annotation.py <sample_name>")
    sys.exit(1)

profile_file = '/home/data/metaphlan/' + sys.argv[1] + '_profile.txt'
output_file = profile_file.rsplit(".", 1)[0] + ".annot"

lines = []
header_line = None

with open(profile_file) as f:
    for line in f:
        if line.startswith('#clade_name'):
            header_line = line.lstrip('#')  # 헤더 줄, '#' 제거
        elif line.startswith('#') or line.strip() == '':
            continue
        else:
            lines.append(line)

if header_line is None:
    print("Error: header line with '#clade_name' not found.")
    sys.exit(1)

# DataFrame 생성
df = pd.read_csv(StringIO(header_line + ''.join(lines)), sep='\t')

# annotation 파일 작성
with open(output_file, 'w', encoding='utf-8') as out:
    out.write('clade_name\tannotation_type\tannotation_value\n')
    for clade in df['clade_name']:
        # 'k__Bacteria|p__Actinobacteria' -> 'Bacteria.Actinobacteria'
        clade_clean = '.'.join([level.split('__')[1] if '__' in level else level for level in clade.split('|')])
        color = '#1f78b4'  # 임의 색상 예시 (파란색)
        out.write(f"{clade_clean}\tcolor\t{color}\n")

