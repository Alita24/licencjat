import pandas as pd
from collections import Counter

files = ['450_nep1.tsv']

interpro_counts = Counter()

for file in files:
    with open(file) as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith("#") and "InterPro" in line: 
                header_line = i
    
    if header_line is None:
        print("⚠️ No header line with 'InterPro' found, skipping this file.")
        continue

    header = lines[header_line].lstrip("#").strip().split("\t")

    df = pd.read_csv(file, sep='\t', comment='#',names=header) 
    
    for ids in df['InterPro'].dropna():
        for ipr in ids.split(','):
            interpro_counts[ipr.strip()] += 1

result = pd.DataFrame(interpro_counts.items(), columns=['InterPro_ID', 'Count'])
result = result.sort_values(by='Count', ascending=False)

print(result.head(20))

result.to_csv("interpro_count_nep1.tsv", sep='\t', index=False)