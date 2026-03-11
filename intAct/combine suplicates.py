import sys
import csv
from collections import defaultdict

def combine_duplicates(input_file, output_file):
    """
    Reads a MITAB 2.5 format file (or similar TSV with ID(s) interactor A/B as columns 0 and 1),
    and combines rows where (idA, idB) are the same (or swapped). Merges all fields with '|'.
    Outputs the combined table.
    """
    pair_to_rows = defaultdict(list)
    with open(input_file, encoding='utf-8') as fin:
        reader = csv.reader(fin, delimiter='\t')
        header = next(reader)
        if header[0].lower().startswith("# id(s) interactor a"):
            # MITAB comment header (keep it)
            fieldnames = header
        else:
            # No MITAB comment header, treat as data
            fieldnames = None
            fin.seek(0)
            reader = csv.reader(fin, delimiter='\t')

        rows = list(reader)

    # Combine rows by interactor pairs (unordered, i.e. (A,B) == (B,A))
    for row in rows:
        if not row or len(row) < 2:
            continue
        pair = tuple(sorted([row[0], row[1]]))
        pair_to_rows[pair].append(row)

    combined_rows = []
    for pair, row_group in pair_to_rows.items():
        # To handle different number of columns, align and fill missing with ""
        ncols = max(len(row) for row in row_group)
        merged = []
        for i in range(ncols):
            merged_vals = set()
            for row in row_group:
                if i < len(row):
                    merged_vals.update(x for x in row[i].split('|') if x)
            merged.append('|'.join(sorted(merged_vals)))
        combined_rows.append(merged)

    # Write the output
    with open(output_file, 'w', encoding='utf-8', newline='') as fout:
        writer = csv.writer(fout, delimiter='\t')
        if fieldnames:
            writer.writerow(fieldnames)
        for row in combined_rows:
            writer.writerow(row)

if __name__ == "__main__":
    name = 'trm5'
    combine_duplicates(f"intAct_{name}.tsv", f"intAct_{name}_combined.tsv")
