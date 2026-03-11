# Script to check if the UniProt ID from the data is in ../intAct/uniprot_ids.txt and write those rows to an output file

def load_uniprot_ids(filepath):
    """Load UniProt IDs from the given text file, skip "-" or empty lines."""
    ids = set()
    with open(filepath, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line == "-":
                continue
            ids.add(line)
    return ids

def filter_rows_with_uniprot(input_tsv, uniprot_ids, output_tsv):
    """Write rows whose UniProt ID (3rd column, tab-separated) is in uniprot_ids set."""
    with open(input_tsv, encoding="utf-8") as f, open(output_tsv, "w", encoding="utf-8") as outf:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            uniprot_id = parts[2]
            if uniprot_id.endswith("-F1"):
                uniprot_id = uniprot_id[:-3]
            # They're sometimes in form P12345-2, should match exactly
            if uniprot_id in uniprot_ids:
                outf.write(line)

if __name__ == "__main__":
    # Example usage, adapt file names as needed:
    uniprot_set = load_uniprot_ids("../intAct/uniprot_ids.txt")
    # Assuming your input file is named 'input.tsv', and saving matched rows to 'matched_in_intAct.tsv'
    filter_rows_with_uniprot("38tys-_trmD.tsv", uniprot_set, "matched_in_intAct.tsv")
