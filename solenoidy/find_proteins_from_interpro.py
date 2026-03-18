import csv

from pathlib import Path

here = Path(__file__).resolve().parent

def find_proteins_with_interpro(tsv_path, id_col=0, name_col=1, ipr_col=3, delimiter="\t"):
    """
    Given a TSV file (csv), returns a list of protein IDs that contain the given InterPro code.
    The InterPro codes are assumed to be in the ipr_col column (0-indexed), separated by '|'.
    Args:
        tsv_path (str or Path): Path to the TSV file.
        interpro_code (str): The InterPro code to search for (e.g., "IPR009091").
        id_col (int): The 0-based column index for protein ID (default: 0).
        name_col (int): The 0-based column index for protein name (default: 1).
        ipr_col (int): The 0-based column index for InterPro codes (default: 3).
        delimiter (str): Delimiter for the TSV file (default: '\t').
    """
    # Dictionary mapping InterPro code -> list of (uniprot_id, name) pairs
    ipr_to_proteins = {}
    with open(tsv_path, "r", encoding="utf-8") as tsvfile:
        reader = csv.reader(tsvfile, delimiter=delimiter)
        for row in reader:
            # Skip comment lines
            if not row or row[0].startswith('#'):
                continue
            if len(row) > max(id_col, name_col, ipr_col):
                uniprot_id = row[id_col].strip()
                name = row[name_col].strip() if name_col < len(row) else ""
                ipr_field = row[ipr_col]
                codes = {code.strip() for code in ipr_field.split("|") if code.strip()}
                print(f'for {uniprot_id}, {name}', end=': ')
                for ipr_code in codes:
                    if not ipr_code:
                        continue
                    print(ipr_code, end=' ')
                    ipr_to_proteins.setdefault(ipr_code, []).append((uniprot_id, name))
                print()
    
    # Go two directories up and create the folder "prot_with_intpr"
    proteins_dir = Path(tsv_path).resolve().parents[1] / "prot_with_intpr"
    proteins_dir.mkdir(parents=True, exist_ok=True)

    for ipr_code, proteins in ipr_to_proteins.items():
        print(f'making file for {ipr_code}')
        if not ipr_code.startswith('IPR'):
            continue
        output_csv = proteins_dir / f"proteins_with_{ipr_code}.csv"
        with open(output_csv, "w", newline='', encoding="utf-8") as outcsv:
            writer = csv.writer(outcsv)
            writer.writerow(["UniProtID", "Name"])  # Header row
            for uniprot_id, name in proteins:
                writer.writerow([uniprot_id, name])

if __name__ == "__main__":
    find_proteins_with_interpro(
        f"{here}/isolated_partners/solenoid_interactor_IPR.csv",
        delimiter=','
    )
