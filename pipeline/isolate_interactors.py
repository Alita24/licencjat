from pathlib import Path
import sys
import pandas as pd

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

def isolate_all_uniprot_ids(group_name, intact_dir):
    """
    Zwraca wszytskie UniProt ID dla podanej grupy białek,
    przeszukując pliki wynikowe IntAct dla tej grupy w zadanym katalogu.
    """
    intact_dir = Path(intact_dir)
    partner_ids = set()

    for tsv_file in intact_dir.glob(f"intact_{group_name}_*.tsv"):
        try:
            df = pd.read_csv(tsv_file, sep="\t")
            # Both columns may contain partners; collect from both interactor columns
            for _, row in df.iterrows():
                partner_ids.add(str(row[0]))
                partner_ids.add(str(row[1]))
        except Exception as e:
            print(f"Failed to read {tsv_file}: {e}")
    return partner_ids

def write_partner_ids_for_group(group_name, intact_dir=None, out_path=None):
    """
    Collects UniProt IDs of all proteins that interact with the given protein group, writes to file.
    Parameters:
        group_name (str): The name of the protein group.
        intact_dir (str or Path): Directory containing intact result tsvs.
        out_path (str or Path): Output filename for partner UniProt IDs.
    Returns:
        set: Set of collected partner UniProt IDs.
    """
    print(f"Starting collection of partner UniProt IDs for group: {group_name}")

    if intact_dir is None:
        here = Path(__file__).resolve().parent
        intact_dir = here / "intact_uniprot_results"
    else:
        intact_dir = Path(intact_dir)

    if out_path is None:
        here = Path(__file__).resolve().parent
        out_path = here / f"{group_name}_partners.txt"
    else:
        out_path = Path(out_path)

    print(f"Collecting all UniProt IDs for group '{group_name}' from files in: {intact_dir}")
    all_uniprot_ids = isolate_all_uniprot_ids(group_name, intact_dir)
    # Filter out UniProt IDs already present in uniprot_ids/{group_name}_uniprot_ids_*.tsv
    here = Path(__file__).resolve().parent
    uniprot_ids_dir = here / "uniprot_ids"
    print(f"Searching for already known group UniProt IDs in directory: {uniprot_ids_dir}")
    group_uniprot_ids = set()
    for f in uniprot_ids_dir.glob(f"{group_name}_uniprot_ids_*.txt"):
        with open(f, "r") as fin:
            group_uniprot_ids.update(line.strip() for line in fin if line.strip())
    partners = all_uniprot_ids - group_uniprot_ids
    print(f"Total partner UniProt IDs (excluding already present): {len(partners)}")
    with open(out_path, "w") as f:
        for pid in sorted(partners):
            f.write(pid + "\n")
    print(f"Partner UniProt IDs for group '{group_name}' written to {out_path}")
    return partners

if __name__ == "__main__":
    # Example usage: provide the group name and directory containing intact results
    # intact_gene_results = Path("pipeline/intact_uniprot_results")  # Change this path as needed
    print('hi')
    # write_partner_ids_for_group('trmd', Path(intact_gene_results))


