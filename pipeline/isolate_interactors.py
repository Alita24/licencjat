from pathlib import Path
import sys
import pandas as pd
from targets import GENE_TO_IPRS

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
            df = pd.read_csv(tsv_file, sep="\t", skiprows=1)
            
            # Both columns may contain partners; collect from both interactor columns
            for _, row in df.iterrows():
                id1, id2 = row.iloc[0], row.iloc[1]
                if str(id1).startswith('uniprotkb'):
                    id1 = str(id1).split(':', 1)[1]
                if str(id2).startswith('uniprotkb'):
                    id2 = str(id2).split(':', 1)[1]
                partner_ids.add(id1)
                partner_ids.add(id2)
        except Exception as e:
            print(f"Failed to read {tsv_file}: {e}")
    return partner_ids

def write_partner_ids_for_group(group_name, intact_dir, out_path):
    """
    Collects UniProt IDs of all proteins that interact with the given protein group, writes to file.
    Parameters:
        group_name (str): The name of the protein group.
        intact_dir (str or Path): Directory containing intact result tsvs, default.
        out_path (str or Path): Output filename for partner UniProt IDs.
    Returns:
        set: Set of collected partner UniProt IDs.
    """
    print(f"Starting collection of partner UniProt IDs for group: {group_name}")


    intact_dir = Path(intact_dir)
    out_path = Path(out_path)
    out_path.mkdir(parents=True, exist_ok=True)

    out_file = out_path / f"{group_name}_partners.txt"

    print(f"Collecting all UniProt IDs for group '{group_name}' from files in: {intact_dir}")
    all_uniprot_ids = isolate_all_uniprot_ids(group_name, intact_dir)
   
    # Filter out UniProt IDs already present in uniprot_ids/{group_name}_uniprot_ids_*.tsv
    
    uniprot_ids_dir = parent_dir / 'pipeline' / "uniprot_ids"
    print(f"Searching for already known group UniProt IDs in directory: {uniprot_ids_dir}")
    
    group_uniprot_ids = set()
    for f in uniprot_ids_dir.glob(f"{group_name}_uniprot_ids_*.txt"):
        with open(f, "r") as fin:
            group_uniprot_ids.update(line.strip() for line in fin if line.strip())
    partners = all_uniprot_ids - group_uniprot_ids
# 
    print(f"Total partner UniProt IDs (excluding already present): {len(partners)}")
    with open(out_file, "w") as f:
        for pid in sorted(partners):
            f.write(pid + "\n")
    print(f"Partner UniProt IDs for group '{group_name}' written to {out_file}")
    return partners

def isolate_partners(intact_dir, output_dir):
    for fam in GENE_TO_IPRS.keys():
        write_partner_ids_for_group(fam, intact_dir=intact_dir, out_path=output_dir)

if __name__ == "__main__":
    # Example usage: provide the group name and directory containing intact results
    intact_gene_results = Path("pipeline/intact_uniprot_results")  # Change this path as needed
    # print('hi')
    write_partner_ids_for_group('trmD', Path(intact_gene_results), Path('pipeline/isolated_partners'))


