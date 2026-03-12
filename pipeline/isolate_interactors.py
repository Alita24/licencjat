from pathlib import Path
import pandas as pd

def isolate_partner_uniprot_ids(group_name, intact_dir):
    """
    For a given group, collect ALL UniProt IDs that interacted with any protein from the group
    (i.e., the partners from the IntAct tsvs: intact_{group}_{uniprot}.tsv).
    Returns a set of partner UniProt IDs (not from the original group).
    """
    intact_dir = Path(intact_dir)
    partner_ids = set()

    for tsv_file in intact_dir.glob(f"intact_{group_name}_*.tsv"):
        try:
            df = pd.read_csv(tsv_file, sep="\t")
            id_cols = [c for c in df.columns if 'uniprot' in c.lower()]
            if len(id_cols) >= 2:
                uniprot_a = id_cols[0]
                uniprot_b = id_cols[1]
            elif 'ID interactor A' in df.columns and 'ID interactor B' in df.columns:
                uniprot_a = 'ID interactor A'
                uniprot_b = 'ID interactor B'
            else:
                continue
            # Here, partner is usually interactor B
            partner_ids.update(str(row[uniprot_b]) for _, row in df.iterrows())
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

    partners = isolate_partner_uniprot_ids(group_name, intact_dir)
    with open(out_path, "w") as f:
        for pid in sorted(partners):
            f.write(pid + "\n")
    print(f"Partner UniProt IDs for group '{group_name}' written to {out_path}")
    return partners

