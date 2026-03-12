import sys
from pathlib import Path

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

from intAct.api import fetch_intact_interactions

def get_uniprot_ids_from_file(file_path):
    """Read UniProt IDs from a file (one per line)."""
    with open(file_path, 'r') as f:
        ids = [line.strip() for line in f if line.strip()]
    return ids

def process_all_uniprot_id_files(ids_dir, out_dir):
    """
    For every file in ids_dir, treat the filename (without '_uniprot_ids.txt') as a gene,
    fetch IntAct interactions for each UniProt ID listed, and save results in out_dir.
    """
    ids_dir = Path(ids_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    for file in ids_dir.glob("*_uniprot_ids.txt"):
        gene_name = file.stem.replace("_uniprot_ids", "")
        uniprot_ids = get_uniprot_ids_from_file(file)
        for uniprot_id in uniprot_ids:
            print(f"Fetching IntAct for {uniprot_id} (gene: {gene_name})")
            out_file = out_dir / f"intact_{gene_name}_{uniprot_id}.tsv"
            try:
                fetch_intact_interactions(uniprot_id, out_tsv_path=out_file)
            except Exception as e:
                print(f"Failed to fetch for {uniprot_id}: {e}")

def main():
    """
    Fetches IntAct interaction data for all UniProt IDs listed in pipeline/alphaknot_uniprot_ids/*.txt,
    and stores the results as TSV files in 'intact_uniprot_results'.
    """
    ids_dir = Path(__file__).resolve().parent / "alphaknot_uniprot_ids"
    out_dir = Path(__file__).resolve().parent / "intact_uniprot_results"
    process_all_uniprot_id_files(ids_dir, out_dir)

if __name__ == "__main__":
    main()