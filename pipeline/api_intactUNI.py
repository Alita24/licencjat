import sys
from pathlib import Path
import pandas as pd

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

from intAct.api import fetch_intact_interactions

def intact_uniprot_pipeline(ids_dir, out_dir):
    """
    For every file in ids_dir, treat the filename (without '_uniprot_ids_COMBINED.txt') as a gene,
    fetch IntAct interactions for each UniProt ID listed, and save results in out_dir.
    """
    ids_dir = Path(ids_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    for file in ids_dir.glob("*_uniprot_ids_COMBINED.txt"):
        gene_name = file.stem.replace("_uniprot_ids_COMBINED", "")

        with open(file, 'r') as f:
            uniprot_ids = [line.strip() for line in f if line.strip()]
        for uniprot_id in uniprot_ids:
            print(f"Fetching IntAct for {uniprot_id} (gene: {gene_name})")
            out_file = out_dir / f"intact_{gene_name}_{uniprot_id}.tsv"
            try:
                fetch_intact_interactions(uniprot_id, out_tsv_path=out_file)
            except Exception as e:
                print(f"Failed to fetch for {uniprot_id}: {e}")
                pass

        all_results = []
        files_to_remove = []
        for file in out_dir.glob(f"intact_{gene_name}_*.tsv"):
            df = pd.read_csv(file, sep="\t")
            all_results.append(df)
            files_to_remove.append(file)
        combined_file = out_dir / f"intact_{gene_name}_UNIPROT.tsv"
        pd.concat(all_results).to_csv(combined_file, sep="\t", index=False)
  
        for file in files_to_remove:
            if file != combined_file:
                file.unlink()


def main():
    ids_dir = Path(__file__).resolve().parent  / "DELETE"
    out_dir = Path(__file__).resolve().parent / "intact_uniprot_results"
    intact_uniprot_pipeline(ids_dir, out_dir)
    print(f"IntAct interactions for all UniProt IDs saved to {out_dir}")

if __name__ == "__main__":
    main()