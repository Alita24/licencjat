import sys
from pathlib import Path
import pandas as pd
from shutil import rmtree

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

from intAct.api import fetch_intact_interactions

def intact_uniprot_pipeline(ids_dir, out_dir, file_pattern="_uniprot_ids"):
    """
    For every file in ids_dir, treat the filename (without '_uniprot_ids.txt') as a gene,
    fetch IntAct interactions for each UniProt ID listed, and save results in out_dir.
    """
    ids_dir = Path(ids_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    
    for file in ids_dir.glob('*'+file_pattern+".txt"):
        fam_name = file.stem.replace(file_pattern, "")

        combined_file = out_dir / f"intact_{fam_name}_UNIPROT.tsv"

        if combined_file.exists():
            print(f"Skipping {fam_name}: combined file already exists.")
            continue

        with open(file, 'r') as f:
            uniprot_ids = [line.strip() for line in f if line.strip()]

        fam_subdir = out_dir / f"intact_tmp_{fam_name}"
        fam_subdir.mkdir(exist_ok=True)

        for uniprot_id in uniprot_ids:
            out_file = fam_subdir / f"intact_{fam_name}_{uniprot_id}.tsv"

            if out_file.exists():
                print(f"Skipping fetch for {uniprot_id} (name: {fam_name}): outfile already exists.")
                continue

            print(f"Fetching IntAct for {uniprot_id} (name: {fam_name})")
            try:
                fetch_intact_interactions(uniprot_id, out_tsv_path=out_file)
            except Exception as e:
                print(f"Failed to fetch for {uniprot_id}: {e}")
                out_file.touch(exist_ok=True)
                pass

            if not out_file.exists():
                out_file.touch(exist_ok=True)

        all_results = []
        for file in fam_subdir.glob(f"intact_{fam_name}_*.tsv"):
            try:
                df = pd.read_csv(file, sep="\t")
            except pd.errors.EmptyDataError:
                df = pd.DataFrame()
            all_results.append(df)
        if all_results:
            pd.concat(all_results).to_csv(combined_file, sep="\t", index=False)
        else:
            combined_file.touch(exist_ok=True)
 


def main():
    ids_dir = Path(__file__).resolve().parent  / "DELETE"
    out_dir = Path(__file__).resolve().parent / "intact_uniprot_results"
    intact_uniprot_pipeline(ids_dir, out_dir)
    print(f"IntAct interactions for all UniProt IDs saved to {out_dir}")

if __name__ == "__main__":
    main()