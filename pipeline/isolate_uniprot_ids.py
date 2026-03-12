import sys
from pathlib import Path

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))


def extract_uniprot_ids_from_file(filepath: Path) -> set[str]:
    """
    Extract unique UniProt accessions from the given file.
    Assumes fields with accessions contain patterns matching [A-NR-Z0-9]{6,10}
    """
    found = set()
    with filepath.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip().startswith("#"):
                continue
            cols = line.strip().split('\t')
            if len(cols) >= 3:
                found.add(cols[2])
    return found

def isolate_uniprot_ids(
    input_dir: Path,
    output_dir: Path,
    file_pattern: str = "*.tsv"
):
    """
    Extract unique UniProt IDs from all files in input_dir matching file_pattern and
    write them to output_dir as '<original_stem>_uniprot_ids.txt'.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    for file in input_dir.glob(file_pattern):
        uniprot_ids = extract_uniprot_ids_from_file(file)
        out_path = output_dir / f"{file.stem}_uniprot_ids.txt"
        with out_path.open("w", encoding="utf-8") as out_f:
            for uid in uniprot_ids:
                if uid.endswith('-F1'):
                    uid = uid[:-3]
                out_f.write(uid + "\n")
        print(f"{file.name}: wrote {len(uniprot_ids)} UniProt IDs to {out_path}")

# Example usage:
# from pathlib import Path
# isolate_uniprot_ids(
#    input_dir=Path("/path/to/alphaknot_results"),
#    output_dir=Path("/path/to/alphaknot_uniprot_ids")
# )
