import sys
from pathlib import Path

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))


def isolate_uniprot_ids(
    input_dir: Path,
    output_dir: Path,
    file_pattern: str = "*.tsv",
    interpro: bool = False,
):
    """
    wyciaga unikalne identyfikatory UniProt z plików w input_dir i zapisuje je do output_dir
    Args:
        input_dir: ścieżka do katalogu, w którym są pliki z danymi
        output_dir: ścieżka do katalogu, w którym będą zapisywane wyniki
        file_pattern: pattern plików do przetworzenia
        interpro: czy przetwarzaj pliki z InterPro
    Returns:
        None
    """
    print(f"rozpoczynam isolating uniprot ids dla interpro {'yes' if interpro else 'no'}...")
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for file in input_dir.glob(file_pattern):
        # .glob iterator of all the files that match the file_pattern
        found_ids = set()
            # set of unique UniProt IDs
        with file.open("r", encoding="utf-8") as f:
            for line in f:
                if line.strip().startswith("#") or line.strip().startswith("protein_accession"):
                    continue
                cols = line.strip().split("\t")
                if len(cols) >= 3:
                    if interpro:
                        found_ids.add(cols[0])
                    else:
                        found_ids.add(cols[2])

        if interpro:
            out_path = output_dir / f"{file.stem}_uniprot_ids_IN.txt"
        else:
            out_path = output_dir / f"{file.stem}_uniprot_ids_GE.txt"
        with out_path.open("w", encoding="utf-8") as out_f:
            for uid in found_ids:
                if uid.endswith("-F1"):
                    uid = uid[:-3]
                out_f.write(uid + "\n")
        print(f"{file.name}: wrote {len(found_ids)} UniProt IDs to {out_path}")
    
    print("isolate_uniprot_ids completed.")

def combine_list(ids_dir):
    """
    polacz oba listy z uniprot ids
    """
    ids_dir = Path(ids_dir)
    in_files = {f.stem.replace("_uniprot_ids_IN", ""): f for f in ids_dir.glob("*_uniprot_ids_IN.txt")}
    ge_files = {f.stem.replace("_uniprot_ids_GE", ""): f for f in ids_dir.glob("*_uniprot_ids_GE.txt")}

    shared_keys = set(in_files).intersection(set(ge_files))
    for key in shared_keys:
        in_path = in_files[key]
        ge_path = ge_files[key]
        combined_path = ids_dir / f"{key}_uniprot_ids_COMBINED.txt"

        ids = set()
        with in_path.open("r", encoding="utf-8") as f:
            ids.update(line.strip() for line in f if line.strip())
        with ge_path.open("r", encoding="utf-8") as f:
            ids.update(line.strip() for line in f if line.strip())

        with combined_path.open("w", encoding="utf-8") as out_f:
            for uid in sorted(ids):
                out_f.write(uid + "\n")
        print(f"Combined {key}: wrote {len(ids)} unique UniProt IDs to {combined_path}")

if __name__ == "__main__":
    isolate_uniprot_ids(
        input_dir=Path("pipeline/DELETE"),
        output_dir=Path("pipeline/uniprot_ids"),
        interpro=True
    )
    combine_list(Path(parent_dir / 'pipeline'/ 'uniprot_ids'))
    # print('oo')
