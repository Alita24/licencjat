import sys
from pathlib import Path

# Ensure parent directory is in sys.path for relative imports.
parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

from alphaknot.api import fetch_alphaknot_for_gene
from pipeline.targets import GENE_TO_IPRS, normalize_gene_name

def fetch_and_save_alphaknot_results(out_dir: Path):
    """
    For each gene/IPR mapping, fetch AlphaKnot results and save to `out_dir`.
    """
    for gene_name, iprs in GENE_TO_IPRS.items():
        canonical_name = normalize_gene_name(gene_name)
        fetch_alphaknot_for_gene(gene_name=canonical_name, iprs=iprs, out_dir=out_dir)

def extract_uniprot_ids_from_file(filepath: Path) -> set[str]:
    """
    Extract unique UniProt accessions from a results file (TSV).
    Ignores comment lines and expects the UniProt accession in column 3.
    """
    ids = set()
    with filepath.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split('\t')
            if len(cols) >= 3:
                ids.add(cols[2])
    return ids

def isolate_uniprot_ids(alphaknot_results_dir: Path, out_dir: Path):
    """
    For each .tsv in `alphaknot_results_dir`, extract UniProt IDs
    and write to `{out_dir}/{inputfile}_uniprot_ids.txt`.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    for file in alphaknot_results_dir.glob("*.tsv"):
        uniprot_ids = extract_uniprot_ids_from_file(file)
        out_path = out_dir / f"{file.stem}_uniprot_ids.txt"
        with out_path.open("w", encoding="utf-8") as out_f:
            for uid in uniprot_ids:
                if uid.endswith('-F1'):
                    uid = uid[:-3]
                out_f.write(uid + "\n")
        print(f"{file.name}: wrote {len(uniprot_ids)} UniProt IDs to {out_path}")

def main():
    """
    Complete AlphaKnot pipeline:
    1. Fetch AlphaKnot results for all defined gene/IPR sets.
    2. Extract UniProt IDs from each result and save for later use.
    """
    this_dir = Path(__file__).resolve().parent
    results_dir = this_dir / "alphaknot_results"
    uniprot_ids_dir = this_dir / "alphaknot_uniprot_ids"

    print("Fetching AlphaKnot data for genes/IPRs...")
    fetch_and_save_alphaknot_results(results_dir)
    print("Done fetching AlphaKnot results.")

    print("Extracting UniProt IDs from AlphaKnot results...")
    isolate_uniprot_ids(results_dir, uniprot_ids_dir)
    print("Done extracting UniProt IDs.")

if __name__ == "__main__":
    main()
