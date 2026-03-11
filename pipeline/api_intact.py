import sys
from pathlib import Path


parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

from typing import Dict

import pandas as pd

from intAct.api import fetch_intact_interactions
from pipeline.targets import GENE_TO_IPRS, normalize_gene_name, get_all_names_for_group


# Katalog na wyniki z IntAct
OUT_DIR = parent_dir / "pipeline" / "intact_results"

def fetch_intact_for_gene(gene_name: str) -> Path:
    """
    Pobiera interakcje z IntAct dla podanej nazwy genu
    i zapisuje je jako plik TSV w OUT_DIR.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    canonical_name = normalize_gene_name(gene_name)
    query_names = get_all_names_for_group(gene_name)

    dfs: list[pd.DataFrame] = []
    for query in query_names:
        print(f"Fetching IntAct interactions for query '{query}' (group: {canonical_name})")
        df = fetch_intact_interactions(query, out_tsv_path=None)
        if not df.empty:
            df = df.copy()
            df["intact_query"] = query
            dfs.append(df)

    if dfs:
        combined = pd.concat(dfs, ignore_index=True).drop_duplicates()
    else:
        combined = pd.DataFrame()

    out_path = OUT_DIR / f"{canonical_name}_intact.tsv"
    combined.to_csv(out_path, sep="\t", index=False)
    print(
        f"{canonical_name}: retrieved {len(combined)} interactions "
        f"combined from queries {query_names} -> {out_path}"
    )
    return out_path


def fetch_all_intact_for_targets() -> Dict[str, Path]:
    """
    Dla wszystkich genów zdefiniowanych w pipeline.targets.GENE_TO_IPRS
    pobiera dane z IntAct na podstawie nazwy genu.

    Zwraca słownik gene_name -> ścieżka do pliku TSV z wynikami.
    """
    results: Dict[str, Path] = {}
    for gene_name in GENE_TO_IPRS.keys():
        results[gene_name] = fetch_intact_for_gene(gene_name)
    return results


def main() -> None:
    fetch_all_intact_for_targets()


if __name__ == "__main__":
    main()
