import sys
from pathlib import Path

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))


import pandas as pd
from interpro.api import get_reviewed_proteins_for_interpro, write_combined_csv
from pipeline.targets import GENE_TO_IPRS

# Katalog na surowe pliki CSV z InterPro dla poszczególnych IPR
RAW_INTERPRO_DIR = parent_dir / "interpro"

# Katalog na zintegrowane TSV per grupa (trmD, trm5, nep1, ...)
OUT_DIR = parent_dir / "pipeline" / "interpro_results"


def download_ipr_csv(interpro_id: str) -> Path:
    """
    Pobiera z InterPro wszystkie białka powiązane z danym IPR
    i zapisuje je jako pojedynczy plik CSV w RAW_INTERPRO_DIR.
    """
    RAW_INTERPRO_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = RAW_INTERPRO_DIR / f"{interpro_id}.csv"

    if csv_path.exists():
        print(f"{interpro_id}: CSV already exists at {csv_path}, skipping download.")
        return csv_path

    proteins = get_reviewed_proteins_for_interpro(interpro_id)
    write_combined_csv(proteins, output_file=str(csv_path))
    return csv_path


def gather_proteins_for_group(group_name: str, ipr_list: list[str]) -> pd.DataFrame:
    """
    Dla danej grupy białek, agreguje wszystkie wpisy z jej określonych IPR.
    Jeśli brakuje lokalnego CSV dla danego IPR, najpierw go pobiera.
    """
    dfs: list[pd.DataFrame] = []

    for ipr in ipr_list:
        ipr_path = RAW_INTERPRO_DIR / f"{ipr}.csv"
        if not ipr_path.exists():
            print(f"{group_name}: {ipr_path} not found, downloading data for {ipr}")
            ipr_path = download_ipr_csv(ipr)

        if ipr_path.exists():
            df = pd.read_csv(ipr_path, comment="#")
            if not df.empty:
                df = df.copy()
                df["protein_group"] = group_name
                dfs.append(df)

    if dfs:
        return pd.concat(dfs, ignore_index=True)
    return pd.DataFrame()


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    for group_name, ipr_list in GENE_TO_IPRS.items():
        df_group = gather_proteins_for_group(group_name, ipr_list)
        if not df_group.empty:
            output_file = OUT_DIR / f"{group_name}_interpro.tsv"
            core_columns = [
                "protein_accession",
                "protein_name",
                "protein_length",
                "gene",
                "organism_taxid",
                "organism_name",
                "organism_fullname",
                "in_alphafold",
                "interpro_accession",
                "entry_type",
                "entry_integrated",
                "interpro_protein_length",
                "fragment_start",
                "fragment_end",
                "fragment_status",
                "protein_group",
            ]
            available_columns = [col for col in core_columns if col in df_group.columns]
            df_group[available_columns].to_csv(output_file, sep="\t", index=False)
            print(f"Wrote {len(df_group)} records for {group_name} to {output_file}")
        else:
            print(f"No entries found for group {group_name}")


if __name__ == "__main__":
    main()