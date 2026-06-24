import sys
from pathlib import Path

parent_dir = Path(__file__).resolve().parents[1]
here = Path(__file__).resolve().parents[0]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))



import pandas as pd
from interpro.api import get_reviewed_proteins_for_interpro, write_combined_csv
from data.targets import GENE_TO_IPRS

def download_ipr_csv(interpro_id: str, out_dir: Path) -> Path:
    """
    Pobiera z InterPro wszystkie białka powiązane z danym IPR
    i zapisuje je jako pojedynczy plik CSV w out_dir.
    Args:
        interpro_id: identyfikator IPR
        out_dir: ścieżka do katalogu, w którym będą zapisywane wyniki
    Returns:
        ścieżka do pliku CSV
    """
    out_dir.mkdir(parents=True, exist_ok=True)
        # out_dir - Path object; parents=True - creates parent directories if they don't exist; exist_ok=True - doesn't raise an error if the directory already exists
    
    csv_path = out_dir / f"{interpro_id}.csv"

    if csv_path.exists():
        print(f"{interpro_id}: CSV already exists at {csv_path}, skipping download.")
        return csv_path

    proteins = get_reviewed_proteins_for_interpro(interpro_id)
    write_combined_csv(proteins, output_file=str(csv_path))
    return csv_path


def gather_proteins_for_group(group_name: str, ipr_list: list[str], out_dir: Path) -> pd.DataFrame:
    """
    zbiera dane z InterPro dla danych genów i IPR
    Args:
        group_name: nazwa grupy białek
        ipr_list: lista identyfikatorów IPR
        out_dir: ścieżka do katalogu, w którym będą zapisywane wyniki
    Returns:
        DataFrame
    """
    dfs: list[pd.DataFrame] = []
    
    for ipr in ipr_list:
        ipr_path = out_dir / f"{ipr}.csv"
        if not ipr_path.exists():
            print(f"{group_name}: {ipr_path} not found, downloading data for {ipr}")
            ipr_path = download_ipr_csv(ipr, out_dir)

        if ipr_path.exists():
            df = pd.read_csv(ipr_path, comment="#")
            if not df.empty:
                df = df.copy()
                df["protein_group"] = group_name
                dfs.append(df)

    if dfs:
        # concatenates all dataframes in the list dfs into a single dataframe
        return pd.concat(dfs, ignore_index=True)
            # ignore_index=True - resets the index of the concatenated dataframe
    return pd.DataFrame()


def interpro_pipeline(
    out_dir: Path = None
):
    """
    funkcja pobiera dane z InterPro dla danych genów i IPR
    i zapisuje je w pliku TSV
    Args:
        out_dir: ścieżka do katalogu, w którym będą zapisywane wyniki
    Returns:
        None
    """
    print("Starting interpro_pipeline...")
    if out_dir is None:
        out_dir = parent_dir /'pipeline' / "interpro_results"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    alphaknot_dir = here / Path("alphaknot_results")
    for group_name, ipr_list in GENE_TO_IPRS.items():
        print(f"\nProcessing group: {group_name} with IPRs: {ipr_list}")

        alphaknot_file = alphaknot_dir / f"{group_name}.tsv"
        print(alphaknot_file)
        if alphaknot_file.exists():
            print(f"Group {group_name} already present in alphaknot_results ({alphaknot_file}), skipping...")
            continue

        df_group = gather_proteins_for_group(group_name, ipr_list, out_dir=out_dir)
        if not df_group.empty:
            output_file = out_dir / f"{group_name}.tsv"
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
            # tworzy df z kolumnami available_columns i zapisuje go do output_file
            print(f"Wrote {len(df_group)} records for {group_name} to {output_file}")
        else:
            print(f"No entries found for group {group_name}")

    print("interpro_pipeline completed.")

if __name__ == "__main__":
    interpro_pipeline(out_dir=Path("pipeline/interpro_results"))
    # print('kk')