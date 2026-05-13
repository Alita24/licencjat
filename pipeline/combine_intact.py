import pandas as pd
from pathlib import Path

import pandas as pd
from pathlib import Path

def _read_intact_file(path: Path) -> pd.DataFrame:
    """
    Reads an IntAct TSV file, standardizes header names, and cleans up empty columns.
    Returns a DataFrame with canonical column names.
    """
    print(f"- Reading IntAct file: {path}")
    if not path.exists():
        print(f"  File {path} does not exist.")
        return pd.DataFrame()

    # Maps various IntAct header variants to canonical names
    canonical_map = {
        "# ID(s) interactor A": "idA",
        "#ID(s) interactor A": "idA",
        "ID(s) interactor A": "idA",
        "ID(s) interactor B": "idB",
        "Alt. ID(s) interactor A": "altIdsA",
        "Alt. ID(s) interactor B": "altIdsB",
        "Alias(es) interactor A": "aliasesA",
        "Alias(es) interactor B": "aliasesB",
        "Interaction detection method(s)": "detectionMethod",
        "Publication 1st author(s)": "firstAuthor",
        "Publication Identifier(s)": "publicationIdentifiers",
        "Taxid interactor A": "taxidA",
        "Taxid interactor B": "taxidB",
        "Interaction type(s)": "interactionType",
        "Source database(s)": "sourceDatabase",
        "Interaction identifier(s)": "interactionIdentifiers",
        "Confidence value(s)": "confidence",
    }

    try:
        # Try reading the file as tab-separated values
        df = pd.read_csv(path, sep="\t")

        # Check for presence of header; handle if header is missing or formatted oddly
        if not df.columns.str.contains("ID\(s\)").any():
            print(f"  Header not detected, retrying read for {path}")
            df = pd.read_csv(path, sep="\t")

        # Normalize headers and apply canonical mapping
        print(f"- Normalizing column headers for {path}")
        normalized_cols = {}
        for col in df.columns:
            # Remove problematic characters and normalize (robust against weird formatting)
            norm = (
                col.replace("#", "")
                   .replace("(", "")
                   .replace(")", "")
                   .replace(" ", "")
                   .replace("-", "")
                   .replace(".", "")
                   .replace("_", "")
                   .lower()
            )
            # Map to canonical or keep as-is if not mapped
            mapped = canonical_map.get(col, 
                       canonical_map.get(norm, col))
            normalized_cols[col] = mapped
        print(f"  Mapped headers: {normalized_cols}")

        df = df.rename(columns=normalized_cols)

        # Drop columns that are entirely empty or blank
        df = df.dropna(axis=1, how="all")
        df = df.loc[:, ~(df == "").all()]

        print(f"- Finished processing {path}")
        return df

    except Exception as e:
        print(f"! Error reading {path}: {e}")
        return pd.DataFrame()

def _deduplicate_interactions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate IntAct interaction DataFrame so each interaction (idA, idB, interactionType) is unique within each file,
    regardless of order (idA-idB vs idB-idA), assuming 'idA' and 'idB' are strings.
    The deduplication is performed within a given file before combining with others.
    """
    print(f"Deduplicating interaction DataFrame of shape {df.shape} ...")
    if df.empty or not all(col in df.columns for col in ["idA", "idB"]):
        print("DataFrame is empty or missing idA/idB columns, skipping deduplication.")
        return df

    # Canonicalize the (idA, idB) pairs so that idA always <= idB lexicographically
    def canonical_pair(row):
        id1 = str(row["idA"])
        id2 = str(row["idB"])
        sorted_pair = sorted([id1, id2])
        return tuple(sorted_pair)

    canon_pairs = df.apply(canonical_pair, axis=1)
    df = df.copy()
    df["__canon_pair__"] = canon_pairs

    # Optionally include interactionType in deduplication if present
    if "interactionType" in df.columns:
        dedup_keys = ["__canon_pair__", "interactionType"]
    else:
        dedup_keys = ["__canon_pair__"]

    df_dedup = df.drop_duplicates(subset=dedup_keys)
    df_dedup = df_dedup.drop(columns="__canon_pair__")
    print(f"Deduplicated DataFrame to shape {df_dedup.shape}")

    return df_dedup

def combine_intact_gene_interpro(uniprot_dir: str | Path, gene_dir: str | Path, output_dir: str | Path):
    """
    Synchronizes protein interaction data from UniProt and Gene-centric sources.
    """
    print("Starting combination of IntAct gene and UniProt tables...")
    uniprot_dir, gene_dir, output_dir = Path(uniprot_dir), Path(gene_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Set output directory: {output_dir}")

    print("Listing genes in UniProt and Gene result folders...")
    uniprot_genes = {f.stem.split("_")[1] for f in uniprot_dir.glob("intact_*_UNIPROT.tsv")}
    gene_genes = {f.stem.split("_")[1] for f in gene_dir.glob("intact_*_GENE.tsv")}
    all_genes = sorted(uniprot_genes | gene_genes)
    print(f"Identified {len(all_genes)} unique gene(s): {all_genes}")

    all_canon_columns = [
        "idA", "idB", "altIdsA", "altIdsB", "aliasesA", "aliasesB",
        "detectionMethod", "firstAuthor", "publicationIdentifiers",
        "taxidA", "taxidB", "interactionType", "sourceDatabase",
        "interactionIdentifiers", "confidence",
    ]

    for gene_name in all_genes:
        print(f"\nProcessing gene: {gene_name}")
        dfs = []
        files = [
            uniprot_dir / f"intact_{gene_name}_UNIPROT.tsv",
            gene_dir / f"intact_{gene_name}_GENE.tsv"
        ]

        for file_path in files:
            print(f"  Checking file: {file_path}")
            if file_path.exists():
                df = _read_intact_file(file_path)
                if not df.empty:
                    print(f"  Deduplicating {file_path} ...")
                    # Individual deduplication
                    df = _deduplicate_interactions(df)
             
                    dfs.append(df)
                else:
                    print(f"  {file_path} is empty after reading.")
            else:
                print(f"  {file_path} does not exist. Skipping.")

        if dfs:
            print(f"  Aligning columns and concatenating dataframes for gene {gene_name} ...")
            dfs_aligned = [df.reindex(columns=all_canon_columns) for df in dfs]
            df_combined = pd.concat(dfs_aligned, ignore_index=True)
            
            print(f"  Deduplicating combined dataframe for gene {gene_name} ...")
            # Global deduplication across both sources
            df_combined = _deduplicate_interactions(df_combined)

            out_path = output_dir / f"intact_{gene_name}_combined.tsv"
            print(f"  Writing combined results to {out_path} ...")
            df_combined.to_csv(out_path, sep="\t", index=False)
        else:
            print(f"  No dataframes to combine for gene {gene_name}, skipping.")

    print("Combination complete.")


if __name__ == "__main__":
    # Usage Example
    combine_intact_gene_interpro(
        uniprot_dir="pipeline/intact_uniprot_results",
        gene_dir="pipeline/intact_gene_results",
        output_dir="pipeline/combined_intact_results"
    )