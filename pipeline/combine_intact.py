import pandas as pd
from pathlib import Path

def _read_intact_file(path: Path) -> pd.DataFrame:
    """
    Reads an IntAct interaction file (TSV), detects and skips comment lines/header variants,
    normalizing column headers to a consistent set between GENE and UNIPROT files.
    Also removes columns that are all NaN or empty.
    """
    if not path.exists():
        return pd.DataFrame()
    
    # Try reading header lines: skip initial '#' or comment-rows
    comment_rows = 0
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                comment_rows += 1
            else:
                break

    # Read with Pandas after skipping comments
    df = pd.read_csv(path, sep="\t", skiprows=comment_rows)
    # Normalize whitespace in column names
    df.columns = [c.strip().lstrip("#").replace(" ", "") for c in df.columns]

    # Map header variants from GENE/UNIPROT to a canonical set
    header_map = {
        "ID(s)interactorA": "idA",
        "ID(s)interactorB": "idB",
        "Alt.ID(s)interactorA": "altIdsA",
        "Alt.ID(s)interactorB": "altIdsB",
        "Alias(es)interactorA": "aliasesA",
        "Alias(es)interactorB": "aliasesB",
        "Interactiondetectionmethod(s)": "detectionMethod",
        "Publication1stauthor(s)": "firstAuthor",
        "PublicationIdentifier(s)": "publicationIdentifiers",
        "TaxidinteractorA": "taxidA",
        "TaxidinteractorB": "taxidB",
        "Interactiontype(s)": "interactionType",
        "Sourcedatabase(s)": "sourceDatabase",
        "Interactionidentifier(s)": "interactionIdentifiers",
        "Confidencevalue(s)": "confidence",
    }
    # Rename columns using the map if present
    df = df.rename(columns={c: header_map.get(c, c) for c in df.columns})
    return df

def _deduplicate_interactions(df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate IntAct interaction DataFrame so each interaction (idA, idB, interactionType) is unique within each file,
    regardless of order (idA-idB vs idB-idA), assuming 'idA' and 'idB' are strings.
    The deduplication is performed within a given file before combining with others.
    """
    if df.empty or not all(col in df.columns for col in ["idA", "idB"]):
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
    return df_dedup

def combine_intact_gene_interpro(uniprot_dir: str | Path, gene_dir: str | Path, output_dir: str | Path):
    """
    Combines IntAct interaction tables from UNIPROT-based (uniprot_dir) and GENE-based (gene_dir) sources,
    applying consistent headers, deduplicates using _deduplicate_interactions, and saves the unified results into output_dir by gene name.
    """
    uniprot_dir = Path(uniprot_dir)
    gene_dir = Path(gene_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for uniprot_file in uniprot_dir.glob("intact_*_UNIPROT.tsv"):
        gene_name = uniprot_file.stem.split("_")[1]
        gene_file = gene_dir / f"intact_{gene_name}_GENE.tsv"
        dfs = []

        # Read both files with header normalization and deduplicate each before combining
        df_uni = _read_intact_file(uniprot_file)
        if not df_uni.empty:
            df_uni = _deduplicate_interactions(df_uni)
            dfs.append(df_uni)
            print(f"Read and deduplicated {len(df_uni)} rows from {uniprot_file}")

        df_gene = _read_intact_file(gene_file)
        if not df_gene.empty:
            df_gene = _deduplicate_interactions(df_gene)
            dfs.append(df_gene)
            print(f"Read and deduplicated {len(df_gene)} rows from {gene_file}")

        # Combine only if at least one is non-empty
        if dfs:
            # Find shared columns: use canonical order for output
            all_canon_columns = [
                "idA", "idB", "altIdsA", "altIdsB", "aliasesA", "aliasesB",
                "detectionMethod", "firstAuthor", "publicationIdentifiers",
                "taxidA", "taxidB", "interactionType", "sourceDatabase",
                "interactionIdentifiers", "confidence",
            ]
            present_cols = [c for c in all_canon_columns if any(c in df.columns for df in dfs)]
            dfs_aligned = [df[[c for c in present_cols if c in df.columns]].copy() for df in dfs]
            # Pad missing columns in each df with NaN, then reindex
            dfs_aligned = [df.reindex(columns=present_cols) for df in dfs_aligned]

            # Combine
            df_combined = pd.concat(dfs_aligned, ignore_index=True)

            # Final deduplication step after combine (in case the formats overlapped)
            df_combined = _deduplicate_interactions(df_combined)

            out_path = output_dir / f"intact_{gene_name}_combined.tsv"
            df_combined.to_csv(out_path, sep="\t", index=False)
            print(f"Wrote combined (deduplicated/unified) file for {gene_name}: {out_path}")
        else:
            print(f"No data found for {gene_name}, skipping.")

if __name__ == "__main__":
    # Usage Example
    combine_intact_gene_interpro(
        uniprot_dir="pipeline/intact_uniprot_results",
        gene_dir="pipeline/intact_gene_results",
        output_dir="pipeline/combined_intact_results"
    )