import requests
import pandas as pd
from io import StringIO

def fetch_intact_interactions(query, out_tsv_path=None):
    """
    Fetch interaction data from IntAct for the specified query.

    Parameters:
        query (str): The query string, e.g. a UniProt ID or gene name.
        out_tsv_path (str or None): If provided, saves results as a TSV file at this path.

    Returns:
        pd.DataFrame: DataFrame of MITAB 2.5 data.
    """
    url = (
        "https://www.ebi.ac.uk/Tools/webservices/psicquic/intact/"
        f"webservices/current/search/query/{query}?format=tab25"
    )
    print(url)
    resp = requests.get(url)
    resp.raise_for_status()

    df = pd.read_csv(StringIO(resp.text), sep="\t", header=None, dtype=str)

    # standard mitab 2.5 column names
    mitab25_cols = [
        "idA",                # 0
        "idB",                # 1
        "altIdsA",            # 2
        "altIdsB",            # 3
        "aliasesA",           # 4
        "aliasesB",           # 5
        "detectionMethod",    # 6
        "firstAuthor",        # 7
        "publicationIdentifiers", # 8
        "taxidA",             # 9
        "taxidB",             # 10
        "interactionType",    # 11
        "sourceDatabase",     # 12
        "interactionIdentifiers", # 13
        "confidence"          # 14
    ]

    num_cols = df.shape[1]
    if num_cols <= len(mitab25_cols):
        df.columns = mitab25_cols[:num_cols]
    else:
        extra = [f"col{i}" for i in range(len(mitab25_cols), num_cols)]
        df.columns = mitab25_cols + extra

    if out_tsv_path:
        df.to_csv(out_tsv_path, sep="\t", index=False)
    return df

if __name__ == "__main__":
    # Example usage: fetch for gene 'Trmd' and save output to 'intact_Trmd.tsv'
    query = "trmd"
    out_file = f"try_intact_{query}.tsv"
    print(f"Fetching IntAct interactions for query: '{query}'")
    df = fetch_intact_interactions(query, out_tsv_path=out_file)
    print(f"Retrieved {len(df)} interactions. Saved to {out_file}")
