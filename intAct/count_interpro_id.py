# -----------------------------------------------------------------------------
# Ten skrypt służy do analizy pliku interakcji białek w formacie MITAB, w którym 
# interaktory są zidentyfikowane przez identyfikatory UniProt. Skrypt mapuje każdy 
# identyfikator UniProt na odpowiadające mu domeny InterPro (używając pliku mapowania) 
# i liczy, jak często każda para domen InterPro występuje razem w interakcjach białek. 
#
# Wynikiem działania skryptu są pliki CSV zestawiające liczbę wystąpień poszczególnych 
# par InterPro, informacje o organizmach i identyfikatory UniProt, które przyczyniły się 
# do danego połączenia. 
#-----------------------------------------------------------------------------

#!/usr/bin/env python3
"""
mitab_to_interpro_pairs.py

Parse MITAB-like interactions, map UniProt -> InterPro and count InterPro-InterPro occurrences.

Inputs (edit filenames below):
- MITAB_FILE: tab-separated interactions (columns 0 and 1 = interactor A and B)
- INTERPRO_MAP: CSV with columns "UniProt_ID", "InterPro_ID", "InterPro_Name", "organism"

Outputs:
- interpro_interactions_counts.csv  (Count first)
- interpro_matrix.tsv               (pivot matrix: rows InterPro_A, cols InterPro_B)
- (optional) heatmap_topN.png       (visualize top N InterPro pairs)
- family-specific CSVs if you provide family InterPro lists
"""

import re
from collections import defaultdict, Counter
from itertools import product
import pandas as pd
import numpy as np

# -------------------------
# Config - set filenames
# -------------------------
gene = 'nep1'

MITAB_FILE = f"intAct_{gene}.tsv"        # your downloaded MITAB-like file
INTERPRO_MAP = f"uniprot_to_interpro_{gene}.csv"     # produced earlier: UniProt_ID,InterPro_ID,InterPro_Name,organism

OUTPUT_COUNTS = f"interpro_{gene}_interactions_counts.csv"

# If you want filtered outputs for families (fill with InterPro IDs)
TRMD_INTERPRO = []   # e.g. ["IPRxxxx"]
TRM5_INTERPRO = []
NEP1_INTERPRO  = []

# Sleep / rate-limits not needed here because we parse local file & mapping.


# -------------------------
# Helpers
# -------------------------
_uniprot_re = re.compile(r"uniprotkb:([A-Z0-9\-]+)", re.IGNORECASE)

def extract_uniprot_ids(field_text):
    """Extract all UniProt accessions from a MITAB field (e.g. 'uniprotkb:Q9W4J5|intact:EBI-...')."""
    if not isinstance(field_text, str):
        return set()
    return set(m.group(1) for m in _uniprot_re.finditer(field_text))


def load_interpro_map(map_csv):
    """
    Load UniProt -> list of (InterPro ID, InterPro Name, organism)
    Expect CSV with headers: UniProt_ID,InterPro_ID,InterPro_Name,organism
    """
    prot_to_ipr = defaultdict(list)
    ipr_names = dict()
    ipr_organism = dict()
    df = pd.read_csv(map_csv, dtype=str).fillna("")
    expected_cols = {"UniProt_ID", "InterPro_ID"}
    if not expected_cols.issubset(set(df.columns)):
        raise ValueError(f"Mapping CSV must contain columns: {expected_cols}. Found: {df.columns.tolist()}")

    for _, r in df.iterrows():
        uid = r["UniProt_ID"].strip()
        ipr = r["InterPro_ID"].strip()
        # InterPro_Name and organism always present, possibly empty string
        name = r["InterPro_Name"].strip() if "InterPro_Name" in r else ""
        organism = r["organism"].strip() if "organism" in r else ""
        if uid and ipr:
            prot_to_ipr[uid].append((ipr, organism))
            if ipr not in ipr_names and name:
                ipr_names[ipr] = name
            # Only keep one representative organism per InterPro, notably the first found
            if ipr not in ipr_organism and organism:
                ipr_organism[ipr] = organism
    return prot_to_ipr, ipr_names, ipr_organism


# -------------------------
# Main pipeline
# -------------------------
def main():
    print("Loading UniProt->InterPro map...")
    prot2ipr, ipr_names, ipr_organism = load_interpro_map(INTERPRO_MAP)
    print(f"Loaded mapping for {len(prot2ipr)} UniProt IDs, {len(ipr_names)} unique InterPro IDs.")

    # counter for domain-domain pairs
    pair_counter = Counter()

    # also keep track of which UniProt pairs contributed to each IPR pair (optional)
    contributors = defaultdict(set)  # key = (iprA, iprB) -> set of "UIDa|UIDb"
    org_contributors = defaultdict(set)  # key = (iprA, iprB) -> set of tuples (orgA, orgB)

    print("Parsing MITAB file and building IPR pairs...")
    with open(MITAB_FILE, "r", encoding="utf-8") as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 2:
                continue

            fieldA = cols[0]
            fieldB = cols[1]

            uidsA = extract_uniprot_ids(fieldA)
            uidsB = extract_uniprot_ids(fieldB)

            if not uidsA or not uidsB:
                # skip interactions where we cannot find UniProt IDs on either side
                continue

            # For all UniProt pairs, map to all InterPro domains (and collect organisms)
            for uidA, uidB in product(uidsA, uidsB):
                iprsA_list = prot2ipr.get(uidA, [])
                iprsB_list = prot2ipr.get(uidB, [])
                if not iprsA_list or not iprsB_list:
                    continue

                for (iprA, orgA), (iprB, orgB) in product(iprsA_list, iprsB_list):
                    key = tuple(sorted((iprA, iprB)))
                    # sort organisms in the same way as key
                    if iprA <= iprB:
                        org_key = (orgA, orgB)
                    else:
                        org_key = (orgB, orgA)
                    pair_counter[key] += 1
                    contributors[key].add(f"{uidA}|{uidB}")
                    org_contributors[key].add(org_key)

            if line_no % 50000 == 0:
                print(f"  processed {line_no} lines...")

    print(f"Total distinct InterPro pairs found: {len(pair_counter)}")

    # Build DataFrame for output
    rows = []
    for (iprA, iprB), count in pair_counter.items():
        nameA = ipr_names.get(iprA, "")
        nameB = ipr_names.get(iprB, "")

        # Build semicolon-separated organism info from org_contributors for (iprA, iprB)
        org_pairs_set = org_contributors.get((iprA, iprB), set())
        # If present, join them, else use empty string
        organism_pairs = ";".join(f"{oA}|{oB}" for oA, oB in sorted(org_pairs_set)) if org_pairs_set else ""

        # optional: join contributor UniProts as semicolon list (could get large)
        contribs = ";".join(sorted(contributors[(iprA, iprB)])) if (iprA, iprB) in contributors else ""
        rows.append({
            "Count": count,
            "InterPro_A": iprA,
            "InterPro_B": iprB,
            "InterPro_A_Name": nameA,
            "InterPro_B_Name": nameB,
            "Organism_Pairs": organism_pairs,
            "UniProt_Contributors": contribs
        })

    df_out = pd.DataFrame(rows)
    if df_out.empty:
        print("No InterPro pairs produced. Exiting.")
        return

    # Sort by Count descending (Count first as user requested earlier)
    df_out = df_out.sort_values("Count", ascending=False)
    df_out.to_csv(OUTPUT_COUNTS, index=False)
    print("Wrote counts to:", OUTPUT_COUNTS)

    # Optionally produce family-specific csvs
    def save_family_filtered(family_iprs, label):
        if not family_iprs:
            return
        mask = df_out["InterPro_A"].isin(family_iprs) | df_out["InterPro_B"].isin(family_iprs)
        df_f = df_out[mask].sort_values("Count", ascending=False)
        if not df_f.empty:
            fname = f"{label}_interactions.csv"
            df_f.to_csv(fname, index=False)
            print("Wrote family filtered:", fname)

    save_family_filtered(TRMD_INTERPRO, "TRMD")
    save_family_filtered(TRM5_INTERPRO, "TRM5")
    save_family_filtered(NEP1_INTERPRO, "NEP1")

    print("DONE.")


if __name__ == "__main__":
    main()
