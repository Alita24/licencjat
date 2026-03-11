import sys
from pathlib import Path


parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))
    
from typing import Dict, Set

import pandas as pd
import re

from pipeline.targets import GENE_TO_IPRS, normalize_gene_name


PROJECT_ROOT = Path(__file__).resolve().parents[1]

# We assume:
# - IntAct MITAB results per family gene are here:
INTACT_RESULTS_DIR = PROJECT_ROOT / "pipeline" / "intact_results"
# - InterPro-derived protein lists per family gene are here (from api_interpro.py):
INTERPRO_PROTEINS_DIR = parent_dir / "pipeline" / "interpro_results"
# - Output lists of interaction partners will go here:
OUT_DIR = parent_dir / "pipeline" / "intact_partners"


UNIPROTKB_RE = re.compile(r"uniprotkb:([A-Z0-9]+)")


def extract_uniprot_ids(field: str | None) -> Set[str]:
    """
    Extract all UniProt accessions of the form 'uniprotkb:XXXXX' from a MITAB field.
    """
    if not field or not isinstance(field, str):
        return set()
    return set(UNIPROTKB_RE.findall(field))


def load_group_proteins(group_name: str) -> Set[str]:
    """
    Load the set of UniProt accessions that belong to a given family group,
    using the InterPro aggregation output (protein_accession column).
    """
    tsv_path = INTERPRO_PROTEINS_DIR / f"{group_name}_interpro.tsv"
    if not tsv_path.exists():
        print(f"[WARN] InterPro file not found for group '{group_name}': {tsv_path}")
        return set()

    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    if "protein_accession" not in df.columns:
        print(f"[WARN] 'protein_accession' column missing in {tsv_path}")
        return set()

    return set(df["protein_accession"].dropna().unique())


def collect_partners_for_group(group_name: str) -> Set[str]:
    """
    For a given family gene group (e.g. trmD, trm5, nep1), read its IntAct results
    and return a set of UniProt IDs that are interaction partners (i.e. the
    proteins on the other side of the interaction, excluding the group proteins
    themselves when we can identify them).
    """
    canonical_name = normalize_gene_name(group_name)
    intact_path = INTACT_RESULTS_DIR / f"{canonical_name}_intact.tsv"
    if not intact_path.exists():
        print(f"[WARN] IntAct results not found for group '{group_name}': {intact_path}")
        return set()

    df = pd.read_csv(intact_path, sep="\t", dtype=str)

    group_proteins = load_group_proteins(canonical_name)
    partners: Set[str] = set()

    for _, row in df.iterrows():
        ids_a = extract_uniprot_ids(row.get("idA"))
        ids_b = extract_uniprot_ids(row.get("idB"))

        all_ids = ids_a | ids_b

        if group_proteins:
            # Separate family proteins from their partners
            family_side = all_ids & group_proteins
            partner_side = all_ids - family_side
        else:
            # Fallback: we don't know which are family proteins, so treat all as partners
            partner_side = all_ids

        partners.update(partner_side)

    return partners


def write_partners_lists() -> Dict[str, Path]:
    """
    Create per-group lists of unique UniProt IDs of interaction partners,
    based on IntAct MITAB results and (optionally) InterPro protein lists.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    result_paths: Dict[str, Path] = {}

    for group_name in GENE_TO_IPRS.keys():
        partners = collect_partners_for_group(group_name)
        out_txt = OUT_DIR / f"{group_name}_interaction_partners.txt"

        with out_txt.open("w", encoding="utf-8") as f:
            for acc in sorted(partners):
                f.write(f"{acc}\n")

        result_paths[group_name] = out_txt
        print(f"{group_name}: wrote {len(partners)} unique partner UniProt IDs to {out_txt}")

    return result_paths


def main() -> None:
    write_partners_lists()


if __name__ == "__main__":
    main()
