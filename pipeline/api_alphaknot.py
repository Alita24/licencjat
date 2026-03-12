import sys
from pathlib import Path

# Ensure project root is in sys.path for imports
parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

from alphaknot.api import fetch_alphaknot_for_gene
from pipeline.targets import GENE_TO_IPRS

def alphaknot_pipeline(out_dir: Path):
    """
    Orchestrates the complete AlphaKnot workflow:
    1. Fetch AlphaKnot data for all gene/IPR sets, storing TSVs in 'alphaknot_results'..
    """
    results_dir = out_dir

    results_dir.mkdir(parents=True, exist_ok=True)

    print("Fetching AlphaKnot data for genes/IPRs and writing results...")
    for gene_name, iprs in GENE_TO_IPRS.items():
        fetch_alphaknot_for_gene(gene_name=gene_name, iprs=iprs, out_dir=results_dir)
    print("Done fetching AlphaKnot results.\n")

if __name__ == "__main__":
    alphaknot_pipeline(out_dir=Path("alphaknot_results"))
