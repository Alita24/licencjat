import sys
from pathlib import Path

# Ensure project root is in sys.path for imports
parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

from alphaknot.api import fetch_alphaknot_for_gene
from targets import GENE_TO_IPRS

def alphaknot_pipeline(out_dir: Path):
    """
    pobiera dane z AlphaKnot dla danych genów i IPR
    i zapisuje je w pliku TSV
    Args:
        out_dir: ścieżka do katalogu, w którym będą zapisywane wyniki
    Returns:
        None
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print("Fetching AlphaKnot data for genes/IPRs and writing results...")
    for gene_name, iprs in GENE_TO_IPRS.items():
        fetch_alphaknot_for_gene(gene_name=gene_name, iprs=iprs, out_dir=out_dir)
    print("Done fetching AlphaKnot results.\n")

if __name__ == "__main__":
    alphaknot_pipeline(out_dir=Path("alphaknot_results"))
    # print('lets go')
