import sys
from pathlib import Path
parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

from alphaknot.api import fetch_alphaknot_for_gene
from pipeline.targets import GENE_TO_IPRS, normalize_gene_name


def main() -> None:
    """
    Pobiera pliki TSV z AlphaKnot dla genów/IPR zdefiniowanych w `pipeline/`
    i zapisuje je w katalogu `alphaknot_results`.
    """
    out_dir = Path(__file__).resolve().parent / "alphaknot_results"
    for gene_name, iprs in GENE_TO_IPRS.items():
        canonical_name = normalize_gene_name(gene_name)
        fetch_alphaknot_for_gene(gene_name=canonical_name, iprs=iprs, out_dir=out_dir)


if __name__ == "__main__":
    main()
