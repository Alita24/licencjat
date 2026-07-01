import sys
from pathlib import Path
import requests

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

from data.targets import GENE_TO_IPRS 

def intact_gene_pipeline(output_dir):
    """
    For each gene family found in pipeline.targets.GENE_TO_IPRS, downloads IntAct MITAB25 results.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    url = "https://www.ebi.ac.uk/intact/ws/graph/export/interaction/list"

    for fam_name in GENE_TO_IPRS.keys():
        out_fn = output_dir / f"intact_{fam_name}_FAM.tsv"
        if out_fn.exists():
            print(f"Skipping {fam_name}: file already exists.")
            continue
        params = {
            "query": fam_name,
            "format": "miTab25"
        }
        print(f"Fetching IntAct data for gene family: {fam_name}")
        try:
            response = requests.post(url, data=params)
            if response.status_code == 200:
            
                with open(out_fn, "wb") as file:
                    file.write(response.content)
                print(f"Downloaded for {fam_name}: {out_fn}")
            else:
                print(f"Failed to download for {fam_name}: status code {response.status_code}")
        except Exception as e:
            print(f"Error fetching {fam_name}: {e}")

if __name__ == "__main__":
    # results will be saved to intact_gene_results/
    intact_gene_pipeline(input_dir=Path(__file__).resolve().parent / "alphaknot_results", output_dir=Path(__file__).resolve().parent / "intact_gene_results")