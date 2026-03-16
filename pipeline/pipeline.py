from pathlib import Path
import pandas as pd

from api_alphaknot import alphaknot_pipeline
from api_interpro import interpro_pipeline

from api_intactUNI import intact_uniprot_pipeline
from api_intactGENE import intact_gene_pipeline

from isolate_uniprot_ids import isolate_uniprot_ids, combine_list
from combine_intact import combine_intact_gene_interpro




def pipeline_main():
    """
    Orchestrates the complete pipeline:
    1. Prepare alphaknot results (manual or custom parser)
    2. Download group proteins from InterPro
    3. Extract UniProt IDs per group from alphaknot (should already be in uniprot_ids/*.txt)
    4. Fetch IntAct interaction data per UniProt ID
    5. Isolate interaction partners from IntAct and annotate with InterPro domains
    """
    here = Path(__file__).resolve().parent
    # 1. Prepare alphaknot results
    alphaknot_pipeline(
        out_dir=here / "alphaknot_results"
    )

    # 2. Download group proteins from InterPro
    interpro_pipeline() 

    # 3. Extract UniProt IDs per group from alphaknot (should already be in uniprot_ids/*.txt)
    isolate_uniprot_ids(
        input_dir=here / "alphaknot_results",
        output_dir=here / "uniprot_ids"
    )
    isolate_uniprot_ids(
        input_dir=here / "interpro_results",
        output_dir=here / "uniprot_ids",
        interpro=True
    )
    combine_list(
        ids_dir="uniprot_ids"
    )
    print("UniProt IDs extracted from alphaknot.")

    # 4. Fetch IntAct interaction data per UniProt ID and via gene name
    print("Fetching IntAct interaction data per UniProt ID and via gene name...")
    intact_uniprot_pipeline(
        ids_dir=here / "uniprot_ids",
        out_dir=here / "intact_uniprot_results"
    )
    
    intact_gene_pipeline(
        output_dir=here / "intact_gene_results"
    )
    print("IntAct interaction data fetched per UniProt ID and via gene name.")

    # 5. Combine IntAct interaction data per UniProt ID and via gene name
    print("Combining IntAct interaction data per UniProt ID and via gene name...")
    combine_intact_gene_interpro(
        uniprot_dir=here / "intact_uniprot_results",
        gene_dir=here / "intact_gene_results",
        output_dir=here / "combined_intact_results"
    )
    print("IntAct interaction data combined per UniProt ID and via gene name.")

    # 6. isolate interactors

if __name__ == "__main__":
    pipeline_main()
    