from pathlib import Path
import pandas as pd
from api_alphaknot import alphaknot_pipeline
from api_interpro import interpro_pipeline
from api_intactUNI import intact_uniprot_pipeline
from isolate_uniprot_ids import isolate_uniprot_ids
from isolate_interactors import isolate_interactors
from annotate_interactors import annotate_interactors

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
    print("Preparing alphaknot results...")
    alphaknot_pipeline(
        out_dir=here / "alphaknot_results"
    )
    print("Alphaknot results prepared.")

    # 2. Download group proteins from InterPro
    print("Downloading group proteins from InterPro...")
    interpro_pipeline() 
    print("Group proteins downloaded from InterPro.")

    # 3. Extract UniProt IDs per group from alphaknot (should already be in uniprot_ids/*.txt)
    print("Extracting UniProt IDs per group from alphaknot...")
    isolate_uniprot_ids(
        input_dir=here / "alphaknot_results",
        output_dir=here / "uniprot_ids"
    )
    isolate_uniprot_ids(
        input_dir=here / "interpro_results",
        output_dir=here / "uniprot_ids"
    )
    print("UniProt IDs extracted from alphaknot.")

    # 4. Fetch IntAct interaction data per UniProt ID
    print("Fetching IntAct interaction data per UniProt ID...")
    intact_uniprot_pipeline()