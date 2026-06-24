from pathlib import Path
import pandas as pd

from api_alphaknot import alphaknot_pipeline
from api_interpro import interpro_pipeline

from api_intactUNI import intact_uniprot_pipeline
from api_intactGENE import intact_gene_pipeline

from isolate_uniprot_ids import isolate_uniprot_ids, combine_list
from combine_intact import combine_intact_gene_interpro

from isolate_interactors import isolate_partners
from anotate_interactors import annotate_pipeline, count_interpro


def pipeline_main():
    """
    Tutaj znajduje się główny pipeline służący do analizy danych bazując na nazwie rodziny oraz kodzie InterPro domeny definiującej tą rodzinę.
    
    SPRAWDŹ!: przed puszczeniem potoku czy w data/targets.py są prawidłowe dane w słowniku

    Przykład komendy w terminalu: 
        python3 pipelineFAM.py
    """
    here = Path(__file__).resolve().parent
    
    # 1. alphaknot results - for knotted proteins
    alphaknot_pipeline(
        out_dir=here / "alphaknot_results"
    )


    # 2. download proteins from InterPro for unknotted analogues 
    interpro_pipeline(
        out_dir=here/ "interpro_results"
    )

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

    # 4. Fetch IntAct interaction data per UniProt ID and via family name
    print("Fetching IntAct interaction data per UniProt ID and via family name...")

    intact_uniprot_pipeline(
        ids_dir=here / "uniprot_ids",
        out_dir=here / "intact_fam_uniprot",
        file_pattern="_uniprot_ids_KNOT"
    )
    intact_uniprot_pipeline(
        ids_dir=here / "uniprot_ids",
        out_dir=here / "intact_fam_uniprot",
        file_pattern="_uniprot_ids_UNKNOT"
    )

    intact_gene_pipeline(
        output_dir=here / "intact_fam_name"
    )
    print("IntAct interaction data fetched per UniProt ID and via family name.")

    # 5. Combine IntAct interaction data per UniProt ID and via family name
    print("Combining IntAct interaction data per UniProt ID and via family name...")

    combine_intact_gene_interpro(
        uniprot_dir=here / "intact_fam_uniprot",
        gene_dir=here / "intact_fam_name",
        output_dir=here / "combined_intact_results"
    )
    print("IntAct interaction data combined per UniProt ID and via gene name.")

    # 6. isolate interactors
    isolate_partners(
        intact_dir=here/'intact_fam_uniprot',
        output_dir=here/'isolated_partners',
        add_publication=True
    )

    #7. uniprot --> interpro domains
    annotate_pipeline(
        input_dir=here/'isolated_partners',
    )

    count_interpro(
        input_dir=here/'isolated_partners',
        output_dir=here/'count_files',
    )
    

if __name__ == "__main__":
    pipeline_main()
    