# Protein interaction analysis pipeline

This directory contains scripts for collecting protein interaction data and identifying InterPro domains that occur frequently among the interaction partners of selected protein groups.

## Project configuration and components

- `../data/targets.py` defines the protein groups and their InterPro (IPR) identifiers.
- `api_alphaknot.py` and `api_interpro.py` retrieve protein data from AlphaKnot and InterPro.
- `isolate_uniprot_ids.py` extracts and compares UniProt accession lists.
- `api_intactUNI.py` and `api_intactGENE.py` retrieve IntAct interactions by UniProt accession or group name.
- `combine_intact.py` combines and deduplicates interaction tables.
- `isolate_interactors.py` extracts interaction partners, excluding proteins already listed for the target group.
- `anotate_interactors.py` retrieves UniProt annotations for partners and counts their InterPro domains.
- `pipelineGENE.py` and `pipelineLIST.py` are the pipeline entry-point scripts.

## Run the current analysis

Use Python 3.10 or newer. From the repository root, create an environment and install the required packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install requests pandas
```

Run the current group-based analysis:

```bash
python pipeline/pipelineGENE.py
```

The script reads existing interaction tables from `pipeline/intact_uniprot_results/KNOTTED/` and UniProt accession lists from `pipeline/uniprot_ids/`. It writes partner lists, retrieves partner annotations from UniProt, and counts InterPro domains. Network access is required for the annotation and domain metadata requests.

The current `pipelineGENE.py` run does **not** refresh the AlphaKnot, InterPro, or IntAct input data. Those download and preparation calls are commented out. It currently processes the `KNOTTED` results; the corresponding `UNKNOTTED` calls are also commented out.

## Results

The repository includes previously downloaded data and intermediate results in these directories:

- `alphaknot_results/` — AlphaKnot TSV files.
- `interpro_results/` — InterPro protein data.
- `uniprot_ids/` — UniProt accession lists used by the analysis.
- `intact_uniprot_results/` and `intact_gene_results/` — IntAct interaction tables.
- `combined_intact_results/` — merged and deduplicated interaction tables.
- `isolated_partners/` — partner lists and InterPro annotations and frequency summaries.

For each group, the main derived files in `isolated_partners/` are:

- `*_partners.csv` — partner UniProt accessions and, when available, publication identifiers.
- `*_interactor_IPR.csv` — partner protein names, organisms, InterPro accessions, and publication identifiers.
- `*_count_interacter_IPR.csv` — InterPro domain counts, names, descriptions, and associated publications.

The `KNOTTED` and `UNKNOTTED` directory labels describe how the input lists were formed in this project. `UNKNOTTED` contains InterPro accessions not found in the AlphaKnot list; absence from AlphaKnot alone does not establish that a protein has no structural knot. Results are provided as tables, and coverage varies between groups; there is no single consolidated report with biological interpretation.

## Current limitations

- `pipelineLIST.py` expects `data/data.tsv`, which is not currently present in the repository.
- Running the annotation stage again can append duplicate rows to existing `*_interactor_IPR.csv` files. Move or clear those generated annotation files before recalculating them.
- The standalone block at the bottom of `api_intactGENE.py` passes an unsupported `input_dir` argument. The function should be called with `output_dir` only, or that entry point should be corrected, before using the script directly.
