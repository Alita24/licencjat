# Protein Interaction and InterPro Domain Analysis

This repository contains scripts and data for analyzing selected protein families. It covers identifying proteins and their interactions in AlphaKnot, InterPro, and IntAct, then finding the InterPro domains most frequently observed among their interaction partners.

## Repository contents

- `data/targets.py` — configuration for 20 protein groups and their associated InterPro (IPR) identifiers. Edit this file to add or change the groups being analyzed.
- `alphaknot/api.py`, `interpro/api.py`, `intAct/api.py` — clients for AlphaKnot, InterPro, and IntAct.
- `pipeline/` — scripts for downloading, processing, and annotating data, along with intermediate and final results.
- `data/63.tsv`, `data/list.tsv`, `data/solenoid.tsv` — additional data tables. `list.tsv` is currently empty - edit if explicitly providing UniProt identifiers input modality.

## Analysis workflow

1. Protein groups and their IPR identifiers are read from `data/targets.py`.
2. AlphaKnot and InterPro provide lists of proteins associated with those groups.
3. UniProt identifiers are collected and used to retrieve interaction data from IntAct. The scripts can also query IntAct by protein group name.
4. Interaction data is combined, and the UniProt identifiers of interaction partners are separated from those of the proteins in the target group.
5. Partners are annotated using UniProt data. The scripts count their InterPro domains and retrieve domain names and descriptions.

The `KNOTTED` and `UNKNOTTED` directory labels reflect how lists were assembled in this project: `KNOTTED` contains identifiers obtained from AlphaKnot, while `UNKNOTTED` contains InterPro identifiers that were not found in the AlphaKnot list. A protein's absence from AlphaKnot does not, by itself, prove that its structure lacks a knot.

## Setup and run

Python 3.10 or newer is required, along with `requests` and `pandas`. From the repository root, create a virtual environment and install the dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install requests pandas
```

Run the currently active analysis stage with:

```bash
python pipeline/pipelineGENE.py
```

This script uses previously downloaded interaction data in `pipeline/intact_uniprot_results/KNOTTED/`, creates partner lists in `pipeline/isolated_partners/KNOTTED/`, retrieves partner annotations from UniProt, and writes InterPro domain counts. An internet connection is required to retrieve annotations.

**This command does not refresh the full dataset.** The AlphaKnot, InterPro, and IntAct download steps in `pipeline/pipelineGENE.py` are commented out. Scripts for individual stages are in `pipeline/`; before refreshing the data, check their configuration and input/output paths. The experimental `pipeline/pipelineLIST.py` script expects `data/data.tsv`, which is not currently included in the repository.

When rerunning annotation, move or remove existing `*_interactor_IPR.csv` files from the output directory first. Annotations are appended to these files, so rerunning can duplicate rows.

## Results

The repository already includes downloaded data and results for some stages. Key output locations include:

- `pipeline/alphaknot_results/` — AlphaKnot query results in TSV format.
- `pipeline/interpro_results/` — downloaded and combined InterPro data.
- `pipeline/intact_uniprot_results/` and `pipeline/intact_gene_results/` — IntAct interaction tables.
- `pipeline/isolated_partners/` — partner lists (`*_partners.csv`), UniProt/InterPro annotations (`*_interactor_IPR.csv`), and domain frequency summaries (`*_count_interacter_IPR.csv`).

The final outputs are per-group tables listing identified interaction partners, publication identifiers, and InterPro domains ranked by frequency. The available results vary by group; the repository does not include a single consolidated report with biological interpretation.

## Key scripts

- `pipeline/api_alphaknot.py` — downloads AlphaKnot data for the groups in `data/targets.py`.
- `pipeline/api_interpro.py` — downloads proteins associated with the configured IPR identifiers.
- `pipeline/isolate_uniprot_ids.py` — extracts, combines, and compares UniProt identifier lists.
- `pipeline/api_intactUNI.py` and `pipeline/api_intactGENE.py` — retrieve IntAct interactions by UniProt identifier or group name.
- `pipeline/combine_intact.py` — combines and deduplicates interaction tables.
- `pipeline/isolate_interactors.py` — extracts interaction partners.
- `pipeline/anotate_interactors.py` — retrieves UniProt annotations and counts InterPro domains.
