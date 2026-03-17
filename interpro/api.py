import requests
import csv
import time

def get_reviewed_proteins_for_interpro(interpro_id):
    """
    Pobiera wszystkie białka powiązane z danym IPR
    Args:
        interpro_id: identyfikator IPR
    Returns:
        list of proteins
    """
    url = f"https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry/interpro/{interpro_id}"
    results = []
    params = {"page_size": 200}

    page = 0

    print(f"Downloading proteins for {interpro_id}...\n")

    while url:
        page += 1

        r = requests.get(url, params=params)
        if r.status_code != 200:
            raise Exception(f"HTTP {r.status_code}: {r.text}")

        d = r.json()
        batch = d.get("results", [])
        results.extend(batch)

        print(f"Page {page} | {len(batch)} proteins | total: {len(results)}")

        url = d.get("next")
        params = {}     # remove params after first request

    print(f"\nFinished. Total proteins: {len(results)}")
    return results


def extract_rows(protein):
    """
    wyciaga dane z jsona i zwraca liste slownikow
    Args:
        protein: json object
    Returns:
        list of dictionaries
    """

    metadata = protein.get("metadata", {})
    interpro_entries = protein.get("entries", [])

    rows = []

    for e in interpro_entries:
        for loc in e.get("entry_protein_locations", []):
            for frag in loc.get("fragments", []):
                rows.append({
                    "protein_accession": metadata.get("accession"),
                    "protein_name": metadata.get("name"),
                    "protein_length": metadata.get("length"),
                    "gene": metadata.get("gene"),
                    "organism_taxid": metadata.get("source_organism", {}).get("taxId"),
                    "organism_name": metadata.get("source_organism", {}).get("scientificName"),
                    "organism_fullname": metadata.get("source_organism", {}).get("fullName"),
                    "in_alphafold": metadata.get("in_alphafold"),

                    "interpro_accession": e.get("accession"),
                    "entry_type": e.get("entry_type"),
                    "entry_integrated": e.get("entry_integrated"),
                    "interpro_protein_length": e.get("protein_length"),

                    "fragment_start": frag.get("start"),
                    "fragment_end": frag.get("end"),
                    "fragment_status": frag.get("dc-status"),
                })

    return rows


def write_combined_csv(all_proteins, output_file="interpro_output.csv"):
    """
    zapisuje dane do csv
    Args:
        all_proteins: list of proteins
        output_file: ścieżka do pliku csv
    Returns:
        None
    """
    all_rows = []
        # list of dictionaries
        # each dictionary is a protein and has the following keys:
        # protein_accession, protein_name, protein_length, gene, organism_taxid, organism_name, organism_fullname, in_alphafold, interpro_accession, entry_type, entry_integrated, interpro_protein_length, fragment_start, fragment_end, fragment_status

    for protein in all_proteins:
        rows = extract_rows(protein)
        all_rows.extend(rows)
    


    if not all_rows:
        print("No annotation rows extracted.")
        return

    fieldnames = list(all_rows[0].keys())

    with open(output_file, "w", encoding="utf-8", newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(all_rows)

    print(f"\nCSV written → {output_file}")
    print(f"Total rows: {len(all_rows)}")


if __name__ == "__main__":
    INTERPRO_ID = "IPR002649"

    proteins = get_reviewed_proteins_for_interpro(INTERPRO_ID)

    write_combined_csv(proteins, f"{INTERPRO_ID}.csv")
