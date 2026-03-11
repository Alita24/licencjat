# ------------------------------------------------------------
# Ten skrypt:
# 1. izoluje  wszystkie unikalne identyfikatory UniProt (uniprotkb:XXX) z pierwszych dwóch kolumn pliku TSV z wynikami IntAct, 
#    generując plik tekstowy z listą tych identyfikatorów.
# 2. Dla każdego z tych uniprot ID pobiera odpowiadające im wpisy InterPro (ID, nazwa i organizm) z API UniProt,
#    a następnie zapisuje te powiązania do pliku CSV: UniProt_ID, InterPro_ID, InterPro_Name, organism.
# Użycie:
#   - Zmień nazwę genu w zmiennej gene.
#   - Umieść plik "intAct_{gene}.tsv" w tym katalogu.
#   - Uruchom ten skrypt; wygeneruje "uniprotkb_ids_{gene}.txt" i "uniprot_to_interpro_{gene}.csv".
# ------------------------------------------------------------

import re
import requests
import csv
import time

gene = 'nep1'

def extract_uniprotkb_ids_from_tsv_first_two_cols(tsv_path):
    """
    Otwiera plik TSV i wyciąga wszystkie unikalne identyfikatory uniprotkb:xxx,
    ale tylko z 1 lub 2 kolumny (indeksowane od zera).
    Zwraca posortowaną listę.
    """
    uniprotkb_pattern = re.compile(r'uniprotkb:([A-Z0-9]+)')
    ids = set()
    with open(tsv_path, encoding="utf-8") as f:
        for line in f:
            # Pomijaj wiersz nagłówka jeśli istnieje
            cols = line.rstrip('\n').split('\t')
            # Działaj tylko jeśli są przynajmniej 2 kolumny
            for col in cols[:2]:
                matches = uniprotkb_pattern.findall(col)
                for match in matches:
                    ids.add(match)
    sorted(ids)
    with open(f"uniprotkb_ids_{gene}.txt", "w", encoding="utf-8") as out_f:
        for id_ in ids:
            out_f.write(f"{id_}\n")



sciezka_do_tsv = f"intAct_{gene}.tsv"
extract_uniprotkb_ids_from_tsv_first_two_cols(sciezka_do_tsv)


INPUT_FILE = f"uniprotkb_ids_{gene}.txt"      # <-- your file with one UniProt ID per line
OUTPUT_FILE = f"uniprot_to_interpro_{gene}.csv"


def get_interpro(uniprot_id):
    """Return InterPro IDs and names for a UniProt ID."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"

    try:
        r = requests.get(url, timeout=20)
        if r.status_code != 200:
            return []

        data = r.json()
        interpro_entries = []

        # Also capture organism name from the UniProt record
        organism = None
        # Try to get the 'organism' field (fall back to scientificName)
        if "organism" in data and isinstance(data["organism"], dict):
            organism = data["organism"].get("scientificName")
        for xref in data.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "InterPro":
                interpro_id = xref.get("id")
                name = None
                props = xref.get("properties", [])
                if props:
                    name = props[0].get("value")
                interpro_entries.append((interpro_id, name, organism))

        return interpro_entries

    except Exception as e:
        print(f"Error with {uniprot_id}: {e}")
        return []


def main():
    with open(INPUT_FILE, "r") as f:
        uniprot_ids = [line.strip() for line in f if line.strip()]

    rows = []

    for uid in uniprot_ids:
        print(f"Fetching: {uid}")
        interpro_list = get_interpro(uid)

        if not interpro_list:
            rows.append([uid, None, None])
        else:
            for ipr, name, orga in interpro_list:
                rows.append([uid, ipr, name, orga])

        time.sleep(0.2)  # be polite to the server

    # write output CSV
    with open(OUTPUT_FILE, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["UniProt_ID", "InterPro_ID", "InterPro_Name", 'organism'])
        writer.writerows(rows)

    print(f"\nDone! Results written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
