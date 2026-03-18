import sys
from pathlib import Path
import csv
import requests
from collections import Counter, defaultdict

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))


def get_uniprot_data(uniprot_id):
    headers = {
        "accept": "application/json"
    }
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    response = requests.get(url, headers=headers)
    if not response.ok:
        response.raise_for_status()
        sys.exit()
    data = response.json()
    return data

def extract_interpro_data(data):
    """
    Extract InterPro cross-references, protein name, organism, and publications from UniProt JSON response.
    Returns (nazwa, organizm, interpro_list, publications).
    """
    interpro_list = []
    for xref in data.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "InterPro":
            interpro_list.append(xref.get("id"))

    nazwa = data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'Brak nazwy')
    organizm = data.get('organism', {}).get('scientificName', 'Brak organizmu')
    return nazwa, organizm, interpro_list


def write_interpro_to_file(output_file, uid, nazwa, organizm, interpro_list, publications=None):
    """
    Write one line of InterPro info for a given UniProt id to the output file.
    """
    with open(output_file, "a", encoding="utf-8") as outf:
        interpro_list = '|'.join(interpro_list) if interpro_list else 'Brak InterPro kodów dla tego białka.'
        if publications:
            outf.write(f"{uid},{nazwa},{organizm},{interpro_list},{publications}\n")
        else:
            outf.write(f"{uid},{nazwa},{organizm},{interpro_list}\n")
def write_failure_to_file(output_file, uid):
    """
    Record a failure to the output file.
    """
    with open(output_file, "a", encoding="utf-8") as outf:
        outf.write(f'{uid}\n')

def uniprot_ids_to_interpro(csv_file, output_file):
    """
    Process all UniProt IDs from a CSV and write InterPro and metadata results to output_file.
    """
    uniprot_ids = {}

    with open(csv_file, encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader, None)  # safely skip header

        for row in reader:
            if not row:
                continue

            row = [x.strip() for x in row]

            if len(row) > 1 and row[0] and row[1]:
                uniprot_ids[row[0]] = row[1]
            elif row[0]:
                uniprot_ids[row[0]] = None


    for uid in uniprot_ids:
        if not uid :
            continue
        
        print(f"\nDane dla {uid}:")
        try:
            data = get_uniprot_data(uid)
            if data:
                nazwa, organizm, interpro_list = extract_interpro_data(data)
                print(uniprot_ids.get(uid))
                write_interpro_to_file(output_file, uid, nazwa, organizm, interpro_list,uniprot_ids.get(uid))
            else:
                print("Nie udało się pobrać danych z UniProt.")
        except Exception:
            write_failure_to_file(output_file, uid)


def process_directory_of_reactors(input_dir):
    """
    Processes all reactors_*.csv files in the given directory,
    writing output interpro_kody_{protein}.csv for each.
    """
    input_dir = Path(input_dir)
    for reactor_file in input_dir.glob("*_partners.csv"):
        protein = reactor_file.stem.replace("_partners", "")
        output_file = input_dir / f"{protein}_interactor_IPR.csv"
        uniprot_ids_to_interpro(reactor_file, output_file)

def annotate_pipeline(uniprot_partner_dir, output_dir=None):
    process_directory_of_reactors(uniprot_partner_dir)

# zliczanie wystepujacych domen interpro i sortowanie po czestosci
def get_interpro_data(interpro_id):
    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/{interpro_id}/"
    try:
        r = requests.get(url, headers={"Accept": "application/json"})
        if r.status_code == 200:
            data = r.json()
            name, short_desc = "", ""
            if "metadata" in data:
                md = data["metadata"]
                if isinstance(md.get("name"), dict):
                    name_val = md["name"].get("name", "")
                    name_short = md["name"].get("short", "")
                    name = name_val if name_val else name_short
                elif isinstance(md.get("name"), str):
                    name = md["name"]
                desc = md.get("description", "")
                short_desc = desc if desc else ""
            else:
                name = data.get("name", "")
                short_desc = data.get("abstract") or data.get("description") or ""

            if isinstance(short_desc, list):
                short_desc = "; ".join(item.get('text', str(item)) if isinstance(item, dict) else str(item) for item in short_desc)
            elif not isinstance(short_desc, str):
                short_desc = str(short_desc if short_desc is not None else "")
            return name, short_desc
        else:
            return "", ""
    except Exception as e:
        return "", ""

def count_interpro(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for reactor_file in input_dir.glob("*_interactor_IPR.csv"):
        protein = reactor_file.stem.replace("_interactor_IPR", "")
        output_file = output_dir / f"{protein}_count_interacter_IPR.csv"
        print(f"Processing file: {reactor_file} for protein: {protein}")
        
        interpro_counts = Counter()
        interpro_pubmeds = defaultdict(set)

        with open(reactor_file, encoding='utf-8', newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) < 4:
                    continue
                interpro_raw = row[3].strip()
                pubmed_raw = row[4].strip()

                if not interpro_raw or 'Brak InterPro kodów dla tego białka.' in interpro_raw:
                    continue
                codes = {c.strip() for c in interpro_raw.split('|') if c.strip()}

                pubmed_ids = set()
                if pubmed_raw and pubmed_raw.lower() != 'brak':
                    pubmed_ids = {p.strip() for p in pubmed_raw.split('|') if p.strip()}
                
                for code in codes:
                    interpro_counts[code]+=1
                    interpro_pubmeds[code].update(pubmed_ids)
        print(f"Writing InterPro domain counts to: {output_file}")
        with open(output_file, 'w', encoding='utf-8', newline='') as out:
            writer = csv.writer(out, delimiter=';')
            writer.writerow(["Count", "InterProID", "Name", "ShortDescription",'pubmeds'])
            for code, count in interpro_counts.most_common():
                print(f"Writing InterPro code: {code} ({count} occurrences)")
                name, short_desc = get_interpro_data(code)
                writer.writerow([count, code,"\n".join(sorted(interpro_pubmeds[code])) ,name, short_desc])
    print("Finished counting InterPro domains for all processed files.")


if __name__ == "__main__":
    input_dir = parent_dir / 'pipeline' / 'isolated_partners'
    # process_directory_of_reactors(input_dir)
    count_interpro(input_dir, input_dir)

