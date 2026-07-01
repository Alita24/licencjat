from pathlib import Path
import sys
import pandas as pd
import csv

parent_dir = Path(__file__).resolve().parents[1]
if str(parent_dir) not in sys.path:
    sys.path.append(str(parent_dir))

from data.targets import GENE_TO_IPRS

def isolate_all_uniprot_ids(group_name, intact_dir, add_publication_col=True):
    intact_dir = Path(intact_dir)
    partner_ids = set()
    partner_pub_dict = {} if add_publication_col else None

    for tsv_file in intact_dir.glob(f"intact_{group_name}_*.tsv"):
        try:
            df = pd.read_csv(tsv_file, sep="\t", skiprows=1)

            pub_col = None
            if add_publication_col:
                for cname in df.columns:
                    if "pubmed" in cname.lower() or "publication" in cname.lower():
                        pub_col = cname
                        break
                if pub_col is None:
                    pub_col = df.columns[-1]

            for _, row in df.iterrows():
                id1, id2 = row.iloc[0], row.iloc[1]
                if pd.isna(id1) or pd.isna(id2):
                    continue
                id1, id2 = str(id1).strip(), str(id2).strip()
                if id1.startswith('uniprotkb'):
                    id1 = str(id1).split(':', 1)[1]
                if id2.startswith('uniprotkb'):
                    id2 = str(id2).split(':', 1)[1]

                if id1.startswith('intact:') and id2.startswith('intact:'):
                    continue

                partner_ids.update([id1, id2])

                if add_publication_col and pub_col:
                    pubs = row.get(pub_col)
                    for pid in [id1, id2]:
                        partner_pub_dict.setdefault(pid, set())
                        if pd.notna(pubs):
                            if isinstance(pubs, str):
                                for pub in pubs.split("|"):
                                    pub = pub.strip()
                                    if pub:
                                        partner_pub_dict[pid].add(pub)
                            else:
                                partner_pub_dict[pid].add(str(pubs).strip())
        except Exception as e:
            print(f"[ISOLATE INTERACTORS] Failed to read {tsv_file}: {e}")
    if add_publication_col:
        partner_pub_dict = {k: v for k, v in partner_pub_dict.items() if len(v)}
        return partner_ids, partner_pub_dict
    else:
        return partner_ids

def write_partner_ids_for_group(group_name, intact_dir, out_path, add_publication_col=True):
    intact_dir = Path(intact_dir)
    out_path = Path(out_path)
    out_path.mkdir(parents=True, exist_ok=True)

    out_file = out_path / f"{group_name}_partners.csv"

    if add_publication_col:
        all_uniprot_ids, pub_dict = isolate_all_uniprot_ids(group_name, intact_dir, add_publication_col=True)
    else:
        all_uniprot_ids = isolate_all_uniprot_ids(group_name, intact_dir)
        pub_dict = None

    uniprot_ids_dir = parent_dir / 'pipeline' / "uniprot_ids"

    group_uniprot_ids = set()
    for f in uniprot_ids_dir.glob(f"{group_name}_uniprot_ids_*.txt"):
        with open(f, "r") as fin:
            group_uniprot_ids.update(line.strip() for line in fin if line.strip())

    partners = all_uniprot_ids - group_uniprot_ids

    with open(out_file, "w", newline='') as f:
        writer = csv.writer(f)
        headers = ["PartnerUniProtID"]
        if add_publication_col:
            headers.append("Publications")
        writer.writerow(headers)

        for pid in sorted(partners):
            row = [pid]
            if add_publication_col:
                pubs = pub_dict.get(pid, set()) if pub_dict else set()
                if isinstance(pubs, set):
                    pubs_joined = "|".join(sorted(pubs))
                else:
                    pubs_joined = str(pubs)
                row.append(pubs_joined)
            writer.writerow(row)

    print(f"[ISOLATE INTERACTORS] {group_name}: written to: {out_file}")

def isolate_partners(intact_dir, output_dir, add_publication=False):
    for fam in GENE_TO_IPRS.keys():
        if '/' in fam:
            fam = fam.replace('/', '-')
        write_partner_ids_for_group(fam, intact_dir=intact_dir, out_path=output_dir, add_publication_col=add_publication)

if __name__ == "__main__":
    intact_gene_results = Path("pipeline/intact_uniprot_results")
    write_partner_ids_for_group('trmD', Path(intact_gene_results), Path('pipeline/isolated_partners'))
