import csv
from collections import defaultdict
gene='trm5'

def group_uniprot_by_interpro(tsv_path):
    interpro_to_uniprots = defaultdict(set)
    with open(tsv_path, newline='') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            if not row:
                continue
            if row[0].startswith("#"):
                continue
            if len(row) < 7:
                continue
            uniprot_id = row[2]
            interpro_field = row[6]
            if not uniprot_id or not interpro_field:
                continue
            interpro_ids = [id_.strip() for id_ in interpro_field.split(",") if id_.strip()]
            for interpro_id in interpro_ids:
                interpro_to_uniprots[interpro_id].add(uniprot_id)
    # Prepare a dict with counts for each interpro
    interpro_to_counts = {interpro: len(uniprots) for interpro, uniprots in interpro_to_uniprots.items()}
    return interpro_to_uniprots, interpro_to_counts

def main():
    tsv_path = f"{gene}.tsv"
    interpro_to_uniprots, interpro_to_counts = group_uniprot_by_interpro(tsv_path)

    with open(f"interpro_{gene}_uniprots.txt", "w") as outfile:
        for interpro in sorted(interpro_to_uniprots.keys()):
            uniprots = interpro_to_uniprots[interpro]
            count = interpro_to_counts[interpro]
            outfile.write(f"{interpro} ({count}): {', '.join(sorted(uniprots))}\n")

if __name__ == "__main__":
    main()
