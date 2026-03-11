#jest takie cos jak batch search na intact, w ktorym mozna podac liste uniprot id i uzyskac wyniki dla kazdego z nich
# https://www.ebi.ac.uk/intact/home#batch-search
# moze pokazywac max 1500 interakcji na raz

def get_uniprot_ids_from_intact(name):
    """
    Returns a set of Uniprot IDs from all intAct .tsv files in the given pattern.
    Assumes the 'Uniprot' column is present in each file (as 3rd column, 0-based index 2).
    Returns:
        Set of unique Uniprot IDs found across files.
    """
    uniprot_ids = set()
    with open(name, "r", encoding="utf-8") as fin:
        print(f"Processing file: {name}")
        for line in fin:
            if line.startswith("#") or not line.strip():
                continue  # skip comments/blank lines
            cols = line.rstrip("\n").split("\t")
            # Only add IDs from 1st and 2nd column if they start with "uniprotkb:", and strip that prefix
            for i in range(2):
                if len(cols) > i and cols[i].startswith("uniprotkb:"):
                    uniprot_ids.add(cols[i][len("uniprotkb:"):])
    return uniprot_ids

if __name__ == "__main__":
    ids = get_uniprot_ids_from_intact('intAct_TrmD.tsv')
    with open("uniprot_ids.txt", "w") as fout:
        for uid in ids:
            fout.write(f"{uid}\n")
        # fout.write("P02829\n") #to wiem, ze na pewno jest w intact
