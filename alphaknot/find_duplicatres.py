# Script to compare two files and find unique rows based on the 2nd and 3rd columns,
# and save these unique entries into separate files.

def get_entries_from_file1(filename):
    """In file1, the relevant data is in the second column (index 1)"""
    s = set()
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            s.add(parts[1])
    return s

def get_entries_from_file2(filename):
    """In file2, the relevant data is in the third column (index 2)"""
    s = set()
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            s.add(parts[2])
    return s

def write_unique_to_file(unique_set, outfilename):
    with open(outfilename, "w", encoding="utf-8") as out:
        for val in sorted(unique_set):
            out.write(f"{val}\n")

def find_and_save_unique(file1, file2, out1="unique_in_1.tsv", out2="unique_in_2.tsv"):
    set1 = get_entries_from_file1(file1)
    set2 = get_entries_from_file2(file2)
    unique1 = set1 - set2
    unique2 = set2 - set1

    print(f"Saving unique (file1 not in file2) {file1} entries to {out1}")
    write_unique_to_file(unique1, out1)

    print(f"Saving unique (file2 not in file1) {file2} entries to {out2}")
    write_unique_to_file(unique2, out2)

if __name__ == "__main__":
    find_and_save_unique('38tys_trmd.tsv', 'web_ak_data_trmD.tsv')
    # Script to download alphaknot data by Uniprot IDs listed in a file

    import requests

    def load_uniprot_ids(filename):
        """Read Uniprot IDs from file (1 per line, skip lines starting with #)"""
        ids = []
        with open(filename, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # Remove -F1 at the end, if present
                if line.endswith("-F1"):
                    line = line[:-3]
                ids.append(line)
        print(len(ids))
        return ids

    def download_alphaknot_data(ids, outfile="alphaknot_data.tsv"):
        """Download alphaknot data for a list of uniprot IDs and save to outfile."""
        if not ids:
            print("No Uniprot IDs given.")
            return
        ids_str = ";".join(ids)
        query_params = [
            f"field=Uniprot&val={ids_str}",
            "conj=AND&field=Category&val=AF4",
            "raw=2",
            "result_cols=Knot_type;Category;Uniprot;Organism;pLDDT_knotcore;Protein_name;InterPro;Gene_name"
        ]
        base_url = "https://alphaknot.cent.uw.edu.pl/browse/?" + "&".join(query_params)
        params = None
        try:
            # print(base_url)
            resp = requests.get(base_url, params=params)
            resp.raise_for_status()
            with open(outfile, "w", encoding="utf-8") as f:
                f.write(resp.text)
            print(f"Saved to {outfile}")
        except Exception as e:
            print(f"Failed to download or save alphaknot data: {e}")

    
    # Process both unique_in_2.tsv and unique_in_1.tsv files to download alphaknot data
    files = [
        ("unique_in_2.tsv", "alphaknot_data_from_2.tsv"),
        ("unique_in_1.tsv", "alphaknot_data_from_1.tsv"),
    ]

    for in_fname, out_fname in files:
        uniprot_ids = load_uniprot_ids(in_fname)
        if uniprot_ids:
            print(f"Downloading data for {len(uniprot_ids)} Uniprot IDs from {in_fname}...")
            # Download in batches of 50 and append to the same output file
            batch_size = 50
            # Write header to output file first
            wrote_header = False
            with open(out_fname, "w", encoding="utf-8") as outf:
                for batch_start in range(0, len(uniprot_ids), batch_size):
                    batch_ids = uniprot_ids[batch_start:batch_start + batch_size]
                    # Download and save to a temp file
                    tmp_outfile = "tmp_alphaknot_batch.tsv"
                    download_alphaknot_data(batch_ids, tmp_outfile)
                    # Read in data, write header only once
                    with open(tmp_outfile, "r", encoding="utf-8") as tf:
                        lines = tf.readlines()
                        if not lines:
                            continue
                        # Write header if not yet written, skip subsequent headers
                        header = lines[0]
                        if not wrote_header:
                            outf.write(header)
                            wrote_header = True
                        # Write non-header lines
                        outf.writelines(lines[1:])
            print(f"Saved all batches to {out_fname}")
        else:
            print(f"No Uniprot IDs found in {in_fname}.")
