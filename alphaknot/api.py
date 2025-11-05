# DOWNLOAD ONLY THE TrmD ENTRIES FROM THE ALPHAKNOT DATABASE
import requests
import pandas as pd
from io import StringIO

url = "http://alphaknot.cent.uw.edu.pl/browse/?cats=&v=&adv=TrmD&organisms=&array=0&raw=1"

response = requests.get(url)
response.raise_for_status()

df = pd.read_csv(StringIO(response.text), sep="\t", engine="python", on_bad_lines="skip")

df.to_csv("alphaknot_trmD.tsv", sep="\t", index=False)


# DOWNLOAD THE WHOLE ALPHAKNOT DATABASE
# import requests
# import gzip
# import shutil
# import os

# # URL for the full AlphaKnot dataset (adjust if needed)
# # Example from their API: all structures
# ALPHAKNOT_URL = "https://alphaknot.cent.uw.edu.pl/all.txt.gz"

# # Local paths
# LOCAL_GZ = "alphaknot_all.txt.gz"
# LOCAL_TSV = "alphaknot_all.tsv"

# def download_database(url, local_gz_path):
#     print(f"Downloading database from {url}...")
#     response = requests.get(url, stream=True)
#     response.raise_for_status()
    
#     with open(local_gz_path, "wb") as f:
#         for chunk in response.iter_content(chunk_size=8192):
#             f.write(chunk)
#     print(f"Downloaded to {local_gz_path}")

# def decompress_gzip(gz_path, out_path):
#     print(f"Decompressing {gz_path} to {out_path}...")
#     with gzip.open(gz_path, 'rb') as f_in:
#         with open(out_path, 'wb') as f_out:
#             shutil.copyfileobj(f_in, f_out)
#     print(f"Decompressed to {out_path}")

# def main():
#     # Download the gzipped file
#     download_database(ALPHAKNOT_URL, LOCAL_GZ)
    
#     # Decompress to TSV
#     decompress_gzip(LOCAL_GZ, LOCAL_TSV)
    
#     print("Database is ready to use.")

# if __name__ == "__main__":
#     main()
