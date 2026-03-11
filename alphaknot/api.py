# DOWNLOAD ONLY THE TrmD
import requests
import pandas as pd
from io import StringIO

gene = "TrmD"
url = (
        f"https://alphaknot.cent.uw.edu.pl/browse?"
        f"field=InterPro&val=IPR002649"
        f"&conj=OR&val=IPR016009"
        f"&conj=OR&val=IPR019230"
        f"&conj=OR&val=IPR023148"
        f"&conj=OR&val=IPR029026"
        f"&conj=OR&val=IPR029028"
        f"&raw=2"
    )
resp = requests.get(url)
resp.raise_for_status()
print(url)
df = pd.read_csv(StringIO(resp.text), sep="\t", on_bad_lines='skip', engine='python')
df = df[df['Category'] != 'ESM1']
df.to_csv(f"alphaknot_all_AF2.tsv", sep="\t", index=False)

