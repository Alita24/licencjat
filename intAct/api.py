import requests
import pandas as pd

url = "https://www.ebi.ac.uk/intact/ws/interaction/findInteractions/trmD"
resp = requests.get(url)
resp.raise_for_status()

data = resp.json()
interactions = data.get("content", [])

df = pd.json_normalize(interactions)

df.to_csv("intact_trmD_interactions.tsv", sep="\t", index=False)
