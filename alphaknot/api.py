import requests
from pathlib import Path
from typing import Iterable

# define the iprs to query
# znalazlam przez uniprot wszystkie iprs dla danego bialka
trmD_iprs = [
    "IPR002649", #tRNA_m1G_MeTrfase_TrmD
    "IPR016009", #tRNA_MeTrfase_TRMD/TRM10-like
    # "IPR019230", #RNA_MeTrfase_C_dom — the C-terminal domain of SPOUT (SpoU/TrmD) methyltransferases
    # "IPR023148", #tRNA_m1G_MeTrfase_C_sf — structural family of the C-terminal region of TrmD-like tRNA (guanine-N1)-methyltransferases
    # "IPR029026", #tRNA (guanine-N1)-methyltransferase, N-terminal domain (in TrmD-like enzymes)
    # "IPR029028" #Alpha/beta knot methyltransferases — the “knot” fold in many RNA methyltransferases (SPOUT family).
]
trm5_iprs = [
    # "IPR030382", #SAM-dependent methyltransferase TRM5 / TYW2-type domain.
    # "IPR029063", #S-adenosyl-L-methionine (SAM)-dependent methyltransferase superfamily (general class of MTases).
    "IPR056743", #TRM5/TYW2-like, methyltransferase domain
    # "IPR056744",
    # "IPR025792", #tRNA (guanine(37)-N(1))-methyltransferase, eukaryotic (TRM5-type)
]
nep1_iprs = [
    # "IPR029028", #Alpha/beta knot methyltransferases — the “knot” fold in many RNA methyltransferases (SPOUT family).
    "IPR005304",
    # "IPR029026", #tRNA (guanine-N1)-methyltransferase, N-terminal domain (in TrmD-like enzymes)
]
# define the iprs to query

# map gene names to their IPR lists so we can iterate automatically
GENE_TO_IPRS = {
    "trmD": trmD_iprs,
    "trm5": trm5_iprs,
    "nep1": nep1_iprs,
}


def build_alphaknot_url(iprs, category: str = "AF4") -> str:
    iprs_str = ";".join(iprs)
    query_params = [
        f"field=InterPro&val={iprs_str}",
        f"conj=AND&field=Category&val={category}&raw=1&"
        "result_cols=Knot_type;Category;Uniprot;Organism;pLDDT_knotcore;Protein_name;InterPro;Gene_name",
    ]
    return "https://alphaknot.cent.uw.edu.pl/browse/?" + "&".join(query_params)


def fetch_alphaknot_for_gene(
    gene_name: str,
    iprs: Iterable[str],
    out_dir: str | Path = ".",
    category: str = "AF4",
) -> Path:
    """
    Zapytaj AlphaKnot dla podanego genu i listy IPR i zapisz jako plik TSV
    Zwraca sciezke do zapisanego pliku TSV
    """
    url = build_alphaknot_url(iprs, category=category)
    session = requests.Session() #wykonywanie zapytań sieciowych
    session.trust_env = False #ignorowanie zmiennych srodowiskowych proxy
    resp = session.get(url, timeout=60)
    resp.raise_for_status() #wyjatek jesli wystapi blad

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    tsv_path = out_path / f"{gene_name}.tsv"
    tsv_path.write_text(resp.text, encoding="utf-8")
    if tsv_path.exists():
        print(f"{gene_name}: TSV file already exists at {tsv_path}, skipping fetch.")
        return tsv_path

    lines = [ln for ln in resp.text.splitlines() if ln.strip() != ""] #liczenie wierszy 
    row_count = max(0, len(lines) - 1)
    print(f"{gene_name}: downloaded {row_count} rows from {url}")
    print(f"Saved to {tsv_path}")

    return tsv_path


def fetch_all_default_genes(out_dir: str | Path = "alphaknot_results") -> dict[str, Path]:
    """
    Pobiera dane dla wszystkich znanych genów w GENE_TO_IPRS.
    Zwraca słownik gene_name -> sciezka do zapisanego pliku TSV.
    """
    results: dict[str, Path] = {}
    for gene_name, iprs in GENE_TO_IPRS.items():
        results[gene_name] = fetch_alphaknot_for_gene(gene_name, iprs, out_dir=out_dir)
    return results


if __name__ == "__main__":
    # Uruchomienie tego pliku bezpośrednio pobierze i zapisze pliki TSV dla wszystkich genów
    fetch_all_default_genes()
