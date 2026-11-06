# Analiza interakcji białek i domen InterPro

Repozytorium zawiera skrypty i dane do analizy białek należących do wybranych rodzin: od wyszukania białek i ich interakcji w bazach AlphaKnot, InterPro i IntAct po wskazanie domen InterPro najczęściej spotykanych u partnerów interakcji.

## Co znajduje się w repozytorium

- `data/targets.py` — konfiguracja 20 grup białek i powiązanych z nimi identyfikatorów InterPro (IPR). To tutaj dodaje się lub zmienia analizowane grupy.
- `alphaknot/api.py`, `interpro/api.py`, `intAct/api.py` — funkcje komunikujące się z bazami AlphaKnot, InterPro i IntAct.
- `pipeline/` — skrypty łączące etapy pobierania, przetwarzania i anotacji oraz zapisane wyniki pośrednie i końcowe.
- `data/63.tsv`, `data/list.tsv`, `data/solenoid.tsv` — dodatkowe tabele danych; `list.tsv` jest obecnie pusty.

## Jak działa analiza

1. Grupy białek i ich identyfikatory IPR są odczytywane z `data/targets.py`.
2. AlphaKnot i InterPro dostarczają listy białek powiązanych z tymi grupami.
3. Identyfikatory UniProt są zbierane i wykorzystywane do pobrania interakcji z IntAct. Skrypty mogą też pobierać interakcje po nazwie grupy.
4. Dane interakcji są łączone, a identyfikatory UniProt partnerów są oddzielane od identyfikatorów białek z analizowanej grupy.
5. Partnerzy są anotowani na podstawie danych UniProt. Skrypt zlicza występowanie przypisanych im domen InterPro i zapisuje ich nazwy oraz opisy.

W katalogach `KNOTTED` i `UNKNOTTED` etykiety odzwierciedlają sposób przygotowania list w tym projekcie: `KNOTTED` zawiera identyfikatory uzyskane z AlphaKnot, a `UNKNOTTED` — identyfikatory z InterPro, których nie znaleziono na liście AlphaKnot. Brak białka w AlphaKnot sam w sobie nie potwierdza, że jego struktura nie zawiera węzła.

## Uruchomienie

Wymagany jest Python 3.10 lub nowszy oraz pakiety `requests` i `pandas`. Z katalogu głównego repozytorium utwórz środowisko i zainstaluj zależności:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install requests pandas
```

Następnie uruchom aktualny etap analizy:

```bash
python pipeline/pipelineGENE.py
```

Skrypt korzysta z zapisanych wcześniej interakcji w `pipeline/intact_uniprot_results/KNOTTED/`, tworzy listy partnerów w `pipeline/isolated_partners/KNOTTED/`, pobiera ich anotacje z UniProt i zapisuje zliczenia domen InterPro. Do pobrania anotacji potrzebne jest połączenie z internetem.

**To polecenie nie odświeża całego zbioru danych.** Kroki pobierania danych z AlphaKnot, InterPro i IntAct są w `pipeline/pipelineGENE.py` zakomentowane. Skrypty do poszczególnych etapów znajdują się w `pipeline/`; przed pełnym odświeżeniem należy sprawdzić ich konfigurację i wskazane katalogi wejściowe oraz wyjściowe. Eksperymentalny skrypt `pipeline/pipelineLIST.py` oczekuje pliku `data/data.tsv`, którego obecnie nie ma w repozytorium.

Przy ponownym uruchomieniu anotacji warto usunąć lub przenieść istniejące pliki `*_interactor_IPR.csv` z katalogu wynikowego: anotacje są dopisywane do tych plików, więc ponowne wykonanie może powielić wiersze.

## Wyniki

W repozytorium znajdują się już dane pobrane i wyniki dla części etapów. Najważniejsze pliki wynikowe to:

- `pipeline/alphaknot_results/` — wyniki zapytań do AlphaKnot w formacie TSV.
- `pipeline/interpro_results/` — pobrane i połączone dane InterPro.
- `pipeline/intact_uniprot_results/` i `pipeline/intact_gene_results/` — tabele interakcji z IntAct.
- `pipeline/isolated_partners/` — listy partnerów (`*_partners.csv`), anotacje UniProt/InterPro (`*_interactor_IPR.csv`) oraz podsumowania częstości domen (`*_count_interacter_IPR.csv`).

Końcowym rezultatem są tabele dla poszczególnych grup białek: lista wykrytych partnerów, ich identyfikatory publikacji oraz ranking domen InterPro według częstości. Zakres dostępnych wyników różni się między grupami; repozytorium nie zawiera jednego zbiorczego raportu z interpretacją biologiczną.

## Najważniejsze skrypty

- `pipeline/api_alphaknot.py` — pobieranie danych z AlphaKnot dla grup z `data/targets.py`.
- `pipeline/api_interpro.py` — pobieranie białek przypisanych do skonfigurowanych identyfikatorów IPR.
- `pipeline/isolate_uniprot_ids.py` — wyodrębnianie, łączenie i porównywanie list UniProt.
- `pipeline/api_intactUNI.py` i `pipeline/api_intactGENE.py` — pobieranie interakcji z IntAct po identyfikatorach UniProt lub nazwie grupy.
- `pipeline/combine_intact.py` — łączenie i deduplikacja tabel interakcji.
- `pipeline/isolate_interactors.py` — wyodrębnianie partnerów interakcji.
- `pipeline/anotate_interactors.py` — pobieranie anotacji UniProt i zliczanie domen InterPro.
