"""
Przechowuje nazwy genów oraz odpowiadające im identyfikatory InterPro (IPR),
aby punkty wejścia w katalogu `pipeline/` można było modyfikować bez ingerencji
w podstawowy pakiet `alphaknot/`.

Obsługujemy też aliasy nazw (synonimy genów), tak aby różne nazwy tej samej
grupy białek (np. nep1/emg1) były mapowane na jedno, kanoniczne oznaczenie.
"""

GENE_TO_IPRS: dict[str, list[str]] = {
    "trmD": [
        "IPR002649",  # tRNA_m1G_MeTrfase_TrmD
        "IPR016009",  # tRNA_MeTrfase_TRMD/TRM10-like
    ],
    "trm5": [
        "IPR056743",  # TRM5/TYW2-like, methyltransferase domain
    ],
    "nep1": [
        "IPR005304",
    ],
}

# Mapuje alternatywne nazwy genów na nazwę kanoniczną używaną w GENE_TO_IPRS.
ALIASES: dict[str, str] = {
    # nep1 jest również opisywany jako EMG1
    "emg1": "nep1",
}


def normalize_gene_name(name: str) -> str:
    """
    Zwraca kanoniczną nazwę genu używaną w GENE_TO_IPRS,
    uwzględniając zdefiniowane aliasy (synonimy).
    """
    key = name.lower()
    return ALIASES.get(key, key)


def get_all_names_for_group(name: str) -> list[str]:
    """
    Zwraca listę wszystkich nazw (kanoniczna + aliasy), które odnoszą się
    do tej samej grupy białek.
    """
    canonical = normalize_gene_name(name)
    names = [canonical]
    # wszystkie aliasy, które wskazują na ten kanoniczny gen
    for alias, target in ALIASES.items():
        if target == canonical and alias not in names:
            names.append(alias)
    return names

