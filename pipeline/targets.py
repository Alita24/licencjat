"""
Przechowuje nazwy genów oraz odpowiadające im identyfikatory InterPro (IPR),
aby punkty wejścia w katalogu `pipeline/` można było modyfikować bez ingerencji
w podstawowy pakiet `alphaknot/`.
"""

GENE_TO_IPRS: dict[str, list[str]] = {
    # "trmD": [
    #     "IPR002649",  # tRNA_m1G_MeTrfase_TrmD
    #     "IPR016009",  # tRNA_MeTrfase_TRMD/TRM10-like
    # ],
    "trm5": [
        "IPR056743",  # TRM5/TYW2-like, methyltransferase domain
    ]
    # "nep1": [
    #     "IPR005304",
    # ],
}


