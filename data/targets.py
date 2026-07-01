"""
Przechowuje nazwy genów oraz odpowiadające im identyfikatory InterPro (IPR),
aby punkty wejścia w katalogu `pipeline/` można było modyfikować bez ingerencji
w podstawowy pakiet `alphaknot/`.
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
# }
    # SPOUT Superfamily Components
    "tRNA_guanine_N1_methyltransferase_N_terminal": [
        "IPR029026"
    ],
    "tRNA_methyltransferase_TRM10_type_domain_superfamily": [
        "IPR038459"
    ],
    
    # ITIH Superfamily Components
    "von_Willebrand_factor_type_A": [
        "IPR002035"
    ],
    "VIT_domain": [
        "IPR013694"
    ],
    
    # ATC/OTCase Superfamily Components
    "Aspartate_ornithine_carbamoyltransferase": [
        "IPR036901",
        "IPR006131",
    ],
    
    # Standalone Membrane / Exchanger Families
    "Sodium_calcium_exchanger_membrane": [
        "IPR004837"
    ],
    
    # AdoMet Synthase Superfamily Components
    "S_adenosylmethionine_synthetase": [
        "IPR022628",
        "IPR002133",
        "IPR022636",
    ],
    "S_adenosylmethionine_synthetase_archaea": [
        "IPR027790",
        "IPR042544",
        "IPR002795",
    ],
    
    # Carbonic Anhydrase Superfamily Components
    "Alpha_carbonic_anhydrase": [
        "IPR036398"
    ],
    "Delta_carbonic_anhydrase": [
        "IPR018883"
    ],
    
    # Transporter & Integrin Families
    "Oligopeptide_transporter": [
        "IPR004813"
    ],
    "Integrin_alpha": [
        "IPR032695",
        "IPR013649",
    ],
    
    # TDD Superfamily Components
    "tRNA_uridine_acp_transferase": [
        "IPR005636"
    ],
    "16S_18S_rRNA_acp_transferase_Tsr3": [
        "IPR007209",
        "IPR022968",
        "IPR007177",
    ],
    "Calcium_activated_potassium_channel_BK_alpha_subunit": [
        "IPR003929"
    ],
    
    # UCH Superfamily Components
    "Ubiquitin_carboxyl_terminal_hydrolase": [
        "IPR001578",
        "IPR036959",
    ],
    
    # Domain of Unknown Function (DUF) & Miscellaneous Families
    "DUF2254_membrane": [
        "IPR018723"
    ],
    "Lantibiotic_dehydratase": [
        "IPR006827"
    ],
    "Ribosomal_protein_L37_S30": [
        "IPR010793"
    ],
    "DUF4253": [
        "IPR025349"
    ]
}


