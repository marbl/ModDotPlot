VERSION = "0.9.4"
COLS = [
    "#query_name",
    "query_start",
    "query_end",
    "reference_name",
    "reference_start",
    "reference_end",
    "perID_by_events",
]

ASCII_ART = """
  __  __           _   _____        _     _____  _       _   
 |  \/  |         | | |  __ \      | |   |  __ \| |     | |  
 | \  / | ___   __| | | |  | | ___ | |_  | |__) | | ___ | |_ 
 | |\/| |/ _ \ / _` | | |  | |/ _ \| __| |  ___/| |/ _ \| __|
 | |  | | (_) | (_| | | |__| | (_) | |_  | |    | | (_) | |_ 
 |_|  |_|\___/ \__,_| |_____/ \___/ \__| |_|    |_|\___/ \__|
"""

DIVERGING_PALETTES = [
    "BrBG_3",
    "BrBG_4",
    "BrBG_5",
    "BrBG_6",
    "BrBG_7",
    "BrBG_8",
    "BrBG_9",
    "BrBG_10",
    "BrBG_11",
    "PRGn_3",
    "PRGn_4",
    "PRGn_5",
    "PRGn_6",
    "PRGn_7",
    "PRGn_8",
    "PRGn_9",
    "PRGn_10",
    "PRGn_11",
    "PiYG_3",
    "PiYG_4",
    "PiYG_5",
    "PiYG_6",
    "PiYG_7",
    "PiYG_8",
    "PiYG_9",
    "PiYG_10",
    "PiYG_11",
    "PuOr_3",
    "PuOr_4",
    "PuOr_5",
    "PuOr_6",
    "PuOr_7",
    "PuOr_8",
    "PuOr_9",
    "PuOr_10",
    "PuOr_11",
    "RdBu_3",
    "RdBu_4",
    "RdBu_5",
    "RdBu_6",
    "RdBu_7",
    "RdBu_8",
    "RdBu_9",
    "RdBu_10",
    "RdBu_11",
    "RdGy_3",
    "RdGy_4",
    "RdGy_5",
    "RdGy_6",
    "RdGy_7",
    "RdGy_8",
    "RdGy_9",
    "RdGy_10",
    "RdGy_11",
    "RdYlBu_3",
    "RdYlBu_4",
    "RdYlBu_5",
    "RdYlBu_6",
    "RdYlBu_7",
    "RdYlBu_8",
    "RdYlBu_9",
    "RdYlBu_10",
    "RdYlBu_11",
    "RdYlGn_3",
    "RdYlGn_4",
    "RdYlGn_5",
    "RdYlGn_6",
    "RdYlGn_7",
    "RdYlGn_8",
    "RdYlGn_9",
    "RdYlGn_10",
    "RdYlGn_11",
    "Spectral_3",
    "Spectral_4",
    "Spectral_5",
    "Spectral_6",
    "Spectral_7",
    "Spectral_8",
    "Spectral_9",
    "Spectral_10",
    "Spectral_11",
]

QUALITATIVE_PALETTES = [
    "Accent_3",
    "Accent_4",
    "Accent_5",
    "Accent_6",
    "Accent_7",
    "Accent_8",
    "Dark2_3",
    "Dark2_4",
    "Dark2_5",
    "Dark2_6",
    "Dark2_7",
    "Dark2_8",
    "Paired_3",
    "Paired_4",
    "Paired_5",
    "Paired_6",
    "Paired_7",
    "Paired_8",
    "Paired_9",
    "Paired_10",
    "Paired_11",
    "Paired_12",
    "Pastel1_3",
    "Pastel1_4",
    "Pastel1_5",
    "Pastel1_6",
    "Pastel1_7",
    "Pastel1_8",
    "Pastel1_9",
    "Pastel2_3",
    "Pastel2_4",
    "Pastel2_5",
    "Pastel2_6",
    "Pastel2_7",
    "Pastel2_8",
    "Set1_3",
    "Set1_4",
    "Set1_5",
    "Set1_6",
    "Set1_7",
    "Set1_8",
    "Set1_9",
    "Set2_3",
    "Set2_4",
    "Set2_5",
    "Set2_6",
    "Set2_7",
    "Set2_8",
    "Set3_3",
    "Set3_4",
    "Set3_5",
    "Set3_6",
    "Set3_7",
    "Set3_8",
    "Set3_9",
    "Set3_10",
    "Set3_11",
    "Set3_12",
]

SEQUENTIAL_PALETTES = [
    "Blues_3",
    "Blues_4",
    "Blues_5",
    "Blues_6",
    "Blues_7",
    "Blues_8",
    "Blues_9",
    "BuGn_3",
    "BuGn_4",
    "BuGn_5",
    "BuGn_6",
    "BuGn_7",
    "BuGn_8",
    "BuGn_9",
    "BuPu_3",
    "BuPu_4",
    "BuPu_5",
    "BuPu_6",
    "BuPu_7",
    "BuPu_8",
    "BuPu_9",
    "GnBu_3",
    "GnBu_4",
    "GnBu_5",
    "GnBu_6",
    "GnBu_7",
    "GnBu_8",
    "GnBu_9",
    "Greens_3",
    "Greens_4",
    "Greens_5",
    "Greens_6",
    "Greens_7",
    "Greens_8",
    "Greens_9",
    "Greys_3",
    "Greys_4",
    "Greys_5",
    "Greys_6",
    "Greys_7",
    "Greys_8",
    "Greys_9",
    "OrRd_3",
    "OrRd_4",
    "OrRd_5",
    "OrRd_6",
    "OrRd_7",
    "OrRd_8",
    "OrRd_9",
    "Oranges_3",
    "Oranges_4",
    "Oranges_5",
    "Oranges_6",
    "Oranges_7",
    "Oranges_8",
    "Oranges_9",
    "PuBu_3",
    "PuBu_4",
    "PuBu_5",
    "PuBu_6",
    "PuBu_7",
    "PuBu_8",
    "PuBu_9",
    "PuBuGn_3",
    "PuBuGn_4",
    "PuBuGn_5",
    "PuBuGn_6",
    "PuBuGn_7",
    "PuBuGn_8",
    "PuBuGn_9",
    "PuRd_3",
    "PuRd_4",
    "PuRd_5",
    "PuRd_6",
    "PuRd_7",
    "PuRd_8",
    "PuRd_9",
    "Purples_3",
    "Purples_4",
    "Purples_5",
    "Purples_6",
    "Purples_7",
    "Purples_8",
    "Purples_9",
    "RdPu_3",
    "RdPu_4",
    "RdPu_5",
    "RdPu_6",
    "RdPu_7",
    "RdPu_8",
    "RdPu_9",
    "Reds_3",
    "Reds_4",
    "Reds_5",
    "Reds_6",
    "Reds_7",
    "Reds_8",
    "Reds_9",
    "YlGn_3",
    "YlGn_4",
    "YlGn_5",
    "YlGn_6",
    "YlGn_7",
    "YlGn_8",
    "YlGn_9",
    "YlGnBu_3",
    "YlGnBu_4",
    "YlGnBu_5",
    "YlGnBu_6",
    "YlGnBu_7",
    "YlGnBu_8",
    "YlGnBu_9",
    "YlOrBr_3",
    "YlOrBr_4",
    "YlOrBr_5",
    "YlOrBr_6",
    "YlOrBr_7",
    "YlOrBr_8",
    "YlOrBr_9",
    "YlOrRd_3",
    "YlOrRd_4",
    "YlOrRd_5",
    "YlOrRd_6",
    "YlOrRd_7",
    "YlOrRd_8",
    "YlOrRd_9",
]
