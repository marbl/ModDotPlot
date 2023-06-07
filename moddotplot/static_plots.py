from plotnine import (
    ggsave,
    ggplot,
    aes,
    geom_histogram,
    geom_polygon,
    scale_color_discrete,
    scale_fill_gradientn,
    element_blank,
    theme,
    xlab,
    ggtitle,
    scale_fill_manual,
    scale_color_cmap,
    coord_cartesian,
    ylab,
    scale_x_continuous,
    xlim,
)
import pandas as pd
import numpy as np
import math
from moddotplot.const import COLS
import itertools
from typing import List
from moddotplot.estimate_identity import binomial_distance, containment


def make_k(vals):
    return format(vals / 1e3, ",")


def make_scale(vals: float) -> str:
    return format(vals / 1e6, ",")


def paired_bed_file(
    window_partitions: dict,
    input_name: str,
    id_threshold: int,
    density: int,
    output: str,
    seq_length: int,
    no_plot: bool,
    k: int
) -> None:
    """
    Given a dictionary of window partitions, creates a BED file containing pairs of
    genomic regions with specified identity threshold.
    Args:
    window_partitions: a dictionary of window partitions containing genomic regions.
    input_name: a string representing the name of the input sequence.
    id_threshold: an integer specifying the identity threshold between 50 and 100.
    density: an integer specifying the density of the genomic regions.
    output: a string representing the name of the output BED file.
    Returns:
    None
    """
    assert id_threshold > 50 and id_threshold < 100

    cols: List[str] = COLS
    bed: List[List[str]] = []
    for w in itertools.combinations_with_replacement(window_partitions, 2):
        if input_name:
            query_name = input_name
        else:
            query_name = "input_sequence"
        query_start, query_end = w[0].split("-")
        reference_start, reference_end = w[1].split("-")
        # TODO: Update with Jaccard when available
        perID = binomial_distance(
            containment(set(window_partitions[w[0]]), set(window_partitions[w[1]])), k
        )
        if perID * 100 >= id_threshold:
            # TODO: Add strand orientation into bed file
            bed.append(
                [
                    query_name,
                    int(query_start),
                    int(query_end),
                    query_name,
                    int(reference_start),
                    int(reference_end),
                    perID * 100,
                ]
            )
    df = pd.DataFrame(bed, columns=cols)
    if not output:
        bedfile_output = input_name + ".bed"
        df.to_csv(bedfile_output, sep="\t")
        print(f"Self identity matrix complete! Saved to {bedfile_output} \n")
        if not no_plot:
            print(f"Creating plots... \n")
            create_plots([df], input_name, input_name)
    else:
        bedfile_output = output + ".bed"
        df.to_csv(bedfile_output, sep="\t")
        print(f"Identity computed! Saved to {bedfile_output} \n")
        if not no_plot:
            print(f"Creating plots... \n")
            create_plots([df], output, input_name)

def paired_bed_file_a_vs_b(
    window_partitions_a: dict,
    window_partitions_b: dict,
    a_name: str,
    b_name: str,
    resolution: int,
    id_threshold: int,
    k: int,
) -> None:

    cols: List[str] = COLS
    bed: List[List[str]] = []
    assert id_threshold > 50 and id_threshold < 100
    for i in window_partitions_a:
        for j in window_partitions_b:
            reference_start, reference_end = i.split("-")
            query_start, query_end = j.split("-")
            perID = binomial_distance(
                containment(set(window_partitions_a[i]), set(window_partitions_b[j])), k
            )
            if perID * 100 >= id_threshold:
            # TODO: Add strand orientation into bed file
                bed.append(
                    [
                        a_name,
                        int(reference_start),
                        int(reference_end),
                        b_name,
                        int(query_start),
                        int(query_end),
                        perID * 100,
                    ]
                )
    df = pd.DataFrame(bed, columns=cols)
    bedfile_output = a_name + "_" + b_name + ".bed"
    df.to_csv(bedfile_output, sep="\t")

def get_colors(sdf, ncolors, is_freq):
    assert ncolors > 2 and ncolors < 12
    bot = math.floor(min(sdf["perID_by_events"]))
    top = 100
    # TODO: Sort by frequency if arg is selected.
    breaks = np.unique(
        np.quantile(sdf["perID_by_events"], np.arange(0, 1.01, 1 / ncolors))
    )
    labels = np.arange(len(breaks) - 1)
    # corner case of only one %id value
    if len(breaks) == 1:
        return pd.factorize([1] * len(sdf["perID_by_events"]))[0]
    return pd.cut(
        sdf["perID_by_events"], bins=breaks, labels=labels, include_lowest=True
    )


def read_bedpe(all_files):
    # Load each file into a list of dataframes using fread function from pandas
    l = [pd.read_csv(file, sep="\t") for file in all_files]
    # Concatenate the dataframes into a single dataframe using rbindlist function from pandas
    df = pd.concat(l)

    # Get colors for each row based on the values in the dataframe
    df["discrete"] = get_colors(df, 11, False)

    # Rename columns if they have different names in the dataframe
    if "#query_name" in df.columns:
        df.rename(
            columns={
                "#query_name": "q",
                "query_start": "q_st",
                "query_end": "q_en",
                "reference_name": "r",
                "reference_start": "r_st",
                "reference_end": "r_en",
            },
            inplace=True,
        )

    # Calculate the window size
    window = max(df["q_en"] - df["q_st"])

    # Calculate the position of the first and second intervals
    df["first_pos"] = df["q_st"] / window
    df["second_pos"] = df["r_st"] / window

    return df


def read_df(pj):
    df = pd.concat(pj)
    # TODO: Make this cmd line arg
    ncolors = 8

    # Get colors for each row based on the values in the dataframe
    df["discrete"] = get_colors(df, ncolors, False)

    # Rename columns if they have different names in the dataframe
    if "#query_name" in df.columns:
        df.rename(
            columns={
                "#query_name": "q",
                "query_start": "q_st",
                "query_end": "q_en",
                "reference_name": "r",
                "reference_start": "r_st",
                "reference_end": "r_en",
            },
            inplace=True,
        )

    # Calculate the window size
    window = max(df["q_en"] - df["q_st"])

    # Calculate the position of the first and second intervals
    df["first_pos"] = df["q_st"] / window
    df["second_pos"] = df["r_st"] / window

    return df


def diamond(row):
    side_length = float(row["window"])
    x = float(row["w"])
    y = float(row["z"])

    base = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]]) * np.sqrt(2) / 2
    trans = (base * side_length) + np.array([x, y])
    df = pd.DataFrame(trans, columns=["w", "z"])
    df["discrete"] = int(row["discrete"])
    df["group"] = int(row["group"])
    return df

def make_dot(sdf, rname=""):
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max())
    window = max(sdf["query_end"] - sdf["query_start"])
    
    p = (
        ggplot(sdf) +
        geom_tile(aes(x="q_st", y="r_st", fill="discrete", height=window, width=window)) +
        theme_cowplot() +
        scale_fill_brewer(palette="Spectral", direction=-1) +
        theme(legend_position="none") +
        scale_x_continuous(labels=make_scale, limits=[0, max_val]) +
        scale_y_continuous(labels=make_scale, limits=[0, max_val]) +
        coord_fixed(ratio=1) +
        facet_grid("r ~ q") +
        labs(x="Genomic position (Mbp)", y="", title=rname)
    )
    
    return p

def make_tri(df_d, title_name):
    p_tri = (
        ggplot(df_d)
        + aes(x="w_new", y="z_new", group="group", fill="discrete")
        + geom_polygon()
        + scale_color_discrete(guide=False)
        + scale_fill_gradientn(  # TODO: Replace this with built in color palettes.
            colors=[
                "#5E4FA2",
                "#3388BD",
                "#66C2A4",
                "#ABDEA4",
                "#E6F698",
                "#FFFFBF",
                "#FEE08B",
                "#FDAE61",
                "#F46D43",
                "#D53E4F",
                "#9F0142",
            ],
            guide=False,
        )
        + theme(
            axis_title_y=element_blank(),
            axis_line_y=element_blank(),
            axis_text_y=element_blank(),
            axis_ticks_major_y=element_blank(),
            axis_ticks_minor_y=element_blank(),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            plot_background=element_blank(),
            panel_background=element_blank()
        )
        + xlab("Genomic Position (Mbp)")
        + ggtitle(title_name)
    )
    return p_tri


def make_hist(sdf):

    bot = np.quantile(sdf["perID_by_events"], q=0.001)
    count = sdf.shape[0]
    extra = ""

    if count > 1e5:
        extra = "\n(thousands)"

    p = (
        ggplot(data=sdf, mapping=aes(x="perID_by_events", fill="discrete"))
        + geom_histogram(bins=300)
        + scale_color_cmap(cmap_name="plasma")
        +
        # Make sure there's enough colors
        scale_fill_manual(
            [
                "#5E4FA2",
                "#3388BD",
                "#66C2A4",
                "#ABDEA4",
                "#E6F698",
                "#FFFFBF",
                "#FEE08B",
                "#FDAE61",
                "#F46D43",
                "#D53E4F",
                "#9F0142",
            ]
        )
        + theme(legend_position="none")
        + coord_cartesian(xlim=(bot, 100))
        + xlab("% identity estimate")
        + ylab("# of estimates{}".format(extra))
    )
    return p


def create_plots(sdf, output, input_sequence):

    df = read_df(sdf)
    sdf = df

    sdf["w"] = sdf["first_pos"] + sdf["second_pos"]
    sdf["z"] = -sdf["first_pos"] + sdf["second_pos"]
    window = max(sdf["q_en"] - sdf["q_st"])
    tri_scale = max(sdf["q_st"]) / max(sdf["w"])
    sdf["window"] = max(sdf["q_en"] - sdf["q_st"]) / tri_scale
    sdf["group"] = np.arange(sdf.shape[0])
    sdf["w_new"] = sdf["w"] * tri_scale
    sdf["z_new"] = sdf["z"] * window
    df_d = pd.concat([diamond(row) for _, row in sdf.iterrows()], ignore_index=True)

    # TODO: Scale based on size of the input genome, cant always assume Mbp is appropriate
    df_d["w_new"] = df_d["w"] * tri_scale / 1000000
    df_d["z_new"] = df_d["z"] * window * 2 / 3 / 1000000
    tri = make_tri(df_d, input_sequence)

    plot_filename = f"{output}.png"
    print("Plots created! \n")

    print(f"Saving plots to {plot_filename}... \n")
    ggsave(tri, width=13, height=5, dpi=600, filename=plot_filename, verbose=False)
    print(f"{plot_filename} saved sucessfully.")

