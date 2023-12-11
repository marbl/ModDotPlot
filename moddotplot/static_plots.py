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
    scale_y_continuous,
    geom_tile,
    coord_fixed,
    facet_grid,
    labs,
    element_line,
    element_text
)
import pandas as pd
import numpy as np
import math
import os
from moddotplot.const import (
    DIVERGING_PALETTES,
    QUALITATIVE_PALETTES,
    SEQUENTIAL_PALETTES,
)
import itertools
from typing import List
from moddotplot.estimate_identity import binomial_distance, containment
from palettable.colorbrewer import qualitative, sequential, diverging
from moddotplot.parse_fasta import printProgressBar

def check_st_en_equality(df):
        unequal_rows = df[(df['q_st'] != df['r_st']) | (df['q_en'] != df['r_en'])]
        unequal_rows.loc[:, ['q_en', 'r_en', 'q_st', 'r_st']] = unequal_rows[['r_en', 'q_en', 'r_st', 'q_st']].values
        
        df = pd.concat([df, unequal_rows], ignore_index=True)
        
        return df

def make_k(vals):
    return format(vals / 1e3, ",")


def make_scale(vals: list) -> str:
    scaled = [number / 1000000 for number in vals]
    return scaled


def paired_bed_file(
    window_partitions: dict,
    input_name: str,
    id_threshold: int,
    output: str,
    no_bed: bool,
    no_plot: bool,
    no_hist: bool,
    palette: str,
    palette_orientation: str,
    width: int,
    height: int,
    dpi: int,
    k: int,
    is_freq: bool,
    xlim: int,
    custom_colors: List,
    custom_breakpoints: List,
) -> None:
    assert 50 < id_threshold < 100
    bed = []
    bed.append(
        (
            "#query_name",
            "query_start",
            "query_end",
            "reference_name",
            "reference_start",
            "reference_end",
            "perID_by_events",
        )
    )
    counter = 0
    n = (len(window_partitions) * len(window_partitions)) / 2 + (
        len(window_partitions) / 2
    )
    progress_thresholds = round(n / 77)
    printProgressBar(0, n, prefix="Progress:", suffix="Complete", length=40)

    for w in itertools.combinations_with_replacement(window_partitions, 2):
        counter += 1
        if counter % progress_thresholds == 0:
            printProgressBar(
                counter, n, prefix="Progress:", suffix="Completed", length=40
            )
        if input_name:
            query_name = input_name
        else:
            query_name = "input_sequence"
        query_start, query_end = w[0].split("-")
        reference_start, reference_end = w[1].split("-")

        perID = binomial_distance(
            containment(set(window_partitions[w[0]]), set(window_partitions[w[1]])), k
        )
        if perID * 100 >= id_threshold:
            bed.append(
                (
                    query_name,
                    int(query_start),
                    int(query_end),
                    query_name,
                    int(reference_start),
                    int(reference_end),
                    perID * 100,
                )
            )
    if not output:
        bedfile_output = input_name + ".bed"
    else:
        if not os.path.exists(output):
            os.makedirs(output)
        bedfile_output = os.path.join(output, input_name + ".bed")

    if not no_bed:
        with open(bedfile_output, "w") as bedfile:
            for row in bed:
                bedfile.write("\t".join(map(str, row)) + "\n")
        print("\n")
        print(f"Self identity matrix complete! Saved to {bedfile_output}\n")

    if not no_plot:
        print(f"Creating plots...\n")
        create_plots(
            [bed],
            output if output else ".",
            input_name,
            input_name,
            palette,
            palette_orientation,
            no_hist,
            width,
            height,
            dpi,
            is_freq,
            xlim,
            custom_colors,
            custom_breakpoints,
        )


def paired_bed_file_a_vs_b(
    window_partitions_a: dict,
    window_partitions_b: dict,
    a_name: str,
    b_name: str,
    resolution: int,
    id_threshold: int,
    no_bed: bool,
    no_plot: bool,
    palette: str,
    palette_orientation: str,
    height: int,
    width: int,
    dpi: int,
    k: int,
    is_freq: bool,
    output: str,
    custom_colors: List,
    custom_breakpoints: List,
    compare_only: bool,
) -> None:
    bed = []
    bed.append(
        (
            "#query_name",
            "query_start",
            "query_end",
            "reference_name",
            "reference_start",
            "reference_end",
            "perID_by_events",
        )
    )
    assert id_threshold > 50 and id_threshold < 100
    counter = 0
    n = len(window_partitions_a) * len(window_partitions_b)
    progress_thresholds = round(n / 77)
    printProgressBar(0, n, prefix="Progress:", suffix="Complete", length=40)
    for i in window_partitions_a:
        for j in window_partitions_b:
            counter += 1
            if counter % progress_thresholds == 0:
                printProgressBar(
                    counter, n, prefix="Progress:", suffix="Completed", length=40
                )
            reference_start, reference_end = i.split("-")
            query_start, query_end = j.split("-")
            perID = binomial_distance(
                containment(set(window_partitions_a[i]), set(window_partitions_b[j])), k
            )
            if perID * 100 >= id_threshold:
                if perID * 100 >= id_threshold:
                    bed.append(
                        (
                            b_name,
                            int(query_start),
                            int(query_end),
                            a_name,
                            int(reference_start),
                            int(reference_end),
                            perID * 100,
                        )
                    )
    bedfile_output = ""
    image_prefix = a_name + "_" + b_name
    if not output:
        bedfile_output = a_name + "_" + b_name + ".bed"
    else:
        if not os.path.exists(output):
            os.makedirs(output)
        bedfile_output = os.path.join(output, a_name + "_" + b_name + ".bed")
        image_prefix = os.path.join(output, image_prefix)

    image_output = a_name + "_" + b_name
    if not no_bed:
        with open(bedfile_output, "w") as bedfile:
            for row in bed:
                bedfile.write("\t".join(map(str, row)) + "\n")
        print("\n")
        print(f"Pairwise matrix complete! Saved to {bedfile_output}\n")

    sdf = read_df(
        [bed], palette, palette_orientation, is_freq, custom_colors, custom_breakpoints
    )
    cdf = sdf.dropna(subset=["discrete"])

    if not no_plot:
        print(f"Creating plots... \n")
        c = make_dot(cdf, image_output, palette, palette_orientation, custom_colors)
        ggsave(
            c,
            width=width,
            height=height,
            dpi=dpi,
            format="png",
            filename=image_prefix + ".png",
            verbose=False,
        )
        ggsave(
            c,
            width=width,
            height=height,
            dpi=dpi,
            format="pdf",
            filename=image_prefix + ".pdf",
            verbose=False,
        )
        if compare_only:
            histy = make_hist(
                sdf, palette, palette_orientation, custom_colors, custom_breakpoints
            )
            ggsave(
                histy,
                width=3,
                height=3,
                dpi=dpi,
                format="pdf",
                filename=image_prefix + "_HIST.pdf",
                verbose=False,
            )
            ggsave(
                histy,
                width=3,
                height=3,
                dpi=dpi,
                format="png",
                filename=image_prefix + "_HIST.png",
                verbose=False,
            )
            print(
                f"{image_prefix}.png, {image_prefix}.pdf, {image_prefix}_HIST.pdf, and {image_prefix}_HIST.png saved sucessfully. \n"
            )
        else:
            print(f"{image_prefix}.png and {image_prefix}.pdf saved sucessfully. \n")


def get_colors(sdf, ncolors, is_freq, custom_breakpoints):
    assert ncolors > 2 and ncolors < 12
    bot = math.floor(min(sdf["perID_by_events"]))
    top = 100
    interval = (top - bot) / ncolors
    breaks = []
    if is_freq:
        breaks = np.unique(
            np.quantile(sdf["perID_by_events"], np.arange(0, 1.01, 1 / ncolors))
        )
    else:
        breaks = [bot + i * interval for i in range(ncolors + 1)]

    if custom_breakpoints:
        breaks = np.asfarray(custom_breakpoints)
    labels = np.arange(len(breaks) - 1)

    # corner case of only one %id value
    if len(breaks) == 1:
        return pd.factorize([1] * len(sdf["perID_by_events"]))[0]
    else:
        tmp = pd.cut(
            sdf["perID_by_events"], bins=breaks, labels=labels, include_lowest=True
        )
        return tmp


def read_df(
    pj, palette, palette_orientation, is_freq, custom_colors, custom_breakpoints
):
    data = pj[0]
    df = pd.DataFrame(data[1:], columns=data[0])
    hexcodes = []
    new_hexcodes = []
    if palette in DIVERGING_PALETTES:
        function_name = getattr(diverging, palette)
        hexcodes = function_name.hex_colors
        if palette_orientation == "+":
            palette_orientation = "-"
        else:
            palette_orientation = "+"
    elif palette in QUALITATIVE_PALETTES:
        function_name = getattr(qualitative, palette)
        hexcodes = function_name.hex_colors
    elif palette in SEQUENTIAL_PALETTES:
        function_name = getattr(sequential, palette)
        hexcodes = function_name.hex_colors
    else:
        function_name = getattr(sequential, "Spectral_11")
        palette_orientation = "-"
        hexcodes = function_name.hex_colors

    if palette_orientation == "-":
        new_hexcodes = hexcodes[::-1]
    else:
        new_hexcodes = hexcodes

    if custom_colors:
        new_hexcodes = custom_colors

    ncolors = len(new_hexcodes)

    # Get colors for each row based on the values in the dataframe
    df["discrete"] = get_colors(df, ncolors, is_freq, custom_breakpoints)

    # Rename columns if they have different names in the dataframe
    if "query_name" in df.columns or "#query_name" in df.columns:
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
    side_length = row["window"] * np.sqrt(2) / 2
    base = np.array([[1, 0], [0, 1], [-1, 0], [0, -1]]) * side_length
    trans = base + np.array([row["w"], row["z"]])
    df = pd.DataFrame(trans, columns=["w", "z"])
    df["discrete"] = int(row["discrete"])
    df["group"] = int(row["group"])
    return df


from plotnine import ggplot, geom_tile, scale_color_discrete, scale_fill_manual, theme, element_blank, element_line, element_text, scale_x_continuous, scale_y_continuous, coord_fixed, facet_grid, labs
import numpy as np

def make_dot(sdf, title_name, palette, palette_orientation, colors):
    hexcodes = []
    new_hexcodes = []
    if palette in DIVERGING_PALETTES:
        function_name = getattr(diverging, palette)
        hexcodes = function_name.hex_colors
        if palette_orientation == "+":
            palette_orientation = "-"
        else:
            palette_orientation = "+"
    elif palette in QUALITATIVE_PALETTES:
        function_name = getattr(qualitative, palette)
        hexcodes = function_name.hex_colors
    elif palette in SEQUENTIAL_PALETTES:
        function_name = getattr(sequential, palette)
        hexcodes = function_name.hex_colors
    else:
        function_name = getattr(sequential, "Spectral_11")
        palette_orientation = "-"
        hexcodes = function_name.hex_colors

    if palette_orientation == "-":
        new_hexcodes = hexcodes[::-1]
    else:
        new_hexcodes = hexcodes
    if colors:
        new_hexcodes = colors
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max())
    window = max(sdf["q_en"] - sdf["q_st"])
    p = (
        ggplot(sdf)
        + geom_tile(
            aes(x="q_st", y="r_st", fill="discrete", height=window, width=window)
        )
        + scale_color_discrete(guide=False)
        + scale_fill_manual(
            values=new_hexcodes,
            guide=False,
        )
        + theme(
            legend_position="none",
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            plot_background=element_blank(),
            panel_background=element_blank(),
            axis_line=element_line(color="black", size=2),  # Adjust axis line size
            axis_text=element_text(family="Helvetica", size=12),  # Change axis text font and size
            axis_ticks_major=element_line(size=6),
        )
        + scale_x_continuous(labels=make_scale, limits=[0, max_val])
        + scale_y_continuous(labels=make_scale, limits=[0, max_val])
        + coord_fixed(ratio=1)
        + facet_grid("r ~ q")
        + labs(x="Genomic position (Mbp)", y="", title=title_name, size=40)
    )

    # Adjust x-axis label size
    p += theme(axis_title_x=element_text(size=18))

    return p



def make_tri(
    df_d,
    title_name,
    palette,
    palette_orientation,
    is_freq,
    custom_colors,
    custom_breakpoints,
    xlim
):
    hexcodes = []
    new_hexcodes = []
    if palette in DIVERGING_PALETTES:
        function_name = getattr(diverging, palette)
        hexcodes = function_name.hex_colors
        if palette_orientation == "+":
            palette_orientation = "-"
        else:
            palette_orientation = "+"
    elif palette in QUALITATIVE_PALETTES:
        function_name = getattr(qualitative, palette)
        hexcodes = function_name.hex_colors
    elif palette in SEQUENTIAL_PALETTES:
        function_name = getattr(sequential, palette)
        hexcodes = function_name.hex_colors
    else:
        function_name = getattr(sequential, "Spectral_11")
        palette_orientation = "-"
        hexcodes = function_name.hex_colors

    if palette_orientation == "-":
        new_hexcodes = hexcodes[::-1]
    else:
        new_hexcodes = hexcodes

    if custom_colors:
        new_hexcodes = custom_colors
    p_tri = (
        ggplot(df_d)
        + aes(x="w_new", y="z_new", group="group", fill="discrete")
        + geom_polygon()
        + scale_color_discrete(guide=False)
        + scale_fill_gradientn(  # TODO: Replace this with built in color palettes.
            colors=new_hexcodes,
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
            panel_background=element_blank(),
        )
        + xlab("Genomic Position (Mbp)")
        + ggtitle(title_name)
    )

    if xlim:
        print(f"Resizing x-axis to {xlim}")
        p_tri += coord_cartesian(xlim=(0, xlim))

    return p_tri


def make_hist(sdf, palette, palette_orientation, custom_colors, custom_breakpoints):
    hexcodes = []
    new_hexcodes = []
    if palette in DIVERGING_PALETTES:
        function_name = getattr(diverging, palette)
        hexcodes = function_name.hex_colors
        if palette_orientation == "+":
            palette_orientation = "-"
        else:
            palette_orientation = "+"
    elif palette in QUALITATIVE_PALETTES:
        function_name = getattr(qualitative, palette)
        hexcodes = function_name.hex_colors
    elif palette in SEQUENTIAL_PALETTES:
        function_name = getattr(sequential, palette)
        hexcodes = function_name.hex_colors
    else:
        function_name = getattr(sequential, "Spectral_11")
        palette_orientation = "-"
        hexcodes = function_name.hex_colors

    if palette_orientation == "-":
        new_hexcodes = hexcodes[::-1]
    else:
        new_hexcodes = hexcodes

    if custom_colors:
        new_hexcodes = custom_colors
    bot = np.quantile(sdf["perID_by_events"], q=0.001)
    count = sdf.shape[0]
    extra = ""

    if count > 1e5:
        extra = "\n(thousands)"

    p = (
        ggplot(data=sdf, mapping=aes(x="perID_by_events", fill="discrete"))
        + geom_histogram(bins=300)
        + scale_color_cmap(cmap_name="plasma")
        + scale_fill_manual(new_hexcodes)
        + theme(legend_position="none")
        + coord_cartesian(xlim=(bot, 100))
        + xlab("% identity estimate")
        + ylab("# of estimates{}".format(extra))
    )
    return p


def create_plots(
    sdf,
    directory,
    output,
    input_sequence,
    palette,
    palette_orientation,
    no_hist,
    width,
    height,
    dpi,
    is_freq,
    xlim,
    custom_colors,
    custom_breakpoints,
):
    df = read_df(
        sdf, palette, palette_orientation, is_freq, custom_colors, custom_breakpoints
    )
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
    tri = make_tri(
        df_d,
        input_sequence,
        palette,
        palette_orientation,
        is_freq,
        custom_colors,
        custom_breakpoints,
        xlim
    )

    plot_filename = f"{directory}/{output}"
    symmetric = check_st_en_equality(sdf)
    c = make_dot(symmetric, plot_filename, palette, palette_orientation, custom_colors)
    print("Plots created! \n")

    print(f"Saving plots to {plot_filename}... \n")
    ggsave(
        c,
        width=width,
        height=width,
        dpi=dpi,
        format="pdf",
        filename=plot_filename + "_FULL.pdf",
        verbose=False,
    )
    ggsave(
        c,
        width=width,
        height=width,
        dpi=dpi,
        format="png",
        filename=plot_filename + "_FULL.png",
        verbose=False,
    )
    ggsave(
        tri,
        width=width,
        height=height,
        dpi=dpi,
        format="pdf",
        filename=plot_filename + "_TRI.pdf",
        verbose=False,
    )
    ggsave(
        tri,
        width=width,
        height=height,
        dpi=dpi,
        format="png",
        filename=plot_filename + "_TRI.png",
        verbose=False,
    )

    if not no_hist:
        histy = make_hist(
            sdf, palette, palette_orientation, custom_colors, custom_breakpoints
        )
        ggsave(
            histy,
            width=3,
            height=3,
            dpi=dpi,
            format="pdf",
            filename=plot_filename + "_HIST.pdf",
            verbose=False,
        )
        ggsave(
            histy,
            width=3,
            height=3,
            dpi=dpi,
            format="png",
            filename=plot_filename + "_HIST.png",
            verbose=False,
        )
        print(
            f"{plot_filename}_TRI.png, {plot_filename}_TRI.pdf, {plot_filename}_FULL.png, {plot_filename}_FULL.png, {plot_filename}_HIST.png and {plot_filename}_HIST.pdf, saved sucessfully. \n"
        )

    else:
        print(
            f"{plot_filename}_TRI.png {plot_filename}_TRI.pdf, {plot_filename}_FULL.png and {plot_filename}_FULL.png saved sucessfully. \n"
        )
