from ast import excepthandler
from plotnine import (
    ggsave,
    ggplot,
    aes,
    geom_histogram,
    scale_color_discrete,
    element_blank,
    theme,
    xlab,
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
    element_text,
    theme_light,
    geom_blank,
    annotate,
    element_rect,
    coord_flip,
    theme_minimal,
    geom_raster,
)
import cairosvg
import pandas as pd
import numpy as np
from PIL import Image
import patchworklib as pw
import math
import os
import xml.etree.ElementTree as ET
import sys
import re
from moddotplot.parse_fasta import printProgressBar

from moddotplot.const import (
    DIVERGING_PALETTES,
    QUALITATIVE_PALETTES,
    SEQUENTIAL_PALETTES,
)
from typing import List
from palettable.colorbrewer import qualitative, sequential, diverging


def is_plot_empty(p):
    # Check if the plot has data or any layers
    return len(p.layers) == 0 and p.data.empty


def check_pascal(single_val, double_val):
    try:
        if len(single_val) == 2:
            assert len(double_val) == 1
        elif len(single_val) == 3:
            assert len(double_val) == 3
        elif len(single_val) == 4:
            assert len(double_val) == 6
        elif len(single_val) == 5:
            assert len(double_val) == 10
        elif len(single_val) == 6:
            assert len(double_val) == 15
        elif len(single_val) == 0:
            assert len(double_val) == (1 or 3 or 6 or 10 or 15)
    except AssertionError as e:
        print(
            f"Missing bed files required to create grid. Please verify all bed files are included."
        )
        sys.exit(8)


def reverse_pascal(double_vals):
    if len(double_vals) == 1:
        return 2
    elif len(double_vals) == 3:
        return 3
    elif len(double_vals) == 6:
        return 4
    elif len(double_vals) == 10:
        return 5
    elif len(double_vals) == 15:
        return 6
    else:
        sys.exit(9)


# Hardcoding for now, I have the formula.... I'm just lazy
def transpose_order(double_vals):
    if len(double_vals) == 1:
        return [0]
    elif len(double_vals) == 3:
        return [0, 1, 2]
    elif len(double_vals) == 6:
        return [0, 1, 3, 2, 4, 5]
    elif len(double_vals) == 10:
        return [0, 1, 4, 2, 5, 7, 3, 6, 8, 9]
    elif len(double_vals) == 15:
        return [0, 1, 5, 2, 6, 9, 3, 7, 10, 12, 4, 8, 11, 13, 14]


def check_st_en_equality(df):
    unequal_rows = df[(df["q_st"] != df["r_st"]) | (df["q_en"] != df["r_en"])]
    unequal_rows.loc[:, ["q_en", "r_en", "q_st", "r_st"]] = unequal_rows[
        ["r_en", "q_en", "r_st", "q_st"]
    ].values

    df = pd.concat([df, unequal_rows], ignore_index=True)

    return df


def make_k(vals):
    return [number / 1000 for number in vals]


def make_m(vals):
    return [number / 1e6 for number in vals]


def make_g(vals):
    return [number / 1e9 for number in vals]


def make_scale(vals: list) -> list:
    scaled = [number for number in vals]
    if scaled[-1] < 200000:
        return make_k(scaled)
    elif scaled[-1] > 200000000:
        return make_g(scaled)
    else:
        return make_m(scaled)


def overlap_axis(rotated_plot, filename, prefix):
    scale_factor = math.sqrt(2) + 0.04
    new_width = int(rotated_plot.width / scale_factor)
    new_height = int(rotated_plot.height / scale_factor)
    resized_rotated_plot = rotated_plot.resize((new_width, new_height), Image.LANCZOS)

    # Step 3: Overlay the resized rotated heatmap onto the original axes

    # Open the original heatmap with axes
    image_with_axes = Image.open(filename)

    # Create a blank image with the same size as the original
    final_image = Image.new("RGBA", image_with_axes.size)

    # Calculate the position to center the resized rotated image within the original plot area
    x_offset = (final_image.width - resized_rotated_plot.width) // 2
    y_offset = (final_image.height - resized_rotated_plot.height) // 2
    y_offset += 2400
    x_offset += 30

    # Paste the original image with axes onto the final image
    final_image.paste(image_with_axes, (0, 0))

    # Paste the resized rotated plot onto the final image
    final_image.paste(resized_rotated_plot, (x_offset, y_offset), resized_rotated_plot)
    width, height = final_image.size
    cropped_image = final_image.crop((0, height // 2.6, width, height))

    # Save or show the final image
    cropped_image.save(f"{prefix}_TRI.png")
    cropped_image.save(f"{prefix}_TRI.pdf", "PDF", resolution=100.0)

    # Remove temp files
    if os.path.exists(filename):
        os.remove(filename)


def get_colors(sdf, ncolors, is_freq, custom_breakpoints):
    assert ncolors > 2 and ncolors < 12
    try:
        bot = math.floor(min(sdf["perID_by_events"]))
    except ValueError:
        bot = 0
    top = 100.0
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


# TODO: Remove pandas dependency
def read_df_from_file(file_path):
    data = pd.read_csv(file_path, delimiter="\t")
    return data


def read_df(
    pj,
    palette,
    palette_orientation,
    is_freq,
    custom_colors,
    custom_breakpoints,
    from_file,
):
    df = ""
    if from_file is not None:
        df = from_file
    else:
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
        print(f"Palette {palette} not found. Defaulting to Spectral_11.\n")
        function_name = getattr(diverging, "Spectral_11")
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
    try:
        window = max(df["q_en"] - df["q_st"])
    except ValueError:
        window = 0

    # Calculate the position of the first and second intervals
    df["first_pos"] = df["q_st"] / window
    df["second_pos"] = df["r_st"] / window

    return df


def generate_breaks(number, min_breaks=6, max_breaks=8):
    # Determine the order of magnitude
    magnitude = 10 ** int(
        np.floor(np.log10(number))
    )  # Base power of 10 (e.g., 10M, 1M, 100K)

    # Find a reasonable step size
    possible_steps = [
        magnitude // d
        for d in [10, 5, 4, 2, 1]
        if (number // (magnitude // d)) in range(min_breaks, max_breaks + 1)
    ]

    step = (
        possible_steps[0] if possible_steps else magnitude // 5
    )  # Default to 1/5th if no exact match

    # Generate breakpoints
    breaks = list(range(0, ((number // step) + 1) * step, step))

    return breaks


def make_dot(
    sdf,
    title_name,
    palette,
    palette_orientation,
    colors,
    breaks,
    num_ticks,
    xlim,
    deraster,
    width,
):
    # Select the color palette
    if hasattr(diverging, palette):
        function_name = getattr(diverging, palette)
    elif hasattr(qualitative, palette):
        function_name = getattr(qualitative, palette)
    elif hasattr(sequential, palette):
        function_name = getattr(sequential, palette)
    else:
        function_name = diverging.Spectral_11  # Default palette
        palette_orientation = "-"

    hexcodes = function_name.hex_colors

    # Adjust palette orientation
    if palette in diverging.__dict__:
        palette_orientation = "-" if palette_orientation == "+" else "+"

    new_hexcodes = hexcodes[::-1] if palette_orientation == "-" else hexcodes
    if colors:
        new_hexcodes = colors  # Override colors if provided
    if not xlim:
        xlim = 0
    # Determine maximum genomic position for scaling
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max(), xlim)

    # If user provides breaks, convert to ints
    if not breaks:
        breaks = generate_breaks(int(max_val))
    else:
        [int(x) for x in breaks]
    xlim = xlim or 0
    # Compute window size (handling exceptions)
    try:
        window = max(sdf["q_en"] - sdf["q_st"])
    except ValueError:  # Empty dataframe case
        return ggplot(aes(x=[], y=[])) + theme_minimal()

    # Determine axis label scale based on genomic position size
    if max_val < 200_000:
        x_label = "Genomic Position (Kbp)"
    elif max_val < 200_000_000:
        x_label = "Genomic Position (Mbp)"
    else:
        x_label = "Genomic Position (Gbp)"

    # Create the plot
    common_theme = theme(
        legend_position="none",
        panel_grid_major=element_blank(),
        panel_grid_minor=element_blank(),
        plot_background=element_blank(),
        panel_background=element_blank(),
        axis_line=element_line(color="black"),
        axis_text=element_text(family=["DejaVu Sans"], size=width),
        axis_ticks_major=element_line(
            size=(width), color="black"
        ),  # Increased tick length
        title=element_blank(),  # Center title
        axis_title_x=element_text(size=(width * 1.2), family=["DejaVu Sans"]),
        strip_background=element_blank(),  # Remove facet strip background
        strip_text=element_text(
            size=(width * 1.2), family=["DejaVu Sans"]
        ),  # Customize facet label text size (optional)
    )

    # Construct the plot arguments
    ggplot_args = (
        ggplot(sdf)
        + scale_color_discrete(guide=False)
        + scale_fill_manual(values=new_hexcodes, guide=False)
        + common_theme
        + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + coord_fixed(ratio=1)
        + facet_grid("r ~ q")
        + labs(x=x_label, y="", title=title_name)
    )

    # Select either geom_raster or geom_tile depending on deraster flag
    p = ggplot_args + (geom_tile if deraster else geom_raster)(
        aes(x="q_st", y="r_st", fill="discrete", height=window, width=window)
    )

    return p


def make_dot_grid(
    sdf,
    title_name,
    palette,
    palette_orientation,
    colors,
    breaks,
    on_diagonal,
    xlim,
    deraster,
    width,
):
    # Select the color palette
    if hasattr(diverging, palette):
        function_name = getattr(diverging, palette)
    elif hasattr(qualitative, palette):
        function_name = getattr(qualitative, palette)
    elif hasattr(sequential, palette):
        function_name = getattr(sequential, palette)
    else:
        function_name = diverging.Spectral_11  # Default palette
        palette_orientation = "-"

    hexcodes = function_name.hex_colors

    # Adjust palette orientation
    if palette in diverging.__dict__:
        palette_orientation = "-" if palette_orientation == "+" else "+"

    new_hexcodes = hexcodes[::-1] if palette_orientation == "-" else hexcodes
    if colors:
        new_hexcodes = colors  # Override colors if provided
    if not xlim:
        xlim = 0
    # Determine maximum genomic position for scaling
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max(), xlim)

    # If user provides breaks, convert to ints
    if not breaks:
        breaks = generate_breaks(int(max_val))
    else:
        [int(x) for x in breaks]
    xlim = xlim or 0
    # Compute window size (handling exceptions)
    try:
        window = max(sdf["q_en"] - sdf["q_st"])
    except ValueError:  # Empty dataframe case
        return ggplot(aes(x=[], y=[])) + theme_minimal()

    # Determine axis label scale based on genomic position size
    if max_val < 200_000:
        x_label = "Genomic Position (Kbp)"
    elif max_val < 200_000_000:
        x_label = "Genomic Position (Mbp)"
    else:
        x_label = "Genomic Position (Gbp)"

    # Create the plot
    common_theme = theme(
        legend_position="none",
        panel_grid_major=element_blank(),
        panel_grid_minor=element_blank(),
        plot_background=element_blank(),
        panel_background=element_blank(),
        axis_line=element_line(color="black"),
        axis_text=element_text(family=["DejaVu Sans"], size=width),
        axis_ticks_major=element_line(
            size=(width), color="black"
        ),  # Increased tick length
        title=element_text(size=(width * 1.2), alpha=0),
        axis_title_x=element_text(size=(width * 1.2), family=["DejaVu Sans"]),
        strip_background=element_blank(),  # Remove facet strip background
        strip_text=element_text(
            size=(width * 1.2), family=["DejaVu Sans"]
        ),  # Customize facet label text size (optional)
    )

    # Construct the plot arguments
    ggplot_args = (
        ggplot(sdf)
        + scale_color_discrete(guide=False)
        + scale_fill_manual(values=new_hexcodes, guide=False)
        + common_theme
        + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + coord_fixed(ratio=1)
        + labs(x="", y="", title="")
    )

    # Select either geom_raster or geom_tile depending on deraster flag
    p = ggplot_args + (geom_tile if deraster else geom_raster)(
        aes(x="q_st", y="r_st", fill="discrete", height=window, width=window)
    )

    return p


def make_dot_final(
    sdf,
    width,
    palette,
    palette_orientation,
    colors,
    breaks,
    xlim,
    transpose=False,
    deraster=False,
):
    if hasattr(diverging, palette):
        function_name = getattr(diverging, palette)
    elif hasattr(qualitative, palette):
        function_name = getattr(qualitative, palette)
    elif hasattr(sequential, palette):
        function_name = getattr(sequential, palette)
    else:
        function_name = diverging.Spectral_11  # Default palette
        palette_orientation = "-"

    hexcodes = function_name.hex_colors

    # Adjust palette orientation
    if palette in diverging.__dict__:
        palette_orientation = "-" if palette_orientation == "+" else "+"

    new_hexcodes = hexcodes[::-1] if palette_orientation == "-" else hexcodes
    if colors:
        new_hexcodes = colors  # Override colors if provided
    if not xlim:
        xlim = 0
    # Determine maximum genomic position for scaling
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max(), xlim)

    # If user provides breaks, convert to ints
    if not breaks:
        breaks = generate_breaks(int(max_val))
    else:
        [int(x) for x in breaks]
    xlim = xlim or 0

    max_val = max(sdf["q_en"].max(), sdf["r_en"].max(), xlim)
    try:
        window = max(sdf["q_en"] - sdf["q_st"])
    except:
        p = (
            ggplot(aes(x=[], y=[]))
            + theme_minimal()
            + theme(
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
            )
        )
        return p

    x_col, y_col = ("r_st", "q_st") if transpose else ("q_st", "r_st")

    if deraster:
        p = (
            ggplot(sdf)
            + geom_tile(
                aes(x=x_col, y=y_col, fill="discrete", height=window, width=window)
            )
            + scale_color_discrete(guide=False)
            + scale_fill_manual(values=new_hexcodes, guide=False)
            + theme(
                legend_position="none",
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                plot_background=element_blank(),
                panel_background=element_blank(),
                axis_line=element_line(color="black"),
                axis_text=element_text(family=["DejaVu Sans"], size=width),
                axis_ticks_major=element_line(),
                title=element_text(family=["Dejavu Sans"]),
            )
            + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + coord_fixed(ratio=1)
            + labs(x=None, y=None, title=None)
        )
    else:
        p = (
            ggplot(sdf)
            + geom_raster(
                aes(x=x_col, y=y_col, fill="discrete", height=window, width=window)
            )
            + scale_color_discrete(guide=False)
            + scale_fill_manual(values=new_hexcodes, guide=False)
            + theme(
                legend_position="none",
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                plot_background=element_blank(),
                panel_background=element_blank(),
                axis_line=element_line(color="black"),
                axis_text=element_text(family=["DejaVu Sans"], size=width),
                axis_ticks_major=element_line(),
                title=element_text(family=["Dejavu Sans"]),
            )
            + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + coord_fixed(ratio=1)
            + labs(x=None, y=None, title=None)
        )

    p += theme(axis_title_x=element_blank(), axis_title_y=element_blank())

    return p


def make_tri(
    sdf,
    title_name,
    palette,
    palette_orientation,
    colors,
    breaks,
    xlim,
    num_ticks,
    deraster,
    width,
):
    # Select the color palette
    if hasattr(diverging, palette):
        function_name = getattr(diverging, palette)
    elif hasattr(qualitative, palette):
        function_name = getattr(qualitative, palette)
    elif hasattr(sequential, palette):
        function_name = getattr(sequential, palette)
    else:
        function_name = diverging.Spectral_11  # Default palette
        palette_orientation = "-"

    hexcodes = function_name.hex_colors

    # Adjust palette orientation
    if palette in diverging.__dict__:
        palette_orientation = "-" if palette_orientation == "+" else "+"

    new_hexcodes = hexcodes[::-1] if palette_orientation == "-" else hexcodes
    if colors:
        new_hexcodes = colors  # Override colors if provided
    if not xlim:
        xlim = 0
    # Determine maximum genomic position for scaling
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max(), xlim)

    # If user provides breaks, convert to ints
    if not breaks:
        breaks = generate_breaks(int(max_val))
    else:
        [int(x) for x in breaks]
    xlim = xlim or 0
    # Compute window size (handling exceptions)
    try:
        window = max(sdf["q_en"] - sdf["q_st"])
    except ValueError:  # Empty dataframe case
        return ggplot(aes(x=[], y=[])) + theme_minimal()

    # Determine axis label scale based on genomic position size
    if max_val < 200_000:
        x_label = "Genomic Position (Kbp)"
    elif max_val < 200_000_000:
        x_label = "Genomic Position (Mbp)"
    else:
        x_label = "Genomic Position (Gbp)"

    if not deraster:
        tri = (
            ggplot(sdf)
            + geom_raster(
                aes(x="q_st", y="r_st", fill="discrete", height=window, width=window),
                alpha=1.0,
            )  # Ensure full opacity
            + scale_fill_manual(values=new_hexcodes, guide=False)
            + scale_color_discrete(guide=False)
            + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + coord_fixed(ratio=1)
            + labs(x=x_label, y="", title=title_name)
            + theme(
                legend_position="none",
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                plot_background=element_blank(),
                panel_background=element_blank(),
                axis_text=element_text(family=["DejaVu Sans"], size=width),
                axis_line_x=element_line(),
                axis_line_y=element_blank(),
                axis_ticks_major_x=element_line(),
                axis_ticks_major_y=element_blank(),
                axis_ticks_major=element_line(size=(width)),
                title=element_text(size=(width * 1.4), hjust=0.5),
                axis_title_x=element_text(size=(width * 1.2), family=["DejaVu Sans"]),
                axis_text_y=element_blank(),
            )
        )
        axis = (
            ggplot(sdf)
            + geom_tile(
                aes(x="q_st", y="r_st", fill="discrete", height=window, width=window),
                alpha=0,
            )
            + scale_color_discrete(guide=False)
            + scale_fill_manual(values=new_hexcodes, guide=False)
            + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + coord_fixed(ratio=1)
            + labs(x="", y="", title=title_name)
            + theme(
                legend_position="none",
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                plot_background=element_blank(),
                panel_background=element_blank(),
                axis_line=element_line(color="black"),
                axis_text=element_text(family=["DejaVu Sans"], size=width),
                axis_ticks_major=element_line(),
                axis_line_x=element_line(),
                axis_line_y=element_blank(),
                axis_ticks_major_x=element_line(),
                axis_ticks_major_y=element_blank(),
                axis_text_x=element_line(),
                axis_text_y=element_blank(),
                plot_title=element_blank(),
                axis_title_x=element_text(size=(width * 1.2), family=["DejaVu Sans"]),
            )
        )
    else:
        tri = (
            ggplot(sdf)
            + geom_tile(
                aes(x="q_st", y="r_st", fill="discrete", height=window, width=window),
                alpha=1.0,
            )  # Ensure full opacity
            + scale_fill_manual(values=new_hexcodes, guide=False)
            + scale_color_discrete(guide=False)
            + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + coord_fixed(ratio=1)
            + labs(x=x_label, y="", title=title_name)
            + theme(
                legend_position="none",
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                plot_background=element_blank(),
                panel_background=element_blank(),
                axis_text=element_text(family=["DejaVu Sans"], size=width),
                axis_line_x=element_line(),
                axis_line_y=element_blank(),
                axis_ticks_major_x=element_line(),
                axis_ticks_major_y=element_blank(),
                axis_ticks_major=element_line(),
                axis_text_y=element_blank(),
                title=element_blank(),
                axis_title_x=element_text(size=(width * 1.2), family=["DejaVu Sans"]),
            )
        )
        axis = (
            ggplot(sdf)
            + geom_tile(
                aes(x="q_st", y="r_st", fill="discrete", height=window, width=window),
                alpha=0,
            )
            + scale_color_discrete(guide=False)
            + scale_fill_manual(values=new_hexcodes, guide=False)
            + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + coord_fixed(ratio=1)
            + labs(x="", y="", title="")
            + theme(
                legend_position="none",
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                plot_background=element_blank(),
                panel_background=element_blank(),
                axis_line=element_line(color="black"),
                axis_text=element_text(family=["DejaVu Sans"]),
                axis_ticks_major=element_line(),
                axis_line_x=element_line(),
                axis_line_y=element_blank(),
                axis_ticks_major_x=element_line(),
                axis_ticks_major_y=element_blank(),
                axis_text_x=element_line(),
                axis_text_y=element_blank(),
                plot_title=element_blank(),
                axis_title_x=element_text(size=(width * 1.2), family=["DejaVu Sans"]),
            )
        )

    return tri, axis


def rotate_vectorized_tri(svg_path, scale_x, scale_y):
    # Define SVG namespace
    ns = {"svg": "http://www.w3.org/2000/svg"}

    # Parse the SVG file
    tree = ET.parse(svg_path)
    root = tree.getroot()

    # Find all <g> elements with id="PolyCollection_1"
    g_elements = root.find(".//svg:g[@id='PolyCollection_1']", namespaces=ns)

    if g_elements is not None:
        # Apply the rotation transform to the group
        scale_factor = 1 / math.sqrt(2)
        transform = f"rotate(45 0 0) translate({scale_x}, {scale_y}) scale({scale_factor}, {scale_factor})"
        g_elements.set("transform", transform)

        # Save the modified SVG

    viewBox = root.get("viewBox")

    if viewBox:
        min_x, min_y, width, height = map(float, viewBox.split())
        new_min_y = min_y + height / 2  # Move down by half the height
        new_height = height / 2  # Reduce height by half
        root.set("viewBox", f"{min_x} {new_min_y} {width} {new_height}")
    else:
        print("No viewBox found. Consider adding one manually.")

    # Hacky, but it works to halve the height
    height_svg = root.get("height")
    if height_svg:
        current_height = re.match(r"(\d*\.?\d+)([a-zA-Z%]*)", height_svg)
        if current_height:
            numeric_height, unit = current_height.groups()
            numeric_height = float(numeric_height)

            root.set("height", f"{numeric_height / 1.8}pt")
    tree.write(svg_path)


def rotate_rasterized_tri(svg_path, shift_x, shift_y):
    # Load the SVG file
    tree = ET.parse(svg_path)
    root = tree.getroot()

    # Namespace handling
    ns = {"svg": "http://www.w3.org/2000/svg"}
    # Find all image elements with base64 embedded data
    for image in root.findall(".//svg:image", ns):
        href = image.get("{http://www.w3.org/1999/xlink}href", "")
        if href.startswith("data:image/png;base64,"):
            # Get the current width and height of the image
            width = float(image.get("width", 0)) / math.sqrt(2)
            height = float(image.get("height", 0)) / math.sqrt(2)
            # Set the new width and height
            image.set("width", str(width))
            image.set("height", str(height))
            # Apply a 270-degree rotation (about the top-left corner of the image)
            transform = image.get("transform", "")
            new_transform = (
                f"rotate(45, 0, 0) translate({shift_x}, {shift_y}) {transform}"
                if transform
                else f"rotate(45, 0, {height}) translate({shift_x}, {shift_y})"
            )
            image.set("transform", new_transform)
    # Update viewbox
    viewBox = root.get("viewBox")

    if viewBox:
        min_x, min_y, width, height = map(float, viewBox.split())
        new_min_y = min_y + height / 2  # Move down by half the height
        new_height = height / 2  # Reduce height by half
        root.set("viewBox", f"{min_x} {new_min_y} {width} {new_height}")
    else:
        print("No viewBox found. Consider adding one manually.")

    # Hacky, but it works to halve the height
    height_svg = root.get("height")
    if height_svg:
        current_height = re.match(r"(\d*\.?\d+)([a-zA-Z%]*)", height_svg)
        if current_height:
            numeric_height, unit = current_height.groups()
            numeric_height = float(numeric_height)

            root.set("height", f"{numeric_height / 1.8}pt")
    else:
        print("Warning: Could not parse height attribute.")

    # Save the modified SVG back to the same file
    tree.write(svg_path)


def make_tri_axis(sdf, title_name, palette, palette_orientation, colors, breaks, xlim):
    if not breaks:
        breaks = True
    else:
        breaks = [float(number) for number in breaks]
    if not xlim:
        xlim = 0
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
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max(), xlim)
    window = max(sdf["q_en"] - sdf["q_st"])
    if max_val < 100000:
        x_label = "Genomic Position (Kbp)"
    elif max_val < 100000000:
        x_label = "Genomic Position (Mbp)"
    else:
        x_label = "Genomic Position (Gbp)"
    p = (
        ggplot(sdf)
        + geom_tile(
            aes(x="q_st", y="r_st", fill="discrete", height=window, width=window),
            alpha=0,
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
            axis_line=element_line(color="black"),  # Adjust axis line size
            axis_text=element_text(
                family=["DejaVu Sans"]
            ),  # Change axis text font and size
            axis_ticks_major=element_line(),
            axis_line_x=element_line(),  # Keep the x-axis line
            axis_line_y=element_blank(),  # Remove the y-axis line
            axis_ticks_major_x=element_line(),  # Keep x-axis ticks
            axis_ticks_major_y=element_blank(),  # Remove y-axis ticks
            axis_text_x=element_line(),  # Keep x-axis text
            axis_text_y=element_blank(),
            plot_title=element_blank(),
        )
        + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + coord_fixed(ratio=1)
        + labs(x="", y="", title=title_name)
    )

    # Adjust x-axis label size
    p += theme(axis_title_x=element_text())

    return p


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
        function_name = getattr(diverging, "Spectral_11")
        palette_orientation = "-"
        hexcodes = function_name.hex_colors

    if palette_orientation == "-":
        new_hexcodes = hexcodes[::-1]
    else:
        new_hexcodes = hexcodes

    if custom_colors:
        new_hexcodes = custom_colors
    try:
        bot = np.quantile(sdf["perID_by_events"], q=0.001)
    except IndexError:
        bot = 0
    count = sdf.shape[0]
    extra = ""

    if count > 1e6:
        extra = "\n(thousands)"

    p = (
        ggplot(data=sdf, mapping=aes(x="perID_by_events", fill="discrete"))
        + geom_histogram(bins=300)
        + scale_color_cmap(cmap_name="plasma")
        + scale_fill_manual(new_hexcodes)
        + theme_light()
        + theme(text=element_text(family=["DejaVu Sans"]))
        + theme(legend_position="none")
        + coord_cartesian(xlim=(bot, 100))
        + xlab("% Identity Estimate")
        + ylab("# of Estimates{}".format(extra))
    )
    return p


def create_grid(
    singles,
    doubles,
    directory,
    palette,
    palette_orientation,
    single_names,
    double_names,
    is_freq,
    xlim,
    custom_colors,
    custom_breakpoints,
    axes_label,
    is_bed,
    width,
    breaks,
    deraster,
    vector_format,
):
    new_index = []
    transpose_index = []
    check_pascal(singles, doubles)
    # Singles can be empty if not selected
    for i in range(len(single_names)):
        for j in range(i + 1, len(single_names)):
            try:
                index = double_names.index([single_names[i], single_names[j]])
                transpose_index.append(0)
            except:
                index = double_names.index([single_names[j], single_names[i]])
                transpose_index.append(1)

            new_index.append(index)

    single_list = []
    double_list = []
    single_heatmap_list = []
    normal_heatmap_list = []
    transpose_heatmap_list = []
    for matrix in singles:
        if is_bed:
            df = read_df(
                None,
                palette,
                palette_orientation,
                is_freq,
                custom_colors,
                custom_breakpoints,
                matrix,
            )
        else:
            df = read_df(
                [matrix],
                palette,
                palette_orientation,
                is_freq,
                custom_colors,
                custom_breakpoints,
                None,
            )
        single_list.append(df)
    for matrix in doubles:
        if is_bed:
            df = read_df(
                None,
                palette,
                palette_orientation,
                is_freq,
                custom_colors,
                custom_breakpoints,
                matrix,
            )
        else:
            df = read_df(
                [matrix],
                palette,
                palette_orientation,
                is_freq,
                custom_colors,
                custom_breakpoints,
                None,
            )
        double_list.append(df)
    # This is the diagonals
    for plot in single_list:
        heatmap = make_dot_final(
            sdf=plot,
            width=width,
            palette=palette,
            palette_orientation=palette_orientation,
            colors=custom_colors,
            breaks=axes_label,
            xlim=xlim,
            transpose=False,
            deraster=deraster,
        )
        single_heatmap_list.append(heatmap)
    for indie in new_index:
        xd = new_index.index(indie)
        # These are the non-diagonal transposes and normals!
        if transpose_index[xd] == 0:
            heatmap = make_dot_final(
                sdf=double_list[indie],
                width=width,
                palette=palette,
                palette_orientation=palette_orientation,
                colors=custom_colors,
                breaks=axes_label,
                xlim=xlim,
                transpose=True,
                deraster=deraster,
            )
            normal_heatmap_list.append(heatmap)
            heatmap_t = make_dot_final(
                double_list[indie],
                width,
                palette,
                palette_orientation,
                custom_colors,
                axes_label,
                xlim,
                transpose=False,
                deraster=deraster,
            )
            transpose_heatmap_list.append(heatmap_t)
        else:
            heatmap = make_dot_final(
                double_list[indie],
                width,
                palette,
                palette_orientation,
                custom_colors,
                axes_label,
                xlim,
                transpose=True,
                deraster=deraster,
            )
            normal_heatmap_list.append(heatmap)
            heatmap_t = make_dot_final(
                double_list[indie],
                width,
                palette,
                palette_orientation,
                custom_colors,
                axes_label,
                xlim,
                transpose=False,
                deraster=deraster,
            )
            transpose_heatmap_list.append(heatmap_t)

    assert len(transpose_heatmap_list) == len(normal_heatmap_list)
    single_length = len(single_heatmap_list)
    if single_length == 0:
        single_length = reverse_pascal(len(normal_heatmap_list))

    normal_counter = 0
    trans_counter = 0
    trans_to_use = transpose_order(normal_heatmap_list)
    start_grid = pw.Brick(figsize=(9, 9))
    n = single_length * single_length

    if n > 9:
        print(f"This might take a while\n...\n")

    printProgressBar(0, n, prefix="Progress:", suffix="Complete", length=40)
    tots = 0
    col_names = pw.Brick(figsize=(width / 4.5, width))
    row_names = pw.Brick(figsize=(width, width / 4.5))
    for i in range(single_length):
        row_grid = pw.Brick(figsize=(width, width))
        for j in range(single_length):
            if i == j:
                if len(single_heatmap_list) == 0:
                    g1 = pw.Brick(figsize=(width, width))
                else:
                    g1 = pw.load_ggplot(single_heatmap_list[i], figsize=(width, width))

            elif i < j:
                g1 = pw.load_ggplot(
                    normal_heatmap_list[normal_counter], figsize=(width, width)
                )
                normal_counter += 1
            elif i > j:
                g1 = pw.load_ggplot(
                    transpose_heatmap_list[trans_to_use[trans_counter]],
                    figsize=(width, width),
                )
                trans_counter += 1
            if j == 0:
                row_grid = g1
            else:
                row_grid = row_grid | g1
            tots += 1
            printProgressBar(tots, n, prefix="Progress:", suffix="Complete", length=40)
        if i == 0:
            start_grid = row_grid
        else:
            start_grid = row_grid / start_grid
    for w in range(single_length):
        p1 = (
            ggplot()
            + geom_blank()
            + annotate(  # Use geom_blank to create a plot with no data
                "text",
                x=0,
                y=0,
                label=single_names[w],
                size=width * 3,
                angle=90,
                ha="center",
                va="center",
            )
            + theme(
                # Center the plot area and make backgrounds transparent
                axis_title_x=element_blank(),
                axis_title_y=element_blank(),
                axis_ticks=element_blank(),
                axis_text=element_blank(),
                plot_background=element_rect(
                    fill="none"
                ),  # Transparent plot background
                panel_background=element_rect(
                    fill="none"
                ),  # Transparent panel background
                panel_grid=element_blank(),
                aspect_ratio=0.5,  # Adjust the aspect ratio for the desired width/height
            )
            + coord_flip()  # Rotate the plot 90 degrees counterclockwise
        )
        p2 = (
            ggplot()
            + geom_blank()
            + annotate(  # Use geom_blank to create a plot with no data
                "text",
                x=0,
                y=0,
                label=single_names[w],
                size=width * 3,
                ha="center",
                va="center",
            )
            + theme(
                # Center the plot area and make backgrounds transparent
                axis_title_x=element_blank(),
                axis_title_y=element_blank(),
                axis_ticks=element_blank(),
                axis_text=element_blank(),
                plot_background=element_rect(
                    fill="none"
                ),  # Transparent plot background
                panel_background=element_rect(
                    fill="none"
                ),  # Transparent panel background
                panel_grid=element_blank(),
                aspect_ratio=0.5,  # Adjust the aspect ratio for the desired width/height
            )
        )
        g1 = pw.load_ggplot(p1, figsize=(width / 4.5, width))
        g2 = pw.load_ggplot(p2, figsize=(width, width / 4.5))

        if w == 0:
            col_names = g1
            row_names = g2
        else:
            col_names = g1 / col_names
            row_names = row_names | g2
    # Create a ghost 2x2
    pghost = (
        ggplot()
        + geom_blank()
        + annotate(  # Use geom_blank to create a plot with no data
            "text", x=0, y=0, label="", size=32, ha="center", va="center"
        )
        + theme(
            # Center the plot area and make backgrounds transparent
            axis_title_x=element_blank(),
            axis_title_y=element_blank(),
            axis_ticks=element_blank(),
            axis_text=element_blank(),
            plot_background=element_rect(fill="none"),  # Transparent plot background
            panel_background=element_rect(fill="none"),  # Transparent panel background
            panel_grid=element_blank(),
            aspect_ratio=0.5,  # Adjust the aspect ratio for the desired width/height
        )
    )
    ghosty = pw.load_ggplot(pghost, figsize=(2, 2))
    col_names = ghosty / col_names
    start_grid = col_names | (row_names / start_grid)
    gridname = f"{single_length}x{single_length}_GRID"
    print(f"\nGrid complete! Saving to {directory}/{gridname}...\n")
    start_grid.savefig(f"{directory}/{gridname}.png")
    start_grid.savefig(f"{directory}/{gridname}.{vector_format}", format=vector_format)
    print(f"Grid saved successfully!\n")


def create_plots(
    sdf,
    directory,
    name_x,
    name_y,
    palette,
    palette_orientation,
    no_hist,
    width,
    dpi,
    is_freq,
    xlim,
    custom_colors,
    custom_breakpoints,
    from_file,
    is_pairwise,
    axes_labels,
    axes_tick_number,
    vector_format,
    deraster,
):
    df = read_df(
        sdf,
        palette,
        palette_orientation,
        is_freq,
        custom_colors,
        custom_breakpoints,
        from_file,
    )
    sdf = df
    plot_filename = f"{directory}/{name_x}"
    if is_pairwise:
        plot_filename = f"{directory}/{name_x}_{name_y}"

    histy = make_hist(
        sdf, palette, palette_orientation, custom_colors, custom_breakpoints
    )

    if is_pairwise:
        heatmap = make_dot(
            sdf,
            plot_filename,
            palette,
            palette_orientation,
            custom_colors,
            axes_labels,
            axes_tick_number,
            xlim,
            deraster,
            width,
        )
        print(f"Creating plots and saving to {plot_filename}...\n")
        ggsave(
            heatmap,
            width=width,
            height=width,
            dpi=dpi,
            format=vector_format,
            filename=f"{plot_filename}_COMPARE.{vector_format}",
            verbose=False,
        )
        ggsave(
            heatmap,
            width=width,
            height=width,
            dpi=dpi,
            format="png",
            filename=f"{plot_filename}_COMPARE.png",
            verbose=False,
        )
        try:
            if not heatmap.data:
                print(
                    f"{plot_filename}_COMPARE.pdf and {plot_filename}_COMPARE.png saved sucessfully. \n"
                )
                return 0
        except ValueError:
            print(
                f"{plot_filename}_COMPARE.pdf and {plot_filename}_COMPARE.png saved sucessfully. \n"
            )
            return 0
        if no_hist:
            print(
                f"{plot_filename}_COMPARE.pdf and {plot_filename}_COMPARE.png saved sucessfully. \n"
            )
        else:
            ggsave(
                histy,
                width=3,
                height=3,
                dpi=dpi,
                format=vector_format,
                filename=f"{plot_filename}_HIST.{vector_format}",
                verbose=False,
            )
            ggsave(
                histy,
                width=3,
                height=3,
                dpi=dpi,
                format="png",
                filename=f"{plot_filename}_HIST.png",
                verbose=False,
            )
            print(
                f"{plot_filename}_COMPARE.pdf, {plot_filename}_COMPARE.png, {plot_filename}_HIST.pdf and {plot_filename}_HIST.png saved sucessfully. \n"
            )
    # Self-identity plots: Output _TRI, _FULL, and _HIST
    else:
        if deraster:
            print(f"Derasterization turned off. This may take a while...\n")
        tri_plot = make_tri(
            sdf,
            plot_filename,
            palette,
            palette_orientation,
            custom_colors,
            axes_labels,
            xlim,
            axes_tick_number,
            deraster,
            width,
        )
        full_plot = make_dot(
            check_st_en_equality(sdf),
            plot_filename,
            palette,
            palette_orientation,
            custom_colors,
            axes_labels,
            axes_tick_number,
            xlim,
            deraster,
            width,
        )
        ggsave(
            full_plot,
            width=width,
            height=width,
            dpi=dpi,
            format=vector_format,
            filename=f"{plot_filename}_FULL.{vector_format}",
            verbose=False,
        )
        ggsave(
            full_plot,
            width=width,
            height=width,
            dpi=dpi,
            format="png",
            filename=f"{plot_filename}_FULL.png",
            verbose=False,
        )
        tri_prefix = f"{plot_filename}_TRI"
        ggsave(
            tri_plot[0],
            width=width,
            height=width,
            dpi=dpi,
            format="svg",
            filename=f"{tri_prefix}.svg",
            verbose=False,
        )

        # These scaling values were determined thorugh much trial and error. Please don't delete :)
        if deraster:
            scaling_values = (46.62 * width, -3.75 * width)
            rotate_vectorized_tri(
                f"{tri_prefix}.svg", scaling_values[0], scaling_values[1]
            )
            try:
                cairosvg.svg2png(
                    url=f"{tri_prefix}.svg", write_to=f"{tri_prefix}.png", dpi=dpi
                )
            except:
                print(f"Error installing cairosvg. Unable to convert svg file. \n")
        else:
            scaling_values = (44.6 * width, -23 * width)
            rotate_rasterized_tri(
                f"{tri_prefix}.svg", scaling_values[0], scaling_values[1]
            )
            try:
                cairosvg.svg2png(
                    url=f"{tri_prefix}.svg", write_to=f"{tri_prefix}.png", dpi=dpi
                )
            except:
                print(f"Error installing cairosvg. Unable to convert svg file. \n")

        # Convert from svg to selected vector format. Ignore error if user has issues with cairosvg.
        try:
            if vector_format != "svg":
                if vector_format == "pdf":
                    cairosvg.svg2pdf(
                        url=f"{tri_prefix}.svg", write_to=f"{tri_prefix}.pdf"
                    )
                    if os.path.exists(f"{tri_prefix}.svg"):
                        os.remove(f"{tri_prefix}.svg")
                elif vector_format == "ps":
                    cairosvg.svg2pdf(
                        url=f"{tri_prefix}.svg", write_to=f"{tri_prefix}.ps"
                    )
                    if os.path.exists(f"{tri_prefix}.svg"):
                        os.remove(f"{tri_prefix}.svg")
        except:
            pass

        if no_hist:
            print(
                f"Triangle plots and full plots for {plot_filename} saved sucessfully. \n"
            )
        else:
            ggsave(
                histy,
                width=3,
                height=3,
                dpi=dpi,
                format=vector_format,
                filename=plot_filename + f"_HIST.{vector_format}",
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
                f"Triangle plots, full plots, and histogram for {plot_filename} saved sucessfully. \n"
            )
