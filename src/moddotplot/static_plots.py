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
import svgutils.transform as sg
import cairosvg
import pandas as pd
import numpy as np
import glob
from PIL import Image
import patchworklib as pw
import math
import os
import xml.etree.ElementTree as ET
import sys
import re
from moddotplot.parse_fasta import printProgressBar
from lxml import etree
from pygenometracks.utilities import get_region
import matplotlib.pyplot as plt
from moddotplot.const import (
    DIVERGING_PALETTES,
    QUALITATIVE_PALETTES,
    SEQUENTIAL_PALETTES,
)
from typing import List
from palettable.colorbrewer import qualitative, sequential, diverging
import logging

# Set log level BEFORE importing pygenometracks
for name in logging.root.manager.loggerDict:
    if name.startswith("pygenometracks"):
        logging.getLogger(name).setLevel(logging.CRITICAL)
        logging.getLogger(name).propagate = False  # Don't pass to root logger

# Also make sure the root logger isn’t outputting debug messages
logging.basicConfig(level=logging.CRITICAL)

from pygenometracks.tracksClass import PlotTracks


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


def generate_ini_file(
    bedfile, ininame, chrom, color_value="bed_rgb", x_axis=True, display="collapsed"
):
    try:
        thing = chrom.split(":")[0]
        sections = [
            "[spacer]",
            "# height of space in cm (optional)",
            "height = 0.5",
            "",
            f"[{thing}]",
            f"file = {bedfile}",
            f"Title=",
            "height = 1",
            f"display = {display}",
            f"color = {color_value}",
            "labels = false",
            "fontsize = 10",
            "file_type = bed",
        ]

        if x_axis:
            sections.insert(0, "[x-axis]")

        ini_content = "\n".join(sections)
        with open(f"{ininame}.ini", "w") as file:
            file.write(ini_content)
        print(f"Successfully generated {ininame}.ini\n")
        return f"{ininame}.ini"
    except Exception as err:
        print(f"Error producing ini file: {err}\n")
        return None


def read_annotation_bed(filepath):
    """Reads a BED file into a Pandas DataFrame and ensures correct formatting."""
    col_names = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
        "thickStart",
        "thickEnd",
        "itemRgb",
    ]  # Include additional fields

    df = pd.read_csv(filepath, sep="\t", comment="#", header=None)

    # Ensure at least three required columns exist
    if df.shape[1] < 3:
        raise ValueError(
            "Invalid BED file: must have at least 3 columns (chrom, start, end)."
        )

    # Rename only the expected columns
    df.columns = col_names[: df.shape[1]]

    # Ensure start and end columns contain valid integers
    if df["start"].isna().any() or df["end"].isna().any():
        raise ValueError(
            "Invalid BED file: 'start' and 'end' columns must be integers and contain no missing values."
        )

    return df


def make_svg_background_transparent(svg_path, output_path=None):
    """
    Makes the background of an SVG file transparent by removing/modifying background fills.

    Args:
        svg_path: Path to input SVG file
        output_path: Path to output SVG file (if None, overwrites input)
    """
    import xml.etree.ElementTree as ET
    import re

    if output_path is None:
        output_path = svg_path

    # Parse the SVG
    tree = ET.parse(svg_path)
    root = tree.getroot()

    # Define SVG namespace
    ns = {"svg": "http://www.w3.org/2000/svg"}

    # Remove background rectangles/paths that cover the entire canvas
    # Get SVG dimensions for comparison
    width = root.get("width", "0")
    height = root.get("height", "0")

    # Extract numeric values
    width_num = float(re.sub(r"[a-zA-Z%]+", "", width)) if width != "0" else 0
    height_num = float(re.sub(r"[a-zA-Z%]+", "", height)) if height != "0" else 0

    # Find and modify background elements
    elements_to_modify = []

    # Check all paths, rectangles, and other elements
    for elem in root.iter():
        if (
            elem.tag.endswith("path")
            or elem.tag.endswith("rect")
            or elem.tag.endswith("polygon")
        ):
            # Check if this element has a background-like fill
            style = elem.get("style", "")
            fill = elem.get("fill", "")

            # Look for background colors (light colors, white, etc.)
            background_colors = [
                "#ffffff",
                "#f0ffff",
                "white",
                "lightblue",
                "lightgray",
                "lightgrey",
            ]

            is_background = False
            current_fill = None

            if "fill:" in style:
                # Extract fill from style
                fill_match = re.search(r"fill:\s*([^;]+)", style)
                if fill_match:
                    current_fill = fill_match.group(1).strip()
            elif fill:
                current_fill = fill

            if current_fill and any(
                bg_color in current_fill.lower() for bg_color in background_colors
            ):
                is_background = True

            # For paths, check if it covers a large area (likely background)
            if elem.tag.endswith("path"):
                d = elem.get("d", "")
                # Simple heuristic: if path starts at 0,0 and covers large area, it's likely background
                if "M 0" in d and current_fill:
                    is_background = True

            # For rectangles, check if it covers the full canvas
            if elem.tag.endswith("rect"):
                x = float(elem.get("x", 0))
                y = float(elem.get("y", 0))
                w = float(elem.get("width", 0))
                h = float(elem.get("height", 0))

                # If rectangle covers most/all of the canvas, it's likely background
                if x <= 1 and y <= 1 and w >= width_num * 0.9 and h >= height_num * 0.9:
                    is_background = True

            if is_background:
                elements_to_modify.append(elem)

    # Modify the background elements
    for elem in elements_to_modify:
        style = elem.get("style", "")

        if "fill:" in style:
            # Replace fill in style
            new_style = re.sub(r"fill:\s*[^;]+", "fill: transparent", style)
            elem.set("style", new_style)
        elif elem.get("fill"):
            # Replace fill attribute
            elem.set("fill", "transparent")

    # Save the modified SVG
    tree.write(output_path, encoding="unicode", xml_declaration=True)


def make_all_svg_backgrounds_transparent(directory):
    """
    Makes all SVG files in a directory have transparent backgrounds.
    """

    svg_files = glob.glob(os.path.join(directory, "*.svg"))

    for svg_file in svg_files:
        make_svg_background_transparent(svg_file)

    print(f"Processed {len(svg_files)} SVG files")


def run_pygenometracks(inifile, region, output_file, width):
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create directory if it doesn't exist

    try:
        trp = PlotTracks(
            inifile,
            width,
            fig_height=None,
            fontsize=10,
            dpi=300,
            track_label_width=0.05,
            plot_regions=region,
            plot_width=width,
        )

        # Extract the chromosome, start, and end from the region
        chrom, start, end = region[0]

        # Call the plot method directly to generate the image
        fig = trp.plot(output_file, chrom, start, end)

        return trp

    except Exception as e:
        if "No valid intervals were found" in str(e):
            print(f"No valid intervals found in BED file for region {region[0]}")
            print(
                "This is expected when the BED file doesn't overlap with the query region."
            )
            return None
        else:
            print(f"Error in run_pygenometracks: {e}")
            raise e


def test_pygenometracks_direct(inifile, chrom, start, end, output_file, width=40):
    """
    Direct test function to call PlotTracks.plot() without any coordinate validation.
    This bypasses the region checking that happens during PlotTracks initialization.

    Args:
        inifile: Path to the tracks configuration file
        chrom: Chromosome name
        start: Start coordinate
        end: End coordinate
        output_file: Output image file path
        width: Figure width in cm
    """

    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Create a dummy region for initialization (this gets overridden in plot())
    dummy_region = [(chrom, 1, 1000)]

    trp = PlotTracks(
        inifile,
        width,
        fig_height=None,
        fontsize=10,
        dpi=300,
        track_label_width=0.05,
        plot_regions=dummy_region,
        plot_width=width,
    )

    # Call plot() directly with your desired coordinates
    fig = trp.plot(output_file, chrom, start, end)

    return fig


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
        np.asarray(custom_breakpoints, dtype=np.float64)
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


def generate_breaks(min_number, max_number, min_breaks=5, max_breaks=9):
    # Determine the order of magnitude
    difference = max_number - min_number

    magnitude = 10 ** int(math.floor(math.log10(difference)))
    threshold = math.ceil(difference / magnitude)

    while threshold > max_breaks:
        magnitude *= 2
        threshold = math.ceil(difference / magnitude)

    while threshold < min_breaks:
        magnitude /= 2
        threshold = math.ceil(difference / magnitude)

    # Round down min_number to the nearest multiple of magnitude
    min_aligned = int(min_number // magnitude * magnitude)

    # Generate breakpoints
    upper_bound = int(min_aligned + (threshold + 1) * magnitude)
    breaks = list(range(min_aligned, upper_bound, int(magnitude)))

    return breaks


def make_dot(
    sdf,
    name_x,
    name_y,
    palette,
    palette_orientation,
    colors,
    breaks,
    num_ticks,
    xlim,
    deraster,
    width,
    is_pairwise,
):
    if is_pairwise:
        title_name = f"Comparative Plot: {name_x} vs {name_y}"
    else:
        title_name = f"Self-Identity Plot: {name_x}"
    title_length = 2 * width
    if len(title_name) > 50:
        title_length = 1.5 * width
    elif len(title_name) > 80:
        title_length = width
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
    min_val = min(sdf["q_st"].min(), sdf["r_st"].min())
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max(), xlim)

    # If user provides breaks, convert to ints
    if not breaks:
        breaks = generate_breaks(int(min_val), int(max_val))
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
        title=element_text(
            family=["DejaVu Sans"], size=title_length, hjust=0.5
        ),  # Center title
        axis_title_x=element_text(size=(width * 1.4), family=["DejaVu Sans"]),
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
        + scale_x_continuous(
            labels=make_scale, limits=[min_val, max_val], breaks=breaks
        )
        + scale_y_continuous(
            labels=make_scale, limits=[min_val, max_val], breaks=breaks
        )
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
    min_val = max(sdf["q_st"].min(), sdf["r_st"].min())
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max(), xlim)

    # If user provides breaks, convert to ints
    if not breaks:
        breaks = generate_breaks(int(min_val), int(max_val))
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
    min_val = min(sdf["q_st"].min(), sdf["r_st"].min())
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max(), xlim)

    # If user provides breaks, convert to ints
    if not breaks:
        breaks = generate_breaks(int(min_val), int(max_val))
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
    min_val = max(sdf["q_st"].min(), sdf["r_st"].min())
    max_val = max(sdf["q_en"].max(), sdf["r_en"].max(), xlim)

    # If user provides breaks, convert to ints
    if not breaks:
        breaks = generate_breaks(int(min_val), int(max_val))
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
            + scale_x_continuous(
                labels=make_scale, limits=[min_val, max_val], breaks=breaks
            )
            + scale_y_continuous(
                labels=make_scale, limits=[min_val, max_val], breaks=breaks
            )
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
                axis_title_x=element_text(size=(width * 1.4), family=["DejaVu Sans"]),
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
            + scale_x_continuous(
                labels=make_scale, limits=[min_val, max_val], breaks=breaks
            )
            + scale_y_continuous(
                labels=make_scale, limits=[min_val, max_val], breaks=breaks
            )
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
            + scale_x_continuous(
                labels=make_scale, limits=[min_val, max_val], breaks=breaks
            )
            + scale_y_continuous(
                labels=make_scale, limits=[min_val, max_val], breaks=breaks
            )
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


from lxml import etree


def append_svg(svg1_path, svg2_path, output_path):
    """Appends SVG2 to the bottom of SVG1 without modifying its width."""
    # Load SVG1 and SVG2
    tree1 = etree.parse(svg1_path)
    root1 = tree1.getroot()

    tree2 = etree.parse(svg2_path)
    root2 = tree2.getroot()

    # Extract width and height of SVG1
    width1 = float(root1.get("width", "0").replace("pt", ""))
    height1 = float(root1.get("height", "0").replace("pt", ""))

    # Extract width and height of SVG2
    width2 = float(root2.get("width", "0").replace("pt", ""))
    height2 = float(root2.get("height", "0").replace("pt", ""))

    # Update SVG1 height to accommodate SVG2
    new_height = height1 + height2
    root1.set("height", f"{new_height}pt")

    # Create a translation group for SVG2 and shift it down
    group = etree.Element("g", attrib={"transform": f"translate(0,{height1 + 400})"})
    for child in root2:
        group.append(child)

    # Append translated SVG2 to SVG1
    root1.append(group)

    # Save the new merged SVG
    tree1.write(output_path, pretty_print=True, xml_declaration=True, encoding="utf-8")


def get_svg_size(svg_path):
    """Helper to extract width and height of an SVG in pt units."""
    import xml.etree.ElementTree as ET

    tree = ET.parse(svg_path)
    root = tree.getroot()
    width = root.get("width")
    height = root.get("height")
    return width, height


def parse_size(size_str):
    """Convert '648pt' or '800px' -> float(648)."""
    return float(re.sub(r"[a-zA-Z]+", "", size_str))


def merge_annotation_tri(svg1_path, svg2_path, output_path, deraster, width):
    """Merges two SVG files into a single SVG file with proper size."""

    w1, h1 = get_svg_size(svg1_path)
    w2, h2 = get_svg_size(svg2_path)

    # Ensure they are floats
    w1, h1 = parse_size(w1), parse_size(h1)
    w2, h2 = parse_size(w2), parse_size(h2)
    # Determine total size
    total_width = max(w1, w2)
    total_height = h1 + h2
    # Create figure
    fig = sg.SVGFigure(f"{total_width}px", f"{total_height}px")

    # Load SVGs
    svg1 = sg.fromfile(svg1_path).getroot()
    svg2 = sg.fromfile(svg2_path).getroot()
    make_svg_background_transparent(svg2_path)

    # Position
    if deraster:
        # Its not perfect for width > 18, but good enough
        adjust_svg1 = (0, -5 * (h1 / 6))
        adjust_svg2 = ((9 - (width / 2)), 5 * (h1 / 6))
        svg1.moveto(adjust_svg1[0], adjust_svg1[1])
        svg2.scale(1.077 + (width / 1000))
        svg2.moveto(adjust_svg2[0], adjust_svg2[1])
    else:
        adjust_svg1 = (0, -5 * (h1 / 6))
        adjust_svg2 = (10 + (width / 4.5), 5 * (h1 / 6))
        svg1.moveto(adjust_svg1[0], adjust_svg1[1])
        svg2.moveto(adjust_svg2[0], adjust_svg2[1])
        scaling_factor = width / 2
        svg2.scale(1.034 + (scaling_factor / 1000))

    # Append and save
    fig.append([svg1, svg2])
    fig.set_size((f"{total_width}px", f"{total_height}px"))
    fig.save(output_path)


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
    annotation,
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

    plot_filename = os.path.join(directory, name_x)

    if is_pairwise:
        plot_filename = os.path.join(directory, f"{name_x}_{name_y}")

    histy = make_hist(
        sdf, palette, palette_orientation, custom_colors, custom_breakpoints
    )

    # Just doing triangle plots for now.
    if annotation:
        print(f"Generating ini file for annotation track:\n")
        iniprefix = plot_filename
        inifile = generate_ini_file(
            bedfile=annotation,
            ininame=iniprefix,
            chrom=name_x,
        )
        min_val = max(sdf["q_st"].min(), sdf["r_st"].min())
        max_val = max(sdf["q_en"].max(), sdf["r_en"].max())
        if not xlim:
            xlim = max_val
        region = [(name_x.split(":")[0], min_val, xlim)]
        if inifile:
            try:
                # Check if the BED file has valid intervals for this region first
                bed_df = read_annotation_bed(annotation)
                chrom_name = name_x.split(":")[0]

                # Filter for the chromosome and region of interest
                valid_intervals = bed_df[
                    (bed_df["chrom"] == chrom_name)
                    & (bed_df["end"] >= min_val)
                    & (bed_df["start"] <= xlim)
                ]

                if valid_intervals.empty:
                    print(
                        f"No valid intervals found in {annotation} for region {chrom_name}:{min_val}-{xlim}.\n"
                    )
                    print("Skipping annotation track generation.\n")
                else:
                    bed_track = run_pygenometracks(
                        inifile=inifile,
                        region=region,
                        output_file=f"{iniprefix}_ANNOTATION_TRACK.svg",
                        width=width * 2.05,
                    )
                    bed_track.plot(
                        f"{iniprefix}_ANNOTATION_TRACK.svg",
                        name_x.split(":")[0],
                        min_val,
                        xlim,
                    )
                    bed_track.plot(
                        f"{iniprefix}_ANNOTATION_TRACK.png",
                        name_x.split(":")[0],
                        min_val,
                        xlim,
                    )
                    print(f"\nAnnotation track saved to {iniprefix}_ANNOTATION_TRACK\n")

            except Exception as e:
                print(f"Error processing annotation file {annotation}: {e}\n")
                print("Skipping annotation track generation.\n")

    if is_pairwise:
        heatmap = make_dot(
            sdf,
            name_x,
            name_y,
            palette,
            palette_orientation,
            custom_colors,
            axes_labels,
            axes_tick_number,
            xlim,
            deraster,
            width,
            True,
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
        if not no_hist:
            ggsave(
                histy,
                width=3,
                height=3,
                dpi=dpi,
                format=vector_format,
                filename=f"{plot_filename}_COMPARE_HIST.{vector_format}",
                verbose=False,
            )
            ggsave(
                histy,
                width=3,
                height=3,
                dpi=dpi,
                format="png",
                filename=f"{plot_filename}_COMPARE_HIST.png",
                verbose=False,
            )
        try:
            if not heatmap.data:
                print(
                    f"{plot_filename} comparative plots and histogram saved sucessfully. \n"
                )
                return 0
        except ValueError:
            print(
                f"{plot_filename} comparative plots and histogram saved sucessfully. \n"
            )
            return 0
        if no_hist:
            print(
                f"{plot_filename}_COMPARE.{vector_format} and {plot_filename}_COMPARE.png saved sucessfully. \n"
            )
    # Self-identity plots: Output _TRI, _FULL, and _HIST
    else:
        if deraster:
            print(
                f"Producing dotplots with derasterization turned off. This may take a while...\n"
            )
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
            name_x,
            name_y,
            palette,
            palette_orientation,
            custom_colors,
            axes_labels,
            axes_tick_number,
            xlim,
            deraster,
            width,
            False,
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
        if annotation:
            anno_prefix = f"{plot_filename}_PRE_ANNOTATED"
            annotated_tri = tri_plot[0] + theme(
                axis_title_x=element_blank(),
                axis_line_x=element_blank(),
                axis_text_x=element_blank(),
                axis_ticks_minor_x=element_blank(),
                axis_ticks=element_blank(),
            )
            ggsave(
                annotated_tri,
                width=width,
                height=width,
                dpi=dpi,
                format="svg",
                filename=f"{anno_prefix}.svg",
                verbose=False,
            )
        # These scaling values were determined thorugh much trial and error. Please don't delete :)
        if deraster:
            scaling_values = (46.62 * width, -3.75 * width)
            rotate_vectorized_tri(
                f"{tri_prefix}.svg", scaling_values[0], scaling_values[1]
            )
            if annotation:
                rotate_vectorized_tri(
                    f"{anno_prefix}.svg", scaling_values[0], scaling_values[1]
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
            if annotation:
                if annotation:
                    rotate_rasterized_tri(
                        f"{anno_prefix}.svg", scaling_values[0], scaling_values[1]
                    )
            try:
                cairosvg.svg2png(
                    url=f"{tri_prefix}.svg", write_to=f"{tri_prefix}.png", dpi=dpi
                )
            except:
                print(f"Error installing cairosvg. Unable to convert svg file. \n")
        if annotation:
            # Only merge if annotation was successfully created
            if os.path.exists(f"{iniprefix}_ANNOTATION_TRACK.svg"):
                make_svg_background_transparent(f"{iniprefix}_ANNOTATION_TRACK.svg")
                merge_annotation_tri(
                    f"{anno_prefix}.svg",
                    f"{iniprefix}_ANNOTATION_TRACK.svg",
                    f"{tri_prefix}_ANNOTATED.svg",
                    deraster,
                    width,
                )
                if os.path.exists(f"{anno_prefix}.svg"):
                    os.remove(f"{anno_prefix}.svg")
                try:
                    if vector_format != "svg":
                        if vector_format == "pdf":
                            cairosvg.svg2pdf(
                                url=f"{tri_prefix}_ANNOTATED.svg",
                                write_to=f"{tri_prefix}_ANNOTATED.pdf",
                            )
                            cairosvg.svg2pdf(
                                url=f"{iniprefix}_ANNOTATION_TRACK.svg",
                                write_to=f"{iniprefix}_ANNOTATION_TRACK.pdf",
                            )
                        elif vector_format == "ps":
                            cairosvg.svg2ps(
                                url=f"{tri_prefix}_ANNOTATED.svg",
                                write_to=f"{tri_prefix}_ANNOTATED.ps",
                            )
                            cairosvg.svg2ps(
                                url=f"{tri_prefix}_ANNOTATION_TRACK.svg",
                                write_to=f"{tri_prefix}_ANNOTATION_TRACK.ps",
                            )
                        if os.path.exists(f"{iniprefix}_ANNOTATION_TRACK.svg"):
                            os.remove(f"{iniprefix}_ANNOTATION_TRACK.svg")
                        if os.path.exists(f"{tri_prefix}_ANNOTATED.svg"):
                            os.remove(f"{tri_prefix}_ANNOTATED.svg")
                except Exception as e:
                    print(f"Error converting annotated SVG: {e}")
            else:
                print("Annotation file not created, skipping merge step.")
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
