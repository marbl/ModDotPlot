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
    theme_void,
    geom_blank,
    annotate,
    element_rect,
    coord_flip,
    theme_minimal,
)
import pandas as pd
import numpy as np
from PIL import Image
import patchworklib as pw
import math
import os
import sys
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


# Hardcoding for now, I have the formula just lazy
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
    if scaled[-1] < 100000:
        return make_k(scaled)
    elif scaled[-1] > 100000000:
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
    try:
        window = max(df["q_en"] - df["q_st"])
    except ValueError:
        window = 0

    # Calculate the position of the first and second intervals
    df["first_pos"] = df["q_st"] / window
    df["second_pos"] = df["r_st"] / window

    return df


def make_dot(sdf, title_name, palette, palette_orientation, colors, breaks, xlim):
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
    try:
        window = max(sdf["q_en"] - sdf["q_st"])
    except:
        p = (
            ggplot(aes(x=[], y=[]))
            + theme_minimal()  # Use a minimal theme with gridlines
            + theme(
                panel_grid_major=element_blank(),  # Remove major gridlines (optional)
                panel_grid_minor=element_blank(),
            )  # Remove minor gridlines (optional)
        )
        return p
    if max_val < 100000:
        x_label = "Genomic Position (Kbp)"
    elif max_val < 100000000:
        x_label = "Genomic Position (Mbp)"
    else:
        x_label = "Genomic Position (Gbp)"
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
            axis_line=element_line(color="black"),  # Adjust axis line size
            axis_text=element_text(
                family=["Dejavu Sans"]
            ),  # Change axis text font and size
            axis_ticks_major=element_line(),
            title=element_text(
                family=["Dejavu Sans"],  # Change title font family
            ),
        )
        + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + coord_fixed(ratio=1)
        + facet_grid("r ~ q")
        + labs(x="", y="", title=title_name)
    )

    # Adjust x-axis label size
    # p += theme(axis_title_x=element_text())

    return p


def make_dot2(sdf, title_name, palette, palette_orientation, colors, breaks, xlim):
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
    try:
        window = max(sdf["q_en"] - sdf["q_st"])
    except:
        p = (
            ggplot(aes(x=[], y=[]))
            + theme_minimal()  # Use a minimal theme with gridlines
            + theme(
                panel_grid_major=element_blank(),  # Remove major gridlines (optional)
                panel_grid_minor=element_blank(),
            )  # Remove minor gridlines (optional)
        )
        return p
    if max_val < 100000:
        x_label = "Genomic Position (Kbp)"
    elif max_val < 100000000:
        x_label = "Genomic Position (Mbp)"
    else:
        x_label = "Genomic Position (Gbp)"
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
            axis_line=element_line(color="black"),  # Adjust axis line size
            axis_text=element_text(
                family=["Dejavu Sans"]
            ),  # Change axis text font and size
            axis_ticks_major=element_line(),
            title=element_text(
                family=["Dejavu Sans"],  # Change title font family
            ),
        )
        + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + coord_fixed(ratio=1)
        + labs(x=None, y=None, title=title_name)
    )

    # Adjust x-axis label size
    p += theme(axis_title_x=element_blank())
    p += theme(axis_title_y=element_blank())
    # p += theme(title=element_blank())

    return p


def make_dot3(sdf, title_name, palette, palette_orientation, colors, breaks, xlim):
    # sdf = sdf.transpose()
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
    try:
        window = max(sdf["q_en"] - sdf["q_st"])
    except:
        p = (
            ggplot(aes(x=[], y=[]))
            + theme_minimal()  # Use a minimal theme with gridlines
            + theme(
                panel_grid_major=element_blank(),  # Remove major gridlines (optional)
                panel_grid_minor=element_blank(),
            )  # Remove minor gridlines (optional)
        )
        return p
    if max_val < 100000:
        x_label = "Genomic Position (Kbp)"
    elif max_val < 100000000:
        x_label = "Genomic Position (Mbp)"
    else:
        x_label = "Genomic Position (Gbp)"
    p = (
        ggplot(sdf)
        + geom_tile(
            aes(x="r_st", y="q_st", fill="discrete", height=window, width=window)
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
                family=["Dejavu Sans"]
            ),  # Change axis text font and size
            axis_ticks_major=element_line(),
            title=element_text(
                family=["Dejavu Sans"],  # Change title font family
            ),
        )
        + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + coord_fixed(ratio=1)
        + labs(x=None, y=None, title=title_name)
    )

    # Adjust x-axis label size
    p += theme(axis_title_x=element_blank())
    p += theme(axis_title_y=element_blank())
    # p += theme(title=element_blank())

    return p


def make_tri(sdf, title_name, palette, palette_orientation, colors, breaks, xlim):
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
            axis_line=element_line(color="black"),  # Adjust axis line size
            axis_text=element_text(
                family=["DejaVu Sans"]
            ),  # Change axis text font and size
            axis_ticks_major=element_line(),
            title=element_text(
                family=["DejaVu Sans"],  # Change title font family
            ),
        )
        + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + coord_fixed(ratio=1)
        + labs(x="", y="", title=title_name)
    )

    # Adjust x-axis label size
    p += theme(axis_title_x=element_text())

    return p


def make_tri2(sdf, title_name, palette, palette_orientation, colors, breaks, xlim):
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
        function_name = getattr(sequential, "Spectral_11")
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
    no_hist,
    single_names,
    double_names,
    is_freq,
    xlim,
    custom_colors,
    custom_breakpoints,
    axes_label,
    is_bed,
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
    for plot in single_list:
        heatmap = make_dot2(
            plot,
            plot["q"].iloc[1],
            palette,
            palette_orientation,
            custom_colors,
            axes_label,
            xlim,
        )
        single_heatmap_list.append(heatmap)
    for indie in new_index:
        xd = new_index.index(indie)
        if transpose_index[xd] == 0:
            heatmap = make_dot2(
                double_list[indie],
                f"{double_names[indie][0]}_{double_names[indie][1]}",
                palette,
                palette_orientation,
                custom_colors,
                axes_label,
                xlim,
            )
            normal_heatmap_list.append(heatmap)
            heatmap_t = make_dot3(
                double_list[indie],
                f"{double_names[indie][1]}_{double_names[indie][0]}",
                palette,
                palette_orientation,
                custom_colors,
                axes_label,
                xlim,
            )
            transpose_heatmap_list.append(heatmap_t)
        else:
            heatmap = make_dot3(
                double_list[indie],
                f"{double_names[indie][0]}_{double_names[indie][1]}",
                palette,
                palette_orientation,
                custom_colors,
                axes_label,
                xlim,
            )
            normal_heatmap_list.append(heatmap)
            heatmap_t = make_dot2(
                double_list[indie],
                f"{double_names[indie][1]}_{double_names[indie][0]}",
                palette,
                palette_orientation,
                custom_colors,
                axes_label,
                xlim,
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

    if n > 4:
        print(f"This might take a while\n...\n")

    printProgressBar(0, n, prefix="Progress:", suffix="Complete", length=40)
    tots = 0
    col_names = pw.Brick(figsize=(2, 9))
    row_names = pw.Brick(figsize=(9, 2))
    for i in range(single_length):
        row_grid = pw.Brick(figsize=(9, 9))
        for j in range(single_length):
            if i == j:
                if len(single_heatmap_list) == 0:
                    g1 = pw.Brick(figsize=(9, 9))
                else:
                    g1 = pw.load_ggplot(single_heatmap_list[i], figsize=(9, 9))

            elif i < j:
                g1 = pw.load_ggplot(normal_heatmap_list[normal_counter], figsize=(9, 9))
                normal_counter += 1
            elif i > j:
                g1 = pw.load_ggplot(
                    transpose_heatmap_list[trans_to_use[trans_counter]], figsize=(9, 9)
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
                size=32,
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
                size=32,
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
        g1 = pw.load_ggplot(p1, figsize=(2, 9))
        g2 = pw.load_ggplot(p2, figsize=(9, 2))

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
    start_grid.savefig(f"{directory}/{gridname}.pdf", format="pdf")
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
            xlim,
        )
        print(f"Creating plots and saving to {plot_filename}...\n")
        ggsave(
            heatmap,
            width=9,
            height=9,
            dpi=dpi,
            format="pdf",
            filename=f"{plot_filename}_COMPARE.pdf",
            verbose=False,
        )
        ggsave(
            heatmap,
            width=9,
            height=9,
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
                format="pdf",
                filename=f"{plot_filename}_HIST.pdf",
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
        tri_plot = make_tri(
            sdf,
            plot_filename,
            palette,
            palette_orientation,
            custom_colors,
            axes_labels,
            xlim,
        )
        tri_plot_axis_only = make_tri2(
            sdf,
            plot_filename,
            palette,
            palette_orientation,
            custom_colors,
            axes_labels,
            xlim,
        )
        full_plot = make_dot(
            check_st_en_equality(sdf),
            plot_filename,
            palette,
            palette_orientation,
            custom_colors,
            axes_labels,
            xlim,
        )
        print(f"Creating plots and saving to {plot_filename}...\n")
        triplot_no_axis = tri_plot + theme(
            axis_text_x=element_blank(),
            axis_text_y=element_blank(),
            axis_title_x=element_blank(),
            axis_title_y=element_blank(),
            axis_line_x=element_blank(),
            axis_line_y=element_blank(),
            axis_ticks_major=element_blank(),
            axis_ticks_minor=element_blank(),
            panel_background=element_blank(),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            plot_title=element_blank(),
        )
        ggsave(
            triplot_no_axis,
            width=9,
            height=9,
            dpi=600,
            format="png",
            filename=f"{plot_filename}_TRI_NOAXIS.png",
            verbose=False,
        )
        ggsave(
            tri_plot_axis_only,
            width=9,
            height=9,
            dpi=600,
            format="png",
            filename=f"{plot_filename}_AXIS.png",
            verbose=False,
        )

        png_no_axes = Image.open(f"{plot_filename}_TRI_NOAXIS.png")
        rotated_png = png_no_axes.rotate(315, expand=True)

        rotated_png.save(f"{plot_filename}_ROTATED_TRI_NOAXIS.png")
        overlap_axis(rotated_png, f"{plot_filename}_AXIS.png", plot_filename)

        if os.path.exists(f"{plot_filename}_ROTATED_TRI_NOAXIS.png"):
            os.remove(f"{plot_filename}_ROTATED_TRI_NOAXIS.png")
        if os.path.exists(f"{plot_filename}_TRI_NOAXIS.png"):
            os.remove(f"{plot_filename}_TRI_NOAXIS.png")

        ggsave(
            full_plot,
            width=9,
            height=9,
            dpi=dpi,
            format="pdf",
            filename=f"{plot_filename}_FULL.pdf",
            verbose=False,
        )
        ggsave(
            full_plot,
            width=9,
            height=9,
            dpi=dpi,
            format="png",
            filename=f"{plot_filename}_FULL.png",
            verbose=False,
        )
        if no_hist:
            print(
                f"{plot_filename}_TRI.png, {plot_filename}_TRI.pdf, {plot_filename}_FULL.png and {plot_filename}_FULL.pdf saved sucessfully. \n"
            )
        else:
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
                f"{plot_filename}_TRI.png, {plot_filename}_TRI.pdf, {plot_filename}_FULL.png, {plot_filename}_FULL.pdf, {plot_filename}_HIST.png and {plot_filename}_HIST.pdf, saved sucessfully. \n"
            )
