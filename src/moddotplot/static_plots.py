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
)
import pandas as pd
import numpy as np
from PIL import Image
import patchworklib as pw
import math
import os

from moddotplot.const import (
    DIVERGING_PALETTES,
    QUALITATIVE_PALETTES,
    SEQUENTIAL_PALETTES,
)
from typing import List
from palettable.colorbrewer import qualitative, sequential, diverging


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
    bot = math.floor(min(sdf["perID_by_events"]))
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
    window = max(df["q_en"] - df["q_st"])

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
    bot = np.quantile(sdf["perID_by_events"], q=0.001)
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
    axes_label,
):
    print(singles)


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
        if no_hist:
            print(f"{plot_filename}_COMPARE.pdf and {plot_filename}_COMPARE.png saved sucessfully. \n")
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
