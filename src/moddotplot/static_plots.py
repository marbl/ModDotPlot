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
    theme_minimal,
)
import pandas as pd
import numpy as np
import math
import os
import re
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgb
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter
from moddotplot.native_render import (
    DEFAULT_FONT_FAMILY,
    FALLBACK_FONT_FAMILY,
    MIN_TEXT_SIZE,
    MIN_TITLE_SIZE,
    clamped_font_size,
    configure_dotplot_axis,
    configure_triangle_axis,
    create_triangle_layout,
    draw_rectangular_tiles,
    draw_triangle_tiles,
    genomic_scale,
    genomic_tick_formatter,
    is_glyph_loading_error,
    save_figure_pair,
    save_with_font_fallback,
    set_figure_font_family,
)
from moddotplot.const import (
    DIRECTION_COLORS,
    DIVERGING_PALETTES,
    QUALITATIVE_PALETTES,
    SEQUENTIAL_PALETTES,
)
from palettable.colorbrewer import qualitative, sequential, diverging
from moddotplot.annotations import (
    DEFAULT_ANNOTATION_COLOR,
    annotation_color as _annotation_color,
    read_annotation_bed,
    visible_annotation_intervals as _visible_annotation_intervals,
)

REGION_SUFFIX_PATTERN = re.compile(r"(?::\d+-\d+)+$")


def _plot_font_theme(family=DEFAULT_FONT_FAMILY):
    """Apply one family to every Plotnine text themeable."""

    font = element_text(family=[family])
    return theme(
        text=font,
        title=element_text(family=[family]),
        axis_text=element_text(family=[family]),
        strip_text=element_text(family=[family]),
        legend_text=element_text(family=[family]),
    )


def _save_plot(plot, **kwargs):
    """Save a Plotnine plot in Helvetica, retrying on glyph-load failure."""

    try:
        ggsave(plot + _plot_font_theme(), **kwargs)
    except RuntimeError as error:
        if not is_glyph_loading_error(error):
            raise
        ggsave(plot + _plot_font_theme(FALLBACK_FONT_FAMILY), **kwargs)


def _draw_and_save_plot_pair(
    plot,
    output_prefix,
    *,
    width,
    height,
    dpi,
    vector_format,
):
    """Build one Plotnine figure and save both raster and vector outputs.

    ``ggsave`` redraws a plot for every requested format. Large tile plots and
    histograms therefore paid their complete scale/layout/rasterization cost
    twice. Drawing once also guarantees that both files contain the same axes,
    labels, and tile realization.
    """

    def draw(family):
        styled = (
            plot
            + _plot_font_theme(family)
            + theme(figure_size=(float(width), float(height)), dpi=int(dpi))
        )
        return styled.draw(show=False)

    try:
        figure = draw(DEFAULT_FONT_FAMILY)
    except RuntimeError as error:
        if not is_glyph_loading_error(error):
            raise
        figure = draw(FALLBACK_FONT_FAMILY)
    try:
        return save_figure_pair(
            figure,
            output_prefix,
            vector_format,
            dpi,
            # Match plotnine/ggsave's requested physical canvas exactly.  A
            # tight bounding box changes both the raster dimensions and plot
            # framing (for example, 3 in at 96 dpi no longer yields 288 px).
            bbox_inches=figure.bbox_inches,
        )
    finally:
        plt.close(figure)


def display_sequence_name(name):
    """Return a sequence name without appended region coordinates."""

    return REGION_SUFFIX_PATTERN.sub("", str(name))


def _fit_grid_sequence_labels(figure, axes):
    """Fit grid headings inside the existing figure canvas."""

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    ratios = []

    for axis in axes[0, :]:
        title = axis.title
        if title.get_text():
            title_box = title.get_window_extent(renderer=renderer)
            axis_box = axis.get_window_extent(renderer=renderer)
            if title_box.width:
                ratios.append(axis_box.width * 0.9 / title_box.width)

    for axis in axes[:, 0]:
        label = axis.yaxis.label
        if label.get_text():
            label_box = label.get_window_extent(renderer=renderer)
            axis_box = axis.get_window_extent(renderer=renderer)
            if label_box.height:
                ratios.append(axis_box.height * 0.9 / label_box.height)

    artists = [axis.title for axis in axes[0, :]] + [
        axis.yaxis.label for axis in axes[:, 0]
    ]
    artists = [artist for artist in artists if artist.get_text()]
    if not ratios or not artists:
        return
    scale = min(1.0, min(ratios))
    if scale >= 1.0:
        return

    for artist in artists:
        artist.set_fontsize(artist.get_fontsize() * scale)


def _resolve_native_colors(palette, palette_orientation, custom_colors=None):
    """Resolve plot colors with the same orientation rules as plotnine paths."""
    if palette in DIVERGING_PALETTES:
        palette_colors = getattr(diverging, palette).hex_colors
        palette_orientation = "-" if palette_orientation == "+" else "+"
    elif palette in QUALITATIVE_PALETTES:
        palette_colors = getattr(qualitative, palette).hex_colors
    elif palette in SEQUENTIAL_PALETTES:
        palette_colors = getattr(sequential, palette).hex_colors
    else:
        palette_colors = diverging.Spectral_11.hex_colors
        palette_orientation = "-"

    colors = palette_colors[::-1] if palette_orientation == "-" else palette_colors
    return list(custom_colors) if custom_colors else list(colors)


DIRECTION_ANI_COLUMN = "direction_ani"


def _direction_ani_style(dataframe):
    """Return data and colors for direction hue plus ANI intensity.

    Direction selects the blue or pink hue. The existing ordered ANI bins
    control saturation: weak matches are pale and the strongest bin reaches
    the base direction color. Data without ANI bins retains the solid legacy
    direction colors, which keeps the low-level rendering API usable.
    """

    if "direction" not in dataframe.columns:
        return dataframe, None, None
    if "discrete" not in dataframe.columns:
        return dataframe, DIRECTION_COLORS, "direction"

    categories = (
        list(dataframe["discrete"].cat.categories)
        if isinstance(dataframe["discrete"].dtype, pd.CategoricalDtype)
        else list(pd.unique(dataframe["discrete"].dropna()))
    )
    if not categories:
        return dataframe, DIRECTION_COLORS, "direction"

    colors = {}
    category_count = len(categories)
    for direction, base_color in DIRECTION_COLORS.items():
        base_rgb = np.asarray(to_rgb(base_color))
        for index, category in enumerate(categories):
            fraction = 1.0 if category_count == 1 else index / (category_count - 1)
            strength = 0.25 + (0.75 * fraction)
            rgb = np.ones(3) + ((base_rgb - np.ones(3)) * strength)
            colors[f"{direction}:{category}"] = to_hex(rgb)

    styled = dataframe.copy()
    styled[DIRECTION_ANI_COLUMN] = [
        f"{direction}:{category}"
        for direction, category in zip(styled["direction"], styled["discrete"])
    ]
    return styled, colors, DIRECTION_ANI_COLUMN


def is_plot_empty(p):
    # Check if the plot has data or any layers
    return len(p.layers) == 0 and p.data.empty


def draw_annotation_track(
    axis,
    bed_df,
    chrom,
    region_start,
    region_end,
    fallback=DEFAULT_ANNOTATION_COLOR,
):
    """Draw a label-free, collapsed BED track and return its interval count."""
    intervals = _visible_annotation_intervals(
        bed_df, chrom, region_start, region_end, fallback
    )
    for interval_start, interval_end, color in intervals:
        axis.add_patch(
            Rectangle(
                (interval_start, 0.2),
                interval_end - interval_start,
                0.6,
                facecolor=color,
                edgecolor=color,
                linewidth=0.25,
            )
        )

    axis.set_xlim(region_start, region_end)
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    axis.xaxis.set_major_formatter(formatter)
    axis.tick_params(axis="x", labelsize=8, length=3)
    axis.spines["left"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["top"].set_visible(False)
    axis.set_facecolor("none")
    return len(intervals)


def render_annotation_track(
    bed_df,
    chrom,
    region_start,
    region_end,
    output_prefix,
    width_cm,
    dpi,
    vector_format="svg",
):
    """Render a collapsed BED track as PNG and the selected vector format."""
    if vector_format not in {"svg", "pdf", "ps"}:
        raise ValueError(f"Unsupported vector format: {vector_format}")

    figure, axis = plt.subplots(figsize=(max(float(width_cm), 2.54) / 2.54, 2.0 / 2.54))
    try:
        interval_count = draw_annotation_track(
            axis, bed_df, chrom, region_start, region_end
        )
        if interval_count == 0:
            return False

        figure.subplots_adjust(left=0.02, right=0.995, bottom=0.34, top=0.96)

        def save_outputs():
            figure.savefig(
                f"{output_prefix}.{vector_format}",
                format=vector_format,
                dpi=dpi,
                transparent=vector_format != "ps",
                facecolor="white" if vector_format == "ps" else "none",
            )
            figure.savefig(
                f"{output_prefix}.png", format="png", dpi=dpi, facecolor="white"
            )

        save_with_font_fallback(figure, save_outputs)
        return True
    finally:
        plt.close(figure)


def check_st_en_equality(df):
    """Complete a self-comparison across its diagonal without duplicate tiles."""

    if df.empty:
        return df.copy()

    coordinate_columns = ["q_st", "q_en", "r_st", "r_en"]
    unequal_rows = df[(df["q_st"] != df["r_st"]) | (df["q_en"] != df["r_en"])].copy()
    if unequal_rows.empty:
        return df.copy()

    mirrored_rows = unequal_rows.copy()
    mirrored_rows.loc[:, coordinate_columns] = unequal_rows[
        ["r_st", "r_en", "q_st", "q_en"]
    ].to_numpy()

    existing_coordinates = pd.MultiIndex.from_frame(df[coordinate_columns])
    mirrored_coordinates = pd.MultiIndex.from_frame(mirrored_rows[coordinate_columns])
    mirrored_rows = mirrored_rows.loc[~mirrored_coordinates.isin(existing_coordinates)]
    return pd.concat([df, mirrored_rows], ignore_index=True)


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


def _dotplot_tiles(mapping, deraster=False, **kwargs):
    """Create tiles without materializing a genomic-coordinate-sized image.

    ``plotnine.geom_raster`` expands sparse coordinates into an RGBA array whose
    dimensions are derived from the smallest coordinate spacing.  A 496 Mb
    sequence plotted in 2 kb windows can therefore request roughly
    248,000-by-248,000 pixels even when only a small fraction of those cells
    contain matches.  ``geom_tile`` draws only the cells present in the input
    dataframe.  Setting ``raster=True`` keeps the default compact, rasterized
    appearance in vector output, while ``--deraster`` leaves the tiles as
    vectors.
    """
    return geom_tile(mapping, raster=not deraster, **kwargs)


def get_colors(sdf, ncolors, is_freq, custom_breakpoints):
    if ncolors < 1:
        raise ValueError("At least one color is required")
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
        breaks = np.asarray(custom_breakpoints, dtype=np.float64)
        if len(breaks) != ncolors + 1:
            raise ValueError(
                "The number of breakpoints must equal the number of colors plus one"
            )
        if not np.all(np.isfinite(breaks)):
            raise ValueError("Breakpoints must contain only finite numbers")
        if np.any(np.diff(breaks) <= 0):
            raise ValueError("Breakpoints must be strictly increasing")
        values = np.asarray(sdf["perID_by_events"], dtype=np.float64)
        if values.size and (
            not np.all(np.isfinite(values))
            or values.min() < breaks[0]
            or values.max() > breaks[-1]
        ):
            raise ValueError(
                "Breakpoints must cover all finite identity values in the plot"
            )
    labels = np.arange(len(breaks) - 1)
    # A dataset containing only 100% identity creates repeated default bin
    # edges; frequency bins likewise collapse to one edge when every value is
    # equal.  Both cases represent a single category and must not reach
    # ``pandas.cut``, which requires unique edges.
    if len(np.unique(breaks)) < 2:
        return pd.Categorical(
            np.zeros(len(sdf["perID_by_events"]), dtype=int), categories=[0]
        )
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

    # Generate only breakpoints that fall inside the requested interval.
    # Matplotlib expands an axis when ``set_ticks`` includes an out-of-range
    # value, which previously turned a ~103 Mb grid into a 125 Mb grid.
    upper_bound = int(min_aligned + (threshold + 1) * magnitude)
    breaks = [
        value
        for value in range(min_aligned, upper_bound, int(magnitude))
        if min_number <= value <= max_number
    ]

    return breaks


def _requested_axis_bounds(requested_limit):
    """Return an exact ``(start, end)`` pair when one was supplied."""

    if isinstance(requested_limit, (tuple, list)):
        if len(requested_limit) != 2:
            raise ValueError("Axis bounds must contain exactly two values")
        start, end = map(float, requested_limit)
        if end <= start:
            raise ValueError("Axis end must be greater than axis start")
        return start, end
    return None


def _data_axis_limits(sdf, requested_limit=None):
    """Resolve plot limits, preserving exact region bounds when provided."""

    requested_bounds = _requested_axis_bounds(requested_limit)
    if requested_bounds is not None:
        return requested_bounds

    minimum = float(min(sdf["q_st"].min(), sdf["r_st"].min()))
    maximum = float(max(sdf["q_en"].max(), sdf["r_en"].max()))
    if requested_limit:
        maximum = max(maximum, float(requested_limit))
    return minimum, maximum


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
    display_x = display_sequence_name(name_x)
    display_y = display_sequence_name(name_y)
    if is_pairwise:
        title_name = f"Comparative Plot: {display_x} vs {display_y}"
    else:
        title_name = f"Self-Identity Plot: {display_x}"
    title_length = 2 * width
    if len(title_name) > 50:
        title_length = 1.5 * width
    elif len(title_name) > 80:
        title_length = width
    sdf, direction_colors, direction_column = _direction_ani_style(sdf)
    direction_coloring = direction_colors is not None
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
    fill_column = direction_column if direction_coloring else "discrete"
    fill_colors = direction_colors if direction_coloring else new_hexcodes
    # Determine the exact genomic interval. A two-value limit is supplied by
    # FASTA mode so blank edge windows do not shrink or extend the plot.
    min_val, max_val = _data_axis_limits(sdf, xlim)

    # If user provides breaks, convert to ints
    if not breaks:
        breaks = generate_breaks(int(min_val), int(max_val))
    else:
        breaks = [int(x) for x in breaks]
    # Compute window size (handling exceptions)
    try:
        window = max(sdf["q_en"] - sdf["q_st"])
    except ValueError:  # Empty dataframe case
        return ggplot(aes(x=[], y=[])) + theme_minimal()

    # Region-qualified names remain in BEDPE data and filenames, but plot
    # headings should show only the underlying FASTA identifier.
    sdf = sdf.copy()
    sdf["q"] = sdf["q"].map(display_sequence_name)
    sdf["r"] = sdf["r"].map(display_sequence_name)

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
        axis_text=element_text(
            family=[DEFAULT_FONT_FAMILY],
            size=clamped_font_size(width, 2.0),
        ),
        axis_ticks_major=element_line(
            size=(width), color="black"
        ),  # Increased tick length
        title=element_text(
            family=[DEFAULT_FONT_FAMILY],
            size=max(MIN_TITLE_SIZE, title_length),
            hjust=0.5,
        ),  # Center title
        axis_title_x=element_text(
            size=clamped_font_size(width, 2.8),
            family=[DEFAULT_FONT_FAMILY],
        ),
        strip_background=element_blank(),  # Remove facet strip background
        strip_text=element_text(
            size=clamped_font_size(width, 1.2), family=[DEFAULT_FONT_FAMILY]
        ),  # Customize facet label text size (optional)
    )

    # Construct the plot arguments
    ggplot_args = (
        ggplot(sdf)
        + scale_color_discrete(guide=None)
        + scale_fill_manual(values=fill_colors, guide=None)
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

    p = ggplot_args + _dotplot_tiles(
        aes(x="q_st", y="r_st", fill=fill_column, height=window, width=window),
        deraster,
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
    title_name = display_sequence_name(title_name)
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
        breaks = [int(x) for x in breaks]
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
        axis_text=element_text(
            family=[DEFAULT_FONT_FAMILY], size=clamped_font_size(width, 1.0)
        ),
        axis_ticks_major=element_line(
            size=(width), color="black"
        ),  # Increased tick length
        title=element_text(
            size=clamped_font_size(width, 1.2, MIN_TITLE_SIZE),
            family=[DEFAULT_FONT_FAMILY],
            alpha=0,
        ),
        axis_title_x=element_text(
            size=clamped_font_size(width, 1.2),
            family=[DEFAULT_FONT_FAMILY],
        ),
        strip_background=element_blank(),  # Remove facet strip background
        strip_text=element_text(
            size=clamped_font_size(width, 1.2), family=[DEFAULT_FONT_FAMILY]
        ),  # Customize facet label text size (optional)
    )

    # Construct the plot arguments
    ggplot_args = (
        ggplot(sdf)
        + scale_color_discrete(guide=None)
        + scale_fill_manual(values=new_hexcodes, guide=None)
        + common_theme
        + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
        + coord_fixed(ratio=1)
        + labs(x="", y="", title="")
    )

    p = ggplot_args + _dotplot_tiles(
        aes(x="q_st", y="r_st", fill="discrete", height=window, width=window),
        deraster,
    )

    return p


def direction_dataframe(
    canonical_matrix,
    forward_matrix,
    window_size,
    name_x,
    name_y,
    self_identity,
    x_offset=0,
    y_offset=0,
):
    """Build plotting records that distinguish forward and reverse matches.

    Canonical k-mers match in either orientation, while forward-only k-mers
    match only same-strand sequence.  A canonical hit missing from the
    forward-only matrix therefore represents a reverse-orientation match.
    """
    canonical_matrix = np.asarray(canonical_matrix, dtype=float)
    forward_matrix = np.asarray(forward_matrix, dtype=float)
    if canonical_matrix.shape != forward_matrix.shape:
        raise ValueError("Canonical and forward matrices must have matching shapes")

    records = []
    for x_index, y_index in np.argwhere(canonical_matrix > 0):
        if self_identity and x_index > y_index:
            continue
        query_start = x_index * window_size + x_offset
        reference_start = y_index * window_size + y_offset
        records.append(
            {
                "q": name_x,
                "q_st": query_start,
                "q_en": query_start + window_size - 1,
                "r": name_y,
                "r_st": reference_start,
                "r_en": reference_start + window_size - 1,
                "direction": (
                    "Forward" if forward_matrix[x_index, y_index] > 0 else "Reverse"
                ),
            }
        )

    dataframe = pd.DataFrame.from_records(
        records,
        columns=["q", "q_st", "q_en", "r", "r_st", "r_en", "direction"],
    )
    dataframe["direction"] = pd.Categorical(
        dataframe["direction"], categories=["Forward", "Reverse"], ordered=True
    )
    return dataframe


def create_direction_plot(
    canonical_matrix,
    forward_matrix,
    window_size,
    directory,
    name_x,
    name_y,
    self_identity,
    width,
    dpi,
    vector_format,
    deraster=False,
    xlim=None,
    axes_labels=None,
    x_offset=0,
    y_offset=0,
):
    """Save a blue/pink plot showing match orientation."""
    dataframe = direction_dataframe(
        canonical_matrix,
        forward_matrix,
        window_size,
        name_x,
        name_y,
        self_identity,
        x_offset,
        y_offset,
    )
    if dataframe.empty:
        print(f"No directional matches found for {name_x} and {name_y}. Skipping.\n")
        return None

    dataframe = dataframe.assign(
        q_position=dataframe["q_st"] + window_size / 2,
        r_position=dataframe["r_st"] + window_size / 2,
    )
    requested_bounds = _requested_axis_bounds(xlim)
    if requested_bounds is not None:
        min_val, max_val = requested_bounds
    else:
        min_val = min(dataframe["q_st"].min(), dataframe["r_st"].min())
        max_val = max(
            dataframe["q_en"].max() + 1,
            dataframe["r_en"].max() + 1,
            xlim or 0,
        )
    breaks = (
        [int(value) for value in axes_labels]
        if axes_labels
        else generate_breaks(int(min_val), int(max_val))
    )
    title = (
        f"Direction Plot: {display_sequence_name(name_x)}"
        if self_identity
        else f"Direction Plot: {display_sequence_name(name_x)} vs "
        f"{display_sequence_name(name_y)}"
    )
    plot = (
        ggplot(dataframe)
        + _dotplot_tiles(
            aes(
                x="q_position",
                y="r_position",
                fill="direction",
                height=window_size,
                width=window_size,
            ),
            deraster,
        )
        + scale_fill_manual(
            values={"Forward": "#2166AC", "Reverse": "#D01C8B"},
            name="Direction",
        )
        + scale_x_continuous(
            labels=make_scale, limits=[min_val, max_val], breaks=breaks
        )
        + scale_y_continuous(
            labels=make_scale, limits=[min_val, max_val], breaks=breaks
        )
        + coord_fixed(ratio=1)
        + labs(
            x="Genomic Position",
            y="",
            title=title,
            caption="Blue: forward   Pink: reverse",
        )
        + theme_light()
        + theme(
            legend_position="none",
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            axis_text=element_text(
                family=[DEFAULT_FONT_FAMILY],
                size=clamped_font_size(width, 1.0),
            ),
            title=element_text(
                family=[DEFAULT_FONT_FAMILY],
                size=clamped_font_size(width, 1.4, MIN_TITLE_SIZE),
                hjust=0.5,
            ),
            axis_title_x=element_text(
                size=clamped_font_size(width, 1.2),
                family=[DEFAULT_FONT_FAMILY],
            ),
        )
    )

    os.makedirs(directory, exist_ok=True)
    filename = (
        f"{name_x}_DIRECTION" if self_identity else f"{name_x}_{name_y}_DIRECTION"
    )
    prefix = os.path.join(directory, filename)
    _save_plot(
        plot,
        width=width,
        height=width,
        dpi=dpi,
        format=vector_format,
        filename=f"{prefix}.{vector_format}",
        verbose=False,
    )
    _save_plot(
        plot,
        width=width,
        height=width,
        dpi=dpi,
        format="png",
        filename=f"{prefix}.png",
        verbose=False,
    )
    print(f"Direction plots saved to {prefix}.png and {prefix}.{vector_format}.\n")
    return plot


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
    if sdf.empty:
        max_val = xlim or 1
        if not breaks:
            breaks = generate_breaks(0, int(max_val))
        else:
            breaks = [int(x) for x in breaks]
        return (
            ggplot(sdf)
            + geom_blank()
            + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + coord_fixed(ratio=1)
            + labs(x=None, y=None, title=None)
            + theme(
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                plot_background=element_blank(),
                panel_background=element_blank(),
                axis_line=element_line(color="black"),
                axis_text=element_text(
                    family=[DEFAULT_FONT_FAMILY],
                    size=clamped_font_size(width, 1.0),
                ),
                axis_ticks_major=element_line(),
                axis_title_x=element_blank(),
                axis_title_y=element_blank(),
            )
        )

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
        breaks = [int(x) for x in breaks]
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
            + _dotplot_tiles(
                aes(x=x_col, y=y_col, fill="discrete", height=window, width=window),
                deraster,
            )
            + scale_color_discrete(guide=None)
            + scale_fill_manual(values=new_hexcodes, guide=None)
            + theme(
                legend_position="none",
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                plot_background=element_blank(),
                panel_background=element_blank(),
                axis_line=element_line(color="black"),
                axis_text=element_text(
                    family=[DEFAULT_FONT_FAMILY],
                    size=clamped_font_size(width, 1.0),
                ),
                axis_ticks_major=element_line(),
                title=element_text(family=[DEFAULT_FONT_FAMILY]),
            )
            + scale_x_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + scale_y_continuous(labels=make_scale, limits=[0, max_val], breaks=breaks)
            + coord_fixed(ratio=1)
            + labs(x=None, y=None, title=None)
        )
    else:
        p = (
            ggplot(sdf)
            + _dotplot_tiles(
                aes(x=x_col, y=y_col, fill="discrete", height=window, width=window),
                deraster,
            )
            + scale_color_discrete(guide=None)
            + scale_fill_manual(values=new_hexcodes, guide=None)
            + theme(
                legend_position="none",
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                plot_background=element_blank(),
                panel_background=element_blank(),
                axis_line=element_line(color="black"),
                axis_text=element_text(
                    family=[DEFAULT_FONT_FAMILY],
                    size=clamped_font_size(width, 1.0),
                ),
                axis_ticks_major=element_line(),
                title=element_text(family=[DEFAULT_FONT_FAMILY]),
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
    title_name = display_sequence_name(title_name)
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
        breaks = [int(x) for x in breaks]
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
            + _dotplot_tiles(
                aes(x="q_st", y="r_st", fill="discrete", height=window, width=window),
                deraster,
                alpha=1.0,
            )  # Ensure full opacity
            + scale_fill_manual(values=new_hexcodes, guide=None)
            + scale_color_discrete(guide=None)
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
                axis_text=element_text(
                    family=[DEFAULT_FONT_FAMILY],
                    size=clamped_font_size(width, 1.0),
                ),
                axis_line_x=element_line(),
                axis_line_y=element_blank(),
                axis_ticks_major_x=element_line(),
                axis_ticks_major_y=element_blank(),
                axis_ticks_major=element_line(size=(width)),
                title=element_text(
                    family=[DEFAULT_FONT_FAMILY],
                    size=clamped_font_size(width, 1.4, MIN_TITLE_SIZE),
                    hjust=0.5,
                ),
                axis_title_x=element_text(
                    size=clamped_font_size(width, 1.4),
                    family=[DEFAULT_FONT_FAMILY],
                ),
                axis_text_y=element_blank(),
            )
        )
        axis = (
            ggplot(sdf)
            + geom_tile(
                aes(x="q_st", y="r_st", fill="discrete", height=window, width=window),
                alpha=0,
            )
            + scale_color_discrete(guide=None)
            + scale_fill_manual(values=new_hexcodes, guide=None)
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
                axis_text=element_text(
                    family=[DEFAULT_FONT_FAMILY],
                    size=clamped_font_size(width, 1.0),
                ),
                axis_ticks_major=element_line(),
                axis_line_x=element_line(),
                axis_line_y=element_blank(),
                axis_ticks_major_x=element_line(),
                axis_ticks_major_y=element_blank(),
                axis_text_x=element_line(),
                axis_text_y=element_blank(),
                plot_title=element_blank(),
                axis_title_x=element_text(
                    size=clamped_font_size(width, 1.2),
                    family=[DEFAULT_FONT_FAMILY],
                ),
            )
        )
    else:
        tri = (
            ggplot(sdf)
            + _dotplot_tiles(
                aes(x="q_st", y="r_st", fill="discrete", height=window, width=window),
                deraster,
                alpha=1.0,
            )  # Ensure full opacity
            + scale_fill_manual(values=new_hexcodes, guide=None)
            + scale_color_discrete(guide=None)
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
                axis_text=element_text(
                    family=[DEFAULT_FONT_FAMILY],
                    size=clamped_font_size(width, 1.0),
                ),
                axis_line_x=element_line(),
                axis_line_y=element_blank(),
                axis_ticks_major_x=element_line(),
                axis_ticks_major_y=element_blank(),
                axis_ticks_major=element_line(),
                axis_text_y=element_blank(),
                title=element_blank(),
                axis_title_x=element_text(
                    size=clamped_font_size(width, 1.2),
                    family=[DEFAULT_FONT_FAMILY],
                ),
            )
        )
        axis = (
            ggplot(sdf)
            + geom_tile(
                aes(x="q_st", y="r_st", fill="discrete", height=window, width=window),
                alpha=0,
            )
            + scale_color_discrete(guide=None)
            + scale_fill_manual(values=new_hexcodes, guide=None)
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
                axis_text=element_text(family=[DEFAULT_FONT_FAMILY]),
                axis_ticks_major=element_line(),
                axis_line_x=element_line(),
                axis_line_y=element_blank(),
                axis_ticks_major_x=element_line(),
                axis_ticks_major_y=element_blank(),
                axis_text_x=element_line(),
                axis_text_y=element_blank(),
                plot_title=element_blank(),
                axis_title_x=element_text(
                    size=clamped_font_size(width, 1.2),
                    family=[DEFAULT_FONT_FAMILY],
                ),
            )
        )

    return tri, axis


def make_tri_axis(sdf, title_name, palette, palette_orientation, colors, breaks, xlim):
    title_name = display_sequence_name(title_name)
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
        + scale_color_discrete(guide=None)
        + scale_fill_manual(
            values=new_hexcodes,
            guide=None,
        )
        + theme(
            legend_position="none",
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            plot_background=element_blank(),
            panel_background=element_blank(),
            axis_line=element_line(color="black"),  # Adjust axis line size
            axis_text=element_text(
                family=[DEFAULT_FONT_FAMILY]
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

    sdf, direction_colors, direction_column = _direction_ani_style(sdf)
    fill_column = direction_column or "discrete"
    fill_colors = direction_colors or new_hexcodes
    p = (
        ggplot(data=sdf, mapping=aes(x="perID_by_events", fill=fill_column))
        + geom_histogram(bins=300)
        + scale_color_cmap(cmap_name="plasma")
        + scale_fill_manual(fill_colors)
        + theme_light()
        + _plot_font_theme()
        + theme(legend_position="none")
        + coord_cartesian(xlim=(bot, 100))
        + xlab("% Identity Estimate")
        + ylab("# of Estimates{}".format(extra))
    )
    return p


def _missing_symmetric_rows(dataframe):
    """Return rows whose coordinate-transposed counterpart is absent.

    FASTA self matrices normally contain only one triangle. Drawing those
    rows a second time with transposed coordinates completes the full plot
    without allocating a second, mirrored dataframe. Loaded BEDPE files may
    already contain both triangles, so the general path checks coordinate
    membership before selecting rows to mirror.
    """

    if dataframe.empty:
        return dataframe.iloc[0:0]

    q_start = dataframe["q_st"]
    q_end = dataframe["q_en"]
    r_start = dataframe["r_st"]
    r_end = dataframe["r_en"]
    off_diagonal = (q_start != r_start) | (q_end != r_end)
    if not off_diagonal.any():
        return dataframe.iloc[0:0]

    # The mirror collection needs coordinates and its resolved color only; do
    # not copy names, identity estimates, or plotting-helper columns from a
    # potentially large dataframe.
    render_columns = ["q_st", "q_en", "r_st", "r_en"]
    render_columns.extend(
        column
        for column in ("discrete", "direction", DIRECTION_ANI_COLUMN)
        if column in dataframe.columns
    )

    # The production FASTA path is strictly triangular. Avoid building two
    # MultiIndexes for that common case; a strict start-coordinate ordering
    # proves that no transposed off-diagonal row can already be present.
    start_difference = (
        q_start.loc[off_diagonal].to_numpy() - r_start.loc[off_diagonal].to_numpy()
    )
    if np.all(start_difference < 0) or np.all(start_difference > 0):
        return dataframe.loc[off_diagonal, render_columns]

    existing = pd.MultiIndex.from_arrays([q_start, q_end, r_start, r_end])
    mirrored = pd.MultiIndex.from_arrays([r_start, r_end, q_start, q_end])
    return dataframe.loc[off_diagonal & ~mirrored.isin(existing), render_columns]


def _full_plot_limits(dataframe, requested_limit):
    """Resolve exact full-plot bounds, including an empty sparse matrix."""

    requested_bounds = _requested_axis_bounds(requested_limit)
    if requested_bounds is not None:
        return requested_bounds
    if not dataframe.empty:
        return _data_axis_limits(dataframe, requested_limit)
    if requested_limit:
        return 0.0, float(requested_limit)
    raise ValueError(
        "Cannot infer full-plot bounds from an empty identity table; "
        "provide explicit axis bounds"
    )


def _build_full_figure(
    sdf,
    name_x,
    name_y,
    palette,
    palette_orientation,
    custom_colors,
    axes_labels,
    xlim,
    deraster,
    width,
    is_pairwise,
):
    """Build a sparse full or comparative dotplot with native Matplotlib.

    The tile geometry intentionally retains the historical Plotnine contract:
    BEDPE start coordinates are tile centers and ``q_en - q_st`` is the tile
    width. The default rasterizes only the sparse tile collections in vector
    output; ``--deraster`` leaves each tile as vector geometry.
    """

    display_x = display_sequence_name(name_x)
    display_y = display_sequence_name(name_y)
    title = (
        f"Comparative Plot: {display_x} vs {display_y}"
        if is_pairwise
        else f"Self-Identity Plot: {display_x}"
    )
    region_start, region_end = _full_plot_limits(sdf, xlim)
    breaks = (
        [float(value) for value in axes_labels]
        if axes_labels
        else generate_breaks(int(region_start), int(region_end))
    )

    styled, direction_colors, direction_column = _direction_ani_style(sdf)
    colors = direction_colors or _resolve_native_colors(
        palette, palette_orientation, custom_colors
    )
    color_column = direction_column or "discrete"

    figure, axis = plt.subplots(figsize=(float(width), float(width)))
    try:
        draw_rectangular_tiles(
            axis,
            styled,
            colors,
            color_column=color_column,
            rasterized=not deraster,
        )
        if not is_pairwise:
            missing_rows = _missing_symmetric_rows(styled)
            if not missing_rows.empty:
                draw_rectangular_tiles(
                    axis,
                    missing_rows,
                    colors,
                    color_column=color_column,
                    transpose=True,
                    rasterized=not deraster,
                )

        configure_dotplot_axis(
            axis,
            region_start,
            region_end,
            breaks=breaks,
        )
        _divisor, unit = genomic_scale(region_end)
        axis.set_xlabel(
            f"Genomic Position ({unit})",
            fontsize=clamped_font_size(width, 2.8),
            fontfamily=DEFAULT_FONT_FAMILY,
        )
        axis.tick_params(
            axis="both",
            labelsize=clamped_font_size(width, 2.0),
            length=max(3.5, float(width)),
            colors="black",
        )
        axis.grid(False)
        axis.set_facecolor("none")
        for spine in axis.spines.values():
            spine.set_color("black")

        # Plotnine's one-cell facet supplies a query label above the panel and
        # a reference label at its right edge. Retain those identifiers while
        # placing the descriptive title independently above them.
        axis.set_title(
            display_x,
            fontsize=clamped_font_size(width, 1.2),
            fontfamily=DEFAULT_FONT_FAMILY,
            pad=5,
        )
        axis.set_ylabel(
            display_y,
            fontsize=clamped_font_size(width, 1.2),
            fontfamily=DEFAULT_FONT_FAMILY,
            rotation=-90,
            labelpad=16,
        )
        axis.yaxis.set_label_position("right")

        title_size = 2.0 * float(width)
        if len(title) > 80:
            title_size = float(width)
        elif len(title) > 50:
            title_size = 1.5 * float(width)
        figure.suptitle(
            title,
            fontsize=max(MIN_TITLE_SIZE, title_size),
            fontfamily=DEFAULT_FONT_FAMILY,
            y=0.975,
        )
        figure.subplots_adjust(
            left=0.14,
            right=0.87,
            bottom=0.14,
            top=0.84,
        )
        set_figure_font_family(figure, DEFAULT_FONT_FAMILY)
    except Exception:
        plt.close(figure)
        raise
    return figure


def _triangle_limits(sdf, xlim):
    requested_bounds = _requested_axis_bounds(xlim)
    if requested_bounds is not None:
        region_start, region_end = requested_bounds
    elif not sdf.empty:
        region_start = max(float(sdf["q_st"].min()), float(sdf["r_st"].min()))
        region_end = max(
            float(sdf["q_en"].max()),
            float(sdf["r_en"].max()),
            float(xlim or 0),
        )
    else:
        raise ValueError(
            "Cannot infer triangle bounds from an empty identity table; "
            "provide explicit axis bounds"
        )
    if region_end <= region_start:
        raise ValueError("Triangle plot end must be greater than its start")
    return region_start, region_end


def _build_triangle_figure(
    sdf,
    title,
    palette,
    palette_orientation,
    custom_colors,
    axes_labels,
    xlim,
    deraster,
    width,
    annotation_df=None,
    annotation_chrom=None,
):
    """Build a triangle, optionally aligned with a BED annotation track."""
    region_start, region_end = _triangle_limits(sdf, xlim)
    breaks = (
        [float(value) for value in axes_labels]
        if axes_labels
        else generate_breaks(int(region_start), int(region_end))
    )
    sdf, direction_colors, direction_column = _direction_ani_style(sdf)
    colors = direction_colors or _resolve_native_colors(
        palette, palette_orientation, custom_colors
    )
    color_column = direction_column or "discrete"
    with_annotation = annotation_df is not None
    layout = create_triangle_layout(width, with_annotation=with_annotation)
    try:
        draw_triangle_tiles(
            layout.triangle_axis,
            sdf,
            colors,
            color_column=color_column,
            rasterized=not deraster,
        )
        configure_triangle_axis(
            layout.triangle_axis,
            region_start,
            region_end,
            breaks=breaks,
            label=not with_annotation,
        )
        if with_annotation:
            # The annotation axis is the sole genomic axis in the combined
            # figure. Remove the triangle baseline and its tick marks.
            layout.triangle_axis.spines["bottom"].set_visible(False)
            layout.triangle_axis.tick_params(axis="x", bottom=False, labelbottom=False)
        layout.triangle_axis.set_title(
            display_sequence_name(title),
            fontsize=clamped_font_size(width, 1.4, MIN_TITLE_SIZE),
            fontfamily=DEFAULT_FONT_FAMILY,
        )

        if with_annotation:
            if not annotation_chrom:
                raise ValueError("An annotation chromosome is required")
            draw_annotation_track(
                layout.annotation_axis,
                annotation_df,
                annotation_chrom,
                region_start,
                region_end,
            )
            layout.annotation_axis.set_xticks(breaks)
            layout.annotation_axis.xaxis.set_major_formatter(
                genomic_tick_formatter(region_end)
            )
            _, unit = genomic_scale(region_end)
            layout.annotation_axis.set_xlabel(f"Genomic Position ({unit})")

        layout.figure.subplots_adjust(
            left=0.08,
            right=0.98,
            bottom=0.14 if not with_annotation else 0.10,
            top=0.90,
            hspace=0.05,
        )
        if with_annotation:
            # ``set_aspect('equal', adjustable='box')`` narrows the triangle
            # axis inside its GridSpec cell to retain 45-degree diagonals. A
            # normal annotation axis keeps the full cell width, so shared data
            # limits alone do not produce physical pixel alignment. Match the
            # BED axis to the triangle's final horizontal bounds.
            triangle_position = layout.triangle_axis.get_position()
            annotation_position = layout.annotation_axis.get_position()
            layout.annotation_axis.set_position(
                [
                    triangle_position.x0,
                    annotation_position.y0,
                    triangle_position.width,
                    annotation_position.height,
                ]
            )
        set_figure_font_family(layout.figure, DEFAULT_FONT_FAMILY)
    except Exception:
        plt.close(layout.figure)
        raise
    return layout.figure


def _ordered_grid_names(single_names, double_names):
    names = []
    for name in single_names:
        if name in names:
            raise ValueError(f"Duplicate self comparison for sequence {name!r}")
        names.append(name)
    for pair in double_names:
        if len(pair) != 2:
            raise ValueError("Each grid comparison must name exactly two sequences")
        for name in pair:
            if name not in names:
                names.append(name)
    if not names:
        raise ValueError("No sequence names were provided for the grid")
    return names


def _read_grid_dataframe(
    matrix,
    is_bed,
    palette,
    palette_orientation,
    is_freq,
    custom_colors,
    custom_breakpoints,
):
    return read_df(
        None if is_bed else [matrix],
        palette,
        palette_orientation,
        is_freq,
        custom_colors,
        custom_breakpoints,
        matrix if is_bed else None,
    )


def _grid_axis_limits(dataframes, requested_limit):
    requested_bounds = _requested_axis_bounds(requested_limit)
    if requested_bounds is not None:
        return requested_bounds
    if requested_limit:
        return 0.0, float(requested_limit)
    minima = [
        float(min(dataframe["q_st"].min(), dataframe["r_st"].min()))
        for dataframe in dataframes
        if not dataframe.empty
    ]
    maxima = [
        float(max(dataframe["q_en"].max(), dataframe["r_en"].max()))
        for dataframe in dataframes
        if not dataframe.empty
    ]
    if not maxima:
        return 0.0, 1.0
    return min(minima), max(maxima)


def _build_grid_figure(
    singles,
    doubles,
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
):
    """Build a native Matplotlib comparison grid and return its axes.

    ``width`` controls the complete square figure, not each panel, so memory
    use does not grow quadratically in pixels as sequences are added.
    """
    if len(singles) != len(single_names):
        raise ValueError("Self-comparison matrices and names must have equal lengths")
    if len(doubles) != len(double_names):
        raise ValueError("Pairwise matrices and names must have equal lengths")

    names = _ordered_grid_names(single_names, double_names)
    # Keep columns in input order and reverse rows so self-comparisons occupy
    # the anti-diagonal: bottom-left to top-right.
    row_names = list(reversed(names))
    single_frames = {}
    for name, matrix in zip(single_names, singles):
        single_frames[name] = _read_grid_dataframe(
            matrix,
            is_bed,
            palette,
            palette_orientation,
            is_freq,
            custom_colors,
            custom_breakpoints,
        )

    pair_frames = {}
    pair_orientations = {}
    for pair, matrix in zip(double_names, doubles):
        query_name, reference_name = pair
        if query_name == reference_name:
            raise ValueError("Pairwise grid comparisons must use two distinct names")
        key = frozenset((query_name, reference_name))
        if key in pair_frames:
            raise ValueError(
                f"Duplicate pairwise comparison for {query_name!r} and {reference_name!r}"
            )
        pair_frames[key] = _read_grid_dataframe(
            matrix,
            is_bed,
            palette,
            palette_orientation,
            is_freq,
            custom_colors,
            custom_breakpoints,
        )
        pair_orientations[key] = (query_name, reference_name)

    missing_pairs = [
        (names[row], names[column])
        for row in range(len(names))
        for column in range(row + 1, len(names))
        if frozenset((names[row], names[column])) not in pair_frames
    ]
    if missing_pairs:
        formatted = ", ".join(f"{left}/{right}" for left, right in missing_pairs)
        raise ValueError(f"Missing pairwise grid comparisons: {formatted}")

    all_frames = list(single_frames.values()) + list(pair_frames.values())
    axis_start, axis_end = _grid_axis_limits(all_frames, xlim)
    _, axis_unit = genomic_scale(axis_end)
    axis_breaks = axes_label or breaks
    if not axis_breaks:
        axis_breaks = generate_breaks(int(axis_start), int(axis_end))
    axis_breaks = [float(value) for value in axis_breaks]
    colors = _resolve_native_colors(palette, palette_orientation, custom_colors)

    grid_size = len(names)
    figure_width = max(float(width), 2.0)
    heading_size = clamped_font_size(figure_width, 1.2, 8.0, 12.0)
    # Numeric genomic labels are intentionally twice the previous size. The
    # inverse grid-size factor keeps larger grids proportionate.
    tick_size = clamped_font_size(
        figure_width,
        3.0 / grid_size,
        MIN_TEXT_SIZE,
        18.0,
    )
    axis_title_size = clamped_font_size(
        figure_width,
        1.6 / grid_size,
        MIN_TEXT_SIZE,
        14.0,
    )
    figure, axes = plt.subplots(
        grid_size,
        grid_size,
        figsize=(figure_width, figure_width),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    try:
        for row, row_name in enumerate(row_names):
            for column, column_name in enumerate(names):
                axis = axes[row, column]
                dataframe = None
                transpose = False
                if row_name == column_name:
                    dataframe = single_frames.get(row_name)
                else:
                    key = frozenset((row_name, column_name))
                    dataframe = pair_frames[key]
                    query_name, reference_name = pair_orientations[key]
                    if (query_name, reference_name) == (column_name, row_name):
                        transpose = False
                    elif (query_name, reference_name) == (row_name, column_name):
                        transpose = True
                    else:
                        raise ValueError(
                            f"Grid comparison names do not match {row_name!r}/{column_name!r}"
                        )

                if dataframe is not None and not dataframe.empty:
                    if row_name == column_name:
                        dataframe = check_st_en_equality(dataframe)
                    (
                        dataframe,
                        direction_colors,
                        direction_column,
                    ) = _direction_ani_style(dataframe)
                    draw_rectangular_tiles(
                        axis,
                        dataframe,
                        direction_colors or colors,
                        color_column=direction_column or "discrete",
                        transpose=transpose,
                        rasterized=not deraster,
                    )
                configure_dotplot_axis(
                    axis,
                    axis_start,
                    axis_end,
                    breaks=axis_breaks,
                    show_x=row == grid_size - 1,
                    show_y=column == 0,
                )
                axis.tick_params(axis="both", labelsize=tick_size)
                axis.grid(False)
                if row == 0:
                    axis.set_title(
                        display_sequence_name(column_name),
                        fontsize=heading_size,
                        fontfamily=DEFAULT_FONT_FAMILY,
                    )
                if column == 0:
                    axis.set_ylabel(
                        display_sequence_name(row_name),
                        fontsize=heading_size,
                        fontfamily=DEFAULT_FONT_FAMILY,
                        labelpad=2,
                    )

        # Keep both shared genomic-axis titles with the bottom-left cell, where
        # both sets of numeric tick labels are visible. Figure-wide titles sit
        # far from that cell and enlarge tightly cropped output canvases.
        bottom_left_axis = axes[-1, 0]
        axis_title = f"Genomic Position ({axis_unit})"
        bottom_left_axis.set_xlabel(
            axis_title,
            fontsize=axis_title_size,
            fontfamily=DEFAULT_FONT_FAMILY,
            labelpad=2,
        )
        bottom_margin = max(0.12, 0.35 / figure_width)
        left_margin_inches = 0.72 + max(0.0, tick_size - MIN_TEXT_SIZE) / 72.0
        left_margin = max(0.14, left_margin_inches / figure_width)
        figure.subplots_adjust(
            left=left_margin,
            right=0.98,
            bottom=bottom_margin,
            top=0.92,
            wspace=0.08,
            hspace=0.08,
        )
        _fit_grid_sequence_labels(figure, axes)
        vertical_title = bottom_left_axis.annotate(
            axis_title,
            xy=(0, 0.5),
            xycoords=bottom_left_axis.yaxis.label,
            xytext=(-4, 0),
            textcoords="offset points",
            ha="center",
            va="center",
            rotation=90,
            rotation_mode="anchor",
            fontsize=axis_title_size,
            fontfamily=DEFAULT_FONT_FAMILY,
            annotation_clip=False,
        )
        vertical_title.set_gid("grid-y-axis-title")
        set_figure_font_family(figure, DEFAULT_FONT_FAMILY)
    except Exception:
        plt.close(figure)
        raise
    return figure, axes


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
    dpi=300,
):
    figure, axes = _build_grid_figure(
        singles=singles,
        doubles=doubles,
        palette=palette,
        palette_orientation=palette_orientation,
        single_names=single_names,
        double_names=double_names,
        is_freq=is_freq,
        xlim=xlim,
        custom_colors=custom_colors,
        custom_breakpoints=custom_breakpoints,
        axes_label=axes_label,
        is_bed=is_bed,
        width=width,
        breaks=breaks,
        deraster=deraster,
    )
    grid_size = axes.shape[0]
    directional = any(
        (
            "direction" in matrix.columns
            if isinstance(matrix, pd.DataFrame)
            else bool(matrix) and "direction" in matrix[0]
        )
        for matrix in [*singles, *doubles]
    )
    grid_label = "DIRECTION_GRID" if directional else "GRID"
    grid_prefix = os.path.join(directory, f"{grid_size}x{grid_size}_{grid_label}")
    print(f"\nGrid complete! Saving to {grid_prefix}...\n")
    try:
        save_figure_pair(
            figure,
            grid_prefix,
            vector_format,
            dpi,
            bbox_inches=figure.bbox_inches,
        )
    finally:
        plt.close(figure)
    print("Grid saved successfully!\n")
    return [f"{grid_prefix}.{vector_format}", f"{grid_prefix}.png"]


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
    os.makedirs(directory, exist_ok=True)
    created_files = []
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
    directional = "direction" in sdf.columns

    plot_filename = os.path.join(directory, name_x)

    if is_pairwise:
        plot_filename = os.path.join(directory, f"{name_x}_{name_y}")

    histy = make_hist(
        sdf, palette, palette_orientation, custom_colors, custom_breakpoints
    )

    annotation_track_created = False
    annotation_bed_df = None
    annotation_chrom = None
    # Just doing triangle plots for now.
    if annotation:
        print("Generating annotation track:\n")
        iniprefix = plot_filename
        requested_bounds = _requested_axis_bounds(xlim)
        if requested_bounds is not None:
            min_val, annotation_end = requested_bounds
        else:
            min_val = max(sdf["q_st"].min(), sdf["r_st"].min())
            max_val = max(sdf["q_en"].max(), sdf["r_en"].max())
            annotation_end = xlim or max_val
        chrom_name = name_x.split(":")[0]
        try:
            bed_df = read_annotation_bed(annotation)
            annotation_track_created = render_annotation_track(
                bed_df=bed_df,
                chrom=chrom_name,
                region_start=min_val,
                region_end=annotation_end,
                output_prefix=f"{iniprefix}_ANNOTATION_TRACK",
                width_cm=width * 2.05,
                dpi=dpi,
                vector_format=vector_format,
            )
            if annotation_track_created:
                created_files.extend(
                    [
                        f"{iniprefix}_ANNOTATION_TRACK.{vector_format}",
                        f"{iniprefix}_ANNOTATION_TRACK.png",
                    ]
                )
                annotation_bed_df = bed_df
                annotation_chrom = chrom_name
                print(f"\nAnnotation track saved to {iniprefix}_ANNOTATION_TRACK\n")
            else:
                print(
                    f"No valid intervals found in {annotation} for region "
                    f"{chrom_name}:{min_val}-{annotation_end}.\n"
                )
                print("Skipping annotation track generation.\n")
        except Exception as e:
            print(f"Error processing annotation file {annotation}: {e}\n")
            print("Skipping annotation track generation.\n")

    if is_pairwise:
        print(f"Creating plots and saving to {plot_filename}...\n")
        full_suffix = "_DIRECTION_FULL" if directional else "_COMPARE"
        hist_suffix = "_DIRECTION_HIST" if directional else "_COMPARE_HIST"
        full_figure = _build_full_figure(
            sdf=sdf,
            name_x=name_x,
            name_y=name_y,
            palette=palette,
            palette_orientation=palette_orientation,
            custom_colors=custom_colors,
            axes_labels=axes_labels,
            xlim=xlim,
            deraster=deraster,
            width=width,
            is_pairwise=True,
        )
        try:
            save_figure_pair(
                full_figure,
                f"{plot_filename}{full_suffix}",
                vector_format,
                dpi,
                bbox_inches=full_figure.bbox_inches,
            )
        finally:
            plt.close(full_figure)
        created_files.extend(
            [
                f"{plot_filename}{full_suffix}.{vector_format}",
                f"{plot_filename}{full_suffix}.png",
            ]
        )
        if not no_hist:
            _draw_and_save_plot_pair(
                histy,
                f"{plot_filename}{hist_suffix}",
                width=3,
                height=3,
                dpi=dpi,
                vector_format=vector_format,
            )
            created_files.extend(
                [
                    f"{plot_filename}{hist_suffix}.{vector_format}",
                    f"{plot_filename}{hist_suffix}.png",
                ]
            )
        if not no_hist:
            print(
                f"{plot_filename} comparative plots and histogram saved sucessfully. \n"
            )
        else:
            print(
                f"{plot_filename}{full_suffix}.{vector_format} and "
                f"{plot_filename}{full_suffix}.png saved sucessfully. \n"
            )
    # Self-identity plots: Output _TRI, _FULL, and _HIST
    else:
        if deraster:
            print(
                f"Producing dotplots with derasterization turned off. This may take a while...\n"
            )
        full_suffix = "_DIRECTION_FULL" if directional else "_FULL"
        tri_suffix = "_DIRECTION_TRI" if directional else "_TRI"
        hist_suffix = "_DIRECTION_HIST" if directional else "_HIST"
        full_figure = _build_full_figure(
            sdf=sdf,
            name_x=name_x,
            name_y=name_y,
            palette=palette,
            palette_orientation=palette_orientation,
            custom_colors=custom_colors,
            axes_labels=axes_labels,
            xlim=xlim,
            deraster=deraster,
            width=width,
            is_pairwise=False,
        )
        try:
            save_figure_pair(
                full_figure,
                f"{plot_filename}{full_suffix}",
                vector_format,
                dpi,
                bbox_inches=full_figure.bbox_inches,
            )
        finally:
            plt.close(full_figure)
        created_files.extend(
            [
                f"{plot_filename}{full_suffix}.{vector_format}",
                f"{plot_filename}{full_suffix}.png",
            ]
        )
        tri_prefix = f"{plot_filename}{tri_suffix}"
        triangle_figure = _build_triangle_figure(
            sdf=sdf,
            title=name_x,
            palette=palette,
            palette_orientation=palette_orientation,
            custom_colors=custom_colors,
            axes_labels=axes_labels,
            xlim=xlim,
            deraster=deraster,
            width=width,
        )
        try:
            save_figure_pair(
                triangle_figure,
                tri_prefix,
                vector_format,
                dpi,
                bbox_inches="tight",
            )
        finally:
            plt.close(triangle_figure)
        created_files.extend([f"{tri_prefix}.{vector_format}", f"{tri_prefix}.png"])

        if annotation_track_created:
            annotated_figure = _build_triangle_figure(
                sdf=sdf,
                title=name_x,
                palette=palette,
                palette_orientation=palette_orientation,
                custom_colors=custom_colors,
                axes_labels=axes_labels,
                xlim=xlim,
                deraster=deraster,
                width=width,
                annotation_df=annotation_bed_df,
                annotation_chrom=annotation_chrom,
            )
            try:
                save_figure_pair(
                    annotated_figure,
                    f"{tri_prefix}_ANNOTATED",
                    vector_format,
                    dpi,
                    bbox_inches="tight",
                )
            finally:
                plt.close(annotated_figure)
            created_files.extend(
                [
                    f"{tri_prefix}_ANNOTATED.{vector_format}",
                    f"{tri_prefix}_ANNOTATED.png",
                ]
            )

        if no_hist:
            print(
                f"Triangle plots and full plots for {plot_filename} saved sucessfully. \n"
            )
        else:
            _draw_and_save_plot_pair(
                histy,
                f"{plot_filename}{hist_suffix}",
                width=3,
                height=3,
                dpi=dpi,
                vector_format=vector_format,
            )
            created_files.extend(
                [
                    plot_filename + f"{hist_suffix}.{vector_format}",
                    plot_filename + f"{hist_suffix}.png",
                ]
            )
            print(
                f"Triangle plots, full plots, and histogram for {plot_filename} saved sucessfully. \n"
            )
    return created_files
