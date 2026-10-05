"""Sparse, native Matplotlib rendering primitives for ModDotPlot.

The functions in this module deliberately know nothing about named palettes or
how percent-identity bins are calculated.  Callers provide a sequence or
mapping of colors after applying their chosen palette.  Keeping those concerns
separate makes the geometry usable by individual dotplots, triangle plots, and
multi-panel grids without creating intermediate SVG files.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.collections import PolyCollection
from matplotlib.figure import Figure
from matplotlib.text import Text
from matplotlib.ticker import FuncFormatter
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd

ColorSource = Union[Sequence[str], Mapping[object, str]]
TickFormatter = Callable[[float, int], str]

DEFAULT_FONT_FAMILY = "Helvetica"
FALLBACK_FONT_FAMILY = "DejaVu Sans"
MIN_TEXT_SIZE = 8.0
MIN_TITLE_SIZE = 10.0


def clamped_font_size(
    width: float,
    multiplier: float,
    minimum: float = MIN_TEXT_SIZE,
    maximum: Optional[float] = None,
) -> float:
    """Scale a font with figure width without allowing unreadable sizes."""

    size = float(width) * float(multiplier)
    if not np.isfinite(size):
        raise ValueError("Calculated font size must be finite")
    size = max(float(minimum), size)
    if maximum is not None:
        size = min(float(maximum), size)
    return size


def set_figure_font_family(figure: Figure, family: str) -> None:
    """Set every existing text artist in a figure to one font family."""

    for artist in figure.findobj(match=Text):
        artist.set_fontfamily(family)


def is_glyph_loading_error(error: BaseException) -> bool:
    """Return whether Matplotlib failed while loading a font glyph."""

    return (
        isinstance(error, RuntimeError) and "failed to load glyph" in str(error).lower()
    )


def save_with_font_fallback(figure: Figure, save: Callable[[], None]) -> None:
    """Save with Helvetica, retrying the complete operation with DejaVu Sans."""

    set_figure_font_family(figure, DEFAULT_FONT_FAMILY)
    try:
        save()
    except RuntimeError as error:
        if not is_glyph_loading_error(error):
            raise
        set_figure_font_family(figure, FALLBACK_FONT_FAMILY)
        save()


@dataclass(frozen=True)
class TriangleLayout:
    """Axes belonging to a native triangle figure.

    ``annotation_axis`` is ``None`` for an unannotated layout.  When present it
    shares its x-axis with ``triangle_axis``, so BED features and transformed
    triangle tiles use the same genomic coordinates.
    """

    figure: Figure
    triangle_axis: Axes
    annotation_axis: Optional[Axes]


def tile_width(dataframe: pd.DataFrame) -> float:
    """Return the legacy tile width, ``max(q_en - q_st)``.

    ModDotPlot historically passes ``q_st`` and ``r_st`` to ``geom_tile`` as
    tile centers.  The native renderer intentionally preserves that behavior
    for the first release of the new plotting pipeline.
    """

    _require_columns(dataframe, ("q_st", "q_en"))
    if dataframe.empty:
        return 0.0
    widths = _numeric_values(dataframe, "q_en") - _numeric_values(dataframe, "q_st")
    if not np.all(np.isfinite(widths)) or np.any(widths <= 0):
        raise ValueError("Tile widths must be finite and greater than zero")
    return float(np.max(widths))


def rectangular_tile_vertices(
    dataframe: pd.DataFrame, *, transpose: bool = False
) -> np.ndarray:
    """Build one square polygon per sparse input row.

    The squares are centered on the start columns and all use the maximum query
    interval width, matching the existing plotnine tile semantics.  Transpose
    swaps query and reference positions without copying or expanding the sparse
    input into a dense genomic matrix.
    """

    _require_columns(dataframe, ("q_st", "q_en", "r_st"))
    if dataframe.empty:
        return np.empty((0, 4, 2), dtype=float)

    window = tile_width(dataframe)
    half_window = window / 2.0
    query = _numeric_values(dataframe, "q_st")
    reference = _numeric_values(dataframe, "r_st")
    if not np.all(np.isfinite(query)) or not np.all(np.isfinite(reference)):
        raise ValueError("Tile centers must contain only finite numbers")
    if transpose:
        query, reference = reference, query

    vertices = np.empty((query.size, 4, 2), dtype=float)
    vertices[:, 0, 0] = query - half_window
    vertices[:, 0, 1] = reference - half_window
    vertices[:, 1, 0] = query + half_window
    vertices[:, 1, 1] = reference - half_window
    vertices[:, 2, 0] = query + half_window
    vertices[:, 2, 1] = reference + half_window
    vertices[:, 3, 0] = query - half_window
    vertices[:, 3, 1] = reference + half_window
    return vertices


def transform_triangle_points(points: np.ndarray) -> np.ndarray:
    """Rotate dotplot coordinates into triangle coordinates.

    For each ``(query, reference)`` point, the result is
    ``((query + reference) / 2, (reference - query) / 2)``.
    """

    values = np.asarray(points, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("Triangle points must be an N-by-2 array")
    if not np.all(np.isfinite(values)):
        raise ValueError("Triangle points must contain only finite numbers")
    query = values[:, 0]
    reference = values[:, 1]
    return np.column_stack(((query + reference) / 2.0, (reference - query) / 2.0))


def triangle_tile_vertices(
    dataframe: pd.DataFrame, *, baseline: float = 0.0
) -> Sequence[np.ndarray]:
    """Transform sparse square tiles and clip them at the triangle baseline."""

    if not np.isfinite(baseline):
        raise ValueError("Triangle baseline must be finite")
    polygons = []
    for rectangle in rectangular_tile_vertices(dataframe):
        transformed = transform_triangle_points(rectangle)
        clipped = _clip_polygon_above_baseline(transformed, float(baseline))
        if len(clipped) >= 3:
            polygons.append(clipped)
    return polygons


def draw_rectangular_tiles(
    axis: Axes,
    dataframe: pd.DataFrame,
    colors: ColorSource,
    *,
    color_column: str = "discrete",
    transpose: bool = False,
    rasterized: bool = True,
    edgecolors: str = "none",
    linewidth: float = 0.0,
    alpha: float = 1.0,
) -> PolyCollection:
    """Draw sparse rectangular tiles and return their ``PolyCollection``."""

    vertices = rectangular_tile_vertices(dataframe, transpose=transpose)
    facecolors = _resolve_facecolors(dataframe, colors, color_column)
    collection = PolyCollection(
        vertices,
        facecolors=facecolors,
        edgecolors=edgecolors,
        linewidths=linewidth,
        rasterized=rasterized,
        alpha=alpha,
    )
    axis.add_collection(collection)
    return collection


def draw_triangle_tiles(
    axis: Axes,
    dataframe: pd.DataFrame,
    colors: ColorSource,
    *,
    color_column: str = "discrete",
    baseline: float = 0.0,
    rasterized: bool = True,
    edgecolors: str = "none",
    linewidth: float = 0.0,
    alpha: float = 1.0,
) -> PolyCollection:
    """Draw transformed sparse tiles, clipped to ``y >= baseline``."""

    vertices = triangle_tile_vertices(dataframe, baseline=baseline)
    # A tile can be wholly below a non-default baseline, so resolve only the
    # corresponding leading colors.  The normal ModDotPlot baseline is zero
    # and its upper-triangle input retains every row.
    facecolors = _resolve_facecolors(dataframe, colors, color_column)
    if len(vertices) != len(facecolors):
        facecolors = _triangle_visible_facecolors(
            dataframe, colors, color_column, baseline
        )
    collection = PolyCollection(
        vertices,
        facecolors=facecolors,
        edgecolors=edgecolors,
        linewidths=linewidth,
        rasterized=rasterized,
        alpha=alpha,
    )
    axis.add_collection(collection)
    return collection


def genomic_scale(maximum: float) -> Tuple[float, str]:
    """Return the divisor and unit used for genomic tick labels."""

    if not np.isfinite(maximum):
        raise ValueError("Genomic axis maximum must be finite")
    maximum = abs(float(maximum))
    if maximum < 200_000:
        return 1_000.0, "Kbp"
    if maximum > 200_000_000:
        return 1_000_000_000.0, "Gbp"
    return 1_000_000.0, "Mbp"


def genomic_tick_formatter(maximum: float) -> FuncFormatter:
    """Create a Matplotlib formatter matching ModDotPlot's genomic scaling."""

    divisor, _ = genomic_scale(maximum)

    def format_tick(value: float, _position: int) -> str:
        scaled = value / divisor
        return f"{scaled:g}"

    return FuncFormatter(format_tick)


def configure_dotplot_axis(
    axis: Axes,
    region_start: float,
    region_end: float,
    *,
    breaks: Optional[Sequence[float]] = None,
    formatter: Optional[TickFormatter] = None,
    show_x: bool = True,
    show_y: bool = True,
    equal_aspect: bool = True,
) -> Axes:
    """Configure limits and scaled ticks for a rectangular dotplot axis."""

    _validate_region(region_start, region_end)
    axis.set_xlim(region_start, region_end)
    axis.set_ylim(region_start, region_end)
    if breaks is not None:
        axis.set_xticks(breaks)
        axis.set_yticks(breaks)
        # Matplotlib deliberately expands view limits to make every fixed tick
        # visible. Re-apply the genomic interval so custom or automatically
        # generated ticks can never add whitespace outside the data bounds.
        axis.set_xlim(region_start, region_end)
        axis.set_ylim(region_start, region_end)
    tick_formatter = formatter or genomic_tick_formatter(region_end)
    axis.xaxis.set_major_formatter(
        FuncFormatter(tick_formatter) if formatter else tick_formatter
    )
    axis.yaxis.set_major_formatter(
        FuncFormatter(tick_formatter) if formatter else tick_formatter
    )
    axis.tick_params(axis="x", labelbottom=show_x, bottom=show_x)
    axis.tick_params(axis="y", labelleft=show_y, left=show_y)
    if equal_aspect:
        axis.set_aspect("equal", adjustable="box")
    return axis


def configure_triangle_axis(
    axis: Axes,
    region_start: float,
    region_end: float,
    *,
    breaks: Optional[Sequence[float]] = None,
    formatter: Optional[TickFormatter] = None,
    baseline: float = 0.0,
    label: bool = True,
) -> Axes:
    """Configure a transformed triangle axis with a genomic x-axis."""

    _validate_region(region_start, region_end)
    if not np.isfinite(baseline):
        raise ValueError("Triangle baseline must be finite")
    axis.set_xlim(region_start, region_end)
    axis.set_ylim(baseline, baseline + (region_end - region_start) / 2.0)
    if breaks is not None:
        axis.set_xticks(breaks)
        axis.set_xlim(region_start, region_end)
    tick_formatter = formatter or genomic_tick_formatter(region_end)
    axis.xaxis.set_major_formatter(
        FuncFormatter(tick_formatter) if formatter else tick_formatter
    )
    axis.set_yticks([])
    axis.set_aspect("equal", adjustable="box")
    axis.spines["left"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["top"].set_visible(False)
    if label:
        _, unit = genomic_scale(region_end)
        axis.set_xlabel(f"Genomic Position ({unit})")
    return axis


def create_triangle_layout(
    width: float,
    *,
    with_annotation: bool = False,
    annotation_height: float = 0.8,
    hspace: float = 0.05,
) -> TriangleLayout:
    """Create triangle axes, optionally with a shared BED annotation axis."""

    width = float(width)
    annotation_height = float(annotation_height)
    if not np.isfinite(width) or width <= 0:
        raise ValueError("Figure width must be finite and greater than zero")
    if not np.isfinite(annotation_height) or annotation_height <= 0:
        raise ValueError("Annotation height must be finite and greater than zero")

    triangle_height = width / 2.0
    extra_height = annotation_height if with_annotation else 0.0
    figure = plt.figure(figsize=(width, triangle_height + extra_height))
    if not with_annotation:
        triangle_axis = figure.add_subplot(1, 1, 1)
        return TriangleLayout(figure, triangle_axis, None)

    grid = figure.add_gridspec(
        2,
        1,
        height_ratios=(triangle_height, annotation_height),
        hspace=hspace,
    )
    triangle_axis = figure.add_subplot(grid[0, 0])
    annotation_axis = figure.add_subplot(grid[1, 0], sharex=triangle_axis)
    triangle_axis.tick_params(axis="x", labelbottom=False)
    return TriangleLayout(figure, triangle_axis, annotation_axis)


def save_figure_pair(
    figure: Figure,
    output_prefix: Union[str, Path],
    vector_format: str,
    dpi: int,
    *,
    transparent: bool = False,
    bbox_inches: Optional[Union[str, Bbox]] = "tight",
) -> Tuple[Path, Path]:
    """Save one figure directly as PNG and SVG, PDF, or PostScript."""

    vector_format = str(vector_format).lower().lstrip(".")
    if vector_format not in {"svg", "pdf", "ps"}:
        raise ValueError("Vector format must be one of: svg, pdf, ps")
    if int(dpi) <= 0:
        raise ValueError("DPI must be greater than zero")

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    # Sequence identifiers frequently contain dots (for example
    # ``PAN010.chr14``), so ``Path.with_suffix`` would discard part of a valid
    # output prefix. Append extensions exactly as the legacy CLI does.
    png_path = Path(f"{prefix}.png")
    vector_path = Path(f"{prefix}.{vector_format}")
    save_options = {
        "dpi": int(dpi),
        "transparent": transparent,
        "bbox_inches": bbox_inches,
    }

    def save_outputs() -> None:
        figure.savefig(png_path, format="png", **save_options)
        figure.savefig(vector_path, format=vector_format, **save_options)

    save_with_font_fallback(figure, save_outputs)
    return png_path, vector_path


def _require_columns(dataframe: pd.DataFrame, columns: Sequence[str]) -> None:
    missing = [column for column in columns if column not in dataframe.columns]
    if missing:
        raise ValueError(f"Missing required tile columns: {', '.join(missing)}")


def _numeric_values(dataframe: pd.DataFrame, column: str) -> np.ndarray:
    try:
        return pd.to_numeric(dataframe[column], errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Tile column '{column}' must contain only numbers") from error


def _resolve_facecolors(
    dataframe: pd.DataFrame, colors: ColorSource, color_column: str
) -> Sequence[str]:
    if dataframe.empty:
        return []
    _require_columns(dataframe, (color_column,))
    values = dataframe[color_column]

    if isinstance(colors, Mapping):
        try:
            return [colors[value] for value in values]
        except KeyError as error:
            raise ValueError(
                f"No color was provided for tile category {error.args[0]!r}"
            ) from error

    if isinstance(colors, str):
        palette = [colors]
    else:
        palette = list(colors)
    if not palette:
        raise ValueError("At least one tile color is required")

    if isinstance(values.dtype, pd.CategoricalDtype):
        categories = list(values.cat.categories)
    else:
        categories = list(pd.unique(values))
    if len(palette) < len(categories):
        raise ValueError(
            "The color sequence must contain at least one color per tile category"
        )
    color_by_category = dict(zip(categories, palette))
    return [color_by_category[value] for value in values]


def _triangle_visible_facecolors(
    dataframe: pd.DataFrame,
    colors: ColorSource,
    color_column: str,
    baseline: float,
) -> Sequence[str]:
    all_colors = _resolve_facecolors(dataframe, colors, color_column)
    rectangles = rectangular_tile_vertices(dataframe)
    return [
        color
        for rectangle, color in zip(rectangles, all_colors)
        if len(
            _clip_polygon_above_baseline(
                transform_triangle_points(rectangle), float(baseline)
            )
        )
        >= 3
    ]


def _clip_polygon_above_baseline(polygon: np.ndarray, baseline: float) -> np.ndarray:
    """Clip a convex polygon against the half-plane ``y >= baseline``."""

    if len(polygon) == 0:
        return np.empty((0, 2), dtype=float)
    output = []
    previous = polygon[-1]
    previous_inside = previous[1] >= baseline
    for current in polygon:
        current_inside = current[1] >= baseline
        if current_inside != previous_inside:
            fraction = (baseline - previous[1]) / (current[1] - previous[1])
            output.append(previous + fraction * (current - previous))
        if current_inside:
            output.append(current)
        previous = current
        previous_inside = current_inside
    return np.asarray(output, dtype=float).reshape((-1, 2))


def _validate_region(region_start: float, region_end: float) -> None:
    if not np.isfinite(region_start) or not np.isfinite(region_end):
        raise ValueError("Genomic axis limits must be finite")
    if region_end <= region_start:
        raise ValueError("Genomic region end must be greater than its start")
