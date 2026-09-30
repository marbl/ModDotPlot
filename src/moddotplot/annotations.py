"""Shared BED annotation parsing and interval-selection helpers."""

import numpy as np
import pandas as pd

DEFAULT_ANNOTATION_COLOR = "#4C72B0"
BED_COLUMNS = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
]


def read_annotation_bed(filepath):
    """Read the BED3-BED9 subset used by ModDotPlot annotations."""

    try:
        dataframe = pd.read_csv(filepath, sep="\t", comment="#", header=None, dtype=str)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=BED_COLUMNS[:3])

    if not 3 <= dataframe.shape[1] <= len(BED_COLUMNS):
        raise ValueError(
            "Invalid BED file: expected between 3 and 9 tab-separated columns."
        )

    dataframe.columns = BED_COLUMNS[: dataframe.shape[1]]
    dataframe["chrom"] = dataframe["chrom"].astype(str)

    for column in ("start", "end"):
        try:
            values = pd.to_numeric(dataframe[column], errors="raise")
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"Invalid BED file: '{column}' must contain only integers."
            ) from error
        if values.isna().any() or not np.all(np.isfinite(values)):
            raise ValueError(
                f"Invalid BED file: '{column}' must contain only finite integers."
            )
        if np.any(values % 1 != 0):
            raise ValueError(
                f"Invalid BED file: '{column}' must contain only integers."
            )
        dataframe[column] = values.astype(np.int64)

    if (dataframe["start"] < 0).any():
        raise ValueError("Invalid BED file: 'start' must be non-negative.")
    if (dataframe["end"] <= dataframe["start"]).any():
        raise ValueError(
            "Invalid BED file: 'end' must be greater than 'start' for every interval."
        )

    return dataframe


def read_annotation_beds(filepaths):
    """Read and combine one or more annotation BED files."""

    frames = [read_annotation_bed(filepath) for filepath in filepaths]
    if not frames:
        return pd.DataFrame(columns=BED_COLUMNS[:3])
    return pd.concat(frames, ignore_index=True, sort=False)


def annotation_color(value, fallback=DEFAULT_ANNOTATION_COLOR):
    """Return an RGB tuple for a BED ``itemRgb`` value."""

    if value is None or pd.isna(value):
        return fallback

    fields = [field.strip() for field in str(value).split(",")]
    if len(fields) != 3:
        return fallback

    try:
        channels = tuple(int(field) for field in fields)
    except ValueError:
        return fallback
    if any(channel < 0 or channel > 255 for channel in channels):
        return fallback
    return tuple(channel / 255 for channel in channels)


def visible_annotation_intervals(
    bed_df, chrom, region_start, region_end, fallback=DEFAULT_ANNOTATION_COLOR
):
    """Select and clip BED intervals to a plotted genomic region."""

    if region_end <= region_start:
        raise ValueError("Annotation region end must be greater than its start.")
    if bed_df.empty:
        return []

    intervals = []
    matching = bed_df[bed_df["chrom"] == str(chrom)]
    has_item_rgb = "itemRgb" in matching.columns
    for row in matching.itertuples(index=False):
        interval_start = int(row.start)
        interval_end = int(row.end)
        clipped_start = max(interval_start, region_start)
        clipped_end = min(interval_end, region_end)
        if clipped_end <= clipped_start:
            continue
        rgb = getattr(row, "itemRgb", None) if has_item_rgb else None
        intervals.append((clipped_start, clipped_end, annotation_color(rgb, fallback)))
    return intervals
