import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import pandas as pd
import pytest

import moddotplot.static_plots as static_plots
from moddotplot.const import DIRECTION_COLORS
from moddotplot.static_plots import (
    DEFAULT_ANNOTATION_COLOR,
    _build_triangle_figure,
    draw_annotation_track,
    read_annotation_bed,
    render_annotation_track,
)


@pytest.mark.parametrize(
    ("record", "expected_columns"),
    [
        ("chr1\t10\t20\n", ["chrom", "start", "end"]),
        (
            "chr1\t10\t20\tfeature\t42\t+\t11\t19\t12,34,56\n",
            [
                "chrom",
                "start",
                "end",
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
            ],
        ),
    ],
)
def test_read_annotation_bed_accepts_bed3_through_bed9(
    tmp_path, record, expected_columns
):
    bed_path = tmp_path / "annotations.bed"
    bed_path.write_text(record)

    dataframe = read_annotation_bed(bed_path)

    assert list(dataframe.columns) == expected_columns
    assert dataframe.loc[0, "chrom"] == "chr1"
    assert dataframe.loc[0, "start"] == 10
    assert dataframe.loc[0, "end"] == 20


@pytest.mark.parametrize(
    "record",
    [
        "chr1\tstart\t20\n",
        "chr1\t10.5\t20\n",
        "chr1\t-1\t20\n",
        "chr1\t20\t20\n",
        "chr1\t21\t20\n",
        "chr1\t1\t2\ta\t0\t+\t1\t2\t0,0,0\textra\n",
    ],
)
def test_read_annotation_bed_rejects_malformed_records(tmp_path, record):
    bed_path = tmp_path / "invalid.bed"
    bed_path.write_text(record)

    with pytest.raises(ValueError, match="Invalid BED file"):
        read_annotation_bed(bed_path)


def test_draw_annotation_track_filters_clips_and_uses_rgb(tmp_path):
    bed_path = tmp_path / "annotations.bed"
    bed_path.write_text(
        "chr1\t0\t12\tleft\t0\t+\t0\t12\t255,0,0\n"
        "chr1\t18\t30\tright\t0\t+\t18\t30\tnot-a-color\n"
        "chr1\t20\t25\toutside\t0\t+\t20\t25\t0,255,0\n"
        "chr2\t10\t20\twrong-chrom\t0\t+\t10\t20\t0,0,255\n"
    )
    dataframe = read_annotation_bed(bed_path)
    figure, axis = plt.subplots()

    interval_count = draw_annotation_track(axis, dataframe, "chr1", 10, 20)

    assert interval_count == 2
    assert axis.get_xlim() == pytest.approx((10, 20))
    assert len(axis.patches) == 2
    assert axis.patches[0].get_x() == 10
    assert axis.patches[0].get_width() == 2
    assert axis.patches[0].get_facecolor() == pytest.approx(to_rgba((1, 0, 0)))
    assert axis.patches[1].get_x() == 18
    assert axis.patches[1].get_width() == 2
    assert axis.patches[1].get_facecolor() == pytest.approx(
        to_rgba(DEFAULT_ANNOTATION_COLOR)
    )
    assert {patch.get_y() for patch in axis.patches} == {0.2}
    plt.close(figure)


def test_render_annotation_track_writes_expected_svg_and_png(tmp_path):
    bed_path = tmp_path / "annotations.bed"
    bed_path.write_text("chr1\t10\t20\n")
    dataframe = read_annotation_bed(bed_path)
    output_prefix = tmp_path / "sample_ANNOTATION_TRACK"

    created = render_annotation_track(
        dataframe,
        "chr1",
        0,
        100,
        output_prefix,
        width_cm=20,
        dpi=72,
    )

    svg_path = tmp_path / "sample_ANNOTATION_TRACK.svg"
    png_path = tmp_path / "sample_ANNOTATION_TRACK.png"
    assert created is True
    assert svg_path.stat().st_size > 0
    assert png_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    ET.parse(svg_path)


def test_triangle_uses_blue_and_pink_direction_colors():
    dataframe = pd.DataFrame(
        {
            "q": ["chr1"] * 4,
            "q_st": [10, 20, 30, 40],
            "q_en": [20, 30, 40, 50],
            "r": ["chr1"] * 4,
            "r_st": [10, 40, 60, 80],
            "r_en": [20, 50, 70, 90],
            "discrete": pd.Categorical([0, 1, 0, 1], categories=[0, 1]),
            "direction": ["Forward", "Forward", "Reverse", "Reverse"],
        }
    )

    figure = _build_triangle_figure(
        sdf=dataframe,
        title="chr1",
        palette="Spectral_11",
        palette_orientation="+",
        custom_colors=["#000000", "#ffffff"],
        axes_labels=[0, 50, 100],
        xlim=(0, 100),
        deraster=True,
        width=4,
    )
    try:
        rendered = {
            tuple(color)
            for collection in figure.axes[0].collections
            for color in collection.get_facecolors()
        }
        assert to_rgba(DIRECTION_COLORS["Forward"]) in rendered
        assert to_rgba(DIRECTION_COLORS["Reverse"]) in rendered
        assert len(rendered) == 4
    finally:
        plt.close(figure)


def test_render_annotation_track_skips_empty_overlap(tmp_path):
    bed_path = tmp_path / "annotations.bed"
    bed_path.write_text("chr2\t10\t20\n")
    dataframe = read_annotation_bed(bed_path)
    output_prefix = tmp_path / "sample_ANNOTATION_TRACK"

    created = render_annotation_track(
        dataframe,
        "chr1",
        0,
        100,
        output_prefix,
        width_cm=20,
        dpi=72,
    )

    assert created is False
    assert not (tmp_path / "sample_ANNOTATION_TRACK.svg").exists()
    assert not (tmp_path / "sample_ANNOTATION_TRACK.png").exists()


def test_annotated_triangle_physically_aligns_axes_and_hides_heatmap_baseline(
    tmp_path,
):
    bed_path = tmp_path / "annotations.bed"
    bed_path.write_text("chr1\t0\t100\n")
    dataframe = pd.DataFrame(
        {
            "q_st": [0],
            "q_en": [100],
            "r_st": [0],
            "r_en": [100],
            "discrete": [0],
        }
    )

    figure = _build_triangle_figure(
        sdf=dataframe,
        title="chr1:1-100",
        palette="Spectral_11",
        palette_orientation="+",
        custom_colors=None,
        axes_labels=None,
        xlim=(0, 100),
        deraster=False,
        width=6,
        annotation_df=read_annotation_bed(bed_path),
        annotation_chrom="chr1",
    )
    try:
        triangle_axis, annotation_axis = figure.axes
        triangle_position = triangle_axis.get_position()
        annotation_position = annotation_axis.get_position()

        assert annotation_position.x0 == pytest.approx(triangle_position.x0)
        assert annotation_position.width == pytest.approx(triangle_position.width)
        assert not triangle_axis.spines["bottom"].get_visible()
        assert not any(
            tick.tick1line.get_visible() for tick in triangle_axis.xaxis.majorTicks
        )
        assert annotation_axis.spines["bottom"].get_visible()
        assert annotation_axis.get_xlabel() == "Genomic Position (Kbp)"
    finally:
        plt.close(figure)


def test_unannotated_triangle_keeps_its_genomic_axis():
    dataframe = pd.DataFrame(
        {
            "q_st": [0],
            "q_en": [100],
            "r_st": [0],
            "r_en": [100],
            "discrete": [0],
        }
    )

    figure = _build_triangle_figure(
        sdf=dataframe,
        title="chr1",
        palette="Spectral_11",
        palette_orientation="+",
        custom_colors=None,
        axes_labels=None,
        xlim=(0, 100),
        deraster=False,
        width=6,
    )
    try:
        triangle_axis = figure.axes[0]
        assert triangle_axis.spines["bottom"].get_visible()
        assert triangle_axis.get_xlabel() == "Genomic Position (Kbp)"
    finally:
        plt.close(figure)


class _FakePlot:
    def __add__(self, _other):
        return self


def _stub_create_plots_dependencies(monkeypatch, *, directional=False):
    dataframe = pd.DataFrame(
        [
            {
                "q": "chr1",
                "q_st": 0,
                "q_en": 100,
                "r": "chr1",
                "r_st": 0,
                "r_en": 100,
                "perID_by_events": 100.0,
                "discrete": 0,
            }
        ]
    )
    if directional:
        dataframe["direction"] = ["Forward"]
    monkeypatch.setattr(static_plots, "read_df", lambda *_args, **_kwargs: dataframe)
    monkeypatch.setattr(
        static_plots, "make_hist", lambda *_args, **_kwargs: _FakePlot()
    )
    monkeypatch.setattr(static_plots, "make_dot", lambda *_args, **_kwargs: _FakePlot())

    def fake_plot_pair(
        _plot,
        output_prefix,
        *,
        width,
        height,
        dpi,
        vector_format,
    ):
        del width, height, dpi
        png = Path(f"{output_prefix}.png")
        vector = Path(f"{output_prefix}.{vector_format}")
        png.write_bytes(b"png")
        vector.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')
        return png, vector

    monkeypatch.setattr(static_plots, "_draw_and_save_plot_pair", fake_plot_pair)


def _run_create_plots(output_dir, annotation, vector_format="svg"):
    return static_plots.create_plots(
        sdf=None,
        directory=str(output_dir),
        name_x="chr1",
        name_y="chr1",
        palette="Spectral_11",
        palette_orientation="+",
        no_hist=True,
        width=4,
        dpi=72,
        is_freq=False,
        xlim=100,
        custom_colors=None,
        custom_breakpoints=None,
        from_file=None,
        is_pairwise=False,
        axes_labels=None,
        axes_tick_number=7,
        vector_format=vector_format,
        deraster=True,
        annotation=str(annotation),
    )


def test_create_plots_creates_and_skips_annotation_artifacts(tmp_path, monkeypatch):
    _stub_create_plots_dependencies(monkeypatch)

    matching_dir = tmp_path / "matching"
    matching_dir.mkdir()
    matching_bed = tmp_path / "matching.bed"
    matching_bed.write_text("chr1\t10\t20\n")
    matching_outputs = _run_create_plots(matching_dir, matching_bed)

    assert (matching_dir / "chr1_ANNOTATION_TRACK.svg").exists()
    assert (matching_dir / "chr1_ANNOTATION_TRACK.png").exists()
    assert (matching_dir / "chr1_TRI_ANNOTATED.svg").exists()
    assert (matching_dir / "chr1_TRI_ANNOTATED.png").exists()
    assert not (matching_dir / "chr1_PRE_ANNOTATED.svg").exists()
    assert str(matching_dir / "chr1_TRI_ANNOTATED.svg") in matching_outputs
    assert str(matching_dir / "chr1_ANNOTATION_TRACK.svg") in matching_outputs

    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    nonmatching_bed = tmp_path / "nonmatching.bed"
    nonmatching_bed.write_text("chr2\t10\t20\n")
    empty_outputs = _run_create_plots(empty_dir, nonmatching_bed)

    assert not (empty_dir / "chr1_ANNOTATION_TRACK.svg").exists()
    assert not (empty_dir / "chr1_ANNOTATION_TRACK.png").exists()
    assert not (empty_dir / "chr1_PRE_ANNOTATED.svg").exists()
    assert not (empty_dir / "chr1_TRI_ANNOTATED.svg").exists()
    assert not (empty_dir / "chr1_TRI_ANNOTATED.png").exists()
    assert not any("ANNOTATION" in output for output in empty_outputs)


def test_create_plots_uses_direction_output_names(tmp_path, monkeypatch):
    _stub_create_plots_dependencies(monkeypatch, directional=True)

    outputs = static_plots.create_plots(
        sdf=None,
        directory=str(tmp_path / "directionality"),
        name_x="chr1",
        name_y="chr1",
        palette="Spectral_11",
        palette_orientation="+",
        no_hist=False,
        width=4,
        dpi=72,
        is_freq=False,
        xlim=100,
        custom_colors=None,
        custom_breakpoints=None,
        from_file=None,
        is_pairwise=False,
        axes_labels=None,
        axes_tick_number=7,
        vector_format="svg",
        deraster=True,
        annotation=None,
    )

    expected_stems = {
        "chr1_DIRECTION_FULL",
        "chr1_DIRECTION_TRI",
        "chr1_DIRECTION_HIST",
    }
    assert {Path(output).stem for output in outputs} == expected_stems
    assert all(Path(output).parent.name == "directionality" for output in outputs)


@pytest.mark.parametrize(
    ("vector_format", "magic"),
    [("svg", b"<svg"), ("pdf", b"%PDF"), ("ps", b"%!PS")],
)
def test_create_plots_writes_native_triangle_and_annotation_formats(
    tmp_path, monkeypatch, vector_format, magic
):
    _stub_create_plots_dependencies(monkeypatch)
    annotation = tmp_path / "annotation.bed"
    annotation.write_text("chr1\t10\t20\n")

    _run_create_plots(tmp_path, annotation, vector_format=vector_format)

    for stem in ("chr1_TRI", "chr1_TRI_ANNOTATED", "chr1_ANNOTATION_TRACK"):
        assert (tmp_path / f"{stem}.png").read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        vector = (tmp_path / f"{stem}.{vector_format}").read_bytes().lstrip()
        if vector_format == "svg":
            assert magic in vector[:1024]
        else:
            assert vector.startswith(magic)

    if vector_format != "svg":
        assert not list(tmp_path.glob("*_TRI*.svg"))
        assert not list(tmp_path.glob("*_ANNOTATION_TRACK.svg"))
