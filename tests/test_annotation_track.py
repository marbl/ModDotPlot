import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import pandas as pd
import pytest

import moddotplot.static_plots as static_plots
from moddotplot.static_plots import (
    DEFAULT_ANNOTATION_COLOR,
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


class _FakePlot:
    def __add__(self, _other):
        return self


def _stub_create_plots_dependencies(monkeypatch):
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
    monkeypatch.setattr(static_plots, "read_df", lambda *_args, **_kwargs: dataframe)
    monkeypatch.setattr(
        static_plots, "make_hist", lambda *_args, **_kwargs: _FakePlot()
    )
    monkeypatch.setattr(static_plots, "make_dot", lambda *_args, **_kwargs: _FakePlot())

    def fake_ggsave(*_args, **kwargs):
        output = Path(kwargs["filename"])
        if output.suffix == ".png":
            output.write_bytes(b"png")
        else:
            output.write_text('<svg xmlns="http://www.w3.org/2000/svg"/>')

    monkeypatch.setattr(static_plots, "ggsave", fake_ggsave)


def _run_create_plots(output_dir, annotation, vector_format="svg"):
    static_plots.create_plots(
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
    _run_create_plots(matching_dir, matching_bed)

    assert (matching_dir / "chr1_ANNOTATION_TRACK.svg").exists()
    assert (matching_dir / "chr1_ANNOTATION_TRACK.png").exists()
    assert (matching_dir / "chr1_TRI_ANNOTATED.svg").exists()
    assert (matching_dir / "chr1_TRI_ANNOTATED.png").exists()
    assert not (matching_dir / "chr1_PRE_ANNOTATED.svg").exists()

    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    nonmatching_bed = tmp_path / "nonmatching.bed"
    nonmatching_bed.write_text("chr2\t10\t20\n")
    _run_create_plots(empty_dir, nonmatching_bed)

    assert not (empty_dir / "chr1_ANNOTATION_TRACK.svg").exists()
    assert not (empty_dir / "chr1_ANNOTATION_TRACK.png").exists()
    assert not (empty_dir / "chr1_PRE_ANNOTATED.svg").exists()
    assert not (empty_dir / "chr1_TRI_ANNOTATED.svg").exists()
    assert not (empty_dir / "chr1_TRI_ANNOTATED.png").exists()


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
