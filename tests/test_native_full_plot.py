import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import numpy as np
import pandas as pd
from PIL import Image
import pytest

import moddotplot.static_plots as static_plots


def _processed_tiles(rows):
    frame = pd.DataFrame.from_records(
        rows,
        columns=[
            "q",
            "q_st",
            "q_en",
            "r",
            "r_st",
            "r_en",
            "perID_by_events",
            "discrete",
        ],
    )
    frame["discrete"] = pd.Categorical(
        frame["discrete"], categories=[0, 1], ordered=True
    )
    return frame


@pytest.mark.parametrize("deraster", [False, True])
def test_native_full_self_plot_mirrors_only_missing_triangle(deraster):
    frame = _processed_tiles(
        [
            ("chr1", 10, 14, "chr1", 10, 14, 100.0, 1),
            ("chr1", 10, 14, "chr1", 20, 24, 95.0, 0),
        ]
    )

    figure = static_plots._build_full_figure(
        sdf=frame,
        name_x="chr1:1-30",
        name_y="chr1:1-30",
        palette="Spectral_11",
        palette_orientation="+",
        custom_colors=["#010203", "#abcdef"],
        axes_labels=[1, 15, 30],
        xlim=(1, 30),
        deraster=deraster,
        width=3,
        is_pairwise=False,
    )
    try:
        axis = figure.axes[0]
        assert [len(collection.get_paths()) for collection in axis.collections] == [
            2,
            1,
        ]
        assert [collection.get_rasterized() for collection in axis.collections] == [
            not deraster,
            not deraster,
        ]
        np.testing.assert_allclose(
            axis.collections[1].get_paths()[0].vertices[:4],
            [[18, 8], [22, 8], [22, 12], [18, 12]],
        )
        np.testing.assert_allclose(
            axis.collections[1].get_facecolors(), [to_rgba("#010203")]
        )
        assert axis.get_xlim() == pytest.approx((1, 30))
        assert axis.get_ylim() == pytest.approx((1, 30))
        assert axis.get_xlabel() == "Genomic Position (Kbp)"
        assert axis.get_title() == "chr1"
        assert axis.get_ylabel() == "chr1"
        assert figure._suptitle.get_text() == "Self-Identity Plot: chr1"
    finally:
        plt.close(figure)


def test_native_full_self_plot_does_not_duplicate_loaded_symmetric_rows():
    frame = _processed_tiles(
        [
            ("chr1", 10, 14, "chr1", 10, 14, 100.0, 1),
            ("chr1", 10, 14, "chr1", 20, 24, 95.0, 0),
            ("chr1", 20, 24, "chr1", 10, 14, 95.0, 0),
        ]
    )

    figure = static_plots._build_full_figure(
        sdf=frame,
        name_x="chr1",
        name_y="chr1",
        palette="Spectral_11",
        palette_orientation="+",
        custom_colors=None,
        axes_labels=None,
        xlim=(1, 30),
        deraster=False,
        width=3,
        is_pairwise=False,
    )
    try:
        assert len(figure.axes[0].collections) == 1
        assert len(figure.axes[0].collections[0].get_paths()) == 3
    finally:
        plt.close(figure)


def test_native_full_plot_supports_empty_sparse_data_with_explicit_bounds():
    frame = _processed_tiles([])

    figure = static_plots._build_full_figure(
        sdf=frame,
        name_x="chrM",
        name_y="chrM",
        palette="Spectral_11",
        palette_orientation="+",
        custom_colors=None,
        axes_labels=None,
        xlim=(1, 16_569),
        deraster=False,
        width=2,
        is_pairwise=False,
    )
    try:
        axis = figure.axes[0]
        assert axis.get_xlim() == pytest.approx((1, 16_569))
        assert axis.get_ylim() == pytest.approx((1, 16_569))
        assert len(axis.collections) == 1
        assert len(axis.collections[0].get_paths()) == 0
    finally:
        plt.close(figure)


@pytest.mark.parametrize("savefig_bbox", [None, "tight"])
def test_create_plots_uses_native_comparative_renderer_and_exact_canvas(
    tmp_path, monkeypatch, savefig_bbox
):
    bed = [
        (
            "#query_name",
            "query_start",
            "query_end",
            "reference_name",
            "reference_start",
            "reference_end",
            "perID_by_events",
        ),
        ("alpha", 1, 100, "beta", 101, 200, 95.0),
    ]

    def fail_plotnine_full(*_args, **_kwargs):
        raise AssertionError("the full/compare path must not call make_dot")

    monkeypatch.setattr(static_plots, "make_dot", fail_plotnine_full)
    monkeypatch.setattr(static_plots, "make_hist", lambda *_args, **_kwargs: None)

    with plt.rc_context({"savefig.bbox": savefig_bbox}):
        outputs = static_plots.create_plots(
            sdf=[bed],
            directory=str(tmp_path),
            name_x="alpha",
            name_y="beta",
            palette="Spectral_11",
            palette_orientation="+",
            no_hist=True,
            width=2,
            dpi=40,
            is_freq=False,
            xlim=(1, 200),
            custom_colors=None,
            custom_breakpoints=None,
            from_file=None,
            is_pairwise=True,
            axes_labels=None,
            axes_tick_number=7,
            vector_format="svg",
            deraster=False,
            annotation=None,
        )

    png = tmp_path / "alpha_beta_COMPARE.png"
    svg = tmp_path / "alpha_beta_COMPARE.svg"
    assert set(outputs) == {str(png), str(svg)}
    assert png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert svg.read_bytes().lstrip().startswith(b"<?xml")
    with Image.open(png) as image:
        assert image.size == (80, 80)
def test_native_comparative_direction_colors_preserve_hue_and_ani_strength():
    frame = _processed_tiles(
        [
            ("query", 10, 14, "reference", 10, 14, 90.0, 0),
            ("query", 20, 24, "reference", 20, 24, 99.0, 1),
        ]
    )
    frame["direction"] = pd.Categorical(
        ["Forward", "Reverse"], categories=["Forward", "Reverse"], ordered=True
    )

    figure = static_plots._build_full_figure(
        sdf=frame,
        name_x="query",
        name_y="reference",
        palette="Spectral_11",
        palette_orientation="+",
        custom_colors=None,
        axes_labels=None,
        xlim=(1, 30),
        deraster=False,
        width=3,
        is_pairwise=True,
    )
    try:
        colors = figure.axes[0].collections[0].get_facecolors()
        assert colors.shape == (2, 4)
        # Forward remains blue-dominant; reverse remains red-dominant. The
        # stronger second ANI bin is darker than the weak first bin.
        assert colors[0, 2] > colors[0, 0]
        assert colors[1, 0] > colors[1, 2]
        assert np.mean(colors[1, :3]) < np.mean(colors[0, :3])
        assert figure.axes[0].get_title() == "query"
        assert figure.axes[0].get_ylabel() == "reference"
    finally:
        plt.close(figure)
