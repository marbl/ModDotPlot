from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import numpy as np
import pandas as pd
import pytest

from moddotplot.native_render import (
    configure_dotplot_axis,
    configure_triangle_axis,
    create_triangle_layout,
    draw_rectangular_tiles,
    draw_triangle_tiles,
    genomic_scale,
    rectangular_tile_vertices,
    save_figure_pair,
    tile_width,
    transform_triangle_points,
    triangle_tile_vertices,
)


def _tiles():
    return pd.DataFrame(
        {
            "q_st": [0, 10],
            "q_en": [4, 14],
            "r_st": [10, 20],
            "r_en": [14, 24],
            "discrete": pd.Categorical([0, 1], categories=[0, 1]),
        }
    )


def test_rectangular_tiles_preserve_center_and_max_query_width_semantics():
    data = _tiles()

    assert tile_width(data) == 4
    vertices = rectangular_tile_vertices(data)
    transposed = rectangular_tile_vertices(data, transpose=True)

    np.testing.assert_allclose(vertices[0], [[-2, 8], [2, 8], [2, 12], [-2, 12]])
    np.testing.assert_allclose(transposed[0], [[8, -2], [12, -2], [12, 2], [8, 2]])


def test_rectangular_collection_is_sparse_colored_and_optionally_rasterized():
    figure, axis = plt.subplots()
    try:
        collection = draw_rectangular_tiles(
            axis,
            _tiles(),
            {0: "#010203", 1: "#abcdef"},
            rasterized=False,
        )

        assert len(collection.get_paths()) == 2
        assert collection.get_rasterized() is False
        np.testing.assert_allclose(
            collection.get_facecolors(),
            [to_rgba("#010203"), to_rgba("#abcdef")],
        )
        assert not axis.images
    finally:
        plt.close(figure)


def test_triangle_transform_and_baseline_clipping():
    points = np.asarray([[2, 6], [4, 4], [8, 2]])
    np.testing.assert_allclose(
        transform_triangle_points(points), [[4, 2], [4, 0], [5, -3]]
    )

    diagonal = pd.DataFrame(
        {
            "q_st": [10],
            "q_en": [14],
            "r_st": [10],
            "r_en": [14],
            "discrete": [0],
        }
    )
    polygon = triangle_tile_vertices(diagonal)[0]
    assert np.min(polygon[:, 1]) == pytest.approx(0)
    assert np.max(polygon[:, 1]) == pytest.approx(2)
    assert np.min(polygon[:, 0]) == pytest.approx(8)
    assert np.max(polygon[:, 0]) == pytest.approx(12)


def test_triangle_collection_omits_tiles_wholly_below_baseline():
    data = pd.DataFrame(
        {
            "q_st": [0, 20],
            "q_en": [4, 24],
            "r_st": [10, 0],
            "r_en": [14, 4],
            "discrete": ["visible", "hidden"],
        }
    )
    figure, axis = plt.subplots()
    try:
        collection = draw_triangle_tiles(
            axis,
            data,
            {"visible": "red", "hidden": "blue"},
            rasterized=True,
        )
        assert len(collection.get_paths()) == 1
        assert collection.get_rasterized() is True
        np.testing.assert_allclose(collection.get_facecolors(), [to_rgba("red")])
        assert np.min(collection.get_paths()[0].vertices[:, 1]) >= 0
    finally:
        plt.close(figure)


def test_axis_configuration_supports_scaled_and_custom_tick_formatters():
    figure, (dot_axis, triangle_axis) = plt.subplots(1, 2)
    try:
        configure_dotplot_axis(
            dot_axis,
            0,
            1_000_000,
            breaks=[0, 500_000, 1_000_000, 1_250_000],
            formatter=lambda value, _position: f"{value / 1000:.0f}k",
        )
        configure_triangle_axis(
            triangle_axis,
            0,
            1_000_000,
            breaks=[0, 1_000_000],
        )

        assert dot_axis.get_xlim() == pytest.approx((0, 1_000_000))
        assert dot_axis.xaxis.get_major_formatter()(500_000, 0) == "500k"
        assert triangle_axis.get_ylim() == pytest.approx((0, 500_000))
        assert triangle_axis.xaxis.get_major_formatter()(1_000_000, 0) == "1"
        assert triangle_axis.get_xlabel() == "Genomic Position (Mbp)"
        assert genomic_scale(100_000) == (1_000.0, "Kbp")
        assert genomic_scale(500_000_000) == (1_000_000_000.0, "Gbp")
    finally:
        plt.close(figure)


def test_annotated_triangle_layout_shares_genomic_x_axis():
    layout = create_triangle_layout(6, with_annotation=True)
    try:
        assert layout.annotation_axis is not None
        assert layout.triangle_axis.get_shared_x_axes().joined(
            layout.triangle_axis, layout.annotation_axis
        )
        assert tuple(layout.figure.get_size_inches()) == pytest.approx((6, 3.8))
    finally:
        plt.close(layout.figure)


@pytest.mark.parametrize(
    ("vector_format", "magic"),
    [("svg", b"<?xml"), ("pdf", b"%PDF"), ("ps", b"%!PS")],
)
def test_save_figure_pair_writes_real_png_and_vector_outputs(
    tmp_path, vector_format, magic
):
    figure, axis = plt.subplots()
    axis.plot([0, 1], [0, 1])
    try:
        png_path, vector_path = save_figure_pair(
            figure, tmp_path / "nested" / "triangle", vector_format, 72
        )
    finally:
        plt.close(figure)

    assert png_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert vector_path.read_bytes().lstrip().startswith(magic)
    assert vector_path == Path(tmp_path / "nested" / f"triangle.{vector_format}")


def test_save_figure_pair_preserves_dots_in_sequence_names(tmp_path):
    figure, axis = plt.subplots()
    axis.plot([0, 1], [1, 0])
    try:
        png_path, vector_path = save_figure_pair(
            figure, tmp_path / "PAN010.chr14_TRI", "svg", 72
        )
    finally:
        plt.close(figure)

    assert png_path.name == "PAN010.chr14_TRI.png"
    assert vector_path.name == "PAN010.chr14_TRI.svg"
    assert png_path.exists()
    assert vector_path.exists()


def test_invalid_geometry_and_export_options_fail_clearly(tmp_path):
    invalid_width = _tiles()
    invalid_width.loc[0, "q_en"] = invalid_width.loc[0, "q_st"]
    with pytest.raises(ValueError, match="greater than zero"):
        rectangular_tile_vertices(invalid_width)

    figure = plt.figure()
    try:
        with pytest.raises(ValueError, match="svg, pdf, ps"):
            save_figure_pair(figure, tmp_path / "plot", "jpeg", 72)
    finally:
        plt.close(figure)
