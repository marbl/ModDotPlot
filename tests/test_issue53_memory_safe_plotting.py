import matplotlib.pyplot as plt
import pandas as pd
from plotnine.geoms.geom_tile import geom_tile

from moddotplot.static_plots import make_dot, make_dot_final, make_dot_grid, make_tri

GENOME_SIZE = 496_000_000
WINDOW_SIZE = 2_000


def _large_sparse_dotplot_data():
    # The adjacent first two coordinates establish a 2 kb resolution while the
    # final coordinate establishes the ~496 Mb extent from issue #53.  A
    # geom_raster layer would try to allocate about 248,000**2 RGBA pixels.
    starts = [0, WINDOW_SIZE, GENOME_SIZE - WINDOW_SIZE]
    return pd.DataFrame(
        {
            "q": ["chr1"] * len(starts),
            "q_st": starts,
            "q_en": [start + WINDOW_SIZE - 1 for start in starts],
            "r": ["chr1"] * len(starts),
            "r_st": starts,
            "r_en": [start + WINDOW_SIZE - 1 for start in starts],
            "discrete": pd.Categorical([0, 1, 2]),
        }
    )


def _make_full_plot(data, deraster=False):
    return make_dot(
        sdf=data,
        name_x="chr1",
        name_y="chr1",
        palette="Spectral_11",
        palette_orientation="+",
        colors=None,
        breaks=[0, GENOME_SIZE // 2, GENOME_SIZE],
        num_ticks=3,
        xlim=GENOME_SIZE,
        deraster=deraster,
        width=2,
        is_pairwise=False,
    )


def _assert_tile_layer(plot, rasterized):
    layer = plot.layers[0]
    assert isinstance(layer.geom, geom_tile)
    assert layer.geom._kwargs["raster"] is rasterized


def test_large_sparse_plot_does_not_build_coordinate_sized_raster():
    plot = _make_full_plot(_large_sparse_dotplot_data())

    _assert_tile_layer(plot, rasterized=True)
    figure = plot.draw(show=False)
    try:
        axis = figure.axes[0]
        assert not axis.images
        assert len(axis.collections) == 1
        assert axis.collections[0].get_rasterized() is True
        assert len(axis.collections[0].get_paths()) == 3
    finally:
        plt.close(figure)


def test_deraster_keeps_memory_safe_tiles_as_vectors():
    plot = _make_full_plot(_large_sparse_dotplot_data(), deraster=True)

    _assert_tile_layer(plot, rasterized=False)
    figure = plot.draw(show=False)
    try:
        assert figure.axes[0].collections[0].get_rasterized() is False
    finally:
        plt.close(figure)


def test_grid_and_triangle_paths_use_the_same_memory_safe_geometry():
    data = _large_sparse_dotplot_data()
    common = {
        "sdf": data,
        "palette": "Spectral_11",
        "palette_orientation": "+",
        "colors": None,
        "breaks": [0, GENOME_SIZE // 2, GENOME_SIZE],
        "xlim": GENOME_SIZE,
        "deraster": False,
        "width": 2,
    }

    grid_cell = make_dot_final(**common)
    grid_plot = make_dot_grid(
        **common,
        title_name="grid",
        on_diagonal=True,
    )
    triangle, _ = make_tri(
        **common,
        title_name="triangle",
        num_ticks=3,
    )

    for plot in (grid_cell, grid_plot, triangle):
        _assert_tile_layer(plot, rasterized=True)
