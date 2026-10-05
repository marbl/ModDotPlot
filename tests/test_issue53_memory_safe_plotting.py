import matplotlib.pyplot as plt
import pandas as pd

from moddotplot.static_plots import make_dot, make_dot_final, make_dot_grid, make_tri

GENOME_SIZE = 496_000_000
WINDOW_SIZE = 2_000


def _large_sparse_dotplot_data():
    # The adjacent first two coordinates establish a 2 kb resolution while the
    # final coordinate establishes the ~496 Mb extent from issue #53.  A
    # A coordinate-sized raster would allocate about 248,000**2 RGBA pixels.
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


def _assert_tile_collection(figure, rasterized):
    axis = figure.axes[0]
    assert not axis.images
    assert axis.collections
    assert axis.collections[0].get_rasterized() is rasterized


def test_large_sparse_plot_does_not_build_coordinate_sized_raster():
    figure = _make_full_plot(_large_sparse_dotplot_data())

    try:
        _assert_tile_collection(figure, rasterized=True)
        axis = figure.axes[0]
        assert len(axis.collections) == 1
        assert len(axis.collections[0].get_paths()) == 3
    finally:
        plt.close(figure)


def test_deraster_keeps_memory_safe_tiles_as_vectors():
    figure = _make_full_plot(_large_sparse_dotplot_data(), deraster=True)

    try:
        _assert_tile_collection(figure, rasterized=False)
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

    for figure in (grid_cell, grid_plot, triangle):
        try:
            _assert_tile_collection(figure, rasterized=True)
        finally:
            plt.close(figure)
