import pandas as pd
import pytest
from matplotlib.colors import to_rgba

from moddotplot.const import DIRECTION_COLORS
from moddotplot.static_plots import (
    display_sequence_name,
    generate_breaks,
    get_colors,
    make_dot,
)


def test_get_colors_uses_string_custom_breakpoints():
    identity_scores = pd.DataFrame(
        {
            "perID_by_events": [
                86.0,
                89.0,
                90.0,
                97.6,
                98.1,
                98.9,
                99.6,
                100.0,
            ]
        }
    )
    custom_breakpoints = [
        "86",
        "90",
        "97.5",
        "97.75",
        "98.0",
        "98.25",
        "98.5",
        "98.75",
        "99.0",
        "99.25",
        "99.5",
        "100.0",
    ]

    bins = get_colors(
        identity_scores,
        ncolors=11,
        is_freq=False,
        custom_breakpoints=custom_breakpoints,
    )

    assert bins.astype(int).tolist() == [0, 0, 0, 2, 4, 7, 10, 10]


def test_make_dot_uses_custom_color_scale():
    custom_colors = ["#010203", "#456789", "#abcdef"]
    plot_data = pd.DataFrame(
        {
            "q": ["query"] * 3,
            "q_st": [0, 10, 20],
            "q_en": [10, 20, 30],
            "r": ["reference"] * 3,
            "r_st": [0, 10, 20],
            "r_en": [10, 20, 30],
            "discrete": pd.Categorical([0, 1, 2], categories=[0, 1, 2]),
        }
    )

    plot = make_dot(
        sdf=plot_data,
        name_x="query",
        name_y="reference",
        palette="Spectral_11",
        palette_orientation="+",
        colors=custom_colors,
        breaks=[0, 10, 20, 30],
        num_ticks=4,
        xlim=30,
        deraster=False,
        width=4,
        is_pairwise=True,
    )

    fill_scale = plot.scales.get_scales("fill")
    assert fill_scale.palette(len(custom_colors)) == custom_colors


def test_make_dot_uses_direction_colors_when_orientation_is_present():
    plot_data = pd.DataFrame(
        {
            "q": ["query"] * 4,
            "q_st": [0, 10, 20, 30],
            "q_en": [10, 20, 30, 40],
            "r": ["reference"] * 4,
            "r_st": [0, 10, 20, 30],
            "r_en": [10, 20, 30, 40],
            "discrete": pd.Categorical([0, 1, 0, 1], categories=[0, 1]),
            "direction": ["Forward", "Forward", "Reverse", "Reverse"],
        }
    )

    plot = make_dot(
        sdf=plot_data,
        name_x="query",
        name_y="reference",
        palette="Spectral_11",
        palette_orientation="+",
        colors=["#000000", "#ffffff"],
        breaks=[0, 10, 20, 30, 40],
        num_ticks=3,
        xlim=40,
        deraster=True,
        width=4,
        is_pairwise=True,
    )

    figure = plot.draw(show=False)
    try:
        rendered = {
            tuple(color)
            for collection in figure.axes[0].collections
            for color in collection.get_facecolors()
        }
        assert to_rgba(DIRECTION_COLORS["Forward"]) in rendered
        assert to_rgba(DIRECTION_COLORS["Reverse"]) in rendered
        assert len(rendered) == 4
        assert to_rgba("#000000") not in rendered
    finally:
        import matplotlib.pyplot as plt

        plt.close(figure)


def test_make_dot_honors_exact_region_bounds():
    plot_data = pd.DataFrame(
        {
            "q": ["query"],
            "q_st": [200],
            "q_en": [300],
            "r": ["reference"],
            "r_st": [200],
            "r_en": [300],
            "discrete": pd.Categorical([0]),
        }
    )

    plot = make_dot(
        sdf=plot_data,
        name_x="query",
        name_y="reference",
        palette="Spectral_11",
        palette_orientation="+",
        colors=None,
        breaks=None,
        num_ticks=4,
        xlim=(101, 400),
        deraster=False,
        width=4,
        is_pairwise=True,
    )

    assert plot.scales.get_scales("x").limits == (101.0, 400.0)
    assert plot.scales.get_scales("y").limits == (101.0, 400.0)


def test_display_names_omit_region_and_full_axes_are_twice_as_large():
    name = "PAN010.chr14.haplotype1.paternal:1-4000000"
    plot_data = pd.DataFrame(
        {
            "q": [name],
            "q_st": [1],
            "q_en": [100],
            "r": [name],
            "r_st": [1],
            "r_en": [100],
            "discrete": pd.Categorical([0]),
        }
    )

    plot = make_dot(
        sdf=plot_data,
        name_x=name,
        name_y=name,
        palette="Spectral_11",
        palette_orientation="+",
        colors=None,
        breaks=[1, 50, 100],
        num_ticks=3,
        xlim=(1, 100),
        deraster=False,
        width=4,
        is_pairwise=False,
    )

    assert display_sequence_name(name) == "PAN010.chr14.haplotype1.paternal"
    assert ":1-4000000" not in plot.labels.title
    assert plot.data["q"].unique().tolist() == ["PAN010.chr14.haplotype1.paternal"]

    figure = plot.draw(show=False)
    try:
        axis = figure.axes[0]
        assert axis.get_xticklabels()[0].get_fontsize() == pytest.approx(8)
        # Plotnine rounds text sizes to whole points: 2 * (width * 1.4) = 11.2.
        assert axis.xaxis.label.get_fontsize() == pytest.approx(11)
    finally:
        import matplotlib.pyplot as plt

        plt.close(figure)


def test_generated_breaks_never_extend_past_sequence_length():
    breaks = generate_breaks(1, 103_156_783)

    assert breaks
    assert all(1 <= value <= 103_156_783 for value in breaks)
    assert breaks[-1] == 100_000_000


def test_get_colors_allows_large_custom_palettes():
    scores = pd.DataFrame({"perID_by_events": [0.0, 50.0, 100.0]})
    breaks = list(range(0, 105, 5))

    bins = get_colors(scores, ncolors=20, is_freq=False, custom_breakpoints=breaks)

    assert bins.astype(int).tolist() == [0, 9, 19]


@pytest.mark.parametrize("is_freq", [False, True])
def test_get_colors_handles_a_single_perfect_identity_bin(is_freq):
    scores = pd.DataFrame({"perID_by_events": [100.0, 100.0, 100.0]})

    bins = get_colors(
        scores,
        ncolors=11,
        is_freq=is_freq,
        custom_breakpoints=None,
    )

    assert bins.tolist() == [0, 0, 0]


@pytest.mark.parametrize(
    ("breakpoints", "message"),
    [
        ([0, 50], "number of breakpoints"),
        ([0, 75, 50, 100], "strictly increasing"),
        ([0, float("nan"), 50, 100], "finite"),
    ],
)
def test_get_colors_rejects_invalid_custom_breakpoints(breakpoints, message):
    scores = pd.DataFrame({"perID_by_events": [0.0, 50.0, 100.0]})

    with pytest.raises(ValueError, match=message):
        get_colors(
            scores,
            ncolors=3,
            is_freq=False,
            custom_breakpoints=breakpoints,
        )


def test_get_colors_rejects_breakpoints_that_do_not_cover_observed_values():
    scores = pd.DataFrame({"perID_by_events": [86.0, 100.0]})

    with pytest.raises(ValueError, match="cover all finite identity values"):
        get_colors(
            scores,
            ncolors=2,
            is_freq=False,
            custom_breakpoints=[86, 90, 95],
        )
