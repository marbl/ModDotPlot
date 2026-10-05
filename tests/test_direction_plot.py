import numpy as np
import pytest

from moddotplot.static_plots import create_direction_plot, direction_dataframe


def test_direction_dataframe_classifies_canonical_only_hits_as_reverse():
    canonical = np.array([[100.0, 95.0], [95.0, 100.0]])
    forward = np.array([[100.0, 0.0], [0.0, 100.0]])

    result = direction_dataframe(
        canonical,
        forward,
        window_size=10,
        name_x="chr1",
        name_y="chr1",
        self_identity=True,
        x_offset=100,
        y_offset=100,
    )

    assert result["direction"].tolist() == ["Forward", "Reverse", "Forward"]
    assert result[["q_st", "r_st"]].values.tolist() == [
        [100, 100],
        [100, 110],
        [110, 110],
    ]


def test_direction_dataframe_requires_equal_shapes():
    with pytest.raises(ValueError, match="matching shapes"):
        direction_dataframe(
            np.ones((2, 2)),
            np.ones((2, 3)),
            10,
            "x",
            "y",
            False,
        )


def test_create_direction_plot_writes_raster_and_vector_outputs(tmp_path):
    canonical = np.array([[100.0, 95.0], [95.0, 100.0]])
    forward = np.array([[100.0, 0.0], [0.0, 100.0]])

    plot = create_direction_plot(
        canonical,
        forward,
        window_size=10,
        directory=tmp_path,
        name_x="chr1",
        name_y="chr1",
        self_identity=True,
        width=1,
        dpi=72,
        vector_format="svg",
    )

    assert plot is not None
    assert (tmp_path / "chr1_DIRECTION.png").stat().st_size > 0
    assert (tmp_path / "chr1_DIRECTION.svg").stat().st_size > 0
