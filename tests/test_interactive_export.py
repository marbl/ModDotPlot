import pytest

from moddotplot.interactive import figure_to_bed


def _figure(z, *, x=(100, 110), y=(200, 210), x_name="chr1", y_name="chr2"):
    return {
        "data": [{"z": z, "x": list(x), "y": list(y)}],
        "layout": {
            "xaxis": {"title": {"text": x_name}},
            "yaxis": {"title": {"text": y_name}},
        },
    }


def test_figure_to_bed_includes_coordinate_offsets():
    rows, filename = figure_to_bed(_figure([[98.0, 0.0], [91.0, 99.0]]))

    assert filename == "chr1-chr2.bedpe"
    assert rows[1][:6] == ("chr1", 100, 109, "chr2", 200, 209)
    assert rows[-1][:6] == ("chr1", 110, 119, "chr2", 210, 219)


def test_figure_to_bed_names_self_identity_export():
    rows, filename = figure_to_bed(
        _figure([[100.0, 95.0], [95.0, 100.0]], x_name="chr1", y_name="chr1")
    )

    assert filename == "chr1.bedpe"
    assert len(rows) == 4  # header plus the upper triangle


def test_figure_to_bed_rejects_missing_window_size():
    with pytest.raises(ValueError, match="positive window size"):
        figure_to_bed(_figure([[100.0]], x=(), y=()))
