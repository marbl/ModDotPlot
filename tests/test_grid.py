import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
import numpy as np
from pathlib import Path
import pytest

import moddotplot.static_plots as static_plots
from moddotplot.const import DIRECTION_COLORS
from moddotplot.native_render import FALLBACK_FONT_FAMILY, set_figure_font_family
from moddotplot.static_plots import _build_grid_figure, create_grid

BED_HEADER = (
    "#query_name",
    "query_start",
    "query_end",
    "reference_name",
    "reference_start",
    "reference_end",
    "perID_by_events",
)


def _bed(matrix, query_name, reference_name, *, self_identity):
    """Build the in-memory BEDPE representation accepted by ``create_grid``."""
    rows = [BED_HEADER]
    for query_index, matrix_row in enumerate(matrix):
        for reference_index, identity in enumerate(matrix_row):
            if self_identity and query_index > reference_index:
                continue
            if identity < 86:
                continue
            query_start = query_index * 10
            reference_start = reference_index * 10
            rows.append(
                (
                    query_name,
                    query_start,
                    query_start + 9,
                    reference_name,
                    reference_start,
                    reference_start + 9,
                    float(identity),
                )
            )
    return rows


def _records(query_name, reference_name, cells):
    """Build BEDPE rows from ``(q_start, r_start, identity)`` cells."""
    return [BED_HEADER] + [
        (
            query_name,
            query_start,
            query_start + 10,
            reference_name,
            reference_start,
            reference_start + 10,
            float(identity),
        )
        for query_start, reference_start, identity in cells
    ]


def _grid_kwargs(
    *,
    singles,
    doubles,
    single_names,
    double_names,
    custom_colors=None,
    custom_breakpoints=None,
    deraster=False,
):
    return {
        "singles": singles,
        "doubles": doubles,
        "palette": "Spectral_11",
        "palette_orientation": "+",
        "single_names": single_names,
        "double_names": double_names,
        "is_freq": False,
        "xlim": 100,
        "custom_colors": custom_colors,
        "custom_breakpoints": custom_breakpoints,
        "axes_label": [0, 50, 100],
        "is_bed": False,
        "width": 1,
        "breaks": [0, 50, 100],
        "deraster": deraster,
    }


def _basic_two_sequence_grid(**overrides):
    names = ["sequence_a", "sequence_b"]
    kwargs = _grid_kwargs(
        singles=[
            _records(names[0], names[0], [(10, 20, 91)]),
            _records(names[1], names[1], [(30, 40, 92)]),
        ],
        doubles=[_records(names[0], names[1], [(20, 70, 95)])],
        single_names=names,
        double_names=[[names[0], names[1]]],
    )
    kwargs.update(overrides)
    return kwargs


def _create_grid(tmp_path, *, vector_format="svg", **kwargs):
    return create_grid(
        directory=tmp_path,
        vector_format=vector_format,
        dpi=72,
        **kwargs,
    )


def _collection_center(axis):
    vertices = [
        path.vertices
        for collection in axis.collections
        for path in collection.get_paths()
        if path.vertices.size
    ]
    assert vertices, "expected the grid cell to contain a plotted collection"
    points = np.concatenate(vertices)
    return (
        (points[:, 0].min() + points[:, 0].max()) / 2,
        (points[:, 1].min() + points[:, 1].max()) / 2,
    )


def _collection_centers(axis):
    return sorted(
        (
            (path.vertices[:, 0].min() + path.vertices[:, 0].max()) / 2,
            (path.vertices[:, 1].min() + path.vertices[:, 1].max()) / 2,
        )
        for collection in axis.collections
        for path in collection.get_paths()
        if path.vertices.size
    )


def _all_artist_colors(axes):
    colors = set()
    for axis in axes.flat:
        for collection in axis.collections:
            colors.update(tuple(color) for color in collection.get_facecolors())
        colors.update(to_rgba(patch.get_facecolor()) for patch in axis.patches)
    return colors


def test_create_grid_handles_empty_pairwise_comparison(tmp_path):
    names = ["sequence_a", "sequence_b", "sequence_c"]
    self_matrix = [[100, 92], [92, 100]]
    pair_matrix = [[91, 94], [96, 99]]
    empty_pair_matrix = [[0, 0], [0, 0]]

    singles = [_bed(self_matrix, name, name, self_identity=True) for name in names]
    doubles = [
        _bed(pair_matrix, "sequence_a", "sequence_b", self_identity=False),
        _bed(pair_matrix, "sequence_a", "sequence_c", self_identity=False),
        _bed(
            empty_pair_matrix,
            "sequence_b",
            "sequence_c",
            self_identity=False,
        ),
    ]

    # An all-zero comparison produces a BED table containing only its header.
    assert len(doubles[-1]) == 1

    _create_grid(
        tmp_path,
        **_grid_kwargs(
            singles=singles,
            doubles=doubles,
            single_names=names,
            double_names=[
                ["sequence_a", "sequence_b"],
                ["sequence_a", "sequence_c"],
                ["sequence_b", "sequence_c"],
            ],
        ),
    )

    for suffix in ("png", "svg"):
        output = tmp_path / f"3x3_GRID.{suffix}"
        assert output.is_file()
        assert output.stat().st_size > 0


@pytest.mark.parametrize(
    ("vector_format", "signature"),
    [("svg", b"<svg"), ("pdf", b"%PDF-"), ("ps", b"%!PS-Adobe")],
)
def test_create_grid_writes_png_and_valid_selected_vector_format(
    tmp_path, vector_format, signature
):
    output_files = _create_grid(
        tmp_path,
        vector_format=vector_format,
        **_basic_two_sequence_grid(),
    )

    png = tmp_path / "2x2_GRID.png"
    vector = tmp_path / f"2x2_GRID.{vector_format}"
    assert output_files == [str(vector), str(png)]
    assert png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    if vector_format == "svg":
        assert signature in vector.read_bytes()[:1024]
    else:
        assert vector.read_bytes().startswith(signature)


def test_create_grid_uses_direction_output_name(tmp_path):
    kwargs = _basic_two_sequence_grid()
    direction_header = (*BED_HEADER, "direction")
    kwargs["singles"] = [
        [direction_header] + [(*row, "Forward") for row in records[1:]]
        for records in kwargs["singles"]
    ]
    kwargs["doubles"] = [
        [direction_header] + [(*row, "Reverse") for row in records[1:]]
        for records in kwargs["doubles"]
    ]

    outputs = _create_grid(tmp_path / "directionality", **kwargs)

    assert {Path(output).name for output in outputs} == {
        "2x2_DIRECTION_GRID.svg",
        "2x2_DIRECTION_GRID.png",
    }


def test_reversed_pair_metadata_orients_asymmetric_coordinates_by_grid_cell():
    # The BED record is B (query) vs A (reference), while the requested grid is
    # ordered A, B across columns and B, A down rows. Pairwise comparisons sit
    # on the other diagonal and must remain exact transposes of one another.
    kwargs = _basic_two_sequence_grid(
        doubles=[_records("sequence_b", "sequence_a", [(20, 70, 95)])],
        double_names=[["sequence_b", "sequence_a"]],
    )

    figure, axes = _build_grid_figure(**kwargs)
    try:
        assert axes.shape == (2, 2)
        # q_st/r_st remain tile centers, matching the legacy geom_tile plots.
        assert _collection_center(axes[0, 0]) == pytest.approx((70, 20))
        assert _collection_center(axes[1, 1]) == pytest.approx((20, 70))
    finally:
        plt.close(figure)


def test_self_comparisons_run_bottom_left_to_top_right():
    names = ["sequence_a", "sequence_b", "sequence_c"]
    kwargs = _grid_kwargs(
        singles=[
            _records(name, name, [(10 + index * 10, 10 + index * 10, 100)])
            for index, name in enumerate(names)
        ],
        doubles=[
            _records(names[left], names[right], [(10, 20, 95)])
            for left in range(len(names))
            for right in range(left + 1, len(names))
        ],
        single_names=names,
        double_names=[
            [names[left], names[right]]
            for left in range(len(names))
            for right in range(left + 1, len(names))
        ],
    )

    figure, axes = _build_grid_figure(**kwargs)
    try:
        assert _collection_center(axes[2, 0]) == pytest.approx((10, 10))
        assert _collection_center(axes[1, 1]) == pytest.approx((20, 20))
        assert _collection_center(axes[0, 2]) == pytest.approx((30, 30))
    finally:
        plt.close(figure)


def test_self_comparison_grid_diagonal_uses_full_symmetric_dotplots():
    kwargs = _basic_two_sequence_grid(
        singles=[
            _records(
                "sequence_a",
                "sequence_a",
                [(10, 30, 91), (50, 50, 100)],
            ),
            _records(
                "sequence_b",
                "sequence_b",
                [(20, 40, 92), (60, 60, 100)],
            ),
        ]
    )

    figure, axes = _build_grid_figure(**kwargs)
    try:
        assert _collection_centers(axes[1, 0]) == [
            (10, 30),
            (30, 10),
            (50, 50),
        ]
        assert _collection_centers(axes[0, 1]) == [
            (20, 40),
            (40, 20),
            (60, 60),
        ]
        assert len(_collection_centers(axes[0, 0])) == 1
        assert len(_collection_centers(axes[1, 1])) == 1
    finally:
        plt.close(figure)


def test_full_self_comparison_input_is_not_mirrored_twice():
    kwargs = _basic_two_sequence_grid(
        singles=[
            _records("sequence_a", "sequence_a", [(10, 30, 91), (30, 10, 91)]),
            _records("sequence_b", "sequence_b", [(20, 40, 92), (40, 20, 92)]),
        ]
    )

    figure, axes = _build_grid_figure(**kwargs)
    try:
        assert _collection_centers(axes[1, 0]) == [(10, 30), (30, 10)]
        assert _collection_centers(axes[0, 1]) == [(20, 40), (40, 20)]
    finally:
        plt.close(figure)


def test_grid_region_names_fit_panels_and_numeric_labels_are_doubled(monkeypatch):
    monkeypatch.setattr(static_plots, "DEFAULT_FONT_FAMILY", FALLBACK_FONT_FAMILY)
    names = [
        "PAN010.chr14.haplotype1.paternal:1-4000000",
        "PAN010.chr14.haplotype2.maternal:1-4000000",
        "PAN027.chr14.paternal:1-4000000",
    ]
    kwargs = _grid_kwargs(
        singles=[_records(name, name, [(10, 10, 100)]) for name in names],
        doubles=[
            _records(names[left], names[right], [(10, 20, 95)])
            for left in range(len(names))
            for right in range(left + 1, len(names))
        ],
        single_names=names,
        double_names=[
            [names[left], names[right]]
            for left in range(len(names))
            for right in range(left + 1, len(names))
        ],
    )
    kwargs["width"] = 6

    figure, axes = _build_grid_figure(**kwargs)
    try:
        renderer = figure.canvas.get_renderer()
        assert [axis.title.get_text() for axis in axes[0]] == [
            "PAN010.chr14.haplotype1.paternal",
            "PAN010.chr14.haplotype2.maternal",
            "PAN027.chr14.paternal",
        ]
        assert all(
            ":1-4000000" not in axis.yaxis.label.get_text() for axis in axes[:, 0]
        )

        for axis in axes[0]:
            title_width = axis.title.get_window_extent(renderer).width
            assert title_width <= axis.get_window_extent(renderer).width * 0.91 + 1
        for axis in axes[:, 0]:
            label_height = axis.yaxis.label.get_window_extent(renderer).height
            assert label_height <= axis.get_window_extent(renderer).height * 0.91 + 1

        assert axes[-1, 0].get_xticklabels()[0].get_fontsize() == pytest.approx(8)
    finally:
        plt.close(figure)


def test_grid_width_controls_total_figure_size():
    kwargs = _basic_two_sequence_grid(width=9)

    figure, _axes = _build_grid_figure(**kwargs)
    try:
        assert tuple(figure.get_size_inches()) == pytest.approx((9, 9))
    finally:
        plt.close(figure)


def test_grid_without_explicit_limits_uses_only_selected_region_bounds():
    names = ["sequence_a:1000001-2000000", "sequence_b:1000001-2000000"]
    kwargs = _grid_kwargs(
        singles=[
            _records(names[0], names[0], [(1_000_001, 1_999_990, 91)]),
            _records(names[1], names[1], [(1_000_001, 1_999_990, 92)]),
        ],
        doubles=[_records(names[0], names[1], [(1_100_000, 1_900_000, 95)])],
        single_names=names,
        double_names=[[names[0], names[1]]],
    )
    kwargs["xlim"] = None
    kwargs["axes_label"] = None
    kwargs["breaks"] = None

    figure, axes = _build_grid_figure(**kwargs)
    try:
        assert axes[0, 0].get_xlim() == pytest.approx((1_000_001, 2_000_000))
        assert axes[1, 1].get_ylim() == pytest.approx((1_000_001, 2_000_000))
    finally:
        plt.close(figure)


def test_grid_exact_bounds_do_not_expand_to_next_nice_tick():
    kwargs = _basic_two_sequence_grid(
        xlim=(1, 103_156_783),
        axes_label=None,
        breaks=None,
    )

    figure, axes = _build_grid_figure(**kwargs)
    try:
        assert axes[0, 0].get_xlim() == pytest.approx((1, 103_156_783))
        assert axes[1, 1].get_ylim() == pytest.approx((1, 103_156_783))
        assert max(axes[0, 0].get_xticks()) <= 103_156_783
    finally:
        plt.close(figure)


@pytest.mark.parametrize(
    ("axis_end", "unit"),
    [(100_000, "Kbp"), (103_156_783, "Mbp"), (500_000_000, "Gbp")],
)
def test_grid_labels_genomic_units_only_on_bottom_left_cell(
    axis_end, unit, monkeypatch
):
    monkeypatch.setattr(static_plots, "DEFAULT_FONT_FAMILY", FALLBACK_FONT_FAMILY)
    kwargs = _basic_two_sequence_grid(
        xlim=(1, axis_end),
        axes_label=None,
        breaks=None,
    )

    figure, axes = _build_grid_figure(**kwargs)
    try:
        expected = f"Genomic Position ({unit})"
        bottom_left_axis = axes[-1, 0]
        assert bottom_left_axis.get_xlabel() == expected
        assert sum(axis.get_xlabel() == expected for axis in axes.flat) == 1
        assert getattr(figure, "_supxlabel", None) is None
        assert getattr(figure, "_supylabel", None) is None
        vertical_titles = [
            text
            for axis in axes.flat
            for text in axis.texts
            if text.get_gid() == "grid-y-axis-title"
        ]
        assert len(vertical_titles) == 1
        assert vertical_titles[0].get_text() == expected
        assert list(axis.get_ylabel() for axis in axes[:, 0]) == [
            "sequence_b",
            "sequence_a",
        ]

        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        figure_bounds = figure.bbox
        cell_bounds = bottom_left_axis.get_window_extent(renderer)
        horizontal_bounds = bottom_left_axis.xaxis.label.get_window_extent(renderer)
        vertical_bounds = vertical_titles[0].get_window_extent(renderer)
        row_label_bounds = bottom_left_axis.yaxis.label.get_window_extent(renderer)

        assert cell_bounds.x0 <= horizontal_bounds.x0 + horizontal_bounds.width / 2
        assert horizontal_bounds.x0 + horizontal_bounds.width / 2 <= cell_bounds.x1
        assert cell_bounds.y0 <= vertical_bounds.y0 + vertical_bounds.height / 2
        assert vertical_bounds.y0 + vertical_bounds.height / 2 <= cell_bounds.y1
        assert vertical_bounds.x1 + 1 <= row_label_bounds.x0
        for bounds in (horizontal_bounds, vertical_bounds):
            assert figure_bounds.contains(bounds.x0, bounds.y0)
            assert figure_bounds.contains(bounds.x1, bounds.y1)
    finally:
        plt.close(figure)


@pytest.mark.parametrize(("width", "expected_width"), [(1, 2), (4, 4)])
def test_grid_axis_labels_do_not_expand_saved_canvas(tmp_path, width, expected_width):
    kwargs = _basic_two_sequence_grid(width=width)

    with plt.rc_context({"savefig.bbox": "tight"}):
        _create_grid(tmp_path, **kwargs)

    image = plt.imread(tmp_path / "2x2_GRID.png")
    assert image.shape[:2] == (expected_width * 72, expected_width * 72)


def test_one_cell_grid_axis_titles_stay_inside_canvas():
    name = "sequence_a"
    kwargs = _grid_kwargs(
        singles=[_records(name, name, [(10, 10, 100)])],
        doubles=[],
        single_names=[name],
        double_names=[],
    )
    kwargs["width"] = 4

    figure, axes = _build_grid_figure(**kwargs)
    try:
        vertical_title = next(
            text for text in axes[0, 0].texts if text.get_gid() == "grid-y-axis-title"
        )
        for family in (None, FALLBACK_FONT_FAMILY):
            if family is not None:
                set_figure_font_family(figure, family)
            for dpi in (72, 100, 300, 600):
                figure.set_dpi(dpi)
                figure.canvas.draw()
                renderer = figure.canvas.get_renderer()
                for title in (axes[0, 0].xaxis.label, vertical_title):
                    bounds = title.get_window_extent(renderer)
                    assert figure.bbox.contains(bounds.x0, bounds.y0)
                    assert figure.bbox.contains(bounds.x1, bounds.y1)
    finally:
        plt.close(figure)


def test_compare_only_grid_derives_sequence_names_from_pair_metadata(tmp_path):
    double_names = [
        ["gamma", "alpha"],
        ["beta", "gamma"],
        ["alpha", "beta"],
    ]
    doubles = [
        _records(query, reference, [(10 + index * 10, 60, 90 + index)])
        for index, (query, reference) in enumerate(double_names)
    ]
    kwargs = _grid_kwargs(
        singles=[],
        doubles=doubles,
        single_names=[],
        double_names=double_names,
    )

    _create_grid(tmp_path, **kwargs)

    assert (tmp_path / "3x3_GRID.png").is_file()
    svg = (tmp_path / "3x3_GRID.svg").read_text(encoding="utf-8")
    for name in ("alpha", "beta", "gamma"):
        assert name in svg


def test_native_grid_supports_more_than_six_sequences():
    names = [f"sequence_{index}" for index in range(7)]
    double_names = [
        [names[left], names[right]]
        for left in range(len(names))
        for right in range(left + 1, len(names))
    ]
    doubles = [
        _records(query, reference, [(10, 20, 95)]) for query, reference in double_names
    ]
    kwargs = _grid_kwargs(
        singles=[],
        doubles=doubles,
        single_names=[],
        double_names=double_names,
    )

    figure, axes = _build_grid_figure(**kwargs)
    try:
        assert axes.shape == (7, 7)
        assert sum(bool(axis.collections) for axis in axes.flat) == 42
    finally:
        plt.close(figure)


@pytest.mark.parametrize(
    ("singles", "doubles", "single_names", "double_names", "message"),
    [
        (
            [
                _records("a", "a", [(0, 0, 100)]),
                _records("b", "b", [(0, 0, 100)]),
                _records("c", "c", [(0, 0, 100)]),
            ],
            [
                _records("a", "b", [(0, 0, 90)]),
                _records("a", "c", [(0, 0, 90)]),
            ],
            ["a", "b", "c"],
            [["a", "b"], ["a", "c"]],
            "missing",
        ),
        (
            [
                _records("a", "a", [(0, 0, 100)]),
                _records("b", "b", [(0, 0, 100)]),
            ],
            [
                _records("a", "b", [(0, 0, 90)]),
                _records("b", "a", [(0, 0, 91)]),
            ],
            ["a", "b"],
            [["a", "b"], ["b", "a"]],
            "duplicate",
        ),
    ],
)
def test_grid_pair_metadata_errors_are_actionable(
    singles, doubles, single_names, double_names, message
):
    kwargs = _grid_kwargs(
        singles=singles,
        doubles=doubles,
        single_names=single_names,
        double_names=double_names,
    )

    with pytest.raises(ValueError, match=f"(?i){message}"):
        _build_grid_figure(**kwargs)


def test_grid_honors_custom_breakpoints_and_colors():
    custom_colors = ["#010203", "#456789", "#abcdef"]
    pair = _records(
        "sequence_a",
        "sequence_b",
        [(10, 10, 10), (30, 30, 50), (50, 50, 90)],
    )
    kwargs = _basic_two_sequence_grid(
        doubles=[pair],
        custom_colors=custom_colors,
        custom_breakpoints=[0, 34, 67, 101],
    )

    figure, axes = _build_grid_figure(**kwargs)
    try:
        figure.canvas.draw()
        rendered_colors = _all_artist_colors(axes)
        for color in custom_colors:
            assert to_rgba(color) in rendered_colors
    finally:
        plt.close(figure)


def test_grid_uses_blue_and_pink_direction_colors_in_every_cell():
    kwargs = _basic_two_sequence_grid()
    direction_header = (*BED_HEADER, "direction")

    def add_directions(records):
        return [direction_header] + [
            (*row, "Forward" if index % 2 == 0 else "Reverse")
            for index, row in enumerate(records[1:])
        ]

    kwargs["singles"] = [add_directions(records) for records in kwargs["singles"]]
    pair = _records(
        "sequence_a",
        "sequence_b",
        [(10, 10, 86), (30, 40, 100), (50, 60, 86), (70, 80, 100)],
    )
    kwargs["doubles"] = [
        [direction_header]
        + [
            (*row, direction)
            for row, direction in zip(
                pair[1:], ["Forward", "Forward", "Reverse", "Reverse"]
            )
        ]
    ]

    figure, axes = _build_grid_figure(**kwargs)
    try:
        rendered_colors = _all_artist_colors(axes)
        assert to_rgba(DIRECTION_COLORS["Forward"]) in rendered_colors
        assert to_rgba(DIRECTION_COLORS["Reverse"]) in rendered_colors
        assert len(rendered_colors) >= 4
    finally:
        plt.close(figure)


@pytest.mark.parametrize(
    ("deraster", "expected_rasterized"), [(False, True), (True, False)]
)
def test_grid_deraster_controls_collection_rasterization(deraster, expected_rasterized):
    figure, axes = _build_grid_figure(**_basic_two_sequence_grid(deraster=deraster))
    try:
        plotted_collections = [
            collection
            for axis in axes.flat
            for collection in axis.collections
            if collection.get_paths()
        ]
        assert plotted_collections
        assert {collection.get_rasterized() for collection in plotted_collections} == {
            expected_rasterized
        }
    finally:
        plt.close(figure)
