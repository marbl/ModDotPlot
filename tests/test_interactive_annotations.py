import numpy as np
import plotly.graph_objs as go
import dash

from moddotplot.annotations import read_annotation_beds
from moddotplot.interactive import (
    add_annotation_tracks,
    interactive_axis_bounds,
    preserve_zoom_ranges,
    run_dash,
)


def _figure():
    figure = go.Figure(data=[go.Heatmap(z=[[100.0]])])
    figure.update_xaxes(title_text="chrX")
    figure.update_yaxes(title_text="chrY")
    return figure


def test_multiple_beds_add_tracks_to_both_comparative_axes(tmp_path):
    x_bed = tmp_path / "x.bed"
    y_bed = tmp_path / "y.bed"
    x_bed.write_text("chrX\t10\t30\tx-feature\t0\t+\t10\t30\t255,0,0\n")
    y_bed.write_text("chrY\t40\t80\ty-feature\t0\t+\t40\t80\t0,0,255\n")
    annotations = read_annotation_beds([x_bed, y_bed])
    metadata = {
        "x_name": "chrX",
        "y_name": "chrY",
        "x_size": 100,
        "y_size": 120,
        "self": False,
    }

    figure = add_annotation_tracks(_figure(), metadata, annotations)

    x_shapes = [
        shape
        for shape in figure.layout.shapes
        if shape.xref == "x" and shape.yref == "paper"
    ]
    y_shapes = [
        shape
        for shape in figure.layout.shapes
        if shape.xref == "paper" and shape.yref == "y"
    ]
    assert len(x_shapes) == 2  # track background plus one interval
    assert len(y_shapes) == 2
    assert (x_shapes[1].x0, x_shapes[1].x1) == (10, 30)
    assert x_shapes[1].fillcolor == "rgb(255,0,0)"
    assert (y_shapes[1].y0, y_shapes[1].y1) == (40, 80)
    assert y_shapes[1].fillcolor == "rgb(0,0,255)"
    assert figure.layout.margin.b >= 140
    assert figure.layout.margin.l >= 140


def test_self_plot_adds_only_x_track_and_clips_to_sequence(tmp_path):
    bed = tmp_path / "annotations.bed"
    bed.write_text("chrX\t90\t130\nchrY\t10\t20\n")
    annotations = read_annotation_beds([bed])
    metadata = {
        "x_name": "chrX",
        "y_name": "chrX",
        "x_size": 100,
        "y_size": 100,
        "self": True,
    }

    figure = add_annotation_tracks(_figure(), metadata, annotations)

    shapes = list(figure.layout.shapes)
    assert len(shapes) == 2
    assert all(shape.xref == "x" and shape.yref == "paper" for shape in shapes)
    assert (shapes[1].x0, shapes[1].x1) == (90, 100)


def test_region_suffixed_header_matches_base_bed_chromosome(tmp_path):
    bed = tmp_path / "annotations.bed"
    bed.write_text("chr14_MATERNAL\t14002000\t14004000\n")
    annotations = read_annotation_beds([bed])
    metadata = {
        "x_name": "chr14_MATERNAL:14000001-18000000",
        "y_name": "chr14_MATERNAL:14000001-18000000",
        "x_size": 4_000_000,
        "y_size": 4_000_000,
        "self": True,
    }

    figure = add_annotation_tracks(_figure(), metadata, annotations)

    shapes = list(figure.layout.shapes)
    assert interactive_axis_bounds(metadata, "x") == (14_000_001, 18_000_000)
    assert len(shapes) == 2
    assert (shapes[0].x0, shapes[0].x1) == (14_000_001, 18_000_000)
    assert (shapes[1].x0, shapes[1].x1) == (14_002_000, 14_004_000)


def test_nonmatching_bed_does_not_add_tracks(tmp_path):
    bed = tmp_path / "annotations.bed"
    bed.write_text("other\t10\t20\n")
    annotations = read_annotation_beds([bed])
    metadata = {
        "x_name": "chrX",
        "y_name": "chrY",
        "x_size": 100,
        "y_size": 100,
        "self": False,
    }

    figure = add_annotation_tracks(_figure(), metadata, annotations)

    assert not figure.layout.shapes


def test_bed_shapes_do_not_override_explicit_zoom_range(tmp_path):
    bed = tmp_path / "annotations.bed"
    bed.write_text("chrX\t0\t1000000\n")
    annotations = read_annotation_beds([bed])
    metadata = {
        "x_name": "chrX",
        "y_name": "chrX",
        "x_size": 1_000_000,
        "y_size": 1_000_000,
        "self": True,
    }
    figure = add_annotation_tracks(_figure(), metadata, annotations)

    preserve_zoom_ranges(
        figure,
        x_range=(200_000, 300_000),
        y_range=(400_000, 500_000),
    )

    assert tuple(figure.layout.xaxis.range) == (200_000, 300_000)
    assert tuple(figure.layout.yaxis.range) == (400_000, 500_000)
    assert figure.layout.xaxis.autorange is False
    assert figure.layout.yaxis.autorange is False
    assert figure.layout.shapes[0].x0 == 0
    assert figure.layout.shapes[0].x1 == 1_000_000


def _find_component(component, component_id):
    if getattr(component, "id", None) == component_id:
        return component
    children = getattr(component, "children", None)
    if children is None:
        return None
    if not isinstance(children, (list, tuple)):
        children = [children]
    for child in children:
        match = _find_component(child, component_id)
        if match is not None:
            return match
    return None


def test_dash_initial_comparative_figure_contains_both_tracks(monkeypatch, tmp_path):
    bed = tmp_path / "annotations.bed"
    bed.write_text("chrX\t10\t30\nchrY\t40\t80\n")
    annotations = read_annotation_beds([bed])
    metadata = [
        {
            "x_name": "chrX",
            "y_name": "chrY",
            "x_size": 100,
            "y_size": 120,
            "self": False,
            "min_window_size": 50,
            "max_window_size": 50,
            "resolution": 2,
            "title": "chrX-chrY",
            "sparsities": [1],
        }
    ]
    apps = []
    monkeypatch.setattr(dash.Dash, "run", lambda app, **_kwargs: apps.append(app))

    run_dash(
        [[np.array([[100.0, 90.0], [90.0, 100.0]])]],
        metadata,
        [[[0, 50, 100], [0, 60, 120]]],
        sparsity=1,
        identity=86.0,
        port_number=8050,
        output_dir=None,
        annotations=annotations,
    )

    graph = _find_component(apps[0].layout, "dotplot")
    shapes = list(graph.figure.layout.shapes)
    assert any(shape.xref == "x" and shape.yref == "paper" for shape in shapes)
    assert any(shape.xref == "paper" and shape.yref == "y" for shape in shapes)
