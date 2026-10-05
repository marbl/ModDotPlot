import sys
from pathlib import Path

import numpy as np

import moddotplot.moddotplot as cli


def _patch_common_static_io(monkeypatch, names, hash_sets, tmp_path):
    monkeypatch.setattr(cli, "isValidFasta", lambda _path: True)
    monkeypatch.setattr(cli, "getInputHeaders", lambda _path: names)

    def read_hashes(
        _path,
        _kmer,
        _quiet,
        forward_only,
        ambiguous=False,
        regions=None,
        record_ids=None,
    ):
        assert ambiguous is False
        assert not regions
        assert record_ids == names
        return hash_sets[forward_only]

    monkeypatch.setattr(cli, "readKmersFromFile", read_hashes)

    def convert_matrix_to_bed(
        matrix,
        window_size,
        _identity,
        name_x,
        name_y,
        self_identity,
        x_offset,
        y_offset,
        *_args,
    ):
        rows = [
            (
                "#query_name",
                "query_start",
                "query_end",
                "reference_name",
                "reference_start",
                "reference_end",
                "perID_by_events",
            )
        ]
        for query_index, reference_index in np.argwhere(matrix > 0):
            if self_identity and query_index > reference_index:
                continue
            query_start = query_index * window_size + x_offset
            reference_start = reference_index * window_size + y_offset
            rows.append(
                (
                    name_x,
                    query_start,
                    query_start + window_size - 1,
                    name_y,
                    reference_start,
                    reference_start + window_size - 1,
                    matrix[query_index, reference_index],
                )
            )
        return rows

    monkeypatch.setattr(cli, "convertMatrixToBed", convert_matrix_to_bed)
    plot_calls = []
    grid_calls = []

    def capture_plot(**kwargs):
        plot_calls.append(kwargs)
        return [str(Path(kwargs["directory"]) / "mock_plot.png")]

    def capture_grid(**kwargs):
        grid_calls.append(kwargs)
        return [str(Path(kwargs["directory"]) / "mock_grid.png")]

    monkeypatch.setattr(cli, "create_plots", capture_plot)
    monkeypatch.setattr(cli, "create_grid", capture_grid)
    monkeypatch.setattr(cli, "read_df_from_file", lambda _path: None)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "moddotplot",
            "static",
            "--fasta",
            "sequence.fa",
            "--resolution",
            "10",
            "--no-bedpe",
            "--plot-direction",
            "--output-dir",
            str(tmp_path),
        ],
    )
    return plot_calls, grid_calls


def test_static_direction_colors_standard_self_plot_rows(monkeypatch, tmp_path):
    canonical_hashes = [[1] * 1000]
    forward_hashes = [[2] * 1000]
    plot_calls, _grid_calls = _patch_common_static_io(
        monkeypatch,
        ["chrA"],
        {False: canonical_hashes, True: forward_hashes},
        tmp_path,
    )
    canonical_matrix = np.array([[100.0, 90.0], [90.0, 100.0]])
    forward_matrix = np.array([[100.0, 0.0], [0.0, 100.0]])

    def create_self(_length, sequence, *_args):
        return canonical_matrix if sequence is canonical_hashes[0] else forward_matrix

    monkeypatch.setattr(cli, "createSelfMatrix", create_self)
    cli.main()

    assert len(plot_calls) == 2
    assert "direction" not in plot_calls[0]["sdf"][0][0]
    assert plot_calls[1]["directory"].endswith("directionality")
    direction_bed = plot_calls[1]["sdf"][0]
    assert direction_bed[0][-1] == "direction"
    assert [row[-1] for row in direction_bed[1:]] == [
        "Forward",
        "Reverse",
        "Forward",
    ]
    assert plot_calls[1]["name_x"] == plot_calls[1]["name_y"] == "chrA"
    direction_summary = tmp_path / "chrA" / "directionality" / "plot_summary.txt"
    assert direction_summary.is_file()
    assert "Plot group" not in direction_summary.read_text()


def test_direction_coloring_with_forward_option_keeps_reverse_hits(
    monkeypatch, tmp_path
):
    canonical_hashes = [[1] * 1000]
    forward_hashes = [[2] * 1000]
    plot_calls, _grid_calls = _patch_common_static_io(
        monkeypatch,
        ["chrA"],
        {False: canonical_hashes, True: forward_hashes},
        tmp_path,
    )
    sys.argv.append("--forward")
    canonical_matrix = np.array([[100.0, 90.0], [90.0, 100.0]])
    forward_matrix = np.array([[100.0, 0.0], [0.0, 100.0]])

    def create_self(_length, sequence, *_args):
        return canonical_matrix if sequence is canonical_hashes[0] else forward_matrix

    monkeypatch.setattr(cli, "createSelfMatrix", create_self)
    cli.main()

    assert len(plot_calls) == 2
    direction_bed = plot_calls[1]["sdf"][0]
    assert [row[-1] for row in direction_bed[1:]] == [
        "Forward",
        "Reverse",
        "Forward",
    ]


def test_static_direction_colors_standard_pairwise_plot_rows(monkeypatch, tmp_path):
    canonical_hashes = [[1] * 1000, [2] * 800]
    forward_hashes = [[3] * 1000, [4] * 800]
    plot_calls, _grid_calls = _patch_common_static_io(
        monkeypatch,
        ["chrA", "chrB"],
        {False: canonical_hashes, True: forward_hashes},
        tmp_path,
    )
    sys.argv.extend(["--compare-only"])
    canonical_matrix = np.array([[90.0]])
    forward_matrix = np.array([[0.0]])

    def create_pair(_y_length, _x_length, y_sequence, _x_sequence, *_args):
        return canonical_matrix if y_sequence is canonical_hashes[1] else forward_matrix

    monkeypatch.setattr(cli, "createPairwiseMatrix", create_pair)
    cli.main()

    assert len(plot_calls) == 2
    assert "direction" not in plot_calls[0]["sdf"][0][0]
    assert plot_calls[1]["directory"].endswith("directionality")
    direction_bed = plot_calls[1]["sdf"][0]
    assert direction_bed[0][-1] == "direction"
    assert [row[-1] for row in direction_bed[1:]] == ["Reverse"]
    assert (plot_calls[1]["name_x"], plot_calls[1]["name_y"]) == ("chrA", "chrB")


def test_grid_only_receives_direction_colored_self_and_pairwise_rows(
    monkeypatch, tmp_path
):
    canonical_hashes = [[1] * 1000, [2] * 800]
    forward_hashes = [[3] * 1000, [4] * 800]
    _plot_calls, grid_calls = _patch_common_static_io(
        monkeypatch,
        ["chrA", "chrB"],
        {False: canonical_hashes, True: forward_hashes},
        tmp_path,
    )
    sys.argv.extend(["--grid-only"])
    monkeypatch.setattr(cli, "ModimizerSketchCache", lambda **_kwargs: None)

    def create_self(_length, sequence, *_args):
        if sequence in canonical_hashes:
            return np.array([[100.0, 95.0], [95.0, 100.0]])
        return np.array([[100.0, 0.0], [0.0, 100.0]])

    def create_pair(_y_length, _x_length, y_sequence, _x_sequence, *_args):
        if y_sequence in canonical_hashes:
            return np.array([[95.0]])
        return np.array([[0.0]])

    monkeypatch.setattr(cli, "createSelfMatrix", create_self)
    monkeypatch.setattr(cli, "createPairwiseMatrix", create_pair)

    cli.main()

    assert len(grid_calls) == 2
    assert all("direction" not in matrix[0] for matrix in grid_calls[0]["singles"])
    assert all("direction" not in matrix[0] for matrix in grid_calls[0]["doubles"])
    assert grid_calls[1]["directory"].endswith("directionality")
    grid_call = grid_calls[1]
    assert all(matrix[0][-1] == "direction" for matrix in grid_call["singles"])
    assert all(matrix[0][-1] == "direction" for matrix in grid_call["doubles"])
    assert any(
        row[-1] == "Reverse"
        for matrix in [*grid_call["singles"], *grid_call["doubles"]]
        for row in matrix[1:]
    )
    assert (tmp_path / "directionality" / "plot_summary.txt").is_file()
