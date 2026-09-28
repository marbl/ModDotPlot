import sys

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
    monkeypatch.setattr(
        cli,
        "convertMatrixToBed",
        lambda *_args, **_kwargs: [["header"], ["value"]],
    )
    monkeypatch.setattr(cli, "create_plots", lambda **_kwargs: None)
    monkeypatch.setattr(cli, "create_grid", lambda **_kwargs: None)
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


def test_static_direction_plot_receives_canonical_and_forward_self_matrices(
    monkeypatch, tmp_path
):
    canonical_hashes = [[1] * 1000]
    forward_hashes = [[2] * 1000]
    _patch_common_static_io(
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
    direction_calls = []
    monkeypatch.setattr(
        cli,
        "create_direction_plot",
        lambda **kwargs: direction_calls.append(kwargs),
    )

    cli.main()

    assert len(direction_calls) == 1
    call = direction_calls[0]
    assert call["canonical_matrix"] is canonical_matrix
    assert call["forward_matrix"] is forward_matrix
    assert call["self_identity"] is True
    assert call["name_x"] == call["name_y"] == "chrA"


def test_static_direction_plot_receives_pairwise_matrices(monkeypatch, tmp_path):
    canonical_hashes = [[1] * 1000, [2] * 800]
    forward_hashes = [[3] * 1000, [4] * 800]
    _patch_common_static_io(
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
    direction_calls = []
    monkeypatch.setattr(
        cli,
        "create_direction_plot",
        lambda **kwargs: direction_calls.append(kwargs),
    )

    cli.main()

    assert len(direction_calls) == 1
    call = direction_calls[0]
    assert call["canonical_matrix"] is canonical_matrix
    assert call["forward_matrix"] is forward_matrix
    assert call["self_identity"] is False
    assert (call["name_x"], call["name_y"]) == ("chrA", "chrB")
