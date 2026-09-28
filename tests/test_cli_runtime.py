import sys

import numpy as np
import pytest

import moddotplot.moddotplot as cli


def _patch_fasta_input(monkeypatch, names, kmers):
    monkeypatch.setattr(cli, "isValidFasta", lambda _path: True)
    monkeypatch.setattr(cli, "getInputHeaders", lambda _path: names)
    monkeypatch.setattr(cli, "readKmersFromFile", lambda *_args: kmers)


def _patch_static_calculation(monkeypatch, plot_calls, pair_calls=None):
    monkeypatch.setattr(cli, "createSelfMatrix", lambda *_args: np.full((1, 1), 100.0))

    def create_pairwise(*args):
        if pair_calls is not None:
            pair_calls.append(args)
        return np.full((1, 1), 95.0)

    monkeypatch.setattr(cli, "createPairwiseMatrix", create_pairwise)
    monkeypatch.setattr(
        cli,
        "convertMatrixToBed",
        lambda *_args, **_kwargs: [["header"], ["value"]],
    )
    monkeypatch.setattr(cli, "create_plots", lambda **kwargs: plot_calls.append(kwargs))


def test_main_without_subcommand_reports_parser_error(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["moddotplot"])

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    captured = capsys.readouterr()
    assert exc_info.value.code == 2
    assert "the following arguments are required: command" in captured.err
    assert "{interactive,static}" in captured.err


@pytest.mark.parametrize("option", ["--colors", "--color"])
def test_static_parser_accepts_color_option_aliases(option):
    args = cli.get_parser().parse_args(
        ["static", "--fasta", "sequence.fa", option, "#010203", "#abcdef"]
    )

    assert args.colors == ["#010203", "#abcdef"]


def test_static_config_prefers_canonical_colors_key():
    args = cli.get_parser().parse_args(["static", "--fasta", "sequence.fa"])

    cli._apply_static_config(
        args,
        {
            "fasta": ["sequence.fa"],
            "colors": ["#canonical"],
            "color": ["#legacy"],
        },
    )

    assert args.colors == ["#canonical"]


def test_static_config_supports_legacy_color_key():
    args = cli.get_parser().parse_args(["static", "--fasta", "sequence.fa"])

    cli._apply_static_config(args, {"fasta": ["sequence.fa"], "color": ["#legacy"]})

    assert args.colors == ["#legacy"]


def test_static_delta_defaults_to_half_window():
    args = cli.get_parser().parse_args(["static", "--fasta", "sequence.fa"])

    assert args.delta == 0.5


def test_region_kmer_slice_uses_one_based_inclusive_base_coordinates():
    kmers = np.arange(980)  # 980 21-mers represent a 1,000-base sequence.

    selected = cli._slice_kmers_for_region(kmers, ("chrA", 101, 400), 21)

    assert len(selected) == 280
    np.testing.assert_array_equal(selected, np.arange(100, 380))


@pytest.mark.parametrize(
    "region",
    [("chrA", 0, 100), ("chrA", 900, 1_001), ("chrA", 10, 20)],
)
def test_region_kmer_slice_rejects_invalid_bounds_or_short_intervals(region):
    with pytest.raises(ValueError):
        cli._slice_kmers_for_region(np.arange(980), region, 21)


def test_interactive_window_uses_longest_sequence_by_length(monkeypatch):
    # The shorter k-mer list is lexicographically greater. This reproduces the
    # old len(max(k_list)) bug while keeping the test computation tiny.
    short_kmers = [9] * 100
    long_kmers = [1] * 1000
    _patch_fasta_input(monkeypatch, ["short", "long"], [short_kmers, long_kmers])
    monkeypatch.setattr(cli, "partitionOverlaps", lambda *_args: [])
    monkeypatch.setattr(cli, "convertToModimizers", lambda *_args: [])
    monkeypatch.setattr(cli, "selfContainmentMatrix", lambda *_args: np.zeros((1, 1)))
    dash_calls = []
    monkeypatch.setattr(cli, "run_dash", lambda *args: dash_calls.append(args))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "moddotplot",
            "interactive",
            "--fasta",
            "sequence.fa",
            "--resolution",
            "10",
            "--quick",
        ],
    )

    cli.main()

    metadata = dash_calls[0][1]
    assert {entry["min_window_size"] for entry in metadata} == {102}
    assert {entry["max_window_size"] for entry in metadata} == {102}


def test_no_bedpe_self_plot_uses_sequence_output_directory(monkeypatch, tmp_path):
    _patch_fasta_input(monkeypatch, ["chrA"], [[1] * 1000])
    plot_calls = []
    _patch_static_calculation(monkeypatch, plot_calls)
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
            "--output-dir",
            str(tmp_path),
        ],
    )

    cli.main()

    expected_directory = tmp_path / "chrA"
    assert expected_directory.is_dir()
    assert plot_calls[0]["directory"] == str(expected_directory)
    assert not list(tmp_path.rglob("*.bedpe"))


def test_unmatched_region_fails_instead_of_silently_using_full_sequences(
    monkeypatch, tmp_path, capsys
):
    larger = [1] * 1000
    smaller = [2] * 800
    _patch_fasta_input(monkeypatch, ["chrA", "chrB"], [larger, smaller])
    plot_calls = []
    pair_calls = []
    _patch_static_calculation(monkeypatch, plot_calls, pair_calls)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "moddotplot",
            "static",
            "--fasta",
            "sequences.fa",
            "--compare-only",
            "--region",
            "missing:1-100",
            "--resolution",
            "10",
            "--no-bedpe",
            "--output-dir",
            str(tmp_path),
        ],
    )

    with pytest.raises(SystemExit) as error:
        cli.main()

    captured = capsys.readouterr()
    assert error.value.code == 2
    assert "does not match any FASTA identifier" in captured.out
    assert not pair_calls
    assert not plot_calls
    assert not list(tmp_path.rglob("*.bedpe"))


def test_compare_only_region_can_target_just_one_sequence(monkeypatch, tmp_path):
    larger = [1] * 1000
    smaller = [2] * 800
    _patch_fasta_input(monkeypatch, ["chrA", "chrB"], [larger, smaller])
    plot_calls = []
    pair_calls = []
    _patch_static_calculation(monkeypatch, plot_calls, pair_calls)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "moddotplot",
            "static",
            "--fasta",
            "sequences.fa",
            "--compare-only",
            "--region",
            "chrB:101-400",
            "--resolution",
            "10",
            "--no-bedpe",
            "--output-dir",
            str(tmp_path),
        ],
    )

    cli.main()

    pair_args = pair_calls[0]
    assert pair_args[0] == 280
    assert len(pair_args[2]) == 280
    assert pair_args[3] is larger
    assert plot_calls[0]["name_x"] == "chrA"
    assert plot_calls[0]["name_y"] == "chrB:101-400"
    assert plot_calls[0]["xlim"] == (1, 1020)
    assert not list(tmp_path.rglob("*.bedpe"))
