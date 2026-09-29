import sys
import shlex
from types import SimpleNamespace

import numpy as np
import pytest

import moddotplot.moddotplot as cli


def _matrix_args(**overrides):
    values = {
        "window": None,
        "resolution": 1000,
        "kmer": 21,
        "modimizer": 1000,
    }
    values.update(overrides)
    return SimpleNamespace(**values)


def test_matrix_config_caps_resolution_at_one_valid_kmer_per_window():
    args = _matrix_args()

    config = cli._matrix_config_for_length(16_549, args)

    assert config.window_size == 21
    assert config.resolution == 789
    assert config.modimizer == 21
    assert config.sparsity == 1
    assert config.expectation == 21
    assert args.modimizer == 1000


def test_matrix_config_preserves_default_nuclear_chromosome_parameters():
    config = cli._matrix_config_for_length(248_387_308, _matrix_args())

    assert config.window_size == 248_388
    assert config.resolution == 1000
    assert config.modimizer == 1000
    assert config.sparsity == 128
    assert config.expectation == 1941


def test_streaming_self_runner_finishes_one_record_before_requesting_next(
    monkeypatch,
):
    events = []

    def records(*_args, **_kwargs):
        events.append("yield:first")
        yield "first", "ACGT", "first"
        events.append("yield:second")
        yield "second", "TGCA", "second"

    def process(**kwargs):
        events.append(f"process:{kwargs['sequence_id']}")

    monkeypatch.setattr(cli, "iter_fasta_records", records)
    monkeypatch.setattr(cli, "_process_static_self_record", process)
    args = SimpleNamespace(output_dir=None)

    cli._run_streaming_static_self(
        args,
        ["input.fa"],
        {"input.fa": ["first", "second"]},
        {},
        object(),
    )

    assert events == [
        "yield:first",
        "process:first",
        "yield:second",
        "process:second",
    ]


def test_streaming_self_runner_submits_only_record_descriptors_to_bounded_pool(
    monkeypatch,
):
    submitted = []

    class FinishedFuture:
        def result(self):
            return None

        def cancel(self):
            return False

    class FakeExecutor:
        def __init__(self, *, max_workers, mp_context):
            assert max_workers == 2
            assert mp_context.get_start_method() == "spawn"

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def submit(self, function, task):
            assert function is cli._process_static_self_task
            submitted.append(task)
            return FinishedFuture()

    monkeypatch.setattr(cli, "ProcessPoolExecutor", FakeExecutor)
    monkeypatch.setattr(cli, "as_completed", lambda futures: list(futures))
    monkeypatch.setattr(cli, "supports_indexed_fasta_access", lambda _path: True)
    monkeypatch.setattr(
        cli,
        "iter_fasta_records",
        lambda *_args, **_kwargs: pytest.fail(
            "the parent must not decode sequences for indexed worker tasks"
        ),
    )
    args = SimpleNamespace(output_dir=None, processes=2)

    cli._run_streaming_static_self(
        args,
        ["indexed.fa"],
        {"indexed.fa": ["chr1", "chr2", "chr3"]},
        {"chr2": ("chr2", 10, 20)},
        SimpleNamespace(command="moddotplot -f indexed.fa"),
    )

    assert [task[2] for task in submitted] == ["chr1", "chr2", "chr3"]
    assert submitted[1][3] == ("chr2", 10, 20)
    assert all(len(task) == 5 for task in submitted)


@pytest.mark.parametrize(
    ("requested", "records", "indexed", "cpus", "expected"),
    [
        (None, 25, True, 10, 2),
        (None, 3, True, 2, 2),
        (4, 2, True, 10, 2),
        (4, 25, False, 10, 1),
        (2, 1, True, 10, 1),
    ],
)
def test_streaming_process_count_is_bounded_and_index_aware(
    monkeypatch, requested, records, indexed, cpus, expected
):
    monkeypatch.setattr(cli.os, "cpu_count", lambda: cpus)
    args = SimpleNamespace(processes=requested)

    assert cli._streaming_process_count(args, records, indexed) == expected


def test_compute_only_auto_process_count_can_use_four_workers(monkeypatch):
    monkeypatch.setattr(cli.os, "cpu_count", lambda: 10)
    args = SimpleNamespace(processes=None, no_plot=True)

    assert cli._streaming_process_count(args, 25, indexed_access=True) == 4


@pytest.mark.parametrize("requested", [0, 5, "many"])
def test_streaming_process_count_rejects_invalid_values(requested):
    with pytest.raises(ValueError, match="1 through 4"):
        cli._streaming_process_count(
            SimpleNamespace(processes=requested), 2, indexed_access=True
        )


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


def test_main_without_arguments_defaults_to_static_parser(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["moddotplot"])

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    captured = capsys.readouterr()
    assert exc_info.value.code == 2
    assert "moddotplot static" in captured.err
    assert "one of the arguments -c/--config -l/--load -f/--fasta is required" in (
        captured.err
    )
    assert "the following arguments are required: command" not in captured.err


def test_static_main_errors_when_all_fasta_inputs_are_unreadable(
    monkeypatch, tmp_path, capsys
):
    missing = tmp_path / "missing.fa"
    monkeypatch.setattr(
        sys,
        "argv",
        ["moddotplot", "--fasta", str(missing), "--no-plot"],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    captured = capsys.readouterr()
    assert exc_info.value.code == 2
    assert "no readable FASTA input files remain" in captured.err


def test_parser_defaults_omitted_subcommand_to_static():
    args = cli.parse_args(["--fasta", "sequence.fa", "--no-plot"])

    assert args.command == "static"
    assert args.fasta == ["sequence.fa"]
    assert args.no_plot


def test_parser_preserves_explicit_interactive_subcommand():
    args = cli.parse_args(["interactive", "--fasta", "sequence.fa"])

    assert args.command == "interactive"


@pytest.mark.parametrize("option", ["-s", "--sequence"])
def test_static_parser_accepts_sequence_selection_aliases(option):
    args = cli.parse_args(["--fasta", "sequence.fa", option, "chr1", "chr2", "--grid"])

    assert args.command == "static"
    assert args.sequence == ["chr1", "chr2"]


def test_static_parser_accepts_bounded_process_request():
    args = cli.parse_args(["--fasta", "sequence.fa", "--processes", "3"])

    assert args.processes == 3


@pytest.mark.parametrize("value", [0, 5, True, "many"])
def test_process_count_validation_rejects_invalid_config_values(value):
    with pytest.raises(ValueError, match="integer from 1 through 4"):
        cli._validated_process_count(value)


def test_static_parser_rejects_process_requests_outside_public_bound():
    with pytest.raises(SystemExit) as exc_info:
        cli.parse_args(["--fasta", "sequence.fa", "--processes", "5"])

    assert exc_info.value.code == 2


def test_columnar_bedpe_writer_matches_legacy_text(tmp_path):
    matrix = np.array([[1.0, 0.91, 0.2], [0.91, 1.0, 0.88], [0.2, 0.88, 1.0]])
    kwargs = dict(
        window_size=10,
        id_threshold=86,
        x_name="chr1",
        y_name="chr1",
        self_identity=True,
        x_offset=1,
        y_offset=1,
        x_end=29,
        y_end=29,
    )
    expected_rows = cli.convertMatrixToBed(matrix, **kwargs)
    output = tmp_path / "matrix.bedpe"

    cli._write_matrix_bedpe(
        output,
        cli.iterMatrixToBedChunks(matrix, max_chunk_cells=2, **kwargs),
    )

    expected = "".join("\t".join(map(str, row)) + "\n" for row in expected_rows)
    assert output.read_text() == expected


def test_interactive_short_s_remains_save_flag():
    args = cli.get_parser().parse_args(["interactive", "--fasta", "sequence.fa", "-s"])

    assert args.save is True


def test_sequence_selection_prefers_exact_names_and_accepts_casefold_fallback():
    selected, names = cli._select_fasta_headers(
        {"genome.fa": ["Chr1", "chr1", "Chr2"]}, ["chr1", "CHR2"]
    )

    assert selected == {"genome.fa": ["chr1", "Chr2"]}
    assert names == ["chr1", "Chr2"]


@pytest.mark.parametrize(
    ("headers", "selectors", "message"),
    [
        ({"genome.fa": ["Chr1"]}, ["chr1", "Chr1"], "already requested"),
        ({"genome.fa": ["Chr1", "chr1"]}, ["CHR1"], "ambiguous"),
    ],
)
def test_sequence_selection_rejects_duplicate_or_ambiguous_requests(
    headers, selectors, message
):
    with pytest.raises(ValueError, match=message):
        cli._select_fasta_headers(headers, selectors)


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


def test_static_config_accepts_sequence_selection():
    args = cli.get_parser().parse_args(["static", "--fasta", "sequence.fa"])

    cli._apply_static_config(
        args,
        {"fasta": ["sequence.fa"], "sequence": ["chr1", "chr2"]},
    )

    assert args.sequence == ["chr1", "chr2"]


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


def test_interactive_load_passes_combined_beds_to_dash(monkeypatch, tmp_path):
    first_bed = tmp_path / "first.bed"
    second_bed = tmp_path / "second.bed"
    first_bed.write_text("chrA\t1010\t1020\n")
    second_bed.write_text("chrB\t30\t40\n")
    metadata = [
        {
            "x_name": "chrA:1001-1100",
            "y_name": "chrB",
            "x_size": 100,
            "y_size": 100,
            "self": False,
            "max_window_size": 50,
            "resolution": 2,
        }
    ]
    monkeypatch.setattr(
        cli, "extractFiles", lambda _path: ([[np.ones((2, 2))]], metadata)
    )
    dash_calls = []
    monkeypatch.setattr(cli, "run_dash", lambda *args: dash_calls.append(args))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "moddotplot",
            "interactive",
            "--load",
            "saved-matrices",
            "--bed",
            str(first_bed),
            str(second_bed),
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    assert exc_info.value.code == 0
    x_axis, y_axis = dash_calls[0][2][0]
    assert (x_axis[0], x_axis[-1]) == (1001, 1100)
    assert (y_axis[0], y_axis[-1]) == (0, 100)
    annotations = dash_calls[0][7]
    assert annotations[["chrom", "start", "end"]].to_dict("records") == [
        {"chrom": "chrA", "start": 1010, "end": 1020},
        {"chrom": "chrB", "start": 30, "end": 40},
    ]


def test_interactive_rejects_invalid_annotation_bed(monkeypatch, tmp_path, capsys):
    bed = tmp_path / "invalid.bed"
    bed.write_text("chrA\tnot-a-coordinate\t20\n")
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "moddotplot",
            "interactive",
            "--load",
            "saved-matrices",
            "--bed",
            str(bed),
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    assert exc_info.value.code == 2
    assert "Error reading annotation BED file(s)" in capsys.readouterr().err


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


def test_static_plot_directory_gets_reproducibility_summary(monkeypatch, tmp_path):
    fasta = tmp_path / "source genome.fa"
    fasta.touch()
    _patch_fasta_input(monkeypatch, ["chrA"], [[1] * 1000])
    monkeypatch.setattr(cli, "createSelfMatrix", lambda *_args: np.full((1, 1), 100.0))
    monkeypatch.setattr(
        cli,
        "convertMatrixToBed",
        lambda *_args, **_kwargs: [["header"], ["value"]],
    )

    def create_plot_files(**kwargs):
        prefix = tmp_path / "chrA:101-400" / "chrA:101-400"
        created = []
        for suffix in ("_FULL.svg", "_FULL.png", "_TRI.svg", "_TRI.png"):
            path = prefix.parent / f"{prefix.name}{suffix}"
            path.touch()
            created.append(str(path))
        return created

    monkeypatch.setattr(cli, "create_plots", create_plot_files)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "moddotplot",
            "-f",
            str(fasta),
            "--region",
            "chrA:101-400",
            "--resolution",
            "10",
            "--no-bedpe",
            "--no-hist",
            "--output-dir",
            str(tmp_path),
        ],
    )

    cli.main()

    summary = (tmp_path / "chrA:101-400" / "plot_summary.txt").read_text()
    assert f"Command: {shlex.join(sys.argv)}" in summary
    assert str(fasta.resolve()) in summary
    assert "Window sizes:\n  - 28 bp" in summary
    assert "Regions:\n  - chrA:101-400" in summary
    assert "BED annotation file: None" in summary
    assert (
        str((tmp_path / "chrA:101-400" / "chrA:101-400_TRI.svg").resolve()) in summary
    )


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
