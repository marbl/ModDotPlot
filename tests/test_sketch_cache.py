import gc
import sys
import weakref

import numpy as np
import pytest

import moddotplot.estimate_identity as estimate_identity
import moddotplot.moddotplot as cli
from moddotplot.estimate_identity import (
    ModimizerSketchCache,
    PreparedModimizerSketches,
    create_pairwise_matrix_from_sketches,
    create_self_matrix_from_sketches,
    createPairwiseMatrix,
    createSelfMatrix,
    prepare_modimizer_sketches,
)


def test_prepared_sketch_matrix_results_match_compatibility_apis():
    first = np.arange(1, 57, dtype=np.uint64)
    second = np.arange(29, 85, dtype=np.uint64)
    parameters = {
        "window_size": 20,
        "sparsity": 2,
        "delta": 0.5,
        "k": 5,
        "ambiguous": True,
        "expectation": 8,
    }

    prepared_first = prepare_modimizer_sketches(len(first), first, **parameters)
    prepared_second = prepare_modimizer_sketches(len(second), second, **parameters)

    expected_self = createSelfMatrix(
        len(first),
        first,
        parameters["window_size"],
        parameters["sparsity"],
        parameters["delta"],
        parameters["k"],
        0,
        parameters["ambiguous"],
        parameters["expectation"],
    )
    actual_self = create_self_matrix_from_sketches(
        prepared_first, parameters["k"], 0, parameters["ambiguous"]
    )
    np.testing.assert_array_equal(actual_self, expected_self)

    expected_pair = createPairwiseMatrix(
        len(first),
        len(second),
        first,
        second,
        parameters["window_size"],
        parameters["sparsity"],
        parameters["delta"],
        parameters["k"],
        0,
        parameters["ambiguous"],
        parameters["expectation"],
    )
    actual_pair = create_pairwise_matrix_from_sketches(
        prepared_first, prepared_second, 0, parameters["k"], True
    )
    np.testing.assert_array_equal(actual_pair, expected_pair)


def test_sketch_cache_reuses_exact_configuration_without_copying(monkeypatch):
    calls = []
    cache_sizes_during_prepare = []
    prepared = PreparedModimizerSketches(core=[{1}], neighbors=[{1, 2}])

    def fake_prepare(*args):
        calls.append(args)
        cache_sizes_during_prepare.append(len(cache._cache))
        return prepared

    monkeypatch.setattr(estimate_identity, "prepare_modimizer_sketches", fake_prepare)
    cache = ModimizerSketchCache(max_entries=2)
    arguments = ("sequence-a", 100, [1, 2, 3], 10, 2, 0.5, 5, False, 5)

    first = cache.get_or_prepare(*arguments)
    second = cache.get_or_prepare(*arguments)

    assert first is prepared
    assert second is first
    assert len(calls) == 1

    cache.get_or_prepare("sequence-a", 100, [1, 2, 3], 10, 2, 0.25, 5, False, 5)
    assert len(calls) == 2

    cache.get_or_prepare("sequence-b", 100, [1, 2, 3], 10, 2, 0.5, 5, False, 5)
    assert len(calls) == 3
    assert cache_sizes_during_prepare == [0, 1, 1]
    assert len(cache._cache) == 2


@pytest.mark.parametrize(
    ("lengths", "sizing_arguments", "expected_preparations"),
    [
        ((1000, 1000), ("--window", "100"), 2),
        # With default-style resolution sizing, the shorter self sketch is an
        # exact pairwise hit while the longer sequence needs the shorter
        # pairwise window: three preparations instead of four.
        ((1000, 900), ("--resolution", "10"), 3),
    ],
)
def test_two_sequence_grid_reuses_and_releases_prepared_sketches(
    monkeypatch, tmp_path, lengths, sizing_arguments, expected_preparations
):
    sequences = [[1] * lengths[0], [2] * lengths[1]]
    monkeypatch.setattr(cli, "isValidFasta", lambda _path: True)
    monkeypatch.setattr(cli, "getInputHeaders", lambda _path: ["chrA", "chrB"])
    monkeypatch.setattr(cli, "readKmersFromFile", lambda *_args: sequences)

    prepare_calls = []
    prepared_references = []

    def fake_prepare(*args):
        prepare_calls.append(args)
        prepared = PreparedModimizerSketches(core=[{1}], neighbors=[{1, 2}])
        prepared_references.append(weakref.ref(prepared))
        return prepared

    monkeypatch.setattr(estimate_identity, "prepare_modimizer_sketches", fake_prepare)
    monkeypatch.setattr(
        cli,
        "create_self_matrix_from_sketches",
        lambda *_args: np.full((1, 1), 100.0),
    )
    monkeypatch.setattr(
        cli,
        "create_pairwise_matrix_from_sketches",
        lambda *_args: np.full((1, 1), 95.0),
    )
    monkeypatch.setattr(
        cli,
        "convertMatrixToBed",
        lambda *_args, **_kwargs: [["header"], ["value"]],
    )

    def assert_sketches_released_before_render(**_kwargs):
        gc.collect()
        assert all(reference() is None for reference in prepared_references)

    monkeypatch.setattr(cli, "create_grid", assert_sketches_released_before_render)
    monkeypatch.setattr(cli, "create_plots", lambda **_kwargs: None)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "moddotplot",
            "static",
            "--fasta",
            "sequences.fa",
            "--grid-only",
            *sizing_arguments,
            "--no-bedpe",
            "--output-dir",
            str(tmp_path),
        ],
    )

    cli.main()

    assert len(prepare_calls) == expected_preparations
    assert {call[1][0] for call in prepare_calls} == {1, 2}
