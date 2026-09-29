import numpy as np
import pytest

from moddotplot import _nthash
from moddotplot.estimate_identity import (
    prepare_modimizer_sketches,
    prepare_sequence_sketches,
)
from moddotplot.parse_fasta import _hash_sequence


def _legacy_sketches(
    sequence, window_size, sparsity, delta, k, ambiguous, expectation, canonical
):
    hashes = _hash_sequence(
        sequence, k, fw_only=not canonical, ambiguous=ambiguous
    )
    return prepare_modimizer_sketches(
        len(hashes),
        hashes,
        window_size,
        sparsity,
        delta,
        k,
        ambiguous,
        expectation,
    )


def _assert_prepared_equal(actual, expected):
    assert len(actual.core) == len(expected.core)
    assert len(actual.neighbors) == len(expected.neighbors)
    for actual_sketch, expected_sketch in zip(actual.core, expected.core):
        np.testing.assert_array_equal(actual_sketch, expected_sketch)
    for actual_sketch, expected_sketch in zip(actual.neighbors, expected.neighbors):
        np.testing.assert_array_equal(actual_sketch, expected_sketch)


@pytest.mark.parametrize("canonical", [False, True])
@pytest.mark.parametrize("ambiguous", [False, True])
@pytest.mark.parametrize("delta", [0, 0.35, 0.5])
@pytest.mark.parametrize("sparsity", [1, 8, 64])
def test_native_sequence_sketches_match_legacy_pipeline(
    canonical, ambiguous, delta, sparsity
):
    rng = np.random.default_rng(20260929)
    sequence = "".join(rng.choice(list("ACGT"), size=713))
    sequence = sequence[:91] + "nRy" + sequence[94:351].lower() + "U" + sequence[352:]
    parameters = dict(
        window_size=73,
        sparsity=sparsity,
        delta=delta,
        k=11,
        ambiguous=ambiguous,
        expectation=91,
        canonical=canonical,
    )

    actual = prepare_sequence_sketches(sequence, **parameters)
    expected = _legacy_sketches(sequence, **parameters)

    _assert_prepared_equal(actual, expected)


@pytest.mark.parametrize(
    ("sequence", "window_size", "k", "expectation"),
    [
        ("A" * 401, 83, 11, 200),
        ("ACGT" * 10, 17, 21, 100),
        ("ACGT", 10, 21, 100),
        ("N" * 101, 31, 7, 100),
    ],
)
def test_native_sequence_sketches_match_adaptive_and_empty_edge_cases(
    sequence, window_size, k, expectation
):
    parameters = dict(
        window_size=window_size,
        sparsity=64,
        delta=0.5,
        k=k,
        ambiguous=False,
        expectation=expectation,
        canonical=True,
    )

    _assert_prepared_equal(
        prepare_sequence_sketches(sequence, **parameters),
        _legacy_sketches(sequence, **parameters),
    )


def test_sequence_sketch_path_does_not_materialize_positional_hashes(monkeypatch):
    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("the chromosome-wide positional hash API was used")

    monkeypatch.setattr(_nthash, "hash_kmers", fail_if_called)

    prepared = prepare_sequence_sketches(
        "ACGT" * 10_000,
        window_size=1_000,
        sparsity=64,
        delta=0.5,
        k=21,
        ambiguous=False,
        expectation=16,
    )

    assert prepared.core
    assert prepared.neighbors
    assert all(sketch.dtype == np.uint64 for sketch in prepared.core)


@pytest.mark.parametrize("sparsity", [0, -1, 3, 12])
def test_native_sequence_sketches_require_power_of_two_sparsity(sparsity):
    with pytest.raises(ValueError, match="positive power of two"):
        prepare_sequence_sketches(
            "ACGTACGT",
            window_size=4,
            sparsity=sparsity,
            delta=0.5,
            k=3,
            ambiguous=False,
            expectation=1,
        )


def test_native_sketch_api_rejects_out_of_range_interval_bounds():
    with pytest.raises(ValueError, match="interval bounds"):
        _nthash.sketch_kmers("ACGT", 3, True, [(0, 3)], 1, 1, False)


def test_native_intersections_read_each_sequence_length_only_once():
    packed = np.asarray([7], dtype=np.uint64).tobytes()

    class ChangingLengthSequence:
        def __init__(self):
            self.length_calls = 0

        def __len__(self):
            self.length_calls += 1
            return self.length_calls

        def __getitem__(self, index):
            if index == 0:
                return packed
            raise IndexError(index)

    left = ChangingLengthSequence()
    result = _nthash.intersection_counts(left, [packed])

    np.testing.assert_array_equal(np.frombuffer(result, dtype=np.int32), [1])
    assert left.length_calls == 1
