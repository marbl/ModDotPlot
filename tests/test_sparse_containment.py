import numpy as np
import pytest

from moddotplot import _nthash
import moddotplot.estimate_identity as estimate_identity
from moddotplot.estimate_identity import (
    _sketch_intersection_counts,
    pairwiseContainmentMatrix,
    prepare_modimizer_sketches,
    selfContainmentMatrix,
)


def _scalar_identity(core_a, core_b, expanded_a, expanded_b, identity, k):
    core_a = set(core_a)
    core_b = set(core_b)
    expanded_a = set(expanded_a)
    expanded_b = set(expanded_b)
    a_to_b = len(core_a & expanded_b) / len(core_a) if core_a else 0.0
    b_to_a = len(core_b & expanded_a) / len(core_b) if core_b else 0.0
    estimated_identity = max(a_to_b, b_to_a) ** (1.0 / k)
    return estimated_identity * 100 if estimated_identity >= identity / 100 else 0.0


def _scalar_pairwise(core_x, core_y, expanded_x, expanded_y, identity, k):
    return np.asarray(
        [
            [
                _scalar_identity(
                    core_x[x], core_y[y], expanded_x[x], expanded_y[y], identity, k
                )
                for x in range(len(core_x))
            ]
            for y in range(len(core_y))
        ],
        dtype=float,
    ).reshape(len(core_y), len(core_x))


def _random_sketches(rng, count, universe, maximum_size):
    sketches = []
    expanded = []
    for _ in range(count):
        size = int(rng.integers(0, maximum_size + 1))
        core = set(rng.choice(universe, size=size, replace=False).tolist())
        additions = set(
            rng.choice(
                universe,
                size=int(rng.integers(0, maximum_size + 1)),
                replace=False,
            ).tolist()
        )
        sketches.append(core)
        expanded.append(core | additions)
    return sketches, expanded


@pytest.mark.parametrize(("identity", "k"), [(0, 1), (75, 5), (86, 21), (99, 31)])
def test_sparse_pairwise_matches_scalar_reference(identity, k):
    rng = np.random.default_rng(911)
    core_x, expanded_x = _random_sketches(rng, 7, 200, 30)
    core_y, expanded_y = _random_sketches(rng, 5, 200, 30)

    expected = _scalar_pairwise(core_x, core_y, expanded_x, expanded_y, identity, k)
    actual = pairwiseContainmentMatrix(
        core_x,
        core_y,
        expanded_x,
        expanded_y,
        identity,
        k,
        supress_progress=True,
    )

    np.testing.assert_allclose(actual, expected)
    assert actual.shape == (5, 7)


@pytest.mark.parametrize("ambiguous", [False, True])
def test_sparse_self_matches_scalar_reference_including_empty_diagonal(ambiguous):
    rng = np.random.default_rng(77)
    core, expanded = _random_sketches(rng, 8, 150, 25)
    core[3] = set()
    expanded[3] = set()
    expected = _scalar_pairwise(core, core, expanded, expanded, 86, 21)
    np.fill_diagonal(expected, 100.0)
    if not ambiguous:
        expected[3, 3] = 0.0

    actual = selfContainmentMatrix(core, expanded, 21, 86, ambiguous)

    np.testing.assert_allclose(actual, expected)


def test_sparse_hash_compression_does_not_alias_out_of_range_python_ints():
    # Casting these values to uint64 would alias -1 and 2**64 - 1. The generic
    # compatibility path must continue to treat them as distinct hashes.
    sketches_a = [{-1}, {2**64 - 1}, {2**64 + 1}]
    sketches_b = [{-1}, {2**64 - 1}, {1}]

    np.testing.assert_array_equal(
        _sketch_intersection_counts(sketches_a, sketches_b),
        np.diag([1, 1, 0]).astype(np.int32),
    )


def test_sparse_counts_use_wide_accumulator():
    # uint8 sparse multiplication silently wraps 300 to 44.
    sketch = set(range(300))
    counts = _sketch_intersection_counts([sketch], [sketch])

    assert counts.dtype == np.int32
    assert counts[0, 0] == 300


def test_native_sorted_uint64_intersections_match_scalar_reference(monkeypatch):
    rng = np.random.default_rng(616)
    sketches_a = [
        np.sort(
            rng.choice(5_000, size=int(rng.integers(0, 250)), replace=False)
        ).astype(np.uint64)
        for _ in range(17)
    ]
    sketches_b = [
        np.sort(
            rng.choice(5_000, size=int(rng.integers(0, 250)), replace=False)
        ).astype(np.uint64)
        for _ in range(13)
    ]
    expected = np.asarray(
        [
            [len(set(left.tolist()) & set(right.tolist())) for right in sketches_b]
            for left in sketches_a
        ],
        dtype=np.int32,
    )

    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("the CSR compatibility path was used")

    monkeypatch.setattr(estimate_identity, "csr_matrix", fail_if_called)

    actual = _sketch_intersection_counts(sketches_a, sketches_b)

    np.testing.assert_array_equal(actual, expected)
    assert actual.dtype == np.int32


def test_native_intersection_merge_handles_hash_shared_by_every_window():
    shared = np.uint64(2**63 + 17)
    sketches_a = [np.array([index, shared], dtype=np.uint64) for index in range(20)]
    sketches_b = [
        np.array([index + 100, shared], dtype=np.uint64) for index in range(30)
    ]
    # Keep the production precondition explicit: arrays are sorted and unique.
    sketches_a = [np.sort(sketch) for sketch in sketches_a]
    sketches_b = [np.sort(sketch) for sketch in sketches_b]

    np.testing.assert_array_equal(
        _sketch_intersection_counts(sketches_a, sketches_b),
        np.ones((20, 30), dtype=np.int32),
    )


def test_native_intersection_retains_ephemeral_sequence_items():
    released = []

    class TrackedBytes(bytes):
        def __new__(cls, payload):
            instance = super().__new__(cls, payload)
            instance.release_events = released
            return instance

        def __del__(self):
            self.release_events.append(True)

    class EphemeralSketches:
        def __init__(self, sketches):
            self.sketches = sketches

        def __len__(self):
            return len(self.sketches)

        def __getitem__(self, index):
            # A sequence implementation is allowed to return a newly-created
            # object for each item. No item may be released while the native
            # function is still materializing or reading the sequence.
            assert not released
            return TrackedBytes(self.sketches[index])

    left = EphemeralSketches(
        [
            np.array([1, 3], dtype=np.uint64).tobytes(),
            np.array([2, 3], dtype=np.uint64).tobytes(),
        ]
    )
    right = EphemeralSketches(
        [
            np.array([3, 4], dtype=np.uint64).tobytes(),
            np.array([1, 2], dtype=np.uint64).tobytes(),
        ]
    )

    packed = _nthash.intersection_counts(left, right)

    np.testing.assert_array_equal(
        np.frombuffer(packed, dtype=np.int32).reshape(2, 2),
        np.array([[1, 1], [1, 1]], dtype=np.int32),
    )
    assert len(released) == 4


def test_unsorted_uint64_arrays_retain_compatibility_path():
    sketches_a = [np.array([9, 1, 5], dtype=np.uint64)]
    sketches_b = [np.array([5, 2, 9], dtype=np.uint64)]

    np.testing.assert_array_equal(
        _sketch_intersection_counts(sketches_a, sketches_b),
        np.array([[2]], dtype=np.int32),
    )


def test_matrix_path_does_not_fall_back_to_per_cell_set_intersections(monkeypatch):
    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("per-cell Python containment was used")

    monkeypatch.setattr(estimate_identity, "containment_neighbors", fail_if_called)
    core = [{index, index + 1} for index in range(250)]
    expanded = [sketch | {index + 2} for index, sketch in enumerate(core)]

    matrix = pairwiseContainmentMatrix(
        core, core, expanded, expanded, 0, 21, supress_progress=True
    )

    assert matrix.shape == (250, 250)


def test_prepared_sketches_use_compact_uint64_arrays():
    hashes = np.arange(20_000, dtype=np.uint64)
    prepared = prepare_modimizer_sketches(
        len(hashes),
        hashes,
        window_size=1_000,
        sparsity=8,
        delta=0.5,
        k=21,
        ambiguous=False,
        expectation=125,
    )

    sketches = prepared.core + prepared.neighbors
    assert all(isinstance(sketch, np.ndarray) for sketch in sketches)
    assert all(sketch.dtype == np.uint64 for sketch in sketches)
    assert sum(sketch.nbytes for sketch in sketches) == 8 * sum(
        len(sketch) for sketch in sketches
    )
