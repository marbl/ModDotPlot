import numpy as np
import pytest

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
