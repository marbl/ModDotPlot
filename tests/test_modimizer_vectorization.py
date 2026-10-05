import numpy as np
import pytest

import moddotplot.estimate_identity as estimate_identity
from moddotplot.estimate_identity import convertToModimizers, populateModimizers


def _scalar_reference(partition, sparsity, expectation):
    """Straightforward reference for the adaptive integer-sparsity algorithm."""

    current_sparsity = int(sparsity)
    while True:
        modimizers = {
            int(kmer)
            for kmer in partition
            if kmer is not None
            and not np.ma.is_masked(kmer)
            and int(kmer) % current_sparsity == 0
        }
        if len(modimizers) >= round(expectation / 2) or current_sparsity == 1:
            return modimizers
        current_sparsity = max(1, current_sparsity // 2)


@pytest.mark.parametrize("sparsity", [1, 2, 8, 64])
@pytest.mark.parametrize("expectation", [0, 1, 12, 200])
def test_vectorized_modimizers_match_scalar_reference(sparsity, expectation):
    rng = np.random.default_rng(441)
    values = rng.integers(0, 2**63, size=4096, dtype=np.uint64)
    values[100:160] = values[:60]
    mask = rng.random(values.size) < 0.13
    partition = np.ma.MaskedArray(values, mask=mask)

    assert populateModimizers(
        partition,
        sparsity=sparsity,
        ambiguous=False,
        expectation=expectation,
        k=21,
    ) == _scalar_reference(partition, sparsity, expectation)


def test_legacy_sequence_skips_none_and_masked_values():
    partition = [1, None, np.ma.masked, 2, np.uint64(4), 4, 7]

    assert populateModimizers(partition, 4, False, 4, 21) == {2, 4}


def test_adaptive_sparsity_remains_integer(monkeypatch):
    observed_sparsities = []
    original = estimate_identity._divisible_hashes

    def record_sparsity(values, sparsity):
        observed_sparsities.append(sparsity)
        return original(values, sparsity)

    monkeypatch.setattr(estimate_identity, "_divisible_hashes", record_sparsity)

    assert populateModimizers([], 8, False, 100, 21) == set()
    assert observed_sparsities == [8, 4, 2, 1]
    assert all(type(value) is int for value in observed_sparsities)


def test_convert_to_modimizers_preserves_partition_order():
    hashes = np.ma.MaskedArray(
        np.arange(24, dtype=np.uint64),
        mask=[False] * 8 + [True] * 4 + [False] * 12,
    )

    sketches = convertToModimizers(
        [hashes[:8], hashes[8:16], hashes[16:]],
        sparsity=4,
        ambiguous=False,
        k=5,
        expectation=1,
    )

    assert sketches == [{0, 4}, {12}, {16, 20}]
    assert all(type(value) is int for sketch in sketches for value in sketch)


def test_numpy_fast_path_does_not_iterate_python_scalars():
    class NonIterableArray(np.ndarray):
        def __iter__(self):
            raise AssertionError("the NumPy fast path must not iterate in Python")

    values = np.arange(100_000, dtype=np.uint64).view(NonIterableArray)

    result = populateModimizers(values, 64, False, 1, 21)

    assert result == set(range(0, 100_000, 64))


def test_nomask_fast_path_does_not_materialize_a_mask(monkeypatch):
    partition = np.ma.MaskedArray(
        np.arange(100_000, dtype=np.uint64), mask=np.ma.nomask
    )

    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("nomask should not be expanded into a boolean array")

    monkeypatch.setattr(estimate_identity.np.ma, "getmaskarray", fail_if_called)

    assert populateModimizers(partition, 64, False, 1, 21) == set(range(0, 100_000, 64))


def test_uniqueness_work_scales_with_selected_candidates(monkeypatch):
    # Filtering must happen before uniqueness. At the production-like sparsity
    # below, only 1/1024 of the input should reach the more expensive sort.
    values = np.arange(2**20, dtype=np.uint64)
    unique_input_sizes = []
    original = estimate_identity.np.unique

    def record_unique_input(selected):
        unique_input_sizes.append(selected.size)
        return original(selected)

    monkeypatch.setattr(estimate_identity.np, "unique", record_unique_input)

    result = populateModimizers(values, 1024, False, 1, 21)

    assert len(result) == 1024
    assert unique_input_sizes == [1024]


@pytest.mark.parametrize("sparsity", [0, -1])
def test_modimizers_reject_nonpositive_sparsity(sparsity):
    with pytest.raises(ValueError, match="positive integer"):
        populateModimizers([1, 2, 3], sparsity, False, 1, 21)
