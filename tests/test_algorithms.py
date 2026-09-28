import numpy as np
import pytest

from moddotplot.estimate_identity import (
    containment_neighbors,
    convertMatrixToBed,
    createSelfMatrix,
    pairwiseContainmentMatrix,
    partitionOverlaps,
    populateModimizers,
)
from moddotplot.parse_fasta import generateKmersFromFasta, printProgressBar


def test_bed_conversion_clamps_partial_windows_to_exact_region_end():
    bed = convertMatrixToBed(
        np.ones((2, 2)),
        window_size=100,
        id_threshold=80,
        x_name="query",
        y_name="reference",
        self_identity=False,
        x_offset=101,
        y_offset=201,
        x_end=250,
        y_end=350,
    )

    assert max(row[2] for row in bed[1:]) == 250
    assert max(row[5] for row in bed[1:]) == 350


def test_populate_modimizers_returns_denser_recursive_fallback():
    result = populateModimizers(
        partition=[1, 2, 3, 4],
        sparsity=4,
        ambiguous=True,
        expectation=4,
        k=1,
    )

    assert result == {2, 4}


def test_populate_modimizers_can_fall_back_all_the_way_to_sparsity_one():
    result = populateModimizers(
        partition=[1, 3, 5],
        sparsity=8,
        ambiguous=True,
        expectation=10,
        k=1,
    )

    assert result == {1, 3, 5}


def test_populate_modimizers_empty_partition_terminates_at_sparsity_one():
    assert populateModimizers([], 8, True, 10, 1) == set()


def test_partition_overlaps_uses_consistent_genomic_window_boundaries():
    # Twenty-six 5-mers represent a 30-base sequence.  Each 10-base window
    # contains six internal 5-mers; boundary-spanning k-mers are excluded.
    kmers = list(range(26))

    assert partitionOverlaps(kmers, win=10, delta=0, seq_len=26, k=5) == [
        list(range(0, 6)),
        list(range(10, 16)),
        list(range(20, 26)),
    ]


def test_partition_overlaps_expands_each_window_in_genomic_coordinates():
    kmers = list(range(26))

    assert partitionOverlaps(kmers, win=10, delta=0.5, seq_len=26, k=5) == [
        list(range(0, 11)),
        list(range(5, 21)),
        list(range(15, 26)),
    ]


def test_partition_overlaps_omits_trailing_fragment_without_a_kmer():
    # Eleven bases produce seven 5-mers.  With a 10-base window, the final
    # one-base fragment cannot produce another partition.
    assert partitionOverlaps(list(range(7)), win=10, delta=0, seq_len=7, k=5) == [
        list(range(6))
    ]


@pytest.mark.parametrize(("win", "k"), [(0, 5), (10, 0)])
def test_partition_overlaps_rejects_nonpositive_sizes(win, k):
    with pytest.raises(ValueError):
        partitionOverlaps([1, 2, 3], win=win, delta=0, seq_len=3, k=k)


def test_containment_uses_matching_flanks_outside_core_windows():
    core_a = {1, 2, 3, 4}
    core_b = {5, 6, 7, 8}

    result = containment_neighbors(
        core_a,
        core_b,
        core_a | core_b,
        core_a | core_b,
        identity=0,
        k=1,
    )

    # Neighbor expansion deliberately permits this match: it is the mechanism
    # used to recover repeats that straddle different partition boundaries.
    # Consequently, callers that need strictly core-local identity must use
    # delta=0 rather than silently expecting the expanded sketches to be
    # ignored.
    assert result == 1.0


def test_delta_half_recovers_repeat_shifted_across_window_boundary():
    # The {1, 2, 3, 4} repeat fills window 0, but its second occurrence starts
    # halfway through window 1 and ends halfway through window 2.  Core-only
    # comparison sees just half of it in window 2 and falls below the cutoff;
    # delta=0.5 expands window 2 far enough to recover the full repeat.
    hashes = [1, 2, 3, 4, 90, 91, 1, 2, 3, 4, 92, 93]

    without_neighbors = createSelfMatrix(len(hashes), hashes, 4, 1, 0, 1, 75, True, 4)
    with_neighbors = createSelfMatrix(len(hashes), hashes, 4, 1, 0.5, 1, 75, True, 4)

    assert without_neighbors[0, 2] == 0.0
    assert with_neighbors[0, 2] == 100.0


def test_neighbor_containment_applies_cutoff_after_both_directions():
    # A -> expanded B is 1/2, while B -> expanded A is 3/4.  The stronger
    # direction must be considered before applying the threshold.
    core_a = {1, 2}
    core_b = {3, 4, 5, 6}
    expanded_a = {1, 2, 3, 4, 5}
    expanded_b = {1}

    assert (
        containment_neighbors(core_a, core_b, expanded_a, expanded_b, identity=75, k=1)
        == 0.75
    )
    assert (
        containment_neighbors(core_a, core_b, expanded_a, expanded_b, identity=76, k=1)
        == 0.0
    )


def test_containment_cutoff_is_independent_of_argument_order():
    larger_sketch = {1, 2, 3, 4}
    contained_sketch = {1, 2}

    forward = containment_neighbors(
        larger_sketch,
        contained_sketch,
        larger_sketch,
        contained_sketch,
        identity=75,
        k=1,
    )
    reverse = containment_neighbors(
        contained_sketch,
        larger_sketch,
        contained_sketch,
        larger_sketch,
        identity=75,
        k=1,
    )

    assert forward == reverse == 1.0


def test_pairwise_containment_matrix_is_rectangular_and_keeps_axis_orientation():
    matrix = pairwiseContainmentMatrix(
        mod_set_x=[{1}, {2}, {3}],
        mod_set_y=[{1}, {3}],
        mod_set_x_neighbors=[{1}, {2}, {3}],
        mod_set_y_neighbors=[{1}, {3}],
        identity=0,
        k=1,
        supress_progress=True,
    )

    np.testing.assert_array_equal(
        matrix,
        np.array(
            [
                [100.0, 0.0, 0.0],
                [0.0, 0.0, 100.0],
            ]
        ),
    )
    assert matrix.shape == (2, 3)


def test_pairwise_containment_matrix_supports_more_rows_than_columns():
    matrix = pairwiseContainmentMatrix(
        mod_set_x=[{2}],
        mod_set_y=[{1}, {2}, {3}],
        mod_set_x_neighbors=[{2}],
        mod_set_y_neighbors=[{1}, {2}, {3}],
        identity=0,
        k=1,
        supress_progress=True,
    )

    np.testing.assert_array_equal(matrix, np.array([[0.0], [100.0], [0.0]]))
    assert matrix.shape == (3, 1)


@pytest.mark.parametrize(
    ("mod_set_x", "mod_set_y", "expected_shape"),
    [
        ([], [{1}, {2}], (2, 0)),
        ([{1}, {2}], [], (0, 2)),
        ([], [], (0, 0)),
    ],
)
def test_pairwise_containment_matrix_preserves_empty_axis_dimensions(
    mod_set_x, mod_set_y, expected_shape
):
    matrix = pairwiseContainmentMatrix(
        mod_set_x=mod_set_x,
        mod_set_y=mod_set_y,
        mod_set_x_neighbors=list(mod_set_x),
        mod_set_y_neighbors=list(mod_set_y),
        identity=0,
        k=1,
        supress_progress=True,
    )

    assert matrix.shape == expected_shape


def test_pairwise_containment_matrix_does_not_hide_misaligned_neighbor_data():
    with pytest.raises(IndexError):
        pairwiseContainmentMatrix(
            mod_set_x=[{1}, {2}],
            mod_set_y=[{1}],
            mod_set_x_neighbors=[{1}],
            mod_set_y_neighbors=[{1}],
            identity=0,
            k=1,
            supress_progress=True,
        )


@pytest.mark.parametrize("sequence", ["", "A", "AC"])
def test_generate_kmers_shorter_than_k_with_progress_returns_empty(sequence, capsys):
    assert list(generateKmersFromFasta(sequence, 3, quiet=False, fw_only=True)) == []
    assert "100.0%" in capsys.readouterr().out


def test_generate_one_kmer_with_progress_does_not_use_zero_modulus(capsys):
    result = list(generateKmersFromFasta("ACG", 3, quiet=False, fw_only=True))

    assert result == [np.uint64(0xB13A5310100F646E)]
    output = capsys.readouterr().out
    assert "100.0%" in output
    assert "Completed" in output


def test_print_progress_bar_accepts_zero_total(capsys):
    printProgressBar(0, 0, prefix="Progress:", suffix="Completed", length=4)

    output = capsys.readouterr().out
    assert "|████|" in output
    assert "100.0%" in output
