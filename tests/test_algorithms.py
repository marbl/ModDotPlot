import numpy as np
import pytest

from moddotplot.estimate_identity import (
    BEDPE_HEADER,
    containment_neighbors,
    convertMatrixToBed,
    convertMatrixToBedDataFrame,
    createSelfMatrix,
    iterMatrixToBedChunks,
    pairwiseContainmentMatrix,
    partitionOverlaps,
    populateModimizers,
)
from moddotplot.parse_fasta import generateKmersFromFasta, printProgressBar


def _scalar_bed_reference(
    matrix,
    window_size,
    id_threshold,
    x_name,
    y_name,
    self_identity,
    x_offset,
    y_offset,
    x_end=None,
    y_end=None,
):
    """Original scalar implementation used as a parity oracle."""

    bed = [BEDPE_HEADER]
    for x in range(matrix.shape[0]):
        for y in range(matrix.shape[1]):
            value = matrix[x, y]
            if self_identity and x > y:
                continue
            if not value >= id_threshold / 100:
                continue
            start_x = x * window_size + x_offset
            end_x = start_x + window_size - 1
            start_y = y * window_size + y_offset
            end_y = start_y + window_size - 1
            if x_end is not None:
                end_x = min(end_x, x_end)
            if y_end is not None:
                end_y = min(end_y, y_end)
            bed.append(
                (
                    x_name,
                    int(start_x),
                    int(end_x),
                    y_name,
                    int(start_y),
                    int(end_y),
                    float(value),
                )
            )
    return bed


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


@pytest.mark.parametrize("self_identity", [False, True])
def test_vectorized_bed_conversion_preserves_row_order_and_values(self_identity):
    matrix = np.array(
        [
            [0.0, 86.5, 0.0, 92.25],
            [87.0, 0.0, 99.0, 0.0],
            [0.0, 91.0, 88.0, 0.0],
        ]
    )
    expected = [
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
    for x in range(matrix.shape[0]):
        for y in range(matrix.shape[1]):
            value = matrix[x, y]
            if (not self_identity or x <= y) and value >= 86 / 100:
                expected.append(
                    (
                        "x",
                        x * 100 + 11,
                        min(x * 100 + 110, 250),
                        "y",
                        y * 100 + 21,
                        min(y * 100 + 120, 350),
                        float(value),
                    )
                )

    assert (
        convertMatrixToBed(
            matrix,
            window_size=100,
            id_threshold=86,
            x_name="x",
            y_name="y",
            self_identity=self_identity,
            x_offset=11,
            y_offset=21,
            x_end=250,
            y_end=350,
        )
        == expected
    )


def test_bed_conversion_returns_only_header_when_no_tiles_pass():
    bed = convertMatrixToBed(
        np.zeros((1_000, 1_000)),
        window_size=100,
        id_threshold=86,
        x_name="x",
        y_name="y",
        self_identity=True,
        x_offset=0,
        y_offset=0,
    )

    assert len(bed) == 1


@pytest.mark.parametrize("self_identity", [False, True])
@pytest.mark.parametrize("max_chunk_cells", [1, 7, 64, 10_000])
def test_chunked_bed_conversion_matches_scalar_reference_randomized(
    self_identity, max_chunk_cells
):
    rng = np.random.default_rng(709)
    # Exercise a non-contiguous view as well as values immediately around the
    # legacy threshold. NaN must remain filtered by the comparison.
    source = rng.uniform(0.0, 1.5, size=(14, 24))
    matrix = source[::2, 1::2]
    matrix[0, :4] = [0.859999, 0.86, 0.860001, np.nan]
    kwargs = dict(
        window_size=37,
        id_threshold=86,
        x_name="query",
        y_name="reference",
        self_identity=self_identity,
        x_offset=13,
        y_offset=29,
        x_end=251,
        y_end=411,
    )

    expected = _scalar_bed_reference(matrix, **kwargs)
    actual = convertMatrixToBed(matrix, max_chunk_cells=max_chunk_cells, **kwargs)

    assert actual == expected


def test_bed_chunks_are_strictly_bounded_across_row_and_column_boundaries():
    matrix = np.arange(30, dtype=float).reshape(3, 10)
    chunks = list(
        iterMatrixToBedChunks(
            matrix,
            window_size=10,
            id_threshold=0,
            x_name="x",
            y_name="y",
            self_identity=False,
            x_offset=0,
            y_offset=0,
            max_chunk_cells=4,
        )
    )

    # A ten-column row must be split because it is wider than the cap. The
    # concatenated chunks still follow exact C order across every boundary.
    assert [len(chunk) for chunk in chunks] == [4, 4, 2] * 3
    assert all(tuple(chunk.columns) == BEDPE_HEADER for chunk in chunks)
    observed = [
        tuple(row)
        for chunk in chunks
        for row in chunk.itertuples(index=False, name=None)
    ]
    assert observed == _scalar_bed_reference(matrix, 10, 0, "x", "y", False, 0, 0)[1:]


def test_bed_dataframe_helper_preserves_columns_clipping_and_float_values():
    matrix = np.array([[0.85, 0.9, 1.25], [0.95, 0.1, 1.5]], dtype=np.float32)
    kwargs = dict(
        window_size=10.5,
        id_threshold=86,
        x_name="chrQ",
        y_name="chrR",
        self_identity=False,
        x_offset=-3.25,
        y_offset=101.75,
        x_end=9.5,
        y_end=119.25,
        max_chunk_cells=2,
    )

    frame = convertMatrixToBedDataFrame(matrix, **kwargs)
    expected = _scalar_bed_reference(
        matrix,
        **{key: value for key, value in kwargs.items() if key != "max_chunk_cells"},
    )

    assert tuple(frame.columns) == BEDPE_HEADER
    assert list(frame.itertuples(index=False, name=None)) == expected[1:]
    assert frame["perID_by_events"].dtype == np.dtype(float)


@pytest.mark.parametrize("shape", [(0, 4), (4, 0), (0, 0)])
def test_empty_bed_dataframe_has_exact_schema(shape):
    frame = convertMatrixToBedDataFrame(
        np.empty(shape), 100, 86, "x", "y", True, 0, 0, max_chunk_cells=1
    )

    assert frame.empty
    assert tuple(frame.columns) == BEDPE_HEADER
    assert convertMatrixToBed(
        np.empty(shape), 100, 86, "x", "y", True, 0, 0, max_chunk_cells=1
    ) == [BEDPE_HEADER]


@pytest.mark.parametrize("max_chunk_cells", [0, -1, 1.5, True])
def test_bed_conversion_rejects_invalid_chunk_bound(max_chunk_cells):
    with pytest.raises((TypeError, ValueError), match="positive integer"):
        convertMatrixToBed(
            np.ones((1, 1)),
            100,
            86,
            "x",
            "y",
            False,
            0,
            0,
            max_chunk_cells=max_chunk_cells,
        )


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
        )


@pytest.mark.parametrize("legacy_progress_setting", [False, True])
def test_pairwise_legacy_progress_argument_is_silent(legacy_progress_setting, capsys):
    pairwiseContainmentMatrix(
        mod_set_x=[{1}],
        mod_set_y=[{1}],
        mod_set_x_neighbors=[{1}],
        mod_set_y_neighbors=[{1}],
        identity=0,
        k=1,
        supress_progress=legacy_progress_setting,
    )

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


@pytest.mark.parametrize("sequence", ["", "A", "AC"])
def test_generate_kmers_shorter_than_k_returns_empty_without_progress(sequence, capsys):
    assert list(generateKmersFromFasta(sequence, 3, quiet=False, fw_only=True)) == []
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_generate_one_kmer_does_not_emit_progress(capsys):
    result = list(generateKmersFromFasta("ACG", 3, quiet=False, fw_only=True))

    assert result == [np.uint64(0xB13A5310100F646E)]
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_legacy_progress_helper_is_a_silent_noop(capsys):
    assert printProgressBar(1, 1, prefix="Progress:", suffix="Completed") is None

    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
