import numpy as np
import pytest

from moddotplot import _nthash
from moddotplot.estimate_identity import partitionOverlaps, populateModimizers
from moddotplot.parse_fasta import (
    HASH_ALGORITHM,
    _hash_sequence,
    generateKmersFromFasta,
)


UPSTREAM_SEQUENCE = "ACATGCATGCA"
UPSTREAM_CANONICAL_K5 = np.array(
    [
        0xF59ECB45F0E22B9C,
        0x38CC00F940AEBDAE,
        0x603A48C5A11C794A,
        0x603A48C5A11C794A,
        0x38CC00F940AEBDAE,
        0x38CC00F940AEBDAE,
        0x603A48C5A11C794A,
    ],
    dtype=np.uint64,
)
REFERENCE_SEQUENCE = "ACGTCGTCAGTCGATGCAGT"
REFERENCE_FORWARD_K5 = np.array(
    [
        12092198163974795865,
        18369466880275892494,
        16642930255818702813,
        18031069624071434367,
        16308855512700874334,
        15484339800684336654,
        9427241265389440194,
        2202942986454218781,
        3417906301597334153,
        17521573918962941951,
        14501826884856700705,
        10960263084059400115,
        3737613185721920956,
        6733294302466194309,
        842074597590751927,
        1487097976980995052,
    ],
    dtype=np.uint64,
)
REFERENCE_CANONICAL_K5 = np.array(
    [
        17280910991238361383,
        12143167639373124245,
        8315064964825903610,
        6347800447189266283,
        600083231416985220,
        8268053965520490328,
        2337193125131410193,
        8232026757601677174,
        15014898482287274366,
        14172688858799321146,
        141383167342963224,
        12915775154031495997,
        5947576043254925263,
        6933934589939054922,
        10998487937625170617,
        11900314154005552483,
    ],
    dtype=np.uint64,
)
IUPAC_COMPLEMENT = str.maketrans(
    "ACGTURYMKSWHBVDNacgturymkswhbvdn",
    "TGCAAYRKMSWDVBHNt gcaayrkmswdvbhn".replace(" ", ""),
)


def _unpack_hashes(raw_hashes):
    return np.frombuffer(raw_hashes, dtype=np.uint64)


def _reverse_complement(sequence):
    return sequence.translate(IUPAC_COMPLEMENT)[::-1]


def test_native_hashes_match_official_nthash2_v2_vectors():
    """Values come from ntHash v2.4.0's upstream tests/tests.cpp."""
    raw_hashes, raw_mask = _nthash.hash_kmers(UPSTREAM_SEQUENCE, 5, True)
    hashes = _unpack_hashes(raw_hashes)

    assert HASH_ALGORITHM == _nthash.ALGORITHM == "ntHash_v2"
    assert _nthash.UPSTREAM_VERSION == "2.4.0"
    assert _nthash.UPSTREAM_COMMIT == "c26bd4572a19de81e30d55042dbd33c1fd21d4b6"
    assert raw_mask == b""
    np.testing.assert_array_equal(hashes, UPSTREAM_CANONICAL_K5)
    # These are the two base canonical hashes explicitly fixed by upstream.
    assert hashes[1] == np.uint64(0x38CC00F940AEBDAE)
    assert hashes[2] == np.uint64(0x603A48C5A11C794A)


@pytest.mark.parametrize(
    ("canonical", "fw_only", "expected"),
    [
        (False, True, REFERENCE_FORWARD_K5),
        (True, False, REFERENCE_CANONICAL_K5),
    ],
)
def test_native_and_public_paths_match_exact_20_base_vectors(
    canonical, fw_only, expected
):
    raw_hashes, raw_mask = _nthash.hash_kmers(REFERENCE_SEQUENCE, 5, canonical)
    public_hashes = _hash_sequence(
        REFERENCE_SEQUENCE, 5, fw_only=fw_only, ambiguous=False
    )

    assert raw_mask == b""
    np.testing.assert_array_equal(_unpack_hashes(raw_hashes), expected)
    np.testing.assert_array_equal(public_hashes.data, expected)


def test_public_hash_helper_preserves_official_values_and_uint64_dtype():
    hashes = _hash_sequence(UPSTREAM_SEQUENCE, 5, fw_only=False, ambiguous=False)

    assert isinstance(hashes, np.ma.MaskedArray)
    assert hashes.dtype == np.dtype(np.uint64)
    assert hashes.shape == (len(UPSTREAM_SEQUENCE) - 5 + 1,)
    assert not np.ma.getmaskarray(hashes).any()
    np.testing.assert_array_equal(hashes.data, UPSTREAM_CANONICAL_K5)
    assert any(int(value) > np.iinfo(np.uint32).max for value in hashes)


@pytest.mark.parametrize("as_bytes", [False, True])
def test_native_canonical_hashes_are_reverse_complement_invariant(as_bytes):
    sequence = "AACCGTACACTGGACTGAGTCT"
    reverse_complement = _reverse_complement(sequence)
    if as_bytes:
        sequence = sequence.encode("ascii")
        reverse_complement = reverse_complement.encode("ascii")

    forward_raw, forward_mask = _nthash.hash_kmers(sequence, 7, True)
    reverse_raw, reverse_mask = _nthash.hash_kmers(reverse_complement, 7, True)

    assert forward_mask == reverse_mask == b""
    np.testing.assert_array_equal(
        _unpack_hashes(forward_raw), _unpack_hashes(reverse_raw)[::-1]
    )


def test_forward_hashes_remain_strand_specific():
    sequence = "AACCGTACACTGGACTGAGTCT"
    reverse_complement = _reverse_complement(sequence)

    forward = _hash_sequence(sequence, 7, fw_only=True).data
    reverse = _hash_sequence(reverse_complement, 7, fw_only=True).data[::-1]

    assert np.any(forward != reverse)


def test_native_mask_marks_every_and_only_ambiguous_window():
    sequence = "ACGTNRYACGT"
    k = 3
    raw_hashes, raw_mask = _nthash.hash_kmers(sequence, k, True)

    assert len(raw_hashes) == (len(sequence) - k + 1) * np.dtype(np.uint64).itemsize
    assert list(raw_mask) == [0, 0, 1, 1, 1, 1, 1, 0, 0]
    assert len(_unpack_hashes(raw_hashes)) == len(raw_mask) == 9


def test_public_ambiguity_flag_masks_or_retains_position_preserving_fallbacks():
    sequence = "ACGTNRYACGT"
    expected_mask = np.array([0, 0, 1, 1, 1, 1, 1, 0, 0], dtype=bool)

    excluded = _hash_sequence(sequence, 3, fw_only=False, ambiguous=False)
    retained = _hash_sequence(sequence, 3, fw_only=False, ambiguous=True)

    assert len(excluded) == len(retained) == len(sequence) - 3 + 1
    np.testing.assert_array_equal(np.ma.getmaskarray(excluded), expected_mask)
    assert not np.ma.getmaskarray(retained).any()
    np.testing.assert_array_equal(
        excluded.data[~expected_mask], retained.data[~expected_mask]
    )
    np.testing.assert_array_equal(
        excluded.data[expected_mask], retained.data[expected_mask]
    )


def test_generator_keeps_ambiguous_window_positions_as_none_by_default():
    sequence = "ACGTNRYACGT"

    excluded = list(
        generateKmersFromFasta(sequence, 3, quiet=True, fw_only=False, ambiguous=False)
    )
    retained = list(
        generateKmersFromFasta(sequence, 3, quiet=True, fw_only=False, ambiguous=True)
    )

    assert len(excluded) == len(retained) == len(sequence) - 3 + 1
    assert [value is None for value in excluded] == [
        False,
        False,
        True,
        True,
        True,
        True,
        True,
        False,
        False,
    ]
    assert all(isinstance(value, int) for value in retained)


def test_ambiguous_canonical_fallback_is_reverse_complement_invariant():
    sequence = "AURYNACGTRYSWKMBDHVNACGTU"
    reverse_complement = _reverse_complement(sequence)

    forward = _hash_sequence(sequence, 5, fw_only=False, ambiguous=True)
    reverse = _hash_sequence(reverse_complement, 5, fw_only=False, ambiguous=True)

    np.testing.assert_array_equal(forward.data, reverse.data[::-1])


@pytest.mark.parametrize("sequence", ["", "A", "AC"])
def test_short_sequences_return_empty_uint64_results(sequence):
    raw_hashes, raw_mask = _nthash.hash_kmers(sequence, 3, True)
    hashes = _hash_sequence(sequence, 3, fw_only=False)

    assert raw_hashes == raw_mask == b""
    assert hashes.shape == (0,)
    assert hashes.dtype == np.dtype(np.uint64)
    assert list(generateKmersFromFasta(sequence, 3, quiet=True, fw_only=False)) == []


@pytest.mark.parametrize("k", [0, -1, 65536])
def test_native_rejects_unsupported_kmer_sizes(k):
    with pytest.raises(ValueError, match="k must be between 1 and 65535"):
        _nthash.hash_kmers("ACGT", k, True)


def test_hashing_is_ascii_case_insensitive_and_treats_u_as_t():
    dna = _hash_sequence("ACGTTACG", 4, fw_only=False)
    lower = _hash_sequence("acgttacg", 4, fw_only=False)
    rna = _hash_sequence("ACGUUACG", 4, fw_only=False)

    np.testing.assert_array_equal(dna.data, lower.data)
    np.testing.assert_array_equal(dna.data, rna.data)


def test_mask_survives_window_partitioning_without_coordinate_collapse():
    hashes = _hash_sequence("AAAAANAAAAA", 3, fw_only=False, ambiguous=False)

    partitions = partitionOverlaps(hashes, win=6, delta=0, seq_len=len(hashes), k=3)

    assert [len(partition) for partition in partitions] == [4, 3]
    assert [np.ma.getmaskarray(partition).tolist() for partition in partitions] == [
        [False, False, False, True],
        [False, False, False],
    ]


def test_modimizer_population_skips_masked_hashes_and_returns_python_ints():
    partition = np.ma.MaskedArray(
        np.array([2, 100, 4, 200, 6], dtype=np.uint64),
        mask=[False, True, False, True, False],
    )

    result = populateModimizers(
        partition, sparsity=2, ambiguous=False, expectation=1, k=3
    )

    assert result == {2, 4, 6}
    assert all(type(value) is int for value in result)


def test_ambiguity_flag_controls_whether_fallbacks_enter_modimizer_sketches():
    excluded = _hash_sequence("NNNNN", 3, fw_only=False, ambiguous=False)
    retained = _hash_sequence("NNNNN", 3, fw_only=False, ambiguous=True)

    excluded_modimizers = populateModimizers(
        excluded, sparsity=1, ambiguous=False, expectation=1, k=3
    )
    retained_modimizers = populateModimizers(
        retained, sparsity=1, ambiguous=True, expectation=1, k=3
    )

    assert excluded_modimizers == set()
    assert len(retained_modimizers) == 1
