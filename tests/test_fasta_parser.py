import gzip
import os
import struct
import zlib

import numpy as np
import pytest
import moddotplot.parse_fasta as fasta_parser

from moddotplot.parse_fasta import (
    _iter_fasta_records,
    extractRegion,
    generateKmersFromFasta,
    getInputHeaders,
    getInputSeqLength,
    iter_fasta_records,
    isValidFasta,
    readKmersFromFile,
    supports_indexed_fasta_access,
)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("chrA:1-400", ("chrA", 1, 400)),
        (
            "PAN010.chr14.haplotype1.paternal:1-4000000",
            ("PAN010.chr14.haplotype1.paternal", 1, 4_000_000),
        ),
        ("sample-name:10-20", ("sample-name", 10, 20)),
        ("sample:1-400:101-300", ("sample", 101, 300)),
    ],
)
def test_extract_region_accepts_realistic_fasta_identifiers(value, expected):
    assert extractRegion(value) == expected


@pytest.mark.parametrize("value", ["sample", "sample:1", "sample:one-two"])
def test_extract_region_rejects_missing_or_malformed_coordinates(value):
    assert extractRegion(value) is None


def _bgzf_block(data):
    compressor = zlib.compressobj(level=6, wbits=-15)
    compressed = compressor.compress(data) + compressor.flush()
    block_size = 18 + len(compressed) + 8
    if block_size > 65_536:
        raise ValueError("test BGZF block is too large")

    header = struct.pack("<BBBBLBBH", 31, 139, 8, 4, 0, 0, 255, 6)
    extra = b"BC" + struct.pack("<HH", 2, block_size - 1)
    footer = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
    return header + extra + compressed + footer


def _indexed_fasta(records, line_bases=8):
    contents = bytearray()
    index_lines = []
    for name, sequence in records:
        contents.extend(f">{name} description\n".encode("ascii"))
        sequence_offset = len(contents)
        for start in range(0, len(sequence), line_bases):
            contents.extend(sequence[start : start + line_bases].encode("ascii"))
            contents.extend(b"\n")
        index_lines.append(
            f"{name}\t{len(sequence)}\t{sequence_offset}\t"
            f"{line_bases}\t{line_bases + 1}\n"
        )
    return bytes(contents), "".join(index_lines)


def _write_bgzf(path, uncompressed, block_size=29):
    blocks = []
    gzip_entries = []
    compressed_offset = 0
    uncompressed_offset = 0
    for start in range(0, len(uncompressed), block_size):
        if blocks:
            gzip_entries.append((compressed_offset, uncompressed_offset))
        data = uncompressed[start : start + block_size]
        block = _bgzf_block(data)
        blocks.append(block)
        compressed_offset += len(block)
        uncompressed_offset += len(data)
    blocks.append(_bgzf_block(b""))
    path.write_bytes(b"".join(blocks))
    gzi = struct.pack("<Q", len(gzip_entries)) + b"".join(
        struct.pack("<QQ", *entry) for entry in gzip_entries
    )
    path.with_name(path.name + ".gzi").write_bytes(gzi)


def test_indexed_access_capability_requires_all_compression_indexes(tmp_path):
    contents, fai = _indexed_fasta([("chr1", "ACGT" * 20)])

    plain = tmp_path / "plain.fa"
    plain.write_bytes(contents)
    assert not supports_indexed_fasta_access(plain)
    plain.with_name(plain.name + ".fai").write_text(fai)
    assert supports_indexed_fasta_access(plain)

    compressed = tmp_path / "indexed.fa.gz"
    _write_bgzf(compressed, contents)
    compressed.with_name(compressed.name + ".fai").write_text(fai)
    assert supports_indexed_fasta_access(compressed)
    compressed.with_name(compressed.name + ".gzi").unlink()
    assert not supports_indexed_fasta_access(compressed)


def test_plain_fasta_supports_wrapping_crlf_case_blanks_and_empty_record(tmp_path):
    fasta = tmp_path / "records.fa"
    fasta.write_bytes(
        b">alpha descriptive header\r\n"
        b"ac\r\n"
        b"\r\n"
        b"gT\r\n"
        b">beta another description\r\n"
        b"NN\r\n"
        b">empty"
    )

    assert list(_iter_fasta_records(fasta)) == [
        ("alpha", "acgT"),
        ("beta", "NN"),
        ("empty", ""),
    ]
    assert getInputHeaders(fasta) == ["alpha", "beta", "empty"]
    assert getInputSeqLength(fasta) == [4, 2, 0]
    assert isValidFasta(fasta) is True
    assert not (tmp_path / "records.fa.fai").exists()


def test_gzip_is_detected_by_content_instead_of_extension(tmp_path):
    fasta = tmp_path / "compressed.data"
    with gzip.open(fasta, "wt", encoding="ascii") as output:
        output.write(">alpha description\nAC\nGT\n>beta\ntt")

    assert list(_iter_fasta_records(fasta)) == [
        ("alpha", "ACGT"),
        ("beta", "tt"),
    ]
    assert getInputHeaders(fasta) == ["alpha", "beta"]
    assert getInputSeqLength(fasta) == [4, 2]


def test_bgzf_concatenated_blocks_are_read_as_one_fasta_stream(tmp_path):
    fasta = tmp_path / "records.fa.bgz"
    fasta.write_bytes(
        _bgzf_block(b">alpha description\nAC")
        + _bgzf_block(b"GT\n>beta\ntt\n")
        + _bgzf_block(b"")
    )

    assert list(_iter_fasta_records(fasta)) == [
        ("alpha", "ACGT"),
        ("beta", "tt"),
    ]


def test_bgzf_fai_supplies_headers_and_lengths_without_decompression(
    tmp_path, monkeypatch
):
    records = [("alpha", "ACGT" * 20), ("beta", "TGCA" * 18)]
    contents, fai = _indexed_fasta(records)
    fasta = tmp_path / "indexed.fa.gz"
    _write_bgzf(fasta, contents)
    fasta.with_name(fasta.name + ".fai").write_text(fai)

    def fail_decompression(*_args, **_kwargs):
        raise AssertionError("fresh compressed .fai metadata should be sufficient")

    monkeypatch.setattr(fasta_parser, "_open_fasta_text", fail_decompression)

    assert getInputHeaders(fasta) == ["alpha", "beta"]
    assert getInputSeqLength(fasta) == [80, 72]


def test_bgzf_gzi_fetches_requested_region_across_blocks(tmp_path, monkeypatch):
    records = [("alpha", "ACGT" * 20), ("beta", "TGCATGCA" * 16)]
    contents, fai = _indexed_fasta(records)
    fasta = tmp_path / "indexed.fa.bgz"
    _write_bgzf(fasta, contents, block_size=23)
    fasta.with_name(fasta.name + ".fai").write_text(fai)

    def fail_streaming(*_args, **_kwargs):
        raise AssertionError("indexed BGZF should not use sequential decompression")

    monkeypatch.setattr(fasta_parser, "_iter_selected_fasta_records", fail_streaming)
    region = ("beta", 7, 91)

    assert list(
        iter_fasta_records(
            fasta,
            regions={"beta": region},
            record_ids=["beta"],
        )
    ) == [("beta", records[1][1][6:91], "beta:7-91")]


@pytest.mark.parametrize("gzi_contents", [None, b"not-a-gzi-index"])
def test_bgzf_without_usable_gzi_keeps_sequential_fallback(
    tmp_path, monkeypatch, gzi_contents
):
    records = [("alpha", "ACGT" * 20), ("beta", "TGCA" * 18)]
    contents, fai = _indexed_fasta(records)
    fasta = tmp_path / "fallback.fa.gz"
    _write_bgzf(fasta, contents)
    fasta.with_name(fasta.name + ".fai").write_text(fai)
    gzi_path = fasta.with_name(fasta.name + ".gzi")
    if gzi_contents is None:
        gzi_path.unlink()
    else:
        gzi_path.write_bytes(gzi_contents)

    def fail_indexed_fetch(*_args, **_kwargs):
        raise AssertionError("unusable .gzi must select the streaming fallback")

    monkeypatch.setattr(fasta_parser, "_fetch_indexed_region", fail_indexed_fetch)

    assert list(iter_fasta_records(fasta, record_ids=["beta"])) == [
        ("beta", records[1][1], "beta")
    ]


def test_regular_gzip_fai_metadata_does_not_force_bgzf_random_access(
    tmp_path, monkeypatch
):
    records = [("alpha", "ACGT" * 20), ("beta", "TGCA" * 18)]
    contents, fai = _indexed_fasta(records)
    fasta = tmp_path / "ordinary.fa.gz"
    with gzip.open(fasta, "wb") as output:
        output.write(contents)
    fasta.with_name(fasta.name + ".fai").write_text(fai)

    def fail_indexed_fetch(*_args, **_kwargs):
        raise AssertionError("ordinary gzip is not BGZF seekable")

    monkeypatch.setattr(fasta_parser, "_fetch_indexed_region", fail_indexed_fetch)

    assert getInputHeaders(fasta) == ["alpha", "beta"]
    assert list(iter_fasta_records(fasta, record_ids=["alpha"])) == [
        ("alpha", records[0][1], "alpha")
    ]


def test_stale_bgzf_gzi_keeps_sequential_fallback(tmp_path, monkeypatch):
    records = [("alpha", "ACGT" * 20)]
    contents, fai = _indexed_fasta(records)
    fasta = tmp_path / "stale.fa.gz"
    _write_bgzf(fasta, contents)
    fasta.with_name(fasta.name + ".fai").write_text(fai)
    os.utime(fasta.with_name(fasta.name + ".gzi"), (1, 1))

    def fail_indexed_fetch(*_args, **_kwargs):
        raise AssertionError("a stale .gzi must not be used for random access")

    monkeypatch.setattr(fasta_parser, "_fetch_indexed_region", fail_indexed_fetch)

    assert list(iter_fasta_records(fasta)) == [("alpha", records[0][1], "alpha")]


def test_public_record_iterator_preserves_requested_order_and_reports_missing(tmp_path):
    fasta = tmp_path / "records.fa"
    fasta.write_text(">alpha\nACGT\n>beta\nTTAA\n")

    assert list(iter_fasta_records(fasta, record_ids=["beta", "alpha"])) == [
        ("beta", "TTAA", "beta"),
        ("alpha", "ACGT", "alpha"),
    ]
    with pytest.raises(ValueError, match="FASTA record.*not found: 'missing'"):
        list(iter_fasta_records(fasta, record_ids=["missing"]))


def test_read_kmers_preserves_record_order_and_public_return_shape(tmp_path):
    fasta = tmp_path / "records.fa"
    fasta.write_text(">alpha description\nACGT\n>beta\nTTAA\n")

    result = readKmersFromFile(
        str(fasta), ksize=3, quiet=True, fw_only=True, ambiguous=False
    )

    assert isinstance(result, list)
    assert len(result) == 2
    assert all(isinstance(record, np.ma.MaskedArray) for record in result)
    assert result[0].tolist() == list(
        generateKmersFromFasta("ACGT", 3, quiet=True, fw_only=True)
    )
    assert result[1].tolist() == list(
        generateKmersFromFasta("TTAA", 3, quiet=True, fw_only=True)
    )


def test_read_kmers_hashes_only_selected_records_from_unindexed_fasta(tmp_path):
    sequences = {
        "Chr1": "ACGTAC",
        "Chr2": "TTAACC",
        "Chr3": "GGGGGG",
    }
    fasta = tmp_path / "selected.fa"
    fasta.write_text(
        "".join(f">{name}\n{sequence}\n" for name, sequence in sequences.items())
    )

    result = readKmersFromFile(
        str(fasta),
        ksize=3,
        quiet=True,
        fw_only=True,
        ambiguous=False,
        record_ids=["Chr1", "Chr2"],
    )

    assert len(result) == 2
    for hashes, sequence in zip(result, (sequences["Chr1"], sequences["Chr2"])):
        assert hashes.tolist() == list(
            generateKmersFromFasta(sequence, 3, quiet=True, fw_only=True)
        )


def test_read_kmers_uses_index_for_selected_record_subset(tmp_path, monkeypatch):
    fasta = tmp_path / "selected-indexed.fa"
    fasta.write_bytes(b">alpha\nACGTAC\n>beta\nTTAACC\n>gamma\nGGGGGG\n")
    (tmp_path / "selected-indexed.fa.fai").write_text(
        "alpha\t6\t7\t6\t7\n" "beta\t6\t20\t6\t7\n" "gamma\t6\t34\t6\t7\n"
    )

    def fail_streaming(*_args, **_kwargs):
        raise AssertionError("indexed subset should not use the streaming fallback")

    monkeypatch.setattr(fasta_parser, "_iter_selected_fasta_records", fail_streaming)
    result = readKmersFromFile(
        str(fasta),
        ksize=3,
        quiet=True,
        fw_only=True,
        ambiguous=False,
        record_ids=["beta", "alpha"],
    )

    assert [hashes.tolist() for hashes in result] == [
        list(generateKmersFromFasta(sequence, 3, quiet=True, fw_only=True))
        for sequence in ("TTAACC", "ACGTAC")
    ]


def test_read_kmers_hashes_only_the_requested_region(tmp_path):
    sequence = "ACGT" * 250
    fasta = tmp_path / "region.fa"
    fasta.write_text(f">sample.with.dots\n{sequence}\n")

    result = readKmersFromFile(
        str(fasta),
        ksize=21,
        quiet=True,
        fw_only=True,
        ambiguous=False,
        regions={"sample.with.dots": ("sample.with.dots", 101, 400)},
    )

    assert len(result[0]) == 280
    assert result[0].tolist() == list(
        generateKmersFromFasta(sequence[100:400], 21, quiet=True, fw_only=True)
    )


def test_indexed_region_uses_random_access_instead_of_record_parser(
    tmp_path, monkeypatch
):
    sequence = "ACGT" * 250
    fasta = tmp_path / "indexed.fa"
    header = b">sample.with.dots\n"
    fasta.write_bytes(header + sequence.encode("ascii") + b"\n")
    (tmp_path / "indexed.fa.fai").write_text(
        f"sample.with.dots\t{len(sequence)}\t{len(header)}\t{len(sequence)}\t{len(sequence) + 1}\n"
    )

    def fail_streaming(*_args, **_kwargs):
        raise AssertionError("indexed FASTA should not use the streaming fallback")

    monkeypatch.setattr(fasta_parser, "_iter_selected_fasta_records", fail_streaming)
    result = readKmersFromFile(
        str(fasta),
        ksize=21,
        quiet=True,
        fw_only=True,
        regions={"sample.with.dots": ("sample.with.dots", 101, 400)},
        record_ids=["sample.with.dots"],
    )

    assert len(result[0]) == 280
    assert result[0].tolist() == list(
        generateKmersFromFasta(sequence[100:400], 21, quiet=True, fw_only=True)
    )


def test_single_record_stream_stops_after_requested_region(tmp_path):
    fasta = tmp_path / "streamed.fa"
    fasta.write_text(">sample\nACGT\nACGT\nBRO KEN\n")

    result = readKmersFromFile(
        str(fasta),
        ksize=3,
        quiet=True,
        fw_only=True,
        regions={"sample": ("sample", 1, 8)},
        record_ids=["sample"],
    )

    assert result[0].tolist() == list(
        generateKmersFromFasta("ACGTACGT", 3, quiet=True, fw_only=True)
    )


def test_header_reader_does_not_assemble_sequence_records(tmp_path, monkeypatch):
    fasta = tmp_path / "headers.fa"
    fasta.write_text(">alpha\nACGT\n>beta\nTTAA\n")

    monkeypatch.setattr(
        fasta_parser,
        "_iter_fasta_records",
        lambda *_args: (_ for _ in ()).throw(AssertionError("full parser used")),
    )

    assert getInputHeaders(fasta) == ["alpha", "beta"]


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        (b"ACGT\n", "sequence data before the first header"),
        (b">   \nACGT\n", "empty header"),
        (b">alpha\nAC GT\n", "whitespace within sequence data"),
        (b">alpha\nAC\n>alpha description\nGT\n", "duplicate sequence identifier"),
        (b"\n\n", "no FASTA records found"),
    ],
)
def test_malformed_fasta_has_a_clear_error(tmp_path, contents, message):
    fasta = tmp_path / "malformed.fa"
    fasta.write_bytes(contents)

    with pytest.raises(ValueError, match=message):
        list(_iter_fasta_records(fasta))


def test_public_header_reader_rejects_duplicate_first_token_ids(tmp_path):
    fasta = tmp_path / "duplicate.fa"
    fasta.write_text(">alpha first\nAC\n>alpha second\nGT\n")

    with pytest.raises(ValueError, match="duplicate sequence identifier 'alpha'"):
        getInputHeaders(fasta)


def test_validation_returns_false_and_reports_malformed_fasta(tmp_path, capsys):
    fasta = tmp_path / "malformed.fa"
    fasta.write_text(">\nACGT\n")

    assert isValidFasta(fasta) is False
    assert "empty header at line 1" in capsys.readouterr().out


def test_validation_preserves_missing_file_exit_code(tmp_path, capsys):
    missing = tmp_path / "missing.fa"

    with pytest.raises(SystemExit) as error:
        isValidFasta(missing)

    assert error.value.code == 5
    assert "Unable to find fasta" in capsys.readouterr().out


def test_validation_preserves_unreadable_compression_exit_code(tmp_path, capsys):
    fasta = tmp_path / "corrupt.fa.gz"
    fasta.write_bytes(b"\x1f\x8bnot-a-valid-gzip-stream")

    with pytest.raises(SystemExit) as error:
        isValidFasta(fasta)

    assert error.value.code == 6
    assert "An error occurred" in capsys.readouterr().out
