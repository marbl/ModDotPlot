import gzip
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
    isValidFasta,
    readKmersFromFile,
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
