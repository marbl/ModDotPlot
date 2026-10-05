from typing import (
    Iterable,
    Iterator,
    List,
    NamedTuple,
    Optional,
    Sequence,
    TextIO,
    Tuple,
)
from bisect import bisect_right
import sys
import os
import pickle
import re
import struct
import numpy as np
import gzip

from moddotplot import _nthash

HASH_ALGORITHM = _nthash.ALGORITHM


class FastaIndexEntry(NamedTuple):
    name: str
    length: int
    offset: int
    line_bases: int
    line_width: int


class BgzfIndexEntry(NamedTuple):
    compressed_offset: int
    uncompressed_offset: int


def _is_gzip(filename: str) -> bool:
    with open(filename, "rb") as probe:
        return probe.read(2) == b"\x1f\x8b"


def _is_bgzf(filename: str) -> bool:
    """Return whether *filename* starts with a BGZF gzip member.

    BGZF is distinguished from ordinary gzip by the ``BC`` extra subfield in
    each member header.  Checking the header prevents an unrelated or stale
    ``.gzi`` file from making a normal gzip stream look seekable.
    """

    try:
        with open(filename, "rb") as compressed:
            fixed_header = compressed.read(12)
            if (
                len(fixed_header) != 12
                or fixed_header[:3] != b"\x1f\x8b\x08"
                or not fixed_header[3] & 0x04
            ):
                return False
            extra_length = struct.unpack_from("<H", fixed_header, 10)[0]
            extra = compressed.read(extra_length)
    except OSError:
        return False

    offset = 0
    while offset + 4 <= len(extra):
        subfield_id = extra[offset : offset + 2]
        subfield_length = struct.unpack_from("<H", extra, offset + 2)[0]
        offset += 4
        subfield_end = offset + subfield_length
        if subfield_end > len(extra):
            return False
        if subfield_id == b"BC" and subfield_length == 2:
            return True
        offset = subfield_end
    return False


def _open_fasta_text(filename: str) -> TextIO:
    """Open a plain, gzip, or BGZF FASTA file as text.

    Compression is detected from the gzip magic bytes rather than the filename
    extension. BGZF is a blocked form of gzip and is decoded by Python's gzip
    reader as a concatenated gzip stream.
    """
    if _is_gzip(filename):
        return gzip.open(filename, "rt", encoding="ascii", newline=None)
    return open(filename, "rt", encoding="ascii", newline=None)


def _read_fasta_index(filename: str) -> Optional[List[FastaIndexEntry]]:
    """Read a fresh samtools-style ``.fai`` index when one is available."""

    index_path = f"{os.fspath(filename)}.fai"
    if not os.path.isfile(index_path):
        return None
    if os.path.getmtime(index_path) < os.path.getmtime(filename):
        return None

    entries = []
    seen_names = set()
    try:
        with open(index_path, "rt", encoding="utf-8") as index:
            for line_number, raw_line in enumerate(index, start=1):
                fields = raw_line.rstrip("\r\n").split("\t")
                if len(fields) < 5:
                    return None
                name = fields[0]
                if not name or name in seen_names:
                    return None
                seen_names.add(name)
                entry = FastaIndexEntry(
                    name,
                    int(fields[1]),
                    int(fields[2]),
                    int(fields[3]),
                    int(fields[4]),
                )
                if entry.length < 0 or entry.offset < 0:
                    return None
                if entry.length and (
                    entry.line_bases <= 0 or entry.line_width < entry.line_bases
                ):
                    return None
                entries.append(entry)
    except (OSError, UnicodeError, ValueError):
        return None
    return entries or None


def _read_bgzf_index(filename: str) -> Optional[List[BgzfIndexEntry]]:
    """Read a fresh samtools-style ``.gzi`` index for a BGZF stream.

    ``.gzi`` stores compressed and uncompressed offsets for every BGZF block
    after the first.  The implicit origin is added here so callers can binary
    search every uncompressed FASTA byte offset, including offsets in block 0.
    Malformed, stale, or unrelated indexes are ignored and the FASTA reader can
    transparently fall back to sequential gzip decompression.
    """

    if not _is_bgzf(filename):
        return None
    index_path = f"{os.fspath(filename)}.gzi"
    if not os.path.isfile(index_path):
        return None
    try:
        if os.path.getmtime(index_path) < os.path.getmtime(filename):
            return None
        index_size = os.path.getsize(index_path)
        compressed_size = os.path.getsize(filename)
        with open(index_path, "rb") as index:
            count_bytes = index.read(8)
            if len(count_bytes) != 8:
                return None
            entry_count = struct.unpack("<Q", count_bytes)[0]
            if index_size != 8 + entry_count * 16:
                return None
            entries = [
                BgzfIndexEntry(*struct.unpack("<QQ", index.read(16)))
                for _ in range(entry_count)
            ]
    except (OSError, OverflowError, struct.error):
        return None

    if not entries or entries[0] != BgzfIndexEntry(0, 0):
        entries.insert(0, BgzfIndexEntry(0, 0))
    previous = entries[0]
    if previous != BgzfIndexEntry(0, 0):
        return None
    for entry in entries[1:]:
        if (
            entry.compressed_offset <= previous.compressed_offset
            or entry.uncompressed_offset <= previous.uncompressed_offset
            or entry.compressed_offset >= compressed_size
        ):
            return None
        previous = entry
    return entries


def supports_indexed_fasta_access(filename: str) -> bool:
    """Return whether individual records can be fetched without a full scan.

    Plain FASTA needs a fresh ``.fai``.  Compressed FASTA additionally needs
    to be BGZF with a usable ``.gzi``.  This predicate is intentionally
    conservative because the chromosome process pool must never make every
    worker decompress an ordinary gzip stream from the beginning.
    """

    if _read_fasta_index(filename) is None:
        return False
    return not _is_gzip(filename) or _read_bgzf_index(filename) is not None


def _iter_fasta_headers(filename: str) -> Iterator[str]:
    """Yield FASTA identifiers without assembling or validating sequences."""

    seen_ids = set()
    found_header = False
    with _open_fasta_text(filename) as fasta:
        for line_number, raw_line in enumerate(fasta, start=1):
            if raw_line.startswith(">"):
                description = raw_line[1:].strip()
                if not description:
                    raise ValueError(
                        f"Invalid FASTA {filename!s}: empty header at line {line_number}"
                    )
                sequence_id = description.split(maxsplit=1)[0]
                if sequence_id in seen_ids:
                    raise ValueError(
                        f"Invalid FASTA {filename!s}: duplicate sequence identifier "
                        f"{sequence_id!r} at line {line_number}"
                    )
                seen_ids.add(sequence_id)
                found_header = True
                yield sequence_id
            elif not found_header and raw_line.strip():
                raise ValueError(
                    f"Invalid FASTA {filename!s}: sequence data before the first "
                    f"header at line {line_number}"
                )

    if not found_header:
        raise ValueError(f"Invalid FASTA {filename!s}: no FASTA records found")


def _read_bgzf_range(
    filename: str,
    bgzf_index: Sequence[BgzfIndexEntry],
    start: int,
    size: int,
) -> bytes:
    """Read an uncompressed byte range from a BGZF stream."""

    if start < 0 or size < 0:
        raise ValueError("BGZF byte ranges must be non-negative")
    uncompressed_offsets = [entry.uncompressed_offset for entry in bgzf_index]
    block_number = bisect_right(uncompressed_offsets, start) - 1
    if block_number < 0:
        raise ValueError("BGZF index does not contain the start of the stream")
    block = bgzf_index[block_number]
    skip = start - block.uncompressed_offset

    with open(filename, "rb") as compressed:
        compressed.seek(block.compressed_offset)
        with gzip.GzipFile(fileobj=compressed, mode="rb") as uncompressed:
            if len(uncompressed.read(skip)) != skip:
                raise ValueError("BGZF index points past the end of the FASTA stream")
            return uncompressed.read(size)


def _fetch_indexed_region(
    filename: str,
    entry: FastaIndexEntry,
    start: int,
    end: int,
    *,
    bgzf_index: Optional[Sequence[BgzfIndexEntry]] = None,
) -> str:
    """Fetch one 1-based inclusive interval directly from an indexed FASTA."""

    if entry.length == 0 and start == 1 and end == 0:
        return ""
    if start < 1 or end < start or end > entry.length:
        raise ValueError(
            f"region {entry.name}:{start}-{end} is outside sequence length "
            f"{entry.length}"
        )
    start_index = start - 1
    end_index = end - 1
    start_byte = (
        entry.offset
        + (start_index // entry.line_bases) * entry.line_width
        + start_index % entry.line_bases
    )
    end_byte = (
        entry.offset
        + (end_index // entry.line_bases) * entry.line_width
        + end_index % entry.line_bases
    )
    byte_count = end_byte - start_byte + 1
    if _is_gzip(filename):
        if bgzf_index is None:
            bgzf_index = _read_bgzf_index(filename)
        if bgzf_index is None:
            raise ValueError(
                f"Compressed FASTA {filename!s} does not have a usable BGZF .gzi index"
            )
        raw_sequence = _read_bgzf_range(filename, bgzf_index, start_byte, byte_count)
    else:
        with open(filename, "rb") as fasta:
            fasta.seek(start_byte)
            raw_sequence = fasta.read(byte_count)

    sequence_bytes = raw_sequence.replace(b"\n", b"").replace(b"\r", b"")
    expected_length = end - start + 1
    if len(sequence_bytes) != expected_length:
        raise ValueError(
            f"FASTA index for {entry.name!r} returned {len(sequence_bytes)} bases; "
            f"expected {expected_length}"
        )
    if re.search(rb"[\t\v\f ]", sequence_bytes):
        raise ValueError(f"Invalid FASTA {filename!s}: whitespace within sequence data")
    return sequence_bytes.decode("ascii")


def _iter_selected_fasta_records(
    filename: str, regions, record_ids=None
) -> Iterator[Tuple[str, str]]:
    """Stream requested FASTA records and intervals in file order."""

    requested_ids = None if record_ids is None else list(record_ids)
    if requested_ids is not None and len(requested_ids) != len(set(requested_ids)):
        raise ValueError("FASTA record identifiers must be unique")
    requested_set = None if requested_ids is None else set(requested_ids)
    if requested_set == set():
        return

    sequence_id = None
    sequence_parts = []
    sequence_position = 0
    selected_region = None
    selection_complete = False
    collect_sequence = False
    seen_ids = set()
    yielded_ids = set()

    def selected_sequence():
        sequence = "".join(sequence_parts)
        if selected_region:
            _name, start, end = selected_region
            expected_length = end - start + 1
            if len(sequence) != expected_length:
                raise ValueError(
                    f"region {sequence_id}:{start}-{end} is outside sequence length "
                    f"{sequence_position}"
                )
        return sequence

    with _open_fasta_text(filename) as fasta:
        for line_number, raw_line in enumerate(fasta, start=1):
            if raw_line.startswith(">"):
                if sequence_id is not None and collect_sequence:
                    yield sequence_id, selected_sequence()
                    yielded_ids.add(sequence_id)
                    if requested_set is not None and yielded_ids == requested_set:
                        return

                description = raw_line[1:].strip()
                if not description:
                    raise ValueError(
                        f"Invalid FASTA {filename!s}: empty header at line {line_number}"
                    )
                sequence_id = description.split(maxsplit=1)[0]
                if sequence_id in seen_ids:
                    raise ValueError(
                        f"Invalid FASTA {filename!s}: duplicate sequence identifier "
                        f"{sequence_id!r} at line {line_number}"
                    )
                seen_ids.add(sequence_id)
                sequence_parts = []
                sequence_position = 0
                collect_sequence = requested_set is None or sequence_id in requested_set
                selected_region = (
                    regions.get(sequence_id) if collect_sequence and regions else None
                )
                selection_complete = False
                continue

            if sequence_id is None:
                if not raw_line.strip():
                    continue
                raise ValueError(
                    f"Invalid FASTA {filename!s}: sequence data before the first "
                    f"header at line {line_number}"
                )
            if selection_complete:
                continue

            line = raw_line.strip()
            if not line:
                continue
            if any(character.isspace() for character in line):
                raise ValueError(
                    f"Invalid FASTA {filename!s}: whitespace within sequence data "
                    f"at line {line_number}"
                )
            if not collect_sequence:
                continue

            line_start = sequence_position + 1
            line_end = sequence_position + len(line)
            if selected_region:
                _name, start, end = selected_region
                overlap_start = max(start, line_start)
                overlap_end = min(end, line_end)
                if overlap_start <= overlap_end:
                    local_start = overlap_start - line_start
                    local_end = overlap_end - line_start + 1
                    sequence_parts.append(line[local_start:local_end])
                sequence_position = line_end
                if sequence_position >= end:
                    selection_complete = True
                    if (
                        requested_set is not None
                        and (yielded_ids | {sequence_id}) == requested_set
                    ):
                        yield sequence_id, selected_sequence()
                        return
            else:
                sequence_parts.append(line)
                sequence_position = line_end

    if sequence_id is None:
        raise ValueError(f"Invalid FASTA {filename!s}: no FASTA records found")
    if collect_sequence:
        yield sequence_id, selected_sequence()


def _sequence_label(sequence_id: str, regions) -> str:
    if regions and sequence_id in regions:
        _name, start, end = regions[sequence_id]
        return f"{sequence_id}:{start}-{end}"
    return sequence_id


def iter_fasta_records(
    filename: str, regions=None, record_ids=None
) -> Iterator[Tuple[str, str, str]]:
    """Yield selected records as ``(identifier, sequence, display_label)``.

    A supplied ``record_ids`` sequence controls output order; with no explicit
    selection, records retain index or FASTA file order. Fresh ``.fai`` indexes
    provide record metadata for every FASTA. Plain FASTA and BGZF inputs with a
    valid ``.gzi`` are fetched directly, while ordinary gzip and unusable BGZF
    indexes retain the sequential decompression fallback.
    """

    requested_ids = None if record_ids is None else list(record_ids)
    if requested_ids is not None and len(requested_ids) != len(set(requested_ids)):
        raise ValueError("FASTA record identifiers must be unique")

    fasta_index = _read_fasta_index(filename)
    indexed_entries = (
        {entry.name: entry for entry in fasta_index} if fasta_index else None
    )
    if indexed_entries is not None:
        selected_ids = list(indexed_entries) if requested_ids is None else requested_ids
        missing_ids = [
            sequence_id
            for sequence_id in selected_ids
            if sequence_id not in indexed_entries
        ]
        if missing_ids:
            formatted = ", ".join(repr(sequence_id) for sequence_id in missing_ids)
            raise ValueError(f"FASTA record(s) not found: {formatted}")

        compressed = _is_gzip(filename)
        bgzf_index = _read_bgzf_index(filename) if compressed else None
        if not compressed or bgzf_index is not None:
            for sequence_id in selected_ids:
                entry = indexed_entries[sequence_id]
                if regions and sequence_id in regions:
                    _name, start, end = regions[sequence_id]
                else:
                    start, end = 1, entry.length
                sequence = _fetch_indexed_region(
                    filename,
                    entry,
                    start,
                    end,
                    bgzf_index=bgzf_index,
                )
                yield sequence_id, sequence, _sequence_label(sequence_id, regions)
            return
    else:
        selected_ids = requested_ids

    streamed_records = _iter_selected_fasta_records(
        filename, regions, record_ids=selected_ids
    )
    if selected_ids is None:
        for sequence_id, sequence in streamed_records:
            yield sequence_id, sequence, _sequence_label(sequence_id, regions)
        return

    # Streaming naturally discovers records in file order. Buffer only records
    # that precede the next explicitly requested identifier so the public API
    # can preserve caller order without requiring a separate header scan.
    pending = {}
    seen_selected = set()
    selected_position = 0
    for sequence_id, sequence in streamed_records:
        pending[sequence_id] = sequence
        seen_selected.add(sequence_id)
        while (
            selected_position < len(selected_ids)
            and selected_ids[selected_position] in pending
        ):
            selected_id = selected_ids[selected_position]
            yield (
                selected_id,
                pending.pop(selected_id),
                _sequence_label(selected_id, regions),
            )
            selected_position += 1

    if selected_position != len(selected_ids):
        missing_ids = [
            sequence_id
            for sequence_id in selected_ids
            if sequence_id not in seen_selected
        ]
        formatted = ", ".join(repr(sequence_id) for sequence_id in missing_ids)
        raise ValueError(f"FASTA record(s) not found: {formatted}")


def _iter_fasta_records(filename: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(identifier, sequence)`` records from a FASTA file.

    Identifiers follow the convention used by ``pysam.FastaFile``: only the
    first whitespace-delimited token after ``>`` is retained. Empty sequence
    records and blank lines are accepted, while malformed content, empty
    identifiers, and duplicate identifiers raise ``ValueError`` with the
    offending line number.
    """
    sequence_id = None
    sequence_parts = []
    seen_ids = set()

    with _open_fasta_text(filename) as fasta:
        for line_number, raw_line in enumerate(fasta, start=1):
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if sequence_id is not None:
                    yield sequence_id, "".join(sequence_parts)

                description = line[1:].strip()
                if not description:
                    raise ValueError(
                        f"Invalid FASTA {filename!s}: empty header at line {line_number}"
                    )

                sequence_id = description.split(maxsplit=1)[0]
                if sequence_id in seen_ids:
                    raise ValueError(
                        f"Invalid FASTA {filename!s}: duplicate sequence identifier "
                        f"{sequence_id!r} at line {line_number}"
                    )
                seen_ids.add(sequence_id)
                sequence_parts = []
                continue

            if sequence_id is None:
                raise ValueError(
                    f"Invalid FASTA {filename!s}: sequence data before the first "
                    f"header at line {line_number}"
                )
            if any(character.isspace() for character in line):
                raise ValueError(
                    f"Invalid FASTA {filename!s}: whitespace within sequence data "
                    f"at line {line_number}"
                )
            sequence_parts.append(line)

    if sequence_id is None:
        raise ValueError(f"Invalid FASTA {filename!s}: no FASTA records found")

    yield sequence_id, "".join(sequence_parts)


def extractRegion(seq_name):
    """Extract chromosome and region from seq_name.

    Supports:
    - Standard format: "chrY:50-3000"
    - Extended format: "HG002_chr13_MATERNAL:1-4000000:1000000-3000000"
      (keeps the last range).
    """
    # FASTA identifiers commonly contain periods, hyphens, and other
    # punctuation.  Treat everything before the first trailing coordinate
    # range as the identifier instead of restricting it to ``\w`` characters.
    # Repeated ranges are retained for compatibility with saved names such as
    # ``sample:1-4000000:1000000-3000000``; the final range wins.
    region_pattern = r"^(.+?)(?::\d+-\d+)+$"

    match = re.fullmatch(region_pattern, seq_name)
    if match:
        chrom = match.group(1)
        # Get all ranges (everything after the chrom)
        ranges = re.findall(r"(\d+)-(\d+)", seq_name)
        if ranges:
            # Take the last range
            lower_bound, upper_bound = map(int, ranges[-1])
            return chrom, lower_bound, upper_bound

    # No match
    return None


def _hash_sequence(
    seq: Sequence[str],
    k: int,
    fw_only: bool,
    ambiguous: bool = False,
) -> np.ma.MaskedArray:
    """Bulk-hash every genomic k-mer start with ntHash2.

    ntHash2 ordinarily skips windows containing non-ACGTU characters.  The
    native wrapper instead returns one hash for every start together with an
    ambiguity mask, allowing downstream window slicing to remain aligned with
    genomic coordinates.  Unless ``ambiguous`` is requested, those fallback
    hashes stay masked and therefore cannot enter a modimizer sketch.
    """
    if k <= 0:
        raise ValueError("k-mer size must be greater than zero")

    n = len(seq)
    total_kmers = max(n - k + 1, 0)

    if isinstance(seq, str):
        sequence = seq.upper().encode("ascii")
    else:
        sequence = bytes(seq).upper()

    hash_buffer, ambiguity_buffer = _nthash.hash_kmers(sequence, k, not fw_only)
    hashes = np.frombuffer(hash_buffer, dtype=np.uint64)
    if len(hashes) != total_kmers:
        raise RuntimeError(
            "ntHash2 returned an unexpected number of hashes: "
            f"expected {total_kmers}, received {len(hashes)}"
        )

    mask = np.ma.nomask
    if ambiguity_buffer:
        ambiguous_windows = np.frombuffer(ambiguity_buffer, dtype=np.uint8)
        if len(ambiguous_windows) != total_kmers:
            raise RuntimeError(
                "ntHash2 returned an ambiguity mask with an unexpected length: "
                f"expected {total_kmers}, received {len(ambiguous_windows)}"
            )
        if not ambiguous:
            mask = ambiguous_windows.astype(bool, copy=False)

    return np.ma.MaskedArray(hashes, mask=mask, copy=False)


def generateKmersFromFasta(
    seq: Sequence[str],
    k: int,
    quiet: bool,
    fw_only: bool,
    ambiguous: bool = False,
) -> Iterable[Optional[int]]:
    """Yield position-preserving ntHash2 values for a sequence.

    The public iterator remains compatible with existing callers.  Ambiguous
    windows are yielded as ``None`` unless ``ambiguous`` is enabled, while the
    FASTA reader below keeps the compact bulk NumPy representation in memory.
    ``quiet`` is retained for call compatibility; hashing no longer emits a
    progress display in either mode.
    """
    hashes = _hash_sequence(seq, k, fw_only, ambiguous)
    for kmer_hash in hashes:
        yield None if np.ma.is_masked(kmer_hash) else int(kmer_hash)


def isValidFasta(file_path):
    try:
        for _sequence_id, _sequence in _iter_fasta_records(file_path):
            pass
        return True
    except FileNotFoundError:
        print(f"Unable to find fasta {file_path}. Check filename and/or directory!\n")
        sys.exit(5)
    except ValueError as error:
        print(f"An error occurred: {error}")
        return False
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(6)


def extractFiles(folder_path):
    # Check to see at least one compressed numpy matrix, and one metadata pickle are included
    metadata = []
    matrices = []
    tmp = []
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)  # Full path to the file
        if filename.endswith(".pkl"):
            with open(file_path, "rb") as f:
                metadata = pickle.load(f)  # Append loaded data to the metadata list

    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if filename.endswith(".npz"):
            pattern = rf"_(\d+)\.npz"  # Using f-string to include the value of i in the regex pattern
            tmp2 = re.split(pattern, filename, maxsplit=1)
            ff = np.load(file_path, allow_pickle=True)
            tmp.append((tmp2[0], tmp2[1], ff))
    sorted_list = sorted(tmp, key=lambda x: (x[0], x[1]))

    unique_lists = {}

    # Iterate over the sorted list
    for item in sorted_list:
        key = item[0]  # Get the first element of the tuple
        if key in unique_lists:
            unique_lists[key].append(item)  # Append the item to the existing list
        else:
            unique_lists[key] = [item]  # Create a new list with the item

    # Convert dictionary values to lists
    result_lists = list(unique_lists.values())
    sorted_result_lists = [
        lst for title in metadata for lst in result_lists if lst[0][0] == title["title"]
    ]
    for unique_list in sorted_result_lists:
        matrices.append([])
        for val in unique_list:
            matrices[-1].append(val[-1]["data"])
    return matrices, metadata


def printProgressBar(
    iteration,
    total,
    prefix="",
    suffix="",
    decimals=1,
    length=100,
    fill="█",
    printEnd="\r",
):
    """Compatibility no-op retained after removal of progress displays."""

    return None


def readKmersFromFile(
    filename: str,
    ksize: int,
    quiet: bool,
    fw_only: bool,
    ambiguous: bool = False,
    regions=None,
    record_ids=None,
) -> List[np.ma.MaskedArray]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    all_kmers = []
    for seq_id, sequence, sequence_label in iter_fasta_records(
        filename, regions=regions, record_ids=record_ids
    ):
        if not quiet:
            print(f"Retrieving k-mers from {sequence_label}.... \n")
        if len(sequence) < ksize:
            if regions and seq_id in regions:
                _name, start, end = regions[seq_id]
                raise ValueError(
                    f"region {seq_id}:{start}-{end} is shorter than k-mer size "
                    f"{ksize}"
                )

        kmers_for_seq = _hash_sequence(sequence, ksize, fw_only, ambiguous)
        all_kmers.append(kmers_for_seq)
        if not quiet:
            print(f"\n{sequence_label} k-mers retrieved! \n")

    return all_kmers


def getInputHeaders(filename: str) -> List[str]:
    fasta_index = _read_fasta_index(filename)
    if fasta_index is not None:
        return [entry.name for entry in fasta_index]
    return list(_iter_fasta_headers(filename))


def getInputSeqLength(filename: str) -> List[int]:
    fasta_index = _read_fasta_index(filename)
    if fasta_index is not None:
        return [entry.length for entry in fasta_index]
    return [len(sequence) for _sequence_id, sequence in _iter_fasta_records(filename)]
