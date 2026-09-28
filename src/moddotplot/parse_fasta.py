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
import sys
import os
import pickle
import re
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


def _is_gzip(filename: str) -> bool:
    with open(filename, "rb") as probe:
        return probe.read(2) == b"\x1f\x8b"


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

    if _is_gzip(filename):
        return None
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


def _fetch_indexed_region(
    filename: str, entry: FastaIndexEntry, start: int, end: int
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
    with open(filename, "rb") as fasta:
        fasta.seek(start_byte)
        raw_sequence = fasta.read(end_byte - start_byte + 1)

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
    filename: str, regions, single_record: bool
) -> Iterator[Tuple[str, str]]:
    """Stream selected intervals, stopping early for a one-record FASTA."""

    sequence_id = None
    sequence_parts = []
    sequence_position = 0
    selected_region = None
    selection_complete = False
    seen_ids = set()

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
                if sequence_id is not None:
                    yield sequence_id, selected_sequence()

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
                selected_region = regions.get(sequence_id) if regions else None
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
                    if single_record:
                        yield sequence_id, selected_sequence()
                        return
            else:
                sequence_parts.append(line)
                sequence_position = line_end

    if sequence_id is None:
        raise ValueError(f"Invalid FASTA {filename!s}: no FASTA records found")
    yield sequence_id, selected_sequence()


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
    """
    total_kmers = max(len(seq) - k + 1, 0)
    if not quiet:
        printProgressBar(
            0, total_kmers, prefix="Progress:", suffix="Complete", length=40
        )

    hashes = _hash_sequence(seq, k, fw_only, ambiguous)
    progress_threshold = max(round(total_kmers / 77), 1)
    for index, kmer_hash in enumerate(hashes):
        if not quiet and index % progress_threshold == 0:
            printProgressBar(
                index,
                total_kmers,
                prefix="Progress:",
                suffix="Complete",
                length=40,
            )

        yield None if np.ma.is_masked(kmer_hash) else int(kmer_hash)

        if not quiet and index == total_kmers - 1:
            printProgressBar(
                total_kmers,
                total_kmers,
                prefix="Progress:",
                suffix="Completed",
                length=40,
            )


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
    if total <= 0:
        percent = f"{100:.{decimals}f}"
        filledLength = length
    else:
        percent = f"{100 * (iteration / total):.{decimals}f}"
        filledLength = int(length * iteration // total)
    bar = [fill] * filledLength + ["-"] * (length - filledLength)
    bar_str = "".join(bar)
    print(f"\r{prefix} |{bar_str}| {percent}% {suffix}", end=printEnd)
    if iteration == total:
        print()


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
    record_ids = (
        list(record_ids) if record_ids is not None else getInputHeaders(filename)
    )
    fasta_index = _read_fasta_index(filename)
    indexed_entries = (
        {entry.name: entry for entry in fasta_index}
        if fasta_index and [entry.name for entry in fasta_index] == record_ids
        else None
    )

    if indexed_entries is not None:

        def indexed_records():
            for seq_id in record_ids:
                entry = indexed_entries[seq_id]
                if regions and seq_id in regions:
                    _name, start, end = regions[seq_id]
                else:
                    start, end = 1, entry.length
                sequence_label = (
                    f"{seq_id}:{start}-{end}"
                    if regions and seq_id in regions
                    else seq_id
                )
                print(f"Retrieving k-mers from {sequence_label}.... \n")
                sequence = _fetch_indexed_region(filename, entry, start, end)
                yield seq_id, sequence, sequence_label

        selected_records = indexed_records()
    else:
        selected_records = (
            (
                seq_id,
                sequence,
                (
                    f"{seq_id}:{regions[seq_id][1]}-{regions[seq_id][2]}"
                    if regions and seq_id in regions
                    else seq_id
                ),
            )
            for seq_id, sequence in _iter_selected_fasta_records(
                filename, regions, single_record=len(record_ids) == 1
            )
        )

    for seq_id, sequence, sequence_label in selected_records:
        if indexed_entries is None:
            print(f"Retrieving k-mers from {sequence_label}.... \n")
        if len(sequence) < ksize:
            if regions and seq_id in regions:
                _name, start, end = regions[seq_id]
                raise ValueError(
                    f"region {seq_id}:{start}-{end} is shorter than k-mer size "
                    f"{ksize}"
                )

        total_kmers = max(len(sequence) - ksize + 1, 0)
        if not quiet:
            printProgressBar(
                0, total_kmers, prefix="Progress:", suffix="Complete", length=40
            )
        kmers_for_seq = _hash_sequence(sequence, ksize, fw_only, ambiguous)
        if not quiet:
            printProgressBar(
                total_kmers,
                total_kmers,
                prefix="Progress:",
                suffix="Completed",
                length=40,
            )
        all_kmers.append(kmers_for_seq)
        print(f"\n{sequence_label} k-mers retrieved! \n")

    return all_kmers


def getInputHeaders(filename: str) -> List[str]:
    fasta_index = _read_fasta_index(filename)
    if fasta_index is not None:
        return [entry.name for entry in fasta_index]
    return list(_iter_fasta_headers(filename))


def getInputSeqLength(filename: str) -> List[int]:
    return [len(sequence) for _sequence_id, sequence in _iter_fasta_records(filename)]
