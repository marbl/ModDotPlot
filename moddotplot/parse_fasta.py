from typing import Generator, List
import pysam
import sys


def custom_hash_fn(h):
    h ^= h >> 33
    h *= 0xFF51AFD7ED558CCD
    h ^= h >> 33
    h *= 0xC4CEB9FE1A85EC53
    h ^= h >> 33
    return h


def is_valid_fasta(file_path):
    try:
        with open(file_path, "r") as file:
            in_sequence = False
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    if in_sequence:
                        print("Fasta formatting error: '>' found within sequence")
                        sys.exit(2)
                    in_sequence = True
                elif in_sequence and not line:
                    # Empty line encountered after the sequence header
                    in_sequence = False  # Exit the sequence mode but continue checking
                elif in_sequence:
                    # Check if the line contains valid sequence characters (modify this condition as needed)
                    # For example, you could add your sequence validation logic here
                    pass
            return in_sequence
    except FileNotFoundError:
        print("Unable to find fasta file. Check filename and/or directory!")
        sys.exit(5)
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        sys.exit(6)


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
    percent = f"{100 * (iteration / total):.{decimals}f}"
    filledLength = int(length * iteration // total)
    bar = [fill] * filledLength + ["-"] * (length - filledLength)
    bar_str = "".join(bar)
    print(f"\r{prefix} |{bar_str}| {percent}% {suffix}", end=printEnd)
    if iteration == total:
        print()


def generate_kmers(sequence: str, k: int) -> Generator[int, None, None]:
    """
    Given a DNA sequence and a value k, generate all k-mers from the sequence.

    Args:
    - sequence (str): a DNA sequence
    - k (int): the length of the k-mers to generate

    Yields:
    - int: the hash value of each k-mer and its reverse complement

    """

    # Calculate the length of the sequence
    n = len(sequence)
    moddh = round(n / 77)

    # Progress bar updates
    printProgressBar(0, n, prefix="Progress:", suffix="Complete", length=40)

    # Calculate a mask to use for masking off any bits outside of the k-mer range
    mask = (1 << (3 * k)) - 1

    # Initialize the first k-mer and its reverse complement
    kmer = 0
    rc_kmer = 0
    for i in range(k):
        # Add the next base to the k-mer and its reverse complement
        kmer = (kmer << 3) | encode_base(sequence[i])
        rc_kmer = (rc_kmer >> 3) | (encode_base(sequence[i], True) ^ 0b111) << (
            3 * (k - 1)
        )

    # Yield the hash value of the first k-mer and its reverse complement
    yield custom_hash_fn(rc_kmer) if custom_hash_fn(rc_kmer) < custom_hash_fn(
        kmer
    ) else custom_hash_fn(kmer)

    # Generate the remaining k-mers
    for i in range(k, n):
        # Update the k-mer and its reverse complement by adding the next base and dropping the leftmost base
        kmer = ((kmer << 3) & mask) | encode_base(sequence[i])
        rc_kmer = (rc_kmer >> 3) | (encode_base(sequence[i], True) ^ 0b111) << (
            3 * (k - 1)
        )

        # Print progress only if
        if i % moddh == 0:
            printProgressBar(i, n, prefix="Progress:", suffix="Completed", length=40)

        # Yield the hash value of the current k-mer and its reverse complement
        yield custom_hash_fn(rc_kmer) if custom_hash_fn(rc_kmer) < custom_hash_fn(
            kmer
        ) else custom_hash_fn(kmer)


def encode_base(base: str, reverse_complement: bool = False) -> int:
    """Encode a DNA base to a numeric representation.

    If `reverse_complement` is True, return the complement of the base instead.
    """
    if not reverse_complement:
        if base == "A":
            return 0b000
        elif base == "C":
            return 0b001
        elif base == "G":
            return 0b010
        elif base == "T":
            return 0b011
        elif base == "N":
            return 0b100
        else:
            return 0b100
    else:
        if base == "T":
            return 0b111
        elif base == "G":
            return 0b110
        elif base == "C":
            return 0b101
        elif base == "A":
            return 0b100
        elif base == "N":
            return 0b011
        else:
            return 0b100


def report_all_kmers(sequence: str, k: int) -> List[int]:
    """
    Given a DNA sequence and an integer k, returns a list of all k-mers found in the sequence.
    """
    all_mers = []
    # Generate all k-mers of length k using a bitmask
    kmers = generate_kmers(sequence, k)

    # Iterate through the k-mers and append them to the list
    for kmer in kmers:
        all_mers.append(kmer)

    return all_mers


def read_kmers_from_file(filename: str, ksize: int) -> List[List[int]]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    all_kmers = []
    seq = pysam.FastaFile(filename)

    for seq_id in seq.references:
        print(f"Retrieving k-mers from {seq_id}.... \n")
        all_kmers.append(report_all_kmers(seq.fetch(seq_id), ksize))
        print(f"{seq_id} k-mers retrieved! \n")

    return all_kmers


def get_input_headers(filename: str) -> List[str]:
    header_list = []
    seq = pysam.FastaFile(filename)
    for seq_id in seq.references:
        header_list.append(seq_id)

    return header_list


def get_input_seq_length(filename: str) -> List[int]:
    length_list = []
    seq = pysam.FastaFile(filename)
    return seq.lengths
