from typing import Iterable, List, Sequence
import pysam
import sys
import mmh3

tab_b = bytes.maketrans(b"ACTG", b"TGAC")

def generate_kmers_from_fasta(seq: Sequence[str], k: int) -> Iterable[int]:
    n = len(seq)
    progress_thresholds = round(n / 77)
    printProgressBar(0, n - k + 1, prefix='Progress:', suffix='Complete', length=40)

    for i in range(n - k + 1):
        if i % progress_thresholds == 0:
            printProgressBar(i, n - k + 1, prefix='Progress:', suffix='Complete', length=40)
        if i == n - k:
            printProgressBar(n - k + 1, n - k + 1, prefix='Progress:', suffix='Completed', length=40)
            print('\n')

        kmer = seq[i:i + k]
        fh = mmh3.hash(kmer)
        
        # Calculate reverse complement hash directly without the need for translation
        rc = mmh3.hash(kmer[::-1].translate(tab_b))
        
        yield fh if fh < rc else rc

def is_valid_fasta(file_path):
    try:
        with open(file_path, "r") as file:
            in_sequence = False
            for line in file:
                line = line.strip()
                if line.startswith(">"):
                    in_sequence = True
                elif in_sequence and not line:
                    # Empty line encountered after the sequence header
                    in_sequence = False
                elif in_sequence:
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


def read_kmers_from_file(filename: str, ksize: int) -> List[List[int]]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    all_kmers = []
    seq = pysam.FastaFile(filename)

    for seq_id in seq.references:
        print(f"Retrieving k-mers from {seq_id}.... \n")
        kmers_for_seq = []
        for kmer_hash in generate_kmers_from_fasta(seq.fetch(seq_id), ksize):
            kmers_for_seq.append(kmer_hash)
        all_kmers.append(kmers_for_seq)
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
