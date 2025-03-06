from enum import unique
from typing import Iterable, List, Sequence
import pysam
import sys
import mmh3
import os
import pickle
import re
import numpy as np
import gzip

tab_b = bytes.maketrans(b"ACTG", b"TGAC")


def generateKmersFromFasta(seq: Sequence[str], k: int, quiet: bool) -> Iterable[int]:
    n = len(seq)
    if not quiet:
        progress_thresholds = round(n / 77)
        printProgressBar(0, n - k + 1, prefix="Progress:", suffix="Complete", length=40)

    for i in range(n - k + 1):
        if not quiet:
            if i % progress_thresholds == 0:
                printProgressBar(
                    i, n - k + 1, prefix="Progress:", suffix="Complete", length=40
                )
            if i == n - k:
                printProgressBar(
                    n - k + 1,
                    n - k + 1,
                    prefix="Progress:",
                    suffix="Completed",
                    length=40,
                )
        # Remove case sensitivity
        kmer = seq[i : i + k].upper()
        fh = mmh3.hash(kmer)

        # Calculate reverse complement hash directly without the need for translation
        rc = mmh3.hash(kmer[::-1].translate(tab_b))

        yield fh if fh < rc else rc


def isValidFasta(file_path):
    try:
        open_func = gzip.open if file_path.endswith(".gz") else open
        with open_func(file_path, "rt") as file:
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
        print(f"Unable to find fasta {file_path}. Check filename and/or directory!\n")
        sys.exit(5)
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
    percent = f"{100 * (iteration / total):.{decimals}f}"
    filledLength = int(length * iteration // total)
    bar = [fill] * filledLength + ["-"] * (length - filledLength)
    bar_str = "".join(bar)
    print(f"\r{prefix} |{bar_str}| {percent}% {suffix}", end=printEnd)
    if iteration == total:
        print()


def readKmersFromFile(filename: str, ksize: int, quiet: bool) -> List[List[int]]:
    """
    Given a filename and an integer k, returns a list of all k-mers found in the sequences in the file.
    """
    all_kmers = []
    seq = pysam.FastaFile(filename)

    for seq_id in seq.references:
        print(f"Retrieving k-mers from {seq_id}.... \n")
        kmers_for_seq = []
        for kmer_hash in generateKmersFromFasta(seq.fetch(seq_id), ksize, quiet):
            kmers_for_seq.append(kmer_hash)
        all_kmers.append(kmers_for_seq)
        print(f"\n{seq_id} k-mers retrieved! \n")

    return all_kmers


def getInputHeaders(filename: str) -> List[str]:
    header_list = []
    try:
        seq = pysam.FastaFile(filename)
    except OSError:
        seq = None
    for seq_id in seq.references:
        header_list.append(seq_id)

    return header_list


def getInputSeqLength(filename: str) -> List[int]:
    seq = pysam.FastaFile(filename)
    return seq.lengths
