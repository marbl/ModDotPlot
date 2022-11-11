#TODO: Replace screed with a faster fasta farser
import screed 
import mmh3
import math

def hash_kmer(kmer):
    # calculate the reverse complement
    orientation = 1
    rc_kmer = screed.rc(kmer)

    # determine whether original k-mer or reverse complement is lesser
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer
        orientation = 0

    # calculate murmurhash using a hash seed of 42
    hashy = hash(canonical_kmer)
    if hashy < 0: hashy += 2**64

    # done
    return [hashy, orientation]

def build_kmers(sequence, ksize):
    kmers = []
    orientation = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer_hash = hash_kmer(sequence[i:i + ksize])
        kmers.append(kmer_hash[0])
        orientation.append(kmer_hash[1])

    return kmers


def read_kmers_from_file(filename, ksize):
    all_kmers = []
    for record in screed.open(filename):
        sequence = record.sequence

        kmers = build_kmers(sequence, ksize)
        all_kmers += kmers

    #print(all_kmers[1][0:10])
    return all_kmers