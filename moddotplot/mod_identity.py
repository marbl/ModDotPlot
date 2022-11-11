#!/usr/bin/env python3
from moddotplot.parse_fasta import read_kmers_from_file
from estimate_identity import binomial_distance, containment, get_mods, partition_windows
import sys
import argparse
import logging
import itertools
import pandas as pd
import os
import glob

def paired_bed_file(window_partitions, input_name):
    cols = [
        '#query_name', 
        'query_start',
        'query_end',
        'reference_name',
        'reference_start',
        'reference_end',
        'perID_by_events',
        'strand'
    ]
    bed = []
    for w in itertools.combinations_with_replacement(window_partitions, 2):
        if input_name:
            query_name = input_name
        else:
            query_name = "input_sequence"
        query_start = w[0].split("-")[0]
        query_end = w[0].split("-")[1]
        reference_start = w[1].split("-")[0]
        reference_end = w[1].split("-")[1]
        perID = binomial_distance(containment(set(window_partitions[w[0]]), set(window_partitions[w[1]])),21)
        if (perID*100 >= 80): 
            # TODO: Incorporate strand data into bedfile
            bed.append([query_name, query_start, query_end, query_name, reference_start, reference_end, perID*100, '+'])
    df = pd.DataFrame(bed, columns=cols)
    bedfile_output = input_name + ".bed"
    df.to_csv(bedfile_output, sep="\t")

def get_args_parse():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__,
    )
    mod_params = parser.add_argument_group(
        "Specific Patient / library highlight commands"
    )
    mod_params.add_argument(
        "-k", 
        default=21,
        help="Kmer size to use. Default: 21"
        )

    mod_params.add_argument(
        "-p",
        "--prefix",
        default=None,
        help="Prefix for input sequence. Defaults to name in the index file.",
    )

    mod_params.add_argument(
        "-o",
        "--output-folder",
        default=None,
        help="Folder for output tables and plots. Defaults to modplot_output.bed",
    )

    mod_params.add_argument(
        "-d",
        '--density',
        default=None,
        help="Mod value. Default is to take d as close to 1Mbp. Eg: d=64 for the Y chromosome.",
        type=int
    )
    mod_params.add_argument(
        "-nc",
        '--non-canonical',
        default=False,
        help="Only consider forward strand when computing kmers. Default = False"
    )

    mod_params.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input sequence. Must be indexed with samtools faidx!"
    )

    mod_params.add_argument(
        "--index",
        default=None,
        help="Input sequence index. Default: input.fai"
    )

    mod_params.add_argument(
        "-w",
        "--window-size",
        default=1000,
        type=int,
        help="Window size to use. Default = 1000"
    )
    args = parser.parse_args()

    return args

def main():
    args = get_args_parse()
    #TODO: Check fasta index file
    print("Parsing k-mers...")
    kmer_list = read_kmers_from_file(args.input, args.k)
    print("K-mers parsed!")
    logging.info("Computing modimizers...")
    mod_list = get_mods(kmer_list, args.density)
    print("Modimizers done!")
    print("Creating coordinates...")
    windows = partition_windows(mod_list)
    print("Coordinates done!")
    print("Computing identity...")
    paired_bed_file(windows, args.prefix)
    print("Bed file created!")

if __name__ == "__main__":
    main()
