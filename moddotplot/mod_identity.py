#!/usr/bin/env python3
from moddotplot.parse_fasta import read_kmers_from_file, get_input_headers
from moddotplot.estimate_identity import binomial_distance, containment, get_mods, partition_windows
#from ...const import ASCII_ART
import sys
import argparse
import logging
import itertools
import pandas as pd
import os
import glob

def paired_bed_file(window_partitions, input_name, id_threshold):

    assert id_threshold > 50 and id_threshold < 100

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
        if (perID*100 >= id_threshold): 
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
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Mod.Plot, A Rapid and Interactive Visualization of Tandem Repeats.",
    )

    parser.add_argument(
        "Input",
        help="Path to input fasta file"
        )

    mod_params = parser.add_argument_group(
        "Mod.Plot commands"
    )

    mod_params.add_argument(
        "-k",
        "--kmer",
        default=21,
        help="k-mer length"
    )

    mod_params.add_argument(
        "-d",
        '--density',
        help="Modimizer density value (default: d = 2/Mbp)",
        default=argparse.SUPPRESS,
        type=int
    )

    mod_params.add_argument(
        "-r",
        "--resolution",
        default=1000,
        type=int,
        help="Dotplot resolution. Recommended > 500 for high quality dot plots"
    )

    mod_params.add_argument(
        "-id",
        "--identity",
        default=80,
        type=int,
        help="Identity cutoff threshold"
    )

    mod_params.add_argument(
        "-o",
        "--dir",
        default=argparse.SUPPRESS,
        help="Folder for output tables and plots. Defaults to current directory",
    )

    mod_params.add_argument(
        "-nc",
        '--non-canonical',
        default=False,
        help="Only consider forward strand when computing k-mers"
    )

    mod_params.add_argument(
        "-i",
        "--interactive",
        default=False,
        help="Launch a interactive Dash application on localhost."
    )

    mod_params.add_argument(
        "--port",
        default="8050",
        type=int,
        help="Port number for launching interactive mode on localhost. Only used in interactive mode"
    )

    args = parser.parse_args()

    return args

def main():
    args = get_args_parse()
    # Validate input sequence, display warning if multiple fasta files are shown:
    print("Parsing k-mers...")
    headers = get_input_headers(args.Input)
    kmer_list = read_kmers_from_file(args.Input, args.kmer)
    print("K-mers parsed!")
    for seq in range(len(headers)):
        print(f" Computing modimizers for {headers[seq]}...")
        #TODO: Function to auto-determine density
        mod_list = get_mods(kmer_list[seq], args.density)
        print("Modimizers done!")
        print("Creating coordinates...")
        windows = partition_windows(mod_list, args.resolution)
        print("Coordinates done!")
        print("Computing identity...")
        paired_bed_file(windows, headers[seq], args.identity)
        print("Bed file created!")

if __name__ == "__main__":
    main()
