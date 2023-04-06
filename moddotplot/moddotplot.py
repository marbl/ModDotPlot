#!/usr/bin/env python3
from moddotplot.parse_fasta import read_kmers_from_file, get_input_headers
from moddotplot.estimate_identity import binomial_distance, containment, get_mods, partition_windows
from moddotplot.interactive import run_dash
from moddotplot.interactive import get_custom_mods, create_custom_coordinates, get_matrix, run_dash 
from moddotplot.const import ASCII_ART
import argparse
import logging
import itertools
import pandas as pd
import os
import math
import glob
from moddotplot.static_plots import paired_bed_file

def get_args_parse():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Mod.Plot, A Rapid and Interactive Visualization of Tandem Repeats.",
    )

    req_params = parser.add_argument_group(
        "Required input"
    )

    req_params.add_argument(
        "-i",
        "--input",
        required=True,
        default=argparse.SUPPRESS,
        help="Path to input fasta file(s)",
        nargs="+",
        )

    dist_params = parser.add_argument_group(
        "Mod.Plot distance matrix commands"
    )

    dist_params.add_argument(
        "-k",
        "--kmer",
        default=21,
        help="k-mer length. Must be < 32"
    )

    dist_params.add_argument(
        "-d",
        '--density',
        help="Modimizer density value",
        default=None,
        type=int
    )

    dist_params.add_argument(
        "-r",
        "--resolution",
        default=1000,
        type=int,
        help="Dotplot resolution"
    )

    dist_params.add_argument(
        "-id",
        "--identity",
        default=80,
        type=int,
        help="Identity cutoff threshold"
    )

    dist_params.add_argument(
        "-o",
        "--output",
        default=argparse.SUPPRESS,
        help="Name for bed file and plots. Defaults to input fasta name",
    )

    dist_params.add_argument(
        "-nc",
        '--non-canonical',
        default=False,
        help="Only consider forward strand when computing k-mers"
    )

    plot_params = parser.add_argument_group(
        "Static plotting commands"
    )

    plot_params.add_argument(
        "--no-bed",
        action='store_true',
        help="Don't output bed file"
    )

    plot_params.add_argument(
        "--no-diag",
        action='store_true',
        help="Don't output diagonal plot"
    )

    plot_params.add_argument(
        "--no-pairwise",
        action='store_true',
        help="Don't output pairwise plot"
    )

    plot_params.add_argument(
        "--num-colors",
        default=11,
        type=int,
        help="Number of colors to map. Must be < 15"
    )

    interactive_params = parser.add_argument_group(
        "Interactive plotting commands"
    )

    interactive_params.add_argument(
        "--interactive",
        action='store_true',
        help="Launch a interactive Dash application on localhost."
    )

    interactive_params.add_argument(
        "--port",
        default="8050",
        type=int,
        help="Port number for launching interactive mode on localhost. Only used in interactive mode"
    )

    logging_params = parser.add_argument_group(
        "Logging options"
    )

    logging_params.add_argument(
        "-q",
        "--quiet",
        action='store_true',
        help="Supress "
    )

    args = parser.parse_args()

    return args

def main():
    print(ASCII_ART)
    args = get_args_parse()
    # Validate input sequence, display warning if multiple fasta files are shown:
    for i in args.input:
        headers = get_input_headers(i)
        kmer_list = read_kmers_from_file(i, args.kmer)

        for seq in range(len(headers)):
            if args.interactive:
                run_dash(kmer_list[seq], args.resolution)
            else:
                print(f"Computing modimizers for {headers[seq]}... \n")
                if not args.density:
                    # Get input sequence length
                    args.density = math.ceil(len(kmer_list[seq])/500000)
                    print(f"Density not provided. Using d = {args.density}. \n")
                mod_list = get_mods(kmer_list[seq], args.density)
                print("Modimizers done! \n")
                print("Creating coordinates...\n")
                windows = partition_windows(mod_list, args.resolution)
                print("Coordinates done! \n")
                print("Computing identity... \n")
                paired_bed_file(windows, headers[seq], args.identity, args.density, None)
                print("Bed file created!")

if __name__ == "__main__":
    main()
