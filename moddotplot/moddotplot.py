#!/usr/bin/env python3
from bz2 import compress
import sys
from moddotplot.parse_fasta import (
    read_kmers_from_file,
    get_input_headers,
    is_valid_fasta,
)

from moddotplot.estimate_identity import (
    get_mods,
    convert_set,
    convert_set_neighbors,
    partition_evenly_spaced_modimizers,
    partition_pairwise_modimizers_different_size,
    convert_to_modimizers,
    self_containment_matrix,
    pairwise_containment_matrix,
    next_power_of_two,
    partition_overlaps
)
from moddotplot.interactive import run_dash
from moddotplot.const import ASCII_ART, VERSION

import argparse
import math
from moddotplot.static_plots import (
    paired_bed_file,
    paired_bed_file_a_vs_b,
    read_df_from_file,
    create_plots
)
import itertools
import json
import numpy as np

def get_parser():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="ModDotPlot: Visualization of Complex Repeat Structures",
    )
    subparsers = parser.add_subparsers(dest='command', help='Choose mode: interactive or static')
    interactive_parser = subparsers.add_parser('interactive', help='Interactive mode commands')
    static_parser = subparsers.add_parser('static', help='Static mode commands')

    # -----------INTERACTIVE MODE SUBCOMMANDS-----------
    interactive_input_group = interactive_parser.add_mutually_exclusive_group(required=True)

    interactive_input_group.add_argument(
        "-f",
        "--fasta",
        default=argparse.SUPPRESS,
        help="Path to input fasta file(s).",
        nargs="+",
    )

    interactive_input_group.add_argument(
        "-l",
        "--load",
        default=None,
        type=str,
        help="Load previously computed hierarchical matrices.",
    )

    # Add a mutually exclusive group for compare and compare only.
    compare_group = interactive_parser.add_mutually_exclusive_group(required=False)
    window_size_group = interactive_parser.add_mutually_exclusive_group(required=False)

    interactive_parser.add_argument(
        "-k", "--kmer", default=21, type=int, help="k-mer length."
    )

    interactive_parser.add_argument(
        "-s",
        "--sparsity",
        help="Modimizer sparsity value at base value. Higher value will reduce the number of modimizers, but will increase performance. Default set to 1 per Mbp of sequence, rounded to the nearest power of two.",
        default=None,
        type=int,
    )

    window_size_group.add_argument(
        "-r",
        "--resolution",
        default=1000,
        type=int,
        help="Dotplot resolution, or the number of intervals to compare against.",
    )

    window_size_group.add_argument(
        "-w",
        "--window",
        default=None,
        type=int,
        help="Window size, or the length in genomic coordinates of each interval. Default is set to (genome length)/(resolution)",
    )

    interactive_parser.add_argument(
        "-id",
        "--identity",
        default=86,
        type=int,
        help="Identity cutoff threshold.",
    )

    interactive_parser.add_argument(
        "-d",
        "--delta",
        default=0.5,
        type=float,
        help="Fraction of neighboring partition to include in identity estimation. Must be between 0 and 1, use > 0.5 is not recommended.",
    )

    interactive_parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Directory name for saving matrices and coordinate logs. Defaults to working directory.",
    )

    compare_group.add_argument(
        "--compare",
        action="store_true",
        help="Create a dotplot with two different sequences (in addition to self-identity plots).",
    )

    compare_group.add_argument(
        "--compare-only",
        action="store_true",
        help="Create a dotplot with two different sequences (skips self-identity plots).",
    )

    #TODO: Change to min:resolution
    interactive_parser.add_argument(
        "--layers",
        default=3,
        type=int,
        help="Number of matrix hierarchy layers to compute when preparing interactive mode.",
    )

    interactive_parser.add_argument(
        "--save",
        action="store_true",
        help="Save hierarchical matrices to file. Only used in interactive mode.",
    )

    interactive_parser.add_argument(
        "--port",
        default="8050",
        type=int,
        help="Port number for launching interactive mode on localhost. Only used in interactive mode.",
    )

    # -----------STATIC MODE SUBCOMMANDS-----------
    static_input_group = static_parser.add_mutually_exclusive_group(required=True)
    static_input_group.add_argument(
        "-c",
        "--config",
        default=None,
        type=str,
        help="Config file to use. Takes precedence over any other competing command line arguments.",
    )

    static_input_group.add_argument(
        "-b",
        "--bed",
        default=argparse.SUPPRESS,
        help="Path to input bed file(s). Exclusively used in static mode.",
        nargs="+",
    )

    static_input_group.add_argument(
        "-f",
        "--fasta",
        default=argparse.SUPPRESS,
        help="Path to input fasta file(s).",
        nargs="+",
    )

    # Add a mutually exclusive group for compare and compare only.
    static_compare_group = static_parser.add_mutually_exclusive_group(required=False)
    static_window_size_group = static_parser.add_mutually_exclusive_group(required=False)

    static_parser.add_argument(
        "-k", "--kmer", default=21, type=int, help="k-mer length."
    )

    static_parser.add_argument(
        "-s",
        "--sparsity",
        help="Modimizer sparsity value. Lower value will result in a more accurate plot, but will reduce runtime and memory performance. Default set to 1 per Mbp of sequence, rounded up.",
        default=None,
        type=int,
    )

    static_window_size_group.add_argument(
        "-r",
        "--resolution",
        default=1000,
        type=int,
        help="Dotplot resolution, or the number of intervals to compare against.",
    )

    static_window_size_group.add_argument(
        "-w",
        "--window",
        default=None,
        type=int,
        help="Window size, or the length in genomic coordinates of each interval. Default is set to (genome length)/(resolution)",
    )

    static_parser.add_argument(
        "-id",
        "--identity",
        default=86,
        type=int,
        help="Identity cutoff threshold.",
    )

    static_parser.add_argument(
        "-d",
        "--delta",
        default=0.5,
        type=float,
        help="Fraction of neighboring partition to include in identity estimation. Must be between 0 and 1, use > 0.5 is not recommended.",
    )

    static_parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Directory name for saving bed files and plots. Defaults to working directory.",
    )

    static_compare_group.add_argument(
        "--compare",
        action="store_true",
        help="Create a dotplot with two different sequences (in addition to self-identity plots).",
    )

    static_compare_group.add_argument(
        "--compare-only",
        action="store_true",
        help="Create a dotplot with two different sequences (skips self-identity plots).",
    )

    static_parser.add_argument(
        "--no-bed", action="store_true", help="Skip output of bed file."
    )

    static_parser.add_argument(
        "--no-plot", action="store_true", help="Skip output of plots."
    )

    static_parser.add_argument(
        "--no-hist", action="store_true", help="Skip output of histogram color legend."
    )

    static_parser.add_argument(
        "--width", default=9, type=float, nargs="+", help="Plot width."
    )

    static_parser.add_argument(
        "--height", default=5, type=float, nargs="+", help="Plot height, for self-identity plots. "
    )

    static_parser.add_argument(
        "--xaxis",
        default=None,
        type=float,
        nargs="+",
        help="Change x axis for self identity plots. Default is length of the sequence, in mbp.",
    )

    static_parser.add_argument(
        "--dpi", default=600, type=int, help="Plot dpi."
    )

    # TODO: Create list of accepted colors.
    static_parser.add_argument(
        "--palette",
        default="Spectral_11",
        help="Select color palette. See RColorBrewer for list of accepted palettes. Will default to Spectral_11 if not used.",
        type=str,
    )

    static_parser.add_argument(
        "--palette-orientation",
        default="+",
        choices=["+", "-"],
        help="Color palette orientation. + for forward, - for reverse.",
        type=str,
    )

    static_parser.add_argument(
        "--colors",
        default=None,
        nargs="+",
        help="Use a custom color palette, entered in either hexcode or rgb format.",
    )

    static_parser.add_argument(
        "--breakpoints",
        default=None,
        nargs="+",
        help="Introduce custom color thresholds. Must be between identity threshold and 100.",
    )

    static_parser.add_argument(
        "--bin-freq",
        action="store_true",
        help="By default, histograms are evenly spaced based on the number of colors and the identity threshold. Select this argument to bin based on the frequency of observed identity values.",
    )

    # TODO: Implement static mode logging options

    return parser


def main():
    print(ASCII_ART)
    print(f"v{VERSION} \n")
    args = get_parser().parse_args()

    # -----------MUTUALLY EXCLUSIVE: INTERACTIVE OR STATIC MODE-----------
    if args.command == "interactive":
        print(f"Running ModDotPlot in interactive mode\n")
        # -----------LOAD MATRICES FOR INTERACTIVE MODE-----------
        if hasattr(args, 'load') and args.load:
            print(f"Loading matrix hierarchy from {args.load}... \n")
            #TODO: Throw error if invalid
            npz_file = np.load(args.load)
            file_list = list(npz_file)
            print(npz_file[file_list[-1]])
            # TODO: Get names & sequence ids from compressed numpy file
            '''run_dash(
                    matrix_metadata[0],
                    matrix_metadata[1],
                    matrix_metadata[3],
                    matrix_metadata[4],
                    args.resolution,
                    args.sparsity,
                    args.kmer,
                    args.identity,
                    args.port,
                    args.palette,
                    args.palette_orientation,
                    args.delta,
                    image_pyramid,
                )'''
            sys.exit(0)
    elif args.command == "static":
        print(f"Running ModDotPlot in static mode\n")
        # -----------CONFIG PARSING-----------
        #TODO: Change to yml file, add radme to config folder
        if args.config:
            with open(args.config, "r") as f:
                config = json.load(f)
                #TODO: Make sure these are mutually exclusive
                args.fasta = config.get("fasta")
                args.bed = config.get("bed")

                # Distance matrix commands
                args.kmer = config.get("kmer", args.kmer)
                args.sparsity = config.get("sparsity", args.sparsity)
                args.resolution = config.get("resolution", args.resolution)
                args.window = config.get("window", args.window)
                args.identity = config.get("identity", args.identity)
                args.delta = config.get("delta", args.delta)
                args.output_dir = config.get("output_dir", args.output_dir)
                args.compare = config.get("compare", args.compare)
                args.compare_only = config.get("compare_only", args.compare_only)

                # Interactive plotting commands
                args.layers = config.get("layers", args.layers)
                args.save = config.get("save", args.save)
                args.load = config.get("load", args.load)
                args.port = config.get("port", args.port)

                # Static plotting commands
                args.static = config.get("static", args.static)
                args.no_bed = config.get("no_bed", args.no_bed)
                args.no_plot = config.get("no_plot", args.no_plot)
                args.no_hist = config.get("no_hist", args.no_hist)
                args.width = config.get("width", args.width)
                args.height = config.get("height", args.height)
                args.xaxis = config.get("xaxis", args.xaxis)
                args.dpi = config.get("dpi", args.dpi)
                args.palette = config.get("palette", args.palette)
                args.palette_orientation = config.get(
                    "palette_orientation", args.palette_orientation
                )
                args.colors = config.get("color", args.colors)
                args.breakpoints = config.get("breakpoints", args.breakpoints)
                args.bin_freq = config.get("bin_freq", args.bin_freq)

            # TODO: Include logging options here

        # -----------INPUT COMMAND VALIDATION-----------
        # TODO: More tests!
        # Validation for breakpoints command
        if args.breakpoints:
            # Check that start value for breakpoints = identity threshold value
            if float(args.breakpoints[0]) != float(args.identity):
                print(
                    f"Identity threshold is {args.identity}, but starting breakpoint is {args.breakpoints[0]}! \n"
                )
                print(
                    f"Please modify identity threshold using --identity {args.breakpoints[0]}.\n"
                )
                #TODO: Chronicle exit codes
                # Exit code 2: breakpoint value != idnetity threshold
                sys.exit(2)

        # -----------BEDFILE INPUT FOR STATIC MODE-----------
        if hasattr(args, 'bed') and args.bed:
            try:
                for bed in args.bed:
                    # If args.bed is provided as input, run static mode directly from the bed file. Skip counting input k-mers.
                    df = read_df_from_file(bed)
                    unique_query_names = df['#query_name'].unique()
                    unique_reference_names = df['reference_name'].unique()
                    if len(unique_query_names) == 1 and len(unique_reference_names) == 1:
                        print(f"Input bed file {bed} read successfully! Creating plots... \n")
                        create_plots(None, args.output_dir if args.output_dir else ".", unique_query_names[0], unique_reference_names[0], args.palette, args.palette_orientation, args.no_hist, 9, 5, args.dpi, args.bin_freq, None, args.colors, args.breakpoints, df)
                    else:
                        print(f"Unresolvable bed file\n")
                        # Exit code 6: Weird columns in bed file...
                        sys.exit(6)
            except Exception as e:
                # Exit code 7: Error getting info from bed file: 
                #TODO: Change to logs
                print(f"Error in bed file: {e}")
                sys.exit(7)

    # -----------INPUT SEQUENCE VALIDATION-----------
    seq_list = []
    for i in args.fasta:
        is_valid_fasta(i)
        headers = get_input_headers(i)
        if len(headers) > 1:
            print(f"File {i} contains multiple fasta entries: \n")
            counter = 1
            for j in headers:
                counter += 1
        for j in headers:
            seq_list.append(j)

    # -----------LOAD SEQUENCES INTO MEMORY-----------
    kmer_list = []
    for i in args.fasta:
        kmer_list.append(read_kmers_from_file(i, args.kmer))
    k_list = [item for sublist in kmer_list for item in sublist]

    # -----------SET SPARSITY VALUE-----------
    if not args.sparsity:
        max_length = max(len(x) for x in k_list)
        args.sparsity = -(-max_length // (500000 * 2))

    # -----------LAUNCH INTERACTIVE MODE-----------
    if args.command == "interactive":
        # Single sequence, can set window length immediately
        if len(seq_list) == 1:
            if args.resolution:
                args.window = math.ceil( (len(k_list[0]) + args.kmer - 1) / args.resolution)
        # Set warning if >2 sequences detected
        elif len(seq_list) > 2:
            print(f"{len(seq_list)} sequences were detected, however interactive mode can only load two sequences at a time.\n")
            print(f"Interactive mode will proceed with {seq_list[0]} and {seq_list[1]}\n")
        # Set warning if <1 sequence detected
        elif len(seq_list) < 1:
            print(f"Error: No sequences detected!")
            sys.exit(5)
        # Set window size based on largest sequence out of 2
        else:
            if args.resolution:
                max_size = max(len(k_list[0]), len(k_list[1]))
                args.window = math.ceil(max_size/args.resolution)

        # Set sparsity to be the closest power of 2
        args.sparsity = next_power_of_two(args.sparsity)
        print(f"Setting top layer sparsity = {args.sparsity}. \n")

        # TODO: assert that args layers is a reasonable number
        image_pyramid_1 = []
        image_pyramid_2 = []
        image_pyramid_combined = []
        
        # -----------BUILD IMAGE PYRAMID FOR SELF MATRICES-----------
        if not args.compare_only:
            for j in range(min(len(seq_list),2)):
                print(f"Building matrix hierarchy for {seq_list[j]} with {args.layers} layers.... \n")
                for i in range(args.layers -1, -1, -1):
                    layer_sparsity = max(args.sparsity // (2**i), 1)
                    layer_window_size = args.window // (2**i)
                    layer_neigh = partition_overlaps(k_list[j], layer_window_size, args.delta, len(k_list[j]), args.kmer)
                    layer_sing = partition_overlaps(k_list[j], layer_window_size, 0, len(k_list[j]), args.kmer)

                    mods_neigh = convert_to_modimizers(layer_neigh, layer_sparsity)
                    mods_sing = convert_to_modimizers(layer_sing, layer_sparsity)
                    print(f"Layer {i+1} using sparsity {args.sparsity // (2**i)}\n")
                    matrix_layer = self_containment_matrix(
                        mods_sing, mods_neigh, args.kmer, args.identity
                    )
                    if j == 0:
                        image_pyramid_1.insert(0, matrix_layer)
                    else:
                        image_pyramid_2.insert(0, matrix_layer)
        # -----------BUILD IMAGE PYRAMID FOR COMPARATIVE MATRICES-----------
        if (args.compare or args.compare_only) and len(seq_list) > 1:
            #Determine which is smaller, which is larger
            larger_name = ""
            smaller_name = ""
            larger_seq = []
            smaller_seq = []
            if len(k_list[0]) > len(k_list[1]):
                larger_name = seq_list[0]
                larger_seq = k_list[0]
                smaller_name = seq_list[1]
                smaller_seq = k_list[1]
            else:
                larger_name = seq_list[1]
                larger_seq = k_list[1]
                smaller_name = seq_list[0]
                smaller_seq = k_list[0]

            for i in range(args.layers -1, -1, -1):
                layer_sparsity = args.sparsity // (2**i)
                layer_window_size = args.window // (2**i)
                larger_neigh = partition_overlaps(larger_seq, layer_window_size, args.delta, len(larger_seq), args.kmer)
                larger_sing = partition_overlaps(larger_seq, layer_window_size, 0, len(larger_seq), args.kmer)
                smaller_neigh = partition_overlaps(smaller_seq, layer_window_size, args.delta, len(smaller_seq), args.kmer)
                smaller_sing = partition_overlaps(smaller_seq, layer_window_size, 0, len(smaller_seq), args.kmer)
                print(len(larger_neigh), len(smaller_neigh))
                larger_mods_neigh = convert_to_modimizers(larger_neigh, layer_sparsity)
                larger_mods_sing = convert_to_modimizers(larger_sing, layer_sparsity)
                smaller_mods_neigh = convert_to_modimizers(smaller_neigh, layer_sparsity)
                smaller_mods_sing = convert_to_modimizers(smaller_sing, layer_sparsity)
                matrix_layer = pairwise_containment_matrix(
                    larger_mods_sing,
                    smaller_mods_sing,
                    larger_mods_neigh,
                    smaller_mods_neigh,
                    args.identity,
                    args.kmer,
                    False
                )
                print(matrix_layer.shape)
                print(matrix_layer)
                image_pyramid_combined.insert(0,matrix_layer)
        if args.compare_only:
            run_dash(
                len(k_list[0]),
                None,
                seq_list[0],
                seq_list[0],
                args.resolution,
                args.sparsity,
                args.kmer,
                args.identity,
                args.port,
                args.delta,
                image_pyramid_combined,
                None,
                None
            )
        else:
            run_dash(
                len(k_list[0]),
                None,
                seq_list[0],
                seq_list[0],
                args.resolution,
                args.sparsity,
                args.kmer,
                args.identity,
                args.port,
                args.delta,
                image_pyramid_1,
                image_pyramid_2,
                image_pyramid_combined
            )

    # -----------SETUP STATIC MODE-----------
    elif args.command == "static":
        # -----------COMPUTE SELF-IDENTITY PLOTS-----------
        if not args.compare_only:
            for i in range(len(seq_list)):
                win = args.window
                res = args.resolution
                if args.window:
                    # Change the resolution of each plot 
                    res = math.ceil(len(k_list[i])/args.window)
                else:
                    win = math.ceil(len(k_list[i])/args.resolution)

                #TODO: Update width and height once self identity plot sizes get fixed:
                xaxis = 0
                width = 0
                height = 0
                if isinstance(args.xaxis, int) or args.xaxis == None:
                    xaxis = args.xaxis
                else:
                    assert len(args.xaxis) == len(seq_list)
                    xaxis = args.xaxis[i]
                if isinstance(args.width, int):
                    width = args.width
                else:
                    assert len(args.width) == len(seq_list)
                    width = args.width[i]
                if isinstance(args.height, int):
                    height = args.height
                else:
                    assert len(args.height) == len(seq_list)
                    height = args.height[i]

                print(f"Computing self identity matrix for {seq_list[i]}... \n")
                #TODO: Logging here
                print(f"\tSparsity value s: {args.sparsity}\n")
                print(f"\tSequence length n: {len(k_list[i]) + args.kmer - 1}\n")
                print(f"\tResolution r: {res}\n")
                print(f"\tWindow size w: {win}\n")
                mod_list = get_mods(k_list[i], args.sparsity, res)
                mod_set_neighbors = convert_set_neighbors(mod_list, args.delta)

                new_windows = partition_evenly_spaced_modimizers(
                    mod_list if args.delta == 0 else mod_set_neighbors,
                    len(k_list[i]) + args.kmer - 1,
                    res
                )
                print(len(new_windows))

                paired_bed_file(
                    new_windows,
                    seq_list[i],
                    args.identity,
                    args.output_dir,
                    args.no_bed,
                    args.no_plot,
                    args.no_hist,
                    args.palette,
                    args.palette_orientation,
                    width,
                    height,
                    args.dpi,
                    args.kmer,
                    args.bin_freq,
                    xaxis,
                    args.colors,
                    args.breakpoints
                )

        # -----------COMPUTE COMAPRATIVE PLOTS-----------
        #TODO: Optimize computations so that largest sequence doesn't need to be redone all the time
        if (args.compare or args.compare_only) and len(seq_list) > 1:
            # Set window size to args.window. Otherwise, set it to n/resolution
            win = args.window
            res = args.resolution
            max_length = max(len(x) for x in k_list)
            if args.window:
                res = math.ceil(max_length/args.window)
            else:
                win = math.ceil(max_length/args.window)
            for i, seq_i in enumerate(seq_list[:-1]):
                for seq_j in seq_list[i + 1:]:
                    print(f"Computing comparative plot {seq_i} vs. {seq_j}... \n")
                    # Determine which seqeunce is smaller, which one is larger:
                    larger_seq, smaller_seq = (seq_i, seq_j) if len(k_list[i]) > len(k_list[seq_list.index(seq_j)]) else (seq_j, seq_i)
                    # Get the smaller/larger ratio. Modify the resolution of the smaller one based on the larger
                    ratio = round(args.resolution * len(k_list[seq_list.index(smaller_seq)]) / len(k_list[seq_list.index(larger_seq)]))
                    print(f"{ratio}")

                    larger_seq_mods = get_mods(k_list[seq_list.index(larger_seq)], args.sparsity, args.resolution)
                    larger_seq_mod_neighbors = convert_set_neighbors(larger_seq_mods, args.delta)
                    smaller_seq_mods = get_mods(k_list[seq_list.index(smaller_seq)], args.sparsity, ratio)
                    smaller_seq_mod_neighbors = convert_set_neighbors(smaller_seq_mods, args.delta)

                    short_dict, large_dict = partition_pairwise_modimizers_different_size(larger_seq_mod_neighbors, smaller_seq_mod_neighbors, ratio, win)
                    paired_bed_file_a_vs_b(
                        short_dict,
                        large_dict,
                        seq_i,
                        seq_j,
                        args.resolution,
                        args.identity,
                        args.no_bed,
                        args.no_plot,
                        args.palette,
                        args.palette_orientation,
                        args.width,
                        args.width,
                        args.dpi,
                        args.kmer,
                        args.bin_freq,
                        args.output_dir,
                        args.colors,
                        args.breakpoints,
                        args.compare_only,
                    )
                


if __name__ == "__main__":
    main()
