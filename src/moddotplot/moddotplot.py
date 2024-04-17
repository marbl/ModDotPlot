#!/usr/bin/env python3
import sys
from moddotplot.parse_fasta import (
    readKmersFromFile,
    getInputHeaders,
    isValidFasta,
    extractFiles,
)

from moddotplot.estimate_identity import (
    convertToModimizers,
    selfContainmentMatrix,
    pairwiseContainmentMatrix,
    convertMatrixToBed,
    createSelfMatrix,
    createPairwiseMatrix,
    partitionOverlaps,
)
from moddotplot.interactive import run_dash
from moddotplot.const import ASCII_ART, VERSION

import argparse
import math
from moddotplot.static_plots import read_df_from_file, create_plots
import json
import numpy as np
import pickle
import os


def get_parser():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="ModDotPlot: Visualization of Complex Repeat Structures",
    )
    subparsers = parser.add_subparsers(
        dest="command", help="Choose mode: interactive or static"
    )
    interactive_parser = subparsers.add_parser(
        "interactive", help="Interactive mode commands"
    )
    static_parser = subparsers.add_parser("static", help="Static mode commands")

    # -----------INTERACTIVE MODE SUBCOMMANDS-----------
    interactive_input_group = interactive_parser.add_mutually_exclusive_group(
        required=True
    )

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

    interactive_parser.add_argument(
        "-k", "--kmer", default=21, type=int, help="k-mer length."
    )

    interactive_parser.add_argument(
        "-m",
        "--modimizer",
        help="Modimizer sketch size. A lower value will reduce the number of modimizers, but will increase performance. Must be less than window length `-w`. ",
        default=1000,
        type=int,
    )

    interactive_parser.add_argument(
        "-r",
        "--resolution",
        default=1000,
        type=int,
        help="Dotplot resolution, or the number of intervals to compare against.",
    )

    interactive_parser.add_argument(
        "-w",
        "--window",
        default=None,
        type=int,
        help="Window size, or the length in genomic coordinates of each interval. Default is set to (genome length)/(resolution)",
    )

    interactive_parser.add_argument(
        "-id",
        "--identity",
        default=86.0,
        type=float,
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
        help="Create a dotplot for each pairwise combination of input sequences (in addition to self-identity plots).",
    )

    compare_group.add_argument(
        "--compare-only",
        action="store_true",
        help="Create a dotplot for each pairwise combination of input sequences (skips self-identity plots).",
    )

    interactive_parser.add_argument(
        "-s",
        "--save",
        action="store_true",
        help="Save hierarchical matrices to file.",
    )

    interactive_parser.add_argument(
        "--port",
        default="8050",
        type=int,
        help="Port number for launching interactive mode on localhost. Only used in interactive mode.",
    )

    interactive_parser.add_argument(
        "--ambiguous",
        action="store_true",
        help="Preserve diagonal when handling strings of ambiguous homopolymers (eg. long runs of N's).",
    )

    interactive_parser.add_argument(
        "-q",
        "--quick",
        action="store_true",
        help="Launch a quick, non-interactive version of interactive mode.",
    )

    interactive_parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Prevent launching dash after saving. Must be used in combination with --save.",
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
    static_window_size_group = static_parser.add_mutually_exclusive_group(
        required=False
    )

    static_parser.add_argument(
        "-k", "--kmer", default=21, type=int, help="k-mer length."
    )

    static_parser.add_argument(
        "-m",
        "--modimizer",
        help="Modimizer sketch size. A lower value will reduce the number of modimizers, but will increase performance. Must be less than window length `-w`. ",
        default=1000,
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
        default=86.0,
        type=float,
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
        "--xaxis",
        default=None,
        type=float,
        nargs="+",
        help="Change x axis for self identity plots. Default is length of the sequence, in mbp.",
    )

    static_parser.add_argument("--dpi", default=600, type=int, help="Plot dpi.")

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

    static_parser.add_argument(
        "--ambiguous",
        action="store_true",
        help="Preserve diagonal when handling strings of ambiguous homopolymers (eg. long runs of N's).",
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
        if hasattr(args, "load") and args.load:
            print(f"Loading matrix hierarchy from {args.load}... \n")
            matrices, metadata = extractFiles(args.load)
            sparsity = math.ceil(
                (metadata[0]["max_window_size"] + 1) / metadata[0]["resolution"]
            )
            sparsity = 2 ** math.floor(math.log2(sparsity))
            axes = []
            for i in range(len(matrices)):
                matrix_axes = []
                counter = 0
                for matrix in matrices[i]:
                    x_axis = [
                        j * round(metadata[i]["x_size"] / matrix.shape[0])
                        for j in range(matrix.shape[0])
                    ]
                    y_axis = [
                        j * round(metadata[i]["y_size"] / matrix.shape[1])
                        for j in range(matrix.shape[1])
                    ]
                    x_axis.append(metadata[i]["x_size"])
                    y_axis.append(metadata[i]["y_size"])
                    matrix_axes.append(x_axis)
                    matrix_axes.append(y_axis)
                axes.append(matrix_axes)
            run_dash(
                matrices,
                metadata,
                axes,
                sparsity,
                args.identity,
                args.port,
                args.output_dir,
            )
            sys.exit(0)
    elif args.command == "static":
        print(f"Running ModDotPlot in static mode\n")
        # -----------CONFIG PARSING-----------
        # TODO: Change to yml file, add readme to config folder
        if args.config:
            with open(args.config, "r") as f:
                config = json.load(f)
                # TODO: Remove args that are interactive only
                args.fasta = config.get("fasta")
                args.bed = config.get("bed")

                # Distance matrix commands
                args.kmer = config.get("kmer", args.kmer)
                args.modimizer = config.get("modimizer", args.modimizer)
                args.resolution = config.get("resolution", args.resolution)
                args.window = config.get("window", args.window)
                args.identity = config.get("identity", args.identity)
                args.delta = config.get("delta", args.delta)
                args.output_dir = config.get("output_dir", args.output_dir)
                args.compare = config.get("compare", args.compare)
                args.compare_only = config.get("compare_only", args.compare_only)

                args.no_bed = config.get("no_bed", args.no_bed)
                args.no_plot = config.get("no_plot", args.no_plot)
                args.no_hist = config.get("no_hist", args.no_hist)
                args.width = config.get("width", args.width)
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
                # TODO: Chronicle exit codes
                # Exit code 2: breakpoint value != identity threshold
                sys.exit(2)

        # -----------BEDFILE INPUT FOR STATIC MODE-----------
        if hasattr(args, "bed") and args.bed:
            try:
                for bed in args.bed:
                    # If args.bed is provided as input, run static mode directly from the bed file. Skip counting input k-mers.
                    df = read_df_from_file(bed)
                    unique_query_names = df["#query_name"].unique()
                    unique_reference_names = df["reference_name"].unique()
                    assert len(unique_query_names) == len(unique_reference_names)
                    # TODO: Change this to allow for multiple seqs in bed file
                    assert len(unique_reference_names) == 1
                    self_id_scores = df[df["#query_name"] == df["reference_name"]]
                    pairwise_id_scores = df[df["#query_name"] != df["reference_name"]]
                    print(
                        f"Input bed file {bed} read successfully! Creating plots... \n"
                    )
                    # Create directory
                    if not args.output_dir:
                        args.output_dir = os.getcwd()
                    if not os.path.exists(args.output_dir):
                        os.makedirs(args.output_dir)
                    if len(self_id_scores) > 1:
                        create_plots(
                            sdf=None,
                            directory=args.output_dir if args.output_dir else ".",
                            name_x=unique_query_names[0],
                            name_y=unique_query_names[0],
                            palette=args.palette,
                            palette_orientation=args.palette_orientation,
                            no_hist=args.no_hist,
                            width=args.width,
                            dpi=args.dpi,
                            is_freq=args.bin_freq,
                            xlim=None,  # TODO: Get xlim working
                            custom_colors=args.colors,
                            custom_breakpoints=args.breakpoints,
                            from_file=df,
                            is_pairwise=False,
                        )
                    # Case 2: Pairwise bed file
                    if len(pairwise_id_scores) > 1:
                        create_plots(
                            sdf=None,
                            directory=args.output_dir if args.output_dir else ".",
                            name_x=unique_query_names[0],
                            name_y=unique_reference_names[0],
                            palette=args.palette,
                            palette_orientation=args.palette_orientation,
                            no_hist=args.no_hist,
                            width=args.width,
                            dpi=args.dpi,
                            is_freq=args.bin_freq,
                            xlim=None,  # TODO: Get xlim working
                            custom_colors=args.colors,
                            custom_breakpoints=args.breakpoints,
                            from_file=df,
                            is_pairwise=True,
                        )
                # Exit once all bed files have been iterated through
                sys.exit(0)
            except Exception as e:
                # Exit code 7: Error getting info from bed file:
                # TODO: Change to logs
                print(f"Error in bed file: {e}")
                sys.exit(7)

    # -----------INPUT SEQUENCE VALIDATION-----------
    seq_list = []
    for i in args.fasta:
        isValidFasta(i)
        headers = getInputHeaders(i)
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
        kmer_list.append(readKmersFromFile(i, args.kmer, False))
    k_list = [item for sublist in kmer_list for item in sublist]

    # Throw error if compare only selected with one sequence.
    if len(k_list) < 2 and args.compare_only:
        print(f"Error: Can't create a comparative plot with only one sequence")
        sys.exit(2)

    # -----------LAUNCH INTERACTIVE MODE-----------
    if args.command == "interactive":
        # Single sequence, can set window length immediately.
        hgi = len(max(k_list))
        hgi = hgi + args.kmer - 1
        min_window_size = 0
        window_lengths = []
        if not args.window:
            if args.quick:
                min_window_size = round((hgi / args.resolution))
                args.window = min_window_size
            else:
                min_window_size = round((hgi / args.resolution) / 2)
                args.window = min_window_size
        else:
            min_window_size = args.window
            if args.window and args.quick:
                print(f"Conflict with `--quick` argument.")
        max_window_size = math.ceil(hgi / args.resolution)
        # If only sequence is too small, throw an error.
        if max_window_size < 10:
            print(
                    f"Error: sequence too small for analysis.\n"
                )
            print(
                f"ModDotPlot requires a minimum window size of 10. Sequences less than 10Kbp will not work with ModDotPlot under normal resolution. We recommend rerunning ModDotPlot with --r {math.ceil(hgi / 10)}.\n"
            )
            sys.exit(0)
        while min_window_size <= max_window_size:
            window_lengths.append(min_window_size)
            min_window_size = min_window_size * 2

        # Set warning if >2 sequences detected
        if len(seq_list) > 2:
            print(
                f"{len(seq_list)} sequences were detected, however interactive mode can only load two sequences at a time.\n"
            )
            print(
                f"Interactive mode will proceed with {seq_list[0]} and {seq_list[1]}\n"
            )
        # Set warning if <1 sequence detected
        elif len(seq_list) < 1:
            print(f"Error: No sequences detected!")
            sys.exit(5)

        # Set sparsity to be the closest power of 2
        sparsities = []
        if window_lengths[0] < 1000:
            sparsities.append(1)
        else:
            sparsities.append(round(window_lengths[0] / args.modimizer))
        if sparsities[0] <= args.modimizer:
            sparsities[0] = 2 ** int(math.log2(sparsities[0]))
        else:
            sparsities[0] = 2 ** (int(math.log2(sparsities[0] - 1)) + 1)
        # expectation = round(win/seq_sparsity)
        for i in range(1, len(window_lengths)):
            if window_lengths[i] > 1000:
                sparsities.append(sparsities[-1] * 2)
            else:
                sparsities.append(1)
        expectation = round(window_lengths[-1] / sparsities[-1])
        matrices = []
        metadata = []
        # -----------BUILD IMAGE PYRAMID FOR SELF MATRICES-----------
        if not args.compare_only:
            for j in range(min(len(seq_list), 2)):
                image_pyramid = []
                if args.quick or len(window_lengths) == 1:
                    print(
                        f"Building 1 self-identity matrix for {seq_list[j]}, using a window size of {window_lengths[0]}.... \n"
                    )
                else:
                    print(
                        f"Building {len(window_lengths)} self-identity matrices for {seq_list[j]}, using a minimum window size of {window_lengths[0]}.... \n"
                    )
                if not args.quick:
                    print(
                        f"Creating base layer using window length {window_lengths[0]}...\n"
                    )

                for i in range(len(window_lengths)):
                    layer_sparsity = sparsities[i]
                    layer_window_size = window_lengths[i]
                    layer_neigh = partitionOverlaps(
                        k_list[j],
                        layer_window_size,
                        args.delta,
                        len(k_list[j]),
                        args.kmer,
                    )
                    layer_sing = partitionOverlaps(
                        k_list[j], layer_window_size, 0, len(k_list[j]), args.kmer
                    )

                    mods_neigh = convertToModimizers(
                        layer_neigh,
                        layer_sparsity,
                        args.ambiguous,
                        args.kmer,
                        expectation,
                    )
                    mods_sing = convertToModimizers(
                        layer_sing,
                        layer_sparsity,
                        args.ambiguous,
                        args.kmer,
                        expectation,
                    )
                    if not args.quick and i > 0:
                        print(f"Layer {i+1} using window length {layer_window_size}\n")
                    matrix_layer = selfContainmentMatrix(
                        mods_sing, mods_neigh, args.kmer, args.identity, args.ambiguous
                    )
                    image_pyramid.insert(0, matrix_layer)
                matrices.append(image_pyramid)
                metadata.append(
                    {
                        "x_name": seq_list[j],
                        "y_name": seq_list[j],
                        "x_size": len(k_list[j]) + args.kmer - 1,
                        "y_size": len(k_list[j]) + args.kmer - 1,
                        "self": True,
                        "min_window_size": window_lengths[0],
                        "max_window_size": window_lengths[-1],
                        "resolution": args.resolution,
                        "kmer_length": args.kmer,
                        "title": f"{seq_list[j]}",
                        "sparsities": sparsities,
                    }
                )
        # -----------BUILD IMAGE PYRAMID FOR COMPARATIVE MATRICES-----------
        if (args.compare or args.compare_only) and len(seq_list) > 1:
            # Determine which is smaller, which is larger
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
            if args.quick:
                print(
                    f"Quickly building pairwise matrices for {seq_list[0]} and {seq_list[1]}, using a window size of {window_lengths[0]}.... \n"
                )
            else:
                print(
                    f"Building pairwise matrices for {seq_list[0]} and {seq_list[1]}, using a minimum window size of {window_lengths[0]}.... \n"
                )
            image_pyramid = []
            for i in range(len(window_lengths)):
                layer_sparsity = sparsities[i]
                layer_window_size = window_lengths[i]
                larger_neigh = partitionOverlaps(
                    larger_seq,
                    layer_window_size,
                    args.delta,
                    len(larger_seq),
                    args.kmer,
                )
                larger_sing = partitionOverlaps(
                    larger_seq, layer_window_size, 0, len(larger_seq), args.kmer
                )
                smaller_neigh = partitionOverlaps(
                    smaller_seq,
                    layer_window_size,
                    args.delta,
                    len(smaller_seq),
                    args.kmer,
                )
                smaller_sing = partitionOverlaps(
                    smaller_seq, layer_window_size, 0, len(smaller_seq), args.kmer
                )

                larger_mods_neigh = convertToModimizers(
                    larger_neigh, layer_sparsity, args.ambiguous, args.kmer, expectation
                )
                larger_mods_sing = convertToModimizers(
                    larger_sing, layer_sparsity, args.ambiguous, args.kmer, expectation
                )
                smaller_mods_neigh = convertToModimizers(
                    smaller_neigh,
                    layer_sparsity,
                    args.ambiguous,
                    args.kmer,
                    expectation,
                )
                smaller_mods_sing = convertToModimizers(
                    smaller_sing, layer_sparsity, args.ambiguous, args.kmer, expectation
                )
                if not args.quick:
                    print(f"Layer {i+1} using window length {layer_window_size}\n")
                matrix_layer = pairwiseContainmentMatrix(
                    larger_mods_sing,
                    smaller_mods_sing,
                    larger_mods_neigh,
                    smaller_mods_neigh,
                    args.identity,
                    args.kmer,
                    False,
                )
                image_pyramid.insert(0, matrix_layer)
            matrices.append(image_pyramid)
            metadata.append(
                {
                    "x_name": larger_name,
                    "y_name": smaller_name,
                    "x_size": len(larger_seq) + args.kmer - 1,
                    "y_size": len(smaller_seq) + args.kmer - 1,
                    "self": False,
                    "min_window_size": window_lengths[0],
                    "max_window_size": window_lengths[-1],
                    "resolution": args.resolution,
                    "kmer_length": args.kmer,
                    "title": f"{larger_name}-{smaller_name}",
                    "sparsities": sparsities,
                }
            )

        if args.save:
            # Check if this value already exists
            if not args.output_dir:
                args.output_dir = os.getcwd()
            folder_path = os.path.join(args.output_dir, "interactive_matrices")
            if not os.path.exists(folder_path):
                print(f"Saving interactive matrices in {folder_path}\n")
                os.makedirs(folder_path)
            else:
                print(f"Saving interactive matrices in {folder_path}\n")
                print(f"{folder_path} already exists, overwriting its contents.\n")

            for i in range(len(metadata)):
                for j in range(len(matrices[i])):
                    if metadata[i]["self"]:
                        tmp = f"{metadata[i]['x_name']}_{j}.npz"
                    else:
                        tmp = f"{metadata[i]['x_name']}-{metadata[i]['y_name']}_{j}.npz"
                    saved_path = os.path.join(folder_path, tmp)
                    np.savez_compressed(saved_path, data=matrices[i][j])
            pickle_path = os.path.join(folder_path, "metadata.pkl")
            # Save the dictionary as a pickle file
            with open(pickle_path, "wb") as f:
                pickle.dump(metadata, f)
            # Check if no plot arg is used
            if args.no_plot:
                print(f"Saved matrices to {folder_path}. Thank you for using ModDotPlot!\n")
                sys.exit(0)


        # Before running dash, change into intervals...
        axes = []
        for matrices_set, meta in zip(matrices, metadata):
            matrix_axes = []
            for matrix in matrices_set:
                x_axis = np.linspace(0, meta["x_size"], matrix.shape[0] + 1)
                y_axis = np.linspace(0, meta["y_size"], matrix.shape[1] + 1)
                matrix_axes.append(x_axis)
                matrix_axes.append(y_axis)
            axes.append(matrix_axes)
        run_dash(
            matrices,
            metadata,
            axes,
            sparsities[0],
            args.identity,
            args.port,
            args.output_dir,
        )

    # -----------SETUP STATIC MODE-----------
    elif args.command == "static":
        # -----------SET SPARSITY VALUE-----------
        sequences = list(zip(seq_list, k_list))
        sequences.sort(key=lambda x: len(x[1]), reverse=True)

        # Create output directory, if doesn't exist:
        if (args.output_dir) and not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
        # -----------COMPUTE SELF-IDENTITY PLOTS-----------
        if not args.compare_only:
            for i in range(len(sequences)):
                seq_length = len(sequences[i][1])
                win = args.window
                res = args.resolution
                if args.window:
                    # Change the resolution of each plot
                    res = math.ceil(seq_length / args.window)
                else:
                    win = math.ceil(seq_length / args.resolution)
                if win < args.modimizer:
                    raise ValueError(
                        "Window size must be greater than or equal to the modimizer sketch size"
                    )

                seq_sparsity = round(win / args.modimizer)
                if seq_sparsity <= args.modimizer:
                    seq_sparsity = 2 ** int(math.log2(seq_sparsity))
                else:
                    seq_sparsity = 2 ** (int(math.log2(seq_sparsity - 1)) + 1)
                expectation = round(win / seq_sparsity)
                xaxis = 0
                width = 0
                if isinstance(args.xaxis, int) or args.xaxis == None:
                    xaxis = args.xaxis
                else:
                    assert len(args.xaxis) == len(sequences)
                    xaxis = args.xaxis[i]
                if isinstance(args.width, int):
                    width = args.width
                else:
                    assert len(args.width) == len(sequences)
                    width = args.width[i]

                print(f"Computing self identity matrix for {sequences[i][0]}... \n")
                # TODO: Logging here
                # print(f"\tSparsity value s: {seq_sparsity}\n")
                print(f"\tSequence length n: {len(k_list[i]) + args.kmer - 1}\n")
                print(f"\tWindow size w: {win}\n")
                print(f"\tModimizer sketch value: {expectation}\n")
                print(f"\tPlot Resolution r: {res}\n")
                self_mat = createSelfMatrix(
                    seq_length,
                    sequences[i][1],
                    win,
                    seq_sparsity,
                    args.delta,
                    args.kmer,
                    args.identity,
                    args.ambiguous,
                    expectation,
                )
                bed = convertMatrixToBed(
                    self_mat, win, args.identity, seq_list[i], seq_list[i], True
                )

                if not args.no_bed:
                    # Log saving bed file
                    if not args.output_dir:
                        bedfile_output = sequences[i][0] + ".bed"
                    else:
                        bedfile_output = os.path.join(
                            args.output_dir, sequences[i][0] + ".bed"
                        )
                    with open(bedfile_output, "w") as bedfile:
                        for row in bed:
                            bedfile.write("\t".join(map(str, row)) + "\n")
                    print(f"Saved bed file to {bedfile_output}\n")

                if not args.no_plot:
                    create_plots(
                        sdf=[bed],
                        directory=args.output_dir if args.output_dir else ".",
                        name_x=sequences[i][0],
                        name_y=sequences[i][0],
                        palette=args.palette,
                        palette_orientation=args.palette_orientation,
                        no_hist=args.no_hist,
                        width=width,
                        dpi=args.dpi,
                        is_freq=args.bin_freq,
                        xlim=xaxis,
                        custom_colors=args.colors,
                        custom_breakpoints=args.breakpoints,
                        from_file=None,
                        is_pairwise=False,
                    )

        # -----------COMPUTE COMAPRATIVE PLOTS-----------
        # TODO: Optimize computations so that largest sequence doesn't need to be redone all the time
        if (args.compare or args.compare_only) and len(sequences) > 1:
            # Set window size to args.window. Otherwise, set it to n/resolution
            win = args.window
            res = args.resolution
            if args.window:
                res = math.ceil(len(sequences[0][1]) / args.window)
            else:
                win = math.ceil(len(sequences[0][1]) / args.resolution)
            if win < args.modimizer:
                raise ValueError(
                    "Window size must be greater than or equal to the modimizer sketch size"
                )

            seq_sparsity = round(win / args.modimizer)
            if seq_sparsity <= args.modimizer:
                seq_sparsity = 2 ** int(math.log2(seq_sparsity))
            else:
                seq_sparsity = 2 ** (int(math.log2(seq_sparsity - 1)) + 1)
            expectation = round(win / seq_sparsity)

            for i in range(len(sequences)):
                larger_seq = sequences[i][1]

                for j in range(i + 1, len(sequences)):
                    smaller_seq = sequences[j][1]
                    larger_length = len(larger_seq)
                    smaller_length = len(smaller_seq)
                    print(
                        f"Computing pairwise identity matrix for {sequences[i][0]} and {sequences[j][0]}... \n"
                    )
                    # TODO: Logging here
                    # print(f"\tSparsity value s: {seq_sparsity}\n")
                    print(
                        f"\tSequence length {sequences[i][0]}: {larger_length + args.kmer - 1}\n"
                    )
                    print(
                        f"\tSequence length {sequences[j][0]}: {smaller_length + args.kmer - 1}\n"
                    )
                    print(f"\tWindow size w: {win}\n")
                    print(f"\tModimizer sketch value: {expectation}\n")
                    print(f"\tPlot Resolution r: {res}\n")

                    pair_mat = createPairwiseMatrix(
                        larger_length,
                        smaller_length,
                        sequences[i][1],
                        sequences[j][1],
                        win,
                        seq_sparsity,
                        args.delta,
                        args.kmer,
                        args.identity,
                        args.ambiguous,
                        expectation,
                    )
                    # Throw error if the matrix is empty
                    if np.all(pair_mat == 0):
                        print(
                            f"The pairwise identity matrix for {sequences[i][0]} and {sequences[j][0]} is empty. Skipping.\n"
                        )
                    else:
                        bed = convertMatrixToBed(
                            pair_mat,
                            win,
                            args.identity,
                            seq_list[i],
                            seq_list[j],
                            False,
                        )

                        if not args.no_bed:
                            # Log saving bed file
                            if not args.output_dir:
                                bedfile_output = sequences[i][0] + ".bed"
                            else:
                                bedfile_output = os.path.join(
                                    args.output_dir,
                                    sequences[i][0] + "_" + sequences[j][0] + ".bed",
                                )
                            with open(bedfile_output, "w") as bedfile:
                                for row in bed:
                                    bedfile.write("\t".join(map(str, row)) + "\n")
                            print(f"Saved bed file to {bedfile_output}\n")

                        if not args.no_plot:
                            create_plots(
                                sdf=[bed],
                                directory=args.output_dir if args.output_dir else ".",
                                name_x=sequences[i][0],
                                name_y=sequences[j][0],
                                palette=args.palette,
                                palette_orientation=args.palette_orientation,
                                no_hist=args.no_hist,
                                width=None,
                                dpi=args.dpi,
                                is_freq=args.bin_freq,
                                xlim=None,
                                custom_colors=args.colors,
                                custom_breakpoints=args.breakpoints,
                                from_file=None,
                                is_pairwise=True,
                            )


if __name__ == "__main__":
    main()
