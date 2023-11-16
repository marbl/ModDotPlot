#!/usr/bin/env python3
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
    self_containment_matrix,
    pairwise_containment_matrix,
    next_power_of_two,
)
from moddotplot.interactive import run_dash
from moddotplot.const import ASCII_ART, VERSION

import argparse
import math
from moddotplot.static_plots import (
    paired_bed_file,
    paired_bed_file_a_vs_b,
)
import itertools
import json


def get_args_parse():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Mod.Plot, Visualization of Complex Repeat Structures.",
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        "-i",
        "--input",
        default=argparse.SUPPRESS,
        help="Path to input files. Accepts fasta file(s).",
        nargs="+",
    )

    group.add_argument(
        "-c",
        "--config",
        default=None,
        type=str,
        help="Config file to use. Takes precedence over any other competing command line arguments.",
    )

    dist_params = parser.add_argument_group("Mod.Plot distance matrix commands")

    # Add a mutually exclusive group for compare and compare only.
    dist_group = parser.add_mutually_exclusive_group(required=False)

    dist_params.add_argument(
        "-k", "--kmer", default=21, type=int, help="k-mer length. Must be < 32."
    )

    dist_params.add_argument(
        "-s",
        "--sparsity",
        help="Modimizer sparsity value. Higher value will reduce the number of modimizers, but will increase performance. Default set to 2 per Mbp of sequence, rounded up to nearest even integer.",
        default=None,
        type=int,
    )

    dist_params.add_argument(
        "-r",
        "--resolution",
        default=1000,
        type=int,
        help="Dotplot resolution, or the number of cells to compare against. Prioritized over window size.",
    )

    # TODO: Implement this functionality
    dist_params.add_argument(
        "-w",
        "--window",
        default=None,
        type=int,
        help="Window size, or how many k-mers partitoned per cell. Used in lieu of resolution",
    )

    dist_params.add_argument(
        "-id",
        "--identity",
        default=86,
        type=int,
        help="Identity cutoff threshold. Set higher for faster computations.",
    )

    dist_params.add_argument(
        "-a",
        "--alpha",
        default=0.2,
        type=float,
        help="Fraction of neighboring k-mers to include in identity estimation.",
    )

    dist_params.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Directory name for bed file, plots, and stored matrices. Defaults to working directory.",
    )

    dist_group.add_argument(
        "--compare",
        action="store_true",
        help="Create a dotplot with two different sequences. Valid for 2+ input fasta sequences.",
    )

    dist_group.add_argument(
        "--compare-only",
        action="store_true",
        help="Create a dotplot with two different sequences, skipping self-identity dotplot.",
    )

    interactive_params = parser.add_argument_group("Interactive plotting commands")

    interactive_params.add_argument(
        "-l",
        "--layers",
        default=3,
        type=int,
        help="Number of matrix hierarchy layers to compute when preparing interactive mode.",
    )

    # TODO: Implement this functionality
    """interactive_params.add_argument(
        "--save",
        action="store_true",
        help="Save hierarchical matrices to file. Only used in interactive mode.",
    )"""

    # TODO: Implement this functionality
    """interactive_params.add_argument(
        "--load",
        default=None,
        type=str,
        help="Load previously computed hierarchical matrices.",
    )"""

    interactive_params.add_argument(
        "--port",
        default="8050",
        type=int,
        help="Port number for launching interactive mode on localhost. Only used in interactive mode.",
    )

    plot_params = parser.add_argument_group("Static plotting commands")

    plot_params.add_argument(
        "--static",
        action="store_true",
        help="Create static bitmap and vector dotplots. Prevents launching of Plotly & Dash. ",
    )

    plot_params.add_argument(
        "--no-bed", action="store_true", help="Don't output bed file."
    )

    plot_params.add_argument(
        "--no-plot", action="store_true", help="Don't output pdf and png image files."
    )

    plot_params.add_argument(
        "--no-hist", action="store_true", help="Don't output histogram color legend."
    )

    plot_params.add_argument(
        "--width", default=9, type=float, nargs="+", help="Change default plot width."
    )

    plot_params.add_argument(
        "--height", default=5, type=float, nargs="+", help="Change default plot height."
    )

    plot_params.add_argument(
        "--xaxis",
        default=None,
        type=float,
        nargs="+",
        help="Change x axis for static plot. Default is length of the sequence, in mbp.",
    )

    plot_params.add_argument(
        "--dpi", default=600, type=int, help="Change default plot dpi."
    )

    # TODO: Create list of accepted colors.
    plot_params.add_argument(
        "--palette",
        default="Spectral_11",
        help="Select color palette. See RColorBrewer for list of accepted palettes. Will default to Spectral_11 if not used.",
        type=str,
    )

    plot_params.add_argument(
        "--palette-orientation",
        default="+",
        choices=["+", "-"],
        help="Color palette orientation. + for forward, - for reverse.",
        type=str,
    )

    plot_params.add_argument(
        "--colors",
        default=None,
        nargs="+",
        help="Use a custom color palette, entered in either hexcode or rgb format.",
    )

    plot_params.add_argument(
        "--breakpoints",
        default=None,
        nargs="+",
        help="Introduce custom color thresholds. Must be between identity threshold and 100.",
    )

    plot_params.add_argument(
        "--bin-freq",
        action="store_true",
        help="By default, histograms are evenly spaced based on the number of colors and the identity threshold. Select this argument to bin based on the frequency of observed identity values.",
    )
    # TODO: implement logging options

    args = parser.parse_args()

    return args


def main():
    print(ASCII_ART)
    print(f"v{VERSION} \n")
    args = get_args_parse()

    # If config file selected, parse those arguments first
    if args.config:
        with open(args.config, "r") as f:
            config = json.load(f)
            args.input = config.get("input")

            # Distance matrix commands
            args.kmer = config.get("kmer", args.kmer)
            args.sparsity = config.get("sparsity", args.sparsity)
            args.resolution = config.get("resolution", args.resolution)
            args.window = config.get("window", args.window)
            args.identity = config.get("identity", args.identity)
            args.alpha = config.get("alpha", args.alpha)
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

    # Tests for specific bugs here
    # TODO: More tests!
    if args.breakpoints:
        # Check that start value fro breakpoints = identity
        if args.breakpoints[0] != args.identity:
            print(
                f"Identity threshold is {args.identity}, but starting breakpoint is {args.breakpoints[0]}! \n"
            )
            print(
                f"Please modify identity threshold using --identity {args.breakpoints[0]}.\n"
            )
            sys.exit(2)

    # Validate input sequence:
    seq_list = []
    for i in args.input:
        is_valid_fasta(i)
        headers = get_input_headers(i)
        if len(headers) > 1:
            print(f"File {i} contains multiple fasta entries: \n")
            counter = 1
            for j in headers:
                print(f"{counter}) {j} \n")
                counter += 1
        for j in headers:
            seq_list.append(j)

    kmer_list = []
    for i in args.input:
        kmer_list.append(read_kmers_from_file(i, args.kmer))
    k_list = [item for sublist in kmer_list for item in sublist]

    # Set sparsity value if not determined
    if not args.sparsity:
        len_list = []
        for x in k_list:
            len_list.append(len(x))
        args.sparsity = math.ceil(max(len_list) / 500000)

    if not args.static:
        args.sparsity = next_power_of_two(args.sparsity)
        print(f"Setting top layer sparsity = {args.sparsity}. \n")
        # Compute base layer for modimizers

        if len(seq_list) > 1 and (args.compare or args.compare_only):
            print(f"Building matrix hierarchy with {args.layers} layers.... \n")

            # TODO: assert that args layers is a reasonable number
            image_pyramid = []
            for i in range(args.layers):
                mod_list_1 = get_mods(
                    k_list[0], args.sparsity // (2**i), args.resolution * (2**i)
                )
                mod_list_2 = get_mods(
                    k_list[1], args.sparsity // (2**i), args.resolution * (2**i)
                )
                mod_set_1 = convert_set(mod_list_1)
                mod_set_2 = convert_set(mod_list_2)
                mod_set_neighbors_1 = convert_set_neighbors(mod_list_1, args.alpha)
                mod_set_neighbors_2 = convert_set_neighbors(mod_list_2, args.alpha)
                print(f"Layer {i+1} using sparsity {args.sparsity // (2**i)}\n")
                xd = pairwise_containment_matrix(
                    mod_set_1,
                    mod_set_2,
                    mod_set_neighbors_1,
                    mod_set_neighbors_2,
                    args.identity,
                    args.kmer,
                    False,
                )
                image_pyramid.append(xd)

            run_dash(
                k_list[0],
                k_list[1],
                seq_list[0],
                seq_list[1],
                args.resolution,
                args.sparsity,
                args.kmer,
                args.identity,
                args.port,
                args.palette,
                args.palette_orientation,
                args.alpha,
                image_pyramid,
            )
        else:
            if len(seq_list) > 1:
                # Only using the first seq!!!
                print(f"Multiple sequences detected, only one currently used.")

            print(f"Building matrix hierarchy with {args.layers} layers.... \n")

            # TODO: assert that args layers is a reasonable number

            image_pyramid = []
            for i in range(args.layers):
                mod_list = get_mods(
                    k_list[0], args.sparsity // (2**i), args.resolution * (2**i)
                )
                mod_set = convert_set(mod_list)
                mod_set_neighbors = convert_set_neighbors(mod_list, args.alpha)
                print(f"Layer {i+1} using sparsity {args.sparsity // (2**i)}\n")
                xd = self_containment_matrix(
                    mod_set, mod_set_neighbors, args.kmer, args.identity
                )
                image_pyramid.append(xd)

            run_dash(
                k_list[0],
                None,
                seq_list[0],
                seq_list[0],
                args.resolution,
                args.sparsity,
                args.kmer,
                args.identity,
                args.port,
                args.palette,
                args.palette_orientation,
                args.alpha,
                image_pyramid,
            )

    # Static plots
    else:
        print(f"Using s = {args.sparsity}. \n")
        if not args.compare_only:
            for i in range(len(seq_list)):
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
                mod_list = get_mods(k_list[i], args.sparsity, args.resolution)
                mod_set_neighbors_hi = convert_set_neighbors(mod_list, args.alpha)
                # TODO: make cleaner
                new_windows = []
                if args.alpha == 0:
                    new_windows = partition_evenly_spaced_modimizers(
                        mod_list, len(k_list[i]) + args.kmer - 1, args.resolution
                    )
                else:
                    new_windows = partition_evenly_spaced_modimizers(
                        mod_set_neighbors_hi,
                        len(k_list[i]) + args.kmer - 1,
                        args.resolution,
                    )
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
                    args.breakpoints,
                )
        # Compare only
        if (args.compare or args.compare_only) and len(seq_list) > 1:
            seq_kmers = {}
            for j in range(len(seq_list)):
                seq_kmers[seq_list[j]] = k_list[j]
            for i in itertools.combinations(seq_kmers, 2):
                print(f"Computing {i[0]} vs. {i[1]}... \n")
                ratio = 0
                if len(seq_kmers[i[0]]) > len(seq_kmers[i[1]]):
                    ratio = round(
                        args.resolution * len(seq_kmers[i[1]]) / len(seq_kmers[i[0]])
                    )
                    first = get_mods(seq_kmers[i[1]], args.sparsity, ratio)
                    second = get_mods(
                        seq_kmers[i[0]][0 : len(seq_kmers[i[1]])], args.sparsity, ratio
                    )
                    third = get_mods(
                        seq_kmers[i[0]][
                            len(seq_kmers[i[1]]) + 1 : len(seq_kmers[i[0]])
                        ],
                        args.identity,
                        ratio,
                    )

                    short_dict = {}
                    large_dict = {}
                    for w in range(ratio):
                        start_site = w * round(len(seq_kmers[i[1]]) / ratio)
                        end_site = (
                            w * round(len(seq_kmers[i[1]]) / ratio)
                            + round(len(seq_kmers[i[1]]) / ratio)
                            - 1
                        )
                        name = f"{start_site}-{end_site}"
                        short_dict[name] = first[w]
                        large_dict[name] = second[w]
                    for x in range(args.resolution - ratio):
                        start_site = (ratio + x) * round(len(seq_kmers[i[1]]) / ratio)
                        end_site = (
                            (ratio + x) * round(len(seq_kmers[i[1]]) / ratio)
                            + round(len(seq_kmers[i[1]]) / ratio)
                            - 1
                        )
                        name = f"{start_site}-{end_site}"
                        large_dict[name] = third[x]

                    paired_bed_file_a_vs_b(
                        short_dict,
                        large_dict,
                        i[0],
                        i[1],
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
                else:
                    ratio = round(
                        args.resolution * len(seq_kmers[i[0]]) / len(seq_kmers[i[1]])
                    )
                    first = get_mods(seq_kmers[i[0]], args.sparsity, ratio)
                    second = get_mods(
                        seq_kmers[i[1]][0 : len(seq_kmers[i[0]])], args.sparsity, ratio
                    )
                    third = []
                    if ratio != args.resolution:
                        third = get_mods(
                            seq_kmers[i[1]][
                                len(seq_kmers[i[0]]) + 1 : len(seq_kmers[i[1]])
                            ],
                            args.identity,
                            args.resolution - ratio,
                        )

                    short_dict = {}
                    large_dict = {}
                    for w in range(ratio):
                        start_site = w * round(len(seq_kmers[i[0]]) / ratio)
                        end_site = (
                            w * round(len(seq_kmers[i[0]]) / ratio)
                            + round(len(seq_kmers[i[0]]) / ratio)
                            - 1
                        )
                        name = f"{start_site}-{end_site}"
                        short_dict[name] = first[w]
                        large_dict[name] = second[w]
                    for x in range(args.resolution - ratio):
                        start_site = (ratio + x) * round(len(seq_kmers[i[0]]) / ratio)
                        end_site = (
                            (ratio + x) * round(len(seq_kmers[i[0]]) / ratio)
                            + round(len(seq_kmers[i[0]]) / ratio)
                            - 1
                        )
                        name = f"{start_site}-{end_site}"
                        large_dict[name] = third[x]

                    paired_bed_file_a_vs_b(
                        short_dict,
                        large_dict,
                        i[0],
                        i[1],
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
