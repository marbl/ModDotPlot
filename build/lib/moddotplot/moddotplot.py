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
)
from moddotplot.interactive import run_dash
from moddotplot.const import ASCII_ART, VERSION

import argparse
import math
from moddotplot.static_plots import (
    paired_bed_file,
    paired_bed_file_a_vs_b,
    create_plots,
)
import itertools
import pandas as pd


def get_args_parse():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Mod.Plot, A Rapid and Interactive Visualization of Tandem Repeats.",
    )

    req_params = parser.add_argument_group("Required input")

    req_params.add_argument(
        "-i",
        "--input",
        required=True,
        default=argparse.SUPPRESS,
        help="Path to input fasta file(s)",
        nargs="+",
    )

    dist_params = parser.add_argument_group("Mod.Plot distance matrix commands")

    dist_params.add_argument(
        "-k", "--kmer", default=21, type=int, help="k-mer length. Must be < 32."
    )

    dist_params.add_argument(
        "-s",
        "--sparsity",
        help="Modimizer sparsity value. Higher value will reduce the number of modimizers, but will increase performance. Default set to 2 per Mbp of sequence, rounded up to nearest even integer)",
        default=None,
        type=int,
    )

    dist_params.add_argument(
        "-r", "--resolution", default=800, type=int, help="Dotplot resolution."
    )

    dist_params.add_argument(
        "-id", "--identity", default=86, type=int, help="Identity cutoff threshold."
    )

    dist_params.add_argument(
        "-a", "--alpha", default=0.4, type=float, help="Fraction of neighboring k-mers to include in identity estimation."
    )

    dist_params.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Name for bed file and plots. Will be set to input fasta file name if not provided.",
    )

    dist_params.add_argument(
        "-nc",
        "--non-canonical",
        default=False,
        help="Only consider forward strand when computing k-mers. Not recommended for most applications.",
    )

    dist_params.add_argument(
        "--compare",
        action="store_true",
        help="Create a dotplot with two different sequences. Valid for 2+ input fasta sequences.",
    )

    dist_params.add_argument(
        "--compare-only",
        action="store_true",
        help="Create a dotplot with two different sequences, skipping self-identity dotplot.",
    )

    interactive_params = parser.add_argument_group("Interactive plotting commands")

    interactive_params.add_argument(
        "--no-interactive",
        action="store_true",
        help="Prevent launching interactive mode: used for producing images only."
    )

    interactive_params.add_argument(
        "-l",
        "--layers",
        default=3,
        type=int,
        help="Number of image pyramid layers to compute when rendering interactive mode."
    )

    interactive_params.add_argument(
        "--save",
        action="store_true",
        help="Save image pyramid to . Only used in interactive mode.",
    )

    interactive_params.add_argument(
        "--load",
        default=None,
        type=str,
        help="Load previously computed image pyramid file. Only used in interactive mode.",
    )

    interactive_params.add_argument(
        "--port",
        default="8050",
        type=int,
        help="Port number for launching interactive mode on localhost. Only used in interactive mode.",
    )

    plot_params = parser.add_argument_group("Static plotting commands")

    plot_params.add_argument(
        "--no-bed", action="store_true", help="Don't output bed file."
    )

    plot_params.add_argument(
        "--no-plot", action="store_true", help="Don't output svg and png image files."
    )

    plot_params.add_argument(
        "--no-hist", action="store_true", help="Don't output histogram color legend."
    )

    plot_params.add_argument(
        "--width", default=9, type=float, help="Change default plot width."
    )

    plot_params.add_argument(
        "--height", default=5, type=float, help="Change default plot height."
    )

    plot_params.add_argument(
        "--xaxis", default=None, type=float, help="Change x axis for static plot."
    )

    plot_params.add_argument(
        "--dpi", default=300, type=int, help="Change default plot dpi."
    )

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
        "-b",
        "--bedfile",
        default=None,
        help="Take bedfile as input instead.",
    )

    plot_params.add_argument(
        "--bin-freq",
        action="store_true",
        help="By default, histograms are evenly spaced based on the number of colors and the identity threshold. Select this argument to bin based on the frequency of observed identity values.",
    )
    # TODO: implement logging options
    """logging_params = parser.add_argument_group("Logging options")

    logging_params.add_argument(
        "-q", "--quiet", action="store_true", help="Supress help text when running."
    )

    logging_params.add_argument(
        "-v", "--verbose", action="store_true", help="Add additional logging info when running."
    )"""

    args = parser.parse_args()

    return args


def main():
    print(ASCII_ART)
    print(f"v{VERSION} \n")
    args = get_args_parse()
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
        print(f"Using s = {args.sparsity}. \n")

    if not args.no_interactive:
        # Compute base layer for modimizers

        if len(seq_list) > 1 and (args.compare or args.compare_only):
            print(f"Building image pyramid with {args.layers} layers.... \n")

            #TODO: assert that args layers is a reasonable number

            image_pyramid = []
            for i in range(args.layers):
                
                mod_list_1 = get_mods(k_list[0], args.sparsity // (2**i), args.resolution * (2**i))
                mod_list_2 = get_mods(k_list[1], args.sparsity // (2**i), args.resolution * (2**i))
                mod_set_1 = convert_set(mod_list_1)
                mod_set_2 = convert_set(mod_list_2)
                mod_set_neighbors_1 = convert_set_neighbors(mod_list_1, args.alpha)
                mod_set_neighbors_2 = convert_set_neighbors(mod_list_2, args.alpha)
                print(f"Layer {i+1} \n")
                xd = pairwise_containment_matrix(mod_set_1, mod_set_2, mod_set_neighbors_1, mod_set_neighbors_2, args.identity, args.kmer, False)
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
                image_pyramid
            )
        else:
            if len(seq_list) > 1:
                # Only using the first seq!!!
                print(f"Multiple sequences detected, only one currently used.")

            print(f"Building image pyramid with {args.layers} layers.... \n")

            #TODO: assert that args layers is a reasonable number

            image_pyramid = []
            for i in range(args.layers):
                mod_list = get_mods(k_list[0], args.sparsity // (2**i), args.resolution * (2**i))
                mod_set = convert_set(mod_list)
                mod_set_neighbors = convert_set_neighbors(mod_list, args.alpha)
                print(f"Layer {i+1} \n")
                xd = self_containment_matrix(mod_set, mod_set_neighbors, args.kmer, args.identity)
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
                image_pyramid
            )

    else:
        if not args.compare_only:
            for i in range(len(seq_list)):
                print(f"Computing self identity matrix for {seq_list[i]}... \n")
                mod_list = get_mods(k_list[i], args.sparsity, args.resolution)
                windows = partition_evenly_spaced_modimizers(
                    mod_list, len(k_list[i]) + args.kmer - 1, args.resolution
                )
                paired_bed_file(
                    windows,
                    seq_list[i],
                    args.identity,
                    args.output_dir,
                    args.no_bed,
                    args.no_plot,
                    args.no_hist,
                    args.palette,
                    args.palette_orientation,
                    args.width,
                    args.height,
                    args.dpi,
                    args.kmer,
                    args.bin_freq,
                    args.xaxis,
                )
        # --------------------
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
                        args.output_dir
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
                    )


if __name__ == "__main__":
    main()
