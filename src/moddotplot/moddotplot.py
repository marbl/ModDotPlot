#!/usr/bin/env python3
import sys
from moddotplot.parse_fasta import (
    HASH_ALGORITHM,
    readKmersFromFile,
    getInputHeaders,
    isValidFasta,
    extractFiles,
    extractRegion,
)

from moddotplot.estimate_identity import (
    convertToModimizers,
    selfContainmentMatrix,
    pairwiseContainmentMatrix,
    convertMatrixToBed,
    convertMatrixToCool,
    createSelfMatrix,
    createPairwiseMatrix,
    create_self_matrix_from_sketches,
    create_pairwise_matrix_from_sketches,
    ModimizerSketchCache,
    partitionOverlaps,
)
from moddotplot.interactive import interactive_axis_bounds, run_dash
from moddotplot.annotations import read_annotation_beds
from moddotplot.const import ASCII_ART, VERSION

import argparse
import math
import json
import numpy as np
import pickle
import os
import shlex

from moddotplot.plot_summary import PlotSummaryWriter


# Static plotting pulls in the Plotnine and Matplotlib stacks. Keep those
# imports behind the static command boundary so ``--help`` and interactive
# mode do not pay their startup cost.
read_df_from_file = None
create_plots = None
create_grid = None

COMMANDS = frozenset({"interactive", "static"})
INTERACTIVE_DEPRECATION_MESSAGE = (
    "Warning: interactive mode is deprecated and maintenance-only. "
    "It remains available, but will not receive new features."
)


def _load_static_plotting():
    global read_df_from_file, create_plots, create_grid

    from moddotplot import static_plots

    if read_df_from_file is None:
        read_df_from_file = static_plots.read_df_from_file
    if create_plots is None:
        create_plots = static_plots.create_plots
    if create_grid is None:
        create_grid = static_plots.create_grid


def get_parser():
    """
    Argument parsing for stand-alone runs.

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="ModDotPlot: Visualization of Tandem Repeats",
    )
    subparsers = parser.add_subparsers(
        dest="command",
        required=False,
        help="Choose mode; static is used when omitted",
    )
    static_parser = subparsers.add_parser(
        "static", help="Static mode commands (default)"
    )
    interactive_parser = subparsers.add_parser(
        "interactive",
        help="Interactive mode commands (deprecated; explicit use only)",
        description=(
            "Deprecated interactive mode. This mode remains available but is "
            "maintenance-only and will not receive new features."
        ),
    )

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
        help="Fraction of each neighboring window included when estimating identity. Default: 0.5.",
    )

    interactive_parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Directory name for saving matrices and coordinate logs. Defaults to working directory.",
    )

    interactive_parser.add_argument(
        "-b",
        "--bed",
        default=None,
        nargs="+",
        help=(
            "BED3-BED9 annotation file(s). Tracks are shown for BED chromosome "
            "names matching the FASTA headers on each matrix axis."
        ),
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
        help="Include k-mer windows containing non-ACGTU IUPAC bases instead of masking them.",
    )

    interactive_parser.add_argument(
        "--forward",
        action="store_true",
        help="Enforce forward only k-mers instead of canonical k-mers. Warning: only use if you want strand-specific output!",
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
        "-l",
        "--load",
        default=argparse.SUPPRESS,
        help="Path to input paired-end bed file(s). Exclusively used in static mode.",
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

    static_parser.add_argument("-b", "--bed", default=None, help="Bed file annotation.")

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
        "--region",
        default=None,
        help="Genomic region to analyze. Syntax is seq_id:start-end (e.g. chr1:100000-200000).",
        nargs="+",
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
        help="Fraction of each neighboring window included when estimating identity. Default: 0.5.",
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
        "--compare-order",
        choices=["sequential", "size"],
        default="sequential",
        help="Order in which sequences appear in the comparative plot. Default is 'sequential': First file on x-axis, second file on y-axis. Another option is 'size': The larger sequence on the x-axis and the smaller on y-axis.",
    )

    static_parser.add_argument(
        "--cooler", action="store_true", help="Output matrix to cooler file."
    )

    static_parser.add_argument(
        "--no-bedpe", action="store_true", help="Skip output of paired-end bed file."
    )

    static_parser.add_argument(
        "--no-plot", action="store_true", help="Skip output of plots."
    )

    static_parser.add_argument(
        "--no-hist", action="store_true", help="Skip output of histogram color legend."
    )

    static_parser.add_argument(
        "--width", default=9, type=float, help="Plot width (also height for _FULL)."
    )

    static_parser.add_argument("--dpi", default=300, type=int, help="Plot dpi.")

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
        "--forward",
        action="store_true",
        help="Enforce forward only k-mers instead of canonical k-mers. Warning: only use if you want strand-specific output!",
    )

    static_parser.add_argument(
        "--plot-direction",
        action="store_true",
        help=(
            "Color matches in _FULL, _TRI, and grid plots by k-mer orientation: "
            "blue for the same orientation and pink for reverse orientation."
        ),
    )

    static_parser.add_argument(
        "--colors",
        "--color",
        dest="colors",
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
        "-a",
        "--axes-limits",
        default=None,
        type=float,
        help="Change x and y axis limits for self identity plots. Default is length of the sequence. Can't be shorter than length of sequence.",
    )

    static_parser.add_argument(
        "-t",
        "--axes-ticks",
        default=None,
        nargs="+",
        type=int,
        help="Tick labels to include in x and y axis for custom plots.",
    )
    # CURRENTLY NOT WORKING
    static_parser.add_argument(
        "--axes-number",
        default=7,
        help="Number of axis ticks labels to include in x and y axis for custom plots, including 0 and seq_length. A minimum of 2 is required, maximum 25.",
    )

    static_parser.add_argument(
        "--bin-freq",
        action="store_true",
        help="By default, histograms are evenly spaced based on the number of colors and the identity threshold. Select this argument to bin based on the frequency of observed identity values.",
    )

    static_parser.add_argument(
        "--ambiguous",
        action="store_true",
        help="Include k-mer windows containing non-ACGTU IUPAC bases instead of masking them.",
    )

    static_parser.add_argument(
        "--grid",
        action="store_true",
        help="Plot comparative plots in an NxN grid like format.",
    )

    static_parser.add_argument(
        "--grid-only",
        action="store_true",
        help="Plot comparative plots in an NxN grid like format, skipping individual plots.",
    )

    static_parser.add_argument(
        "--vector",
        choices=["svg", "pdf", "ps"],
        default="svg",
        help="Output format for vector format.",
    )

    static_parser.add_argument(
        "--deraster",
        action="store_true",
        help="De-rasterize dotplot in vector format. Note this can lead to large image sizes and make it unusable in image editing software.",
    )

    return parser


def _arguments_with_default_command(arguments=None):
    """Return CLI arguments with static mode inserted when no mode is given."""

    normalized = list(sys.argv[1:] if arguments is None else arguments)
    if normalized[:1] and normalized[0] in ("-h", "--help"):
        return normalized
    if not normalized or normalized[0] not in COMMANDS:
        normalized.insert(0, "static")
    return normalized


def parse_args(arguments=None):
    """Parse command-line arguments, defaulting omitted subcommands to static."""

    return get_parser().parse_args(_arguments_with_default_command(arguments))


def _apply_static_config(args, config):
    """Apply static-mode JSON configuration values to parsed arguments."""
    # TODO: Remove args that are interactive only
    args.fasta = config.get("fasta")
    args.load = config.get("load")
    args.bed = config.get("bed")

    # Distance matrix commands
    args.kmer = config.get("kmer", args.kmer)
    args.modimizer = config.get("modimizer", args.modimizer)
    args.resolution = config.get("resolution", args.resolution)
    args.window = config.get("window", args.window)
    args.region = config.get("region", args.region)
    args.identity = config.get("identity", args.identity)
    args.delta = config.get("delta", args.delta)
    args.output_dir = config.get("output_dir", args.output_dir)
    args.compare = config.get("compare", args.compare)
    args.compare_only = config.get("compare_only", args.compare_only)
    args.compare_order = config.get("compare_order", args.compare_order)

    args.cooler = config.get("cooler", args.cooler)
    args.no_bedpe = config.get("no_bedpe", args.no_bedpe)
    args.no_plot = config.get("no_plot", args.no_plot)
    args.no_hist = config.get("no_hist", args.no_hist)
    args.width = config.get("width", args.width)
    args.axes_limits = config.get("axes_limits", args.axes_limits)
    args.dpi = config.get("dpi", args.dpi)
    args.palette = config.get("palette", args.palette)
    args.palette_orientation = config.get(
        "palette_orientation", args.palette_orientation
    )
    args.colors = config.get("colors", config.get("color", args.colors))
    args.axes_ticks = config.get("axes_ticks", args.axes_ticks)
    args.axes_number = config.get("axes_number", args.axes_number)
    args.breakpoints = config.get("breakpoints", args.breakpoints)
    args.bin_freq = config.get("bin_freq", args.bin_freq)
    args.forward = config.get("forward", args.forward)
    args.plot_direction = config.get("plot_direction", args.plot_direction)
    args.ambiguous = config.get("ambiguous", args.ambiguous)
    args.grid = config.get("grid", args.grid)
    args.grid_only = config.get("grid_only", args.grid_only)
    args.vector = config.get("vector", args.vector)
    args.deraster = config.get("deraster", args.deraster)

    return args


def _parse_region_arguments(region_arguments, sequence_names):
    """Validate CLI regions and index them by their exact FASTA identifier."""

    if not region_arguments:
        return {}
    if isinstance(region_arguments, str):
        region_arguments = [region_arguments]

    available_names = {
        parsed[0] if (parsed := extractRegion(name)) else name
        for name in sequence_names
    }
    regions = {}
    for value in region_arguments:
        parsed = extractRegion(value)
        if not parsed:
            raise ValueError(f"invalid region {value!r}; expected FASTA_ID:start-end")
        sequence_name, start, end = parsed
        if start < 1 or end < start:
            raise ValueError(
                f"invalid region {value!r}; coordinates must satisfy 1 <= start <= end"
            )
        if sequence_name not in available_names:
            raise ValueError(f"region {value!r} does not match any FASTA identifier")
        if sequence_name in regions:
            raise ValueError(
                f"multiple regions were provided for FASTA identifier {sequence_name!r}"
            )
        regions[sequence_name] = parsed
    return regions


def _slice_kmers_for_region(kmers, region, kmer_size, sequence_start=1):
    """Return the exact k-mer slice for a 1-based inclusive base interval."""

    _sequence_name, start, end = region
    sequence_base_length = len(kmers) + kmer_size - 1
    sequence_end = sequence_start + sequence_base_length - 1
    if start < sequence_start or end > sequence_end:
        raise ValueError(
            f"region {start}-{end} is outside the available interval "
            f"{sequence_start}-{sequence_end}"
        )
    region_base_length = end - start + 1
    if region_base_length < kmer_size:
        raise ValueError(
            f"region length {region_base_length} is shorter than k-mer size {kmer_size}"
        )

    # A base interval [start, end] contains k-mers beginning at genomic
    # positions start through end - k + 1, inclusive. Translate those positions
    # into the available sequence's zero-based coordinates.
    local_start = start - sequence_start
    local_stop = end - sequence_start - kmer_size + 2
    selected = kmers[local_start:local_stop]
    expected_count = region_base_length - kmer_size + 1
    if len(selected) != expected_count:
        raise ValueError(
            f"region produced {len(selected)} k-mers; expected {expected_count}"
        )
    return selected


def _bedpe_window_sizes(dataframe):
    """Extract the window sizes represented by a loaded BEDPE dataframe."""

    for start, end in (("query_start", "query_end"), ("q_st", "q_en")):
        if start in dataframe and end in dataframe:
            sizes = dataframe[end] - dataframe[start]
            return sorted({int(size) for size in sizes if size > 0})
    return []


def _regions_from_names(names):
    """Return normalized region strings embedded in sequence names."""

    regions = []
    for name in names:
        parsed = extractRegion(name)
        if parsed:
            region = f"{parsed[0]}:{parsed[1]}-{parsed[2]}"
            if region not in regions:
                regions.append(region)
    return regions


def _annotate_bed_directions(
    bed,
    canonical_matrix,
    forward_matrix,
    window_size,
    x_offset,
    y_offset,
):
    """Add forward/reverse orientation labels to canonical BEDPE rows."""

    canonical_matrix = np.asarray(canonical_matrix)
    forward_matrix = np.asarray(forward_matrix)
    if canonical_matrix.shape != forward_matrix.shape:
        raise ValueError("Canonical and forward matrices must have matching shapes")
    if not bed:
        return bed

    header = list(bed[0])
    try:
        query_start_index = header.index("query_start")
        reference_start_index = header.index("reference_start")
    except ValueError as error:
        raise ValueError(
            "BEDPE data is missing direction coordinate columns"
        ) from error

    annotated = [tuple([*header, "direction"])]
    for row in bed[1:]:
        query_index = round((int(row[query_start_index]) - x_offset) / window_size)
        reference_index = round(
            (int(row[reference_start_index]) - y_offset) / window_size
        )
        try:
            direction = (
                "Forward"
                if forward_matrix[query_index, reference_index] > 0
                else "Reverse"
            )
        except IndexError as error:
            raise ValueError(
                "BEDPE coordinates fall outside direction matrices"
            ) from error
        annotated.append(tuple([*row, direction]))
    return annotated


def main():
    print(ASCII_ART)
    print(f"v{VERSION} \n")
    args = parse_args()
    summary_writer = PlotSummaryWriter(shlex.join(sys.argv))
    annotation_df = None
    if args.command == "interactive" and args.bed:
        try:
            annotation_df = read_annotation_beds(args.bed)
        except (OSError, ValueError) as error:
            print(f"Error reading annotation BED file(s): {error}", file=sys.stderr)
            sys.exit(2)
    if args.command == "static":
        _load_static_plotting()
    # -----------MUTUALLY EXCLUSIVE: INTERACTIVE OR STATIC MODE-----------
    if args.command == "interactive":
        print(INTERACTIVE_DEPRECATION_MESSAGE, file=sys.stderr)
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
                for matrix in matrices[i]:
                    x_start, x_end = interactive_axis_bounds(metadata[i], "x")
                    y_start, y_end = interactive_axis_bounds(metadata[i], "y")
                    x_axis = [
                        value
                        for value in np.linspace(x_start, x_end, matrix.shape[0] + 1)
                    ]
                    y_axis = [
                        value
                        for value in np.linspace(y_start, y_end, matrix.shape[1] + 1)
                    ]
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
                annotation_df,
            )
            sys.exit(0)
    elif args.command == "static":
        print(f"Running ModDotPlot in static mode\n")
        # -----------CONFIG PARSING-----------
        # TODO: Change to yml file, add readme to config folder
        if args.config:
            with open(args.config, "r") as f:
                config = json.load(f)
                _apply_static_config(args, config)

        # -----------INPUT COMMAND VALIDATION-----------
        if args.plot_direction and getattr(args, "load", None):
            print(
                "Error: --plot-direction requires FASTA input because strand "
                "orientation cannot be recovered from a BEDPE file.\n"
            )
            sys.exit(2)

        # TODO: More tests!
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
        if hasattr(args, "load") and args.load:
            if args.grid or args.grid_only:
                single_vals = []
                double_vals = []
                single_val_name = []
                double_val_name = []
                xlim_val_grid = 0
                loaded_window_sizes = []
                loaded_regions = []
            for bed in args.load:
                # If args.load is provided as input, run static mode directly from the paired-end bed file. Skip counting input k-mers.
                df = read_df_from_file(bed)

                unique_query_names = df["#query_name"].unique()
                unique_reference_names = df["reference_name"].unique()
                bed_window_sizes = _bedpe_window_sizes(df)
                bed_regions = _regions_from_names(
                    [*unique_query_names, *unique_reference_names]
                )
                if args.grid or args.grid_only:
                    loaded_window_sizes.extend(bed_window_sizes)
                    loaded_regions.extend(bed_regions)
                assert len(unique_query_names) == len(unique_reference_names)
                assert len(unique_reference_names) == 1
                self_id_scores = df[df["#query_name"] == df["reference_name"]]
                pairwise_id_scores = df[df["#query_name"] != df["reference_name"]]
                if not args.grid_only:
                    print(f"Input bed file {bed} read successfully! Creating plots: \n")
                # Create directory
                if not args.output_dir:
                    args.output_dir = os.getcwd()
                if not os.path.exists(args.output_dir):
                    os.makedirs(args.output_dir)
                if len(self_id_scores) > 1:
                    if not args.grid_only:
                        plot_files = create_plots(
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
                            xlim=args.axes_limits,
                            custom_colors=args.colors,
                            custom_breakpoints=args.breakpoints,
                            from_file=df,
                            is_pairwise=False,
                            axes_labels=args.axes_ticks,
                            axes_tick_number=args.axes_number,
                            vector_format=args.vector,
                            deraster=args.deraster,
                            annotation=args.bed,
                        )
                        summary_writer.add(
                            args.output_dir,
                            plot_files or [],
                            window_sizes=bed_window_sizes,
                            regions=bed_regions,
                            bed_file=args.bed,
                            bedpe_inputs=[bed],
                        )
                    if args.grid or args.grid_only:
                        single_vals.append(df)
                        single_val_name.append(unique_query_names[0])
                        try:
                            max_so_far = max(
                                df["query_end"].max(), df["reference_end"].max()
                            )
                        except KeyError:
                            max_so_far = max(df["q_en"].max(), df["r_en"].max())
                        xlim_val_grid = (
                            max_so_far if max_so_far > xlim_val_grid else xlim_val_grid
                        )
                # Case 2: Pairwise bed file
                if len(pairwise_id_scores) > 1:
                    if not args.grid_only:
                        # Potentially sort
                        plot_files = create_plots(
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
                            xlim=args.axes_limits,
                            custom_colors=args.colors,
                            custom_breakpoints=args.breakpoints,
                            from_file=df,
                            is_pairwise=True,
                            axes_labels=args.axes_ticks,
                            axes_tick_number=args.axes_number,
                            vector_format=args.vector,
                            deraster=args.deraster,
                            annotation=args.bed,
                        )
                        summary_writer.add(
                            args.output_dir,
                            plot_files or [],
                            window_sizes=bed_window_sizes,
                            regions=bed_regions,
                            bed_file=args.bed,
                            bedpe_inputs=[bed],
                        )
                    if args.grid or args.grid_only:
                        double_vals.append(df)
                        double_val_name.append(
                            [unique_query_names[0], unique_reference_names[0]]
                        )
            # Exit once all bed files have been iterated through
            if args.grid or args.grid_only:
                # Determine the number of sequences
                if args.axes_limits:
                    xlim_val_grid = args.axes_limits
                print(
                    f"Creating a {len(single_val_name)}x{len(single_val_name)} grid.\n"
                )
                plot_files = create_grid(
                    singles=single_vals,
                    doubles=double_vals,
                    directory=args.output_dir if args.output_dir else ".",
                    palette=args.palette,
                    palette_orientation=args.palette_orientation,
                    single_names=single_val_name,
                    double_names=double_val_name,
                    is_freq=args.bin_freq,
                    xlim=xlim_val_grid,
                    custom_colors=args.colors,
                    custom_breakpoints=args.breakpoints,
                    axes_label=args.axes_ticks,
                    is_bed=True,
                    width=args.width,
                    breaks=args.axes_ticks,
                    deraster=args.deraster,
                    vector_format=args.vector,
                    dpi=args.dpi,
                )
                summary_writer.add(
                    args.output_dir,
                    plot_files or [],
                    window_sizes=loaded_window_sizes,
                    regions=loaded_regions,
                    bed_file=args.bed,
                    bedpe_inputs=args.load,
                )
            sys.exit(0)

    # -----------INPUT SEQUENCE VALIDATION-----------
    seq_list = []
    fasta_list = args.fasta.copy()
    fasta_headers = {}
    for i in args.fasta:
        try:
            headers = getInputHeaders(i)
            fasta_headers[i] = headers

            if len(headers) > 1:
                print(f"File {i} contains multiple fasta entries.\n")

            seq_list.extend(headers)  # Add all headers to seq_list

        except Exception as e:
            print(
                f"\nUnable to open {i}. Please check it is correctly formatted or compressed...\n"
            )
            fasta_list.remove(i)

    fasta_source_by_name = {}
    for fasta_path, headers in fasta_headers.items():
        for header in headers:
            parsed_header = extractRegion(header)
            base_header = parsed_header[0] if parsed_header else header
            fasta_source_by_name.setdefault(base_header, fasta_path)

    try:
        region_by_name = _parse_region_arguments(
            getattr(args, "region", None), seq_list
        )
    except ValueError as error:
        print(f"Error: {error}.\n")
        sys.exit(2)

    # -----------LOAD SEQUENCES INTO MEMORY-----------
    kmer_list = []
    for i in fasta_list:
        if args.forward:
            kmer_list.append(
                readKmersFromFile(
                    i,
                    args.kmer,
                    False,
                    True,
                    args.ambiguous,
                    region_by_name,
                    fasta_headers[i],
                )
            )
        else:
            kmer_list.append(
                readKmersFromFile(
                    i,
                    args.kmer,
                    False,
                    False,
                    args.ambiguous,
                    region_by_name,
                    fasta_headers[i],
                )
            )
    k_list = [item for sublist in kmer_list for item in sublist]

    # Direction plots need both canonical and forward-only hashes.  Load the
    # opposite representation only when requested so normal runs retain their
    # existing memory footprint.
    direction_k_list = None
    if args.command == "static" and args.plot_direction:
        alternate_kmer_list = [
            readKmersFromFile(
                path,
                args.kmer,
                False,
                not args.forward,
                args.ambiguous,
                region_by_name,
                fasta_headers[path],
            )
            for path in fasta_list
        ]
        direction_k_list = [item for sublist in alternate_kmer_list for item in sublist]
        if len(direction_k_list) != len(k_list):
            raise ValueError(
                "Canonical and forward-only FASTA parsing produced different sequence counts"
            )
    # Throw error if compare only selected with one sequence.
    if len(k_list) < 2 and args.compare_only:
        print(
            f"Error: Can't create a comparative plot with only one sequence. Please re-run without --compare-only."
        )
        sys.exit(2)

    # -----------LAUNCH INTERACTIVE MODE-----------
    if args.command == "interactive":
        # Use the longest sequence to size the shared interactive image pyramid.
        hgi = max(len(kmers) for kmers in k_list)
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
            print(f"Error: sequence too small for analysis.\n")
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
                        "hash_algorithm": HASH_ALGORITHM,
                        "format_version": 2,
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
                    "hash_algorithm": HASH_ALGORITHM,
                    "format_version": 2,
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
                print(
                    f"Saved matrices to {folder_path}. Thank you for using ModDotPlot!\n"
                )
                sys.exit(0)

        # Before running dash, change into intervals...
        axes = []
        for matrices_set, meta in zip(matrices, metadata):
            matrix_axes = []
            x_start, x_end = interactive_axis_bounds(meta, "x")
            y_start, y_end = interactive_axis_bounds(meta, "y")
            for matrix in matrices_set:
                x_axis = np.linspace(x_start, x_end, matrix.shape[0] + 1)
                y_axis = np.linspace(y_start, y_end, matrix.shape[1] + 1)
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
            annotation_df,
        )

    # -----------SETUP STATIC MODE-----------
    elif args.command == "static":
        # -----------SET SPARSITY VALUE-----------
        direction_rendering = args.plot_direction and (
            not args.no_plot or args.grid or args.grid_only
        )
        if args.grid or args.grid_only:
            grid_val_singles = []
            grid_val_single_names = []
            grid_window_sizes = []
            if direction_rendering:
                direction_grid_val_singles = []
                direction_grid_val_doubles = []
        if direction_k_list is None:
            new_sequences = list(zip(seq_list, k_list))
        else:
            new_sequences = list(zip(seq_list, k_list, direction_k_list))
        if args.compare_order == "size":
            sequences = sorted(new_sequences, key=lambda seq: len(seq[1]), reverse=True)
        else:
            sequences = new_sequences

        # Record exact base-coordinate bounds independently of sparse hits.
        # Renderers must not infer these bounds from BEDPE rows: thresholding
        # can remove edge windows, and a partial final window can extend past
        # the selected interval.
        selected_intervals = []
        for sequence in sequences:
            sequence_name = sequence[0]
            header_range = extractRegion(sequence_name)
            base_name = header_range[0] if header_range else sequence_name
            selected_range = region_by_name.get(base_name)
            if selected_range:
                interval_start, interval_end = selected_range[1:]
            else:
                interval_start = int(header_range[1]) if header_range else 1
                interval_end = interval_start + len(sequence[1]) + args.kmer - 2
            selected_intervals.append((interval_start, interval_end))
        grid_axis_bounds = (
            min(start for start, _end in selected_intervals),
            max(end for _start, end in selected_intervals),
        )
        sketch_cache = (
            ModimizerSketchCache(max_entries=2) if args.grid or args.grid_only else None
        )
        if len(sequences) > 6 and (args.grid or args.grid_only):
            print(
                f"Creating a large {len(sequences)}x{len(sequences)} grid; "
                "rendering may take additional time and memory.\n"
            )

        # Create output directory, if doesn't exist:
        if (args.output_dir) and not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir, exist_ok=True)
        # -----------COMPUTE SELF-IDENTITY PLOTS-----------
        if not args.compare_only:
            for i in range(len(sequences)):
                sequence_name = sequences[i][0]
                header_range = extractRegion(sequence_name)
                base_name = header_range[0] if header_range else sequence_name
                sequence_start = int(header_range[1]) if header_range else 1
                seq_range = region_by_name.get(base_name)
                matrix_sequence = sequences[i][1]
                alternate_sequence = sequences[i][2] if args.plot_direction else None

                if seq_range:
                    selected_sequence_start = seq_range[1]
                    try:
                        matrix_sequence = _slice_kmers_for_region(
                            matrix_sequence,
                            seq_range,
                            args.kmer,
                            sequence_start=selected_sequence_start,
                        )
                        if args.plot_direction:
                            alternate_sequence = _slice_kmers_for_region(
                                alternate_sequence,
                                seq_range,
                                args.kmer,
                                sequence_start=selected_sequence_start,
                            )
                    except ValueError as error:
                        print(f"Error: invalid region for {base_name}: {error}.\n")
                        sys.exit(2)
                    _, seq_start_pos, subseq_end_pos = seq_range
                    seq_name = f"{base_name}:{seq_start_pos}-{subseq_end_pos}"
                    print(f"Using region {seq_name}\n")
                else:
                    seq_start_pos = sequence_start
                    subseq_end_pos = (
                        seq_start_pos + len(matrix_sequence) + args.kmer - 2
                    )
                    seq_name = sequence_name

                plot_axis_bounds = args.axes_limits or (
                    seq_start_pos,
                    subseq_end_pos,
                )

                seq_length = len(matrix_sequence)
                win = args.window
                res = args.resolution
                if args.window:
                    # Change the resolution of each plot
                    res = math.ceil(seq_length / args.window)
                else:
                    win = math.ceil(seq_length / args.resolution)

                if win < args.modimizer:
                    args.modimizer = win
                if win < 10:
                    print(f"Error: sequence too small for analysis.\n")
                    print(
                        f"ModDotPlot requires a minimum window size of 10. Sequences less than 10Kbp will not work with ModDotPlot under normal resolution. We recommend rerunning ModDotPlot with --r {math.ceil(seq_length / 10)}.\n"
                    )
                    sys.exit(0)

                seq_sparsity = round(win / args.modimizer)
                if seq_sparsity <= args.modimizer:
                    seq_sparsity = 2 ** int(math.log2(seq_sparsity))
                else:
                    seq_sparsity = 2 ** (int(math.log2(seq_sparsity - 1)) + 1)
                expectation = round(win / seq_sparsity)

                print(f"Computing self identity matrix for {seq_name}... \n")
                # TODO: Logging here
                # print(f"\tSparsity value s: {seq_sparsity}\n")
                print(f"\tSequence length n: {seq_length + args.kmer - 1}\n")
                print(f"\tWindow size w: {win}\n")
                print(f"\tModimizer sketch size: {expectation}\n")
                print(f"\tPlot Resolution r: {res}\n")

                if sketch_cache is None:
                    self_mat = createSelfMatrix(
                        seq_length,
                        matrix_sequence,
                        win,
                        seq_sparsity,
                        args.delta,
                        args.kmer,
                        args.identity,
                        args.ambiguous,
                        expectation,
                    )
                else:
                    source_region = (seq_range[1], seq_range[2]) if seq_range else None
                    prepared_self = sketch_cache.get_or_prepare(
                        (i, source_region),
                        seq_length,
                        matrix_sequence,
                        win,
                        seq_sparsity,
                        args.delta,
                        args.kmer,
                        args.ambiguous,
                        expectation,
                    )
                    self_mat = create_self_matrix_from_sketches(
                        prepared_self, args.kmer, args.identity, args.ambiguous
                    )
                    # The cache owns the reusable reference. Keeping this loop
                    # local alive can pin an evicted sketch until all grid
                    # calculations finish.
                    del prepared_self
                direction_self_mat = None
                if direction_rendering:
                    direction_self_mat = createSelfMatrix(
                        seq_length,
                        alternate_sequence,
                        win,
                        seq_sparsity,
                        args.delta,
                        args.kmer,
                        args.identity,
                        args.ambiguous,
                        expectation,
                    )
                bed = convertMatrixToBed(
                    self_mat,
                    win,
                    args.identity,
                    seq_name,
                    seq_name,
                    True,
                    seq_start_pos,
                    seq_start_pos,
                    subseq_end_pos,
                    subseq_end_pos,
                )
                plot_bed = bed
                if direction_rendering:
                    if args.forward:
                        canonical_matrix = direction_self_mat
                        forward_matrix = self_mat
                        canonical_bed = convertMatrixToBed(
                            canonical_matrix,
                            win,
                            args.identity,
                            seq_name,
                            seq_name,
                            True,
                            seq_start_pos,
                            seq_start_pos,
                            subseq_end_pos,
                            subseq_end_pos,
                        )
                    else:
                        canonical_matrix = self_mat
                        forward_matrix = direction_self_mat
                        canonical_bed = bed
                    plot_bed = _annotate_bed_directions(
                        canonical_bed,
                        canonical_matrix,
                        forward_matrix,
                        win,
                        seq_start_pos,
                        seq_start_pos,
                    )
                if args.grid or args.grid_only:
                    grid_val_singles.append(bed)
                    grid_val_single_names.append(seq_name)
                    grid_window_sizes.append(win)
                    if direction_rendering:
                        direction_grid_val_singles.append(plot_bed)

                if args.cooler:
                    try:
                        cooler_path = "."
                        if not args.output_dir:
                            cooler_path = os.path.join(cooler_path, seq_name)
                        else:
                            cooler_path = os.path.join(args.output_dir, seq_name)
                        os.makedirs(cooler_path, exist_ok=True)
                        cooler_output = os.path.join(cooler_path, seq_name + ".cooler")
                        convertMatrixToCool(
                            matrix=self_mat,
                            window_size=win,
                            id_threshold=args.identity,
                            x_name=seq_name,
                            y_name=seq_name,
                            self_identity=True,
                            x_offset=seq_start_pos,
                            y_offset=seq_start_pos,
                            chromsizes=seq_length,
                            output_cool=cooler_output,
                        )
                        print(
                            f"Saved self-identity matrix as a cooler file to {cooler_output}\n"
                        )
                    except Exception as e:
                        print(f"Error creating cooler file: {e}")

                bedpe_path = os.path.join(args.output_dir or ".", seq_name)
                if (not args.no_bedpe) or ((not args.no_plot) and (not args.grid_only)):
                    os.makedirs(bedpe_path, exist_ok=True)

                if not args.no_bedpe:
                    # Log saving bed file
                    bedfile_output = os.path.join(bedpe_path, seq_name + ".bedpe")

                    with open(bedfile_output, "w") as bedfile:
                        for row in bed:
                            bedfile.write("\t".join(map(str, row)) + "\n")
                    print(
                        f"Saved self-identity matrix as a paired-end bed file to {bedfile_output}\n"
                    )

                if (not args.no_plot) and (not args.grid_only):
                    plot_files = create_plots(
                        sdf=[bed],
                        directory=bedpe_path,
                        name_x=seq_name,
                        name_y=seq_name,
                        palette=args.palette,
                        palette_orientation=args.palette_orientation,
                        no_hist=args.no_hist,
                        width=args.width,
                        dpi=args.dpi,
                        is_freq=args.bin_freq,
                        xlim=plot_axis_bounds,
                        custom_colors=args.colors,
                        custom_breakpoints=args.breakpoints,
                        from_file=None,
                        is_pairwise=False,
                        axes_labels=args.axes_ticks,
                        axes_tick_number=args.axes_number,
                        vector_format=args.vector,
                        deraster=args.deraster,
                        annotation=args.bed,
                    )
                    self_region = (
                        [f"{base_name}:{seq_range[1]}-{seq_range[2]}"]
                        if seq_range
                        else _regions_from_names([sequence_name])
                    )
                    summary_writer.add(
                        bedpe_path,
                        plot_files or [],
                        fasta_files=[fasta_source_by_name[base_name]],
                        window_sizes=[win],
                        regions=self_region,
                        bed_file=args.bed,
                    )
                    if direction_rendering:
                        direction_directory = os.path.join(bedpe_path, "directionality")
                        direction_files = create_plots(
                            sdf=[plot_bed],
                            directory=direction_directory,
                            name_x=seq_name,
                            name_y=seq_name,
                            palette=args.palette,
                            palette_orientation=args.palette_orientation,
                            no_hist=args.no_hist,
                            width=args.width,
                            dpi=args.dpi,
                            is_freq=args.bin_freq,
                            xlim=plot_axis_bounds,
                            custom_colors=args.colors,
                            custom_breakpoints=args.breakpoints,
                            from_file=None,
                            is_pairwise=False,
                            axes_labels=args.axes_ticks,
                            axes_tick_number=args.axes_number,
                            vector_format=args.vector,
                            deraster=args.deraster,
                            annotation=None,
                        )
                        summary_writer.add(
                            direction_directory,
                            direction_files or [],
                            fasta_files=[fasta_source_by_name[base_name]],
                            window_sizes=[win],
                            regions=self_region,
                            bed_file=args.bed,
                        )

        # -----------COMPUTE COMPARATIVE PLOTS-----------
        # TODO: Optimize computations so that largest sequence doesn't need to be redone all the time
        if (args.compare or args.compare_only or args.grid or args.grid_only) and len(
            sequences
        ) > 1:
            # Set window size to args.window. Otherwise, set it to n/resolution

            if args.grid or args.grid_only:
                grid_val_doubles = []
                grid_val_double_names = []
                xlim_val_grid = args.axes_limits or grid_axis_bounds

            for i in range(len(sequences)):
                for j in range(i + 1, len(sequences)):
                    # Larger = x, smaller = y. This is pre-sorted earlier.
                    larger_seq = sequences[i][1]
                    smaller_seq = sequences[j][1]
                    larger_direction_seq = (
                        sequences[i][2] if args.plot_direction else None
                    )
                    smaller_direction_seq = (
                        sequences[j][2] if args.plot_direction else None
                    )
                    larger_sequence_name = sequences[i][0]
                    smaller_sequence_name = sequences[j][0]
                    larger_header_range = extractRegion(larger_sequence_name)
                    smaller_header_range = extractRegion(smaller_sequence_name)
                    larger_base_name = (
                        larger_header_range[0]
                        if larger_header_range
                        else larger_sequence_name
                    )
                    smaller_base_name = (
                        smaller_header_range[0]
                        if smaller_header_range
                        else smaller_sequence_name
                    )
                    larger_sequence_start = (
                        int(larger_header_range[1]) if larger_header_range else 1
                    )
                    smaller_sequence_start = (
                        int(smaller_header_range[1]) if smaller_header_range else 1
                    )
                    larger_seq_range = region_by_name.get(larger_base_name)
                    smaller_seq_range = region_by_name.get(smaller_base_name)
                    larger_subseq = larger_seq
                    smaller_subseq = smaller_seq
                    larger_direction_subseq = larger_direction_seq
                    smaller_direction_subseq = smaller_direction_seq
                    larger_seq_start_pos = larger_sequence_start
                    smaller_seq_start_pos = smaller_sequence_start
                    larger_seq_end_pos = (
                        larger_seq_start_pos + len(larger_subseq) + args.kmer - 2
                    )
                    smaller_seq_end_pos = (
                        smaller_seq_start_pos + len(smaller_subseq) + args.kmer - 2
                    )
                    larger_seq_name = larger_sequence_name
                    smaller_seq_name = smaller_sequence_name

                    try:
                        if larger_seq_range:
                            selected_larger_start = larger_seq_range[1]
                            larger_subseq = _slice_kmers_for_region(
                                larger_seq,
                                larger_seq_range,
                                args.kmer,
                                sequence_start=selected_larger_start,
                            )
                            if args.plot_direction:
                                larger_direction_subseq = _slice_kmers_for_region(
                                    larger_direction_seq,
                                    larger_seq_range,
                                    args.kmer,
                                    sequence_start=selected_larger_start,
                                )
                            _, larger_seq_start_pos, larger_end = larger_seq_range
                            larger_seq_end_pos = larger_end
                            larger_seq_name = f"{larger_base_name}:{larger_seq_start_pos}-{larger_end}"
                            print(f"Using region {larger_seq_name}\n")

                        if smaller_seq_range:
                            selected_smaller_start = smaller_seq_range[1]
                            smaller_subseq = _slice_kmers_for_region(
                                smaller_seq,
                                smaller_seq_range,
                                args.kmer,
                                sequence_start=selected_smaller_start,
                            )
                            if args.plot_direction:
                                smaller_direction_subseq = _slice_kmers_for_region(
                                    smaller_direction_seq,
                                    smaller_seq_range,
                                    args.kmer,
                                    sequence_start=selected_smaller_start,
                                )
                            _, smaller_seq_start_pos, smaller_end = smaller_seq_range
                            smaller_seq_end_pos = smaller_end
                            smaller_seq_name = f"{smaller_base_name}:{smaller_seq_start_pos}-{smaller_end}"
                            print(f"Using region {smaller_seq_name}\n")
                    except ValueError as error:
                        print(f"Error: invalid comparison region: {error}.\n")
                        sys.exit(2)

                    larger_length = len(larger_subseq)
                    smaller_length = len(smaller_subseq)
                    pair_axis_bounds = args.axes_limits or (
                        min(larger_seq_start_pos, smaller_seq_start_pos),
                        max(larger_seq_end_pos, smaller_seq_end_pos),
                    )

                    win = args.window
                    res = args.resolution
                    if args.window:
                        res = math.ceil(smaller_length / args.window)
                    else:
                        win = math.ceil(smaller_length / args.resolution)
                    if win < args.modimizer:
                        args.modimizer = win

                    seq_sparsity = round(win / args.modimizer)
                    if seq_sparsity <= args.modimizer:
                        seq_sparsity = 2 ** int(math.log2(seq_sparsity))
                    else:
                        seq_sparsity = 2 ** (int(math.log2(seq_sparsity - 1)) + 1)
                    expectation = round(win / seq_sparsity)
                    print(
                        f"Computing pairwise identity matrix for {larger_seq_name} and {smaller_seq_name}... \n"
                    )
                    # TODO: Logging here
                    print(
                        f"\tSequence length {larger_seq_name}: {larger_length + args.kmer - 1}\n"
                    )
                    print(
                        f"\tSequence length {smaller_seq_name}: {smaller_length + args.kmer - 1}\n"
                    )
                    print(f"\tWindow size w: {win}\n")
                    print(f"\tModimizer sketch size: {expectation}\n")
                    print(f"\tPlot Resolution r: {res}\n")

                    if sketch_cache is None:
                        pair_mat = createPairwiseMatrix(
                            smaller_length,
                            larger_length,
                            smaller_subseq,
                            larger_subseq,
                            win,
                            seq_sparsity,
                            args.delta,
                            args.kmer,
                            args.identity,
                            args.ambiguous,
                            expectation,
                        )
                    else:
                        smaller_source_region = (
                            (smaller_seq_range[1], smaller_seq_range[2])
                            if smaller_seq_range
                            else None
                        )
                        larger_source_region = (
                            (larger_seq_range[1], larger_seq_range[2])
                            if larger_seq_range
                            else None
                        )
                        prepared_smaller = sketch_cache.get_or_prepare(
                            (j, smaller_source_region),
                            smaller_length,
                            smaller_subseq,
                            win,
                            seq_sparsity,
                            args.delta,
                            args.kmer,
                            args.ambiguous,
                            expectation,
                        )
                        prepared_larger = sketch_cache.get_or_prepare(
                            (i, larger_source_region),
                            larger_length,
                            larger_subseq,
                            win,
                            seq_sparsity,
                            args.delta,
                            args.kmer,
                            args.ambiguous,
                            expectation,
                        )
                        pair_mat = create_pairwise_matrix_from_sketches(
                            prepared_smaller,
                            prepared_larger,
                            args.identity,
                            args.kmer,
                        )
                        # Avoid retaining entries after the bounded cache
                        # evicts or clears them.
                        del prepared_smaller, prepared_larger
                    direction_pair_mat = None
                    if direction_rendering:
                        direction_pair_mat = createPairwiseMatrix(
                            smaller_length,
                            larger_length,
                            smaller_direction_subseq,
                            larger_direction_subseq,
                            win,
                            seq_sparsity,
                            args.delta,
                            args.kmer,
                            args.identity,
                            args.ambiguous,
                            expectation,
                        )
                        canonical_pair_mat = (
                            direction_pair_mat if args.forward else pair_mat
                        )
                    else:
                        canonical_pair_mat = pair_mat
                    # Throw error if the matrix is empty
                    if np.all(canonical_pair_mat == 0) and not (
                        args.grid or args.grid_only
                    ):
                        print(
                            f"The pairwise identity matrix for {sequences[i][0]} and {sequences[j][0]} is empty. Skipping.\n"
                        )
                    else:
                        if args.cooler:
                            try:
                                cooler_path = "."
                                if not args.output_dir:
                                    cooler_path = os.path.join(
                                        cooler_path,
                                        f"{larger_seq_name}_{smaller_seq_name}",
                                    )
                                else:
                                    cooler_path = os.path.join(
                                        args.output_dir,
                                        f"{larger_seq_name}_{smaller_seq_name}",
                                    )
                                os.makedirs(cooler_path, exist_ok=True)
                                cooler_output = os.path.join(
                                    cooler_path,
                                    f"{larger_seq_name}_{smaller_seq_name}.cooler",
                                )
                                convertMatrixToCool(
                                    matrix=pair_mat,
                                    window_size=win,
                                    id_threshold=args.identity,
                                    x_name=larger_seq_name,
                                    y_name=smaller_seq_name,
                                    self_identity=False,
                                    x_offset=larger_seq_start_pos,
                                    y_offset=smaller_seq_start_pos,
                                    chromsizes=larger_length,
                                    output_cool=cooler_output,
                                )
                                print(
                                    f"Saved comparative matrix as a cooler file to {cooler_output}\n"
                                )
                            except Exception as e:
                                print(f"Error creating pairwise cooler file: {e}")
                        bed = convertMatrixToBed(
                            pair_mat,
                            win,
                            args.identity,
                            # check if this is correct
                            larger_seq_name,
                            smaller_seq_name,
                            False,
                            larger_seq_start_pos,
                            smaller_seq_start_pos,
                            larger_seq_end_pos,
                            smaller_seq_end_pos,
                        )
                        plot_bed = bed
                        if direction_rendering:
                            if args.forward:
                                canonical_matrix = direction_pair_mat
                                forward_matrix = pair_mat
                                canonical_bed = convertMatrixToBed(
                                    canonical_matrix,
                                    win,
                                    args.identity,
                                    larger_seq_name,
                                    smaller_seq_name,
                                    False,
                                    larger_seq_start_pos,
                                    smaller_seq_start_pos,
                                    larger_seq_end_pos,
                                    smaller_seq_end_pos,
                                )
                            else:
                                canonical_matrix = pair_mat
                                forward_matrix = direction_pair_mat
                                canonical_bed = bed
                            plot_bed = _annotate_bed_directions(
                                canonical_bed,
                                canonical_matrix,
                                forward_matrix,
                                win,
                                larger_seq_start_pos,
                                smaller_seq_start_pos,
                            )
                        if args.grid or args.grid_only:
                            grid_val_doubles.append(bed)
                            grid_val_double_names.append(
                                [larger_seq_name, smaller_seq_name]
                            )
                            grid_window_sizes.append(win)
                            if direction_rendering:
                                direction_grid_val_doubles.append(plot_bed)
                        bedfile_prefix = larger_seq_name + "_" + smaller_seq_name
                        bedpe_path = os.path.join(
                            args.output_dir or ".", bedfile_prefix
                        )
                        if (not args.no_bedpe) or (
                            (not args.no_plot) and (not args.grid_only)
                        ):
                            os.makedirs(bedpe_path, exist_ok=True)

                        if not args.no_bedpe:
                            # Log saving bed file
                            bedfile_output = os.path.join(
                                bedpe_path, bedfile_prefix + "_COMPARE.bedpe"
                            )
                            with open(bedfile_output, "w") as bedfile:
                                for row in bed:
                                    bedfile.write("\t".join(map(str, row)) + "\n")
                            print(
                                f"Saved comparative matrix as a paired-end bed file to {bedfile_output}\n"
                            )

                        if (not args.no_plot) and (not args.grid_only):
                            plot_files = create_plots(
                                sdf=[bed],
                                directory=bedpe_path,
                                name_x=larger_seq_name,
                                name_y=smaller_seq_name,
                                palette=args.palette,
                                palette_orientation=args.palette_orientation,
                                no_hist=args.no_hist,
                                width=args.width,
                                dpi=args.dpi,
                                is_freq=args.bin_freq,
                                xlim=pair_axis_bounds,
                                custom_colors=args.colors,
                                custom_breakpoints=args.breakpoints,
                                from_file=None,
                                is_pairwise=True,
                                axes_labels=args.axes_ticks,
                                axes_tick_number=args.axes_number,
                                vector_format=args.vector,
                                deraster=args.deraster,
                                annotation=args.bed,
                            )
                            pair_regions = []
                            if larger_seq_range:
                                pair_regions.append(
                                    f"{larger_base_name}:{larger_seq_range[1]}-"
                                    f"{larger_seq_range[2]}"
                                )
                            else:
                                pair_regions.extend(
                                    _regions_from_names([larger_sequence_name])
                                )
                            if smaller_seq_range:
                                pair_regions.append(
                                    f"{smaller_base_name}:{smaller_seq_range[1]}-"
                                    f"{smaller_seq_range[2]}"
                                )
                            else:
                                pair_regions.extend(
                                    _regions_from_names([smaller_sequence_name])
                                )
                            pair_fasta_files = [
                                fasta_source_by_name[name]
                                for name in (larger_base_name, smaller_base_name)
                            ]
                            summary_writer.add(
                                bedpe_path,
                                plot_files or [],
                                fasta_files=pair_fasta_files,
                                window_sizes=[win],
                                regions=pair_regions,
                                bed_file=args.bed,
                            )
                            if direction_rendering:
                                direction_directory = os.path.join(
                                    bedpe_path, "directionality"
                                )
                                direction_files = create_plots(
                                    sdf=[plot_bed],
                                    directory=direction_directory,
                                    name_x=larger_seq_name,
                                    name_y=smaller_seq_name,
                                    palette=args.palette,
                                    palette_orientation=args.palette_orientation,
                                    no_hist=args.no_hist,
                                    width=args.width,
                                    dpi=args.dpi,
                                    is_freq=args.bin_freq,
                                    xlim=pair_axis_bounds,
                                    custom_colors=args.colors,
                                    custom_breakpoints=args.breakpoints,
                                    from_file=None,
                                    is_pairwise=True,
                                    axes_labels=args.axes_ticks,
                                    axes_tick_number=args.axes_number,
                                    vector_format=args.vector,
                                    deraster=args.deraster,
                                    annotation=None,
                                )
                                summary_writer.add(
                                    direction_directory,
                                    direction_files or [],
                                    fasta_files=pair_fasta_files,
                                    window_sizes=[win],
                                    regions=pair_regions,
                                    bed_file=args.bed,
                                )

            if sketch_cache is not None:
                sketch_cache.clear()

            if args.grid or args.grid_only:
                if args.axes_limits:
                    xlim_val_grid = args.axes_limits
                print(f"Creating a {len(sequences)}x{len(sequences)} grid.\n")
                plot_files = create_grid(
                    singles=grid_val_singles,
                    doubles=grid_val_doubles,
                    directory=args.output_dir if args.output_dir else ".",
                    palette=args.palette,
                    palette_orientation=args.palette_orientation,
                    single_names=grid_val_single_names,
                    double_names=grid_val_double_names,
                    is_freq=args.bin_freq,
                    xlim=xlim_val_grid,
                    custom_colors=args.colors,
                    custom_breakpoints=args.breakpoints,
                    axes_label=args.axes_ticks,
                    is_bed=False,
                    width=args.width,
                    breaks=args.axes_ticks,
                    deraster=args.deraster,
                    vector_format=args.vector,
                    dpi=args.dpi,
                )
                grid_directory = args.output_dir if args.output_dir else "."
                grid_regions = [
                    f"{name}:{region[1]}-{region[2]}"
                    for name, region in region_by_name.items()
                ]
                for embedded_region in _regions_from_names(
                    [sequence[0] for sequence in sequences]
                ):
                    if embedded_region not in grid_regions:
                        grid_regions.append(embedded_region)
                summary_writer.add(
                    grid_directory,
                    plot_files or [],
                    fasta_files=fasta_list,
                    window_sizes=grid_window_sizes,
                    regions=grid_regions,
                    bed_file=args.bed,
                )
                if direction_rendering:
                    direction_directory = os.path.join(grid_directory, "directionality")
                    direction_files = create_grid(
                        singles=direction_grid_val_singles,
                        doubles=direction_grid_val_doubles,
                        directory=direction_directory,
                        palette=args.palette,
                        palette_orientation=args.palette_orientation,
                        single_names=grid_val_single_names,
                        double_names=grid_val_double_names,
                        is_freq=args.bin_freq,
                        xlim=xlim_val_grid,
                        custom_colors=args.colors,
                        custom_breakpoints=args.breakpoints,
                        axes_label=args.axes_ticks,
                        is_bed=False,
                        width=args.width,
                        breaks=args.axes_ticks,
                        deraster=args.deraster,
                        vector_format=args.vector,
                        dpi=args.dpi,
                    )
                    summary_writer.add(
                        direction_directory,
                        direction_files or [],
                        fasta_files=fasta_list,
                        window_sizes=grid_window_sizes,
                        regions=grid_regions,
                        bed_file=args.bed,
                    )


if __name__ == "__main__":
    main()
