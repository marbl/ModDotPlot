![](images/logo.png)
---
[![PyPI](https://img.shields.io/pypi/v/ModDotPlot?color=blue&label=PyPI)](https://pypi.org/project/ModDotPlot/)
[![CI](https://github.com/marbl/ModDotPlot/actions/workflows/ci.yml/badge.svg)](https://github.com/marbl/ModDotPlot/actions/workflows/ci.yml)

- [](#)
- [Cite](#cite)
- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
  - [Static Mode](#static-mode)
  - [Interactive Mode](#interactive-mode)
  - [Standard arguments](#standard-arguments)
  - [Static Mode Commands](#static-mode-commands)
    - [Input/Output \& Formatting Commands](#inputoutput--formatting-commands)
    - [Plot Customization Commands](#plot-customization-commands)
  - [Sample run - Static Plots](#sample-run---static-plots)
    - [Using a config file](#using-a-config-file)
    - [Adding custom bed file annotations](#adding-custom-bed-file-annotations)
    - [Comparing two sequences](#comparing-two-sequences)
  - [Interactive Mode Commands](#interactive-mode-commands)
  - [Sample run - Interactive Mode](#sample-run---interactive-mode)
  - [Sample run - Port Forwarding](#sample-run---port-forwarding)
- [Questions](#questions)
- [Known Issues](#known-issues)


## Cite

Alexander P Sweeten, Michael C Schatz, Adam M Phillippy, ModDotPlot—rapid and interactive visualization of tandem repeats, Bioinformatics, Volume 40, Issue 8, August 2024, btae493, [https://doi.org/10.1093/bioinformatics/btae493](https://doi.org/10.1093/bioinformatics/btae493)

If you use ModDotPlot for your research, please cite our software!

---

## About

_ModDotPlot_ is a dot plot visualization tool designed to be used at scale, both for smaller sequences and whole genomes. _ModDotPlot_ is the spiritual successor to [StainedGlass](https://mrvollger.github.io/StainedGlass/). The core algorithm breaks an input sequence down into intervals of sketched *k*-mers called **mod**imizers. This enables the rapid approximation of the Average Nucleotide Identity between combinations of intervals! 

Version 1.0.0 uses a bundled [ntHash2](https://github.com/BirolLab/ntHash) implementation for k-mer hashing, replacing the previous `mmh3` runtime dependency. Hash values and exact sketches therefore differ from pre-1.0 releases; regenerate data instead of mixing sketches produced by the two algorithms. Previously saved interactive matrices remain loadable because they contain completed matrices rather than raw hashes.

FASTA parsing and static BED annotation rendering are also built into ModDotPlot in version 1.0.0, replacing the previous `pysam` and `pyGenomeTracks` runtime dependencies. Plain FASTA, gzip-compressed FASTA, and BGZF-compressed FASTA inputs remain supported, and static annotations produce both PNG and the selected SVG, PDF, or PostScript vector format.

Static triangle plots, annotation layouts, and multi-sequence grids are now composed directly with Matplotlib. This replaces the previous `CairoSVG`, `svgutils`, and `patchworklib` image-conversion and SVG-composition dependencies while retaining raster and vector output formats.

![](images/demo.gif)

If you're interested in learning more about _ModDotPlot_ and how to visualize tandem repeats, we have an in-depth [YouTube video tutorial](https://www.youtube.com/watch?v=_7sQaljB_ys&t=2321s&pp=ygUXYWxleCBzd2VldGVuIG1vZGRvdHBsb3Q%3D) hosted by the [BioDiversity Genomics Academy](https://thebgacademy.org).

--- 

## Installation

_ModDotPlot_ can be installed by running `pip install moddotplot`. Version 1.0.0 supports Python 3.11 through 3.14 and uses the current Matplotlib 3.11 and Plotnine 0.15 release lines. Alternatively, you can download the current release from GitHub by using:

```
git clone https://github.com/marbl/ModDotPlot.git
cd ModDotPlot
```

Although optional, it's recommended to setup a virtual environment before using _ModDotPlot_:

```
python -m venv venv
source venv/bin/activate
```

Once activated, you can install the required dependencies:

```
python -m pip install .
```

Finally, confirm that the installation was installed correctly and that your version is up to date by running `moddotplot -h`:
```       
  __  __           _   _____        _     _____  _       _   
 |  \/  |         | | |  __ \      | |   |  __ \| |     | |  
 | \  / | ___   __| | | |  | | ___ | |_  | |__) | | ___ | |_ 
 | |\/| |/ _ \ / _` | | |  | |/ _ \| __| |  ___/| |/ _ \| __|
 | |  | | (_) | (_| | | |__| | (_) | |_  | |    | | (_) | |_ 
 |_|  |_|\___/ \__,_| |_____/ \___/ \__| |_|    |_|\___/ \__|

 v1.0.0

usage: moddotplot [-h] [{static,interactive}] ...

ModDotPlot: Visualization of Tandem Repeats

positional arguments:
  {static,interactive}  Choose mode; static is used when omitted
    static              Static mode commands (default)
    interactive         Interactive mode commands (deprecated; explicit use only)

options:
  -h, --help            show this help message and exit
  ```

Note that running `moddotplot -h` might take a while at first! This is because the Python interpreter is compiling source code into the __pycache__ directory. Subsequent runs will use the pre-compiled code and load much faster!

--- 

## Usage

_ModDotPlot_ runs in `static` mode by default. The `static` subcommand remains
available for compatibility and clarity, so these forms are equivalent:

```
moddotplot -f sequence.fa <ARGS>
moddotplot static -f sequence.fa <ARGS>
```

### Static Mode

```
moddotplot static <ARGS>
```

Running _ModDotPlot_ in static mode quickly create plots under the specified output directory `-o`. By default, running _ModDotPlot_ in static mode this will produce the following files:

- A paired-end bed file `.bedpe`, containing intervals alongside their corresponding identity estimates.
- A self-identity dotplot for each sequence, as both an upper triangle matrix `_TRI` and full matrix `_FULL` representation.
- A histogram of identity values for each sequence.

![](images/moddotplot_output.png)

Plots and histograms are output as both rasterized `.png` images and vector graphics (default: `.svg`). [Plotnine](https://plotnine.org/) provides the primary plotting interface, while Matplotlib directly renders triangle plots, annotation layouts, multi-sequence grids, and each requested output format. Grid axes state their genomic unit (Kbp, Mbp, or Gbp). Plot text uses Helvetica by default with an automatic DejaVu Sans fallback if Helvetica cannot render a glyph.

Every directory containing generated static plots also receives a `plot_summary.txt` reproducibility record. It lists the creation time, absolute plot and input paths, window size, any selected region or annotation BED file, and the exact command used for the run.

_ModDotPlot_ supports highly customizable plotting features in static mode. See [static mode commands](#static-mode-commands) for a complete list of features.


### Interactive Mode

```
moddotplot interactive <ARGS>
```

Interactive mode is deprecated and maintenance-only. It remains available, but
will not receive new features. It runs only when the `interactive` subcommand is
explicitly provided.

Running _ModDotPlot_ in interactive mode will launch a [Dash application](https://plotly.com/dash/) on your machine's localhost. Open any web browser and go to `http://127.0.0.1:<PORT_NUMBER>` to view the interactive plot (this should happen automatically, but depending on your environment you might need to copy and paste this URL into your web browser). Running `Ctrl+C` on the command line will exit the Dash application. The default port number used by Dash is `8050`, but this can be customized using the `--port` command (see [interactive mode commands](#interactive-mode-commands) for further info, and [Sample run - Port Forwarding](#sample-run---port-forwarding) for tips on running interactive mode on an HPC environment).

--- 

### Standard arguments

The following arguments are the same in both interactive and static mode:

`-f / --fasta <file>`

Fasta files to input. Multifasta files are accepted. Interactive mode will only support a maximum of two sequences at a time.

`-b / --bed <.bed file(s)>`

Input BED3-BED9 annotation file used for dotplot annotation (this is not the paired-end BEDPE file produced by ModDotPlot). The BED chromosome field must match a FASTA header, excluding any trailing `:start-end` region suffix. Static mode accepts one BED file and produces an annotation track plus annotated triangle output. Interactive mode accepts one or more BED files, combines their matching intervals, and displays a collapsed track beneath the x axis. Comparative interactive plots also display a track beside the y axis when that sequence has matching annotations. BED `itemRgb` colors are used when present.

`-k / --kmer <int>`

K-mer size to use. This should be large enough to distinguish unique k-mers with enough specificity, but not too large that sensitivity is removed. Default: 21.

`-o / --output-dir <string>`

Name of output directory for bed file & plots. Default is current working directory.

`-id / --identity <int>`

Minimum sequence identity cutoff threshold. Default is 86. While it is possible to go as low as 50% sequence identity, anything below 80% is not recommended. 

`--delta <float>`

Each partition includes a fraction of the adjacent windows' k-mers when estimating identity. This recovers repetitive matches that straddle different window boundaries. The default is 0.5, and the accepted range is between 0 and 1; values greater than 0.5 are not recommended. Set this to 0 only when strictly core-local comparisons are desired.

`-m / --modimizer <int>`

Modimizer sketch size. Must be lower than window size `w`. A lower sketch size means less k-mers to compare (and faster runtime), at the expense of lower accuracy. Recommended to be kept >= 1000.

`--forward <bool>`

Use forward k-mers only, instead of the default of canonical k-mers. Warning: this will give strand specific output.

`-r / --resolution <int>`

Dotplot resolution. This corresponds to the number of windows each input sequence is partitioned into. Default is 1000. Overrides the `--window` parameter.

`--compare <bool>`

If set when 2 or more sequences are input into ModDotPlot, this will show an A vs. B style plot, in addition to a self-identity plot. Note that interactive mode currently only supports a maximum of two sequences. If more than two sequences are input, only the first two will be shown.

`--compare-only <bool>`

If set when 2 or more sequences are input into ModDotPlot, this will show an A vs. B style plot, without showing self-identity plots.

`--ambiguous <bool>`

By default, every k-mer window containing a non-ACGTU character is excluded from identity estimation without changing its genomic position. This produces gaps through regions containing ambiguous IUPAC bases. To include deterministic hashes for those windows, set the `--ambiguous` flag in either interactive or static mode.

--- 

### Static Mode Commands

#### Input/Output & Formatting Commands

`-l / --load <.bedpe file>`

Create a plot from a previously computed pairwise bed file. Skips Average Nucleotide Identity computation. Used instead of `-f/--fasta`. Will only accept paired-end bed files produced by ModDotPlot. 

`-c / --config <.json file>`

Run moddotplot static with a config file instead of command line args. Example syntax in `config/config.json`. Recommended when creating a really customized plot. Used instead of -f/--fasta.

`--cooler <bool>`

If set, will output a matrix as a cooler file for each input sequence, in addition to a bedpe file.

`--no-bedpe <bool>`

Skip output of bed file.

`--no-hist <bool>`

Skip output of histogram legend.

`--no-plot <bool>`

Save .bedpe to file, but skip rendering of plots.

`--width <float>`

Adjust the output figure width. For a grid this is the width of the complete grid, not each cell. Default is 9 inches.

`--dpi <int>`

Image resolution in dots per inch (not to be confused with dotplot resolution). Default is `300`.

`--vector <str>` 

Vectorized image format to output to. Must be one of ["svg", "pdf", "ps"]. Default: `svg`

`--deraster <bool>`

By default, vectorized outputs rasterize the actual dotplot (not the axis). This is done to save space, as a high-resolution dotplot can be extremely space inefficient and prevent use of image manipulation software. This plot rasterization can be removed using this flag. 

#### Plot Customization Commands

`-w / --window <int>`

Window size. Unlike interactive mode, only one matrix will be created, so this represents the *only* window size. Default is set to `n/1000` (eg. 3000bp for a 3Mbp sequence). 

`--region <list of strs>`

Plot only the requested 1-based, inclusive range for each named sequence. Syntax is `FASTA_ID:start-end`; the identifier must exactly match the FASTA header's first whitespace-delimited token. Supply one value per sequence when every grid row and column should be cropped, for example `--region sample.hap1:1-4000000 sample.hap2:1-4000000`. Region limits apply to self plots, pairwise plots, BEDPE coordinates, and grid axes.

`--palette <str>`

List of accepted palettes can be found [here](https://jiffyclub.github.io/palettable/colorbrewer/). Palettes are segregated into 3 types: _Diverging_, _Qualitative_, and _Sequential_. Syntax is the name of the palette, followed by an underscore and the number of colors, eg. `OrRd_8`. Default is  `Spectral_11`.

`--breakpoints <list of ints>`

Add custom identity threshold breakpoints. Note that the number of breakpoints must be equal to the number of colors + 1, otherwise an error will occur. 

`--palette-orientation <bool>`

Flip sequential order of color palette. Set to `-` by default for divergent palettes. 

`--colors <list of hexcodes>` (legacy alias: `--color`)

List of custom colors in hexcode format can be entered sequentially, mapped from low to high identity. 

`--plot-direction <bool>`

With FASTA input, retain the standard ANI-colored plots and additionally create a `directionality` subfolder. Direction plots use blue for same-orientation matches and pink for reverse-orientation matches, with darker shades representing stronger ANI. Self-comparisons are named `_DIRECTION_FULL`, `_DIRECTION_TRI`, and `_DIRECTION_HIST`; a requested grid is named `_DIRECTION_GRID`. This option reads each input in both canonical and forward-only modes and cannot be reconstructed from a loaded BEDPE file.

`--grid <bool>`

Create a square grid containing every self comparison on the diagonal and every pairwise comparison off the diagonal. The grid is rendered as one Matplotlib figure and supports three or more input sequences, although large grids become visually dense.

`--grid-only <bool>`

Create the comparison grid without writing the individual dotplots.

`-t / --axes-ticks <list of ints>`

Custom tickmarks for x and y axis. Values outside of the `--axes-limits` will not be shown. 

`-a / --axes-limits <int>`

Change axis limits for x and y axis. Useful when comparing multiple plots, allowing them to stay in scale. 

`--bin-freq <bool>`

By default, histograms are evenly spaced based on the number of colors and the identity threshold. Select this argument to bin based on the frequency of observed identity values.

### Sample run - Static Plots

#### Using a config file

When running _ModDotPlot_ to produce static plots, it is recommended to use a config file. The config file is provided in JSON, and accepts the same syntax as the command line arguments shown above. Here is an sample run using a centromeric sequence of _Arabadopsis thaliana_:

```
$ cat config/config.json

{
    "identity": 90,
    "palette": "Blues_7",
    "breakpoints": [
        90,
        91,
        92,
        93,
        96,
        98,
        99,
        100
    ],
    "output_dir": "Arabadopsis",
    "fasta": [
        "sequences/Arabadopsis_chr1_centromere.fa"
    ]
}
```

```
$ moddotplot static -c config/config.json               
  __  __           _   _____        _     _____  _       _   
 |  \/  |         | | |  __ \      | |   |  __ \| |     | |  
 | \  / | ___   __| | | |  | | ___ | |_  | |__) | | ___ | |_ 
 | |\/| |/ _ \ / _` | | |  | |/ _ \| __| |  ___/| |/ _ \| __|
 | |  | | (_) | (_| | | |__| | (_) | |_  | |    | | (_) | |_ 
 |_|  |_|\___/ \__,_| |_____/ \___/ \__| |_|    |_|\___/ \__|

Running ModDotPlot in static mode

Retrieving k-mers from Chr1:14000001-18000000.... 

Progress: |████████████████████████████████████████| 100.0% Completed

Chr1:14000001-18000000 k-mers retrieved! 

Computing self identity matrix for Chr1:14000001-18000000... 

        Sequence length n: 4000000

        Window size w: 4000

        Modimizer sketch size: 1000

        Plot Resolution r: 1000

Progress: |████████████████████████████████████████| 100.0% Completed


Saved self-identity matrix as a paired-end bed file to Arabadopsis/Chr1:14000001-18000000/Chr1:14000001-18000000.bedpe

Triangle plots, full plots, and histogram for Arabadopsis/Chr1:14000001-18000000/Chr1:14000001-18000000 saved sucessfully.
```
![](images/Chr1:14000001-18000000_FULL.png)

Using `samtools faidx` will result in a genomic range being added to a fasta file's header (eg. in the above sequence, the header is Chr1:14000001-18000000). _ModDotPlot_ will parse this syntax to add the appropriate axis.

#### Adding custom bed file annotations

If providing a custom BED3-BED9 annotation file using `--bed/-b`, _ModDotPlot_ will output additional files:

- A collapsed annotation track `_ANNOTATION_TRACK` in PNG and the selected SVG, PDF, or PostScript vector format. Interval colors use the BED `itemRgb` value in column 9 when present, with a default color for BED3-BED8 records or invalid RGB values.
- The annotation track overlayed with a self-identity dotplot `_ANNOTATED` for each sequence present in the annotation track.

```
$ moddotplot static -f sequences/HG002_chr13_MATERNAL:1-4000000.fa -b config/hg002v1.1.cenSatv2.0.bed
  __  __           _   _____        _     _____  _       _   
 |  \/  |         | | |  __ \      | |   |  __ \| |     | |  
 | \  / | ___   __| | | |  | | ___ | |_  | |__) | | ___ | |_ 
 | |\/| |/ _ \ / _` | | |  | |/ _ \| __| |  ___/| |/ _ \| __|
 | |  | | (_) | (_| | | |__| | (_) | |_  | |    | | (_) | |_ 
 |_|  |_|\___/ \__,_| |_____/ \___/ \__| |_|    |_|\___/ \__|

Running ModDotPlot in static mode

...

Annotation track saved to chr13_MATERNAL:1-4000000/chr13_MATERNAL:1-4000000_ANNOTATION_TRACK

Triangle plots, full plots, and histogram for chr13_MATERNAL:1-4000000/chr13_MATERNAL:1-4000000 saved sucessfully. 

```
![](images/chr13_MATERNAL:1-4000000_TRI_ANNOTATED.png)


#### Comparing two sequences

![](images/moddotplot_comparative.png)

ModDotPlot can produce an a vs. b style dotplot for each pairwise combination of input sequences. Use the `--compare` command line argument to include these plots. When running `--compare` in interactive mode, a dropdown menu will appear, allowing the user to switch between self-identity and pairwise plots. Note that a maximum of two sequences are allowed in interactive mode. If you want to skip the creation of self-identity plots, you can use `--compare-only`:

```
moddotplot static -f sequences/*_MATERNAL*.fa --compare-only
```

![](images/chr13_MATERNAL:1-4000000_chr14_MATERNAL:1-4000000_COMPARE.png)

--- 

### Interactive Mode Commands

`-b / --bed <.bed file> [<.bed file> ...]`

Add one or more BED3-BED9 annotation files. A self-identity plot shows the
matching track beneath its x axis. A comparative plot shows independent x- and
y-axis tracks when BED chromosome names match both FASTA headers. A FASTA
header such as `chr14_MATERNAL:1-4000000` matches BED chromosome
`chr14_MATERNAL`, and the interactive axes retain those genomic coordinates.
For example:

```
moddotplot interactive -f sample1.fa sample2.fa --compare \
    --bed sample1.bed sample2.bed
```

`--port <int>`

Port to display ModDotPlot on. Default is 8050, this can be changed to any accepted port. 

`-w / --window <int>`

Minimum window size. By default, interactive mode sets a minimum window size based on the sequence length `n/2000` (eg. a 3Mbp sequence will have a 1500bp window). The maximum window size will always be set to `n/1000` (3000bp under the same example). This means that 2 matrices will be created.

`-q / --quick <bool>`

This will automatically run interactive mode with a minimum window size equal to the maximum window size (`n/1000`). This will result in a quick launch, however the resolution of the plot will not improve upon zooming in.

`-s / --save <bool>`

Save the matrices produced in interactive mode. By default, a folder called `interactive_matrices` will be saved in `--output_dir`, containing each matrix in compressed NumPy format, as well as metadata for each matrix in a pickle. Modifying the files in `interactive_matrices` will cause errors when attempting to load them in the future.

`--no-plot <bool>`

Save .bedpe to file, but skip rendering of plots. Must be used with `--save`.

`-l / --load <directory>`

Load previously saved matrices. Used instead of `-f/--fasta`.


--- 

### Sample run - Interactive Mode

```
$ moddotplot interactive -f sequences/Chr1_cen.fa     

  __  __           _   _____        _     _____  _       _   
 |  \/  |         | | |  __ \      | |   |  __ \| |     | |  
 | \  / | ___   __| | | |  | | ___ | |_  | |__) | | ___ | |_ 
 | |\/| |/ _ \ / _` | | |  | |/ _ \| __| |  ___/| |/ _ \| __|
 | |  | | (_) | (_| | | |__| | (_) | |_  | |    | | (_) | |_ 
 |_|  |_|\___/ \__,_| |_____/ \___/ \__| |_|    |_|\___/ \__|

Running ModDotPlot in interactive mode

Retrieving k-mers from Chr1:14000000-18000000.... 

Progress: |████████████████████████████████████████| 100.0% Completed

Chr1:14000000-18000000 k-mers retrieved! 

Building self-identity matrices for Chr1:14000000-18000000, using a minimum window size of 2000.... 

Layer 1 using window length 2000

Progress: |████████████████████████████████████████| 100.0% Completed


Layer 2 using window length 4000

Progress: |████████████████████████████████████████| 100.0% Completed


ModDotPlot interactive mode is successfully running on http://127.0.0.1:8050/ 

Dash is running on http://127.0.0.1:8050/
```

![](images/chr1_screenshot.png)

The plotly plot can be navigated using the zoom (magnifying glass) and pan (hand) icons. The plot can be reset by double-clicking or selecting the home button. The identity threshold can be modified by seelcting the slider. Colors can be readjusted according to the same gradient based on the new identity levels. 

### Sample run - Port Forwarding

Running interactive mode on an HPC environment can be accomplished through the use of port forwarding. On your remote server, run ModDotPlot as normal:

```
moddotplot interactive -f INPUT_FASTA_FILE(S) --port HPC_PORT_NUMBER
```

Then on your local machine, set up a port forwarding tunnel:

```
ssh -N -f -L <LOCAL_PORT_NUMBER>:127.0.0.1:<HPC_PORT_NUMBER> HPC@LOGIN.CREDENTIALS
```

You should now be able to view interactive mode using `http://127.0.0.1:<LOCAL_PORT_NUMBER>`. Note that your own HPC environment may have specific instructions and/or restrictions for setting up port forwarding.

VSCode now has automatic port forwarding built into the terminal menu. See [VSCode documentation](https://code.visualstudio.com/docs/editor/port-forwarding) for further details 

![](images/portforwarding.png)

--- 

## Questions

For bug reports or general usage questions, please raise a GitHub issue, or email alex ~dot~ sweeten ~at~ nih ~dot~ gov

--- 

## Known Issues

- Mac users might encounter the following unexpected command line output: `/bin/sh: lscpu: command not found`. This is a known issue with Plotnine, the Python plotting library used by ModDotPlot. This can be safely ignored.

- The error ` UserWarning: h5py is running against HDF5 1.xx.x when it was built against 1.xx.x, this may cause problems` is due to the h5py library used by cooler having conflicting versions in the dependency tree. This can also be safely ignored, but if you want to remove this message run `pip uninstall -y h5py` `pip install --no-binary=h5py h5py`
