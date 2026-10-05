![](images/logo.png)

---
[![PyPI](https://img.shields.io/pypi/v/ModDotPlot?color=blue&label=PyPI)](https://pypi.org/project/ModDotPlot/)
[![CI](https://github.com/marbl/ModDotPlot/actions/workflows/ci.yml/badge.svg)](https://github.com/marbl/ModDotPlot/actions/workflows/ci.yml)

- [Cite](#cite)
- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
  - [Command Line Arguments](#command-line-arguments)
    - [General Options](#general-options)
    - [Input Options](#input-options)
    - [Analysis Options](#analysis-options)
    - [Output Options](#output-options)
    - [Plot Formatting Options](#plot-formatting-options)
    - [Plot Customization Options](#plot-customization-options)
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

![](images/demo.gif)

If you're interested in learning more about _ModDotPlot_ and how to visualize tandem repeats, we have an in-depth [YouTube video tutorial](https://www.youtube.com/watch?v=_7sQaljB_ys&t=2321s&pp=ygUXYWxleCBzd2VldGVuIG1vZGRvdHBsb3Q%3D) hosted by the [BioDiversity Genomics Academy](https://thebgacademy.org).

--- 

## Installation

_ModDotPlot_ can be installed by running `pip install moddotplot`. Alternatively, you can download the current release from GitHub by using:

```
git clone https://github.com/marbl/ModDotPlot.git
cd ModDotPlot
```

Although optional, it's recommended to set up a virtual environment before using _ModDotPlot_:

```
python -m venv venv
source venv/bin/activate
```

Once activated, install either the base package for static plotting:

```
python -m pip install .
```

If desired, add the deprecated interactive plotting dependencies, Cooler export support, or both only when needed:

```
python -m pip install ".[interactive]"
python -m pip install ".[cooler]"
python -m pip install ".[interactive,cooler]"
```

Finally, confirm that the installation completed correctly and that your version is up to date by running `moddotplot -h`:
```       
  __  __           _   _____        _     _____  _       _   
 |  \/  |         | | |  __ \      | |   |  __ \| |     | |  
 | \  / | ___   __| | | |  | | ___ | |_  | |__) | | ___ | |_ 
 | |\/| |/ _ \ / _` | | |  | |/ _ \| __| |  ___/| |/ _ \| __|
 | |  | | (_) | (_| | | |__| | (_) | |_  | |    | | (_) | |_ 
 |_|  |_|\___/ \__,_| |_____/ \___/ \__| |_|    |_|\___/ \__|

 v1.0.0

usage: moddotplot [-h] [--quiet] [{static,interactive}] ...

ModDotPlot: Visualization of Tandem Repeats

positional arguments:
  {static,interactive}  Choose mode; static is used when omitted
    static              Static mode commands (default)
    interactive         Interactive mode commands (deprecated; explicit use only)

options:
  -h, --help            show this help message and exit
  --quiet               suppress all console output, including warnings and errors
  ```

Note that running `moddotplot -h` might take a while at first! This is because the Python interpreter is compiling source code into the __pycache__ directory. Subsequent runs will use the pre-compiled code and load much faster!

--- 

## Usage

_ModDotPlot_ runs in static mode by default. The explicit `static` subcommand
is retained for clarity and compatibility, so these commands are equivalent:

```bash
moddotplot -f sequence.fa <ARGS>
moddotplot static -f sequence.fa <ARGS>
```

A standard run writes a BEDPE identity table, full and triangle dotplots,
and an identity histogram beneath the selected output directory. Plots and
histograms are rendered with Matplotlib as both PNG and the selected vector
format (SVG by default). Every plot directory also receives a
`plot_summary.txt` reproducibility record containing input paths, plotting
parameters, and the command used for the run.

![](images/moddotplot_output.png)

The deprecated interactive application remains available only through the
explicit `interactive` subcommand. Its commands are documented under
[Interactive Mode Commands](#interactive-mode-commands).

### Command Line Arguments

The default/static command accepts every option in the following two sections.
Options marked **shared** are also accepted by the deprecated interactive
command; interactive-specific behavior is summarized separately below.
Flags such as `--grid` and `--no-plot` are switches and do not take a
`true` or `false` value.

#### General Options

| Argument | Description |
| --- | --- |
| `-h, --help` | Show help for the root command or selected subcommand and exit. |
| `--quiet` | Suppress console output, including warnings and errors. The exit status still reports success or failure. May appear before or after the subcommand. |

#### Input Options

| Argument | Description |
| --- | --- |
| `-f, --fasta FILE [FILE ...]` | Read one or more FASTA, gzip-compressed FASTA, or BGZF FASTA files. Static mode analyzes every record unless `--sequence` limits the selection. Mutually exclusive with static `--load`. |
| `-l, --load BEDPE [BEDPE ...]` | Plot one or more ModDotPlot BEDPE files without recomputing identity. Mutually exclusive with `--fasta`, but may be combined with `--config`; explicit input paths override input paths in the config. Interactive `--load` has different behavior described below. |
| `-c, --config JSON` | Load command settings from a JSON config. Explicit FASTA, BEDPE, config, and output paths are resolved independently; explicit input/output paths take precedence, while config values supply other settings. |
| `-s, --sequence ID [ID ...]` | Analyze only named records from FASTA input. IDs match the first whitespace-delimited FASTA token, preferring exact matches and then unambiguous case-insensitive matches. |
| `--region ID:START-END [ID:START-END ...]` | Analyze 1-based inclusive regions. Region limits propagate to identity computation, BEDPE coordinates, plots, and grid axes. |
| `--pairs FILE` | Restrict `--compare` or `--compare-only` to pairs in a two-column text file. Blank lines and `#` comments are ignored. Requires indexed FASTA input; duplicate, self, unknown, and ambiguous pairs are errors. |
| `-b, --bed BED` | Add one BED3–BED9 annotation file. BED chromosome names must match FASTA identifiers after any `:start-end` suffix is removed. Valid `itemRgb` values are retained. Interactive mode accepts multiple BED files. |

#### Analysis Options

| Argument | Description |
| --- | --- |
| `--compare` | Add pairwise comparisons while retaining self-comparisons. Static mode considers all selected pairs unless `--pairs` restricts them. |
| `--compare-only` | Produce pairwise comparisons without self-comparisons. Mutually exclusive with `--compare`. |
| `-k, --kmer INT` | Set k-mer length. Default: `21`. |
| `-m, --modimizer INT` | Set the modimizer sketch target. Must be smaller than the window length. Smaller values are faster but reduce accuracy. Default: `1000`. |
| `-r, --resolution INT` | Set the approximate number of sequence windows. Mutually exclusive with `--window` in static mode. Default: `1000`. |
| `-w, --window INT` | Set window length in base pairs. When omitted, it is inferred from sequence length and resolution. Mutually exclusive with `--resolution` in static mode. |
| `-id, --identity FLOAT` | Set the minimum estimated identity percentage. Default: `86.0`. Values below 80 are generally not recommended. |
| `-d, --delta FLOAT` | Include this fraction of each neighboring window when estimating identity. Accepted range: 0–1. Default: `0.5`. |
| `--forward` | Hash forward k-mers only instead of canonical k-mers, producing strand-specific output. |
| `--ambiguous` | Include deterministic hashes for windows containing non-ACGTU IUPAC bases instead of masking those windows. |
| `--processes N` | Use 1–4 independent chromosome or comparison-group workers. When omitted, indexed multi-record input uses a bounded automatic worker count; unindexed and ordinary gzip input remain sequential. |
| `--memory-limit GIB` | Set an aggregate GiB memory budget for comparison workers. When omitted, available memory is used when the platform exposes it. |
| `--sketch-cache DIRECTORY` | Persist compact prepared sketches for reuse in later indexed comparative runs. Entries are keyed by input identity, record/region, strand mode, and sketch parameters. |

#### Output Options

| Argument | Description |
| --- | --- |
| `-o, --output-dir DIRECTORY` | Set the output directory. Static mode writes BEDPE, plots, and summaries; interactive mode writes saved matrices and coordinate logs. Default: the current directory. |
| `--cooler` | Write Cooler matrices in addition to BEDPE. Install with `pip install "ModDotPlot[cooler]"` or `pip install ".[cooler]"`. |
| `--no-bedpe` | Skip BEDPE output. |
| `--no-plot` | Skip all plot rendering in static mode. In interactive mode, prevent Dash from launching; must be combined with `--save`. |
| `--no-hist` | Skip identity histograms. |

#### Plot Formatting Options

| Argument | Description |
| --- | --- |
| `--grid` | Render selected self- and pairwise comparisons in a single square grid, in addition to individual plots. |
| `--grid-only` | Render only the comparison grid and skip individual plots. |
| `--compare-order {sequential,size}` | Choose comparative axis order. `sequential` preserves input order; `size` places the larger sequence on the x-axis. Default: `sequential`. |
| `-a, --axes-limits FLOAT` | Set common x/y axis limits for self-identity plots. The value cannot be shorter than the sequence. |
| `-t, --axes-ticks INT [INT ...]` | Set explicit x/y tick positions. Ticks outside the visible limits are omitted. |
| `--axes-number VALUE` | Retained for configuration compatibility as the requested number of axis ticks; currently unused by the Matplotlib renderer. Default: `7`. |
| `--width FLOAT` | Set plot width in inches. For grids, this is the total grid width, not the width of each cell. Default: `9`. |
| `--dpi INT` | Set raster resolution in dots per inch. Default: `300`. |
| `--vector {svg,pdf,ps}` | Select the vector output format. Default: `svg`. |
| `--deraster` | Keep dotplot tiles as vector geometry instead of rasterizing them inside vector output. This can produce very large files. |

#### Plot Customization Options

| Argument | Description |
| --- | --- |
| `--palette NAME_COUNT` | Select an exact discrete [ColorBrewer](https://colorbrewer2.org/) palette, such as `OrRd_8`. Default: `Spectral_11`. |
| `--palette-orientation {+,-}` | Select forward or reversed palette order. Diverging palettes retain ModDotPlot's historical orientation convention. Default: `+`. |
| `--colors COLOR [COLOR ...]`, `--color ...` | Supply a custom low-to-high color sequence in hexadecimal or RGB form. `--color` is a legacy alias. |
| `--breakpoints VALUE [VALUE ...]` | Supply custom identity thresholds between the identity cutoff and 100. The number of breakpoints must equal the number of colors plus one. |
| `--bin-freq` | Derive identity color bins from the observed value distribution instead of evenly spacing them between the identity cutoff and 100. |
| `--plot-direction` | Compute strand direction and color matches blue for the same orientation and pink for reverse orientation, with ANI represented by shade intensity. Available only with FASTA input. |
---

### Sample run - Static Plots

#### Using a config file

When running _ModDotPlot_ to produce static plots, it is recommended to use a config file. The config file is provided in JSON, and accepts the same syntax as the command line arguments shown above. Here is a sample run using a centromeric sequence of _Arabidopsis thaliana_:

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
    "output_dir": "Arabidopsis",
    "fasta": [
        "sequences/Arabidopsis_chr1_centromere.fa"
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

Chr1:14000001-18000000 k-mers retrieved! 

Computing self identity matrix for Chr1:14000001-18000000... 

        Sequence length n: 4000000

        Window size w: 4000

        Modimizer sketch size: 1000

        Plot Resolution r: 1000

Saved self-identity matrix as a paired-end bed file to Arabidopsis/Chr1:14000001-18000000/Chr1:14000001-18000000.bedpe

Triangle plots, full plots, and histogram for Arabidopsis/Chr1:14000001-18000000/Chr1:14000001-18000000 saved successfully.
```
![](images/Chr1:14000001-18000000_FULL.png)

Using `samtools faidx` will result in a genomic range being added to a FASTA file's header (e.g., in the above sequence, the header is Chr1:14000001-18000000). _ModDotPlot_ will parse this syntax to add the appropriate axis.

#### Adding custom bed file annotations

If providing a custom BED3-BED9 annotation file using `--bed/-b`, _ModDotPlot_ will output additional files:

- A collapsed annotation track `_ANNOTATION_TRACK` in PNG and the selected SVG, PDF, or PostScript vector format. Interval colors use the BED `itemRgb` value in column 9 when present, with a default color for BED3-BED8 records or invalid RGB values.
- The annotation track overlaid with a self-identity dotplot `_ANNOTATED` for each sequence present in the annotation track.

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

Triangle plots, full plots, and histogram for chr13_MATERNAL:1-4000000/chr13_MATERNAL:1-4000000 saved successfully.

```
![](images/chr13_MATERNAL:1-4000000_TRI_ANNOTATED.png)


#### Comparing two sequences

![](images/moddotplot_comparative.png)

ModDotPlot can produce an a vs. b style dotplot for each pairwise combination of input sequences. Use the `--compare` command line argument to include these plots. When running `--compare` in interactive mode, a dropdown menu will appear, allowing the user to switch between self-identity and pairwise plots. Note that a maximum of two sequences are allowed in interactive mode. If you want to skip the creation of self-identity plots, you can use `--compare-only`:

```
moddotplot static -f sequences/*_MATERNAL*.fa --compare-only
```

For a diploid multi-record assembly, a pair manifest avoids comparing every
chromosome and unplaced contig against every other record:

```text
# homologs.tsv
chr1_mat_hsa1 chr1_pat_hsa1
chr2_mat_hsa3 chr2_pat_hsa3
```

```bash
moddotplot static -f diploid.fa.gz --compare-only \
  --pairs homologs.tsv --processes 4 --memory-limit 32 \
  --sketch-cache .moddotplot-sketches
```

For BGZF-compressed FASTA, both `.fai` and `.gzi` indexes must accompany the
input. Indexed comparative mode fetches one record at a time, sketches it
directly, and releases pair-local data before continuing; it does not retain a
`uint64` positional hash for every base in the genome.

![](images/chr13_MATERNAL:1-4000000_chr14_MATERNAL:1-4000000_COMPARE.png)

--- 

### Interactive Mode Commands

**Deprecated:** interactive mode is maintenance-only and will not receive new
features. New browser-based work should use
[ModDotPlot Browser](https://marbl.github.io/ModDotPlot-Browser/). The legacy
Dash application remains available through the explicit subcommand:

```bash
moddotplot interactive <ARGS>
```

Install its optional dependencies with
`pip install "ModDotPlot[interactive]"`, or
`pip install ".[interactive]"` from a source checkout. The application
listens on `http://127.0.0.1:8050` by default and exits when the process is
stopped with `Ctrl+C`.

Interactive mode accepts the following arguments. Some share names with the
default/static command but have interactive-specific behavior.

| Argument | Description |
| --- | --- |
| `--quiet` | Suppress all console output, including warnings and errors. It may appear before or after `interactive`. |
| `-f, --fasta FILE [FILE ...]` | Read FASTA input and compute an interactive matrix hierarchy. Mutually exclusive with `--load`; interactive displays support at most two sequences. |
| `-l, --load DIRECTORY` | Load a previously saved `interactive_matrices` directory containing compressed matrices and `metadata.pkl`. Mutually exclusive with `--fasta`. This is not the static BEDPE loader. |
| `-b, --bed BED [BED ...]` | Add one or more BED3-BED9 annotation files. Tracks are aligned to matching FASTA identifiers on each matrix axis. |
| `-o, --output-dir DIRECTORY` | Set the directory used for saved matrices and coordinate logs. Default: current directory. |
| `-k, --kmer INT` | Set k-mer length. Default: `21`. |
| `-m, --modimizer INT` | Set the modimizer sketch target. Default: `1000`. |
| `-r, --resolution INT` | Set interactive dotplot resolution. Default: `1000`. |
| `-w, --window INT` | Set the minimum interactive window length. When omitted it is inferred from sequence length and resolution. |
| `-id, --identity FLOAT` | Set the minimum estimated identity percentage. Default: `86.0`. |
| `-d, --delta FLOAT` | Include this fraction of neighboring windows during identity estimation. Default: `0.5`. |
| `--compare` | Add a pairwise comparison while retaining self-comparisons. |
| `--compare-only` | Produce the pairwise comparison without self-comparisons. Mutually exclusive with `--compare`. |
| `--ambiguous` | Include deterministic hashes for windows containing non-ACGTU IUPAC bases. |
| `--forward` | Hash forward k-mers only instead of canonical k-mers. |
| `-s, --save` | Save the matrix hierarchy under `OUTPUT_DIR/interactive_matrices` as compressed NumPy arrays plus `metadata.pkl`. |
| `--port INT` | Set the localhost port used by Dash. Default: `8050`. |
| `-q, --quick` | Build a single matrix layer instead of the normal hierarchy for a faster launch without progressively finer zoom resolution. |
| `--no-plot` | Save matrices without launching Dash. Must be combined with `--save`. |

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

Chr1:14000000-18000000 k-mers retrieved! 

Building self-identity matrices for Chr1:14000000-18000000, using a minimum window size of 2000.... 

Layer 1 using window length 2000

Layer 2 using window length 4000

ModDotPlot interactive mode is successfully running on http://127.0.0.1:8050/ 

Dash is running on http://127.0.0.1:8050/
```

![](images/chr1_screenshot.png)

The Plotly plot can be navigated using the zoom (magnifying glass) and pan (hand) icons. The plot can be reset by double-clicking or selecting the home button. The identity threshold can be modified by selecting the slider. Colors can be readjusted according to the same gradient based on the new identity levels.

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

VS Code now has automatic port forwarding built into the terminal menu. See [VS Code documentation](https://code.visualstudio.com/docs/editor/port-forwarding) for further details.

![](images/portforwarding.png)

--- 

## Questions

For bug reports or general usage questions, please raise a GitHub issue, or email alex ~dot~ sweeten ~at~ nih ~dot~ gov

--- 

## Known Issues

- When the optional `ModDotPlot[cooler]` dependencies are installed, Cooler may report `UserWarning: h5py is running against HDF5 1.xx.x when it was built against 1.xx.x, this may cause problems`. This can be safely ignored. To remove the warning, reinstall h5py against the local HDF5 library with `pip uninstall -y h5py` followed by `pip install --no-binary=h5py h5py`.
