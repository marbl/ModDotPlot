![](images/logo.png)

- [About](#about)
- [Installation](#installation)
- [Usage](#usage)
  - [Interactive Mode](#interactive-mode)
  - [Static Mode](#static-mode)
  - [Standard arguments](#standard-arguments)
  - [Interactive Mode Commands](#interactive-mode-commands)
  - [Static Mode Commands](#static-mode-commands)
  - [Sample run - Interactive Mode](#sample-run---interactive-mode)
  - [Sample run - static plots](#sample-run---static-plots)
  - [Sample run - comparing two sequences](#sample-run---comparing-two-sequences)
- [Questions](#questions)
- [Known Issues](#known-issues)
- [Cite](#cite)

## About

ModDotPlot is a novel dot plot visualization tool, similar to [StainedGlass](https://mrvollger.github.io/StainedGlass/). ModDotPlot utilizes modimizers to compute the [Containment Index](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1841-x) between pairwise combinations of genomic intervals, and rapidly approximates their Average Nucleotide Identity. This significantly reduces the computational time required to produce these plots, enough to view multiple layers of resolution in real time!

![](images/demo.gif)

--- 

## Installation

```
git clone https://github.com/marbl/ModDotPlot.git
cd ModDotPlot
```

Although optional, it's recommended to setup a virtual environment before using ModDotPlot:

```
python -m venv venv
source venv/bin/activate
```

Once activated, you can install the required dependencies:

```
python -m pip install .
```

Finally, confirm that the installation was installed correctly by running `moddotplot -h`:
```       
  __  __           _   _____        _     _____  _       _   
 |  \/  |         | | |  __ \      | |   |  __ \| |     | |  
 | \  / | ___   __| | | |  | | ___ | |_  | |__) | | ___ | |_ 
 | |\/| |/ _ \ / _` | | |  | |/ _ \| __| |  ___/| |/ _ \| __|
 | |  | | (_) | (_| | | |__| | (_) | |_  | |    | | (_) | |_ 
 |_|  |_|\___/ \__,_| |_____/ \___/ \__| |_|    |_|\___/ \__|

v0.8.0 

usage: moddotplot [-h] {interactive,static} ...

ModDotPlot: Visualization of Complex Repeat Structures

positional arguments:
  {interactive,static}  Choose mode: interactive or static
    interactive         Interactive mode commands
    static              Static mode commands

options:
  -h, --help            show this help message and exit`
  ```
--- 

## Usage

ModDotPlot must be run either in `interactive` mode, or `static` mode:

### Interactive Mode

```
moddotplot interactive -f FASTA_FILE(S)
```

This will launch a [Dash application](https://plotly.com/dash/) on your machine's localhost. Open any web browser and go to `http://127.0.0.1:<PORT_NUMBER>` to view the interactive plot. Running `Ctrl+C` on the command line will exit the Dash application. The default port number used by Dash is `8050`, but this can be customized using the `--port` command (see [interactive mode commands](#interactive-mode-commands) for further info).

### Static Mode

```
moddotplot static -f FASTA_FILE(S) -o OUTPUT_DIRECTORY
```

This skips running Dash and quickly creates plots under the output directory using [plotnine](https://plotnine.readthedocs.io/en/v0.12.4/). By default, running ModDotPlot in static mode this will produce the following files:

- A paired-end bed file, containing intervals alongside their corresponding identity estimates
- A self-identity dotplot for each sequence
- A histogram of identity values for each sequence.
  
See [static mode commands](#static-mode-commands) for further info.

--- 

### Standard arguments

The following arguments are the same in both interactive and static mode:

`-k / --kmer <int>`

K-mer size to use. This should be large enough to distinguish unique k-mers with enough specificity, but not too large that sensitivity is removed. Default: 21.

`-o / --output-dir <string>`

Name of output directory for bed file & plots. Default is current working directory.

`-id / --identity <int>`

Minimum sequence identity cutoff threshold. Default is 86. While it is possible to go as low as 50% sequence identity, anything below 80% is not recommended. 

`--delta <float>`

Each partition takes into account a fraction of its neighboring partitions k-mers. This is to avoid sub-optimal identity scores when partitons don't overlap identically. Default is 0.5, and the accepted range is between 0 and 1. Anything greater than 0.5 is not recommended.

`-s / --sparsity <int>`

A higher sparsity value means less k-mers to compare, at the expense of lower accuracy. Modimizers are selected `0 mod s`, an inverse of the selected k-mer density. The default is `s = 2/Mbp` of sequence for fast performance without compromising accuracy. For example, on a human Y chromosome ~ 62Mbp, ModDotPlot will set `s = 124`. Interactive mode will automatically round up to the nearest power of 2, in order to make computations easier. 


`-r / --resolution <int>`

Dotplot resolution. This corresponds to the number of windows each input sequence is partitioned into. Default is 1000. Overrides the `--window` parameter.

`--compare <bool>`

If set when 2 or more sequences are input into ModDotPlot, this will show an a vs. b style plot, in addition to a self-identity plot. Note that interactive mode currently only supports a maximum of two sequences. If more than two sequences are input, only the first two will be shown.

`--compare-only <bool>`

If set when 2 or more sequences are input into ModDotPlot, this will show an a vs. b style plot, without showing self-identity plots.

`--ambiguous <bool>`
By default, k-mers that are homopolymers of ambiguous IUPAC codes (eg. NNNNNNNNNNN’s) are excluded from identity estimation. This results in gaps along the central diagonal for these regions.  If desired, these can be kept by setting the `—-ambiguous` flag in both interactive and static mode.  

--- 

### Interactive Mode Commands

`--port <int>`

Port to display ModDotPlot on. Default is 8050, this can be changed to any accepted port. 

Running interactive mode on an HPC environment can be accomplished through the use of port forwarding. On your remote server, run ModDotPlot as normal:

```
moddotplot -i INPUT_FASTA_FILE(S) --port HPC_PORT_NUMBER
```

Then on your local machine, set up a port forwarding tunnel:

```
ssh -N -f -L LOCAL_PORT_NUMBER:127.0.0.1:HPC_PORT_NUMBER HPC@LOGIN.CREDENTIALS
```

You should now be able to view interactive mode using `http://127.0.0.1:<LOCAL_PORT_NUMBER>`. Note that your own HPC environment may have specific instructions and/or restrictions for setting up port forwarding.

`-w / --window <int>`

Minimum window size. By default, interactive mode sets a minimum window size based on the sequence length `n/2000` (eg. a 3Mbp sequence will have a 1500bp window). The maximum window size will always be set to `w/1000` (3000bp under the same example). This means that 2 matrices will be created.

`-s / --save <bool>`

Save the matrices produced in interactive mode. By default, a folder called `interactive_matrices` will be saved in `--output_dir`, containing each matrix in compressed NumPy format, as well as metadata for each matrix in a pickle. Modifying the files in `interactive_matrices` will cause errors when attempting to load them in the future.

`-l / --load <directory>`

Load previously saved matrices. Used instead of `-f/--fasta`

--- 

### Static Mode Commands

`-c / --config`
Run `moddotplot static` with a config file, rather than (sample syntax). Recommended when creating a really customized plot. Used instead of `-f/--fasta`.

`-b / --bed`
Create a plot from a previously computed pairwise bed file. Skips Average Nucleotide Identity computation. Used instead of `-f/--fasta`. 

`-w / --window <int>`

Window size. Unlike interactive mode, only one matrix will be created, so this represents the *only* window size. Default is set to `n/1000` (eg. 3000bp for a 3Mbp sequence). 

`--no-bed <bool>`

Skip output of bed file.

`--no-plot <bool>`

Skip output of pdf and png image files.

`--no-hist <bool>`

Skip output of histogram legend.

`--width <int>`

Adjust width of self dot plots. Default is 9 inches.

`--dpi`

Image resolution in dots per inch (not to be confused with dotplot resolution). Default is 600.

`--palette`

List of accepted palettes can be found [here](https://jiffyclub.github.io/palettable/colorbrewer/). The syntax is to have the name of the palette, followed by an underscore with the number of colors, eg. `OrRd_8`. Default is `Spectral_11`.

`--palette-orientation <bool>`

Flip sequential order of color palette. Set to `-` by default for divergent palettes. 

`--breakpoints <list>`

Add custom identity threshold breakpoints. Note that the number of breakpoints must be equal to the number of colors + 1.

`--bin-freq <bool>`

By default, histograms are evenly spaced based on the number of colors and the identity threshold. Select this argument to bin based on the frequency of observed identity values.

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

v0.8.0 

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

### Sample run - static plots

When running ModDotPlot to produce static plots, it is recommended to use a config file, especially when creating extremely customized plots. The config file is provided in JSON, and accepts the same syntax as the command line arguments shown above. 

```
$ cat config/config.json

{
    "identity": 90,
    "sparsity": 10,
    "palette": "OrRd_7",
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
    "output_dir": "Chr1_cen_plots",
    "fasta": [
        "../sequences/Chr1_cen.fa"
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

v0.8.0 

Running ModDotPlot in static mode

Retrieving k-mers from Chr1:14M-18M.... 

Progress: |████████████████████████████████████████| 100.0% Completed

Chr1:14M-18M k-mers retrieved! 

Computing self identity matrix for Chr1:14M-18M... 

        Sparsity value s: 10

        Sequence length n: 4000001

        Window size w: 4000

Progress: |████████████████████████████████████████| 100.0% Completed


Saved bed file to Chr1_cen_plots/Chr1:14M-18M.bed

Plots created! Saving to Chr1_cen_plots/Chr1:14M-18M...

Chr1_cen_plots/Chr1:14M-18M_TRI.png, Chr1_cen_plots/Chr1:14M-18M_TRI.pdf, Chr1_cen_plots/Chr1:14M-18M_FULL.png, Chr1_cen_plots/Chr1:14M-18M_FULL.png, Chr1_cen_plots/Chr1:14M-18M_HIST.png and Chr1_cen_plots/Chr1:14M-18M_HIST.pdf, saved sucessfully. 
```
![](images/Chr1:14M-18M_FULL.png)

--- 


### Sample run - comparing two sequences

ModDotPlot can produce an a vs. b style dotplot for each pairwise combination of input sequences. Use the `--compare` command line argument to include these plots. When running `--compare` in interactive mode, a dropdown menu will appear, allowing the user to switch between self-identity and pairwise plots. Note that a maximum of two sequences are allowed in interactive mode. If you want to skip the creation of self-identity plots, you can use `--compare-only`:

```
moddotplot interactive -f sequences/chr14_segment.fa sequences/chr21_segment.fa --compare-only
```

--- 

## Questions

For bug reports or general usage questions, please raise a GitHub issue, or email alex ~dot~ sweeten ~at~ nih ~dot~ gov

--- 

## Known Issues

Mac users might encounter the following unexpected command line output:

`/bin/sh: lscpu: command not found`

This is a known issue with Plotnine, the Python plotting library used by ModDotPlot. This can be safely ignored.


---


## Cite

Publication in progress!
