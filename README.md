## Mod.Plot

Mod.Plot is still early in development. We welcome your feedback and look forward to introducing even faster static plots, as well as a Graphical User Interface!

Mod.Plot is a novel dot plot visualization tool used to view tandem repeats, similar to StainedGlass. Our approaches uses sketching and modimizers to significantly reduce the computational time to produce these plots.

View our [Biological Data Science poster](https://docs.google.com/presentation/d/1SR833K-a2alIDtXVuyav_y33SVZ-0N3x/edit?usp=sharing&ouid=116747761671966787462&rtpof=true&sd=true)!

![](images/chrY_levels.png)

## Installation

`git clone https://github.com/marbl/ModDotPlot.git`

`python setup.py install`

## Usage

`cd moddotplot`

`python mod_identity.py -i INPUT [-k KMER_SIZE (default: 21)] [-d DENSITY (default: genome size / 1000000)] [-p PREFIX] [-o OUTPUT_FOLDER]  [-nc NON_CANONICAL] [--index INDEX_FILE] [-w WINDOW_SIZE]`

This produces a paired end bedfile. 

`RScript heatmap.R -b BEDFILE -p PREFIX`

Example run:

```
$ python moddotplot/identity.py -i sequences/CP086569.fa -p CP086569 -d 4

Parsing k-mers...
K-mers parsed!
Modimizers done!
Creating coordinates...
Coordinates done!
Computing identity...
Bed file created!

RScript moddotplot/heatmap.R -b CP086569.bed -p CP086569
```

For good performance, we recommend using a modimizer density d <= 1Mbp of sequence. We also recommend rounding to a scalable number for heirarchical adjustments. For example, plotting a human Y chromosome should use d = 64.
