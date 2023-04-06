import sys
from typing import List, Tuple
sys.path.append(".")
from const import COLS
import pandas as pd
import itertools
from moddotplot.estimate_identity import (
    binomial_distance,
    containment,
)

def make_scale(vals: float) -> str:
    return format(vals / 1e6, ",")


def make_k(vals: float) -> str:
    return format(vals / 1e3, ",")


def paired_bed_file(window_partitions: dict, input_name: str, id_threshold: int,
                     density: int, output: str) -> None:
    """
    Given a dictionary of window partitions, creates a BED file containing pairs of
    genomic regions with specified identity threshold.

    Args:
    window_partitions: a dictionary of window partitions containing genomic regions.
    input_name: a string representing the name of the input sequence.
    id_threshold: an integer specifying the identity threshold between 50 and 100.
    density: an integer specifying the density of the genomic regions.
    output: a string representing the name of the output BED file.

    Returns:
    None
    """
    assert id_threshold > 50 and id_threshold < 100

    cols: List[str] = COLS
    bed: List[List[str]] = []
    for w in itertools.combinations_with_replacement(window_partitions, 2):
        if input_name:
            query_name = input_name
        else:
            query_name = "input_sequence"
        query_start, query_end = w[0].split("-")
        reference_start, reference_end = w[1].split("-")
        perID = binomial_distance(
            containment(set(window_partitions[w[0]]), set(window_partitions[w[1]])), 21
        )
        if perID * 100 >= id_threshold:
            # TODO: Add strand orientation into bed file
            bed.append(
                [
                    query_name,
                    int(query_start) * density,
                    int(query_end) * density,
                    query_name,
                    int(reference_start) * density,
                    int(reference_end) * density,
                    perID * 100,
                ]
            )
    df = pd.DataFrame(bed, columns=cols)
    if not output:
        bedfile_output = input_name + ".bed"
        df.to_csv(bedfile_output, sep="\t")
    else:
        bedfile_output = output + ".bed"
        df.to_csv(bedfile_output, sep="\t")
