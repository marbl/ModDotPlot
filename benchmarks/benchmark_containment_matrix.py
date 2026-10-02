#!/usr/bin/env python3
"""Benchmark the exact sparse containment-matrix implementation.

The generated sketches model the default ``delta=0.5`` layout: each expanded
window contains its core plus half of each adjacent core.  This exercises the
same dimensions and sketch cardinalities as a roughly 100 Mb sequence at the
default 1,000-bin resolution without allocating a genome-sized hash array.

Example::

    python benchmarks/benchmark_containment_matrix.py
    python benchmarks/benchmark_containment_matrix.py --resolution 2000
"""

from __future__ import annotations

import argparse
import resource
import sys
import time
from pathlib import Path

import numpy as np

try:
    from moddotplot.estimate_identity import (
        pairwiseContainmentMatrix,
        selfContainmentMatrix,
    )
except ModuleNotFoundError:  # Permit running from an uninstalled source tree.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from moddotplot.estimate_identity import (
        pairwiseContainmentMatrix,
        selfContainmentMatrix,
    )


def build_sketches(resolution: int, sketch_size: int, seed: int):
    rng = np.random.default_rng(seed)
    core = [
        np.unique(
            rng.integers(0, np.iinfo(np.uint64).max, sketch_size, dtype=np.uint64)
        )
        for _ in range(resolution)
    ]
    expanded = []
    halfway = sketch_size // 2
    for index, sketch in enumerate(core):
        pieces = [sketch]
        if index:
            pieces.append(core[index - 1][halfway:])
        if index + 1 < resolution:
            pieces.append(core[index + 1][:halfway])
        expanded.append(np.unique(np.concatenate(pieces)))
    return core, expanded


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, default=1000)
    parser.add_argument("--sketch-size", type=int, default=1612)
    parser.add_argument("--identity", type=float, default=86)
    parser.add_argument("--kmer", type=int, default=21)
    parser.add_argument("--seed", type=int, default=519)
    args = parser.parse_args(argv)

    core, expanded = build_sketches(args.resolution, args.sketch_size, args.seed)
    compact_bytes = sum(sketch.nbytes for sketch in core + expanded)

    started = time.perf_counter()
    self_matrix = selfContainmentMatrix(
        core, expanded, args.kmer, args.identity, ambiguous=False
    )
    self_seconds = time.perf_counter() - started

    started = time.perf_counter()
    pair_matrix = pairwiseContainmentMatrix(
        core,
        core,
        expanded,
        expanded,
        args.identity,
        args.kmer,
    )
    pair_seconds = time.perf_counter() - started

    print(f"resolution: {args.resolution}")
    print(f"core hashes/window: {args.sketch_size}")
    print(f"expanded hashes/window: about {args.sketch_size * 2}")
    print(f"prepared sketch memory: {compact_bytes / 2**20:.2f} MiB")
    print(f"self matrix: {self_seconds:.3f} s ({self_matrix.shape})")
    print(f"pair matrix: {pair_seconds:.3f} s ({pair_matrix.shape})")
    print("two-self-plus-pair estimate: " f"{2 * self_seconds + pair_seconds:.3f} s")
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux and other supported Unix platforms report KiB.
    peak_rss_mib = peak_rss / (2**20 if sys.platform == "darwin" else 2**10)
    print(f"peak process RSS: {peak_rss_mib:.2f} MiB")


if __name__ == "__main__":
    main()
