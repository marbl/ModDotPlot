#!/usr/bin/env python3
"""Benchmark prepared-sketch reuse for a two-sequence static grid.

The benchmark isolates the work before matrix comparison.  With a fixed window,
the uncached workflow prepares both sequences for their self matrices and then
prepares both again for the pairwise matrix (four preparations).  The cache
performs the same access pattern with two preparations and two exact hits.

Example::

    python benchmarks/benchmark_sketch_reuse.py --length 1000000 --repeats 5
"""

from __future__ import annotations

import argparse
import gc
import math
import statistics
import sys
import time
from pathlib import Path
from typing import Callable, Dict, List, Sequence

import numpy as np

try:
    from moddotplot.estimate_identity import (
        ModimizerSketchCache,
        prepare_modimizer_sketches,
    )
except ModuleNotFoundError:  # Permit running from an uninstalled source tree.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from moddotplot.estimate_identity import (
        ModimizerSketchCache,
        prepare_modimizer_sketches,
    )


def _time(function: Callable[[], None]) -> float:
    gc.collect()
    started = time.perf_counter()
    function()
    return time.perf_counter() - started


def _configuration(length: int, resolution: int, modimizer: int):
    window = math.ceil(length / resolution)
    effective_modimizer = min(window, modimizer)
    raw_sparsity = round(window / effective_modimizer)
    if raw_sparsity <= effective_modimizer:
        sparsity = 2 ** int(math.log2(raw_sparsity))
    else:
        sparsity = 2 ** (int(math.log2(raw_sparsity - 1)) + 1)
    return window, sparsity, round(window / sparsity)


def benchmark(
    length: int,
    resolution: int,
    modimizer: int,
    delta: float,
    kmer: int,
    repeats: int,
) -> None:
    rng = np.random.default_rng(56)
    sequences = [
        rng.integers(0, np.iinfo(np.uint64).max, length, dtype=np.uint64)
        for _ in range(2)
    ]
    window, sparsity, expectation = _configuration(length, resolution, modimizer)

    def prepare(sequence):
        return prepare_modimizer_sketches(
            length,
            sequence,
            window,
            sparsity,
            delta,
            kmer,
            False,
            expectation,
        )

    def uncached():
        for index in (0, 1, 1, 0):
            prepare(sequences[index])

    def cached():
        cache = ModimizerSketchCache(max_entries=2)
        for index in (0, 1, 1, 0):
            cache.get_or_prepare(
                index,
                length,
                sequences[index],
                window,
                sparsity,
                delta,
                kmer,
                False,
                expectation,
            )
        cache.clear()

    # Warm NumPy dispatch and allocators before recording samples.
    prepare(sequences[0][: min(length, window)])
    samples: Dict[str, List[float]] = {"uncached": [], "cached": []}
    for repeat in range(repeats):
        order: Sequence[str] = (
            ("uncached", "cached") if repeat % 2 == 0 else ("cached", "uncached")
        )
        for name in order:
            samples[name].append(_time(uncached if name == "uncached" else cached))

    uncached_median = statistics.median(samples["uncached"])
    cached_median = statistics.median(samples["cached"])
    print(
        f"{length:,} hashes/sequence; window={window:,}; delta={delta}; "
        f"resolution={resolution}; repeats={repeats}"
    )
    print(f"uncached (4 preparations): {uncached_median:.6f} s")
    print(f"cached   (2 preparations): {cached_median:.6f} s")
    print(f"preparation speedup:       {uncached_median / cached_median:.2f}x")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, default=100_000)
    parser.add_argument("--resolution", type=int, default=100)
    parser.add_argument("--modimizer", type=int, default=100)
    parser.add_argument("--delta", type=float, default=0.5)
    parser.add_argument("--kmer", type=int, default=21)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args(argv)
    if min(args.length, args.resolution, args.modimizer, args.kmer, args.repeats) <= 0:
        parser.error(
            "length, resolution, modimizer, kmer, and repeats must be positive"
        )
    if args.delta < 0:
        parser.error("delta must be non-negative")
    benchmark(
        args.length,
        args.resolution,
        args.modimizer,
        args.delta,
        args.kmer,
        args.repeats,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
