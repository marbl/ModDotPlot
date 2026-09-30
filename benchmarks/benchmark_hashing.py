#!/usr/bin/env python3
"""Compare ModDotPlot's ntHash2 path with the removed mmh3 implementation.

``mmh3`` is deliberately optional.  When it is installed, this utility
recreates the legacy ModDotPlot loop for an apples-to-apples migration
benchmark.  Otherwise it reports ntHash2 throughput on its own.

Examples::

    python benchmarks/benchmark_hashing.py --length 1000000 --repeats 7
    python benchmarks/benchmark_hashing.py --fasta sequence.fa --kmer 21
    python benchmarks/benchmark_hashing.py --json results.json
"""

from __future__ import annotations

import argparse
import gc
import gzip
import json
import random
import statistics
import sys
import time
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from moddotplot.parse_fasta import _hash_sequence
except ModuleNotFoundError:  # Permit running from an uninstalled source tree.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from moddotplot.parse_fasta import _hash_sequence


DNA_ALPHABET = "ACGT"
REVERSE_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _load_mmh3():
    """Return the optional legacy module without making it a dependency."""
    try:
        import mmh3  # type: ignore[import-not-found]
    except ImportError:
        return None
    return mmh3


def _legacy_mmh3_hashes(
    sequence: str, k: int, canonical: bool, mmh3_module
) -> List[int]:
    """Reproduce ModDotPlot's pre-ntHash2 per-k-mer implementation."""
    result = []
    for start in range(max(len(sequence) - k + 1, 0)):
        kmer = sequence[start : start + k].upper()
        forward = mmh3_module.hash(kmer)
        if canonical:
            reverse = mmh3_module.hash(kmer[::-1].translate(REVERSE_COMPLEMENT))
            result.append(min(forward, reverse))
        else:
            result.append(forward)
    return result


def _moddotplot_nthash2_hashes(sequence: str, k: int, canonical: bool):
    """Exercise the compact batch hashing path used by ModDotPlot's CLI."""
    return _hash_sequence(
        sequence,
        k,
        fw_only=not canonical,
        ambiguous=False,
    )


def _read_first_fasta(path: Path) -> Tuple[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    name: Optional[str] = None
    chunks: List[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    break
                name = line[1:].split()[0] or path.name
            elif name is not None:
                chunks.append(line)
    if name is None:
        raise ValueError(f"No FASTA record found in {path}")
    return name, "".join(chunks)


def _time_call(function: Callable[[], Sequence[int]]) -> Tuple[float, int, int]:
    gc.collect()
    gc.disable()
    try:
        started = time.perf_counter_ns()
        result = function()
        elapsed = (time.perf_counter_ns() - started) / 1e9
    finally:
        gc.enable()
    count = len(result)
    # Touch the result without adding a full O(n) checksum to the timing.
    raw_result = getattr(result, "data", result)
    checksum = 0 if count == 0 else int(raw_result[0]) ^ int(raw_result[-1])
    return elapsed, count, checksum


def _summary(samples: Iterable[float], count: int) -> Dict[str, float]:
    values = list(samples)
    median = statistics.median(values)
    return {
        "minimum_seconds": min(values),
        "median_seconds": median,
        "mean_seconds": statistics.mean(values),
        "stdev_seconds": statistics.stdev(values) if len(values) > 1 else 0.0,
        "maximum_seconds": max(values),
        "median_million_hashes_per_second": count / median / 1_000_000,
    }


def benchmark(sequence: str, k: int, repeats: int, seed: int) -> Dict[str, object]:
    mmh3_module = _load_mmh3()
    rng = random.Random(seed)
    report: Dict[str, object] = {
        "sequence_length": len(sequence),
        "kmer_length": k,
        "repeats": repeats,
        "legacy_mmh3_available": mmh3_module is not None,
        "modes": {},
    }

    # Warm native code, imports, and allocators before recording samples.
    warm_sequence = sequence[: max(k, min(len(sequence), 10_000))]
    _moddotplot_nthash2_hashes(warm_sequence, k, canonical=True)
    if mmh3_module is not None:
        _legacy_mmh3_hashes(warm_sequence, k, True, mmh3_module)

    for mode, canonical in (("forward", False), ("canonical", True)):
        implementations: Dict[str, Callable[[], Sequence[int]]] = {
            "nthash2": lambda canonical=canonical: _moddotplot_nthash2_hashes(
                sequence, k, canonical
            )
        }
        if mmh3_module is not None:
            implementations[
                "legacy_mmh3"
            ] = lambda canonical=canonical: _legacy_mmh3_hashes(
                sequence, k, canonical, mmh3_module
            )

        raw: Dict[str, List[float]] = {name: [] for name in implementations}
        counts: Dict[str, int] = {}
        checksums: Dict[str, int] = {}
        for _ in range(repeats):
            order = list(implementations)
            rng.shuffle(order)
            for name in order:
                elapsed, count, checksum = _time_call(implementations[name])
                raw[name].append(elapsed)
                counts[name] = count
                checksums[name] = checksum

        expected_count = max(len(sequence) - k + 1, 0)
        if any(count != expected_count for count in counts.values()):
            raise RuntimeError(
                f"Unexpected k-mer cardinality: expected {expected_count}, got {counts}"
            )

        mode_report: Dict[str, object] = {
            "hash_count": expected_count,
            "implementations": {
                name: {
                    **_summary(samples, counts[name]),
                    "raw_seconds": samples,
                    "checksum": checksums[name],
                }
                for name, samples in raw.items()
            },
        }
        if mmh3_module is not None:
            summaries = mode_report["implementations"]
            mode_report["speedup_over_legacy_mmh3"] = (
                summaries["legacy_mmh3"]["median_seconds"]
                / summaries["nthash2"]["median_seconds"]
            )
        report["modes"][mode] = mode_report

    return report


def _print_report(label: str, report: Dict[str, object]) -> None:
    print(
        f"Input: {label}; {report['sequence_length']:,} bases; "
        f"k={report['kmer_length']}; {report['repeats']} repeats"
    )
    print("ntHash2 result: compact NumPy uint64 array; legacy result: Python int list")
    if not report["legacy_mmh3_available"]:
        print("Legacy comparison: skipped (optional mmh3 is not installed)")
    print(
        f"{'mode':<10} {'implementation':<14} {'median (s)':>12} "
        f"{'Mhash/s':>10} {'speedup':>10}"
    )
    for mode, mode_report in report["modes"].items():
        speedup = mode_report.get("speedup_over_legacy_mmh3")
        for implementation, stats in mode_report["implementations"].items():
            shown_speedup = (
                f"{speedup:.2f}x"
                if implementation == "nthash2" and speedup is not None
                else "-"
            )
            print(
                f"{mode:<10} {implementation:<14} "
                f"{stats['median_seconds']:>12.6f} "
                f"{stats['median_million_hashes_per_second']:>10.2f} "
                f"{shown_speedup:>10}"
            )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--fasta", type=Path, help="benchmark the first FASTA record")
    source.add_argument(
        "--length", type=int, default=1_000_000, help="synthetic sequence length"
    )
    parser.add_argument("--kmer", type=int, default=21)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--seed", type=int, default=20260927)
    parser.add_argument("--json", type=Path, help="also write raw results as JSON")
    args = parser.parse_args(argv)

    if args.kmer <= 0:
        parser.error("--kmer must be positive")
    if args.length < 0:
        parser.error("--length cannot be negative")
    if args.repeats <= 0:
        parser.error("--repeats must be positive")

    if args.fasta:
        label, sequence = _read_first_fasta(args.fasta)
    else:
        label = f"synthetic(seed={args.seed})"
        sequence = "".join(
            random.Random(args.seed).choices(DNA_ALPHABET, k=args.length)
        )

    report = benchmark(sequence, args.kmer, args.repeats, args.seed)
    _print_report(label, report)
    if args.json:
        args.json.write_text(json.dumps(report, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
