import importlib.util
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_PATH = PROJECT_ROOT / "benchmarks" / "benchmark_hashing.py"


def _load_benchmark_module():
    spec = importlib.util.spec_from_file_location(
        "moddotplot_hash_benchmark", BENCHMARK_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_benchmark_runs_without_optional_legacy_dependency(monkeypatch):
    benchmark = _load_benchmark_module()
    monkeypatch.setattr(benchmark, "_load_mmh3", lambda: None)

    report = benchmark.benchmark("ACGT" * 50, k=5, repeats=1, seed=7)

    assert report["legacy_mmh3_available"] is False
    for mode in ("forward", "canonical"):
        mode_report = report["modes"][mode]
        assert mode_report["hash_count"] == 196
        assert set(mode_report["implementations"]) == {"nthash2"}
        assert "speedup_over_legacy_mmh3" not in mode_report


def test_benchmark_reports_mmh3_only_when_optional_module_is_available(monkeypatch):
    benchmark = _load_benchmark_module()

    class FakeMmh3:
        @staticmethod
        def hash(value):
            return sum(map(ord, value))

    monkeypatch.setattr(benchmark, "_load_mmh3", lambda: FakeMmh3())

    report = benchmark.benchmark("ACGT" * 50, k=5, repeats=1, seed=7)

    assert report["legacy_mmh3_available"] is True
    for mode in ("forward", "canonical"):
        mode_report = report["modes"][mode]
        assert set(mode_report["implementations"]) == {
            "nthash2",
            "legacy_mmh3",
        }
        assert mode_report["speedup_over_legacy_mmh3"] > 0


def test_benchmark_cli_clearly_reports_skipped_legacy_comparison(monkeypatch, capsys):
    benchmark = _load_benchmark_module()
    monkeypatch.setattr(benchmark, "_load_mmh3", lambda: None)

    assert benchmark.main(["--length", "100", "--kmer", "5", "--repeats", "1"]) == 0

    output = capsys.readouterr().out
    assert "Legacy comparison: skipped (optional mmh3 is not installed)" in output
    assert "forward" in output
    assert "canonical" in output
    assert "nthash2" in output
