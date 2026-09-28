import os
from pathlib import Path
import subprocess
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _run_module(*arguments):
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(PROJECT_ROOT / "src")
    return subprocess.run(
        [sys.executable, "-m", "moddotplot", *arguments],
        cwd=PROJECT_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def test_help_does_not_import_static_rendering_stack():
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys; import moddotplot.moddotplot; "
                "assert 'moddotplot.static_plots' not in sys.modules"
            ),
        ],
        cwd=PROJECT_ROOT,
        env={**os.environ, "PYTHONPATH": str(PROJECT_ROOT / "src")},
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr


def test_module_help_is_available():
    result = _run_module("--help")

    assert result.returncode == 0
    assert "{static,interactive}" in result.stdout
    assert "static is used when omitted" in result.stdout
    assert "Static mode commands (default)" in result.stdout
    assert "Interactive mode commands (deprecated; explicit use" in result.stdout
    assert "only)" in result.stdout


def test_module_without_arguments_defaults_to_static_parser():
    result = _run_module()

    assert result.returncode == 2
    assert "static:" in result.stderr
    assert "one of the arguments -c/--config -l/--load -f/--fasta is required" in (
        result.stderr
    )
    assert "the following arguments are required: command" not in result.stderr
