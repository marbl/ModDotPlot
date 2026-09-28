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
    assert "{interactive,static}" in result.stdout


def test_module_without_subcommand_has_clean_usage_error():
    result = _run_module()

    assert result.returncode == 2
    assert "the following arguments are required: command" in result.stderr
