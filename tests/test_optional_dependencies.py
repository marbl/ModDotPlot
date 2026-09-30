import builtins
import sys

import numpy as np
import pytest

import moddotplot.interactive as interactive
import moddotplot.moddotplot as cli
from moddotplot.estimate_identity import convertMatrixToCool, require_cooler_dependency
from moddotplot.optional_dependencies import OptionalDependencyError


INSTALL_HINT = 'python -m pip install "ModDotPlot[interactive]"'


def _block_import(monkeypatch, missing_package):
    real_import = builtins.__import__

    def import_without_optional_package(
        name, globals=None, locals=None, fromlist=(), level=0
    ):
        if name == missing_package or name.startswith(f"{missing_package}."):
            raise ModuleNotFoundError(
                f"No module named '{missing_package}'", name=missing_package
            )
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", import_without_optional_package)


def test_missing_cooler_reports_the_interactive_extra_install_command(monkeypatch):
    _block_import(monkeypatch, "cooler")

    with pytest.raises(OptionalDependencyError) as exc_info:
        require_cooler_dependency()

    assert "Cooler export requires" in str(exc_info.value)
    assert INSTALL_HINT in str(exc_info.value)
    assert isinstance(exc_info.value.__cause__, ModuleNotFoundError)


@pytest.mark.parametrize("missing_package", ["dash", "plotly"])
def test_missing_interactive_package_reports_the_extra_install_command(
    monkeypatch, missing_package
):
    _block_import(monkeypatch, missing_package)

    with pytest.raises(OptionalDependencyError) as exc_info:
        interactive.require_interactive_dependencies()

    assert "Interactive mode requires" in str(exc_info.value)
    assert INSTALL_HINT in str(exc_info.value)
    assert isinstance(exc_info.value.__cause__, ModuleNotFoundError)


def test_static_cooler_request_exits_with_actionable_error(monkeypatch, capsys):
    def missing_cooler():
        raise OptionalDependencyError("Cooler export")

    monkeypatch.setattr(cli, "require_cooler_dependency", missing_cooler)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "moddotplot",
            "static",
            "--fasta",
            "sequence.fa",
            "--cooler",
            "--no-plot",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    assert exc_info.value.code == 2
    assert INSTALL_HINT in capsys.readouterr().err


def test_cooler_extra_writes_a_readable_comparative_matrix(tmp_path):
    cooler = pytest.importorskip("cooler")
    output = tmp_path / "comparison.cool"

    convertMatrixToCool(
        matrix=np.asarray([[1.0, 0.9], [0.9, 1.0]]),
        window_size=10,
        id_threshold=80,
        x_name="chr1",
        y_name="chr2",
        self_identity=False,
        x_offset=0,
        y_offset=0,
        chromsizes={"chr1": 20, "chr2": 20},
        output_cool=str(output),
    )

    matrix = cooler.Cooler(str(output))
    assert tuple(matrix.shape) == (4, 4)
    assert matrix.info["nnz"] == 4


def test_interactive_request_exits_with_actionable_error(monkeypatch, capsys):
    def missing_interactive_dependencies():
        raise OptionalDependencyError("Interactive mode")

    monkeypatch.setattr(cli, "_load_interactive_plotting", lambda: None)
    monkeypatch.setattr(
        cli, "require_interactive_dependencies", missing_interactive_dependencies
    )
    monkeypatch.setattr(
        sys,
        "argv",
        ["moddotplot", "interactive", "--fasta", "sequence.fa"],
    )

    with pytest.raises(SystemExit) as exc_info:
        cli.main()

    assert exc_info.value.code == 2
    assert INSTALL_HINT in capsys.readouterr().err
