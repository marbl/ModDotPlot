from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover
    # Retained for source-tree tooling that may run outside supported Python.
    import tomli as tomllib

from moddotplot.const import VERSION

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def project_metadata():
    with (PROJECT_ROOT / "pyproject.toml").open("rb") as pyproject:
        return tomllib.load(pyproject)["project"]


def test_runtime_and_distribution_versions_match():
    assert project_metadata()["version"] == VERSION


def test_declared_python_floor_matches_documentation():
    assert project_metadata()["requires-python"] == ">=3.11,<3.15"
    readme = (PROJECT_ROOT / "README.md").read_text()
    assert "supports Python 3.11 through 3.14" in readme


def test_plotnine_supports_declared_python_floor():
    dependencies = project_metadata()["dependencies"]
    assert "plotnine>=0.15.8,<0.16" in dependencies
    assert "matplotlib>=3.11.2" in dependencies


def test_mmh3_is_not_a_runtime_dependency():
    dependencies = project_metadata()["dependencies"]
    normalized_names = {
        dependency.split(";", 1)[0]
        .split("[", 1)[0]
        .split("=", 1)[0]
        .split("<", 1)[0]
        .split(">", 1)[0]
        .strip()
        .lower()
        for dependency in dependencies
    }

    assert "mmh3" not in normalized_names


def test_replaced_compiled_and_genome_track_dependencies_are_not_runtime_dependencies():
    dependencies = project_metadata()["dependencies"]
    normalized_names = {
        dependency.split(";", 1)[0]
        .split("[", 1)[0]
        .split("=", 1)[0]
        .split("<", 1)[0]
        .split(">", 1)[0]
        .strip()
        .lower()
        for dependency in dependencies
    }

    assert "pysam" not in normalized_names
    assert "pygenometracks" not in normalized_names
    assert "matplotlib" in normalized_names


def test_svg_composition_dependencies_are_not_runtime_dependencies():
    dependencies = project_metadata()["dependencies"]
    normalized_names = {
        dependency.split(";", 1)[0]
        .split("[", 1)[0]
        .split("=", 1)[0]
        .split("<", 1)[0]
        .split(">", 1)[0]
        .strip()
        .lower()
        for dependency in dependencies
    }

    assert "cairosvg" not in normalized_names
    assert "svgutils" not in normalized_names
    assert "patchworklib" not in normalized_names
    assert "matplotlib" in normalized_names


def test_ci_covers_every_supported_python_minor():
    workflow = (PROJECT_ROOT / ".github/workflows/ci.yml").read_text()
    for minor in range(11, 15):
        assert f'          - "3.{minor}"' in workflow
    assert '          - "3.10"' not in workflow


def test_release_workflow_is_tag_gated_and_uses_trusted_publishing():
    workflow = (PROJECT_ROOT / ".github/workflows/publish-to-pypi.yml").read_text()
    setup_config = (PROJECT_ROOT / "setup.cfg").read_text()

    assert '      - "v*"' in workflow
    assert "Verify tag matches the package version" in workflow
    assert "id-token: write" in workflow
    assert "pypa/gh-action-pypi-publish@release/v1" in workflow
    assert "pypa/cibuildwheel@v4.2.0" in workflow
    assert 'CIBW_BUILD: "cp311-*"' in workflow
    assert "py_limited_api = cp38" in setup_config
    for runner in ("ubuntu-latest", "macos-15-intel", "windows-latest"):
        assert f"          - {runner}" in workflow
    assert "      - sdist" in workflow
    assert "      - wheels" in workflow
