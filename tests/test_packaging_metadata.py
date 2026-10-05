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


def test_declared_python_floor_matches_classifiers():
    metadata = project_metadata()
    assert metadata["requires-python"] == ">=3.10,<3.15"
    for minor in range(10, 15):
        assert f"Programming Language :: Python :: 3.{minor}" in metadata["classifiers"]


def test_core_dependencies_are_the_required_numeric_and_plotting_stack():
    dependencies = project_metadata()["dependencies"]
    names = {
        dependency.split(";", 1)[0]
        .split("=", 1)[0]
        .split("<", 1)[0]
        .split(">", 1)[0]
        .strip()
        .lower()
        for dependency in dependencies
    }

    assert names == {"numpy", "pandas", "scipy", "matplotlib"}
    assert "matplotlib>=3.10.9; python_version < '3.11'" in dependencies
    assert "matplotlib>=3.11.2; python_version >= '3.11'" in dependencies


def test_interactive_dependencies_are_not_installed_with_the_core_package():
    metadata = project_metadata()
    core_dependencies = metadata["dependencies"]
    core_names = {
        dependency.split(";", 1)[0]
        .split("[", 1)[0]
        .split("=", 1)[0]
        .split("<", 1)[0]
        .split(">", 1)[0]
        .strip()
        .lower()
        for dependency in core_dependencies
    }

    assert core_names.isdisjoint(
        {"cooler", "dash", "plotly", "pillow", "plotnine", "palettable", "setproctitle"}
    )
    assert set(metadata["optional-dependencies"]["interactive"]) == {
        "dash>=2.9",
        "plotly",
    }
    assert metadata["optional-dependencies"]["cooler"] == ["cooler"]


def test_python_310_test_dependencies_include_tomli():
    test_dependencies = project_metadata()["optional-dependencies"]["test"]
    assert "tomli>=1.1; python_version < '3.11'" in test_dependencies


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


def test_colorbrewer_license_is_included_in_distributions():
    setup_config = (PROJECT_ROOT / "setup.cfg").read_text()
    assert "THIRD_PARTY_LICENSES.md" in setup_config


def test_ci_covers_every_supported_python_minor():
    workflow = (PROJECT_ROOT / ".github/workflows/ci.yml").read_text()
    for minor in range(10, 15):
        assert f'          - "3.{minor}"' in workflow


def test_package_ci_compares_python_specifiers_semantically():
    workflow = (PROJECT_ROOT / ".github/workflows/ci.yml").read_text()

    assert "SpecifierSet(actual['Requires-Python'])" in workflow
    assert "SpecifierSet(expected['requires-python'])" in workflow


def test_full_test_workflows_install_all_optional_test_dependencies():
    ci_workflow = (PROJECT_ROOT / ".github/workflows/ci.yml").read_text()
    release_workflow = (
        PROJECT_ROOT / ".github/workflows/publish-to-pypi.yml"
    ).read_text()

    assert (
        'python -m pip install --editable ".[test,interactive,cooler]"' in ci_workflow
    )
    assert 'python -m pip install ".[test,interactive,cooler]"' in release_workflow


def test_release_workflow_is_tag_gated_and_uses_trusted_publishing():
    workflow = (PROJECT_ROOT / ".github/workflows/publish-to-pypi.yml").read_text()
    setup_config = (PROJECT_ROOT / "setup.cfg").read_text()

    assert '      - "v*"' in workflow
    assert "Verify tag matches the package version" in workflow
    assert "id-token: write" in workflow
    assert "pypa/gh-action-pypi-publish@release/v1" in workflow
    assert "pypa/cibuildwheel@v4.2.0" in workflow
    assert 'CIBW_BUILD: "cp310-*"' in workflow
    assert "py_limited_api = cp38" in setup_config
    for runner in ("ubuntu-latest", "macos-15-intel", "windows-latest"):
        assert f"          - {runner}" in workflow
    assert "      - sdist" in workflow
    assert "      - wheels" in workflow
