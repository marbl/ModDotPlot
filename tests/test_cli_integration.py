import gzip
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _run_cli(*arguments):
    environment = os.environ.copy()
    environment.update(
        {
            "PYTHONPATH": str(PROJECT_ROOT / "src"),
            "MPLBACKEND": "Agg",
        }
    )
    return subprocess.run(
        [sys.executable, "-m", "moddotplot", *map(str, arguments)],
        cwd=PROJECT_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def _write_multifasta(path):
    path.write_text(
        ">alpha\n"
        + "ACGT" * 300
        + "\n>beta\n"
        + "ACGT" * 275
        + "\n>gamma\n"
        + "ACGT" * 250
        + "\n"
    )


def test_static_cli_computes_all_self_and_pairwise_outputs(tmp_path):
    fasta = tmp_path / "three.fa"
    output = tmp_path / "static"
    _write_multifasta(fasta)

    result = _run_cli(
        "static",
        "--fasta",
        fasta,
        "--window",
        100,
        "--modimizer",
        10,
        "--identity",
        80,
        "--compare",
        "--no-plot",
        "--output-dir",
        output,
    )

    assert result.returncode == 0, result.stderr + result.stdout
    bedpe_files = sorted(path.relative_to(output) for path in output.rglob("*.bedpe"))
    assert bedpe_files == [
        Path("alpha/alpha.bedpe"),
        Path("alpha_beta/alpha_beta_COMPARE.bedpe"),
        Path("alpha_gamma/alpha_gamma_COMPARE.bedpe"),
        Path("beta/beta.bedpe"),
        Path("beta_gamma/beta_gamma_COMPARE.bedpe"),
        Path("gamma/gamma.bedpe"),
    ]
    assert all(path.stat().st_size > 0 for path in output.rglob("*.bedpe"))


def test_static_grid_regions_with_dotted_headers_crop_every_output(tmp_path):
    names = [
        "PAN010.chr14.haplotype1.paternal",
        "PAN010.chr14.haplotype2.maternal",
        "PAN027.chr14.paternal",
    ]
    fasta = tmp_path / "three-dotted.fa"
    fasta.write_text(
        "".join(f">{name}\n{'ACGT' * 300}\n" for name in names),
        encoding="ascii",
    )
    output = tmp_path / "regions"

    result = _run_cli(
        "static",
        "--grid",
        "--fasta",
        fasta,
        "--region",
        *(f"{name}:1-400" for name in names),
        "--window",
        50,
        "--modimizer",
        10,
        "--identity",
        80,
        "--no-plot",
        "--output-dir",
        output,
    )

    assert result.returncode == 0, result.stderr + result.stdout
    assert result.stdout.count("Sequence length n: 400") == 3
    assert "Sequence length n: 1200" not in result.stdout
    bedpe_files = list(output.rglob("*.bedpe"))
    assert len(bedpe_files) == 6
    for bedpe in bedpe_files:
        rows = [line.split("\t") for line in bedpe.read_text().splitlines()[1:]]
        assert rows
        coordinates = [int(row[index]) for row in rows for index in (1, 2, 4, 5)]
        assert min(coordinates) >= 1
        assert max(coordinates) <= 400

    assert (output / "3x3_GRID.png").is_file()
    assert (output / "3x3_GRID.svg").is_file()


def test_static_cli_reads_gzip_and_renders_bed_annotations(tmp_path):
    fasta = tmp_path / "one.fa.gz"
    with gzip.open(fasta, "wt", encoding="ascii") as stream:
        stream.write(">alpha description\n" + "ACGT" * 300 + "\n")

    annotation = tmp_path / "annotations.bed"
    annotation.write_text("alpha\t100\t300\tfeature\t0\t+\t100\t300\t12,34,56\n")
    output = tmp_path / "annotated"

    result = _run_cli(
        "static",
        "--fasta",
        fasta,
        "--bed",
        annotation,
        "--window",
        100,
        "--modimizer",
        10,
        "--identity",
        80,
        "--no-hist",
        "--output-dir",
        output,
    )

    assert result.returncode == 0, result.stderr + result.stdout
    sequence_output = output / "alpha"
    expected = [
        sequence_output / "alpha_ANNOTATION_TRACK.svg",
        sequence_output / "alpha_ANNOTATION_TRACK.png",
        sequence_output / "alpha_TRI_ANNOTATED.svg",
        sequence_output / "alpha_TRI_ANNOTATED.png",
    ]
    assert all(path.stat().st_size > 0 for path in expected)
    ET.parse(sequence_output / "alpha_ANNOTATION_TRACK.svg")
    ET.parse(sequence_output / "alpha_TRI_ANNOTATED.svg")
    assert not list(sequence_output.glob("*.ini"))
    assert not (tmp_path / "one.fa.gz.fai").exists()


def test_interactive_cli_forward_mode_saves_matrix_without_launching_server(tmp_path):
    fasta = tmp_path / "one.fa"
    fasta.write_text(">alpha\n" + "ACGT" * 300 + "\n")
    output = tmp_path / "interactive"

    result = _run_cli(
        "interactive",
        "--fasta",
        fasta,
        "--window",
        100,
        "--resolution",
        10,
        "--modimizer",
        10,
        "--quick",
        "--forward",
        "--save",
        "--no-plot",
        "--output-dir",
        output,
    )

    assert result.returncode == 0, result.stderr + result.stdout
    saved = output / "interactive_matrices"
    assert (saved / "alpha_0.npz").is_file()
    assert (saved / "metadata.pkl").is_file()
    assert "Saved matrices" in result.stdout
