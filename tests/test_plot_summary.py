from datetime import datetime, timezone

from moddotplot.plot_summary import PlotSummaryWriter


def test_plot_summary_records_reproducibility_metadata_and_merges_groups(tmp_path):
    first_plot = tmp_path / "chrA_FULL.png"
    second_plot = tmp_path / "chrA_TRI.svg"
    fasta = tmp_path / "source genome.fa"
    annotation = tmp_path / "features.bed"
    bedpe = tmp_path / "matrix.bedpe"
    for path in (first_plot, second_plot, fasta, annotation, bedpe):
        path.touch()

    writer = PlotSummaryWriter(
        "moddotplot -f 'source genome.fa' --region chrA:1-4000000",
        created_at=datetime(2026, 9, 28, 12, 30, tzinfo=timezone.utc),
    )
    writer.add(
        tmp_path,
        [first_plot],
        fasta_files=[fasta],
        window_sizes=[4000],
        regions=["chrA:1-4000000"],
        bed_file=annotation,
    )
    writer.add(
        tmp_path,
        [second_plot],
        window_sizes=[4000],
        bedpe_inputs=[bedpe],
    )

    summary = (tmp_path / "plot_summary.txt").read_text(encoding="utf-8")
    assert "Created: 2026-09-28T12:30:00+00:00" in summary
    assert (
        "Command: moddotplot -f 'source genome.fa' --region chrA:1-4000000" in summary
    )
    assert str(first_plot.resolve()) in summary
    assert str(second_plot.resolve()) in summary
    assert str(fasta.resolve()) in summary
    assert "4000 bp" in summary
    assert "chrA:1-4000000" in summary
    assert f"BED annotation file: {annotation.resolve()}" in summary
    assert str(bedpe.resolve()) in summary
    assert "Plot group" not in summary
    assert "Created: 2026-09-28T12:30:00+00:00\n\nPlot files:" in summary
    assert "Plot files:" in summary
    assert "\n\nFASTA files:" in summary
    assert "\n\nWindow sizes:" in summary
    assert "\n\nRegions:" in summary
    assert "\n\nBED annotation file:" in summary


def test_plot_summary_ignores_empty_plot_groups(tmp_path):
    writer = PlotSummaryWriter("moddotplot --help")

    result = writer.add(tmp_path, [], fasta_files=["source.fa"])

    assert result is None
    assert not (tmp_path / "plot_summary.txt").exists()
