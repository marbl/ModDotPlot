"""Reproducibility summaries for static plot output directories."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional


def _absolute_paths(paths: Optional[Iterable[str]]) -> list[str]:
    """Return stable, de-duplicated absolute paths in input order."""

    resolved = []
    for path in paths or ():
        absolute = str(Path(path).expanduser().resolve())
        if absolute not in resolved:
            resolved.append(absolute)
    return resolved


class PlotSummaryWriter:
    """Collect and write one human-readable summary per plot directory."""

    filename = "plot_summary.txt"

    def __init__(self, command: str, created_at: Optional[datetime] = None):
        self.command = command
        self.created_at = created_at
        self._records: dict[str, list[dict[str, object]]] = {}

    def add(
        self,
        directory: str,
        plot_files: Iterable[str],
        *,
        fasta_files: Optional[Iterable[str]] = None,
        window_sizes: Optional[Iterable[int]] = None,
        regions: Optional[Iterable[str]] = None,
        bed_file: Optional[str] = None,
        bedpe_inputs: Optional[Iterable[str]] = None,
    ) -> Optional[Path]:
        """Add a plot record and immediately refresh its directory summary.

        Empty plot records are ignored. This keeps mocked or skipped renderers
        from claiming that files were produced when they were not.
        """

        plots = _absolute_paths(plot_files)
        if not plots:
            return None

        output_directory = str(Path(directory).expanduser().resolve())
        record = {
            "created_at": (self.created_at or datetime.now().astimezone()).isoformat(
                timespec="seconds"
            ),
            "plot_files": plots,
            "fasta_files": _absolute_paths(fasta_files),
            "window_sizes": list(
                dict.fromkeys(int(size) for size in window_sizes or ())
            ),
            "regions": list(dict.fromkeys(regions or ())),
            "bed_file": _absolute_paths([bed_file])[0] if bed_file else None,
            "bedpe_inputs": _absolute_paths(bedpe_inputs),
        }
        records = self._records.setdefault(output_directory, [])
        if record not in records:
            records.append(record)
        return self._write(output_directory)

    def _write(self, directory: str) -> Path:
        output_directory = Path(directory)
        output_directory.mkdir(parents=True, exist_ok=True)
        summary_path = output_directory / self.filename
        lines = [
            "ModDotPlot Plot Summary",
            "",
            f"Command: {self.command}",
            "",
            f"Output directory: {output_directory}",
        ]

        for record in self._records[directory]:
            lines.extend(["", f"Created: {record['created_at']}"])
            self._append_list(lines, "Plot files", record["plot_files"])
            self._append_list(
                lines,
                "FASTA files",
                record["fasta_files"],
                empty="None (plot loaded from BEDPE)",
            )
            window_sizes = [f"{size} bp" for size in record["window_sizes"]]
            self._append_list(lines, "Window sizes", window_sizes, empty="Unknown")
            self._append_list(
                lines,
                "Regions",
                record["regions"],
                empty="None (whole sequence)",
            )
            lines.extend(["", f"BED annotation file: {record['bed_file'] or 'None'}"])
            if record["bedpe_inputs"]:
                self._append_list(lines, "Input BEDPE files", record["bedpe_inputs"])

        summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        return summary_path

    @staticmethod
    def _append_list(lines, label, values, *, empty=None):
        lines.extend(["", f"{label}:"])
        if values:
            lines.extend(f"  - {value}" for value in values)
        elif empty is not None:
            lines.append(f"  - {empty}")
