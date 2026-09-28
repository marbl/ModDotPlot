import matplotlib.pyplot as plt
from plotnine import ggplot

import moddotplot.static_plots as static_plots
from moddotplot.native_render import (
    DEFAULT_FONT_FAMILY,
    FALLBACK_FONT_FAMILY,
    MIN_TEXT_SIZE,
    MIN_TITLE_SIZE,
    clamped_font_size,
    save_figure_pair,
)


def test_width_scaled_fonts_have_readable_minimums():
    assert clamped_font_size(1, 1.0) == MIN_TEXT_SIZE
    assert clamped_font_size(1, 1.4, MIN_TITLE_SIZE) == MIN_TITLE_SIZE
    assert clamped_font_size(9, 1.4, MIN_TITLE_SIZE) == 12.6


def test_native_outputs_default_to_helvetica(tmp_path):
    figure, axis = plt.subplots()
    title = axis.set_title("Helvetica title")
    try:
        save_figure_pair(figure, tmp_path / "helvetica", "svg", 72)
        assert title.get_fontfamily() == [DEFAULT_FONT_FAMILY]
        assert (tmp_path / "helvetica.png").stat().st_size > 0
        assert (tmp_path / "helvetica.svg").stat().st_size > 0
    finally:
        plt.close(figure)


def test_native_outputs_retry_with_dejavu_on_glyph_failure(tmp_path, monkeypatch):
    figure, axis = plt.subplots()
    title = axis.set_title("Fallback title")
    original_savefig = figure.savefig
    attempted_families = []

    def fail_for_helvetica(*args, **kwargs):
        family = title.get_fontfamily()[0]
        attempted_families.append(family)
        if family == DEFAULT_FONT_FAMILY:
            raise RuntimeError("failed to load glyph")
        return original_savefig(*args, **kwargs)

    monkeypatch.setattr(figure, "savefig", fail_for_helvetica)
    try:
        save_figure_pair(figure, tmp_path / "fallback", "svg", 72)
        assert attempted_families == [
            DEFAULT_FONT_FAMILY,
            FALLBACK_FONT_FAMILY,
            FALLBACK_FONT_FAMILY,
        ]
        assert title.get_fontfamily() == [FALLBACK_FONT_FAMILY]
        assert (tmp_path / "fallback.png").stat().st_size > 0
        assert (tmp_path / "fallback.svg").stat().st_size > 0
    finally:
        plt.close(figure)


def test_plotnine_outputs_retry_with_dejavu_on_glyph_failure(monkeypatch):
    attempted_families = []

    def fail_for_helvetica(plot, **_kwargs):
        family = plot.theme.getp(("text", "family"))[0]
        attempted_families.append(family)
        if family == DEFAULT_FONT_FAMILY:
            raise RuntimeError("failed to load glyph")

    monkeypatch.setattr(static_plots, "ggsave", fail_for_helvetica)
    static_plots._save_plot(ggplot(), filename="unused.png")

    assert attempted_families == [DEFAULT_FONT_FAMILY, FALLBACK_FONT_FAMILY]
