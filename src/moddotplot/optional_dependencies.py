"""Shared errors and installation guidance for optional features."""

INTERACTIVE_INSTALL_COMMAND = 'python -m pip install "ModDotPlot[interactive]"'


class OptionalDependencyError(ImportError):
    """Raised when a requested feature is missing its optional dependencies."""

    def __init__(self, feature):
        super().__init__(
            f"{feature} requires the optional interactive dependencies. "
            f"Install them with: {INTERACTIVE_INSTALL_COMMAND}"
        )
