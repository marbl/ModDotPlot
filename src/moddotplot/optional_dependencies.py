"""Shared errors and installation guidance for optional features."""

INSTALL_COMMANDS = {
    "Cooler export": 'python -m pip install "ModDotPlot[cooler]"',
    "Interactive mode": 'python -m pip install "ModDotPlot[interactive]"',
}


class OptionalDependencyError(ImportError):
    """Raised when a requested feature is missing its optional dependencies."""

    def __init__(self, feature):
        install_command = INSTALL_COMMANDS.get(
            feature, 'python -m pip install "ModDotPlot[interactive]"'
        )
        super().__init__(
            f"{feature} requires optional dependencies. "
            f"Install them with: {install_command}"
        )
