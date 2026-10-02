#!/usr/bin/env python3
from contextlib import redirect_stderr, redirect_stdout
import os
import sys


def _run(arguments):
    import setproctitle

    from moddotplot.moddotplot import main as run_moddotplot

    setproctitle.setproctitle("ModDotPlot")
    return run_moddotplot(arguments)


def main(arguments=None):
    """Load the CLI under quiet redirection to cover import-time output."""

    raw_arguments = list(sys.argv[1:] if arguments is None else arguments)
    if "--quiet" not in raw_arguments:
        return _run(raw_arguments)

    with open(os.devnull, "w", encoding="utf-8") as sink:
        with redirect_stdout(sink), redirect_stderr(sink):
            return _run(raw_arguments)


if __name__ == "__main__":
    sys.exit(main())
