#!/usr/bin/env python3
import sys
from moddotplot.estimate_identity import *
from moddotplot.moddotplot import main
from moddotplot.parse_fasta import *
import setproctitle

if __name__ == "__main__":
    # Keeping execution behind the standard guard makes this module safe to
    # import in spawn-based chromosome worker processes and as a console-script
    # entry point.
    setproctitle.setproctitle("ModDotPlot")
    sys.exit(main())
