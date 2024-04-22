#!/usr/bin/env python3
import sys
from moddotplot.estimate_identity import *
from moddotplot.moddotplot import main
from moddotplot.parse_fasta import *
import setproctitle

# Set the process title to a custom name
setproctitle.setproctitle("ModDotPlot")

sys.exit(main())
