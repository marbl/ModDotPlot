import os
from setuptools import setup, find_packages

from moddotplot.const import VERSION

SCRIPT = "./moddotplot/__main__.py"
print("Setup of version " + VERSION + " of: " + os.path.abspath(SCRIPT))

setup(
    name="moddotplot",
    version=VERSION,
    url="https://github.com/marbl/ModDotPlot",
    packages=find_packages(),
    install_requires=[
        "pysam",
        "pandas",
        "plotly",
        "matplotlib",
        "dash",
        "plotnine",
        "patchworklib"
    ],
    python_requires=">=3.7",
    author_email="asweeten@cs.jhu.edu",
    entry_points={
        "console_scripts": [
            "moddotplot = moddotplot.__main__:main",
        ],
    },
)
