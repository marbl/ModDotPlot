from setuptools import setup

import os
import sys
from setuptools import setup, find_packages

from const import VERSION

SCRIPT = "./moddotplot/__main__.py"
print("Setup of version " + VERSION + " of: " + os.path.abspath(SCRIPT))

TEST_REQS = ["mock", "pytest", "doctest"]
setup(
    name="moddotplot",
    version=VERSION,
    url="https://github.com/marbl/ModDotPlot",
    packages=find_packages(),
    install_requires=[
        "screed",
        "pandas",
    ],
    python_requires=">=3",
    author_email="asweeten@nih.gov",
    entry_points={
        'console_scripts': [
            'moddotplot = moddotplot.__main__:main',
        ],
    },
    test_suite="nose.collector",
    tests_require=TEST_REQS,
)