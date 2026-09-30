import sys

from setuptools import Extension, setup

if sys.platform == "win32":
    compile_args = ["/O2", "/std:c++17"]
else:
    compile_args = ["-O3", "-std=c++17"]


setup(
    exclude_package_data={"moddotplot": ["*.cpp"]},
    ext_modules=[
        Extension(
            "moddotplot._nthash",
            sources=["src/moddotplot/_nthash.cpp"],
            define_macros=[("Py_LIMITED_API", "0x03080000")],
            extra_compile_args=compile_args,
            language="c++",
            py_limited_api=True,
        )
    ],
)
