from glob import glob
from xml.etree.ElementInclude import include
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules =  [Pybind11Extension(
    "DPM",
    sorted(glob("src/*.cpp")),
    include_dirs=["include/Cell/","/usr/include/c++/11/"],
    language="c++"),
]

setup(ext_modules=ext_modules, cmdclass={"build_ext":build_ext})