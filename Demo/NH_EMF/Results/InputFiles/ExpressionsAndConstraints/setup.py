#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy
from setuptools import setup, Extension
import os


CurrentDirectiory = os.path.dirname(os.path.abspath(__file__))

setup(
    ext_modules=[Extension(CurrentDirectiory + "/_RelaxMat", [CurrentDirectiory + "/_Rates.c", CurrentDirectiory + "/Rates.c"])],
    include_dirs=numpy.get_include(),
)

setup(
    ext_modules=[Extension(CurrentDirectiory + "/Rates/_R1calculation", [CurrentDirectiory + "/_Rates.c", CurrentDirectiory + "/Rates.c"])],
    include_dirs=numpy.get_include(),
)

setup(
    ext_modules=[Extension(CurrentDirectiory + "/Rates/_R2calculation", [CurrentDirectiory + "/_Rates.c", CurrentDirectiory + "/Rates.c"])],
    include_dirs=numpy.get_include(),
)

setup(
    ext_modules=[Extension(CurrentDirectiory + "/Rates/_Sigmacalculation", [CurrentDirectiory + "/_Rates.c", CurrentDirectiory + "/Rates.c"])],
    include_dirs=numpy.get_include(),
)



#run python setup.py build_ext --inplace