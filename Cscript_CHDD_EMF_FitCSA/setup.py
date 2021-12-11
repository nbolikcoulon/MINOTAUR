#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
import numpy.distutils.misc_util
import os


CurrentDirectiory = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(CurrentDirectiory + "/Rates"):
    os.makedirs(CurrentDirectiory + "/Rates")

setup(
    ext_modules=[Extension(CurrentDirectiory + "/_RelaxMat", [CurrentDirectiory + "/_Rates.c", CurrentDirectiory + "/Rates.c"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)

setup(
    ext_modules=[Extension(CurrentDirectiory + "/Rates/_R2calculation", [CurrentDirectiory + "/_Rates.c", CurrentDirectiory + "/Rates.c"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)

setup(
    ext_modules=[Extension(CurrentDirectiory + "/Rates/_R1calculation", [CurrentDirectiory + "/_Rates.c", CurrentDirectiory + "/Rates.c"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)

setup(
    ext_modules=[Extension(CurrentDirectiory + "/Rates/_Sigmacalculation", [CurrentDirectiory + "/_Rates.c", CurrentDirectiory + "/Rates.c"])],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)



#run python setup.py build_ext --inplace