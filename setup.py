#!/usr/bin/env python3

from setuptools import setup, find_packages
import os

setup(
    name="hobby_spline",
    version="0.0.1",
    description="Fit a cubic bezier curve in 2D and 3D space using Hobby's Algorithm.",
    author="cuttlefisch",
    author_email="",
    packages=find_packages(),
    test_suite="nose.collector",
    tests_require=["nose==1.3.7", "setuptools<72"],
    python_requires="==3.12",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)
