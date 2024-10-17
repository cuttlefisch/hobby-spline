#! /usr/bin/env python
import os
import sys

import difflib
import unittest
from hobby_spline.splines import hobby_points
from euclid3 import Point2, Point3, Vector2, Vector3
from math import pi

SEGMENTS = 8


class DiffOutput(unittest.TestCase):

    def assertEqual(self, first, second, msg=None):
        """
        Override assertEqual and print(a context diff if msg=None)
        """
        # Test if both are strings, in Python 2 & 3
        string_types = str if sys.version_info[0] == 3 else basestring

        if isinstance(first, string_types) and isinstance(second, string_types):
            if not msg:
                msg = "Strings are not equal:\n" + "".join(
                    difflib.unified_diff(
                        [first], [second], fromfile="actual", tofile="expected"
                    )
                )
        return super(DiffOutput, self).assertEqual(first, second, msg=msg)


class TestSplines(DiffOutput):
    def setUp(self):
        self.points = [
            Point3(0, 0),
            Point3(1, 1),
            Point3(2, 1),
        ]
        self.points_raw = [
            (0, 0),
            (1, 1),
            (2, 1),
        ]
        self.bezier_controls = [
            Point3(0, 0),
            Point3(1, 1),
            Point3(2, 1),
            Point3(2, -1),
        ]
        self.bezier_controls_raw = [(0, 0), (1, 1), (2, 1), (2, -1)]
        self.hobby_omega = 0.69
        self.subdivisions = 2

    def assertPointsListsEqual(self, a, b):
        str_list = lambda x: list(str(v) for v in x)
        self.assertEqual(str_list(a), str_list(b))

    def test_hobby_points(self):
        expected = [
            Point3(0.00, 0.00, 0.00),
            Point3(0.07, 0.50, 0.00),
            Point3(0.52, 0.82, 0.00),
            Point3(1.00, 1.00, 0.00),
            Point3(1.33, 1.12, 0.00),
            Point3(1.71, 1.19, 0.00),
            Point3(2.00, 1.00, 0.00),
            Point3(2.63, 0.60, 0.00),
            Point3(2.35, -0.28, 0.00),
            Point3(2.00, -1.00, 0.00),
        ]

        actual = hobby_points(self.bezier_controls, self.hobby_omega, close_loop=False)
        self.assertPointsListsEqual(expected, actual)

    def test_hobby_points_raw(self):
        expected = [
            Point3(0.00, 0.00, 0.00),
            Point3(0.07, 0.50, 0.00),
            Point3(0.52, 0.82, 0.00),
            Point3(1.00, 1.00, 0.00),
            Point3(1.33, 1.12, 0.00),
            Point3(1.71, 1.19, 0.00),
            Point3(2.00, 1.00, 0.00),
            Point3(2.63, 0.60, 0.00),
            Point3(2.35, -0.28, 0.00),
            Point3(2.00, -1.00, 0.00),
        ]
        actual = hobby_points(
            self.bezier_controls_raw, self.hobby_omega, close_loop=False
        )
        self.assertPointsListsEqual(expected, actual)

    def test_hobby_points_3d(self):
        controls_3d = [
            Point3(-2, -1, 0),
            Point3(-0.5, -0.5, 1),
            Point3(0.5, 0.5, 1),
            Point3(2, 1, 0),
        ]
        expected = [
            Point3(-2.00, -1.00, 0.00),
            Point3(-1.75, -1.03, 0.62),
            Point3(-1.04, -0.90, 0.86),
            Point3(-0.50, -0.50, 1.00),
            Point3(-0.12, -0.22, 1.10),
            Point3(0.11, 0.24, 1.13),
            Point3(0.50, 0.50, 1.00),
            Point3(1.01, 0.84, 0.83),
            Point3(1.55, 0.88, 0.44),
            Point3(2.00, 1.00, 0.00),
        ]

        actual = hobby_points(controls_3d, self.hobby_omega, close_loop=False)
        self.assertPointsListsEqual(expected, actual)
