from math import copysign, pow, cos, sin, pi, atan2, acos
from itertools import pairwise

from euclid3 import Vector2, Vector3, Point2, Point3, Matrix4

from typing import Sequence, Tuple, Union, List, cast

# ================
# = HOBBY CURVES =
# ================
# Based on this post
# https://www.jakelow.com/blog/hobby-curves
# Comments mostly transferred to respective positions
# This code implements Hobby's algorithm for fitting a cubic Bézier curve onto
# a sequence of points in two dimensions.

# Hobby's algorithm was devised by John D. Hobby and Donald Knuth for use in
# the METAFONT program. Their original paper on the algorithm, titled “Smooth,
# Easy to Compute Interpolating Splines”, was published in 'Discrete and
# Computational Geometry' vol. 1 in 1986. A copy can be found at this URL:
# https://link.springer.com/content/pdf/10.1007/BF02187690.pdf


# Based on the paper "Typographers, programmers and mathematicians, or the case
# of an æsthetically pleasing interpolation" by Bogusław Jackowski, published in
# TUGboat, Volume 34 (2013), No. 2, and available at the following URL:
# https://tug.org/TUGboat/tb34-2/tb107jackowski.pdf
def hobby_points(
    points: List[Point3], omega: float = 0, close_loop: bool = False
) -> List[Point3]:
    """Hobby's algorithm: given a set of points, fit a Bézier spline to them.

    The chosen splines tend to have pleasing, rounded shapes.

    Parameters:
    - points: an array of points as [x, y] pairs
    - omega: a number between 0 and 1 (inclusive); controls how much curl
        there will be at the endpoints of the curve

    Returns: an array of points as [x, y] pairs that define a Bézier spline.

    The output array will have 3n - 2 points where n is the number of input points.
    The output will contain every point in the input (these become knots in the
    Bézier spline), interspersed with pairs of new points which define positions
    of handle points on the spline. All points are in the same coordinate space.
    """
    # n is defined such that the points can be numbered P[0] ... P[n].
    # such that there are a total of n-1 points
    # Ensure points are 3d
    if not isinstance(points[0], Point3):
        points = [Point3(*point) for point in points]

    # Track edge control points to smooth things out
    points.append(points[0])
    points.insert(0, points[-2])

    n = len(points) - 1
    if n < 2:
        raise ValueError(
            "Solving for Hobby points only possible for two or more points."
        )
    # chords[i] is the vector from P[i] to P[i+1].
    # d[i] is the length of the ith chord.
    # Chords from first to last point
    chords = [
        Vector3(*(next_point - point).xyz) for point, next_point in pairwise(points)
    ]
    d = [chord.magnitude() for chord in chords]
    if min(d) <= 0:
        # cannot support successive points being the same
        raise ValueError("No chord can be zero-length.")

    # gamma[i] is the signed turning angle at P[i], i.e. the angle between
    # the chords from P[i-1] to P[i] and from P[i] to P[i+1].
    gamma = [chord.angle(next_chord) for chord, next_chord in pairwise(chords)]
    gamma.insert(0, gamma[-1])
    gamma.append(gamma[1])

    normals = [
        vec3_perp_norm(chord, next_chord) for chord, next_chord in pairwise(chords)
    ]
    normals.insert(0, normals[0])
    normals.append(normals[-1])

    # Set up the system of linear equations (Jackowski, formula 38).
    # We're representing this system as a tridiagonal matrix, because
    # we can solve such a system in O(n) time using the Thomas algorithm.
    #
    # Here, A, B, and C are the matrix diagonals and D is the right-hand side.
    # See Wikipedia for a more detailed explanation:
    # https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    A = [0.0]
    B = [2.0 + omega]
    C = [2.0 * omega + 1]
    D = [-1.0 * C[0] * gamma[1]]

    for i in range(1, n):
        A.append(1.0 / d[i - 1])
        B.append((2.0 * d[i - 1] + 2 * d[i]) / (d[i - 1] * d[i]))
        C.append(1.0 / d[i])
        D.append(
            (-1.0 * (2 * gamma[i] * d[i] + gamma[i + 1] * d[i - 1])) / (d[i - 1] * d[i])
        )
    A.append(2.0 * omega + 1)
    B.append(2.0 + omega)
    C.append(0)
    D.append(0)

    alpha = thomas(A, B, C, D)

    # Use alpha (the chord angle) and gamma (the turning angle of the chord
    # polyline) to solve for beta at each point (beta is like alpha, but for
    # the chord and handle vector arriving at P[i] rather than leaving from it).
    beta = []
    for i in range(n - 1):
        beta.append(-1 * gamma[i + 1] - alpha[i + 1])
    beta.append(-1 * alpha[n])

    # Now we have the angles between the handle vector and the chord
    # both arrriving at and leaving from each point, we can solve for the
    # positions of the handle (control) points themselves.
    c0 = []
    c1 = []
    for i in range(n):
        # Compute the magnitudes of the handle vectors at this point
        a = (rho(alpha[i], beta[i]) * d[i]) / 3.0
        b = (rho(beta[i], alpha[i]) * d[i]) / 3.0

        # Use the magnitutes, chords, and turning angles to find
        # the positions of the control points in the global coord space.
        c0_vec = chords[i].rotate_around(normals[i], alpha[i])
        c1_vec = chords[i].rotate_around(normals[i + 1], -1 * beta[i])

        c0_new_point = points[i] + c0_vec.normalize() * a
        c1_new_point = points[i + 1] - c1_vec.normalize() * b
        c0.append(c0_new_point)
        c1.append(c1_new_point)

    # Finally gather up and return the spline points (knots & control points)
    # as a List of Point2s
    res_controls = [(points[i], c0[i], c1[i]) for i in range(1, n - 1)]
    if close_loop:
        # Append the curve from the last input point
        # to the first input point,
        res_controls.append(
            (points[n - 1], c0[n - 1], c1[0]),
        )
        res_controls.append((points[1],))
    else:
        # Append the last input point only
        res_controls.append((points[n - 1],))
    return [Point3(*p.xyz) for p in flatten(res_controls)]


# Once the angles alpha and beta have been computed for each knot (which
# determine the direction from the knot to each of its neighboring handles),
# this function is used to compute the lengths of the vectors from the knot to
# those handles. Combining the length and angle together lets us solve for the
# handle positions.
#
# The exact choice of function is somewhat arbitrary. The aim is to return
# handle lengths that produce a Bézier curve which is a good approximation of a
# circular arc for points near the knot.
#
# Hobby and Knuth both proposed multiple candidate functions. This code uses
# the function from Jackowski formula 28, due to its simplicity. For other
# choices see Jackowski, section 5.
def rho(alpha: float, beta: float) -> float:
    """ "Velocity function" to compute the length of handles for a bezier spline."""
    c = 2.0 / 3
    return 2.0 / (1 + c * cos(beta) + (1 - c) * cos(alpha))


def thomas(A: Sequence, B: Sequence, C: Sequence, D: Sequence) -> Sequence:
    """Thomas algorithm: Solve a system of equations encoded in a tridiagonal matrix.

    A, B, and C are diagonals of the matrix. B is the main diagonal.
    D is the vector on the right-hand-side of the equation.

    Both B and D will have n elements. The arrays A and C will have
    length n as well, but each has one fewer element than B (the values
    A[0] and C[n-1] are undefined)."""

    # Note: n is defined here so that B[n] is valid, i.e. we are solving
    # a system of n+1 equations.
    n = len(B) - 1

    # Step 1: forward sweep to eliminate A[I] from each Equation
    # allocate lists for modified C and D coefficients, c_prime, d_prime
    c_prime = [C[0] / B[0]]
    d_prime = [D[0] / B[0]]

    for i in range(1, n + 1):
        denom = B[i] - c_prime[i - 1] * A[i]
        c_prime.append(C[i] / denom)
        d_prime.append((D[i] - d_prime[i - 1] * A[i]) / denom)

    # Step 2: back substitution to solve for x
    X = []
    X.append(d_prime[n])
    # Work back to front, solving for x[i] using x[-1]
    for i in range(n - 1, -1, -1):
        X.append(d_prime[i] - c_prime[i] * X[-1])
    return X[::-1]


# ===========
# = HELPERS =
# ===========
def vec3_perp_norm(v: Vector3, other: Vector3) -> Vector3:
    """Return the normalized cross product of vectors V and OTHER.

    Used to fetch a normalized vector perpendicular to
    both V and OTHER about which to rotate when forming smooth curves.
    """
    return v.cross(other).normalized()


def flatten(matr: Sequence[Sequence], keep_none_values: bool = True) -> List:
    return [item for row in matr for item in row if item or keep_none_values]
