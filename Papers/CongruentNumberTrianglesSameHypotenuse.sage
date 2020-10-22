"""
CongruentNumberTrianglesSameHypotenuse.sage

written by
David Lowry-Duda
ICERM & Brown University
<david@lowryduda.com>

--------------------------------------------------------------------------------

# Introduction

The purpose of this software is to provide experimental evidence for the
following conjecture: there are no non-similar primitive integer right triangles
with the same squarefree part of the area and hypotenuse. "Non-similar" in this
context really means that we don't allow one to consider the two triangles (a,
b, c) and (b, a, c) as distinct.

In equations, we conjecture that there is no non-obvious solution to the pair of
equations

    a^2 + b^2 = A^2 + B^2
    ab u^2 = AB v^2

regarded over P^5(Q). The first equation comes from a^2 + b^2 = c^2 = C^2 = A^2
+ B^2, governing the right-ness of the two triangles (a, b, c) and (A, B, C).
The second equation states that the areas are equal up to multiplication by
squares. Considering these projectively eliminates the primitivity requirement.

This is not the optimal sort of equations for understanding this. 
Contact Lowry-Duda and Hassett for more on this.


## Background

This conjecture arose out of the work of Hulse, Kuan, Lowry-Duda, and Walker in
"A shifted sum for the congruent number problem" (in press at the Ramanujan
Journal, 2019-2020).


# Description of technique

To each right triangle (a, b, c) with rational coordinates and (squarefree) area
t, there is a corresponding point on the elliptic curve E_t: Y^2 = X^3 - t^2 X.
This correspondence is given by

    (a, b, c : a^2 + b^2 = c^2; ab/2 = t) <-----> E_t(Q): Y^2 = X^3 - t^2 X
                                  [a, b, c] --> [tb/(c-a), 2t^2/(c-a)]
      [(X^2 - t^2)/Y, 2tX/Y, (X^2 + t^2)/Y] <-- [X, Y]

Given a particular congruent number t, we find generators for E_t(Q) and
generate several points by adding them together. We check the corresponding
hypotenuses under this identification and see that none are the same.

In general finding generators is hard. We take advantage of the LMFDB and
perform this computation for all elliptic curves of this form within the LMFDB.

Additional are embedded within the program itself below.

--------------------------------------------------------------------------------

# LICENSE Information

Copyright David Lowry-Duda (2020)
Licensed under GPLv3.0.
"""
import itertools
import sys


def curve_pt_to_triangle(pt, t):
    """
    Given a (projective) point on the elliptic curve Y^2 = X^3 - t^2 X, return
    the corresponding right triangle.

    Inputs:
      pt: a projective point on E
      t:  the defining t from the elliptic curve

    Outputs:
      A tuple (a, b, c) of rationals
    """
    x, y = pt.xy()
    return (
        (x*x - t*t) / y,
        2 * t * x   / y,
        (x*x + t*t) / y
    )


def scale_to_primitive_triangle(triangle):
    """
    Given a triangle, scale it to become primitive and positive.

    Inputs:
      triangle: a tuple (a, b, c) of rationals.

    Outputs:
      A tuple (A, B, C) or integers
    """
    a, b, c = triangle
    if a < 0:
        a = -a
    if b < 0:
        b = -b
    if c < 0:
        c = -c
    g = gcd((a, b, c))
    return (a/g, b/g, c/g)


def curve_pt_to_primitive_triangle(pt, t):
    """
    Given a (projective) point on the elliptic curve Y^2 = X^3 - t^2 X, return
    the corresponding primitive right triangle (a, b, c) with a <= b < c.

    Inputs:
      pt: a projective point on E
      t:  the defining t from the elliptic crve

    Outputs:
      A tuple (a, b, c) of positive integers with a <= b < c.
    """
    a, b, c = scale_to_primitive_triangle(curve_pt_to_triangle(pt, t))
    if a > b:
        a, b = b, a
    return (a, b, c)


class SeenList:
    """
    A tracker for curves and hypotenuses thus far observed.
    """
    def __init__(self, t):
        """
        Establish the seen list for a particular curve t. Returns the method to
        check for inclusion in a particular t.

        Input:
          t: the defining t from the elliptic curve
        """
        self.t = t
        self.seen_triangles = set()
        self.seen_hypotenuses = set()

    def check_triangle(self, triangle):
        """
        Determines if two non-identical triangles with the same hypotenuse are
        given to this function. If they are, raise a ValueError.

        Input:
          triangle: a (primitive, oriented, integer) right triangle (a, b, c)
                    with a <= b < c
        """
        if triangle in self.seen_triangles:
            return
        self.seen_triangles.add(triangle)
        hypotenuse = triangle[2]
        if hypotenuse in self.seen_hypotenuses:
            print("Hypotenuse {} appears twice for t={}".format(hypotenuse, self.t))
            for oldtriangle in self.seen_triangles:
                if oldtriangle[2] == hypotenuse:
                    print(oldtriangle)
            print(triangle)
            print("EXCEPTION EXCEPTION EXCEPTION XXX")
            raise KeyError("Whoa!")
        self.seen_hypotenuses.add(hypotenuse)

    def check_from_pt(self, pt):
        """
        Given pt on elliptic curve corresponding to `t`, check if corresponding
        triangle is in seen list.

        First converts pt to triangle, and then calls `check_triangle`.
        """
        self.check_triangle(curve_pt_to_primitive_triangle(pt, self.t))


def find_Et_in_lmfdb(t):
    """
    Query LMFDB for data on the curve Y^2 = X^3 - t^2 X.

    Input:
      t: the defining t of the elliptic curve

    Output:
      a json dictionary of data from the LMFDB.
    """
    # Weierstrass model: Y^2 + a1 XY + a2 Y = X^3 + a2 X^2 + a4 X + a6
    # Here, we have Y^2 = X^3 - t^2 X, so the Weierstrass coefficients,
    # called the `ainvs`, are [a1, a2, a3, a4, a5] = [0, 0, 0, -t^2, 0].
    # This is what we lookup in the LMFDB.
    from lmfdb import db
    query = {'ainvs': [0, 0, 0, -t*t, 0]}
    curve = db.ec_curves.lucky(query=query)
    if curve is None:
        raise KeyError("No curve for {} in the LMFDB.".format(t))
    return curve


def gens_for_Et(t):
    """
    Query the LMFDB to find generators for the elliptic curve Y^2 = X^3 - t^2 X.
    If the rank is 0 (so that there are no generators), returns None. Otherwise
    returns a list of the generators.

    Input:
      t: the defining t of the elliptic curve

    Output:
      a list of tuples of generators (given projectively), or None if there are
      no generators.
    """
    curve = find_Et_in_lmfdb(t)
    return gens_for_ET_from_json(curve)


def gens_for_Et_from_json(jsondict):
    """
    Helper function for gens_for_Et.
    """
    if jsondict['rank'] == 0:
        return None
    return list(map(_genstr_to_list, jsondict['gens']))

def _genstr_to_list(genstring):
    """
    Convert a generator string as returned by the LMFDB to a projective point on
    an elliptic curve.

    Input:
      genstring: An LMFDB-style string representing a generator

    Output:
      a tuple containing the generator, given projectively
    """
    # The LMFDB will store generators as strings of the form u'(-3:9,1)'.
    # This
    #  - removes the parentheses (with the [1:-1] slicing) --> u'-3:9:1'
    #  - splits on ':' --> (u'-3', u'9', u'1')
    #  - maps each value to a rational (integer) and returns it as a tuple
    #    --> (-3, 9, 1)
    return tuple(map(QQ, genstring[1:-1].split(':')))


def is_squarefree(n):
    """
    Returns True if n is squarefree, False otherwise.

    This is naive and works only for numbers less than 10**6.
    """
    for i in range(2, 1001):
        if n % (i*i) == 0:
            return False
    return True


statistics = dict()
def run_lmfdb_trial_on_t(t, box_size=10, verbose=True):
    """
    Experimentally verify the conjecture for t. This is the driver class for an
    individual number t.

    This looks up the generators for E_t: Y^2 = X^3 - t^2 X, takes all linear
    combinations of them with coefficients less than `box_size`, and checks that
    all resulting triangles are unique.

    This also collects statistics for each congruent number, for fun.
    """
    if verbose:
        print("\nBeginning trials on t={}".format(t))

    # Only squarefree numbers should be checked.
    if not is_squarefree(t):
        raise ValueError("{} is not squarefree.".format(t))

    seenlist = SeenList(t)
    curve = find_Et_in_lmfdb(t)

    # Extra statistics
    stat = dict()
    stat['rank'] = curve['rank']

    if stat['rank'] == 0:
        if verbose:
            print("{} is not congruent.".format(t))
        statistics[t] = stat
        return

    if verbose:
        print("{} is congruent.".format(t))

    gens = gens_for_Et_from_json(curve)

    # Construct curve from short Weierstrass equation in sage. We will
    # be adding points on the curve, and sage will do this for us.
    E = EllipticCurve([-t*t, 0])

    if verbose:
        print("{} generators found. Expecting {}.".format(len(gens), stat['rank']))
        print("Beginning counting...")

    for array in itertools.product(
        range(-box_size, box_size + 1), repeat=stat['rank']
    ):
        # Multiplying by 0 gives the identity, and should be omitted
        if any(index == 0 for index in array):
            continue

        pt = E.point((0, 1, 0))   # The identity point in the group
        for mult, gen in zip(array, gens):
            pt += mult * E.point(gen)
            seenlist.check_from_pt(pt)

    num_triangles = (2*box_size) ** stat['rank']
    stat['num_triangles_checked'] = num_triangles

    largest_hypot = max(seenlist.seen_hypotenuses)
    stat['largest_hypot_log'] = n(log(largest_hypot))
    if verbose:
        print("{} triangles checked.".format(num_triangles))
        print("Largest hypotenuse: e^{}".format(stat['largest_hypot_log']))
        print("No duplicates found.\n")
    stat['method'] = 'lmfdb'
    statistics[t] = stat
    del seenlist


def run_raw_sage_trial_on_t(t, box_size=10, verbose=True):
    """
    Experimentally verify the conjecture for t. This is the driver class for an
    individual number t, but this class doesn't use the LMFDB. It computes parts
    explicitly.

    This looks up the generators for E_t: Y^2 = X^3 - t^2 X, takes all linear
    combinations of them with coefficients less than `box_size`, and checks that
    all resulting triangles are unique.

    This also collects statistics for each congruent number, for fun.
    """
    if verbose:
        print("\nBeginning raw sage trials on t={}".format(t))

    # Only squarefree numbers should be checked.
    if not is_squarefree(t):
        raise ValueError("{} is not squarefree.".format(t))

    seenlist = SeenList(t)
    curve = EllipticCurve([-t*t, 0])
    E = curve
    usepari = False

    # Extra statistics
    stat = dict()
    try:
        stat['rank'] = curve.rank()
    except RuntimeError:
        print("Trying to use pari to determine rank...")
        pE = pari(curve)
        try:
            stat['rank'] = pE.ellanalyticrank()[0]
            usepari = True
        except PariError:
            print("Sage/pari failed to determine rank.")
            stat['sagefail'] = -1
            statistics[t] = stat
            raise RuntimeError("Sage/pari failed to determine rank.")

    if stat['rank'] == 0:
        if verbose:
            print("{} is not congruent.".format(t))
        statistics[t] = stat
        return

    if verbose:
        print("{} is congruent.".format(t))

    try:
        if usepari:
            gen = pari(curve).ellheegner()
            gens = [curve.point(list(gen))]
        else:
            gens = curve.gens()
    except (RuntimeError, PariError):
        print("Sage/pari is unable to compute the generators of the curve.")
        stat['sagefail'] = -1
        statistics[t] = stat
        raise RunTimeError("Sage/pari failed to make generators.")

    if verbose:
        print("{} generators found. Expecting {}.".format(len(gens), stat['rank']))
        print("Beginning counting...")

    count = 0
    for array in itertools.product(
        range(-box_size, box_size + 1), repeat=stat['rank']
    ):
        array = list(array)
        # We can guarantee first index is nonnegative
        if array[0] < 0:
            continue

        # Multiplying by 0 gives the identity, and should be omitted
        if all(index == 0 for index in array):
            continue

        # Prevent going too large in multidimensional boxes
        if sum(map(abs, array)) > box_size / (stat['rank']**2):
            continue

        count += 1
        pt = E.point((0, 1, 0))   # The identity point in the group
        for mult, gen in zip(array, gens):
            pt += mult * E.point(gen)
        seenlist.check_from_pt(pt)

    num_triangles = count
    stat['num_triangles_checked'] = num_triangles

    min_hypot = min(seenlist.seen_hypotenuses)
    largest_hypot = max(seenlist.seen_hypotenuses)
    stat['largest_hypot_log'] = n(log(largest_hypot))
    stat['smallest_hypot'] = min_hypot
    if verbose:
        print("{} triangles checked.".format(num_triangles))
        print("Largest hypotenuse: e^{}".format(stat['largest_hypot_log']))
        print("No duplicates found.")
    stat['method'] = 'sage'
    statistics[t] = stat
    del seenlist


def driver(lower_bound=5, upper_bound=20, box_size=10, method='sage'):
    """
    Main entry point. Perform tests on each t up to `upper_bound`, using a box
    of side-length `box_size` for multiples.

    Input:
      upper_bound: the largest (integer) t to test
      box_size: the number of (integer) multiples of each generator to consider
      method: 'sage' or 'lmfdb' -- which method to use to find gens
    """
    import time
    # 5 is the smallest congruent number, so begin there
    for t in range(lower_bound, upper_bound+1):
        try:
            start = time.time()
            if method == 'sage':
                run_raw_sage_trial_on_t(t, box_size=box_size, verbose=True)
            elif method == 'lmfdb':
                run_lmfdb_trial_on_t(t, box_size=box_size, verbose=True)
            else:
                raise KeyError("Incorrect method requested.")
            duration = time.time() - start
            print("{} seconds elapsed.".format(duration))
        except ValueError:    # number is not squarefree
            print("{} is not squarefree. Skipping...".format(t))
        except RuntimeError:  # sage is unable to proceed
            print("\n\nFAILED ON {}\n\n".format(t))
            continue
        sys.stdout.flush()


if __name__ == "__main__":
    # In practice, adjust upper_bound and box_size
    driver(lower_bound=100, upper_bound=200, box_size=75, method='sage')
    print("\n\n\n")
    print("Statistics")
    print(statistics)
