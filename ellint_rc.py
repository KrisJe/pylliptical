# http://arxiv.org/pdf/math/9409227v1.pdf

from __future__ import division
import math

_tol = 1e-40
_maxiter = 100

def ellint_rc(x, y):

    A0 = (x +2 * y) / 3
    A = A0
    A4 = A
    r = _tol
    Q = pow( 3 *r, -1 / 8) * abs(A0 - x)

    xm = x
    ym = y
    niter = 1

    while niter < _maxiter:

        sx = math.sqrt(xm)
        sy = math.sqrt(ym)
        lm = 2 * sx * sy + ym
        A = (A + lm) / 4
        xm = (xm + lm) / 4
        ym = (ym + lm) / 4
        A4 = A4 * 4
        if Q < abs(A4):
            break
        niter += 1

    s = (y - A0) / A4

    val = 1 - 3 / 10 * s**2 + 1 / 7 * s**3 + 3 / 8 * s**4 + 9 / 22 * s**5 + 159 / 208 * s**6 + 9 / 8 * s**7

    return val / math.sqrt(A)

import cmath

def c_ellint_rc(x, y):

    A0 = (x + 2 * y) / 3
    A = A0
    A4 = A
    r = _tol
    Q = pow(3 *r, -1 / 8) * abs(A0 - x)

    xm = x
    ym = y
    niter = 1

    while niter < _maxiter:

        sx = cmath.sqrt(xm)
        sy = cmath.sqrt(ym)
        lm = 2 * sx * sy + ym
        A = (A + lm) / 4
        xm = (xm + lm) / 4
        ym = (ym + lm) / 4
        A4 = A4 * 4
        if Q < abs(A4):
            break
        niter += 1

    s = (y - A0) / A4

    val = 1 - 3 / 10 * s**2 + 1 / 7 * s**3 + 3 / 8 * s**4 + 9 / 22 * s**5 + 159 / 208 * s**6 + 9 / 8 * s**7

    return val / cmath.sqrt(A)

if __name__ == '__main__':

    print "Real"
    print "ellint_rc(0, 0.25)"
    print ellint_rc(0, 0.25)
    print "ellint_rc(9 / 4, 2)"
    print ellint_rc(9 / 4, 2)

    print "Complex"
    print "c_ellint_rc(0, 0.25)"
    print c_ellint_rc(0, 0.25)
    print "c_ellint_rc(9 / 4, 2)"
    print c_ellint_rc(9 / 4, 2)
    print "c_ellint_rc(1 / 4, -2)"
    print c_ellint_rc(1 / 4, -2)
    print "c_ellint_rc(0, i)"
    print c_ellint_rc(0, 1j)

