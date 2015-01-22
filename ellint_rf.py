#http://arxiv.org/pdf/math/9409227v1.pdf

import math

_tol = 1e-40
_maxiter = 100

def ellint_rf(x, y, z):

    A0 = (x+y+z) / 3.0
    A = A0
    r = _tol
    Q = pow( 3 *r, -1 / 6.0) * max(abs(A -x), abs(A-y), abs(A-z))
    xm = x
    ym = y
    zm = z
    niter = 0

    while niter < _maxiter:

        sx = math.sqrt(xm)
        sy = math.sqrt(ym)
        sz = math.sqrt(zm)
        lm = sx * sy + sx * sz + sy * sz
        A = (A + lm) / 4.0
        xm = (xm + lm) / 4.0
        ym = (ym + lm) / 4.0
        zm = (zm + lm) / 4.0
        A4 = A * pow(4, niter)
        if Q < abs(A4):
            break
        niter += 1

    xx = (A0 - x) / A4
    yy = (A0 - y) / A4
    zz = -xx - yy

    E2 = xx * yy - zz**2
    E3 = xx * yy * zz

    val = 1 - E2 / 10.0 + E3 / 14.0 + E2**2 / 24.0 - 3.0 / 44.0 * E2 * E3

    return val / math.sqrt(A)

import cmath

def c_ellint_rf(x, y, z):

    A0 = (x+y+z) / 3.0
    A = A0
    r = _tol
    Q = pow( 3 *r, -1 / 6.0) * max(abs(A -x), abs(A-y), abs(A-z))
    xm = x
    ym = y
    zm = z
    maxiter = _maxiter
    niter = 0

    while niter < maxiter:

        sx = cmath.sqrt(xm)
        sy = cmath.sqrt(ym)
        sz = cmath.sqrt(zm)
        lm = sx * sy + sx * sz + sy * sz
        A = (A + lm) / 4.0
        xm = (xm + lm) / 4.0
        ym = (ym + lm) / 4.0
        zm = (zm + lm) / 4.0
        A4 = A * pow(4, niter)
        if Q < abs(A4):
            break
        niter += 1

    xx = (A0 - x) / A4
    yy = (A0 - y) / A4
    zz = -xx - yy

    E2 = xx * yy - zz**2
    E3 = xx * yy * zz

    val = 1 - E2 / 10.0 + E3 / 14.0 + E2**2 / 24.0 - 3.0 / 44.0 * E2 * E3

    return val / cmath.sqrt(A)

if __name__ == '__main__':

    print "Real"
    print "ellint_rf(2, 3, 4)"
    print ellint_rf(2, 3, 4)
    print "ellint_rf(1, 2, 0)"
    print ellint_rf(1, 2, 0)
    print "ellint_rf(0.5, 1, 0)"
    print ellint_rf(0.5, 1, 0)

    print "Complex"
    print "ellint_rf(2, 3, 4)"
    print c_ellint_rf(2, 3, 4)
    print "ellint_rf(1, 2, 0)"
    print c_ellint_rf(1, 2, 0)
    print "ellint_rf(i, -i, 0)"
    print c_ellint_rf(1j, -1j, 0)
    print "ellint_rf(0.5, 1, 0)"
    print c_ellint_rf(0.5, 1, 0)
    print "ellint_rf(-1+1j, 1j, 0)"
    print c_ellint_rf(-1+1j, 1j, 0)
