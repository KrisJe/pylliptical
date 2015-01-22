#http://arxiv.org/pdf/math/9409227v1.pdf
from __future__ import division

import math

_tol = 1e-60
_maxiter = 1000

from ellint_rc import ellint_rc

def ellint_rj(x, y, z, p):

    A0 = (x+y+z+2*p) / 5
    delt = (p-x)*(p-y)*(p-z)
    A = A0
    r = _tol
    Q = pow(r/4, -1 / 6) * max(abs(A -x), abs(A-y), abs(A-z), abs(A-p))
    xm = x
    ym = y
    zm = z
    pm = p
    niter = 0
    accum = 0

    while niter < _maxiter:

        sx = math.sqrt(xm)
        sy = math.sqrt(ym)
        sz = math.sqrt(zm)
        sp = math.sqrt(pm)
        lm = sx * sy + sx * sz + sy * sz
        A = (A + lm) / 4
        dm = (sp+sx)*(sp+sy)*(sp+sz)
        em = pow(4, -3 * niter) * delt / dm**2
        xm = (xm + lm) / 4
        ym = (ym + lm) / 4
        zm = (zm + lm) / 4
        pm = (pm + lm) / 4
        A4 = A * pow(4, niter)
        if Q < abs(A4):
            break
        accum += pow(4, -niter) * ellint_rc(1, 1 + em) / dm
        niter += 1

    xx = (A0 - x) / A4
    yy = (A0 - y) / A4
    zz = (A0 - z) / A4
    pp = (-xx - yy - zz) / 2

    E2 = xx * yy + xx * zz + yy * zz - 3 * pp**2
    E3 = xx * yy * zz + 3 * E2 * pp + 4 * pp**3
    E4 = (2 *xx * yy * zz + E2 * pp + 3 * pp**3) * pp
    E5 = xx * yy * zz * pp**3
    E3 = xx * yy * zz

    val = 1 - 3 / 14 * E2 + E3 / 6 + 9 * E2**2 / 88- 3 / 22 * E4
    val += -9 / 52 * E2 * E3 + 3 / 26 * E5

    return pow(4, -niter) * pow(A, -3 / 2) * val + 6 * accum

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
    print "ellint_rj(0, 1, 2, 3)"
    print ellint_rj(0, 1, 2, 3)
    print "ellint_rj(2, 3,4 ,5)"
    print ellint_rj(2, 3,4 ,5)
