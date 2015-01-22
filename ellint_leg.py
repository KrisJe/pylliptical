from __future__ import division
import math
from ellint_rf import ellint_rf
from ellint_rj import ellint_rj

def ellint_pi(p,k,n):
    c = 1 / math.sin(p)**2
    x = c - 1
    y = c - k**2
    z = c
    p = c + n
    return ellint_rf(x, y, z) - n / 3 * ellint_rj(x, y, z, p)

if __name__ == '__main__':
    n = 1.3
    phi  = 1.1
    k = math.sqrt(0.2)
    print "ellint_pi(phi, k, n)"
    print ellint_pi(phi, k, n)
