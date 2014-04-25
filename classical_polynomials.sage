load('hypergeometric.sage')

def wilson(n, x, a, b, c, d):
    if n < 0:
        return 0

    leading = poch([a + b, a + c, a + d], n)
    poly = hs(
        [-n, n + a + b + c + d - 1, a + I*x, a - I*x],
        [a + b, a + c, a + d], 1)
    
    return leading*poly

def racah(n, x, a, b, c, d):
    maxdg = -min(filter(lambda elm : -elm in NN, [a + 1, b + d + 1, c + 1]))
    leading = poch([a + b, a + c, a + d], n)    

    poly = 0

    if n < 0 or n > maxdg:
        return poly

    if n != maxdg:
        poly = hs([-n, n + a + b + 1, -x, x + c + d + 1],
            [a + 1, b + d + 1, c + 1], 1)
    else:
        if maxdg == -a - 1:
            for k in range(n+1):
                poly += poch([n + a + b + 1, -x, x + c + d + 1], k) \
                    / poch([1, b + d + 1, c + 1], k)
        elif maxdg == -b - d - 1:
            for k in range(n+1):
                poly += poch([n + a + b + 1, -x, x + c + d + 1], k) \
                    / poch([1, a + 1, c + 1], k)
        else:
            for k in range(n+1):
                poly += poch([n + a + b + 1, -x, x + c + d + 1], k) \
                    / poch([1, a + 1, b + d + 1], k)

    return leading*poly

def continuous_dual_hahn(n, x, a, b, c):
    if n < 0:
        return 0

    leading = poch([a + b, a + c], n)
    poly = hs([-n, a + I*x, a - I*x], [a + b, a + c], 1)

    return leading*poly

def continuous_hahn(n, x, a, b, c, d):
    if n < 0:
        return 0

    leading = I**n * poch([a + c, a + d], n) / poch(1, n)
    poly = hs([-n, n + a + b + c + d - 1, a + I*x], [a + c, a + d], 1)

    return leading*poly

def hahn(n, x, a, b, N):
    if n > N or n < 0:
        return 0

    poly = 0
    if n != N:
        poly = hs([-n, n + a + b + 1, -x], [a + 1, -N], 1)
    else:
        for k in range(n+1):
            poly += poch([n + a + b + 1 - x], k) / poch([1, a + 1], k)

    return poly
 
def dual_hahn(n, x, c, d, N):
    if n > N or n < 0:
        return 0

    poly = 0
    if n != N:
        poly = hs([-n, -x, x + c + d + 1], [c + 1, -N], 1)
    else:
        for k in range(n+1):
            poly += poch([-x, x + c + d + 1], k) / poch([1, c + 1], k)

    return poly

def meixner_pollaczek(n, x, l, p):
    if n < 0:
        return 0
    
    leading = poch(2*l, n) / poch(1, n) * exp(I*n*p)
    poly = hs([-n, l + I*x], 2*l, 1 - exp(-2*I*p))

    return leading*poly

def jacobi(n, x, a, b):
    if n < 0:
        return 0

    leading = poch(a + 1, n) / poch(1, n)
    poly = hs([-n, n + a + b + 1], [a + 1], (1 - x)/2)

    return leading*poly

# Special cases of Jacobi

def gegenbauer(n, x, l):
    if n < 0:
        return 0

    leading = poch(2*l, n) / poch(1, n)
    poly = hs([-n, n + 2*l], [l + 1/2], (1 - x)/2)

    return leading*poly

ultraspherical = gegenbauer

def chebyshev1(n, x):
    if n < 0:
        return 0

    return hs([-n, n], [1/2], (1 - x)/2)

def chebyshev2(n, x):
    if n < 0:
        return 0

    return (n + 1) * hs([-n, n+2], [3/2], (1 - x)/2)

def legendre(n, x):
    if n < 0:
        return 0

    return hs([-n, n+1], [1], (1 - x)/2)

spherical = legendre

def pseudo_jacobi(n, x, v, N):
    if n < 0 or n > N:
        return 0

    leading = (-2*I)**n * poch(-N + I*v, n) / poch(n - 2*N - 1, n)
    poly = hs([-n, n - 2*N - 1], [-N + I*v], (1 - I*x)/2)

    return leading*poly

def meixner(n, x, b, c):
    if n < 0:
        return 0

    return hs([-n, -x], b, 1 - 1/c)

def krawtchouk(n, x, p, N):
    if n < 0 or n > N:
        return 0
    
    poly = 0
    if n != N:
        poly = hs([-n, -x], [-N], 1/p)
    else:
        for k in range(n + 1):
            poly += poch([-x], k) / poch([1], k) * (1 / p)**k

    return poly

def laguerre(n, x, a):
    if n < 0:
        return 0

    leading = poch(a + 1, n) / poch(1, n)
    poly = hs([-n], [a + 1], x)

    return leading*poly

def bessel(n, x, a):
    if n < 0:
        return 0

    return hs([-n, n + a + 1], [], -x/2)

def charlier(n, x, a):
    if n < 0:
        return 0

    return hs([-n, -x], -1/a)

def hermite(n, x):
    if n < 0:
        return 0

    leading = (2*x)**n
    poly = hs([-n/2, -(n-1)/2], [], -1/x**2)

    return leading*poly
