load('basichypergeometric.sage')

def askey_wilson(n, z, a, b, c, d, q):
    lc = qpoch([a*b, a*c, a*d], q, n) / a**n
    poly = bhs(
        [q**(-n), a*b*c*d*q**(n-1), a*z, a*z**(-1)],
        [a*b, a*c, a*d], q, q)

    return lc*poly

def qracah(n, x, a, b, c, d, q):
    pass

def continuous_dual_qhahn(n, z, a, b, c, q):
    lc = qpoch([a*b, a*c], q, n) / a**n
    poly = bhs([q**(-n), a*z, a*z**(-1)], [a*b, a*c], q, q)

    return lc*poly

def continuous_qhahn(n, phi, theta, a, b, c, q):
    pass

def big_qjacobi(n, x, a, b, c, q):
    return bhs([q**(-n), a*b*q**(n+1), x], [a*q, c*q], q, q)

def big_q_legendre(n, c, q):
    return bhs([q**(-n), q**(n+1), x], [q, c*q], q, q)

def qhahn(n, x, a, b, N, q):
    pass

def dual_qhahn(n, x, c, d, N, q):
    pass

def al_salam_chihara(n, z, a, b, q):
    lc = qpoch(a*b, q, n) / a**n
    poly = bhs([q**(-n), a*z, a*z**(-1)], [a*b, 0], q, q)

    return lc*bhs

def qmeixner_pollaczek(n, x, a, q):
    pass

def continuous_qjacobi(n, z, a, b, q):
    lc = qpoch(q, q, n) / qpoch(q**(a + 1), q, n)
    poly = bhs(
        [q**(-n), q**(n + a + b + 1), q**(1/2*a + 1/4)*z, q**(1/2*a + 1/4)*z**(-1)],
        [q**(a + 1), -q**(1/2*(a + b + 1)), -q**(1/2*(a + b + 2))], q, q)
    
    return lc*poly

def continuous_qultraspherical(n, z, b, q):
    ret = 0
    for k in range(n+1):
        ret += qpoch(b, q, k) * qpoch(b, q, n-k) \
            / (qpoch(q, q, k) * qpoch(q, q, n-k)) * z**(n - 2*k)

    return ret

rogers = continous_qultraspherical

def continous_qlegendre(n, z, q):
    return bhs(
        [q**(-n), q**(n+1), q**(1/4)*z, q**(1/4)*z**(-1)],
        [q, -q**(1/2), -q], q, q)

def big_qlaguerre(n, x, a, b, q):
    return bhs([q**(-n), 0, x], [a*q, b*q], q, q)

def little_qjacobi(n, x, a, b, q):
    return bhs([q**(-n), a*b*q**(n+1)], [a*q], q, q*x)

def little_qlegendre(n, x, q):
    return bhs([q**(-n), q**(n+1)], [q], q, q*x)

def qmeixner(n, x, b, c, q):
    return bhs([q**(-n), q**(-x)], [b*q], q, -q**(n+1)/c)

def quantum_qkrawtchouk(n, x, p, N, q):
    pass

def qkrawtchouk(n, x, p, N, q):
    pass

def affine_qkrawtchouk(n, x, p, N, q):
    pass

def dual_q_krawtchouk(n, x, c, N, q):
    pass

def continuous_big_qhermite(n, z, a, q):
    lc = a**(-n)
    poly = bhs([q**(-n), a*z, a*z**(-1)], [0, 0], q, q)

    return lc*poly

def continous_qlaguerre(n, z, a, q):
    lc = qpoch(q**(a + 1), q, n) / qpoch(q, q, n)
    poly = bhs(
        [q**(-n), q**(1/2*a + 1/4)*z, q**(1/2*a + 1/4)*z**(-1)],
        [q**(a + 1), 0], q, q)

    return lc*poly

def little_qlaguerre(n, a, q):
    return bhs([q**(-n), 0], [a*q], q, q*x)

wall = little_qlaguerre

def qlaguerre(n, x, a, q):
    lc = qpoch(q**(a + 1), q, n) / qpoch(q, q, n)
    poly = bhs([q**(-n)], [q**(a + 1)], q, -q**(n + a + 1)*x)

    return lc*poly

def qbessel(n, x, a, q):
    return bhs([q**(-n), -a*q**n], [0], q, q*x)

def qcharlier(n, x, a, q):
    return bhs([q**(-n), q**(-x)], [0], q, -q**(n+1)/a)

def al_salam_carlitz1(n, x, a, q):
    lc = (-a)**n * q**(binomial(n, 2))
    poly = bhs([q**(-n), x**(-1)], [0], a, q*x/a)

    return lc*poly

def al_salam_carlitz2(n, x, a, q):
    lc = (-a)**n * q**(-binomial(n, 2))
    poly = bhs([q**(-n), x], [], q, q**n/a)

    return lc*poly

def continous_q_hermite(n, z, q):
    return z*bhs([q**(-n), 0], [], q, q**n * z**(-2))

def stieltjes_wigert(n, x, q):
    lc = 1 / qpoch(q, q, n)
    poly = bhs([q**(-n)], [0], q, -q**(n+1)*x)

    return lc*poly

def discrete_q_hermite1(n, x, q):
    return q**binomial(n, 2) * bhs([q**(-n), x**(-1)], [0], q, -q*x)

def discrete_q_hermite2(n, x, q):
    lc = (I)**(-n) * q**(-binomial(n, 2))
    poly = bhs([q**(-n), I*x], [], q, -q**n)

    return lc*poly
