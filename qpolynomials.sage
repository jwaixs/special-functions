#from basichypergeometric import *
load('basichypergeometric.sage')

def replace_z_x(polynomial):
    q, z, x = var('q z x')

    ret_poly = 0

    while polynomial.full_simplify() != 0:
        coef = polynomial.coefficients(z)
        coef = [ c for c in coef if c[0] != 0 ]

        constant = True
        for c, dg in coef:
            if dg != 0 and c != 0:
                constant = False        

        if constant:
            return ret_poly + polynomial.full_simplify()

        first, last = coef[0], coef[-1]

        if not (first[0] == last[0] and first[1] == -last[1]):
            return None
        
        ret_poly += first[0]*(2*x)**(last[1])
        polynomial -= first[0]*(z + z**(-1))**(last[1])

    return ret_poly

def askey_wilson(n, z, a, b, c, d, q):
    lc = qpoch([a*b, a*c, a*d], q, n) / a**n
    poly = bhs(
        [q**(-n), a*b*c*d*q**(n-1), a*z, a*z**(-1)],
        [a*b, a*c, a*d], q, q) 
    return lc*poly

class Askey_Wilson_polynomials():
    def __init__(self, n, x, a, b, c, d, q):
        self.n = n
        self.x = x
        self.q = q
        self.param = [a, b, c, d]

    def monic_three_term(self):
        n, q, x, a, b, c, d = [self.n, self.q, self.x] + self.param

        nom = (1 - a*b*q**n)*(1 - a*c*q**n)*(1 - a*d*q**n)*(1 - a*b*c*d*q**(n-1))
        denom = a*(1 - a*b*c*d*q**(2*n-1))*(1 - a*b*c*d*q**(2*n))
        An = nom/denom

        nom = a*(1 - q**n)*(1 - b*c*q**(n-1))\
            *(1 - b*d*q**(n-1))*(1 - c*d*q**(n-1))
        denom = (1 - a*b*c*d*q**(2*n-2))*(1 - a*b*c*d*q**(2*n-1))
        Cn = nom/denom

        rec1 = 1
        rec2 = 1/2*(a + a**(-1) - An - Cn)
        rec3 = 1/4*An.substitute({n : n-1})*Cn

        return rec1, rec2, rec3

#    def three_term(self):
#        n, q, x, a, b, c, d = [self.n, self.q, self.x] + self.param
#
#        An = (1 - a*b*q**n)*(1 - a*c*q**n)*(1 - a*d*q**n)*\
#            (1 - a*b*c*d*q**(n-1)) / \
#            (a*(1 - a*b*c*d*q**(2*n-1))*(1 - a*b*c*d*q**(2*n)))
#        Cn = a*(1 - q**n)*(1 - b*c*q**(n-1))*(1 - b*d*q**(n-1))*\
#            (1 - c*d*q**(n-1)) / \
#            ((1 - a*b*c*d*q**(2*n-2))*(1 - a*b*c*d*q**(2*n-1)))
#        Bn = (a + a**(-1) - An - Cn)
#
#        return An, Bn, Cn

    def three_term(self, n=None):
        if n == None:
            n = self.n
        q, x, a, b, c, d = [self.q, self.x] + self.param
        
        An = (1 - a*b*c*d*q**(n-1)) \
            / ((1 - a*b*c*d*q**(2*n-1))*(1 - a*b*c*d*q**(2*n)))
        Cn = (1 - q**n)*(1 - a*b*q**(n-1))*(1 - a*c*q**(n-1))* \
            (1 - a*d*q**(n-1))*(1 - b*c*q**(n-1))*(1 - b*d*q**(n-1))* \
            (1 - c*d*q**(n-1)) / \
            ((1 - a*b*c*d*q**(2*n-2))*(1 - a*b*c*d*q**(2*n-1)))
        s1, s2 = a + b + c + d, a**(-1) + b**(-1) + c**(-1) + d**(-1)
        Bn = q**(n-1)*((1 + a*b*c*d*q**(2*n-1))*(s1*q + s2*a*b*c*d) \
            - q**(n-1)*(1 + q)*a*b*c*d*(s1 + s2*q)) / \
            ((1 - a*b*c*d*q**(2*n-2))*(1 - a*b*c*d*q**(2*n)))

        return An, Bn, Cn

    def upto(self, endn):
        result = [0, 1]
        q, x, a, b, c, d = [self.q, self.x] + self.param

        for k in range(2, endn+1):
            an, bn, cn = self.three_term(n=k)
            result.append(
                (2*x - bn)*an**(-1)*result[k-1] - an**(-1)*cn*result[k-2]
            )

        return result
            
    
def little_q_Jacobi_polynomials(n, x, a, b, q):
    bhs = BasicHypergeometricSeries([q**(-n), a*b*q**(n+1)], [a*q], q, q*x)
    return bhs.evaluate()


def dual_q_Hahn_polynomials(n, x, c, d, N, q):
    if n != N:
        l1 = [q**(-n), q**(-x), c*d*q**(x + 1)]
        l2 = [c*q, q**(-N)]
        result = BasicHypergeometricSeries(l1, l2, q, q).evaluate()
    else:
        # For n = N we must understand the dual q-Hahn polynomials as
        # polynomials by continuity.
        result = 0
        for k in range(N + 1):
            qpoch1 = qPochhammerSymbol([q**(-x), c*d*q**(x + 1)], q, k)
            if qpoch1.evaluate() == 0:
                break
            qpoch2 = qPochhammerSymbol([c*q, q], q, k)
            result += qpoch1.evaluate() / qpoch2.evaluate() * q**k

    return result       

def rogers_polynomials(n, x, b, q, inx=True):
    z = var('z')

    ret = 0
    for k in range(n+1):
        ret += qpoch(b, q, k) * qpoch(b, q, n-k) \
            / (qpoch(q, q, k) * qpoch(q, q, n-k)) * z**(n - 2*k)

    if inx:
        return replace_z_x(ret)
    else:
        return ret

q_ultraspherical = rogers_polynomials

def continuous_q_jacobi(n, x, a, b, q):
    z = var('z')
    lc = qpoch(q**(a + 1), q, n) / qpoch(q, q, n)
    poly = bhs(
        [q**(-n), q**(n + a + b + 1), q**(1/2*a + 1/4)*z, 
            q**(1/2*a + 1/4)*z**(-1)],
        [q**(a + 1), -q**(1/2*(a + b + 1)), q**(1/2*(a + b + 2))], q, q)

    return lc*poly
    #return lc*replace_z_x(poly)
