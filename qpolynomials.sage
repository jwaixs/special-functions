from basichypergeometric import BasicHypergeometricSeries, qPochhammerSymbol


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
