from basichypergeometric import BasicHypergeometricSeries, qPochhammerSymbol


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
       
