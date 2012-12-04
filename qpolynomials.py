from basichypergeometric import BasicHypergeometricSeries

def little_q_Jacobi_polynomials(n, x, a, b, q):
    bhs = BasicHypergeometricSeries([q**(-n), a*b*q**(n+1)], [a*q], q, q*x)
    return bhs.evaluate()
