from basichypergeometric import qPochhammerSymbol
from qpolynomials import dual_q_Hahn_polynomials

def coef(l1, l2, l, p, n, algorithm=None):
    q = var('q')

    if p >= 0:
        result1 = q**(2*n*(l1+l2) + n)

        qpoch1 = qPochhammerSymbol([q**(-4*l2), q**(-4*l1 + 2*p)], q**2, n)
        qpoch2 = qPochhammerSymbol([q**2, q**(2*p + 2)], q**2, n)
        result2 = sqrt(qpoch1.evaluate() / qpoch2.evaluate())
        
        if p <= 2*l1 - 2*l2:
            result3 = dual_q_Hahn_polynomials( \
                n, l1 + l2 - l, q**(-4*l1 + 2*p - 2), \
                q**(-4*l2 - 2*p - 2), 2*l2, q**2 \
            )
        else:
            result3 = dual_q_Hahn_polynomials( \
                n, l1 + l2 - n, q**(-4*l2 - 2), \
                q**(-4*l1 - 2), 2**l1 - 2*p, q**2
            )

    else:
        return coef(l2, l1, l, -p, n)

    if algorithm == 'mathematica':
        return mathematica(result1*result2*result3).FullSimplify().sage() 
    return result1*result2*result3

def normalized_cgc_p(l1, l2, l, p):
    q = var('q')
    cgc = []

    for n in range(2*l2 + 1):
        cg = coef(l1, l2, l, p, n, algorithm='mathematica')
        cgc.append((l1, l2, l, l1-p-n, l2-n, l1-l2-p), cg)

    cgc_sum = sum([ elm[1]**2 for elm in cgc ])
    cgc2 = [ (elm[0], mathematica(sqrt(1/cgc_sum)).FullSimplify().sage()) \
        for elm in cgc ]

    return cgc2
    
    
