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
        cgc.append(((l1, l2, l, l1-p-n, l2-n, l1-l2-p), cg))

    cgc_sum = sum([ elm[1]**2 for elm in cgc ])
    cgc2 = [ (elm[0], mathematica(sqrt(1/cgc_sum)*elm[1]).Simplify([0 < q, q < 1]).sage()) \
        for elm in cgc ]

    return cgc2
    
def normalized_cgc(l1, l2, l):
    q = var('q')
    cgc = {}

    for p in range(2*l + 1):
        cgc2 = normalized_cgc_p(l1, l2, l, p)
        for key, elm in cgc2:
            cgc[key] = elm

    return cgc 

def inverse_cgc(cgc):
    import copy
    cgc2 = copy.deepcopy(cgc)
    for key in cgc:
        cgc2[(key[1], key[0], key[2], key[4], key[3], -key[5])] = cgc[key]
    return cgc2

def print_normalized_cgc(l1, l2, l):
    cgc = normalized_cgc(l1, l2, l)
    for key in cgc:
        print key, "\t=\t", cgc[key]

#cgc1 = inverse_cgc(normalized_cgc(1, 1/2, 1/2))
#cgc2 = normalized_cgc(1/2, 1/2, 0)
#cgc3 = inverse_cgc(normalized_cgc(1/2, 1/2, 1))
#
#import itertools
#
#cgc_sum = 0
#
#for j1, i1, i2, n2 in itertools.product([-1/2, 1/2], repeat=4):
#    for j2, n1 in itertools.product([-1, 0, 1], repeat=2):
#        if j1 - j2 == 1/2 and i1 - i2 == 0 and j1 - i1 == n1 and j2 - i2 == n2 and n1 - n2 == 1/2:
#            c1, c2, c3, c4, c5, cgc_sum = cgc1[(1/2, 1, 1/2, j1, j2, 1/2)], cgc2[(1/2, 1/2, 0, i1, i2, 0)], cgc3[(1/2, 1/2, 1, j1, i1, n1)], cgc1[(1, 1/2, 1/2, j2, i2, n2)], cgc1[(1, 1/2, 1/2, n1, n2, 1/2)], cgc_sum + c1*c2*c3*c4*c5
#
#print mathematica(cgc_sum).Simplify().sage()

def new_normalized_cgc(l1, l2, l, p, n):
    q = var('q')

    sign = (-1)**(l1+l2-l)
    
    qpoch1 = qPochhammerSymbol(q**(-4*l2 - 2*p), q**2, 2*l2).evaluate()
    qpoch2 = qPochhammerSymbol([q**(-4*l1 + 2*p), q**(-4*l1 - 4*l2 - 2), \
        q**(-4*l2)], q**2, l1 + l2 - l).evaluate()
    qpoch3 = qPochhammerSymbol(q**(-4*l1 - 4*l2), q**2, 2*l2).evaluate()
    qpoch4 = qPochhammerSymbol([q**2, q**(-4*l1), q**(-4*l2 - 2*p)], q**2, \
        l1 + l2 - l).evaluate()
    line1 = sqrt(qpoch1 * qpoch2 / (qpoch3 * qpoch4))

    line2 = sqrt((1 - q**(-4*l - 2)) / (1 - q**(-4*l1 - 4*l2 - 2))) \
        * q**((l1+l2-l)*(2*l1+2*l2-p)+l2*(2*p-4*l1)-binomial(l1 + l2 - l, 2))

    qpoch5 = qPochhammerSymbol([q**(-4*l2), q**(-4*l1 + 2*p)], \
        q**2, n).evaluate()
    qpoch6 = qPochhammerSymbol([q**2, q**(2*p + 2)], q**2, n).evaluate()
    poly = dual_q_Hahn_polynomials(n, l1 + l2 - l, q**(-4*l1 + 2*p - 2), \
        q**(-4*l2 - 2*p - 2), 2*l2, q**2)
    line3 = q**((2*l1 + 2*l2 + 1)*n) * sqrt(qpoch5 / qpoch6) * poly

    return line1 * line2 * line3
