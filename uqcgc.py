from basichypergeometric import qPochhammerSymbol
from qpolynomials import dual_q_Hahn_polynomials

def cgc(l1, l2, l, p, n):
    q = var('q')
    
    qpoch1 = qPochhammerSymbol(q**(-4*l2 - 2*p), q**2, 2*l2).evaluate()
    qpoch2 = qPochhammerSymbol([q**(-4*l1 + 2*p), q**(-4*l1 - 4*l2 - 2), \
        q**(-4*l2)], q**2, l1 + l2 - l).evaluate()
    qpoch3 = qPochhammerSymbol(q**(-4*l1 - 4*l2), q**2, 2*l2).evaluate()
    qpoch4 = qPochhammerSymbol([q**2, q**(-4*l1), q**(-4*l2 - 2*p)], q**2, \
        l1 + l2 - l).evaluate()

    sign = (-1)**(l1+l2-l)
    line12 = sqrt(sign * qpoch1 * qpoch2 / (qpoch3 * qpoch4) \
        * (1 - q**(-4*l - 2)) / (1 - q**(-4*l1 - 4*l2 - 2))) \
        * q**((l1+l2-l)*(2*l1+2*l2-p)+l2*(2*p-4*l1)-binomial(l1 + l2 - l, 2))

    #print mathematica(line12).Simplify(q>0).sage()

    qpoch5 = qPochhammerSymbol([q**(-4*l2), q**(-4*l1 + 2*p)], \
        q**2, n).evaluate()
    qpoch6 = qPochhammerSymbol([q**2, q**(2*p + 2)], q**2, n).evaluate()
    poly = dual_q_Hahn_polynomials(n, l1 + l2 - l, q**(-4*l1 + 2*p - 2), \
        q**(-4*l2 - 2*p - 2), 2*l2, q**2)
    line3 = q**((2*l1 + 2*l2 + 1)*n) * sqrt(qpoch5 / qpoch6) * poly

    #print mathematica(line3).Simplify(q>0).sage()

    return line12 * line3

def cgc_test(l, n):
    q = var('q')
    return (-1)**n * q**(4*l - 2*n) \
        * sqrt((1 - q**4)*(1 - q**(2*n+2))/((1 - q**(4*l+2))*(1 - q**(4*l+4))))

def cgc_test2(l, n):
    q = var('q')
    return (-1)**n * q**(2*l - n) \
        * sqrt((1 - q**4)*(1 - q**(4*l - 2*n + 2))/((1 - q**(4*l+2))*(1 - q**(4*l+4))))
