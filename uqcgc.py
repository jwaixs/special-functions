from basichypergeometric import qPochhammerSymbol
from qpolynomials import dual_q_Hahn_polynomials

def cgc(l1, l2, l, p, n):
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
