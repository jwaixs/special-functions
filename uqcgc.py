from basichypergeometric import qPochhammerSymbol
from qpolynomials import dual_q_Hahn_polynomials

def qpochhammer(a, qq, n):
    if n == 0:
        return 1
    return (1 - a*qq**n)*qpochhammer(a, qq, n-1)

def basichypergeometric32(a, b, c, d, e, qq, z, n=0):
    nom = product(map(lambda elm : qpochhammer(elm, qq, n), [a, b, c]))
    denom = product(map(lambda elm : qpochhammer(elm, qq, n), [d, e, qq]))
    
    if nom == 0:
        return 0
    return nom/denom*z**n + basichypergeometric32(a, b, c, d, e, qq, z, n+1)
    
def cgc(l1, l2, l, p, n):
    '''Clebsch-Gordan coefficients for case 2'''
    q = var('q')

    sign = (-1)**(l1+l2-l)

    if -2*l2 <= p and p < 0:
        # Case 1
        qpoch1 = qPochhammerSymbol(q**(-4*l2), q**2, 2*l2+p).evaluate()
        qpoch2 = qPochhammerSymbol(
            [q**(-4*l1), q**(-4*l1-4*l2-2), q**(-4*l2-2*p)],
            q**2, l1 + l2 - l).evaluate()
        qpoch3 = qPochhammerSymbol(
            q**(-4*l1 - 4*l2),
            q**2, 2*l2 + p).evaluate()
        qpoch4 = qPochhammerSymbol(
            [q**2, q**(-4*l1+2*p), q**(-4*l2)],
            q**2, l1 + l2 - l).evaluate()
        line1 = sqrt(qpoch1 * qpoch2 / (qpoch3 * qpoch4))

        line2 = sqrt((1 - q**(-2*l1-2*l2-2*l-2)) \
            / (1 - q**(-4*l1 - 4*l2 - 2))) \
            * q**(-2*l1*(2*l2+p) + (l1+l2-l)*(2*l1+2*l2-p) + l2*(2*p-4*l1) \
            - binomial(l1 + l2 - l, 2))

        qpoch5 = qPochhammerSymbol(
            [q**(-4*l2-2*p), q**(-4*l1)],
            q**2, n).evaluate()
        qpoch6 = qPochhammerSymbol(
            [q**2, q**(-2*p + 2)], q**2, n).evaluate()
        poly = basichypergeometric32(
            q**(-2*n), q**(-2*l1-2*l2+2*l), q**(-2*l1-2*l2-2*l-2),
            q**(-4*l1), q**(-4*l2-2*p), q**2, q**2)
        line3 = q**((2*l1 + 2*l2 + 1)*n) * sqrt(qpoch5 / qpoch6) * poly
    
    elif 0 <= p and p <= 2*l1 - 2*l2:
        # Case 2
        qpoch1 = qPochhammerSymbol(q**(-4*l2 - 2*p), q**2, 2*l2).evaluate()
        qpoch2 = qPochhammerSymbol(
            [q**(-4*l1 + 2*p), q**(-4*l1 - 4*l2 - 2), q**(-4*l2)], 
            q**2, 
            l1 + l2 - l).evaluate()
        qpoch3 = qPochhammerSymbol(
            q**(-4*l1 - 4*l2), 
            q**2, 
            2*l2).evaluate()
        qpoch4 = qPochhammerSymbol(
            [q**2, q**(-4*l1), q**(-4*l2 - 2*p)], 
            q**2,
            l1 + l2 - l).evaluate()
        line1 = sqrt(qpoch1 * qpoch2 / (qpoch3 * qpoch4))

        line2 = sqrt((1 - q**(-4*l - 2)) / (1 - q**(-4*l1 - 4*l2 - 2))) \
            * q**((l1+l2-l)*(2*l1+2*l2-p)+l2*(2*p-4*l1) \
            - binomial(l1 + l2 - l, 2))

        qpoch5 = qPochhammerSymbol(
            [q**(-4*l2), q**(-4*l1 + 2*p)],
            q**2, 
            n).evaluate()
        qpoch6 = qPochhammerSymbol(
            [q**2, q**(2*p + 2)], q**2, n).evaluate()
        poly = dual_q_Hahn_polynomials(
            n, 
            l1 + l2 - l, 
            q**(-4*l1 + 2*p - 2),
            q**(-4*l2 - 2*p - 2), 
            2*l2, 
            q**2)
        line3 = q**((2*l1 + 2*l2 + 1)*n) * sqrt(qpoch5 / qpoch6) * poly
    elif 2*l1 - 2*l2 < p and p <= 2*l1:
        # Case 3
        qpoch1 = qPochhammerSymbol(q**(-4*l1), q**2, 2*l1+p).evaluate()
        qpoch2 = qPochhammerSymbol(
            [q**(-4*l2), q**(-4*l1-4*l2-2), q**(-4*l1-2*p)],
            q**2, l1 + l2 - l).evaluate()
        qpoch3 = qPochhammerSymbol(
            q**(-4*l1 - 4*l2),
            q**2, 2*l1 + p).evaluate()
        qpoch4 = qPochhammerSymbol(
            [q**2, q**(-2*l1-4*l2-2*p), q**(-4*l1)],
            q**2, l1 + l2 - l).evaluate()
        line1 = sqrt(qpoch1 * qpoch2 / (qpoch3 * qpoch4))

        line2 = sqrt((1 - q**(-2*l1-2*l2-2*l-2)) \
            / (1 - q**(-4*l1 - 4*l2 - 2))) \
            * q**(-2*l2*(2*l1-p) + p*(l1+l2-l) - binomial(l1 + l2 - l, 2))

        qpoch5 = qPochhammerSymbol(
            [q**(-4*l2), q**(-4*l1+2*p)],
            q**2, n).evaluate()
        qpoch6 = qPochhammerSymbol(
            [q**2, q**(2*p + 2)], q**2, n).evaluate()
        poly = basichypergeometric32(
            q**(-2*n), q**(-2*l1-2*l2+2*l), q**(-2*l1-2*l2-2*l-2),
            q**(-4*l1+2*p), q**(-4*l2), q**2, q**2)
        line3 = q**((2*l1 + 2*l2 + 1)*n) * sqrt(qpoch5 / qpoch6) * poly
    else:
        return 0

    return line12 * line3

def cgc_test(l, n):
    q = var('q')
    return (-1)**n * q**(4*l - 2*n) \
        * sqrt((1 - q**4)*(1 - q**(2*n+2))/((1 - q**(4*l+2))*(1 - q**(4*l+4))))

def cgc_test2(l, n):
    q = var('q')
    return (-1)**n * q**(2*l - n) \
        * sqrt((1 - q**4)*(1 - q**(4*l - 2*n + 2))/((1 - q**(4*l+2))*(1 - q**(4*l+4))))
