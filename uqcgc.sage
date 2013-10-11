from basichypergeometric import qPochhammerSymbol
from qpolynomials import dual_q_Hahn_polynomials

def product(ll):
    if len(ll) == 0:
        return 1
    return ll.pop()*product(ll)

def qpochhammer(a, qq, n):
    if n == 0:
        return 1
    return (1 - a*qq**(n-1))*qpochhammer(a, qq, n-1)

def basichypergeometric32(a, b, c, d, e, qq, z, n=0):
    if n == 0:
        return 1 + basichypergeometric32(a, b, c, d, e, qq, z, 1)

    nom = product(map(lambda elm : qpochhammer(elm, qq, n), [a, b, c]))
    denom = product(map(lambda elm : qpochhammer(elm, qq, n), [d, e, qq]))
 
    if nom == 0:
        return 0
    return nom/denom*z**n + basichypergeometric32(a, b, c, d, e, qq, z, n+1)

def cgc_or(l1, l2, l, i, j, k):
    if i - j != k:
        return 0

    p = l1 - l2 - k
    if -2*l2 <= p and p < 0:
        # Case 1
        n = l1 - i
    if 0 <= p:
        # Case 2 and 3
        n = l2 - j

    print p, n

    return cgc(l1, l2, l, p, n)
    
def cgc(l1, l2, l, p, n):
    '''Clebsch-Gordan coefficients for case 2'''
    q = var('q')

    sign = (-1)**(l1+l2-l)

    if -2*l2 <= p and p < 0:
        # Case 1
        print "Case 1"
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
            * q**(-2*l1*(2*l2+p) + (l1+l2-l)*(2*l1+2*l2+p) \
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
   
        print line1, line2, line3
 
    elif 0 <= p:
        # Case 2 and case 3
        print "Case 2"
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
        print "Case 3"
        qpoch1 = qPochhammerSymbol(q**(-4*l1), q**2, 2*l1-p).evaluate()
        qpoch2 = qPochhammerSymbol(
            [q**(-4*l2), q**(-4*l1-4*l2-2), q**(-4*l1-2*p)],
            q**2, l1 + l2 - l).evaluate()
        qpoch3 = qPochhammerSymbol(
            q**(-4*l1 - 4*l2),
            q**2, 2*l1 - p).evaluate()
        qpoch4 = qPochhammerSymbol(
            [q**2, q**(-2*l1-4*l2-2*p), q**(-4*l1)],
            q**2, l1 + l2 - l).evaluate()
        line1 = sqrt(qpoch1 * qpoch2 / (qpoch3 * qpoch4))

        line2 = sqrt((1 - q**(-2*l1-2*l2-2*l-2)) \
            / (1 - q**(-4*l1 - 4*l2 - 2))) \
            * q**(-2*l2*(2*l1-p) + p*(l1+l2-l) - binomial(l1 + l2 - l, 2))

        print line1, line2, l1, l2, l, sqrt((1 - q**(-2*l1-2*l2-2*l-2))/(1-q**(-4*l1-4*l2-2)))

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

    return line1 * line2 * line3
