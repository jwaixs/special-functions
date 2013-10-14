from basichypergeometric import qpoch, bhs
from qpolynomials import dual_q_Hahn_polynomials

def product(ll):
    if len(ll) == 0:
        return 1
    return ll.pop()*product(ll)

def qpochhammer(a, qq, n):
    if type(a) == list:
        return product(map(lambda elm : qpochhammer(elm, qq, n)), a)
    if n == 0:
        return 1
    return (1 - a*qq**(n-1))*qpochhammer(a, qq, n-1)

def basichypergeometric32(la, lb, qq, z, n=0):
    if n == 0:
        return 1 + basichypergeometric32(la, lb, qq, z, 1)

    nom = qpochhammer(la, qq, n)
    denom = qpochhammer(lb + [qq], qq, n)
 
    if nom == 0:
        return 0
    return nom/denom*z**n + basichypergeometric32(a, b, c, d, e, qq, z, n+1)

def cgc(l1, l2, l, i, j, k):
    if i - j != k or i < -l1 or i > l1 or j < -l2 or j > l2:
        return 0

    if l2 > l1:
        return cgc(l2, l1, l, j, i, -k)

    p = l1 - l2 - k
    if -2*l2 <= p and p < 0:
        # Case 1
        n = l1 - i
    elif 0 <= p:
        # Case 2 and 3
        n = l2 - j

    q = var('q')

    sign = (-1)**(l1+l2-l)

    if -2*l2 <= p and p < 0:
        # Case 1
        qpoch1 = qpoch(q**(-4*l2), q**2, 2*l2+p)
        qpoch2 = qpoch( 
            [q**(-4*l1), q**(-4*l1-4*l2-2), q**(-4*l2-2*p)],
            q**2, l1 + l2 - l)
        qpoch3 = qpoch(q**(-4*l1 - 4*l2), q**2, 2*l2 + p)
        qpoch4 = qpoch(
            [q**2, q**(-4*l1+2*p), q**(-4*l2)],
            q**2, l1 + l2 - l)
        line1 = sqrt(qpoch1 * qpoch2 / (qpoch3 * qpoch4))

        line2 = sqrt((1 - q**(-2*l1-2*l2-2*l-2)) \
            / (1 - q**(-4*l1 - 4*l2 - 2))) \
            * q**(-2*l1*(2*l2+p) + (l1+l2-l)*(2*l1+2*l2+p) \
            - binomial(l1 + l2 - l, 2))

        qpoch5 = qpoch(
            [q**(-4*l2-2*p), q**(-4*l1)],
            q**2, n)
        qpoch6 = qpoch(
            [q**2, q**(-2*p + 2)], q**2, n)
        poly = bhs(
            [q**(-2*n), q**(-2*l1-2*l2+2*l), q**(-2*l1-2*l2-2*l-2)],
            [q**(-4*l1), q**(-4*l2-2*p)], q**2, q**2)
        line3 = q**((2*l1 + 2*l2 + 1)*n) * sqrt(qpoch5 / qpoch6) * poly
   
 
    elif 0 <= p and p <= 2*l1 - 2*l2:
        # Case 2 and case 3
        qpoch1 = qpoch(q**(-4*l2 - 2*p), q**2, 2*l2)
        qpoch2 = qpoch(
            [q**(-4*l1 + 2*p), q**(-4*l1 - 4*l2 - 2), q**(-4*l2)], 
            q**2, 
            l1 + l2 - l)
        qpoch3 = qpoch(q**(-4*l1 - 4*l2), q**2, 2*l2)
        qpoch4 = qpoch(
            [q**2, q**(-4*l1), q**(-4*l2 - 2*p)], 
            q**2,
            l1 + l2 - l)
        line1 = sqrt(qpoch1 * qpoch2 / (qpoch3 * qpoch4))


        line2 = sqrt((1 - q**(-4*l - 2)) / (1 - q**(-4*l1 - 4*l2 - 2))) \
            * q**((l1+l2-l)*(2*l1+2*l2-p)+l2*(2*p-4*l1) \
            - binomial(l1 + l2 - l, 2))

        qpoch5 = qpoch([q**(-4*l2), q**(-4*l1 + 2*p)], q**2, n)
        qpoch6 = qpoch([q**2, q**(2*p + 2)], q**2, n)
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
        qpoch1 = qpoch(q**(-4*l1), q**2, 2*l1-p)
        qpoch2 = qpoch(
            [q**(-4*l2), q**(-4*l1-4*l2-2), q**(-4*l1-2*p)],
            q**2, l1 + l2 - l)
        qpoch3 = qpoch(q**(-4*l1 - 4*l2), q**2, 2*l1 - p)
        qpoch4 = qpoch(
            [q**2, q**(-2*l1-4*l2-2*p), q**(-4*l1)],
            q**2, l1 + l2 - l)
        line1 = sqrt(qpoch1 * qpoch2 / (qpoch3 * qpoch4))

        line2 = sqrt((1 - q**(-2*l1-2*l2-2*l-2)) \
            / (1 - q**(-4*l1 - 4*l2 - 2))) \
            * q**(-2*l2*(2*l1-p) + p*(l1+l2-l) - binomial(l1 + l2 - l, 2))

        qpoch5 = qpoch([q**(-4*l2), q**(-4*l1+2*p)], q**2, n)
        qpoch6 = qpoch([q**2, q**(2*p + 2)], q**2, n)
        poly = bhs(
            [q**(-2*n), q**(-2*l1-2*l2+2*l), q**(-2*l1-2*l2-2*l-2)],
            [q**(-4*l1+2*p), q**(-4*l2)], q**2, q**2)
        line3 = q**((2*l1 + 2*l2 + 1)*n) * sqrt(qpoch5 / qpoch6) * poly
    else:
        return 0

    return line1 * line2 * line3

def phi_vec(l, l1, l2):
    z = var('z')
    ret = []

    for n in srange(-l, l+1):
        elm = 0
        for n1 in srange(-l1, l1+1):
            n2 = n1 - n
            elm += cgc(l1, l2, l, n1, n2, n)**2 * q**(-n1-n2) * z**(-n1-n2)
        ret.append(elm)

    return ret

def Phi0(l):
    lmatrix = []
    
    for k in srange(-l, l+1):
        lmatrix.append(phi_vec(l, (l+k)/2, (l-k)/2))

    return matrix(lmatrix).transpose()
