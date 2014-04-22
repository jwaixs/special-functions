load('qpolynomials.sage')

q = var('q')
b = var('b')
x = var('x')

def L(x, l):
    ret = matrix(SR, 2*l+1, 2*l+1)
    
    for m in range(2*l + 1):
        for k in range(m + 1):
            ret[m, k] = qpoch(q**2, q**2, m) * qpoch(q**2, q**2, 2*k + 1) \
                / (qpoch(q**2, q**2, m+k+1) * qpoch(q**2, q**2, k)) \
                * rogers_polynomials(m - k, x, q**(2*k+2), q**2)

    return ret

def Linv(x, l):
    '''...A guess...'''
    ret = matrix(SR, 2*l+1, 2*l+1)

    for m in range(2*l + 1):
        for k in range(m + 1):
            ret[m, k] = rogers_polynomials(m - k, x, q**(-2*m), q**2) \
                * qpoch(q**2, q**2, m) * qpoch(q**2, q**2, m + k) \
                / (qpoch(q**2, q**2, 2*m) * qpoch(q**2, q**2, k)) \
                * q**(2*m*(m - k))

    return ret

def scalardiv(M1, M2, M, N):
    newm = matrix(SR, M, N)

    for i in range(M):
        for j in range(N):
            if M1[i,j] != 0:
                newm[i,j] = (M1[i,j] / M2[i,j]).full_simplify().factor()

    return newm

def qultratest(m, n):
    ret = 0
    x = var('x')
    q = var('q')

    for k in range(0, m-n+1):
        i = k + n
        ret += qpoch(q**2, q**2, 2*k + 2*n + 1) \
            * qpoch(q**2, q**2, k + 2*n) \
            / qpoch(q**2, q**2, m + n + k + 1) \
            / qpoch(q**2, q**2, 2*k + 2*n) \
            * q**(2*(k + n)*k) \
            * rogers_polynomials(m - n - k, x, q**(2*n + 2*k + 2), q**2) \
            * rogers_polynomials(k, x, q**(-2*n - 2*k), q**2)

    return ret.full_simplify()

def qultratest2(m, n):
    ret = 0
    x = var('x')
    q = var('q')

    for k in range(0, m - n + 1):
        ret += (1 - q**(4*k + 4*n + 2)) / qpoch(q**(4*n + 2*k + 2), q**2, m - n + 1) * rogers_polynomials(m - n - k, x, q**(2*n + 2*k + 2), q**2) * rogers_polynomials(k, x, q**(-2*n - 2*k), q**2) * q**(2*k*(k + n))

    return ret.full_simplify()

def frontterm(m, n, k):
    return (1 - q**(4*k + 4*n + 2))/(qpoch(q**(4*n + 2*k + 2), q**2, m - n + 1)) * q**(2*k*(k + n))

def L_general(x, l, b):
    ret = matrix(SR, 2*l+1, 2*l+1)
    
    for m in range(2*l + 1):
        for k in range(m + 1):
            ret[m, k] = qpoch(q**2, q**2, m) / qpoch(q**2, q**2, k) \
                * 1 / qpoch(q**(4*k) * b**2, q**2, m - k) \
                * q_ultraspherical(m - k, x, q**(2*k)*b, q**2)

    return ret

def Linv_general(x, l, b):
    '''...A guess...'''
    ret = matrix(SR, 2*l+1, 2*l+1)

    for m in range(2*l + 1):
        for k in range(m + 1):
            ret[m, k] = qpoch(q**2, q**2, m) \
                / qpoch(q**2, q**2, k) \
                / qpoch(q**(2*m + 2*k - 2) * b**2, q**2, m - k) \
                * q_ultraspherical(m - k, x, q**(2 - 2*m) * b**(-1), q**2) \
                * b**(m - k) * q**(2*(m - 1)*(m - k))

    return ret
