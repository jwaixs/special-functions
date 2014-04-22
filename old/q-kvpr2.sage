load('qpolynomials.sage')
load('basichypergeometric.sage')
load('matrixfuncs.sage')

x = var('x')
z = var('z')
q = var('q')


## Replacing functions

def replace_z_x(polynomial):
    q, z, x = var('q z x')

    ret_poly = 0

    while polynomial.full_simplify() != 0:
        coef = polynomial.coefficients(z)
        coef = [ c for c in coef if c[0] != 0 ]

        constant = True
        for c, dg in coef:
            if dg != 0 and c != 0:
                constant = False        

        if constant:
            return ret_poly + polynomial.full_simplify()

        first, last = coef[0], coef[-1]

        if not (first[0] == last[0] and first[1] == -last[1]):
            return None
        
        ret_poly += first[0]*(2*x)**(last[1])
        polynomial -= first[0]*(z + z**(-1))**(last[1])

    return ret_poly

def replace_x_z(polynomial):
    return polynomial.substitute({x : (z + z**(-1))/2})

## LDU decomposition of W

def c(k, l):
    return (1 - q**(4*k + 2)) * qpoch(q**2, q**2, 2*l + k + 1) \
        * qpoch(q**2, q**2, 2*l - k) * qpoch(q**2, q**2, k)**4 \
        / (qpoch(q**2, q**2, 2*k + 1)**2 * qpoch(q**2, q**2, 2*l)**2)

def T(x, l):
    matT = matrix(SR, 2*l+1, 2*l+1)

    for k in range(2*l + 1):
        kweight = replace_z_x(qpoch([z**2, z**(-2)], q**2, k+1))
        matT[k,k] = q**(-2*l) * c(k,l) / 4  * kweight / (1 - x**2)

    return matT

def L(x, l):
    matL = matrix(SR, 2*l+1, 2*l+1)

    for m in range(2*l+1):
        for k in range(m+1):
            matL[m, k] = q**(m - k) * qpoch(q**2, q**2, m) * qpoch(q**2, q**2, 2*k + 1) \
                / (qpoch(q**2, q**2, m+k+1) * qpoch(q**2, q**2, k)) \
                * q_ultraspherical(m - k, x, q**(2*k + 2), q**2)

    return matL

def W(x, l):
    matL = L(x, l)
    matT = T(x, l)

    return matL * matT * transpose(matL)


## Orthogonality of the polynomials

def innerprod(P, Q, W):
    Mret = transpose(P)*W*Q

    return Mret.apply_map(lambda elm : integral(sqrt(1 - x**2)*elm, x, -1, 1))


## Construction of the monic orthogonal polynomials

def monic_polynomials(l, N):
    Pret = [identity_matrix(SR, 2*l+1)]
    weight = SimplifyMatrix(W(x, l))

    for n in range(1, N+1):
        newP = x**n * identity_matrix(SR, 2*l+1)
        for m in range(n):
            Am = -innerprod(Pret[m], Pret[m], weight)**(-1) \
                * innerprod(Pret[m], x**n*identity_matrix(SR, 2*l+1), weight)
            newP += Am * Pret[m]
        Pret.append(newP)

    return Pret

def R_polynomials(l, mpoly):
    Lmat = transpose(L(x, l))
    return map(lambda poly : SimplifyMatrix(Lmat*poly), mpoly)


## q-difference equations for the monic orthogonal polynomials

def M1(l):
    Mret = matrix(SR, 2*l+1, 2*l+1)
    
    for i in range(2*l):
        Mret[i,i+1] = -q**(2 - i) * (1 - q**(2*i + 2)) / (1 - q**2)**2 \
            * z / (1 - z**2)

    for i in range(2*l + 1):
        Mret[i,i] = q**(2 - i) / (1 - q**2)**2 \
            * (1 - q**(2*i + 2)*z**2) / (1 - z**2)

    return Mret

def N1(l):
    Nret = matrix(SR, 2*l+1, 2*l+1)

    for i in range(2*l):
        Nret[i,i+1] = q**(2 - i) * (1 - q**(2*i + 2)) / (1 - q**2)**2 \
            * z / (1 - z**2)
    
    for i in range(2*l + 1):
        Nret[i,i] = q**(2 - i) / (1 - q**2)**2 \
            * (1 - q**(2*i + 2)*z**(-2)) / (1 - z**(-2))

    return Nret

def M2(l):
    return antimatrix(SR, 2*l+1) * M1(l) * antimatrix(SR, 2*l+1)

def N2(l):
    return antimatrix(SR, 2*l+1) * N1(l) * antimatrix(SR, 2*l+1)

def D1(l, poly):
    p1 = M1(l)*poly.apply_map(replace_x_z).substitute({z : q*z})
    p2 = N1(l)*poly.apply_map(replace_x_z).substitute({z : q**(-1)*z})

    return SimplifyMatrix(p1 + p2).apply_map(replace_z_x)

def D2(l, poly):
    p1 = M2(l)*poly.apply_map(replace_x_z).substitute({z : q*z})
    p2 = N2(l)*poly.apply_map(replace_x_z).substitute({z : q**(-1)*z})

    return (p1 + p2).apply_map(replace_z_x)

## Guesses

def Rguess(n, l):
    Rret = matrix(SR, 2*l+1, 2*l+1)

    for i in range(2*l+1):
        for j in range(2*l+1):
            Rret[i,j] = q_ultraspherical(n+j-i, x, q**(2*i + 2), q**2)

    return Rret
