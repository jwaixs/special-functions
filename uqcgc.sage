from basichypergeometric import qpoch, bhs
from qpolynomials import dual_q_Hahn_polynomials
import sys

def cgc(l1, l2, l, i, j, k):
    if i - j != k or i < -l1 or i > l1 or j < -l2 or j > l2 or k < -l or k > l:
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
        #print 'case 1'
        qpoch1 = qpoch(q**(-4*l2), q**2, 2*l2+p)
        qpoch2 = qpoch( 
            [q**(-4*l1), q**(-4*l1-4*l2-2), q**(-4*l2-2*p)],
            q**2, l1 + l2 - l)
        qpoch3 = qpoch(q**(-4*l1 - 4*l2), q**2, 2*l2 + p)
        qpoch4 = qpoch(
            [q**2, q**(-4*l1+2*p), q**(-4*l2)],
            q**2, l1 + l2 - l)
        line1 = sqrt((-1)**(l1 + l2 - l) * qpoch1 * qpoch2 / (qpoch3 * qpoch4))

        line2 = sqrt((1 - q**(-4*l-2)) \
            / (1 - q**(-4*l1 - 4*l2 - 2))) \
            * q**(-2*l1*(2*l2+p) + (l1+l2-l)*(2*l1+2*l2+p) \
            - binomial(l1 + l2 - l, 2))

        qpoch5 = qpoch(
            [q**(-4*l2-2*p), q**(-4*l1)],
            q**2, n)
        qpoch6 = qpoch(
            [q**2, q**(-2*p + 2)], q**2, n)
        #poly = bhs(
        #    [q**(-2*n), q**(-2*l1-2*l2+2*l), q**(-2*l1-2*l2-2*l-2)],
        #    [q**(-4*l1), q**(-4*l2-2*p)], q**2, q**2)
        poly = dual_q_Hahn_polynomials(
            n, l1+l2-l, q**(-4*l1-2), q**(-4*l2-2), 2*l2+p, q**2
        )
        line3 = q**((2*l1 + 2*l2 + 1)*n) * sqrt(qpoch5 / qpoch6) * poly
   
 
    elif 0 <= p and p <= 2*l1 - 2*l2:
        # Case 2
        #print 'case 2'
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
        line1 = sqrt((-1)**(l1 + l2 - l) * qpoch1 * qpoch2 / (qpoch3 * qpoch4))


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
        #print 'case 3'
        qpoch1 = qpoch(q**(-4*l1), q**2, 2*l1-p)
        qpoch2 = qpoch(
            [q**(-4*l2), q**(-4*l1-4*l2-2), q**(-4*l1+2*p)],
            q**2, l1 + l2 - l)
        qpoch3 = qpoch(q**(-4*l1 - 4*l2), q**2, 2*l1 - p)
        qpoch4 = qpoch(
            [q**2, q**(-4*l2-2*p), q**(-4*l1)],
            q**2, l1 + l2 - l)
        line1 = sqrt((-1)**(l1 + l2 - l) * qpoch1 * qpoch2 / (qpoch3 * qpoch4))

        line2 = sqrt((1 - q**(-4*l-2)) / (1 - q**(-4*l1 - 4*l2 - 2))) \
            * q**(-2*l2*(2*l1-p) + (2*l1+2*l2-p)*(l1+l2-l) - binomial(l1 + l2 - l, 2))

        qpoch5 = qpoch([q**(-4*l2), q**(-4*l1+2*p)], q**2, n)
        qpoch6 = qpoch([q**2, q**(2*p + 2)], q**2, n)
        poly = dual_q_Hahn_polynomials(
            n, l1+l2-l, q**(-4*l2-2), q**(-4*l1-2), 2*l1-p, q**2
        )
        line3 = q**((2*l1 + 2*l2 + 1)*n) * sqrt(qpoch5 / qpoch6) * poly
        
    else:
        return 0

    return line1 * line2 * line3

def test_cgc(l1, l2):
    from itertools import product

    # test 1
    for n1, m1 in product(srange(-l1, l1+1), repeat=2):
        for n2, m2 in product(srange(-l2, l2+1), repeat=2):
            s = 0
            for l in srange(abs(l1-l2), l1+l2+1):
                for n in srange(-l, l+1):
                    s += cgc(l1, l2, l, n1, n2, n) \
                        * cgc(l1, l2, l, m1, m2, n)
            if bool(s != kronecker_delta(n1, m1)*kronecker_delta(n2, m2)):
                for l in srange(abs(l1-l2), l1+l2+1):
                    for n in srange(-l, l+1):
                        ret = cgc(l1, l2, l, n1, n2, n) * cgc(l1, l2, l, m1, m2, n)
                        if ret != 0:
                            print l1, l2, l, n1, n2, n
                            print l1, l2, l, m1, m2, n
                            print cgc(l1, l2, l, n1, n2, n) * cgc(l1, l2, l, m1, m2, n)
                return False, 1, n1, n2, m1, m2, s
            sys.stdout.write('.')
            sys.stdout.flush()

    # test 2
    for l, ll in product(srange(abs(l1-l2), l1+l2+1), repeat=2):
        for n, m in product(srange(-l, l+1), srange(-ll, ll+1)):
            s = 0
            for n1, n2 in product(srange(-l1, l1+1), srange(-l2, l2+1)):
                s += cgc(l1, l2, l, n1, n2, n)*cgc(l1, l2, ll, n1, n2, m)
            if bool(s != kronecker_delta(l, ll)*kronecker_delta(n, m)):
                return False, 2, l, ll, n, m, s
            sys.stdout.write('.')
            sys.stdout.flush()

    return True

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

def Phi0(l, simplify=False):
    lmatrix = []
    
    for k in srange(-l, l+1):
        lmatrix.append(phi_vec(l, (l+k)/2, (l-k)/2))

    if simplify:
        return matrix(lmatrix).transpose().apply_map(
            lambda elm : maple(elm).simplify().sage())
    return matrix(lmatrix).transpose()

def phi_star_vec(l, l1, l2):
    z = var('z')
    ret = []

    for n in srange(-l, l+1):
        elm = 0
        for n1 in srange(-l1, l1+1):
            n2 = n1 - n
            elm += cgc(l1, l2, l, n1, n2, n)**2 * q**(n1+n2) * z**(n1+n2)
        ret.append(elm)

    return ret

def Phi0_star(l, simplify=False):
    lmatrix = []
    
    for k in srange(-l, l+1):
        lmatrix.append(phi_star_vec(l, (l+k)/2, (l-k)/2))

    if simplify:
        return matrix(lmatrix).transpose().apply_map(
            lambda elm : maple(elm).simplify().sage())
    return matrix(lmatrix).transpose()

# Explicit computations
def f(z):
    q = var('q')
    return (1 - q**2*z**2) / ((1 - z**2) * (q - q**(-1))**2)

def cc(l, n):
    '''Squared!'''
    q = var('q')

    return (q**(n - l) - q**(l - n))*(q**(-l - n - 1) - q**(l + n + 1)) \
        / (q**(-1) - q)**2

def bb(l, n):
    '''Squared!'''
    q = var('q')

    return (q**(-l + n - 1) - q**(l - n + 1))*(q**(-l - n) - q**(l + n)) \
        / (q**(-1) - q)**2

def matN(l, i):
    q, z = var('q z')
    
    NN = matrix(SR, 2*l+1, 2*l+1)

    if i == 1:
        for n in srange(-l, l+1):
            NN[n+l, n+l] = q**n * f(z**(-1)) \
                + q*z*q**(-n) * (cc(l,n) + bb(l,n)) \
		/ ((1 - z**2)*(z - z**(-1)))
        for n in srange(-l, l):
            NN[n+l, n+l+1] = -q**2*q**(-n-1) * cc(l, n) \
                / ((1 - z**2)*(z - z**(-1)))
        for n in srange(-l+1, l+1):
            NN[n+l, n+l-1] = -q * z**2 * q**(-n) * bb(l, n) \
                / ((1 - z**2)*(z - z**(-1)))
    elif i == 2:
        for n in srange(-l, l+1):
            NN[n+l, n+l] = q**(-n) * f(z**(-1)) \
                + z*q**(1+n) * (cc(l, n) + bb(l, n)) \
		/ ((1 - z**2)*(z - z**(-1)))
        for n in srange(-l, l-1):
            NN[n+l, n+l+1] = -q**(1 + n) * z**2 * cc(l, n) \
                / ((1 - z**2)*(z - z**(-1)))
        for n in srange(-l+1, l):
            NN[n+l, n+l] = -q**(1 + n) * bb(l, n) \
                / ((1 - z**2)*(z - z**(-1)))

    return NN

def matM(l, i):
    q, z = var('q z')

    MM = matrix(SR, 2*l+1, 2*l+1)

    if i == 1:
        for n in srange(-l, l+1):
            MM[n+l, n+l] = q**(-n) * f(z)
    elif i == 2:
        for n in srange(-l, l+1):
            MM[n+l, n+l] = q**(n) * f(z)

    return MM

def MapleSimplifyMatrix(mat):
    # First apply maxima
    mat = mat.apply_map(lambda elm : elm.simplify())
    # Then Maple
    return mat.apply_map(lambda elm : maple(elm).simplify().sage())

def SageSimplifyMatrix(mat, do_fac=True):
    def fsimplify(elm):
        elm = elm.full_simplify()
        if do_fac and elm != 0:
            elm = elm.factor()
        return elm

    return mat.apply_map(fsimplify)

def matE(l):
    Y = matrix(SR, 2*l+1, 2*l+1)

    if l in ZZ:
        for i in srange(l):
            Y[i,i] = 1
            Y[i, 2*l+1-i-1] = 1
            Y[2*l+1-i-1, i] = -1
            Y[2*l+1-i-1, 2*l+1-i-1] = 1
        Y[l,l] = sqrt(2)
    else:
        for i in srange(l+1/2):
            Y[i,i] = 1
            Y[i, 2*l+1-i-1] = 1
            Y[2*l+1-i-1, i] = -1
            Y[2*l+1-i-1, 2*l+1-i-1] = 1

    return 1/sqrt(2)*Y

def sage2maple(m):
	print 'Matrix(['
	
	height = len(m.columns())
	width = len(m.rows())

	for i in range(height):
		print '[',
		for j in range(width):
			if j != width - 1:
				print '%s,' % m[i,j], 
			else:
				print m[i,j],
		print '],'

	print ']);'

def replace_z_x(polynomial):
    q, z, x = var('q z x')

    ret_poly = 0

    while polynomial != 0:
        coef = polynomial.coefficients(z)

        if len(coef) == 1:
            if coef[0][1] != 0:
                return None
            return expand(ret_poly + polynomial)

        first, last = coef[0], coef[-1]

        if not (first[0] == last[0] and first[1] == -last[1]):
            return None
        
        ret_poly += first[0]*(2*x)**(last[1])
        polynomial -= first[0]*(z + z**(-1))**(last[1])

    return ret_poly

def Weight(l, simplify_algorithm='sage', replace=True):
    phi0 = Phi0(l)
    phi0_star = Phi0_star(l)
    weight = phi0_star.substitute({z : q**(-2)*z})*phi0
    assume(q > 0)
    assume(z > 0)
    weight = weight.apply_map(lambda elm : elm.simplify())

    if simplify_algorithm == 'maple':
        weight = MapleSimplifyMatrix(weight)
    elif simplify_algorithm == 'sage':
        weight = SageSimplifyMatrix(weight)

    if replace:
        weight = weight.apply_map(lambda elm : replace_z_x(elm))

    if simplify_algorithm == 'maple':
        weight = MapleSimplifyMatrix(weight)
    elif simplify_algorithm == 'sage':
        weight = SageSimplifyMatrix(weight)

    return weight

def Triangulize(M):
    l = (len(M.rows()) - 1)/2
    E = matE(l)
    return E*M*E.transpose()

def WeightsFromTo(m, n, file_name='weights.pkl', NUM_CORES=12):
    ret = {}
    
    @parallel(NUM_CORES)
    def par_weight(l):
        return Weight(l, replace=False)

    gen = par_weight(range(m, n))
    for rec in gen:
        print 'Computing weight', rec[0][0][0]/2
        ret[rec[0][0][0]] = rec[1]
        save(rec, 'weights')

    return ret

def baseChebychev_U(poly):
    deg = poly.degree(x)
    ret = []

    while poly != 0:
       pass 


# Three term recurrence relation

def aa(l1, l2, m1, m2, l):
    ret = 0

    for j1 in srange(-l1, l1+1):
        for j2 in srange(-l2, l2+1):
            j = j1 - j2
            for i1 in [-1/2, 1/2]:
                i2 = i1
                n1 = j1 - i1
                n2 = j2 - i2
                if n1 - n2 == l and n1 in srange(-m1, m1+1) and n2 in srange(-m2, m2+1):
                    ret += cgc(l1, l2, l, j1, j2, j) \
                        * cgc(1/2, 1/2, 0, i1, i2, 0) \
                        * cgc(l1, 1/2, m1, j1, i1, n1) \
                        * cgc(l2, 1/2, m2, j2, i2, n2) \
                        * cgc(m1, m2, l, n1, n2, l)

    return ret

def zeta(l, d, k):
    return (d + l + k)/2, (d + l - k)/2

def Ad(l, d):
    ret = matrix(SR, 2*l+1, 2*l+1)

    for k in srange(-l, l+1):
        l1, l2 = zeta(l, d, k)
        m1, m2 = zeta(l, d+1, k)
        ret[l+k, l+k] = aa(l1, l2, m1, m2, l)**2

    return ret

def Bd(l, d):
    ret = matrix(SR, 2*l+1, 2*l+1)

    for k in srange(-l, l):
        l1, l2 = zeta(l, d, k)
        m1, m2 = zeta(l, d, k+1)
        ret[l+k, l+k+1] = aa(l1, l2, m1, m2, l)**2

    for k in srange(-l+1, l+1):
        l1, l2 = zeta(l, d, k)
        m1, m2 = zeta(l, d, k-1)
        ret[l+k, l+k-1] = aa(l1, l2, m1, m2, l)**2

    return ret

def Cd(l, d):
    ret = matrix(SR, 2*l+1, 2*l+1)

    for k in srange(-l, l+1):
        l1, l2 = zeta(l, d, k)
        m1, m2 = zeta(l, d-1, k)
        ret[l+k, l+k] = aa(l1, l2, m1, m2, l)**2

    return ret

def try_poly_3():
    print 'Compute P0'
    P0 = identity_matrix(SR, 3)

    print 'Compute P1'
    A0, B0, C0 = Ad(1, 0), Bd(1, 0), Cd(1, 0)
    P1 = (A0**(-1) * B0 - x * A0**(-1)) * P0

    print 'Compute P2'
    A1, B1, C1 = Ad(1, 1), Bd(1, 1), Cd(1, 1)
    P2 = (A1**(-1) * B1 - x * A1**(-1)) * P1 + A1**(-1) * C1 * P0

    print 'Compute P3'
    A2, B2, C2 = Ad(1, 2), Bd(1, 2), Cd(1, 2)
    P3 = (A2**(-1) * B2 - x * A2**(-1)) * P2 + A2**(-1) * C2 * P1

    return P0, P1, P2, P3

def try_poly_2():
    print 'Compute P0'
    P0 = identity_matrix(SR, 2)

    print 'Compute p1'
    A0, B0, C0 = Ad(1/2, 0), Bd(1/2, 0), Cd(1/2, 0)
    P1 = (A0**(-1) * B0 - x * A0**(-1)) * P0

    return P0, P1

def Ea3(l1, l2):
    q, z = var('q z')

    M = matrix(SR, [
        [cgc(0, 1, 1, 0, 1, -1)**2 * q**(-1)*z**(-1), cgc(1/2, 1/2, 1, -1/2, 1/2, -1)**2, cgc(1, 0, 1, -1, 0, -1)**2*q*z],
        [cgc(0, 1, 1, 0, 0, 0)**2, cgc(1/2, 1/2, 1, 1/2, 1/2, 0)**2 * q**(-1)*z**(-1) + cgc(1/2, 1/2, 1, -1/2, -1/2, 0)**2 * q*z, cgc(1, 0, 1, 0, 0, 0)**2],
        [cgc(0, 1, 1, 0, -1, 1)**2 * q*z, cgc(1/2, 1/2, 1, 1/2, -1/2, 1)**2, cgc(1, 0, 1, 1, 0, 1)**2 * q**(-1)*z**(-1)]])
   
    a1, a2, a3 = 0, 0, 0

    for n1 in srange(-l1, l1+1):
        for n2 in srange(-l2, l2+1):
            a1 += cgc(l1, l2, 1, n1, n2, -1)**2 * q**(-n1-n2)*z**(-n1-n2)
            a2 += cgc(l1, l2, 1, n1, n2, 0)**2 * q**(-n1-n2)*z**(-n1-n2)
            a3 += cgc(l1, l2, 1, n1, n2, 1)**2 * q**(-n1-n2)*z**(-n1-n2)

    N = matrix(SR, [[a1], [a2], [a3]])
    
    return M**(-1)*N

def Qd3(d):
    ret = matrix(SR, 3, 3)

    l1, l2 = zeta(1, d, -1)
    vec = Ea3(l1, l2)
    q1, q2, q3 = vec[0,0], vec[1,0], vec[2,0]
    ret[0,0] = q1
    ret[1,0] = q2
    ret[2,0] = q3

    l1, l2 = zeta(1, d, 0)
    vec = Ea3(l1, l2)
    q1, q2, q3 = vec[0,0], vec[1,0], vec[2,0]
    ret[0,1] = q1
    ret[1,1] = q2
    ret[2,1] = q3

    l1, l2 = zeta(1, d, 1)
    vec = Ea3(l1, l2)
    q1, q2, q3 = vec[0,0], vec[1,0], vec[2,0]
    ret[0,2] = q1
    ret[1,2] = q2
    ret[2,2] = q3

    return ret

