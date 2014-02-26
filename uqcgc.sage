load('basichypergeometric.sage')
from qpolynomials import dual_q_Hahn_polynomials
from qfunctions import qbinomial
import sys

q = var('q')

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

        #print line1

        line2 = sqrt((1 - q**(-4*l-2)) \
            / (1 - q**(-4*l1 - 4*l2 - 2))) \
            * q**(-2*l1*(2*l2+p) + (l1+l2-l)*(2*l1+2*l2+p) \
            - binomial(l1 + l2 - l, 2))

        #print line2

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
        #print line3 
 
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

def bottom_cgc(l1, l2, l, i, j, m):
    k = 2*l1
    p = l1 - l2 - m
    q = var('q') 

    #print
    if -2*l2 <= p and p < 0:
        #print '1'
        line1 = qpoch(q**(2*k - 4*l), q**2, l - m)
        line2 = qpoch([q**(2*m - 2*l), q**(-2*k)], q**2, k/2 - i)
        line3 = qpoch(q**(-4*l), q**2, l - m)
        line4 = qpoch([q**2, q**(-2*k + 2*l + 2*m + 2)], q**2, k/2 - i) 
        #qpower = q**(-(2*i - k)*(k - 2*l)*(2*l + 1)*(l + m))
        qpower = q**(-(2*i - k)*(2*l + 1) - 2*k*(l - m))

        #print line1
        #print line2
        #print line3
        #print line4
        #print qpower

        return qpower * line1 * line2 / (line3 * line4)

    elif 0 <= p and p <= 2*l1 - 2*l2:
        #print '2'
        #line1 = qpoch(q**(2*m - 2*l), q**2, 2*l - k)
        line1 = qpoch(q**(-2*k), q**2, l + m)
        line2 = qpoch([q**(2*k - 4*l), q**(-2*l - 2*m)], q**2, m + l - k/2 - i)
        #line3 = qpoch(q**(-4*l), q**2, 2*l - k)
        line3 = qpoch(q**(-4*l), q**2, l + m)
        line4 = qpoch([q**2, q**(2*k - 2*l - 2*m + 2)], q**2, m + l - k/2 - i)
        qpower = q**(-(2*i + k - 2*l - 2*m)*(2*l + 1) + 2*(k - 2*l)*(l + m))
        #qpower = q**(2*(k - 2*l)*(k - 2*l) + (2*l + 1)*(2*l + 2*m - k - 2*i))

        return qpower * line1 * line2 / (line3 * line4)

    elif 2*l1 - 2*l2 < p and p <= 2*l1:
        #print '3'
        line1 = qpoch(q**(-2*k), q**2, l + m)
        line2 = qpoch([q**(2*k - 4*l), q**(-2*l - 2*m)], q**2, m + l - k/2 - i)
        line3 = qpoch(q**(-4*l), q**2, l + m)
        line4 = qpoch([q**2, q**(2*k - 2*l - 2*m + 2)], q**2, m + l - k/2 - i)

        #qpower = q**(2*(l + m)*(k - 2*l) + (1 + 2*l)*(2*l + 2*m - k - 2*i))
        qpower = q**(-(2*i + k - 2*l - 2*m)*(2*l + 1) + 2*(k - 2*l)*(l + m))

        return qpower * line1 * line2 / (line3 * line4)

    else:
        return 0

def bottom_cgc2(l1, l2, l, i, j, m):
    k = 2*l1 - l
    p = l1 - l2 - m
    q = var('q') 

    #print
    if -2*l2 <= p and p < 0:
        #print '1'
        line1 = qpoch(q**(2*k - 2*l), q**2, l - m)
        line2 = qpoch([q**(2*m - 2*l), q**(-2*k - 2*l)], q**2, k/2 + l/2 - i)
        line3 = qpoch(q**(-4*l), q**2, l - m)
        line4 = qpoch([q**2, q**(-2*k + 2*m + 2)], q**2, k/2 + l/2 - i) 
        qpower = q**(-(2*i - k - l)*(2*l + 1) - (2*k + 2*l)*(l - m))

        #print line1
        #print line2
        #print line3
        #print line4
        #print qpower

        return qpower * line1 * line2 / (line3 * line4)

    elif 0 <= p and p <= 2*l1:
        #print '2'
        line1 = qpoch(q**(-2*k - 2*l), q**2, l + m)
        line2 = qpoch([q**(2*k - 2*l), q**(-2*l - 2*m)], q**2, m + l/2 - k/2 - i)
        line3 = qpoch(q**(-4*l), q**2, l + m)
        line4 = qpoch([q**2, q**(2*k - 2*m + 2)], q**2, m + l/2 - k/2 - i)
        qpower = q**(-(2*i + k - l - 2*m)*(2*l + 1) + (2*k - 2*l)*(l + m))

        return qpower * line1 * line2 / (line3 * line4)

    else:
        return 0

def bottom_cgc3(l1, l2, l, i, j, m):
    k = 2*l1 - l
    qpower = q**(2*(j - l2)*(i - l1))
    
    #qbin1 = qbinomial(l - m, l2 + j, q**2)
    #qbin2 = qbinomial(l + m, l1 + i, q**2)
    #qbin3 = qbinomial(2*l, l - k, q**2)
   
    qbin1 = qbinomial(l - k, l2 + j, q**2)
    qbin2 = qbinomial(l + k, l1 + i, q**2)
    qbin3 = qbinomial(2*l, l - m, q**2)

    return qpower * qbin1 * qbin2 / qbin3

def bottom_cgc_guess(l1, l2, l, i, j, k):
    m = 2*l1 - l

    qbin1 = qbinomial(l + m, i + l1, q**2)
    qbin2 = qbinomial(l - m, j + l2, q**2)
    qbin3 = qbinomial(2*l, l - k, q**2)

    return qbin1 * qbin2 / qbin3

def test_cgc_bottom(l):
    q = var('q')
    assume(q > 0)

    for k in srange(2*l+1):
        l1, l2 = k/2, l - k/2
        for i in srange(-l1, l1+1):
            for j in srange(-l2, l2+1):
                m = i - j
                c1 = try_full_simplify(cgc(l1, l2, l, i, j, m)**2) 
                c2 = try_full_simplify(bottom_cgc(l1, l2, l, i, j, m))
                if bool(c1 != c2):
                    print l1, l2, l, i, j, m 
                    print c1
                    print c2
                    return False
                sys.stdout.write('.')
                sys.stdout.flush()
    return True

def test_cgc_bottom2(l):
    q = var('q')
    assume(q > 0)

    for k in srange(2*l+1):
        l1, l2 = k/2, l - k/2
        for i in srange(-l1, l1+1):
            for j in srange(-l2, l2+1):
                m = i - j
                c1 = try_full_simplify(cgc(l1, l2, l, i, j, m)**2) 
                c2 = try_full_simplify(bottom_cgc2(l1, l2, l, i, j, m))
                if bool(c1 != c2):
                    print l1, l2, l, i, j, m 
                    print c1
                    print c2
                    return False
                sys.stdout.write('.')
                sys.stdout.flush()
    return True

def test_cgc_bottom3(l):
    q = var('q')
    n = var('n')
    assume(q > 0)

    for k in srange(2*l+1):
        l1, l2 = k/2, l - k/2
        for i in srange(-l1, l1+1):
            for j in srange(-l2, l2+1):
                m = i - j
                c1 = try_full_simplify(cgc(l1, l2, l, i, j, m)**2)
                c2 = try_full_simplify(bottom_cgc_guess(l1, l2, l, i, j, m))
                #ret = solve([try_full_simplify(c1 / c2) == q**n], n, solution_dict=True)[0][n]
                #print '    poly.substitute({l1 : %s, l2 : %s, l : %s, i : %s, j : %s, k : %s}) == %s,' % (l1, l2, l, i, j, k, ret)
                print try_full_simplify(c1 / c2),
                sys.stdout.flush()

def test_cgc_bottom4(l):
    q = var('q')
    assume(q > 0)

    for k in srange(2*l+1):
        l1, l2 = k/2, l - k/2
        for i in srange(-l1, l1+1):
            for j in srange(-l2, l2+1):
                m = i - j
                c1 = try_full_simplify(cgc(l1, l2, l, i, j, m)**2) 
                c2 = try_full_simplify(bottom_cgc3(l1, l2, l, i, j, m))
                if bool(c1 != c2):
                    print l1, l2, l, i, j, m 
                    print c1
                    print c2
                    return False
                sys.stdout.write('.')
                sys.stdout.flush()
    return True


def test_cgc(l1, l2):
    from itertools import product

    saved_cgc = {}

    # test 1
    for n1, m1 in product(srange(-l1, l1+1), repeat=2):
        for n2, m2 in product(srange(-l2, l2+1), repeat=2):
            s = 0
            for l in srange(abs(l1-l2), l1+l2+1):
                for n in srange(-l, l+1):

                    if (l1, l2, l, n1, n2, n) in saved_cgc.keys():
                        cgc1 = saved_cgc[(l1, l2, l, n1, n2, n)]
                    else:
                        cgc1 = cgc(l1, l2, l, n1, n2, n)
                        saved_cgc[(l1, l2, l, n1, n2, n)] = cgc1

                    if (l1, l2, l, m1, m2, n) in saved_cgc.keys():
                        cgc2 = saved_cgc[(l1, l2, l, m1, m2, n)]
                    else:
                        cgc2 = cgc(l1, l2, l, m1, m2, n)
                        saved_cgc[(l1, l2, l, m1, m2, n)] = cgc2                        
                    s += cgc1 * cgc2

            if bool(s != kronecker_delta(n1, m1)*kronecker_delta(n2, m2)):
                for l in srange(abs(l1-l2), l1+l2+1):
                    for n in srange(-l, l+1):

                        if (l1, l2, l, n1, n2, n) in saved_cgc.keys():
                            cgc1 = saved_cgc[(l1, l2, l, n1, n2, n)]
                        else:
                            cgc1 = cgc(l1, l2, l, n1, n2, n)
                            saved_cgc[(l1, l2, l, n1, n2, n)] = cgc1

                        if (l1, l2, l, m1, m2, n) in saved_cgc.keys():
                            cgc2 = saved_cgc[(l1, l2, l, m1, m2, n)]
                        else:
                            cgc2 = cgc(l1, l2, l, m1, m2, n)
                            saved_cgc[(l1, l2, l, m1, m2, n)] = cgc2                        

                        ret = cgc1 * cgc2
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
                if (l1, l2, l, n1, n2, n) in saved_cgc.keys():
                    cgc1 = saved_cgc[(l1, l2, l, n1, n2, n)]
                else:
                    cgc1 = cgc(l1, l2, l, n1, n2, n)
                    saved_cgc[(l1, l2, l, n1, n2, n)] = cgc1

                if (l1, l2, ll, m1, m2, n) in saved_cgc.keys():
                    cgc2 = saved_cgc[(l1, l2, ll, n1, n2, m)]
                else:
                    cgc2 = cgc(l1, l2, ll, n1, n2, m)
                    saved_cgc[(l1, l2, ll, n1, n2, m)] = cgc2                        
                s += cgc1 * cgc2
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

def dg_base_chebyshev_U(poly, dg):
    x = var('x')

    if dg > poly.degree(x):
        return 0, poly

    if poly.degree(x) == 0:
        if dg == 0:
            return poly, 0
        else:
            return 0, poly

    if dg == 0:
        return poly.coefficient(x, 0), poly - poly.coefficient(x, 0)

    ch_U = chebyshev_U(dg, x).simplify_full()
    lc_U = ch_U.leading_coefficient(x)
    lc_p = poly.full_simplify().leading_coefficient(x)
    
    return (lc_p/lc_U).simplify_full(), (poly - lc_p/lc_U*ch_U).simplify_full()

def base_U(poly):
    x = var('x')
    dg = poly.degree(x)
    ret = [0]*(dg+1)
    res = poly

    for n in range(dg, -1, -1):
        ret[n], res = dg_base_chebyshev_U(res, n)

    return ret

def try_full_simplify(f):
    try:
        return f.full_simplify().factor()
    except:
        return f

def explicit_weight_elm2(l, t, m, n):
    k = 2*l + 2*t - n - m
    s = (2*l - n + m - k)/2
    coef = q**(-k + 2*s*(2*l - s - k)) 
    line1 = (1 - q**(4*l + 2)) / (1 - q**(2*n + 2))
    line2 = qpoch(q**2, q**2, 2*l - m) \
        * qpoch(q**2, q**2, m) \
        / qpoch(q**2, q**2, 2*l)
    line3 = qpoch(q**(4*l - 2*n), q**(-2), m - t) \
        / qpoch(q**(2*n + 4), q**2, m - t)
    line4 = qpoch(q**(4*l + 4 - 2*t), q**2, t) \
        / qpoch(q**2, q**2, t)

    return coef*line1*line2*line3*line4

def explicit_weight_elm(l, t, m, n):
    ret = (-1)**(m - t) * q**(2*(2*l + 1)*m - m^2 - (4*l + 3)*t + t^2 - 2*l + n) \
        * (1 - q**(4*l + 2)) / (1 - q**(2*n + 2)) \
        * qpoch(q**2, q**2, 2*l - m) \
        * qpoch(q**2, q**2, m) \
        / qpoch(q**2, q**2, 2*l) \
        * qpoch(q**(2*n - 4*l), q**2, m - t) \
        / qpoch(q**(2*n + 4), q**2, m - t) \
        * qpoch(q**(4*l + 4 - 2*t), q**2, t) \
        / qpoch(q**2, q**2, t)

    return ret

def test_weights():
    weights = load('weights3.sobj')
    for i in range(5): 
        print 'Test weight %s' % (i/2)
        if not bool(SageSimplifyMatrix(weights[i]) == SageSimplifyMatrix(explicit_weight(i/2))):
            print 'Weight %s failed!' % (i/2)
            return False
    return True

def matrix_base_U(mat, n):
    def get_base_U_n(elm, n):
        try:
            return base_U(elm)[n]
        except:
            return 0

    return mat.apply_map(lambda elm : get_base_U_n(elm, n))

def explicit_weight(l):
    x = var('x')
    W = matrix(SR, 2*l + 1, 2*l + 1)

    for m in range(2*l + 1):
        for n in range(m+1):
            for t in range(n+1):
                W[m, n] += explicit_weight_elm(l, t, m, n) \
                    * chebyshev_U(n + m - 2*t, x)
                W[n, m] = W[m, n]

    return W

def test_all():
    print 'test functions:'
    print 'test1'
    if not test_tests(test0, test1):
        return False
    print 'test2'
    if not test_tests(test1, test2):
        return False
    print 'test3'
    if not test_tests(test2, test3):
        return False
    print 'test4'
    if not test_tests(test3, test4):
        return False
    
    print 'test weights'
    if not test_weights():
        return False

    return True


def f1(l, i, j, k):
    line1 = qpoch(q**2, q**2, k) * qpoch(q**2, q**2, 2*l - k) \
        * qpoch(q**2, q**2, l - j) * qpoch(q**2, q**2, l + j)
    line2 = qpoch(q**2, q**2, i + k/2) * qpoch(q**2, q**2, k/2 - i) \
        * qpoch(q**2, q**2, j + l - k/2) * qpoch(q**2, q**2, l - k/2 - j) \
        * qpoch(q**2, q**2, 2*l)

    return line1 / line2

def f2(l, i, j, k):
    l1 = k/2
    l2 = l - k/2
    m = i - j
    p = k - l - m
    
    line1 = qpoch(q**(2*k - 4*l), q**2, l - m) * qpoch([q**(2*m - 2*l), q**(-2*k)], q**2, k/2 - i)
    line2 = qpoch(q**(-4*l), q**2, l - m) * qpoch([q**2, q**(-2*k + 2*l + 2*m + 2)], q**2, k/2 - i)

    return line1 / line2

def test0(l, p, qq):
    z = var('z')

    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret = 0
    for j1 in srange(-l1, l1+1):
        for j2 in srange(-l2, l2+1):
            for i1 in srange(-m1, m1+1):
                for i2 in srange(-m2, m2+1):
                    for j in srange(-l, l+1):
                        ret += q**(-j1 - j2 - i1 - i2) \
                            * cgc(l1, l2, l, j1, j2, j)**2 \
                            * cgc(m1, m2, l, i1, i2, j)**2 \
                            * z**(-i1 - i2 + j1 + j2)
    
    return ret
    

def test1(l, p, qq):
    z = var('z')

    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret = 0
    for j1 in srange(-l1, l1+1):
        for i1 in srange(-m1, m1+1):
            for j in srange(-l, l+1):
                ret += q**(-2*i1 - 2*j1 + 2*j) \
                    * cgc(l1, l2, l, j1, j1 - j, j)**2 \
                    * cgc(m1, m2, l, i1, i1 - j, j)**2 * z**(2*j1 - 2*i1)

    return ret

def test15(l, p, qq):
    z = var('z')

    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret = 0
    for j1 in srange(-l1, l1+1):
        for i1 in srange(-m1, m1+1):
            for j in srange(-l, l+1):
                ret += q**(-2*i1 + 2*j1 + 2*j) \
                    * cgc(l1, l2, l, -j1, -j1 - j, j)**2 \
                    * cgc(m1, m2, l, i1, i1 - j, j)**2 * z**(-2*j1 - 2*i1)

    return ret

def test2(l, p, qq):
    z = var('z')

    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret = 0
    for j1 in srange(-l1, l1+1):
        for i1 in srange(-m1, m1+1):
            for j in srange(-l, l+1):
                ret += q**(2*(-i1 + j1 + j) + 2*(j1 + l1)*(j1 + j + l2) \
                        - 2*(j1 - l1)*(j1 + j - l2)) \
                    * cgc(l1, l2, l, j1, j1 + j, -j)**2 \
                    * cgc(m1, m2, l, i1, i1 - j, j)**2 * z**(2*(-j1 - i1))

    return ret

def test_tests(t1, t2):
    for l in srange(1/2, 5/2):
        print l, -1/2, -1/2, bool(t1(l, -1/2, -1/2) == t2(l, -1/2, -1/2))
        if not bool(t1(l, -1/2, -1/2) == t2(l, -1/2, -1/2)):
            return False
    for l in srange(3/2, 7/2):
        print l, -1/2, -3/2, bool(t1(l, -1/2, -3/2) == t2(l, -1/2, -3/2))
        if not bool(t1(l, -1/2, -3/2) == t2(l, -1/2, -3/2)):
            return False
    for l in range(2, 4):
        print l, -1, -2, bool(t1(l, -1, -2) == t2(l, -1, -2))
        if not bool(t1(l, -1, -2) == t2(l, -1, -2)):
            return False
    return True

def test3(l, p, qq):
    z = var('z')

    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret = 0
    for j1 in srange(-l1, l1+1):
        for i1 in srange(-m1, m1+1):
            for j in srange(-l, l+1):
                if cgc(l1, l2, l, j1, j1 + j, -j) == 0 or cgc(m1, m2, l, i1, i1 - j, j) == 0:
                    continue
                ret += q**(2*(-i1 + j1 + j) + 2*(j1 + l1)*(j1 + j + l2) \
                        - 2*(j1 - l1)*(j1 + j - l2) ) \
                    * q**(2*(j1 - l1)*(j1 + j - l2)) \
                    * qbinomial(2*l1, l1 + j1, q**2) \
                    * qbinomial(2*l2, l2 + j1 + j, q**2) \
                    / qbinomial(2*l, l + j, q**2) \
                    * q**(2*(i1 - m1)*(i1 - j - m2)) \
                    * qbinomial(2*m1, m1 + i1, q**2) \
                    * qbinomial(2*m2, m2 + i1 - j, q**2) \
                    / qbinomial(2*l, l - j, q**2) * z**(2*(-j1 - i1))

    return ret

def test4(l, p, qq):
    z = var('z')

    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for j in srange(-l1, l1+1):
        for i in srange(-m1, m1+1):
            ret2 = 0
            for k in srange(max(-j - l2, i - m2), min(-j + l2, i + m2) + 1):
                ret2 += q**(2*k + 2*(j + l1)*k - 2*(i - m1)*k) \
                    * qbinomial(2*l2, l2 - j - k, q**2) \
                    * qbinomial(2*m2, m2 - i + k, q**2) \
                    / qbinomial(2*l, l - k, q**2)**2
            ret1 += q**(2*(-i + j) + 2*(j + l1)*(j + l2) + 2*(i - m1)*(i - m2)) \
                * qbinomial(2*l1, l1 + j, q**2) \
                * qbinomial(2*m1, m1 + i, q**2) \
                * ret2 * z**(-2*(i + j))

    return ret1

def funcF(l, i, j, p, qq):
    '''F from the definition'''
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret = 0
    for k in srange(max(-j - l2, i - m2), min(-j + l2, i + m2) + 1):
        ret += q**(2*(1 + j + l1)*k - 2*(i - m1)*k) \
            * qbinomial(2*l2, l2 - j - k, q**2) \
            * qbinomial(2*m2, m2 - i + k, q**2) \
            / qbinomial(2*l, l - k, q**2)**2
    ret *= q**(2*(-i + j) + 2*(j + l1)*(j + l2) + 2*(i - m1)*(i - m2)) \
        * qbinomial(2*l1, l1 + j, q**2) \
        * qbinomial(2*m1, m1 + i, q**2) 

    return ret

def test5(l, p, qq):
    z = var('z')

     
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for r in srange(-(l + (p + qq)/2), l + (p + qq)/2 + 1):
        ret2 = 0
        for i in srange(max(-m1, r - l1), min(m1, r + l1) + 1):
            ret2 += funcF(l, i, r - i, p, qq)
        ret1 += ret2*z**(-2*r)

    return ret1

def funcd1(l, r, p, qq):
    '''replace j -> r - i'''
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret = 0
    for i in srange(max(-m1, r - l1), min(m1, r + l1) + 1):
        ret += funcF(l, i, r - i, p, qq)
    
    return ret

def test6(l, p, qq):
    z = var('z')

     
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret = 0
    for r in srange(-(l + (p + qq)/2), l + (p + qq)/2 + 1):
        ret += funcd1(l, r, p, qq)*z**(-2*r)

    return ret

def test7(l, p, qq):
    z = var('z')

     
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2
    
    ret = 0
    for s in srange(-l - p, l + p + 1):
        ret += funcd1(l, s + (p - qq)/2, p, qq)*z**(-2*(s + (p - qq)/2))

    return ret

def funcd2(l, s, p, qq):
    '''Substitution r(s) -> s + (p - q)/2'''
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for i in srange(max(-m1, s + (p - qq)/2 - l1), min(m1, s + (p - qq)/2 + l1) + 1):
        ret1 += funcF(l, i, s + (p - qq)/2 - i, p, qq)

    #for n in range(l + qq - s + 1):
    #    ret2 = 0
    #    for k in srange(n + s - l, n - p + 1):
    #        ret2 += q**(2*k*(l + p - n) - 2*k*(n + s - l - qq)) \
    #            * qbinomial(l - p, l - n + k, q**2) \
    #            * qbinomial(l - q, n + s - k - qq, q**2) \
    #            / qbinomial(2*l, l - k, q**2)**2 
    #    ret1 += q**(2*(s + (p - q)/2) + 2*(l + p - n)*(l - n) + 2*(n + s)*(n + s - l)) \
    #        * qbinomial(l + p, l + p - n, q**2) \
    #        * qbinomial(l + q, n + s, q**2) * ret2

    return ret1

def test8(l, p, qq):
    z = var('z')

     
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2
    
    ret = 0
    for s in srange(-l - p, l + p + 1):
        ret += funcd2(l, s, p, qq)*z**(-2*(s + (p - qq)/2))

    return ret


def funcd3(l, s, p, qq):
    '''Writing out F in qbinomials.'''
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2


    ret1 = 0
    for n in srange(l + qq - s + 1):
        #ret1 += funcF(l, n + s - (l + qq)/2, (l + p)/2 - n, p, qq)
        ret2 = 0
        i, j = n + s - (l + qq)/2, (l + p)/2 - n
        for k in srange(n + s - l, n - p + 1):
            ret2 += q**(2*(1 + l + p - n)*k - 2*(n + s - l - qq)*k) \
                * qbinomial(2*l2, n - k - p, q**2) \
                * qbinomial(2*m2, l - n - s + k, q**2) \
                / qbinomial(2*l, l - k, q**2)**2
        ret2 *= q**(2*(l - 2*n + (p + qq)/2 - s) + 2*(l + p - n)*(l - n) + 2*(n + s - l - qq)*(n + s - l)) \
            * qbinomial(2*l1, l + p - n, q**2) \
            * qbinomial(2*m1, n + s, q**2) 
        ret1 += ret2

    return ret1

def funcd4(l, s, p, qq):
    '''Substituting k(m) = m + n + s - l.'''
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2


    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        i, j = n + s - (l + qq)/2, (l + p)/2 - n
        for m in srange(0, l - s - p + 1):
            k = m + n + s - l
            ret2 += q**(2*m*(2*l - 2*n + p + qq - s + 1)) \
                * qbinomial(l - p, m + s, q**2) \
                * qbinomial(l - qq, m, q**2) \
                / qbinomial(2*l, m + n + s, q**2)**2
        ret2 *= q**(-2*n*(s + 1)) \
            * qbinomial(2*l1, n, q**2) \
            * qbinomial(2*m1, n + s, q**2) 
        ret1 += ret2

    return q**(2*(l + p)*s + p + qq)*ret1

def funcd5(l, s, p, qq):
    '''Writing out all qbinomials in m-sum and move a minus sign.'''
    # Good test case: 7/2, 1, -1/2, -1/2
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        i, j = n + s - (l + qq)/2, (l + p)/2 - n
        for m in srange(0, l - s - p + 1):
            ret2 += q**(2*m*(2*l - 2*n + p + qq - s + 1)) \
                * q**((2*l - 2*p)*(m + s) - 2*binomial(m + s, 2)) \
                * qpoch(q**(2*p - 2*l), q**2, m + s) \
                / qpoch(q**2, q**2, m + s) \
                * q**((2*l - 2*qq)*m - 2*binomial(m, 2)) \
                * qpoch(q**(2*qq - 2*l), q**2, m) \
                / qpoch(q**2, q**2, m) \
                * q**(-2*(4*l)*(m + n + s) + 4*binomial(m + n + s, 2)) \
                * qpoch(q**2, q**2, m + n + s)**2 \
                / qpoch(q**(-4*l), q**2, m + n + s)**2
        ret2 *= q**(-2*n*(s + 1)) \
            * qbinomial(2*l1, n, q**2) \
            * qbinomial(2*m1, n + s, q**2) 
        ret1 += ret2

    return (-1)**s * q**(2*(l + p)*s + p + qq)*ret1

def funcd6(l, s, p, qq):
    '''Rewrite the m-sum to a basic hypergeometric series and group all q-powers.'''
    # Good test case: 7/2, 1, -1/2, -1/2
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        i, j = n + s - (l + qq)/2, (l + p)/2 - n
        #for m in srange(0, l - s - p + 1):
        #    ret2 += q**(2*m) * qpoch(q**(2*p - 2*l + 2*s), q**2, m) \
        #        / qpoch(q**(2 + 2*s), q**2, m) \
        #        * qpoch(q**(2*qq - 2*l), q**2, m) \
        #        / qpoch(q**2, q**2, m) \
        #        * qpoch(q**(2 + 2*n + 2*s), q**2, m)**2 \
        #        / qpoch(q**(-4*l + 2*n + 2*s), q**2, m)**2
        ret2 += bhs([
                q**(2*p - 2*l + 2*s),   
                q**(2*qq - 2*l),
                q**(2*n + 2*s + 2),
                q**(2*n + 2*s + 2)
            ], [
                q**(2*s + 2),
                q**(2*s + 2*n - 4*l),
                q**(2*s + 2*n - 4*l)
            ], q**2, q**2)
        ret2 *= q**(-2*n*(s + 1) - 2*(4*l + 1)*n + 2*n**2 - (6*l - 4*n + 2*p + 1)*s + s**2) \
            * qpoch(q**(2*p - 2*l), q**2, s) \
            / qpoch(q**2, q**2, s) \
            * qpoch(q**2, q**2, n + s)**2 \
            / qpoch(q**(-4*l), q**2, n + s)**2 \
            * qbinomial(2*l1, n, q**2) \
            * qbinomial(2*m1, n + s, q**2) 
        ret1 += ret2

    return (-1)**s * q**(2*(l + p)*s + p + qq)*ret1

def funcd7(l, s, p, qq):
    '''Simplifying the q-powers.'''
    # Good test case: 7/2, 1, -1/2, -1/2
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        ret2 += bhs([
                q**(2*p - 2*l + 2*s),   
                q**(2*qq - 2*l),
                q**(2*n + 2*s + 2),
                q**(2*n + 2*s + 2)
            ], [
                q**(2*s + 2),
                q**(2*s + 2*n - 4*l),
                q**(2*s + 2*n - 4*l)
            ], q**2, q**2)
        ret2 *= q**(2*n**2 + 2*n*(-4*l + s - 2))  \
            * qpoch(q**(2*p - 2*l), q**2, s) \
            / qpoch(q**2, q**2, s) \
            * qpoch(q**2, q**2, n + s)**2 \
            / qpoch(q**(-4*l), q**2, n + s)**2 \
            * qbinomial(2*l1, n, q**2) \
            * qbinomial(2*m1, n + s, q**2) 
        ret1 += ret2

    return (-1)**s * q**(-4*l*s + s**2 + p + qq - s) * ret1

def funcd8(l, s, p, qq):
    '''Substitute n = -n + l + qq - s.'''
    # Good test case: 7/2, 1, -1/2, -1/2
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        ret2 += bhs([
                q**(2*p - 2*l + 2*s),   
                q**(2*qq - 2*l),
                q**(2*l + 2*qq - 2*n + 2),
                q**(2*l + 2*qq - 2*n + 2)
            ], [
                q**(2*s + 2),
                q**(-2*l - 2*n + 2*qq),
                q**(-2*l - 2*n + 2*qq)
            ], q**2, q**2)
        ret2 *= q**(-2*(3*l + n - qq + 2)*(l - n + qq - s))
        ret2 *= qpoch(q**(2*p - 2*l), q**2, s) \
            / qpoch(q**2, q**2, s) \
            * qpoch(q**2, q**2, -n + l + qq)**2 \
            / qpoch(q**(-4*l), q**2, -n + l + qq)**2 \
            * qbinomial(2*l1, -n + l + qq - s, q**2) \
            * qbinomial(2*m1, - n + l + qq, q**2) 
        ret1 += ret2

    return (-1)**s * q**(-4*l*s + s**2 + p + qq - s) * ret1

def test_funcd(f1, f2):
    for i in range(3):
        print 4, i, bool(f1(4, i, -1, -2) == f2(4, i, -1, -2))
    for i in range(3):
        print 7/2, i, bool(f1(7/2, i, -1/2, -3/2) == f2(7/2, i, -1/2, -3/2))
    for i in range(2):
        print 3, i, bool(f1(3, i, -1, -2) == f2(3, i, -1, -2))

def funcd9(l, s, p, qq):
    '''Expand all qpochhammer sybols and try to group n's.'''
    # Good test case: 7/2, 2, -1/2, -3/2
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        ret2 += bhs([
                q**(2*p - 2*l + 2*s),   
                q**(2*qq - 2*l),
                q**(2*l + 2*qq - 2*n + 2),
                q**(2*l + 2*qq - 2*n + 2)
            ], [
                q**(2*s + 2),
                q**(-2*l - 2*n + 2*qq),
                q**(-2*l - 2*n + 2*qq)
            ], q**2, q**2)
        ret2 *= q**(-2*(3*l + n - qq + 2)*(l - n + qq - s)) * (-1)**(-n + l + qq - s) * q**((2*l + 2*p)*(-n + l + qq - s) - 2*binomial(-n + l + qq - s, 2)) 
        #ret2 *= qpoch(q**(2*p - 2*l), q**2, s) \
        #    / qpoch(q**2, q**2, s) \
        #    * qpoch(q**2, q**2, l + qq)**2 * qpoch(q**(2 + 2*l + 2*qq), q**2, -n)**2 \
        #    / qpoch(q**(-4*l), q**2, l + qq)**2 / qpoch(q**(-2*l + 2*qq), q**2, -n)**2 \
        #    * qpoch(q**(2*p + 2*n - 2*qq + 2*s + 2), q**2, l + qq - s) * qpoch(q**(2*p + 2*n + 2 + 2*l), q**2, -n)\
        #    / qpoch(q**2, q**2, l + qq - s) / qpoch(q**(2 + 2*l + 2*qq - 2*s), q**2, -n) \
        #    * qpoch(q**2, q**2, l + qq) * qpoch(q**(2*n + 2), q**2, -n) \
        #    / qpoch(q**2, q**2, l + qq) / qpoch(q**(2 + 2*l + 2*qq), q**2, -n)
        ret2 *= qpoch(q**(2*p - 2*l), q**2, s) \
            / qpoch(q**2, q**2, s) \
            * qpoch(q**2, q**2, l + qq)**2 * qpoch(q**(2 + 2*l + 2*qq), q**2, -n)**2 \
            / qpoch(q**(-4*l), q**2, l + qq)**2 / qpoch(q**(-2*l + 2*qq), q**2, -n)**2 \
            * qpoch(q**(-2*l - 2*p), q**2, -n + l + qq - s) \
            / qpoch(q**2, q**2, l + qq - s) / qpoch(q**(2 + 2*l + 2*qq - 2*s), q**2, -n) \
            * qpoch(q**2, q**2, l + qq) / qpoch(q**2, q**2, n) \
            / qpoch(q**2, q**2, l + qq) / qpoch(q**(2 + 2*l + 2*qq), q**2, -n)
        ret1 += ret2

    return (-1)**s * q**(-4*l*s + s^2 + p + qq - s) * ret1

def funcd10(l, s, p, qq):
    '''Rewrite -n's to n's.'''
    # Good test case: 7/2, 2, -1/2, -3/2
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        ret2 += bhs([
                q**(2*p - 2*l + 2*s),   
                q**(2*qq - 2*l),
                q**(2*l + 2*qq - 2*n + 2),
                q**(2*l + 2*qq - 2*n + 2)
            ], [
                q**(2*s + 2),
                q**(-2*l - 2*n + 2*qq),
                q**(-2*l - 2*n + 2*qq)
            ], q**2, q**2)
        ret2 *= q**(-2*(3*l + n - qq + 2)*(l - n + qq - s)) * (-1)**(-n + l + qq - s) * q**((2*l + 2*p)*(-n + l + qq - s) - 2*binomial(-n + l + qq - s, 2)) 
        ret2 *= q**(2*n*(-2*l - 2*qq) + 4*binomial(n, 2)) / qpoch(q**(-2*l - 2*qq), q**2, n)**2 \
            / q**(2*n*(2 + 2*l - 2*qq) + 4*binomial(n, 2)) * qpoch(q**(2 + 2*l - 2*qq), q**2, n)**2 \
            * (-1)**n * q**((2 + 2*p - 2*qq + 2*s)*n + 2*binomial(n, 2)) / qpoch(q**(2 + 2*p - 2*qq + 2*s), q**2, n) \
            / q**(n*(-2*l - 2*qq + 2*s) + 2*binomial(n, 2)) * (-1)**n * qpoch(q**(-2*l - 2*qq + 2*s), q**2, n) \
            / qpoch(q**2, q**2, n) \
            / q**(-n*(2*l + 2*qq) + 2*binomial(n, 2)) * (-1)**n * qpoch(q**(-2*l - 2*qq), q**2, n) 
        ret1 += ret2

    ret1 *= qpoch(q**(2*p - 2*l), q**2, s) \
        /  qpoch(q**2, q**2, s) \
        * qpoch(q**2, q**2, l + qq)**2 \
        / qpoch(q**(-4*l), q**2, l + qq)**2 \
        * qpoch(q**(-2*l - 2*p), q**2, l + qq - s) \
        / qpoch(q**2, q**2, l + qq - s) \
        * qpoch(q**2, q**2, l + qq) \
        / qpoch(q**2, q**2, l + qq)

    return (-1)**s * q**(-4*l*s + s^2 + p + qq - s) * ret1

def funcd11(l, s, p, qq):
    '''Simplifying the q-powers, group qpochhammersymbols.'''
    # Good test case: 7/2, 2, -1/2, -3/2
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        ret2 += bhs([
                q**(2*p - 2*l + 2*s),   
                q**(2*qq - 2*l),
                q**(2*l + 2*qq - 2*n + 2),
                q**(2*l + 2*qq - 2*n + 2)
            ], [
                q**(2*s + 2),
                q**(-2*l - 2*n + 2*qq),
                q**(-2*l - 2*n + 2*qq)
            ], q**2, q**2)
        #ret2 *= (-1)**(l + qq - s)
        #ret2 *= q**(-2*(3*l + n - qq + 2)*(l - n + qq - s) + (2*l + 2*p)*(-n + l + qq - s) - 2*binomial(-n + l + qq - s, 2) + 2*n*(-2*l - 2*qq) + 4*binomial(n, 2) - 2*n*(2 + 2*l - 2*qq) - 4*binomial(n, 2) + (2 + 2*p - 2*qq + 2*s)*n + 2*binomial(n, 2) - n*(-2*l - 2*qq + 2*s) - 2*binomial(n, 2) + n*(2*l + 2*qq) - 2*binomial(n, 2))
        ret2 *= q**(2*n)
        ret2 *= 1 / qpoch(q**(-2*l - 2*qq), q**2, n)**2 \
            * qpoch(q**(2 + 2*l - 2*qq), q**2, n)**2 \
            / qpoch(q**(2 + 2*p - 2*qq + 2*s), q**2, n) \
            * qpoch(q**(-2*l - 2*qq + 2*s), q**2, n) \
            / qpoch(q**2, q**2, n) \
            * qpoch(q**(-2*l - 2*qq), q**2, n) 
        ret1 += ret2

    ret1 *= qpoch(q**(2*p - 2*l), q**2, s) \
        /  qpoch(q**2, q**2, s) \
        * qpoch(q**2, q**2, l + qq)**2 \
        / qpoch(q**(-4*l), q**2, l + qq)**2 \
        * qpoch(q**(-2*l - 2*p), q**2, l + qq - s) \
        / qpoch(q**2, q**2, l + qq - s) \
        * qpoch(q**2, q**2, l + qq) \
        / qpoch(q**2, q**2, l + qq)

    ret1 *= (-1)**(l + qq) * q**(-5*l^2 + (2*l + 1)*p - 2*(2*l - p + 1)*qq + qq^2 + 2*(l - p + 1)*s - 3*l)
    return ret1

def funce1(l, s, p, qq):
    '''From the definition of sum over t.'''
    ret = 0
    r = s + (p - qq)/2

    for t in range(l + qq - s + 1):
        ret += explicit_weight_elm(l, t, l + p, l + qq)

    return ret

def funce2(l, s, p, qq):
    '''Writing out the definition and replace r(s) = s + (p - qq)/2.'''
    r = s + (p - qq)/2
    ret = 0

    for t in range(l + qq - s + 1):
        m, n = l + p, l + qq
        k = 2*l + 2*t - n - m
        ret += q**(2*l**2 + 2*l*p + 2*l*qq + 2*p*qq - 4*l*t - 2*p*t - 2*qq*t + 2*t**2 + p + qq - 2*t) \
            * (1 - q**(4*l + 2)) / (1 - q**(2*l + 2*qq + 2)) \
            * qpoch(q**2, q**2, l - p) \
            * qpoch(q**2, q**2, l + p) \
            / qpoch(q**2, q**2, 2*l) \
            * qpoch(q**(2*l - 2*qq), q**(-2), l + p - t) \
            / qpoch(q**(2*l + 2*qq + 4), q**2, l + p - t) \
            * qpoch(q**(4*l + 4 - 2*t), q**2, t) \
            / qpoch(q**2, q**2, t)

    return ret

def funce3(l, s, p, qq):
    '''Replace t to n and get out all non-independant n variables.'''
    r = s + (p - qq)/2
    ret = 0

    for n in range(l + qq - s + 1):
        ret += q**(n*(-4*l - 2*p - 2*qq - 2) + 2*n**2) \
            * qpoch(q**(2*l - 2*qq), q**(-2), l + p - n) \
            / qpoch(q**(2*l + 2*qq + 4), q**2, l + p - n) \
            * qpoch(q**(4*l + 4 - 2*n), q**2, n) \
            / qpoch(q**2, q**2, n)
    
    ret *= q**(2*l**2 + 2*l*p + 2*l*qq + 2*p*qq + p + qq) \
        * (1 - q**(4*l + 2)) / (1 - q**(2*l + 2*qq + 2)) \
        * qpoch(q**2, q**2, l - p) \
        * qpoch(q**2, q**2, l + p) \
        / qpoch(q**2, q**2, 2*l)

    return ret

def func_lhs(l, s, p, qq):
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        ret2 += bhs([
                q**(2*p - 2*l + 2*s),   
                q**(2*qq - 2*l),
                q**(2*l + 2*qq - 2*n + 2),
                q**(2*l + 2*qq - 2*n + 2)
            ], [
                q**(2*s + 2),
                q**(-2*l - 2*n + 2*qq),
                q**(-2*l - 2*n + 2*qq)
            ], q**2, q**2)
        ret2 *= q**(2*n)
        ret2 *= qpoch(q**(2 + 2*l - 2*qq), q**2, n)**2 \
            * qpoch(q**(-2*l - 2*qq + 2*s), q**2, n) \
            / qpoch(q**(-2*l - 2*qq), q**2, n) \
            / qpoch(q**(2 + 2*p - 2*qq + 2*s), q**2, n) \
            / qpoch(q**2, q**2, n)
        ret1 += ret2

    return ret1

def func_rhs1(l, s, p, qq):
    '''Move all undependant variables from left to right.'''
    ret = 0

    for n in range(l + qq - s + 1):
        ret += q**(n*(-4*l - 2*p - 2*qq - 2) + 2*n**2) \
            * qpoch(q**(2*l - 2*qq), q**(-2), l + p - n) \
            / qpoch(q**(2*l + 2*qq + 4), q**2, l + p - n) \
            * qpoch(q**(4*l + 4 - 2*n), q**2, n) \
            / qpoch(q**2, q**2, n)
    
    ret *= q**(2*l**2 + 2*l*p + 2*l*qq + 2*p*qq + p + qq) \
        * (1 - q**(4*l + 2)) / (1 - q**(2*l + 2*qq + 2)) \
        * qpoch(q**2, q**2, l - p) \
        * qpoch(q**2, q**2, l + p) \
        / qpoch(q**2, q**2, 2*l)
    
    ret *= 1 / qpoch(q**(2*p - 2*l), q**2, s) \
        *  qpoch(q**2, q**2, s) \
        / qpoch(q**2, q**2, l + qq)**2 \
        * qpoch(q**(-4*l), q**2, l + qq)**2 \
        / qpoch(q**(-2*l - 2*p), q**2, l + qq - s) \
        * qpoch(q**2, q**2, l + qq - s) \
        / qpoch(q**2, q**2, l + qq) \
        * qpoch(q**2, q**2, l + qq)

    ret *= (-1)**(l + qq) / q**(-5*l^2 + (2*l + 1)*p - 2*(2*l - p + 1)*qq + qq^2 + 2*(l - p + 1)*s - 3*l)
    
    return ret

def func_rhs2(l, s, p, qq):
    '''Simplifying.'''
    ret = 0

    for n in range(l + qq - s + 1):
        ret += q**(n*(-4*l - 2*p - 2*qq - 2) + 2*n**2) \
            * qpoch(q**(2*l - 2*qq), q**(-2), l + p - n) \
            / qpoch(q**(2*l + 2*qq + 4), q**2, l + p - n) \
            * qpoch(q**(4*l + 4 - 2*n), q**2, n) \
            / qpoch(q**2, q**2, n)
    
    ret *= q**(7*l**2 + 3*(2*l + 1)*qq - qq**2 - 2*(l - p + 1)*s + 3*l) \
        * (1 - q**(4*l + 2)) / (1 - q**(2*l + 2*qq + 2)) \
        * qpoch(q**2, q**2, l - p) \
        * qpoch(q**2, q**2, l + p) \
        * qpoch(q**2, q**2, s) \
        * qpoch(q**2, q**2, l + qq) \
        * qpoch(q**(-4*l), q**2, l + qq)**2 \
        * qpoch(q**2, q**2, l + qq - s) \
        / qpoch(q**2, q**2, l + qq)**2 \
        / qpoch(q**(-2*l - 2*p), q**2, l + qq - s) \
        / qpoch(q**2, q**2, l + qq) \
        / qpoch(q**2, q**2, 2*l) \
        / qpoch(q**(2*p - 2*l), q**2, s)

    ret *= (-1)**(l + qq)
    
    return ret

def func_rhs3(l, s, p, qq):
    '''Rewriting everything to (q^2)_n.'''
    ret = 0

    for n in range(l + qq - s + 1):
        ret += q**(n*(-4*l - 2*p - 2*qq - 2) + 2*n**2) \
            * qpoch(q**(2*l - 2*qq), q**(-2), l + p - n) \
            / qpoch(q**(2*l + 2*qq + 4), q**2, l + p - n) \
            * qpoch(q**(4*l + 4 - 2*n), q**2, n) \
            / qpoch(q**2, q**2, n)
    
    ret *= q**(7*l**2 + 3*(2*l + 1)*qq - qq**2 - 2*(l - p + 1)*s + 3*l + 4*binomial(l + qq, 2) - 4*(2*l)*(l + qq) - 2*binomial(l + qq - s, 2) + 2*(l + p)*(l + qq - s) - 2*binomial(s, 2) + 2*(l - p)*s) \
        * (1 - q**(4*l + 2)) / (1 - q**(2*l + 2*qq + 2)) \
        * qpoch(q**2, q**2, l - p) \
        * qpoch(q**2, q**2, l + p) \
        * qpoch(q**2, q**2, s) \
        * qpoch(q**2, q**2, l + qq) \
        * qpoch(q**2, q**2, 2*l)**2 \
        * qpoch(q**2, q**2, l - p - s) \
        * qpoch(q**2, q**2, p - qq + s) \
        * qpoch(q**2, q**2, l + qq - s) \
        / qpoch(q**2, q**2, l - qq)**2 \
        / qpoch(q**2, q**2, l + qq)**2 \
        / qpoch(q**2, q**2, l + p) \
        / qpoch(q**2, q**2, l + qq) \
        / qpoch(q**2, q**2, 2*l) \
        / qpoch(q**2, q**2, l - p) \
    
    return ret

def func_rhs4(l, s, p, qq):
    '''Simplifying q.'''
    ret = 0

    for n in range(l + qq - s + 1):
        ret += q**(n*(-4*l - 2*p - 2*qq - 2) + 2*n**2) \
            * qpoch(q**(2*l - 2*qq), q**(-2), l + p - n) \
            / qpoch(q**(2*l + 2*qq + 4), q**2, l + p - n) \
            * qpoch(q**(4*l + 4 - 2*n), q**2, n) \
            / qpoch(q**2, q**2, n)
    
    ret *= q**(2*(l + p + s + 1)*(l + qq - s)) \
        * qpoch(q**2, q**2, s) \
        * qpoch(q**2, q**2, 2*l + 1) \
        * qpoch(q**2, q**2, l - p - s) \
        * qpoch(q**2, q**2, p - qq + s) \
        * qpoch(q**2, q**2, l + qq - s) \
        / qpoch(q**2, q**2, l - qq)**2 \
        / qpoch(q**2, q**2, l + qq) \
        / qpoch(q**2, q**2, l + qq + 1)
    
    return ret

def func_rhs5(l, s, p, qq):
    '''Replace to qbinomials.'''
    ret = 0

    for n in range(l + qq - s + 1):
        ret += q**(n*(-4*l - 2*p - 2*qq - 2) + 2*n**2) \
            * qpoch(q**(2*l - 2*qq), q**(-2), l + p - n) \
            / qpoch(q**(2*l + 2*qq + 4), q**2, l + p - n) \
            * qpoch(q**(4*l + 4 - 2*n), q**2, n) \
            / qpoch(q**2, q**2, n)
    
    ret *= q**(2*(l + p + s + 1)*(l + qq - s)) \
        * qbinomial(2*l + 1, l - qq, q**2) \
        / qbinomial(l - qq, l - p - s, q**2) \
        / qbinomial(l + qq, s, q**2)
 
    return ret

def func_rhs6(l, s, p, qq):
    '''Simplifying the n-sum.'''
    ret = 0

    for n in range(l + qq - s + 1):
        ret += (-1)**(l + p - n) \
            * q**( l**2 - (4*l + 3)*n + n**2 - p**2 - 2*(l + p)*qq + l + p  ) \
            * qpoch(q**(2*qq - 2*l), q**2, l + p - n) \
            * qpoch(q**(4*l + 4 - 2*n), q**2, n) \
            / qpoch(q**(2*l + 2*qq + 4), q**2, l + p - n) \
            / qpoch(q**2, q**2, n)
    
    ret *= q**(2*(l + p + s + 1)*(l + qq - s)) \
        * qbinomial(2*l + 1, l - qq, q**2) \
        / qbinomial(l - qq, l - p - s, q**2) \
        / qbinomial(l + qq, s, q**2)
 
    return ret

def func_rhs65(l, s, p, qq):
    '''Taking q-powers together.'''
    ret = 0

    for n in range(l + qq - s + 1):
        ret += (-1)**(l + p - n) \
            * q**( n**2 - n*(4*l + 3)  ) \
            * qpoch(q**(2*qq - 2*l), q**2, l + p - n) \
            * qpoch(q**(4*l + 4 - 2*n), q**2, n) \
            / qpoch(q**(2*l + 2*qq + 4), q**2, l + p - n) \
            / qpoch(q**2, q**2, n)
    
    ret *= q**( 3*l**2 + (2*l + 1)*p - p**2 - 2*(p - qq + 1)*s - 2*s^2 + 3*l + 2*qq  ) \
        * qbinomial(2*l + 1, l - qq, q**2) \
        / qbinomial(l - qq, l - p - s, q**2) \
        / qbinomial(l + qq, s, q**2)
 
    return ret

def func_rhs7(l, s, p, qq):
    '''Simplifying the n-sum even more.'''
    ret = 0

    for n in range(l + qq - s + 1):
        ret += (-1)**(l + p) \
            * qpoch(q**(2*qq + 2*p), q**2, -n) \
            * qpoch(q**(-4*l - 2), q**2, n) \
            / qpoch(q**(4*l + 2*p + 2*qq + 4), q**2, -n) \
            / qpoch(q**2, q**2, n)
    
    ret *= q**(2*(l + p + s + 1)*(l + qq - s) + l**2 - p**2 - 2*(l + p)*qq + l + p) \
        * qpoch(q**(2*qq - 2*l), q**2, l + p) \
        / qpoch(q**(2*l + 2*qq + 4), q**2, l + p) \
        * qbinomial(2*l + 1, l - qq, q**2) \
        / qbinomial(l - qq, l - p - s, q**2) \
        / qbinomial(l + qq, s, q**2)
 
    return ret


# Proving a basic hypergeometric equation

def test_hyper(f1, f2):
    for N in range(1, 8):
        for b in range(N):
            for a in range(N-b, b+1):
                for c in range(N - b):
                    if a <= N and b <= N and a <= b and N <= a + b and c <= N - b:
                        print '%i, %i, %i, %i' % (N, a, b, c)
                        try:
                            if not bool(f1(N, a, b, c) == f2(N, a, b, c)):
                                return False
                        except ZeroDivisionError:
                            print 'division by zero'
    
    # Some random tests without zero divisors
    #tests = [(3, 2, 0, 0), (2, 0, 0, 1), (3, 1, 0, 1), (4, 1, 0, 0), (4, 1, 0, 1), (4, 1, 2, 1)]
    #for N, b, a, c in tests:
    #    print '%i, %i, %i, %i' % (N, a, b, c) 
    #    if not bool(f1(N, a, b, c) == f2(N, a, b, c)):
    #        return False
        

    return True

def hyper1(N, a, b, c):
    ret1 = 0

    for n in range(c+1):
        ret1 += qpoch([q**(-c), q**(b + 1), q**(b + 1)], q, n) \
            / qpoch([q, q**(N - a - c + 1), q**(b - N)], q, n) \
            * bhs([
                q**(-b), q**(N - a - b - c), q**(N - b - n + 1), q**(N - b-n+1)
            ], [
                q**(N - b - c + 1), q**(-b - n), q**(-b - n)
            ], q, q)

    return ret1

def hyper2(N, a, b, c):
    ret1 = 0

    for n in range(c+1):
        ret2 = 0
        for k in range(b + 1):

            ret3 = 0
            for j in range(b + 1):
                ret3 += qpoch([q**(-k), q**(-N-1)], q, j) \
                    / qpoch([q, q**(-b - n)], q, j) \
                    * q**((N - b - n + k + 1)*j)

            ret2 += qpoch([q**(-b), q**(N - a - b - c), q**(N - b - n + 1)], q, k) \
                / qpoch([q, q**(N - b - c + 1), q**(-b - n)], q, k) \
                * q**k \
                * ret3
        
        ret1 += qpoch([q**(-c), q**(b + 1), q**(b + 1)], q, n) \
            / qpoch([q, q**(N - a - c + 1), q**(b - N)], q, n) \
            * ret2

    return ret1

def hyper3(N, a, b, c):
    ret1 = 0

    for n in range(c+1):
        ret2 = 0

        for j in range(b + 1):
            for k in range(b - j + 1):
                ret2 += qpoch([q**(-b), q**(N - a - b - c), q**(N - b - n + 1)], q, k+j) \
                    / qpoch([q**(N - b - c + 1), q**(-b - n)], q, k + j) \
                    * q**(k + j) \
                    * qpoch(q**(-N - 1), q, j) \
                    / qpoch([q, q**(-b - n)], q, j) \
                    / qpoch(q, q, k) \
                    * (-1)**j \
                    * q**(binomial(j, 2) - j*(k + j) + (N - b - n + k + j + 1)*j)
 
        ret1 += qpoch([q**(-c), q**(b + 1), q**(b + 1)], q, n) \
            / qpoch([q, q**(N - a - c + 1), q**(b - N)], q, n) \
            * ret2

    return ret1

def hyper4(N, a, b, c):
    '''Rewrite the k-sum to a 3_phi_2'''
    ret1 = 0

    for n in range(c+1):
        ret2 = 0

        for j in range(b + 1):
            ret2 += (-1)**j * q**((N - b - n + 3/2)*j + j**2/2) \
                * qpoch([q**(-b), q**(N - a - b - c), q**(N - b - n + 1)], q, j) \
                / qpoch([q**(N - b - c + 1), q**(-b - n)], q, j) \
                * qpoch(q**(-N - 1), q, j) \
                / qpoch([q, q**(-b - n)], q, j) \
                * bhs(
                    [q**(j - b), q**(j + N - a - b - c), q**(j + N - b - n + 1)],
                    [q**(j + N - b - c + 1), q**(j - b - n)], q, q)
 
        ret1 += qpoch([q**(-c), q**(b + 1), q**(b + 1)], q, n) \
            / qpoch([q, q**(N - a - c + 1), q**(b - N)], q, n) \
            * ret2

    return ret1

def hyper5(N, a, b, c):
    '''Apply transformation III.11 on the 3_phi_2'''
    ret1 = 0

    for n in range(c+1):
        ret2 = 0

        for j in range(b + 1):
            ret2 += (-1)**j * q**((N - b - n + 3/2)*j + j**2/2) \
                * qpoch([q**(-b), q**(N - a - b - c), q**(N - b - n + 1)], q, j) \
                / qpoch([q**(N - b - c + 1), q**(-b - n)], q, j) \
                * qpoch(q**(-N - 1), q, j) \
                / qpoch([q, q**(-b - n)], q, j) \
                * q**((N - a - b + j - n)*(b - j)) \
                * qpoch(q**(a - N), q, b - j) \
                / qpoch(q**(j - b - n), q, b - j) \
                * bhs(
                    [q**(j - b), q**(a + 1), q**(n - c)],
                    [q**(j + N - b - c + 1), q**(a - N)], q, q) 
 
        ret1 += qpoch([q**(-c), q**(b + 1), q**(b + 1)], q, n) \
            / qpoch([q, q**(N - a - c + 1), q**(b - N)], q, n) \
            * ret2

    return ret1

def hyper6(N, a, b, c):
    '''Take everything together.'''
    ret1 = 0

    for n in range(c+1):
        for j in range(b + 1):
            for k in range(b - j + 1):
                ret1 += (-1)**j * q**((N - b - n + 3/2)*j + j**2/2 + k + (N - a - b + j - n)*(b - j)) \
                    * qpoch([q**(-b), q**(N - a - b - c), q**(N - b - n + 1), q**(-N - 1)], q, j) \
                    / qpoch([q, q**(N - b - c + 1), q**(-b - n), q**(-b - n)], q, j) \
                    * qpoch(q**(a - N), q, b - j) \
                    / qpoch(q**(j - b - n), q, b - j) \
                    * qpoch([q**(j - b), q**(a + 1), q**(n - c)], q, k) \
                    / qpoch([q, q**(j + N - b - c + 1), q**(a - N)], q, k) \
                    * qpoch([q**(-c), q**(b + 1), q**(b + 1)], q, n) \
                    / qpoch([q, q**(N - a - c + 1), q**(b - N)], q, n) \

    return ret1

def hyper7(N, a, b, c):
    '''Take all n-sums together.'''
    ret1 = 0

    for j in range(b + 1):
        for k in range(b - j + 1):
            ret2 = 0
            for n in range(c+1):
                ret2 += q**(-b*n) \
                    * qpoch([q**(N - b - n + 1) ], q, j) \
                    / qpoch([q**(-b - n), q**(-b - n)], q, j) \
                    / qpoch(q**(j - b - n), q, b - j) \
                    * qpoch(q**(n - c), q, k) \
                    * qpoch([q**(-c), q**(b + 1), q**(b + 1)], q, n) \
                    / qpoch([q, q**(N - a - c + 1), q**(b - N)], q, n) 
            ret1 += (-1)**j * ret2 \
                * q**(N*b - a*b - b^2 + 1/2*(2*a + 2*b + 3)*j - 1/2*j**2 + k) \
                * qpoch([q**(-b), q**(N - a - b - c), q**(-N-1)], q, j) \
                / qpoch([q, q**(N - b - c + 1)], q, j) \
                * qpoch(q**(a - N), q, b - j) \
                * qpoch([q**(j - b), q**(a + 1)], q, k) \
                / qpoch([q, q**(j + N - b - c + 1), q**(a - N)], q, k) \
        
    return ret1

def hyper8(N, a, b, c):
    '''All qpochhammers to n.'''
    ret1 = 0

    for j in range(b + 1):
        for k in range(b - j + 1):
            ret2 = 0
            for n in range(c+1):
                ret2 += q**(-b*n) \
                    * qpoch(q**(N - b + 1), q, j) * qpoch(q**(b - N), q, n) / qpoch(q**(b - N - j), q, n) * q**(-n*j) \
                    / qpoch([q**(-b - n), q**(-b - n)], q, j) \
                    * qpoch(q, q, n) / qpoch(q**(j - b), q, b - j) / qpoch(q**(1 + b - j), q, n) * q**(n*(b - j)) \
                    * qpoch(q**(-c), q, k) * qpoch(q**(k - c), q, n) / qpoch(q**(-c), q, n) \
                    * qpoch([q**(-c), q**(b + 1), q**(b + 1)], q, n) \
                    / qpoch([q, q**(N - a - c + 1), q**(b - N)], q, n) 
            ret1 += (-1)**j * ret2 \
                * q**(N*b - a*b - b^2 + 1/2*(2*a + 2*b + 3)*j - 1/2*j**2 + k) \
                * qpoch([q**(-b), q**(N - a - b - c), q**(-N-1)], q, j) \
                / qpoch([q, q**(N - b - c + 1)], q, j) \
                * qpoch(q**(a - N), q, b - j) \
                * qpoch([q**(j - b), q**(a + 1)], q, k) \
                / qpoch([q, q**(j + N - b - c + 1), q**(a - N)], q, k) \
        
    return ret1

def hyper9(N, a, b, c):
    '''Simplifying a lot.'''
    ret1 = 0

    for j in range(b + 1):
        for k in range(b - j + 1):
            ret2 = 0
            for n in range(c-k+1):
                ret2 += qpoch([q**(k - c), q**(1 - j + b)], q, n) \
                    / qpoch([q**(b - N - j), q**(N - a - c + 1)], q, n) 
            ret1 += (-1)**j * ret2 \
                * qpoch(q**(N - b + 1), q, j) \
                / qpoch(q**(-b), q, j)  \
                / qpoch(q**(- b), q, b) \
                * qpoch(q**(-c), q, k) \
                * q**(N*b - a*b - b^2 + 1/2*(2*a + 2*b + 3)*j - 1/2*j**2 + k) \
                * qpoch([q**(-b), q**(N - a - b - c), q**(-N-1)], q, j) \
                / qpoch([q, q**(N - b - c + 1)], q, j) \
                * qpoch(q**(a - N), q, b - j) \
                * qpoch([q**(j - b), q**(a + 1)], q, k) \
                / qpoch([q, q**(j + N - b - c + 1), q**(a - N)], q, k) \
        
    return ret1

def funcA1(l, s, p, qq):
    '''Start functions for a q-analogue of proposition A1.'''
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        i, j = n + s - (l + qq)/2, (l + p)/2 - n
        for m in srange(0, l - s - p + 1):
            k = m + n + s - l
            ret2 += q**(2*m*(2*l - 2*n + p + qq - s + 1)) \
                * qbinomial(l - p, m + s, q**2) \
                * qbinomial(l - qq, m, q**2) \
                / qbinomial(2*l, m + n + s, q**2)**2
        ret2 *= q**(-2*n*(s + 1)) \
            * qbinomial(2*l1, n, q**2) \
            * qbinomial(2*m1, n + s, q**2) 
        ret1 += ret2

    return q**(2*(l + p)*s + p + qq)*ret1


def funcA2(l, s, p, qq):
    '''Substitute M = l - p - s - m.'''
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        i, j = n + s - (l + qq)/2, (l + p)/2 - n
        for M in srange(0, l - s - p + 1):
            ret2 += q**(M*(-4*l + 4*n - 2*p - 2*qq + 2*s - 2)) \
                * qbinomial(l - p, l - p - M, q**2) \
                * qbinomial(l - qq, l - p - s - M, q**2) \
                / qbinomial(2*l, l - p - M + n, q**2)**2
        ret2 *= q**(n*(-4*l + 4*p + 2*s - 2)) \
            * qbinomial(2*l1, n, q**2) \
            * qbinomial(2*m1, n + s, q**2) 
        ret1 += ret2

    return q**(2*(l + p)*s + p + qq + 2*(2*l + p + qq - s + 1)*(l - p - s))*ret1

def funcA3(l, s, p, qq):
    '''Rewriting M-sum to qpochhammers.'''
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        for M in srange(0, l - s - p + 1):
            ret2 += q**(M*(-4*l + 4*n - 2*p - 2*qq + 2*s - 2)) \
                * qpoch(q**2, q**2, l - p - M + n) \
                * qpoch(q**2, q**2, l - p - M + n) \
                * qpoch(q**2, q**2, l + p + M - n) \
                * qpoch(q**2, q**2, l + p + M - n) \
                / qpoch(q**2, q**2, M) \
                / qpoch(q**2, q**2, l - p - M) \
                / qpoch(q**2, q**2, l - p - s - M) \
                / qpoch(q**2, q**2, p - qq + s + M)
        ret2 *= q**(n*(-4*l + 4*p + 2*s - 2)) \
            * qpoch(q**2, q**2, l - p) \
            * qpoch(q**2, q**2, l - qq) \
            / qpoch(q**2, q**2, 2*l) \
            / qpoch(q**2, q**2, 2*l) \
            * qbinomial(2*l1, n, q**2) \
            * qbinomial(2*m1, n + s, q**2) 
        ret1 += ret2

    return q**(2*(l + p)*s + p + qq + 2*(2*l + p + qq - s + 1)*(l - p - s))*ret1

def funcA4(l, s, p, qq):
    '''Simplifying qpochhammers.'''
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        for M in srange(0, l - s - p + 1):
            B = qpoch(q**(2*l - 2*p - 2*M + 2), q**2, n) \
                * qpoch(q**(2*p - 2*qq + 2*s + 2*M + 2), q**2, l + qq - n - s)
            ret2 += q**(M*(-4*l + 2*n - 2*p - 2*qq - 2)) \
                * qpoch(q**2, q**2, l - p + n) \
                * qpoch(q**2, q**2, l + p - n) \
                / qpoch(q**2, q**2, l - p - s) \
                * qpoch(q**(2*l + 2*p - 2*n + 2), q**2, M) \
                * qpoch(q**(-2*l + 2*p + 2*s), q**2, M) \
                / qpoch(q**2, q**2, M) \
                / qpoch(q**(-2*l + 2*p - 2*n), q**2, M) \
                * B
        ret2 *= q**(n*(-4*l + 4*p + 2*s - 2)) \
            * qpoch(q**2, q**2, l - p) \
            * qpoch(q**2, q**2, l - qq) \
            / qpoch(q**2, q**2, 2*l) \
            / qpoch(q**2, q**2, 2*l) \
            * qbinomial(2*l1, n, q**2) \
            * qbinomial(2*m1, n + s, q**2) 
        ret1 += ret2

    return q**(2*(l + p)*s + p + qq + 2*(2*l + p + qq - s + 1)*(l - p - s))*ret1

def genA(upto):
    z = var('z')
    P = function('P', z)
    #Aret = [P(1)]

    #for k in range(1, upto+1):
    #    Asum = -P(q**(-2*k))
    #    for l in range(k):
    #        print k, l, (-1)**l * qpoch(q**(-2*k), q**2, l) * Aret[l]
    #        Asum += (-1)**l * qpoch(q**(-2*k), q**2, l) * Aret[l]
    #    Asum *= (-1)**k / qpoch(q**(-2*k), q**2, k)
    #    Aret.append(Asum)

    p0, p1, p2, p3, p4, p5, p6, p7, p8 = var('p0 p1 p2 p3 p4 p5 p6 p7 p8')
    a0 = p0
    a1 = 1 / qpoch(q**(-2), q**2, 1) * (a0 - p1)
    a2 = - 1 / qpoch(q**(-4), q**2, 2) * (a0 - a1*qpoch(q**(-4), q**2, 1) - p2)
    a3 = 1 / qpoch(q**(-6), q**2, 3) * (a0 - a1*qpoch(q**(-6), q**2, 1) + a2*qpoch(q**(-6), q**2, 2) - p3)
    a4 = - 1 / qpoch(q**(-8), q**2, 4) * (a0 - a1*qpoch(q**(-8), q**2, 1) + a2*qpoch(q**(-8), q**2, 2) - a3*qpoch(q**(-8), q**2, 3) - p4)

    return [a0, a1, a2, a3, a4]

def qdiff(func):
    return (lambda z : z**(-1) * (func(z) - func(z*q**(-2))))

def funcA5(l, s, p, qq):
    '''Simplifying qpochhammers.'''
    l1 = (l + p)/2
    l2 = (l - p)/2
    m1 = (l + qq)/2
    m2 = (l - qq)/2

    ret1 = 0
    for n in srange(l + qq - s + 1):
        ret2 = 0
        for M in srange(0, l - s - p + 1):
            B = lambda z : qpoch(q**(2*l - 2*p - 2*z + 2), q**2, n) \
                * qpoch(q**(2*p - 2*qq + 2*s + 2*z + 2), q**2, l + qq - n - s)
            
            A = B
            ret3 = A(1)
            for t in range(1, l + qq - s + 1):
                A = qdiff(A)
                ret3 += (-1)**t * A(1) / (q**(2*(t - 1)) * qpoch(q**(-2*t), q**2, t)) * qpoch(q**(-2*M), q**2, t)

            ret2 += q**(M*(-4*l + 2*n - 2*p - 2*qq - 2)) \
                * qpoch(q**2, q**2, l - p + n) \
                * qpoch(q**2, q**2, l + p - n) \
                / qpoch(q**2, q**2, l - p - s) \
                * qpoch(q**(2*l + 2*p - 2*n + 2), q**2, M) \
                * qpoch(q**(-2*l + 2*p + 2*s), q**2, M) \
                / qpoch(q**2, q**2, M) \
                / qpoch(q**(-2*l + 2*p - 2*n), q**2, M) \
                * ret3
        ret2 *= q**(n*(-4*l + 4*p + 2*s - 2)) \
            * qpoch(q**2, q**2, l - p) \
            * qpoch(q**2, q**2, l - qq) \
            / qpoch(q**2, q**2, 2*l) \
            / qpoch(q**2, q**2, 2*l) \
            * qbinomial(2*l1, n, q**2) \
            * qbinomial(2*m1, n + s, q**2) 
        ret1 += ret2

    return q**(2*(l + p)*s + p + qq + 2*(2*l + p + qq - s + 1)*(l - p - s))*ret1

