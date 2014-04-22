from copy import deepcopy

def lpq(M, pp, qq):
    n = len(M)
    if n >= pp >= 1 and n >= qq >= 1:
        return M[n - qq][pp-1] - pp
    return None


def ajp2(M, pp, jj):
    prod1 = 1
    for i in range(1, pp+2):
        prod1 *= lpq(M, i, pp+1) - lpq(M, jj, pp)
   
    prod2 = 1
    for i in range(1, pp):
        prod2 *= lpq(M, i, pp-1) - lpq(M, jj, pp) - 1

    prod3 = 1
    for i in range(1, pp+1):
        if i != jj:
            prod3 *= (lpq(M, i, pp) - lpq(M, jj, pp)) \
                * (lpq(M, i, pp) - lpq(M, jj, pp) - 1) 
    
    return (abs(prod1*prod2/prod3))**(1/2)


def bjp2(M, pp, jj):
    prod1 = 1
    for i in range(1, pp+2):
        prod1 *= lpq(M, i, pp+1) - lpq(M, jj, pp) + 1
   
    prod2 = 1
    for i in range(1, pp):
        prod2 *= lpq(M, i, pp-1) - lpq(M, jj, pp)

    prod3 = 1
    for i in range(1, pp+1):
        if i != jj:
            prod3 *= (lpq(M, i, pp) - lpq(M, jj, pp)) \
                * (lpq(M, i, pp) - lpq(M, jj, pp) + 1) 
    
    return (abs(prod1*prod2/prod3))**(1/2)


def mshift(N, pp, jj):
    N = deepcopy(N)
    n = len(N)
    N[n - pp][abs(jj) - 1] += sign(jj)
    return N


def getm(M, ii, jj):
    n = len(M)
    return M[n - jj][ii - 1]


def bjp(M, kk, jj):
    prod1 = 1
    for i in range(1, kk+2):
        prod1 *= getm(M, i, kk+1) - getm(M, jj, kk) - i + jj + 1

    prod2 = 1
    for i in range(1, kk):
        prod2 *= getm(M, i, kk-1) - getm(M, jj, kk) - i + jj

    prod3 = 1
    for i in range(1, kk+1):
        if i != jj:
            prod3 *= (getm(M, i, kk) - getm(M, jj, kk) - i + jj + 1) \
                * (getm(M, i, kk) - getm(M, jj, kk) - i + jj)

    return (-prod1 * prod2 / prod3)**(1/2)


def ajp(M, kk, jj):
    prod1 = 1
    for i in range(1, kk+2):
        prod1 *= getm(M, i, kk+1) - getm(M, jj, kk) - i + jj

    prod2 = 1
    for i in range(1, kk):
        prod2 *= getm(M, i, kk-1) - getm(M, jj, kk) - i + jj - 1

    prod3 = 1
    for i in range(1, kk+1):
        if i != jj:
            prod3 *= (getm(M, i, kk) - getm(M, jj, kk) - i + jj) \
                * (getm(M, i, kk) - getm(M, jj, kk) - i + jj - 1)

    return (-prod1 * prod2 / prod3)**(1/2)


m11, m12, m22, m13, m23, m33 = var('m11 m12 m22 m13 m23 m33')

assume(m13, 'integer')
assume(m23, 'integer')
assume(m33, 'integer')
assume(m12, 'integer')
assume(m22, 'integer')
assume(m11, 'integer')
assume(m13 >= m12)
assume(m12 >= m23)
assume(m23 >= m22)
assume(m22 >= m33)
assume(m12 >= m11)
assume(m11 >= m22)

m, l, r = var('m l r')
m13 = m + l
m23 = m + r
m33 = -r
m12 = m + l
m22 = m

M = [[m13, m23, m33], [m12, m22], [m11]]
test_vec = [ (1, M) ]

def A(m11, m12, m22, m13, m23, m33):
    nom = (m13 - m12) * (m12 - m23 + 1) * (m12 - m33 + 2) * (m11 - m22 + 1)
    denom = (m12 - m22 + 1) * (m12 - m22 + 2)
    return (nom / denom)**(1/2)

def B(m11, m12, m22, m13, m23, m33):
    nom = (m13 - m22 + 1) * (m23 - m22) * (m22 - m33 + 1) * (m12 - m11)
    denom = (m12 - m22 + 1) * (m12 - m22)
    return - (nom / denom)**(1/2)

def C(m11, m12, m22, m13, m23, m33):
    nom = (m11 - m22) * (m13 - m12 + 1) * (m12 - m23) * (m12 - m33 + 1)
    denom = (m12 - m22) * (m12 - m22 + 1)
    return (nom / denom)**(1/2)

def D(m11, m12, m22, m13, m23, m33):
    nom = (m12 - m11 + 1) * (m13 - m22 + 2) * (m23 - m22 + 1) * (m22 - m33)
    denom = (m12 - m22 + 1) * (m12 - m22 + 2)
    return -(nom / denom)**(1/2)

def total_func2():
    p1 = C(m11+1, m12+1, m22, m13, m23, m33) * A(m11, m12, m22, m13, m23, m33)
    p2 = D(m11+1, m12, m22+1, m13, m23, m33) * B(m11, m12, m22, m13, m23, m33)
    p3 = A(m11-1, m12-1, m22, m13, m23, m33) * C(m11, m12, m22, m13, m23, m33)
    p4 = B(m11-1, m12, m22-1, m13, m23, m33) * D(m11, m12, m22, m13, m23, m33)

    return -(p1 + p2 + p3 + p4)


def try_full_simplify(exp):
    try:
        return exp.full_simplify().factor()
    except:
        return exp


def simplify_vector(vec):
    ret = []
    done = []

    for s1, v1 in vec:
        if not v1 in done:
            scalar = 0
            for s2, v2 in vec:
                if v2 == v1:
                    scalar += s2
            ret.append((try_full_simplify(scalar), v1))
            done.append(v1)

    return ret


def E11(vec):
    ret = []

    for s, v in vec:
        ret.append((s * getm(v, 1, 1), v))

    return ret

def E22(vec):
    ret = []

    for s, v in vec:
        ret.append((s * (getm(v, 1, 2) + getm(v, 2, 2) - getm(v, 1, 1)), v))

    return ret

def E33(vec):
    ret = []

    for s, v in vec:
        scalar = getm(v, 1, 3) + getm(v, 2, 3) + getm(v, 3, 3) \
            - getm(v, 1, 2) - getm(v, 2, 2)
        ret.append((s * scalar, v))

    return ret

def E12(vec):
    ret = []

    for s, v in vec:
        scalar = s * ajp(v, 1, 1)
        ret.append((scalar, mshift(v, 1, 1)))

    return simplify_vector(ret)

def E21(vec):
    ret = []

    for s, v in vec:
        scalar = s * bjp(v, 1, 1)
        ret.append((scalar, mshift(v, 1, -1)))

    return ret

def E23(vec):
    ret = []

    for s, v in vec:
        s1 = s * ajp(v, 2, 1)
        s2 = s * ajp(v, 2, 2)
        ret.append((s1, mshift(v, 2, 1)))
        ret.append((s2, mshift(v, 2, 2)))

    return simplify_vector(ret)

def E32(vec):
    ret = []

    for s, v in vec:
        s1 = s * bjp(v, 2, 1)
        s2 = s * bjp(v, 2, 2)
        ret.append((s1, mshift(v, 2, -1)))
        ret.append((s2, mshift(v, 2, -2)))

    return simplify_vector(ret)

def E13(vec):
    ret = []

    for s, v in vec:
        s1 = s * (ajp(v, 2, 1) * ajp(mshift(v, 2, 1), 1, 1) \
            - ajp(v, 1, 1) * ajp(mshift(v, 1, 1), 2, 1))
        s2 = s * (ajp(v, 2, 2) * ajp(mshift(v, 2, 2), 1, 1) \
            - ajp(v, 1, 1) * ajp(mshift(v, 1, 1), 2, 2))
        ret.append((s1, mshift(mshift(v, 1, 1), 2, 1)))
        ret.append((s2, mshift(mshift(v, 1, 1), 2, 2)))

    return simplify_vector(ret)

def E31(vec):
    ret = []

    for s, v in vec:
        s1 = s * (bjp(v, 1, 1) * bjp(mshift(v, 1, -1), 2, 1) \
            - bjp(v, 2, 1) * bjp(mshift(v, 2, -1), 1, 1))
        s2 = s * (bjp(v, 1, 1) * bjp(mshift(v, 1, -1), 2, 2) \
            - bjp(v, 2, 2) * bjp(mshift(v, 2, -2), 1, 1))
        ret.append((s1, mshift(mshift(v, 1, -1), 2, -1)))
        ret.append((s2, mshift(mshift(v, 1, -1), 2, -2)))

    return simplify_vector(ret)

def E13E31(vec):
    ret = []

    for s, v in E13(vec):
        ret.append((s, v))

    for s, v in E31(vec):
        ret.append((-s, v))

    return simplify_vector(ret)

def E13E31_only_vectors(vec):
    ret = []

    for v in vec:
        m1 = mshift(mshift(v, 1, 1), 2, 1)
        if not m1 in ret:
            ret.append(m1)

        m2 = mshift(mshift(v, 1, 1), 2, 2)
        if not m2 in ret:
            ret.append(m2)

        m3 = mshift(mshift(v, 1, -1), 2, -1)
        if not m3 in ret:
            ret.append(m3)

        m4 = mshift(mshift(v, 1, -1), 2, -2)
        if not m4 in ret:
            ret.append(m4)


    return ret
