load('qpolynomials.sage')

q = var('q')
x = var('x')
a = var('a')

def scalardiv(M1, M2, M, N):
    newm = matrix(SR, M, N)

    for i in range(M):
        for j in range(N):
            if M1[i,j] != 0:
                newm[i,j] = (M1[i,j] / M2[i,j]).full_simplify().factor()

    return newm

def L(N):
    ret = matrix(SR, N, N)

    for m in range(N):
        for n in range(m+1):
            ret[m,n] = replace_z_x(continuous_q_jacobi(m - n, x, a + n, a + n, q))

    return ret

def Linv(N):
    ret = matrix(SR, N, N)

    for m in range(N):
        for n in range(m+1):
            ret[m, n] = replace_z_x(continuous_q_jacobi(m - n, x, a - n, a - n, q))

    return ret
