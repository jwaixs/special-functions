def SimplifyMatrix(mat, do_fac=True):
    def fsimplify(elm):
        elm = elm.full_simplify()
        if do_fac and elm != 0:
            elm = elm.factor()
        return elm

    return mat.apply_map(fsimplify)


def scalardiv(M1, M2, M, N):
    newm = matrix(SR, M, N)

    for i in range(M):
        for j in range(N):
            if M1[i,j] != 0:
                newm[i,j] = (M1[i,j] / M2[i,j]).full_simplify().factor()

    return newm

def antimatrix(F, N):
    newm = matrix(F, N, N)

    for i in range(N):
        newm[i, N-i-1] = 1

    return newm
