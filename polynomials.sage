from hypergeometric import poch

def jacobi(n, x, a, b):
    ret = 0
    for k in range(n+1):
        ret += poch(-n, k) * poch(n + a + b + 1, k) * poch(a + k + 1, n - k) \
            / poch(1, k) * ((1 - x)/2)**k

    return ret / poch(1, n)

def test_ck(n, m, a, b):
    ret = 0

    for k in range(n, m+1):
        ret += jacobi(m-k, x, a+k, b+k) \
            * (a - b)/(n + a) * jacobi(k - n, x, -a - k, - b - k - 1)
            #((k+b)/(n+a) * jacobi(k-n, x, -a-k, -b-k) \
            #    + (a-b)/(n+a) * jacobi(k-n, x, -a-k, -b-k-1))

    return ret
