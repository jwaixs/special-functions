import itertools

m, n, t, l = var('m n t l')

coefs = [ a1*a2 for (a1, a2) in itertools.combinations_with_replacement([m, n, t, l], 2) ] 
a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21 = var('a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 a13 a14 a15 a16 a17 a18 a19 a20 a21')
free_param = [a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21]

poly = a1*m^2 + a10*l^2 + a2*m*n + a3*m*t + a4*l*m + a5*n^2 + a6*n*t + a7*l*n + a8*t^2 + a9*l*t + a12*m*n*t + a13*l*m*n + a14*l*m*t + a15*l*n*t + a16*m*n*l*t + a17*m + a18*n + a19*t + a20*l + a21

q = var('q')
equations = [
    poly.substitute({l : 0, m : 0, n : 0, t : 0}) == 1,
    poly.substitute({l : 1/2, m : 0, n : 0, t : 0}) == q**(-1),
    poly.substitute({l : 1/2, m : 1, n : 1, t : 1}) == q**(-1),
    poly.substitute({l : 1/2, m : 0, n : 1, t : 0}) == 1,
    poly.substitute({l : 1, m : 0, n : 0, t : 0}) == q**(-2),
    poly.substitute({l : 1, m : 1, n : 1, t : 1}) == q**(-2),
    poly.substitute({l : 1, m : 2, n : 2, t : 2}) == q**(-2),
    poly.substitute({l : 1, m : 0, n : 1, t : 0}) == q**(-1),
    poly.substitute({l : 1, m : 1, n : 2, t : 1}) == q**(-1),
    poly.substitute({l : 1, m : 0, n : 2, t : 0}) == 1,
    poly.substitute({l : 1, m : 1, n : 1, t : 0}) == q**2,
    poly.substitute({l : 3/2, m : 0, n : 0, t : 0}) == q**(-3),
    poly.substitute({l : 3/2, m : 1, n : 1, t : 1}) == q**(-3),
    poly.substitute({l : 3/2, m : 2, n : 2, t : 2}) == q**(-3),
    poly.substitute({l : 3/2, m : 3, n : 3, t : 3}) == q**(-3),
    poly.substitute({l : 3/2, m : 0, n : 1, t : 0}) == q**(-2),
    poly.substitute({l : 3/2, m : 1, n : 2, t : 1}) == q**(-2),
    poly.substitute({l : 3/2, m : 2, n : 3, t : 2}) == q**(-2),
    poly.substitute({l : 3/2, m : 0, n : 2, t : 0}) == q**(-1),
    poly.substitute({l : 3/2, m : 1, n : 3, t : 1}) == q**(-1),
    poly.substitute({l : 3/2, m : 1, n : 1, t : 0}) == q,
    poly.substitute({l : 3/2, m : 2, n : 2, t : 1}) == q,
]

def used_params(fp):
    ret = []
    
    for elm in fp:
        for v in elm.variables():
            if v != q and not v in ret:
                ret.append(v)

    return ret

print solve(equations, used_params(equations))
