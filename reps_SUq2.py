#from qpolynomials import little_q_Jacobi_polynomials
#from qfunctions import qbinomial

load('qpolynomials.py')
load('qfunctions.py')

#F.<q> = FunctionField(AA)
#A.<a,b,c,d> = FreeAlgebra(F, 4)

def clnm(l, n, m, q):
    print l, n, m
    return qbinomial(l+n, n-m, q**2)
    return qbinomial(l-m, n-m, q**2)**(1/2) \
        * qbinomial(l+n, n-m, q**2)**(1/2) \
        * q**(-(n-m)*(l-n))

def reps(l, n, m, q):
    global a, b, c, d


    k = max(abs(m), abs(n))

    if n >= m and m >= -n:
        p = little_q_Jacobi_polynomials(l-k, -q**(-1)*b*c, q**(2*abs(n-m)), \
            q**(2*abs(n+m)), q**2)
        return clnm(l, n, m, q)*d**(n+m)*c**(n-m)*p
    elif m >= n and n >= -m:
        p = little_q_Jacobi_polynomials(l-k, -q**(-1)*b*c, q**(2*abs(n-m)), \
            q**(2*abs(n+m)), q**2)
        return clnm(l, m, n, q)*d**(n+m)*b**(m-n)*p
    elif -n >= m and m >= n:
        p = little_q_Jacobi_polynomials(l-k, -q**(2*m+2*n-1)*b*c, \
            q**(2*abs(n-m)), q**(2*abs(n+m)), q**2)
        print clnm(l, -n, -m, q)
        return clnm(l, -n, -m, q)*b**(m-n)*a**(-m-n)*p
    elif -m >= n and n >= m:
        p = little_q_Jacobi_polynomials(l-k, -q**(2*m+2*n-1)*b*c, \
            q**(2*abs(n-m)), q**(2*abs(n+m)), q**2)
        return clnm(l, -m, -n, q)*c**(n-m)*a**(-m-n)*p
        
        
