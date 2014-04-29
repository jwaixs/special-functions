'''complextools.sage - Created by Noud Aldenhoven'''

def replace_z_x(polynomial, z, x):
    '''replace_z_x(polynomial, z, x): substitute (z + z^(-1))/2 |-> x in polynomial.'''

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


def replace_x_z(polynomial, x, z):
    '''replace_x_z(polynomial, x, z): substitute x |-> (z + z^(-1))/2 in polynomial.'''

    return polynomial.substitute({x : (z + z**(-1))/2})

def askey_wilson_op(func, z):
    '''askey_wilson_op(func, z): apply the Askey-Wilson operator on a function in z.'''

    return (func.substitute({z : q**(1/2)*z}) - func.substitute({z : q**(-1/2)*z})) \
        * 1/(q**(1/2) - q**(-1/2)) * 2/(z - z**(-1))

def average_op(func, z):
    '''average_op(func, z): apply the averaging operator on the function in z.'''

    return 1/2*(func.substitute({z : q**(1/2)*z}) + func.substitute({z : q**(-1/2)*z}))
