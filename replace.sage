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
