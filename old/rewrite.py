"""
rewrite.py - A simple rewrite module for monomials

A very simple script to do basic formula rewrite actions on (non-commutative)
monomials over an arbitrary field. As far as I know there is no other module in 
sage to do this.

AUTHORS::
    - Noud Aldenhoven, initial file (2013)

EXAMPLES::

    Here we implement an algebraic structure for the vectorspace of the
    quantised algebra of polynomials on SU(2)::

    sage: A = FunctionField(AA, 'q')
        sage: q = A.gen()
        sage: B = FreeAlgebra(A, 4, 'abcd')
        sage: a, b, c, d = B.gens()
        sage: rewrite_rules = {b*d: q*d*b, b*c: c*b, c*d: q*d*c, d*a: 1 + 1/q*b*c, a*b: q*b*a, a*c: q*c*a, a*d: 1 + q*b*c}
        sage: rewrite_formula(a**2*d*c, rewrite_rules)
        q*c*a + q^4*c^2*b*a
        sage: rewrite_formula(b*d, rewrite_rules)     
        q*d*b
        sage: rewrite_formula(b*d + b*c + a*d, rewrite_rules)
        1 + (q+1)*c*b + q*d*b
        sage: rewrite_formula(1/(q+2)*b*d + q**(-1)*a*b, rewrite_rules)
        b*a + (q/(q+2))*d*b

"""

import re
from sage.gsl.interpolation import spline
from string import join
from sage.calculus.functional import expand

def expand_term(term):
    old_term = re.split('(\\^|\\*)', term)

    new_term = []
    idx = 0
    
    while idx < len(old_term):
        if idx < len(old_term) - 2 and old_term[idx + 1] == '^':
            for i in xrange(int(old_term[idx + 2]) - 1):
                new_term.append(old_term[idx])
                new_term.append('*')
            new_term.append(old_term[idx])
            idx += 3            
        else:
            new_term.append(old_term[idx])
            idx += 1

    return join(new_term, '')
                
                

def find_replace_term(term, old, new, count=1):
    expanded_term = expand_term(term)
    expanded_old  = expand_term(old)
    expanded_new  = expand(new)

    return expanded_term.replace(expanded_old, \
        '(' + expanded_new + ')', count)


def rewrite_monomial(monomial, rewrite):
    for key in rewrite:
        new_monomial = find_replace_term(monomial, str(key), str(rewrite[key]))
        new_monomial = str(eval(new_monomial.replace('^', '**')))
        if new_monomial != monomial:
            return new_monomial

    return monomial

def reverse_spline(formula):
    newformula = ''

    for (coefficient, monomial) in formula:
        newformula += '+(' + str(coefficient) + ')*(' + str(monomial) + ')'

    newformula = eval(newformula.replace('^', '**'))

    return newformula

def rewrite_monomials(formula, rewrite):
    newformula = map(\
        lambda (c, m): (c, rewrite_monomial(str(m), rewrite)), formula)
    newformula = reverse_spline(newformula)

    return spline(newformula)


def rewrite_formula(formula, rewrite):
    formula1 = spline(formula)
    formula2 = rewrite_monomials(formula1, rewrite)

    while reverse_spline(formula1) != reverse_spline(formula2):
        formula1 = formula2
        formula2 = rewrite_monomials(formula2, rewrite)

    return reverse_spline(formula2)
