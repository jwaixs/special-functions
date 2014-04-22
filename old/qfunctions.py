from basichypergeometric import qPochhammerSymbol
from sage.symbolic.expression import Expression

def qbinomial(n, k, q):
    nominator = qPochhammerSymbol(q**n, q**(-1), k).evaluate()
    denominator = qPochhammerSymbol(q, q, k).evaluate()
    return nominator/denominator

def qdifference(f, x, q):
    """qdifference - q-Difference Equation

    INPUT::
        - ``f'' - a function f(x).
        - ``x'' - a variable to substitute.
        - ``q'' - a q-parameter.

    OUTPUT::
        - returns the q-difference of the function f(x)

    EXAMPLES::

            sage: q, x = var('q x')
            sage: f = qdifference(x**2 + x + 1, x, q)
            sage: f
            (q^2*x^2 + q*x - x^2 - x)/((q - 1)*x)
            sage: limit(f, q=1)
            2*x + 1

        Using different symbolic functions:

            sage: y = var('y')
            sage: qdifference(cos(x*y**2), y, q)
            (cos(q^2*x*y^2) - cos(x*y^2))/((q - 1)*y)
    
    """

    if type(x) == Expression and x in f.variables():
        return (f - f.subs({x : q*x})) / ((1 - q) * x)
    else:
        return None
