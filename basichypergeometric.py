r"""
Basic Hypergeometric Series

A module for computing basic hypergeometric series.

AUTHORS:
    - Noud Aldenhoven - initial file (2012)

REFERENCES:
    - Basic Hypergeometric Series - Gasper, Rahman
"""

from sage.misc.misc import prod
from sage.functions.other import factorial
from sage.rings.arith import binomial

class qPochhammerSymbol():
    r"""
    qPochhammer symbol

    INPUT::
        - ``list_a``    - a list of ring elements.
        - ``q``         - the q-parameter.
        - ``n``         - zero or a positive integer.

    OUTPUT::
        - a qPochhammerSymbol object.

    EXAMPLES::
        
        Let us start with defining qPochhammer symbols::

            sage: qpoch1 = qPochhammerSymbol(5, 1/2, 3)
            sage: qpoch1 
            (5;1/2)_3
            
        with a unknown q-parameter and one unknown variable::

            sage: q, a = var('q a')
            sage: qpoch2 = qPochhammerSymbol(a, q, 3)
            sage: qpoch2
            (a;q)_3

        with an unknown variable label::

            sage: n = var('n')
            sage: qpoch3 = qPochhammerSymbol(a, q, n)
            sage: qpoch3
            (a;q)_n

        with more unknown variables::

            sage: b = var('b')
            sage: qpoch4 = qPochhammerSymbol([a, b], q, 3) 
            sage: qpoch4
            (a,b;q)_3

        We evaluate the qPochhammer symbols as follow::

            sage: qpoch1.evaluate()
            -3/2
            sage: qpoch2.evaluate()
            -(a - 1)*(a*q - 1)*(a*q^2 - 1)
            sage: qpoch3.evaluate()
            0
            sage: qpoch4.evaluate()
            (b - 1)*(a - 1)*(b*q - 1)*(a*q - 1)*(b*q^2 - 1)*(a*q^2 - 1)

    """
    def __init__(self, list_a, q, n):
        """
        See the docstring form :meth:`qPochhammerSymbol`.

        EXAMPLES::

            sage: q, a, b, n = var('q a b n')
            sage: qPochhammerSymbol(5, 1/2, 3)
            (5;1/2)_3
            sage: qPochhammerSymbol(a, q, 3)
            (a;q)_3
            sage: qPochhammerSymbol(a, q, n)
            (a;q)_n
            sage: qPochhammerSymbol([a, b], q, n)
            (a,b;q)_n

        """
        self.list_a = list_a if type(list_a) == list else [list_a]
        self.q = q
        self.n = n

    def _latex_(self):
        """
        EXAMPLES::

            sage: q, a, b, n = var('q a b n')
            sage: latex(qPochhammerSymbol(5, 1/2, 3))
            (5;1/2)_{3}
            sage: latex(qPochhammerSymbol(a, q, 3))
            (a;q)_{3}
            sage: latex(qPochhammerSymbol(a, q, n))
            (a;q)_{n}
            sage: latex(qPochhammerSymbol([a, b], q, n))
            (a,b;q)_{n}
    
        """
        return "(%s;%s)_{%s}" % (','.join(map(str, self.list_a)), str(self.q), \
                                 str(self.n))

    def __repr__(self):
        """
        EXAMPLES::

            sage: q, a, b, n = var('q, a b n')
            sage: qPochhammerSymbol(5, 1/2, 3)
            (5;1/2)_3
            sage: qPochhammerSymbol(a, q, 3)
            (a;q)_3
            sage: qPochhammerSymbol(a, q, n)
            (a;q)_n
            sage: qPochhammerSymbol([a, b], q, n)
            (a,b;q)_n
        
        """
        return "(%s;%s)_%s" % (','.join(map(str, self.list_a)), str(self.q), \
                               str(self.n))

    def _eval_single(self, a, q):
        """
        Private function to compute a single Pochhammer symbol.

        EXAMPLES::

            sage: q, a, n = var('q a n')
            sage: qPochhammerSymbol(5, 1/2, 3).evaluate()
            -3/2
            sage: qPochhammerSymbol(a, 3).evaluate()
            -(a - 1)*(a*q - 1)*(a*q^2 - 1)
            sage: qPochhammerSymbol(a, n).evaluate()
            0
        
        """
        ret = 1

        for i in xrange(self.n):
            ret *= 1 - a*q**i

        return ret

    def evaluate(self):
        """
        Evaluate qPochhammer symbols.

        EXAMPLES::

            sage: q, a, b, n = var('1 a b n')
            sage: qPochhammerSymbol(5, 1/2, 3).evaluate()
            -3/2
            sage: qPochhammerSymbol(a, q, 3).evaluate()
            -(a - 1)*(a*q - 1)*(a*q^2 - 1)
            sage: qPochhammerSymbol(a, q, n).evaluate() # Not implemented yet
            0
            sage: qPochhammerSymbol([a, b], q, 3).evaluate()
            (b - 1)*(a - 1)*(b*q - 1)*(a*q - 1)*(b*q^2 - 1)*(a*q^2 - 1)
        
        """
        if self.n >= 0:
            return prod([ self._eval_single(elm, self.q) \
                            for elm in self.list_a ])
        
        ## TODO:: Implement fully symbolic qPochhammer symbol
        else:
            return 0

class BasicHypergeometricSeries():
    """
    An implementation of the Hypergeometric series.

    INPUT::
        - `list_a`  - list of nominators.
        - `list_b`  - list of denominators.
        - `z`       - the power of the Hypergeometric series.

    OUTPUT::
        - An object HypergeometriSeries which can be used the compute or
          approximate the Hypergeometric series.

    EXAMPLES::

        This example shows the basic properties of a summable Hypergeometric 
        series::

            sage: a, b, z = var('a b z')
            sage: hs = HypergeometricSeries([a, -3], [b], z)
            sage: hs
            (2)_F_(1)(a,-3;b;z)
            sage: latex(hs)
            {}_{2}F_{1}\biggl[\genfrac..{0pt}{}{a,-3}{b};z\biggr]
            sage: hs.evaluate()
            -(a + 1)*(a + 2)*a*z^3/((b + 1)*(b + 2)*b) + 3*(a + 1)*a*z^2/((b + 1)*b) - 3*a*z/b + 1           
 
            sage: hs[3]
            -(a + 1)*(a + 2)*a*z^3/((b + 1)*(b + 2)*b)
            sage: hs[4]
            0
            sage: hs[1:3]
            [-3*a*z/b, 3*(a + 1)*a*z^2/((b + 1)*b)]
            sage: hs[0:]
            [1, -3*a*z/b, 3*(a + 1)*a*z^2/((b + 1)*b), -(a + 1)*(a + 2)*a*z^3/((b + 1)*(b + 2)*b)]            
 
            sage: for coef in hs: print coef
            1
            -3*a*z/b
            3*(a + 1)*a*z^2/((b + 1)*b)
            -(a + 1)*(a + 2)*a*z^3/((b + 1)*(b + 2)*b)

        We can for example approximate the exponential, sine and cosine using 
        Hypergeometric series::

            sage: z = var('z')

            sage: hs_exp = HypergeometricSeries([], [], z)
            sage: sum(hs_exp[:5])
            1/24*z^4 + 1/6*z^3 + 1/2*z^2 + z + 1
            sage: exp(z).taylor(z, 0, 4)
            1/24*z^4 + 1/6*z^3 + 1/2*z^2 + z + 1
            sage: bool(sum(hs_exp[:5]) == exp(z).taylor(z, 0, 4))
            True

            sage: hs_cos = HypergeometricSeries([], [1/2], -1/4*z**2)
            sage: sum(hs_cos[:5])
            1/40320*z^8 - 1/720*z^6 + 1/24*z^4 - 1/2*z^2 + 1
            sage: cos(z).taylor(z, 0, 8)
            1/40320*z^8 - 1/720*z^6 + 1/24*z^4 - 1/2*z^2 + 1
            sage: bool(cos(z).taylor(z, 0, 8) == sum(hs_cos[:5]))
            True
            
            sage: hs_sin = HypergeometricSeries([], [3/2], -1/4*z**2)
            sage: z*sum(hs_sin[:5])
            1/362880*(z^8 - 72*z^6 + 3024*z^4 - 60480*z^2 + 362880)*z
            sage: sin(z).taylor(z, 0, 9)
            1/362880*z^9 - 1/5040*z^7 + 1/120*z^5 - 1/6*z^3 + z
            sage: bool(sin(z).taylor(z, 0, 9) == z*sum(hs_sin[:5]))
            True

    """ 
    def __init__(self, list_a, list_b, q, z):
        """
        See the docstring form :meth:`HypergeometricSeries`.

        EXAMPLES::

            sage: a, b, c, z = var('a b c z')
            sage: HypergeometricSeries([a, b], [c], z)
            (2)_F_(1)(a,b;c;z)

        """
        self.list_a = list_a
        self.list_b = list_b
        self.q = q
        self.z = z

    def _latex_(self):
        """
        EXAMPLES::

            sage: a, b, c, z = var('a b c z')
            sage: latex(HypergeometricSeries([a, b], [c], z))
            {}_{2}F_{1}\biggl[\genfrac..{0pt}{}{a,b}{c};z\biggr]

        """
        return "{}_{%i}\phi_{%i}\\biggl[\genfrac..{0pt}{}{%s}{%s};%s,%s\\biggr]" \
            % (len(self.list_a), len(self.list_b), \
            ','.join(map(str, self.list_a)), \
            ','.join(map(str, self.list_b)), \
            str(self.q), str(self.z))

    def __repr__(self):
        """
        EXAMPLES::

            sage: a, b, c, z = var('a b c z')
            sage: HypergeometricSeries([a, b], [c], z)
            (2)_F_(1)(a,b;c;z)

        """
        return "(%i)_phi_(%i)(%s;%s;%s,%s)" % (len(self.list_a), \
            len(self.list_b), \
            ','.join(map(str, self.list_a)), \
            ','.join(map(str, self.list_b)), \
            str(self.q), str(self.z))

    def __getitem__(self, key):
        """
        EXAMPLES::

            sage: a, b, c, z = var('a b c z')
            sage: hs = HypergeometricSeries([a, b], [c], z)
            sage: for i in range(4): print hs[i]
            1
            a*b*z/c
            1/2*(b + 1)*(a + 1)*a*b*z^2/((c + 1)*c)
            1/6*(b + 1)*(b + 2)*(a + 1)*(a + 2)*a*b*z^3/((c + 1)*(c + 2)*c)

        """
        if key >= 0:
            nominator = qPochhammerSymbol(self.list_a, self.q, key).evaluate()
            denominator = qPochhammerSymbol( \
                self.list_b, \
                self.q, \
                key \
            ).evaluate() * qPochhammerSymbol(self.q, self.q, key).evaluate()
            return nominator / denominator \
                * ((-1)**key * self.q**(binomial(key,2))) \
                * self.z**key
        else:
            return 0

    def __getslice__(self, i, j):
        """
        EXAMPLES::

            sage: a, b, c, z = var('a b c z')
            sage: hs = HypergeometricSeries([a, b], [c], z)
            sage: hs[0:3]
            [1, a*b*z/c, 1/2*(b + 1)*(a + 1)*a*b*z^2/((c + 1)*c)]
            sage: hs[1:3]
            [a*b*z/c, 1/2*(b + 1)*(a + 1)*a*b*z^2/((c + 1)*c)]
            sage: hs[:2]
            [1, a*b*z/c]
            sage: hs[3:2]
            []

        """
        ret = []
        k = i
        
        while k < j:
            if self[k] == 0:
                return ret
            ret.append(self[k])
            k += 1

        return ret
            

    def __iter__(self):
        """
        EXAMPLES::

            sage: a, b, z = var('a b z')
            sage: hs = HypergeometricSeries([a, -3], [b], z)
            sage: for coef in hs: print coef
            1
            -3*a*z/b
            3*(a + 1)*a*z^2/((b + 1)*b)
            -(a + 1)*(a + 2)*a*z^3/((b + 1)*(b + 2)*b)

        """
        i = 0
        coef = self[i]

        while coef != 0:
            yield coef # Yield! \q.0/
            i += 1
            coef = self[i]

    def evaluate(self):
        """
        EXAMPLES::

            sage: a, b, z = var('a b z')
            sage: hs = HypergeometricSeries([a, -3], [b], z)
            sage: hs.evaluate()
            -(a + 1)*(a + 2)*a*z^3/((b + 1)*(b + 2)*b) + 3*(a + 1)*a*z^2/((b + 1)*b) - 3*a*z/b + 1

        """
        return sum(list(self))
