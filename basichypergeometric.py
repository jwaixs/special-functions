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
from sage.rings.infinity import PlusInfinity

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
        return "(%s;%s)_%s" % (','.join(map(str, self.list_a)), str(self.q), \
                               str(self.n))

    def _eval_single(self, a, q):
        """
        Private function to compute a single Pochhammer symbol.

        EXAMPLES::

            sage: q, a, n = var('q a n')
            sage: qPochhammerSymbol(5, 1/2, 3).evaluate()
            -3/2
            sage: qPochhammerSymbol(a, q, 3).evaluate()
            -(a - 1)*(a*q - 1)*(a*q^2 - 1)
            sage: qPochhammerSymbol(a, q, n).evaluate()
            0
            sage: qPochhammerSymbol(q**(-3), q, oo).evaluate()
            -(1/q^3 - 1)*(1/q^2 - 1)*(1/q - 1)
        
        """
        ret = 1

        if type(self.n) == PlusInfinity:
            i, elm = 0, 1 - a
       
            while elm != 0:
                ret *= elm
                i += 1
                elm = 1 - a*q**i
        else:
            for i in xrange(self.n):
                ret *= 1 - a*q**i

        return ret

    def evaluate(self):
        """
        Evaluate qPochhammer symbols.

        EXAMPLES::

            sage: q, a, b, n = var('q a b n')
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
    An implementation of the Basic Hypergeometric Series.

    INPUT::
        - `list_a`  - list of nominators.
        - `list_b`  - list of denominators.
        - `q`       - the q-parameter
        - `z`       - the power of the Hypergeometric series.

    OUTPUT::
        - An object BasicHypergeometriSeries which can be used the compute or
          approximate the Hypergeometric series.

    EXAMPLES::

        This example shows the basic properties of a summable Basic 
        Hypergeometric series::

            sage: a, b, q, z = var('a b q z')
            sage: bhs = BasicHypergeometricSeries([a, q**(-3)], [b], q, z)
            sage: bhs
            (2)_phi_(1)(a,q^(-3);b;q,z)
            sage: bhs.evaluate()
            (1/q^3 - 1)*(1/q^2 - 1)*(1/q - 1)*(a - 1)*(a*q - 1)*(a*q^2 - 1)*z^3/((q - 1)*(b - 1)*(q^2 - 1)*(b*q - 1)*(q^3 - 1)*(b*q^2 - 1)) + (1/q^3 - 1)*(1/q^2 - 1)*(a - 1)*(a*q - 1)*z^2/((q - 1)*(b - 1)*(q^2 - 1)*(b*q - 1)) + (1/q^3 - 1)*(a - 1)*z/((q - 1)*(b - 1)) + 1            
            sage: bhs[3]
            (1/q^3 - 1)*(1/q^2 - 1)*(1/q - 1)*(a - 1)*(a*q - 1)*(a*q^2 - 1)*z^3/((q - 1)*(b - 1)*(q^2 - 1)*(b*q - 1)*(q^3 - 1)*(b*q^2 - 1))
            sage: bhs[4]
            0
            sage: bhs[1:3]
            [(1/q^3 - 1)*(a - 1)*z/((q - 1)*(b - 1)), (1/q^3 - 1)*(1/q^2 - 1)*(a - 1)*(a*q - 1)*z^2/((q - 1)*(b - 1)*(q^2 - 1)*(b*q - 1))]            
            sage: for coef in bhs: print coef
            1
            (1/q^3 - 1)*(a - 1)*z/((q - 1)*(b - 1))
            (1/q^3 - 1)*(1/q^2 - 1)*(a - 1)*(a*q - 1)*z^2/((q - 1)*(b - 1)*(q^2 - 1)*(b*q - 1))
            (1/q^3 - 1)*(1/q^2 - 1)*(1/q - 1)*(a - 1)*(a*q - 1)*(a*q^2 - 1)*z^3/((q - 1)*(b - 1)*(q^2 - 1)*(b*q - 1)*(q^3 - 1)*(b*q^2 - 1))

    """ 
    def __init__(self, list_a, list_b, q, z):
        """
        See the docstring form :meth:`HypergeometricSeries`.

        EXAMPLES::

            sage: a, b, c, q, z = var('a b c q z')
            sage: BasicHypergeometricSeries([a, b], [c], q, z)
            (2)_phi_(1)(a,b;c;q,z)

        TESTS::

            Remove duplicates in the nominator and denominator::

                sage: q, z = var('q z')
                sage: BasicHypergeometricSeries([q**(-2), q**(-3), q**(-4), q**(-3)], [q**(-4), q**(-3), q**(-5)], q, z)
                (2)_phi_(1)(q^(-2),q^(-3);q^(-5);q,z)

            Test the q-binomial theorem::

                sage: bhs1 = BasicHypergeometricSeries([q**(-5)], [], q, z)
                sage: bhs2 = BasicHypergeometricSeries([q**(-5)], [], q, q*z)
                sage: bool(bhs1.evaluate() == bhs2.evaluate()*(1 - q**(-5)*z)/(1 - z))
                True

            Test Heine's transformation formulas::

                

        """
        self.list_a = list_a
        self.list_b = list_b

        # remove duplicates in list_a \cap list_b to avoid dividing by zero
        temp_a = self.list_a[:]
        for elm in temp_a:
            if elm in self.list_b:
                self.list_a.remove(elm)
                self.list_b.remove(elm)

        self.q = q
        self.z = z

    def _latex_(self):
        """
        EXAMPLES::

            sage: a, b, c, q, z = var('a b c q z')
            sage: latex(BasicHypergeometricSeries([a, b], [c], q, z))
            {}_{2}\phi_{1}\biggl[\genfrac..{0pt}{}{a,b}{c};q,z\biggr]

        """
        return "{}_{%i}\phi_{%i}\\biggl[\genfrac..{0pt}{}{%s}{%s};%s,%s\\biggr]" \
            % (len(self.list_a), len(self.list_b), \
            ','.join(map(str, self.list_a)), \
            ','.join(map(str, self.list_b)), \
            str(self.q), str(self.z))

    def __repr__(self):
        """
        EXAMPLES::

            sage: a, b, c, q, z = var('a b c q z')
            sage: BasicHypergeometricSeries([a, b], [c], q, z)
            (2)_phi_(1)(a,b;c;q,z)

        """
        return "(%i)_phi_(%i)(%s;%s;%s,%s)" % (len(self.list_a), \
            len(self.list_b), \
            ','.join(map(str, self.list_a)), \
            ','.join(map(str, self.list_b)), \
            str(self.q), str(self.z))

    def __getitem__(self, key):
        """
        EXAMPLES::

            sage: a, b, c, q, z = var('a b c q z')
            sage: bhs = BasicHypergeometricSeries([a, b], [c], q, z)
            sage: for i in range(4): print bhs[i]
            1
            (b - 1)*(a - 1)*z/((q - 1)*(c - 1))
            (b - 1)*(a - 1)*(b*q - 1)*(a*q - 1)*z^2/((q - 1)*(c - 1)*(q^2 - 1)*(c*q - 1))
            (b - 1)*(a - 1)*(b*q - 1)*(a*q - 1)*(b*q^2 - 1)*(a*q^2 - 1)*z^3/((q - 1)*(c - 1)*(q^2 - 1)*(c*q - 1)*(q^3 - 1)*(c*q^2 - 1))
        
        """
        if key >= 0:
            j, k = len(self.list_a), len(self.list_b)
            nominator = qPochhammerSymbol(self.list_a, self.q, key).evaluate()
            denominator = qPochhammerSymbol( \
                self.list_b, \
                self.q, \
                key \
            ).evaluate() * qPochhammerSymbol(self.q, self.q, key).evaluate()
            return nominator / denominator \
                * ((-1)**key * self.q**(binomial(key,2)))**(1 + k - j) \
                * self.z**key
        else:
            return 0

    def __getslice__(self, i, j):
        """
        EXAMPLES::

            sage: a, b, c, q, z = var('a b c q z')
            sage: bhs = BasicHypergeometricSeries([a, b], [c], q, z)
            sage: bhs[0:3]
            [1, (b - 1)*(a - 1)*z/((q - 1)*(c - 1)), (b - 1)*(a - 1)*(b*q - 1)*(a*q - 1)*z^2/((q - 1)*(c - 1)*(q^2 - 1)*(c*q - 1))]
            sage: bhs[1:3]
            [(b - 1)*(a - 1)*z/((q - 1)*(c - 1)), (b - 1)*(a - 1)*(b*q - 1)*(a*q - 1)*z^2/((q - 1)*(c - 1)*(q^2 - 1)*(c*q - 1))]
            sage: bhs[:2]
            [1, (b - 1)*(a - 1)*z/((q - 1)*(c - 1))]
            sage: bhs[3:2]
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

            sage: a, b, q, z = var('a b q z')
            sage: bhs = BasicHypergeometricSeries([a, q**(-3)], [b], q, z)
            sage: for coef in bhs: print coef
            1
            (1/q^3 - 1)*(a - 1)*z/((q - 1)*(b - 1))
            (1/q^3 - 1)*(1/q^2 - 1)*(a - 1)*(a*q - 1)*z^2/((q - 1)*(b - 1)*(q^2 - 1)*(b*q - 1))
            (1/q^3 - 1)*(1/q^2 - 1)*(1/q - 1)*(a - 1)*(a*q - 1)*(a*q^2 - 1)*z^3/((q - 1)*(b - 1)*(q^2 - 1)*(b*q - 1)*(q^3 - 1)*(b*q^2 - 1))
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

            sage: a, b, q, z = var('a b q z')
            sage: bhs = BasicHypergeometricSeries([a, q**(-3)], [b], q, z)
            sage: bhs.evaluate()
            (1/q^3 - 1)*(1/q^2 - 1)*(1/q - 1)*(a - 1)*(a*q - 1)*(a*q^2 - 1)*z^3/((q - 1)*(b - 1)*(q^2 - 1)*(b*q - 1)*(q^3 - 1)*(b*q^2 - 1)) + (1/q^3 - 1)*(1/q^2 - 1)*(a - 1)*(a*q - 1)*z^2/((q - 1)*(b - 1)*(q^2 - 1)*(b*q - 1)) + (1/q^3 - 1)*(a - 1)*z/((q - 1)*(b - 1)) + 1

        """
        return sum(list(self))
