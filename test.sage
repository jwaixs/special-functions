from sage.symbolic.expression import Expression

class new_exp(Expression):
    def __init__(self, SR, x=0, y=0):
        self.x, self.y = x, y
        Expression.__init__(self, SR, x*y)

    def test(self):
        return 'bla'


class qPochhammerSymbol(Expression):
    def __init__(self, SR, list_a, q, n):
        self.list_a = list_a if type(list_a) == list else [list_a]
        self.q = q
        self.n = n
        Expression.__init__(self, SR, self.evaluate())

    def _latex_(self):
        return "(%s;%s)_{%s}" % (','.join(map(str, self.list_a)), str(self.q), \
                                 str(self.n))

    def __repr__(self):
        return "(%s;%s)_%s" % (','.join(map(str, self.list_a)), str(self.q), \
                               str(self.n))

    def _eval_single(self, a, q):
        ret = 1

        for i in xrange(self.n):
            ret *= 1 - a*q**i

        return ret

    def evaluate(self):
        if self.n >= 0:
            return prod([ self._eval_single(elm, self.q) \
                            for elm in self.list_a ])
        
        ## TODO:: Implement fully symbolic qPochhammer symbol
        else:
            return 0


class BasicHypergeometricSeries(Expression):
    def __init__(self, SR, list_a, list_b, q, z):
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

        Expression.__init__(self, SR, self.evaluate())

    def _latex_(self):
        return "{}_{%i}\phi_{%i}\\biggl[\genfrac..{0pt}{}{%s}{%s};%s,%s\\biggr]" \
            % (len(self.list_a), len(self.list_b), \
            ','.join(map(str, self.list_a)), \
            ','.join(map(str, self.list_b)), \
            str(self.q), str(self.z))

    def __repr__(self):
        return "(%i)_phi_(%i)(%s;%s;%s,%s)" % (len(self.list_a), \
            len(self.list_b), \
            ','.join(map(str, self.list_a)), \
            ','.join(map(str, self.list_b)), \
            str(self.q), str(self.z))

    def __getitem__(self, key):
        if key >= 0:
            j, k = len(self.list_a), len(self.list_b)

            nominator = qPochhammerSymbol(SR, self.list_a, self.q, key)
            denominator = qPochhammerSymbol(SR, self.list_b, self.q, key) \
                * qPochhammerSymbol(SR, self.q, self.q, key)

            return nominator / denominator \
                * ((-1)**key * self.q**(binomial(key,2)))**(1 + k - j) \
                * self.z**key
        else:
            return 0

    def __getslice__(self, i, j):
        ret = []
        k = i
        
        while k < j:
            if self[k] == 0:
                return ret
            ret.append(self[k])
            k += 1

        return ret
            

    def __iter__(self):
        i = 0
        coef = self[i]

        while coef != 0:
            yield coef # Yield! \q.0/
            i += 1
            coef = self[i]

    def evaluate(self):
        rsum, index = 0, 0

        while self[index] != 0:
            rsum += self[index]
            index += 1
        
        return rsum

class Askey_Wilson(Expression):
    def __init__(self, SR, n, z, a, b, c, d, q):
        self.n = n
        self.z = z
        self.q = q
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.param = [a, b, c, d]

        Expression.__init__(self, SR, self.evaluate())

    def __repr__(self):
        return 'p_%i(%s;%s,%s,%s,%s|%s)' % (
            self.n, self.z, self.a, self.b, self.c, self.d, self.q
        )

    def evaluate(self):
        n, q, z, a, b, c, d = [self.n, self.q, self.z] + self.param

        lc = qPochhammerSymbol(SR, [a*b, a*c, a*d], q, n) / a**n
        poly = BasicHypergeometricSeries(SR,
            [q**(-n), a*b*c*d*q**(n-1), a*z, a*z**(-1)],
            [a*b, a*c, a*d], q, q) 
        return lc*poly

    
