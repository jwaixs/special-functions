from sage.symbolic.expression import Expression

load('qhyper.sage')

class AskeyWilson(Expression):
    def __init__(self, n, z, a, b, c, d, q):
        self.n = n
        self.z = z
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.q = q
        
        if self.n in ZZ:
            Expression.__init__(self, SR, self.evaluate())
        else:
            Expression.__init__(self, SR, function(
                "p_(%s)(%s;%s;%s;%s;%s|%s)" % \
                (self.n, self.z, self.a, self.b, self.c, self.d, self.q))

    def _latex_(self):
        return "p_{%s}(%s;%s,%s,%s,%s|%s)" % \
            (self.n, self.z, self.a, self.b, self.c, self.d, self.q)

    def __repr__(self):
        return "p_(%s)(%s;%s,%s,%s,%s|%s)" % \
            (self.n, self.z, self.a, self.b, self.c, self.d, self.q)

    def evaluate(self):
        n, q, z, a, b, c, d = self.n, self.q, self.z, self.a, self.b, self.c, \
            self.d

        lc = qPoch([a*b, a*c, a*d], q, n) / a**n
        poly = BHS([q**(-n), a*b*c*d*q**(n-1), a*z, a*z**(-1)],
            [a*b, a*c, a*d], q, q) 
        return lc*poly

    def weight(self):
        q, z, a, b, c, d = self.n, self.q, self.z, self.a, self.b, self.c, \
            self.d
        
        h = lambda alpha : qPoch([alpha*z, alpha*z**(-1)], q, infty)
        
        return h(1) * h(-1) * h(q**(1/2)) * h(-q**(1/2)) \
            / (h(a) * h(b) * h(c) * h(d))

class qRacah(Expression):
    def __init__(self, n, mux, a, b, c, d, q):
        self.n = n
        self.mux = mux
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.q = q
        
        if self.n in ZZ:
            Expression.__init__(self, SR, self.evaluate())
        else:
            Expression.__init__(self, SR, function(
                "p_(%s)(mu(%s);%s;%s;%s;%s|%s)" % \
                (self.n, self.mux, self.a, self.b, self.c, self.d, self.q))

    def _latex_(self):
        return "p_{%s}(\\mu(%s);%s,%s,%s,%s|%s)" % \
            (self.n, self.mux, self.a, self.b, self.c, self.d, self.q)

    def __repr__(self):
        return "p_(%s)(%s;%s,%s,%s,%s|%s)" % \
            (self.n, self.mux, self.a, self.b, self.c, self.d, self.q)

    def evaluate(self):
        n, q, mux, a, b, c, d = self.n, self.q, self.mux, self.a, self.b, \
            self.c, self.d

        poly = BHS([q**(-n), a*b*q**(n+1), q**(-mux), c*d*q**(mux+1)],
            [a*q, b*d*q, c*q], q, q) 
        return poly
