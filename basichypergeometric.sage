class qPochhammerSymbol():

    def __init__(self, x, q, n):
        self.x = x if type(x) is list else [x]
        self.q = q
        self.n = n

    def _latex_(self):
        return "(%s;%s)_{%s}" % (str(self.x), str(self.q), str(self.n))

    def evaluate(self):
        result = 1 

        for elm in self.x:
            result *= prod([ 1 - elm*self.q**i for i in xrange(self.n) ])

        return result    

class BasicHypergeometricSeries():

    def __init__(self, nominator, denominator, q, z):
        self.nominator = nominator
        self.denominator = denominator
        self.q = q
        self.z = z

    def _latex_(self):
        return "{}_{%i}\\phi_{%i}" % (len(self.nominator), \
                                      len(self.denominator))

    def evaluate(self):
        result, n = 0, 0
        nom, denom, q, z = self.nominator, self.denominator, self.q, self.z

        while True:
            coefficient = qPochhammerSymbol(nom, q, n).evaluate() \
                / qPochhammerSymbol(denom + [q], q, n).evaluate()
            
            if coefficient == 0:
                return result
            
            result += coefficient*z**n
            n += 1
