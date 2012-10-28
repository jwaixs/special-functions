class PochhammerSymbol():

    def __init__(self, x, n):
        self.x = x
        self.n = n

    def _latex_(self):
        return "(%s)_{%s}" % (str(self.x), str(self.n))

    def _repr_(self):
        return "(%s)_{%s}" % (str(self.x), str(self.n))

    def evaluate(self):
        result = 1
        
        if type(self.x) is list:
            for elm in self.x:
                result *= prod([ elm-i for i in xrange(self.n) ])
        else:
            result = prod([ self.x-i for i in xrange(self.n) ])
        
        return result


class HypergeometricSeries():

    def __init__(self, nominator, denominator, z):
        self.nominator = nominator if type(nominator) is list \
                                    else [nominator]
        self.denominator = denominator if type(denominator) is list \
                                    else [denominator]
        self.z = z

    def _latex_(self):
        return "%i_F_%i(%s, %s; %s)" % (len(self.nominator), \
                                        len(self.denominator), \
                                        str(self.nominator), \
                                        str(self.denominator), \
                                        str(self.z))
            
    def evaluate(self):
        result, n = 0, 0

        while True:
            coefficient = PochhammerSymbol(self.nominator, n).evaluate() \
                / PochhammerSymbol(self.denominator + [n], n).evaluate()
            if coefficient == 0:
                break
            result += coefficient*self.z**n
            n += 1

        return result
