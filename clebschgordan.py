from basichypergeometric import qPochhammerSymbol, BasicHypergeometricSeries

def coef(l1, l2, l, p, n):
    q = var('q')

    if p >= 0:
        result1 = q**(2*n*(l1+l2) + n)

        qpoch1 = qPochhammerSymbol([q**(-4*l2), q**(-4*l1 + 2*p)], q**2, n)
        qpoch2 = qPochhammerSymbol([q**2, q**(2*p + 2)], q**2, n)

        result2 = sqrt(qpoch1.evaluate(), qpohc2.evaluate())
        
        list1 = [q**(-2*n), q**(-2*l1 - 2*l2 + 2*l), q**(-2*l1 - 2*l2 - 2*l - 2)]
        list2 = [q**(-4*l1 + 2*p), q**(-4*l2)]

        result3 = BasicHypergeometricSeries(list1, list2, q**2, q**2).evaluate()

    else:
        result1 = q**(2*n*(l1 + l2) + n)
        
        qpoch1 = qPochhammerSymbol([q**(-4*l1), q**(-4*l2 - 2*p)], q**2, n)
        qpoch2 = qPochhammerSymbol([q**2, q**(-2*p + 2)], q**2, n)

        result2 = sqrt(qpoch1.evaluate() / qpoch2.evaluate())

        list1 = [q**(-2*n), q**(-2*l1 - 2*l2 + 2*l), \
            q**(-2*l1 - 2*l2 - 2*l - 2)]
        list2 = [q**(-4*l2 - 2*p), q**(-4*l1)]

        result3 = BasicHypergeometricSeries( \
            list1, list2, q**2, q**2).evaluate()


    return mathematica(result1*result2*result3).FullSimplify().sage() 
