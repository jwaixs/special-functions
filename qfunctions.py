from basichypergeometric import qPochhammerSymbol

def qbinomial(n, k, q):
    nominator = qPochhammerSymbol(q**n, q**(-1), k).evaluate()
    denominator = qPochhammerSymbol(q, q, k).evaluate()
    return nominator/denominator
