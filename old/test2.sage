from sage.symbolic.expression import Expression
from sage.symbolic.function import BuiltinFunction
from sage.structure.coerce import parent as sage_structure_coerce_parent

class Function_qpoch(BuiltinFunction):
    def __init__(self):
        BuiltinFunction.__init__(self, 'qpoch', nargs=3)

    def __call__(self, *args, **kwds):
        if len(args) == 3:
            return BuiltinFunction.__call__(\
                self, SR._force_pyobject(args[0]), args[1], args[2], **kwds)
        else:
            raise TypeError('qpoch takes three arguments.')

    def _eval_(self, a, q, n):
        if isinstance(n, Integer):
            return self._evalf_(a, q, n)

    def _evalf_(self, a, q, n):
        ret = 1

        if isinstance(a, list):
            for elm in a:
                ret *= self._eval_(elm, q, n)
        else:
            if n >= 0:
                for i in range(n):
                    ret *= (1 - a*q**i)
            else:
                return 1/(self._eval_(a*q**n, q, -n))

        return ret

qpoch = Function_qpoch()


class Function_bhs(BuiltinFunction):
    def __init__(self):
        BuiltinFunction.__init__(self, 'bhs', nargs=4)

    def __call__(self, *args, **kwds):
        if len(args) == 4:
            return BuiltinFunction.__call__(self, \
                SR._force_pyobject(args[0]), SR._force_pyobject(args[1]), \
                args[2], args[3], **kwds)

    def _eval_(self, nom, denom, q, z):
        for elm in nom:
            coefs = elm.coefficients(q)
            if len(coefs) == 1 and coefs[0][0] in ZZ:
                if coefs[0][1] in ZZ and coefs[0][1] < 0:
                    return self._evalf_(nom, denom, q, z)
        return None

    def _evalf_(self, nom, denom, q, z):
        ret = 1
        _nom = 1
        k = 1
        
        while _nom != 0:
            _nom = qpoch(nom, q, k)
            _denom = qpoch(denom, q, k)
            ret += _nom / _denom * z**(k)
            k += 1

        return ret

bhs = Function_bhs()
