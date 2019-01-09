import numpy
import copy
from numpy.polynomial import Polynomial as P
from fractions import Fraction as F
from collections import OrderedDict
from itertools import islice
from logolinear.utils import nloop


def PolyRep(poly,indices=None,var='x'):
    if indices is None:
        indices = list(range(len(poly)))
    items = []
    for x, i in zip(poly,indices):
        if not x:
            continue
        if i==0:
            items.append('%s' % (x))
        elif x==1:
            items.append('%s^%s' % (var, i))
        elif x==-1:
            items.append('-%s^%s' % (var, i))
        else:
            items.append('%s %s^%s' % (x, var, i))
    result = ' + '.join(items)
    result = result.replace('+ -', '- ')
    result = result.replace('^1', '')
    return result

class LogoLinear(object):
    def __init__(self):
        self._arr = OrderedDict()
        self.var = 't'

    def __call__(self, t, order=None):
        if order is None:
            return sum(xp(numpy.log(t))*t**p for p, xp in self._arr.items())
        else:
            return sum(xp(numpy.log(t))*t**p for p, xp in self._arr.items() if p<=order)

    def __getitem__(self, key):
        if not isinstance(key, slice):
            return self._arr[F(key)]
        else:
            x = LogoLinear()
            if key.start is None and key.stop is None:
                for k in self.j():
                    x[k] = self[k]
            elif key.start is None:
                for k in self.j():
                    if k < key.stop:
                        x[k] = self[k]
            elif key.stop is None:
                for k in self.j():
                    if key.start <= k:
                        x[k] = self[k]
            else:
                for k in self.j():
                    if key.start <= k < key.stop:
                        x[k] = self[k]
            return x

    def __setitem__(self, key, value):
        if not isinstance(value, P):
            value = P(value)
        self._arr[F(key)] = value

    def j(self):
        return self._arr.keys()

    def array(self):
        return numpy.array(self._arr.values())

    def __repr__(self):
        poly = ['(%s)' % PolyRep(xp.coef, var='log%s' % self.var) for xp in self._arr.values()]
        poly = ['' if i=='()' else i for i in poly]
        return PolyRep(poly, self._arr.keys(), var=self.var)

    def __mul__(self, other):
        x = LogoLinear()
        if isinstance(other, LogoLinear):
            for p in self.j():
                for q in other.j():
                    try:
                        x[p+q] += self[p] * other[q]
                    except KeyError:
                        x[p+q] = self[p] * other[q]
        else:
            for key in self._arr.keys():
                x[key] =  self[key] * other
        return x

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return self * -1

    @classmethod
    def t_to_the(cls, p):
        x = cls()
        x[p] = 1
        return x




def deriv(logol):
    ll = LogoLinear()
    for p in logol.j():
        ll[p-1] = p*logol[p] + logol[p].deriv()
    return ll

def inttj(x, j):
    if j == -1:
        return x.integ()
    if len(x) == 1:
        return x/(j+1)
    else:
        return (x - inttj(x.deriv(),j))/(j+1)

def integ(logol):
    ll = LogoLinear()
    ll[0] = 0
    for p in logol.j():
        ll[p+1] = inttj(logol[p],p)
    return ll


def BellC(x, j):
    """ Defines the complete ordinary Bell polynomial recursively. """
    indices = [i for i in x.j() if F(i) <= j]
    indices = [(p,q) for p, q in nloop(indices, 2) if p!=0 and p+q==j]
    if indices:
        s = sum(p/j * BellC(x, q) * x[p] for p, q in indices if x[p]!=P(0))
        if s==0:
            return P(0.)
        else:
            return s
    else:
        return P(1)

def exp(logol):
    ll = LogoLinear()
    for p in logol.j():
        ll[p] = BellC(logol,p)
    return ll
