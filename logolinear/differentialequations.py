import numpy
from numpy import sqrt
from logolinear import LogoLinear, BellC
from logolinear.utils import nloop
from numpy.polynomial import Polynomial as P
from fractions import Fraction as F

class DifferentialEquation(object):

    def invA(self, j, s):
        return numpy.array([
            [3*j*j-j-4., 0, 3.*j, -s*2*sqrt(6)/3],
            [0, 3*j*j-j-4., -s*3*sqrt(6), 3*j-1.],
            [0, 0, 3.*j*j, -s*2*sqrt(6)*j/3.],
            [0, 0, -s*3*sqrt(6)*j, (3*j-1.)*j]
            ])/(j*(3*j*j-j-4.))

    def F(self, y, j):
        N, phi, h, v = y
        indices = v.j()
        return numpy.array([
            P(0),
            P(0), 
            self.V(phi, j)/3. - sum(h[p]*h[q]+v[p]*v[q]/3. for p, q in nloop(indices, 2) if p!=j!=q and p+q==j),
            -self.dVdphi(phi, j) - 3*sum(h[p]*v[q] for p, q in nloop(indices, 2) if p!=j!=q and p+q==j),
            ])

    def set_kinetic_terms(self, y, N_p, phi_p, s):
        N, p, h, v = y
        N[0] = [N_p, 1./3]
        p[0] = [phi_p,s*numpy.sqrt(2./3)]
        h[0] = 1./3
        v[0] = numpy.sqrt(2/3)

    def set_curvature_terms(self, y, b, s):
        N, p, h, v = y
        N['4/3'] = -9./14 * b
        p['4/3'] = s*27*sqrt(6)/56 * b
        h['4/3'] = -6./7 * b
        v['4/3'] = s*9*sqrt(6)/14 * b

    def recurse(self, y, j, s):
        # Load F array
        Fj = self.F(y, j)
        Fjk = numpy.zeros((max(len(f) for f in Fj),4))
        for k, f in enumerate(Fj):
            Fjk[:len(f),k] = f.coef

        # Compute the xjk via recursion
        xjk = numpy.zeros((len(Fjk)+1,4))
        for k in reversed(range(1,len(xjk))):
            xjk[k-1] = self.invA(j, s).dot(Fjk[k-1] - k * xjk[k])

        # Load into the polynomials
        for k, _ in enumerate(y):
            y[k][j] = xjk.T[k]

    def compute(self, N_p=0, phi_p=-23, s=+1, b=0, order=16):
        s = numpy.sign(s)
        y = [LogoLinear() for _ in range(4)]
        indices = [F(2*k,3) for k in range(order)]

        for j in indices:
            N, phi, h, v = y
            if j == 0:
                self.set_kinetic_terms(y, N_p, phi_p, s)
            elif j == F('4/3'):
                self.set_curvature_terms(y, b, s)
            else:
                self.recurse(y, j, s)
        return y


class FreeField(DifferentialEquation):

    def __init__(self, m):
        self.m = m

    def V(self, phi, j):
        indices = phi.j()
        ans = P(0)
        if j>=2:
            ans = ans + 1./2 * self.m**2 * sum( phi[p] * phi[q] for p, q in nloop(indices, 2) if p+q==j-2)
        return ans

    def dVdphi(self, phi, j):
        ans = P(0)
        if j>=2:
            ans = ans + self.m**2 *  phi[j-2]
        return ans

class Polynomial(DifferentialEquation):

    def __init__(self, m, l, L):
        self.m = m
        self.l = l
        self.L = L

    def V(self, phi, j):
        indices = phi.j()
        ans = P(0)
        if j==2:
            ans = ans + self.L
        if j>=2:
            ans = ans + 1./2 * self.m**2 * sum( phi[p] * phi[q] for p, q in nloop(indices, 2) if p+q==j-2)
            ans = ans + 1./24 * self.l * sum( phi[p] * phi[q] * phi[r] * phi[s] for p, q, r, s in nloop(indices, 4) if p+q==j-2) 
        return ans

    def dVdphi(self, phi, j):
        ans = P(0)
        if j>=2:
            ans = ans + self.m**2 *  phi[j-2]
            ans = ans + self.l/6. * sum( phi[p] * phi[q] * phi[r]  for p, q, r in nloop(indices, 3) if p+q==j-2)
        return ans

class Starobinsky(DifferentialEquation):

    def __init__(self, L):
        self.L = L

    def V(self, phi, j):
        ans = P(0)
        Phi_p = numpy.exp(-phi[0].coef[0]*sqrt(2./3)) 
        s = F(2*int(numpy.sign(phi[0].coef[1])),3)
        if j==2:
            ans = ans + 1
        if j>=2-s:
            ans = ans -2*Phi_p    * BellC(-sqrt(2./3)*phi,j-(2-s))
        if j>=2-2*s:
            ans = ans +  Phi_p**2 * BellC(-2*sqrt(2./3)*phi,j-(2-2*s))
        return ans * self.L**4

    def dVdphi(self, phi, j):
        ans = P(0)
        Phi_p = numpy.exp(-phi[0].coef[0]*sqrt(2./3)) 
        s = F(2*int(numpy.sign(phi[0].coef[1])),3)
        if j>=2-s:
            ans = ans +2*sqrt(2./3)*Phi_p    * BellC(-sqrt(2./3)*phi,j-(2-s))
        if j>=2-2*s:
            ans = ans -2*sqrt(2./3)*Phi_p**2 * BellC(-2*sqrt(2./3)*phi,j-(2-2*s))
        return ans * self.L**4

    def set_curvature_terms(self, y, b, s):
        N, p, h, v = y
        if s >0:
            Phi_p = numpy.exp(-p[0].coef[0]*sqrt(2./3)) 
            N['4/3'] = -9.*b/14. - 9 * self.L**4*Phi_p/4 - 27*self.L**8*Phi_p**4/50
            p['4/3'] = 27*sqrt(6)/56 * ( b + 2*self.L**4*Phi_p/9 + 49*self.L**8*Phi_p**4/300)
            h['4/3'] = -6.*b/7. - 6 * self.L**4*Phi_p/7 - 18*self.L**8*Phi_p**4/25 
            v['4/3'] = 9*sqrt(6)/14 * ( b + 2*self.L**4*Phi_p/9 + 49*self.L**8*Phi_p**4/300)
        else:
            N['4/3'] = -9./14 * b
            p['4/3'] = s*27*sqrt(6)/56 * b
            h['4/3'] = -6./7 * b 
            v['4/3'] = s*9*sqrt(6)/14 * b
