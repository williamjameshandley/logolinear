#!/usr/bin/env python3

import numpy
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import logolinear
from logolinear.differentialequations import FreeField, Starobinsky
from utils import fig_format


# Derivative for odeint
def f(y, t):
    N, phi, H, v, eta = y
    return [H, v, -H**2 - (v**2 - V_(phi))/3, -3*v*H - dVdphi_(phi), numpy.exp(-N)]

# Starobinsky potential
L = numpy.sqrt(1e-5)
def V_(phi):
    return L**4*(numpy.exp(-numpy.sqrt(2/3)*phi)-1)**2

def dVdphi_(phi):
    return 2*L**4*(numpy.exp(-numpy.sqrt(2/3)*phi)-1)*-numpy.sqrt(2/3)*numpy.exp(-numpy.sqrt(2/3)*phi)

# Starobinsky potential plot
fig, ax = plt.subplots()
x = numpy.linspace(-4.2,14.9,1000)
ax.plot(x,V_(x)/L**4,'k-',linewidth=1)
ax.set_xlim(x[0],x[-1])
ax.set_ylim(0,2)
ax.set_xticks([-5,0,5,10])
ax.set_xticklabels([r'$-5m_\mathrm{p}$',r'$0$',r'$5m_\mathrm{p}$',r'$10m_\mathrm{p}$'])
ax.set_yticks([0,1,2])
ax.set_yticklabels([r'$0$',r'$\Lambda^4$',r'$2\Lambda^4$'])
ax.set_xlabel(r'$\phi$')
ax.set_ylabel(r'$V$')
ax.text(0.35, 0.8, r'$V=\Lambda^4(e^{-\sqrt{2/3}\:\phi/m_\mathrm{p}}-1)^2$', fontsize='10', transform=ax.transAxes)

fig.tight_layout()
fig.savefig('star.pdf')


# Starobinsky solutions plot
t = numpy.logspace(0,8,1000)
de = Starobinsky(L)

N_pos, phi_pos, _, _ = de.compute(phi_p=-4.2,order=15,s=+1) 
H_pos, v_pos = logolinear.deriv(N_pos), logolinear.deriv(phi_pos)
y0 = [N_pos(t[0]), phi_pos(t[0]), H_pos(t[0]), v_pos(t[0]), t[0]]
N_pos_, phi_pos_, H_pos_, v_pos_, eta_pos_ = odeint(f,y0,t).T

N_neg, phi_neg, _, _ = de.compute(phi_p=14.9,order=15,s=-1)
H_neg, v_neg = logolinear.deriv(N_neg), logolinear.deriv(phi_neg)
y0 = [N_neg(t[0]), phi_neg(t[0]), H_neg(t[0]), v_neg(t[0]), t[0]]
N_neg_, phi_neg_, H_neg_, v_neg_, eta_neg_ = odeint(f,y0,t).T

fig, ax = plt.subplots()
ax.set_xscale('log')
ax.plot(t, phi_pos_,'k:',linewidth=1)
ax.plot(t, phi_neg_,'k-',linewidth=1)
ax.set_ylim(-4.2,14.9)
ax.set_xlim(t[0],t[-1])
ax.set_xlabel('$t$')
ax.set_ylabel(r'$\phi/m_\mathrm{p}$')

ax.legend(handles=ax.lines,labels=[r'$\phi^+$',r'$\phi^-$'])

fig.tight_layout()
fig.savefig('phi_star.pdf')


# Polynomial potential
m = 1e-5
def V_(phi):
    return 1/2 * m**2 * phi**2

def dVdphi_(phi):
    return  m**2 * phi

# Error plot
fig, ax = plt.subplots()
t = numpy.linspace(10000,22000,1000)
ax.set_ylim(1e-5,1)
ax.set_xlim(t[0],t[-1])
ax.set_yscale('log')

de = FreeField(m)
order = 50
N, phi, _, _ = de.compute(order=order*3//2)
H = logolinear.deriv(N)
v = logolinear.deriv(phi)

t = numpy.linspace(t[0],t[-1],1000)
y0 = [N(t[0]), phi(t[0]), H(t[0]), v(t[0]), t[0]]
N_, phi_, H_, v_, eta_ = odeint(f,y0,t).T

orders = numpy.arange(0,order,2)
colors = numpy.linspace(0.0,0.7,len(orders))
for order, color in zip(orders, colors):
    ax.plot(t,abs((phi(t,order)-phi_)/phi_), color=str(color), linewidth=0.5)

handles = ax.lines[0::5]
labels = [r'$\sim\mathcal{O}(t^{%i})$'% o for o in orders[0::5]]
ax.legend(handles=handles,labels=labels)
ax.set_xlabel(r'$t$')
ax.set_ylabel('Absolute error in $\phi$')

fig.tight_layout()
fig.savefig('phi_error.pdf')


# N plot
fig, ax = plt.subplots()

de = FreeField(m)
N, phi, _, _ = de.compute(b=2e-6, phi_p = -23)
H = logolinear.deriv(N)
v = logolinear.deriv(phi)

t = numpy.logspace(0,8,1000)
order = 0
y0 = [N(t[0],order), phi(t[0],order), H(t[0],order), v(t[0],order), t[0]]
N_, phi_, H_, v_, eta_ = odeint(f,y0,t).T
ax.plot(t, N_, 'k:',linewidth=1)

t = numpy.logspace(0,8,1000)
order = 2
y0 = [N(t[0],order), phi(t[0],order), H(t[0],order), v(t[0],order), t[0]]
N_, phi_, H_, v_, eta_ = odeint(f,y0,t).T
ax.plot(t, N_, 'k-',linewidth=1)

ax.set_xscale('log')
ax.set_ylim(0,70)
ax.set_xlim(t[0],t[-1])
ax.set_xlabel('$t$')
ax.set_ylabel('e-folds $N$')
ax.legend(handles=ax.lines,labels=[r'$\sim\mathcal{O}(t^0)$', r'$\sim\mathcal{O}(t^{4/3})$'])

fig.tight_layout()
fig.savefig('N.pdf')
