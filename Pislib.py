#!/usr/bin/python

# calculates Sigma and Pi from

# http://inspirehep.net/record/1703751

# using symbolic algebra

# Wed Feb 12 15:14:12 CST 2020 Jeff and Mark

# scalarpotential class supplies a variety of useful (self) functions.

from sympy import assoc_legendre
from sympy import cos,sin,Abs,factorial
from sympy import sqrt
from sympy import trigsimp
from sympy import Derivative
from sympy import lambdify

import sympy as sym

class scalarpotential:
    def __init__(self,ell=2,m=0):

        self.ell=ell
        self.m=m
        
        x=sym.Symbol('x')
        y=sym.Symbol('y')
        z=sym.Symbol('z')

        r=sym.Symbol('r',real=True,positive=True)
        theta=sym.Symbol('theta',real=True,positive=True)
        phi=sym.Symbol('phi',real=True,positive=True)

        def c(ell,m):
            if(m>=0):
                return (factorial(ell-1)*(-2)**Abs(m))/factorial(ell+Abs(m))*cos(m*phi)
            else:
                return (factorial(ell-1)*(-2)**Abs(m))/factorial(ell+Abs(m))*sin(Abs(m)*phi)

        legendre_ell=ell+1
        
        Sigma=c(legendre_ell,m)*r**legendre_ell*assoc_legendre(legendre_ell,Abs(m),cos(theta))
        Sigma=sym.simplify(Sigma)
        Sigma=Sigma.subs({Abs(sin(theta)):sin(theta)})
        Sigma=Sigma.expand(trig=True)
        self.Sigma_spherical=Sigma
        Sigma=Sigma.subs({r*cos(theta):z,
                          r*sin(theta)*sin(phi):y,
                          r*sin(theta)*cos(phi):x})
        Sigma=Sigma.subs({r**2*sin(theta)**2:x**2+y**2})
        Sigma=Sigma.subs({r**2:x**2+y**2+z**2})
        Sigma=sym.simplify(Sigma.expand())
        self.Sigma=Sigma

        Pix=Derivative(Sigma,x)
        Piy=Derivative(Sigma,y)
        Piz=Derivative(Sigma,z)

        Pix=Pix.doit()
        Piy=Piy.doit()
        Piz=Piz.doit()

        Pix=sym.simplify(Pix.expand())
        Piy=sym.simplify(Piy.expand())
        Piz=sym.simplify(Piz.expand())

        self.Pix=Pix
        self.Piy=Piy
        self.Piz=Piz
        
        self.fPix=lambdify([x,y,z],Pix)
        self.fPiy=lambdify([x,y,z],Piy)
        self.fPiz=lambdify([x,y,z],Piz)

