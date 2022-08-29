#Functions used in multiple scalar coil (SC) scripts collected together for uniform implementation.
#
# created: 2020/04/21 by M. McCrea
from scipy.optimize import curve_fit

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages #for creating multipage PDF

from datetime import datetime

from matplotlib.pyplot import cm

#user functions and libraries
from Pislib import *

from sympy import integrate
from sympy import pi
import sympy as sym
from sympy import lambdify

def CreateArbGlm(l, m, Glm):
  '''
  function for creating symbolic combination of arbitrary harmonic terms.
  l, m and Glm are equal length ordered lists of values giving the expansion coefficient Glm for the matching l and m arrays.
  '''
  Bx = 0
  By = 0
  Bz = 0
  for (il,im,iGlm) in zip(np.atleast_1d(l),np.atleast_1d(m),np.atleast_1d(Glm)):
    gen = scalarpotential(il,im)
    Bx += iGlm*gen.Pix
    By += iGlm*gen.Piy
    Bz += iGlm*gen.Piz
  
  return Bx,By,Bz

def CreateGlmToOrder(l,Glm):
  '''
  function for creating symbolic combination of arbitrary harmonic terms up to order l for all m.
  l and Glm are equal length ordered lists of values giving the expansion coefficient Glm for the matching l.
  '''
  Bx = 0
  By = 0
  Bz = 0
  
  l_list = []
  m_list = []
  
  for il in range(l+1):
    for im in np.mgrid[-il-1:il+2]:
      #print ("il = " , il , " im = " , im)
      l_list.append(il)
      m_list.append(im)
  
  for (il,im,iGlm) in zip(np.atleast_1d(l_list),np.atleast_1d(m_list),np.atleast_1d(Glm)):
      gen = scalarpotential(il,im)
      Bx += iGlm*gen.Pix
      By += iGlm*gen.Piy
      Bz += iGlm*gen.Piz
  return Bx,By,Bz


def CreateBListToOrder(l):
  '''
  function for creating symbolic combination of arbitrary harmonic terms up to order l for all m.
  l and Glm are equal length ordered lists of values giving the expansion coefficient Glm for the matching l.
  '''
  Bx = []
  By = []
  Bz = []
  
  l_list = []
  m_list = []
  
  for il in range(l+1):
    for im in np.mgrid[-il-1:il+2]:
      #print ("il = " , il , " im = " , im)
      l_list.append(il)
      m_list.append(im)
      
  x=sym.Symbol('x')
  y=sym.Symbol('y')
  z=sym.Symbol('z')


  for (il,im) in zip(np.atleast_1d(l_list),np.atleast_1d(m_list)):
      gen = scalarpotential(il,im)
      
      tBx = lambdify([x,y,z],gen.Pix,"numpy")
      tBy = lambdify([x,y,z],gen.Piy,"numpy")
      tBz = lambdify([x,y,z],gen.Piz,"numpy")
      
      Bx.append(tBx)
      By.append(tBy)
      Bz.append(tBz)
      
  return np.array(Bx),np.array(By),np.array(Bz)

def CreateArbGlmLambdafied(l, m, Glm):
  '''
  identical harmonic functions as CreateArbGlm(), but lambdified for faster execution
  '''
  Bx,By,Bz = CreateArbGlm(l,m,Glm)
  
  x=sym.Symbol('x')
  y=sym.Symbol('y')
  z=sym.Symbol('z')
  
  Bx = lambdify([x,y,z],Bx,"numpy")
  By = lambdify([x,y,z],By,"numpy")
  Bz = lambdify([x,y,z],Bz,"numpy")
  
  return Bx,By,Bz


def CellIntegral(f,cell_rad=0.18, cell_btm=0.05, cell_top=0.2005):
  '''
  symbolically integrates sympy functions over upper cell volume
  distances are meteres
  '''
  x=sym.Symbol('x')
  y=sym.Symbol('y')
  z=sym.Symbol('z')
  
  intz = integrate(f,(z,cell_btm,cell_top))
  intx = integrate(intz,(x, -(cell_rad**2-y**2)**0.5, (cell_rad**2-y**2)**0.5))
  inty = integrate(intx,(y, -cell_rad, cell_rad))
  # print("intx = " , intx)
  # print("intz = " , intz)
  # print("inty = " , inty)
  # print("inty = " , inty.evalf())
  return intx,inty,intz


def PlotArbGlm(l,m,Glm,x_min = -0.2, x_max = 0.2, x_steps = 51j, columns = 3, **plot_kwargs):
  rows = int(len(Glm) / columns) + 1
  print('rows=',rows)
  
  #overlayed plot
  figOver,axOver =plt.subplots(nrows = 1, ncols = 1, figsize=(9,9))
  title = 'Magnetic Field Fit Components on Z-Axis\n Created:%s'%(datetime.today().strftime('%Y-%m-%d %H:%M:%S'))
  print("SC000 plot title = " , title)
  figOver.suptitle(title)
  axOver.set_xlabel("Z-xis Position")
  axOver.set_ylabel("B_z (T)")
  
  #combined plot
  figComb,axComb =plt.subplots(nrows = 1, ncols = 1, figsize=(9,9))
  title = 'Magnetic Components on Z-Axis\n Created:%s'%(datetime.today().strftime('%Y-%m-%d %H:%M:%S'))
  print("SC000 plot title = " , title)
  figComb.suptitle(title)
  
  #initializing PDF for output:
  pdfout = PdfPages('foo.pdf')
  
  #field calculation positions
  points1d=np.mgrid[x_min:x_max:x_steps]
  zeros = np.zeros_like(points1d)
  
  #combined harmonics plot:
  Bx, By, Bz = CreateArbGlmLambdafied(l, m, Glm)
  
  y = Bx(zeros,zeros,points1d)
  if np.isscalar(y):#deal with case of constant function returning scalar
    y = np.full_like(zeros, y)
  axComb.plot(points1d,y,label="Bx(0,0,z)",color="black")
  y = By(zeros,zeros,points1d)
  '''
  if np.isscalar(y):#deal with case of constant function returning scalar
    y = np.full_like(zeros, y)
  axComb.plot(points1d,y,label="By(0,0,z)",color="green")
  y = Bz(zeros,zeros,points1d)
  if np.isscalar(y):#deal with case of constant function returning scalar
    y = np.full_like(zeros, y)
  axComb.plot(points1d,y,label="Bz(0,0,z)",color="blue")
  '''
  axComb.legend()
  axComb.set_xlabel("Axis Position")
  axComb.set_ylabel("B Field Value (T)")
  
  pdfout.savefig(figComb)

  color=iter(cm.rainbow(np.linspace(0,1,len(Glm))))#list of non-repeating colors of equal to number of lines to plot
  
  #going over all l, m, and Glm and creating a scalarpotential function to plot along the z-axis individually.
  for (il, im, iGlm) in zip(l, m, Glm):
    func = scalarpotential(il,im)#note: adding Glm to PisLib.py
    iBz = func.fPiz
    sBz = iGlm*func.Piz
    title = "Bz=%s"%str(sBz) #str() get printable string describing the function
    
    figt,axt =plt.subplots(nrows = 1, ncols = 1, figsize=(9,9))
    axt.set_title(title)
    axt.set_xlabel("Z-xis Position")
    axt.set_ylabel("B_z (T)")
    c=next(color) #change color for each line
    y = iBz(zeros,zeros,points1d)
    #deal with case of constant function returning scalar
    if np.isscalar(y):
      y = np.full_like(zeros, y)
    axOver.plot(points1d,y,label=title ,color=c) #plot overlayed
    axt.plot(points1d,y,label=title ,color=c) #plot individual
    pdfout.savefig(figt)
  
  
  pdfout.savefig(figOver)
  # plt.show()
  pdfout.close()

  return 1