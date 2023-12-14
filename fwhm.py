import numpy as np
import math
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

def main(x,y,offset=0):
  n = len(x)
  mean = (sum(x*y)/(sum(y)))
  offs = np.average(y)
  sigma = (sum((x-mean)**2.)/(n-1.))**0.5
  if offset == 1:
    popt,pcov = curve_fit(gauss2,x,y,p0=[np.amax(y),mean,sigma,offs])
  else:
    popt,pcov = curve_fit(gauss,x,y,p0=[np.amax(y),mean,sigma])
    '''
  print(" ")
  print("###########")
  print("Fit (Gauss) Peak Position : " + "%4.1f" % popt[1])
  print("Width                     : " + "%4.2f" % np.abs(popt[2] * 2*np.sqrt(2*np.log(2))))
  print("Amplitude                 : " + "%4.4f" % popt[0])
  print("Center of Mass            : " + "%4.1f" % mean)
  '''
  if offset == 1:
    print("Offset                    : " + "%4.4f" % popt[3])
    return gauss2(x, *popt), popt, pcov
  else:
    return gauss(x, *popt), popt, pcov

def fit(x,y,func,*params):
  popt,pcov = curve_fit(func,x,y,*params)
  return func(x, *popt), popt, pcov

def fit2(x,y,func,*params):
  popt,pcov = curve_fit(func,x,y,*params)
  return popt, pcov

def gauss(x,a,x0,sigma):
  """gaussian"""
  return a*exp(-(x-x0)**2./(2.*sigma**2.))

def gauss2(x,a,x0,sigma,const):
  """gaussian with offset"""
  return a*exp(-(x-x0)**2./(2.*sigma**2.)) + const

def convolve(x,y,func,params,resolution):
  """convolution of something, not sure if ever tested :/"""
  yachse = func((np.asarray(range(-10.*resolution,10.*resolution))/resolution),params)
  res = np.convolve(y,yachse)
  return res

def PV(x,mu,a,x0,width):
  """pseudo-voigt function"""
  return a*(((mu)*exp(-1.*math.log(2.,exp(1))*(x-x0)**2./(width**2.)))+((1.-mu)/(1.+(((x-x0)**2.)/(width**2.)))))


def PV2(x,mu,a,x0,width,const):
  """pseudo-voigt function with offset"""
  return a*(((mu)*exp(-1.*math.log(2.,exp(1))*(x-x0)**2./(width**2.)))+((1.-mu)/(1.+(((x-x0)**2.)/(width**2.))))) + const

def broaden(y,broaden):
  """broadening with a gaussian"""
  return np.convolve(y,gauss((np.asarray(range(-1000,1000))/100.),.004/broaden,0,broaden))[1000:-999]

def decay(x,amp,base,lifetime):
  return np.exp(-1*(x) / lifetime) * heaviside(x) * amp + base

def decay2(x,amp,base,lifetime,t0):
  """fails quite often"""
  return np.exp(-1*(x-t0) / lifetime) * heaviside(x,t0) * amp + base

def heaviside(x, x0 = 0):
  """basic heaviside function"""
  return np.piecewise(x, [x < x0, x >= x0], [0, 1])

###### GEORG 2018-02-13:
def fitgeorg(x,y,func,params,boundscontainer=None):
  print("bin in routine fitgeorg in fwhm")
  try:
    popt,pcov = curve_fit(func,x,y,params,sigma=None, absolute_sigma=True, check_finite=True,bounds=boundscontainer)
  except RuntimeError:
    print("Error - curve_fit failed")
  return func(x, *popt), popt, pcov

def PVgeorg(x,mu,a,x0,width,const):
  return a*(((mu)*exp(-1.*math.log(2.,exp(1))*(x-x0)**2./(width**2.)))+((1.-mu)/(1.+(((x-x0)**2.)/(width**2.)))))
