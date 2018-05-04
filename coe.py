import matplotlib.pyplot as plt
import tifffile as TF
import numpy as np
import glob
from scipy import asarray as ar, exp
from scipy.optimize import curve_fit
import scipy.stats
from pylab import *
import copy
import math

class sample:
#defines the class sample, which will later hold a list of all scans and scan relevant values, for example dispersion / zero position
  def __init__(self, crystals, path, roi, elist = None, samplename=None):
    self.name = samplename
    self.crystals = crystals
    self.path = path
    self.roi = roi
    self.elist = elist
#    print self.roi
  def getfilelist(self):
    filelist = glob.glob(self.path + "/*.tif")
    self.filelist = filelist
  def sumtif(self):
    sum = np.zeros([len(TF.imread(self.filelist[0])),len(TF.imread(self.filelist[0]).T)])
    i = 0
    for x in self.filelist:
      i += 1
      sum += TF.imread(x)
      if (i % 50 == 0):
        print i
    self.sum = sum
  def getspec(self):
    speclist = []
    for i in range(self.crystals):
      spec = copy.copy(self.sum[self.roi[i][3]])
      for j in range(self.roi[i][2],self.roi[i][3]):
        spec += self.sum[j]
      speclist.append(spec[self.roi[i][0]:self.roi[i][1]])
    self.speclist = speclist
  def getdisp(self):
#    self.disp = [[.4,6843]]
    self.disp = getdisp(self.crystals,self.roi,self.elist,self.path) #returns a,b of a linear fit
  def applydisp(self):
#    e_axis = np.zeros([self.crystals,len(self.speclist[0])])
    e_axis = np.zeros([self.crystals,487])   # each crystal gets own e axis
    for i in range(self.crystals):
      e_axis[i] = ar(range(487))*self.disp[i][0]+self.disp[i][1]
    self.e_axis = e_axis
  def plotsingle(self,norm=0,fit=0,sep=0):
    for i in range(self.crystals):
      if (sep == 1):
        figure()
      if (norm == 1):
        plt.plot(self.e_axis[i][self.roi[i][0]:self.roi[i][1]],self.normalizedspec[i])
        if (fit == 1):
          plt.plot(self.e_axis[i][self.roi[i][0]:self.roi[i][1]], tri_norm(self.e_axis[i][self.roi[i][0]:self.roi[i][1]], self.normfit[i].T[0][0], self.normfit[i].T[0][1], self.normfit[i].T[0][2], self.normfit[i].T[0][3], self.normfit[i].T[0][4], self.normfit[i].T[0][5], self.normfit[i].T[0][6], self.normfit[i].T[0][7], self.normfit[i].T[0][8], self.normfit[i].T[0][9]))
        if (fit == 2):
          plt.plot(self.e_axis[i][self.roi[i][0]:self.roi[i][1]], bi_ps(self.e_axis[i][self.roi[i][0]:self.roi[i][1]], self.pvfit[i].T[0][0], self.pvfit[i].T[0][1], self.pvfit[i].T[0][2], self.pvfit[i].T[0][3], self.pvfit[i].T[0][4], self.pvfit[i].T[0][5], self.pvfit[i].T[0][6]))
      else:
        plt.plot(self.e_axis[i][self.roi[i][0]:self.roi[i][1]],self.speclist[i])
    plt.show()
  def normalize(self,height,minenergy,maxenergy):
    normalizedspec = copy.copy(self.speclist)
    for i in range(self.crystals):
      normalizer = 0.
      for j in range(len(self.speclist[i])):
        if (self.e_axis[i][j] > minenergy) and (self.e_axis[i][j] < maxenergy):
          normalizer += self.speclist[i][j]
      normalizedspec[i] = self.speclist[i] / normalizer
    self.normalizedspec = normalizedspec
  def fitnorm(self):
    self.normfit = fitonnorm(self.e_axis,self.normalizedspec)
  def fitvoigt(self):
    self.pvfit = fitonbips(self.e_axis,self.normalizedspec)
  def sumspec(self, minenergy, maxenergy, steps):
    #interpolstates and sums up the single spectra
    interpar = np.zeros([2,steps])
    interpar[0] = np.linspace(minenergy, maxenergy, steps)
    for i in range(self.crystals):
      interpar[1] += np.interp(interpar[0],self.e_axis[i][self.roi[i][0]:self.roi[i][1]],self.speclist[i])
    self.summedspec = interpar

def readconfig(cfgfile):
  file = open(cfgfile,"r")
  elist = []
  elmo = 0
  for lines in file:
    a = lines[:-1].split(" ")
    if a[0] == "crystals":
      crystals = int(a[1])
    if a[0] == "roi":
      roi = []
      roiadd = [0,0,0,0]
      for i in range(len(a)-1):
        for j in range(4):
          roiadd[j] = int(a[i+1][1:-1].split(",")[j])
#          print j
#        print roiadd
        roi.append(copy.copy(roiadd))
    if elmo == 1:
      elist.append([a[0],int(a[1])])
    if a[0] == "[elastic]":
      elmo = 1
  file.close()
#  print roi
  return elist,crystals,roi

# 2016-09-26:
# das ist manuels rocksolid, aber gauss only:
def fwhm(x,y,amp=0,mean=0,width=1):
  x = ar(x)
  y = ar(y)
  n = len(x)
  if mean == 0:
    mean = x[np.argmax(y)]
  if width == 0:
    width = 1.
  if amp == 0:
    amp = np.amax(y)
  offset = np.average(y)
  def gaus(x,a,x0,width,const):
    return a*exp(-(x-x0)**2/(2*width**2)) + const
  popt,pcov = curve_fit(gaus,x,y,p0=[amp,mean,width,offset])
  return popt[1],(2*math.log(2,math.e))**.5*2*popt[2],popt[0]


def getdisp(crystals,roi,elist,path):
  disparray = np.zeros([len(elist),crystals,3])
  for i in range(len(elist)):
    elaspic = TF.imread(path + r"/elastic" + "/" + elist[i][0])
#    TF.imshow(elaspic)
#    plt.show()
    for j in range(crystals):
      spec = copy.copy(elaspic[roi[j][3]])
      for k in range(roi[j][2],roi[j][3]):
#        print 'Crystal number' + str(j)
        spec += elaspic[k]
      x, width, amp = fwhm(range(roi[j][0],roi[j][1]),spec[roi[j][0]:roi[j][1]])
#      print x, width, amp
#      plt.figure()
#      plt.plot(range(roi[j][0],roi[j][1]),spec[roi[j][0]:roi[j][1]])
#      plt.show()
      disparray[i][j][0] = elist[i][1]
      disparray[i][j][1] = x
      disparray[i][j][2] = amp
  dispersion = np.zeros([crystals,2])
  for i in range(crystals):
    print disparray.T[0][i]
    print disparray.T[1][i]
    print disparray.T[2][i]
    fit = np.polyfit(disparray.T[1][i],disparray.T[0][i],1) # linear fit
    print fit[0]
    dispersion[i][0] = fit[0]
    dispersion[i][1] = fit[1]
  return dispersion

def tri_norm(x, m1, m2, m3, s1, s2, s3, k1, k2, k3, c1):
  return k1*scipy.stats.norm.pdf(x, loc=m1 ,scale=s1) + k2*scipy.stats.norm.pdf(x, loc=m2 ,scale=s2) + k3*scipy.stats.norm.pdf(x, loc=m3 ,scale=s3) + c1

def voigt(x,a,x0,width):
  return a*((.8)*exp(-1.*math.log(2.,exp(1))*(x-x0)**2./(width**2.)))+(.2/(1.+(((x-x0)**2.)/(width**2.))))

def bi_ps(x,x01,x02,width1,width2,a1,a2,const):
  return voigt(x,a1,x01,width1) + voigt(x,a2,x02,width2) + const
  
def fitonnorm(xlist,ylist):
  print "hi"
  params = [6950, 6915, 6930, 1, 1, 1, 1, 1, 1, 0.01]
  fitarray = np.zeros([len(xlist),len(params),2])
  for spec in range(len(xlist)):
    fitted_params,covarmatrix = scipy.optimize.curve_fit(tri_norm,xlist[spec][0:400], ylist[spec], p0=params)
    for i in range(len(fitted_params)):
      print fitted_params[i], covarmatrix[i][i]
      fitarray[spec][i][0] = fitted_params[i]
      fitarray[spec][i][1] = covarmatrix[i][i]
  return fitarray

def fitonbips(xlist,ylist):
  print "bye"
  params = [6915, 6930, 1, 1, 1, 1, 0.01]
  fitarray = np.zeros([len(xlist),len(params),2])
  for spec in range(len(xlist)):
    fitted_params,covarmatrix = scipy.optimize.curve_fit(bi_ps,xlist[spec][0:400], ylist[spec], p0=params)
    for i in range(len(fitted_params)):
      print fitted_params[i], covarmatrix[i][i]
      fitarray[spec][i][0] = fitted_params[i]
      fitarray[spec][i][1] = covarmatrix[i][i]
  return fitarray

class spectrum:
  def __init__(self,x_axis,y_axis):
    self.x_axis = ar(x_axis)
    self.y_axis = ar(y_axis)
  def __add__(self, new):
    if (isinstance(new,spectrum) == True):
      if ((len(self.x_axis) == len(new.x_axis)) and ((self.x_axis-new.x_axis).any() == False)):
        return spectrum(self.x_axis,self.y_axis + new.y_axis)
      else:
        print "wrong axis"
        return 0.
    elif (type(new) == float or type(new) == int):
      return spectrum(self.x_axis,self.y_axis + new)
    else:
      print "wrong type"
      return 0.
  def __mul__(self, new):
    if (type(new) == float or type(new) == int):
      return spectrum(self.x_axis,self.y_axis * new)
    else:
      print "wrong type"
      return 0.
  def plot(self, newfigure=0):
    if newfigure == 1:
      figure()
    plt.plot(self.x_axis,self.y_axis)
    plt.show()
