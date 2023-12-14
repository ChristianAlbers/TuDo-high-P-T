import matplotlib.pyplot as plt
import numpy as np
from scipy import asarray as ar
from scipy.interpolate import InterpolatedUnivariateSpline
import itertools
import heapq
import copy
import fwhm
import os
#import statistics
#plt.ion()

"""
Defines the class "spectrum", a basic spectrum is an object which connects 2 arrays (interpreted as x and y-axis) to perform basic tasks.
In some cases there is a third axis, the statistical error. So far, there was no task where this became useful, so for now it is orphaned.
To create a spectrum use:
 i)   a = spectrum(listx,listy)            | listx and listy will be converted to numpy nd arrays (usually float64, sometimes int64)
		| listx and listy need to have the same lengths
 ii)  a = rsff("xy_filename.txt")        | col #0 = x-axis (optinal), col #1 (optinal) = y-axis
 iii) a = null()                    | empty spectrum ([],[]), a + null() = null() + a = a
 iv)  a = Sy(listy)                | x-axis = 0,...,n, where n is number of elements in y axis

List of operators     | a X b     a=(s)pect, b = (s)pect or (n)umber:
 &            | s & s
 +             | s + s, s + n
 -            | s - s, s - n
 *            | s * s, s * n
 /            | s / s, s / n
 >>            | s >> n
 <<            | s << n

List of functions returning a new spectra    | b = a.abc():
 ()
 b()
 b2()
 base()
 bg()
 bgr()
 c()
 cms()
 cms2()
 dev() == d()
 f()
 fit()
 fit_more()
 gauss()
 i()
 ln()
 log()
 measured()
 n()
 noise()
 pfit()
 rebin()
 reverse() == r()
 s()
 scale()
 sort()

List of other functions:
 ave()
 max()
 min()
 p()
 save()
 sum()
 swap()
"""

class spectrum:
	def __init__(self,x=[],y=[],background=[],error=[]):
		"""make spectrum from x and y list (and maybe errors)"""
		if x == [] and len(y) != 0:
			self.x = ar(range(len(y)))
			self.y = ar(y)
			#self.error = abs(ar(y))**.5
			return None
		if len(x) != len(y):
			print("\n------   not same length of x and y!   ------ \n")
			return None
		self.x = ar(x).astype(float)
		self.y = ar(y).astype(float)
		if len(error) == 0:
			self.error=np.zeros(len(y))
		else:
			self.error = ar(error).astype(float)
		if len(background) == 0:
			self.background=np.zeros(len(y))
		else:
			self.background = ar(background).astype(float)
	def __and__(self,new):
		"""putting all special operations in here, as this seems to be used so rarely, that it will not be used in two different meanings at once,
		replace how ever you need this..."""
		#1) """special sort for summing, does not interpolate, just adding extra points to the plot"""
		"""new_x_axis = ar(list(self.x) + list(new.x))
		new_y_axis = ar(list(self.y) + list(new.y))
		return (spectrum(new_x_axis,new_y_axis).sort())"""
		#2) cubic spline addition with adding all the data points on the x-axis but here we do cubic spline instead of linear interpolation
		"""defines + operator for addition"""
		if (isinstance(new,spectrum) == True):
			if len(new.x) == 0:
				return self
			elif len(self.x) == 0:
				return new
			new_x_axis = ar(list(self.x) + list(new.x))
			new_y_axis = ar(list(self.y) + list(new.y))
			return (spectrum(new_x_axis,new_y_axis).sort())
		elif (type(new) == float or type(new) == int or type(new) == np.float64 or type(new) == np.int or type(new) == np.int64):
			return spectrum(self.x,self.y + new, self.background, self.error)
	def __call__(self, min = 0, max = 0):
		"""defines the () function to cut the spectrum"""
		#extend to add values on side?
		if min == max:
			return spectrum(self.x,self.y,self.background,self.error)
		else:
			 return spectrum(self.x[np.searchsorted(self.x, min, side="left"):np.searchsorted(self.x, max, side="right")],self.y[np.searchsorted(self.x, min, side="left"):np.searchsorted(self.x, max, side="right")],self.background[np.searchsorted(self.x, min, side="left"):np.searchsorted(self.x, max, side="right")],self.error[np.searchsorted(self.x, min, side="left"):np.searchsorted(self.x, max, side="right")])
	def __repr__(self):
		"""defines how spectra are represented by just typing the instance name"""
		return "Spectrum with " + str(len(self.x)) + " datapoints"
		#return str(self.x) + "\n" + str(self.y)
	def __add__(self,new):
		"""defines + operator for addition"""
		if (isinstance(new,spectrum) == True):
			if len(new.x) == 0:
				return self
			elif len(self.x) == 0:
				return new#
			new_x=new.x
			return spectrum(new_x, np.interp(new_x, self.x, self.y, 0, 0) + np.interp(new_x, new.x, new.y, 0, 0), np.interp(new_x, self.x, self.background, 0, 0) + np.interp(new_x, new.x, new.background, 0, 0),np.sqrt(np.interp(new_x, self.x, self.error, 0, 0)**2 + np.interp(new_x, new.x, new.error, 0, 0)**2))
		elif (type(new) == float or type(new) == int or type(new) == np.float64 or type(new) == np.int or type(new) == np.int64):
			return spectrum(self.x,self.y + new, self.background + new, self.error)
	def __sub__(self, new):
		"""defines - operator for subtraction"""
		if (isinstance(new,spectrum) == True):
			new = spectrum(new.x,-1*new.y,new.background,new.error)
			#new=-1*new
			return self + new
		elif (type(new) == float or type(new) == int or type(new) == np.float64 or type(new) == np.int or type(new) == np.int64):
			return spectrum(self.x,self.y - new, self.background - new, self.error)
	def __mul__(self, new):
		"""defines * operator for multiplication"""
		if (isinstance(new,spectrum) == True):
			new_x = ar(list(k for k, v in itertools.groupby(heapq.merge(self.x, new.x))))
			new_y = np.zeros([len(new_x)])
			selfy = np.interp(new_x, self.x, self.y, 0, 0)
			newy = np.interp(new_x, new.x, new.y, 0, 0)
			for i in range(len(new_x)):
				new_y[i] = selfy[i] * newy[i]
			return spectrum(new_x, new_y)
		elif (type(new) == float or type(new) == int or type(new) == np.float64 or type(new) == np.int or type(new) == np.int64):
			return spectrum(self.x, self.y * new,self.background * new, self.error * new)
	def __truediv__(self, new):
		"""defines / operator for division"""
		if (isinstance(new,spectrum) == True):
			new_x = ar(list(k for k, v in itertools.groupby(heapq.merge(self.x, new.x))))
			new_y = np.zeros([len(new_x)])
			selfy = np.interp(new_x, self.x, self.y, 0, 0)
			newy = np.interp(new_x, new.x, new.y, 0, 0)
			for i in range(len(new_x)):
				new_y[i] = selfy[i] / newy[i]
			return spectrum(new_x, new_y)
		elif (type(new) == float or type(new) == int or type(new) == np.float64 or type(new) == np.int or type(new) == np.int64):
			return spectrum(self.x, self.y / new, self.background / new , self.error / new)
	def __lshift__(self, new):
		"""defines the << operator to shift x-axis to smaller values"""
		return spectrum(self.x - new, self.y, self.background, self.error)
	def __rshift__(self, new):
		"""defines the >> operator to shift x-axis to bigger values"""
		return spectrum(self.x + new,self.y, self.background, self.error)
	def sort(self):
		"""sorts the x values (y accordingly)"""
		order = np.argsort(self.x)
		return spectrum(self.x[order],self.y[order], self.background[order], self.error[order])
	def p(self, newfigure=0, newlabel = "", style = "", ng = 0,log=0):
		"""plots the spectrum"""
		if newfigure == 1:
			plt.figure()
		if ng == 1:
			plt.errorbar(self.x,self.y,self.error, label = newlabel)
			#plt.legend()
		else:
			plt.plot(self.x,self.y,label = newlabel)
	def swap(self):
		"""swaps x and y of self instance (1)"""
		temp = copy.deepcopy(self.x)
		self.x = copy.deepcopy(self.y)
		self.y = copy.deepcopy(temp)
	def s(self):
		"""returns new spectrum with x and y swapped"""
		return spectrum(self.y,self.x,self.background,self.error)
	def sum(self, min = 0, max = 0):
		"""sums up the y axis"""
		summe = 0.
		if min == max:
			for i in self.y:
				summe += i
		else:
			for i in range(len(self.y)):
				if (self.x[i] >= min) and (self.x[i] <= max):
					summe += self.y[i]
		return summe
	def ave(self, min = 0, max = 0):
		"""integrates"""
		return self.i(min, max)
	def avg(self, min = 0, max = 0):
		"""integrates, same as ave()"""
		return self.i(min, max)
	def i(self, min = 0, max = 0):
		"""integrates"""
		#maybe interpolate/extrapolate to real min / max for more accurate calculation, maybe not
		average = 0.
		if min == max:
			for i in range(len(self.y)-1):
				average += (self.y[i] + self.y[i+1]) * (.5*(self.x[i+1] - self.x[i]))
		else:
			for i in range(len(self.y)-1):
				if (self.x[i] >= min) and (self.x[i] <= max):
					average += (self.y[i] + self.y[i+1]) * (.5*(self.x[i+1] - self.x[i]))
		return float(average)
	def F(self,const=0):
		"""returns antiderivative (such an ugly word for Stammfunktion)"""
		intarray = np.zeros([2,len(self.x)-1])
		for a in range(len(self.x)-1):
			intarray[0][a] = (self.x[a]+self.x[a+1])/2.
			if a == 0:
				intarray[1][a] = (self.y[a+1]+self.y[a])*((self.x[a+1]-self.x[a]))/2. + const
			else:
				intarray[1][a] = (self.y[a+1]+self.y[a])*((self.x[a+1]-self.x[a]))/2. + intarray[1][a-1]
		return spectrum(intarray[0],intarray[1])
	def max(self):
		"""returns the maximum"""
		return np.max(self.y)
	def min(self):
		"""returns the maximum"""
		return np.min(self.y)
	def n(self, min = 0, max = 0):
		"""normalize to area"""
		return spectrum(self.x,self.y/self.ave(min, max),self.background/self.ave(min, max),self.error/self.ave(min, max))
	def n_minus(self, min = 0, max = 0):
		"""normalize to negative area if x axis is wrong"""
		return spectrum(self.x,self.y/(-self.ave(min, max)),self.background/self.ave(min, max),self.error/self.ave(min, max))
	def rebin(self,rebin_val):
		"""redefine x and y axis
		can be done with new x-axis or a certain number of data points (int/float)"""
		if (type(rebin_val) == np.ndarray):
			return spectrum(rebin_val, np.interp(rebin_val, self.x, self.y, 0, 0))
		elif (type(rebin_val) == float or type(rebin_val) == int):
			new_x = np.linspace(min(self.x),max(self.x),rebin_val - 1)
			return spectrum(new_x, np.interp(new_x, self.x, self.y, 0, 0), np.interp(new_x, self.x, self.background, 0, 0), np.interp(new_x, self.x, self.error, 0, 0))
	def reverse(self):
		"""see r()"""
		return spectrum(self.x[::-1],self.y[::-1],self.background[::-1], self.error[::-1])
	def r(self):
		"""change order of x and y entries"""
		return spectrum(self.x[::-1],self.y[::-1],self.background[::-1], self.error[::-1])
	def gauss(self,offset=0,fullinfo = 0):
		"""fit gaussian, with or without offset"""
		fit_result, fit_params, fit_error = fwhm.main(self.x,self.y,offset)
		fit_spectrum = spectrum(self.x,fit_result)
		if fullinfo != 0:
			fit_spectrum.fit_p = fit_params
			fit_spectrum.fit_e = fit_error
			fit_spectrum.fit_func = offset
			fit_spectrum.fit_flag = 3
		return fit_spectrum
	def fit(self,func,params,fullinfo = 0):
		"""fit any function"""
		fit_result, fit_params, fit_error = fwhm.fit(self.x,self.y,func,params)
		fit_spectrum = spectrum(self.x,fit_result)
		if fullinfo != 0:
			fit_spectrum.fit_p = fit_params
			fit_spectrum.fit_e = fit_error
			fit_spectrum.fit_func = func
			fit_spectrum.fit_flag = 2
		return fit_spectrum
	def pfit(self,order=0,fullinfo = 0):
		"""fit polynom"""
		fit_params = np.polyfit(self.x,self.y,order)
		fit_result = np.zeros([len(self.x)])
		for i in range(len(fit_params)):
			fit_result += (self.x**i) * fit_params[len(fit_params)-(i+1):len(fit_params)-i]
		fit_spectrum = spectrum(self.x,fit_result)
		if fullinfo != 0:
			fit_spectrum.fit_p = fit_params
			fit_spectrum.fit_func = order
			fit_spectrum.fit_flag = 1
		return fit_spectrum
	def dev(self):
		"""see d()"""
		return self.d()
	def d(self):
		"""get derivative"""
		derivarray = np.zeros([2,len(self.x)-1])
		for a in range(len(self.x)-1):
			derivarray[0][a] = (self.x[a]+self.x[a+1])/2.
			derivarray[1][a] = (self.y[a+1]-self.y[a])/((self.x[a+1]-self.x[a]))
		return spectrum(derivarray[0],derivarray[1])
	def save(self,savefilename,force=False):
		"""save as x/y file"""
		if os.path.exists(savefilename) == True and force == False:
			overwritefile = input("output file exists, overwrite? (y/n)")
			if overwritefile != "y" and overwritefile != "Y":
				print("aborted")
				return
		specfile = open(savefilename,"w")
		for i in range(len(self.x)):
			specfile.write(str(self.x[i]) + " " + str(self.y[i]) + " " + str(self.background[i]) + " " + str(self.error[i]))
			specfile.write("\n")
		specfile.close()
	def b(self, broaden):
		"""convolute with gaussian
			 Useage: assuming equal x-spacing: fwhm(b(x)) = x-space * x * 100 * sqrt(2*ln(2))
			 So for a 1 step / eV broadening apply b(b_x) with b_x = 1. / (100 * (a.x[1] - a.x[0]) * 2.3548)
		"""
		#guassian is just defined for 1000 points, so too big FWHM values may result in strange results
		return spectrum(self.x,fwhm.broaden(self.y,broaden))
	def b2(self, broaden):
		"""convolute with gaussian with fwhm = 1 unit on x-axis
			 Only works for spectra with equal spaced x-axis (!)
		"""
		return spectrum(self.x,fwhm.broaden(self.y,broaden / (100. * np.abs(self.x[1] - self.x[0]) * 2 * np.sqrt(2 * np.log(2)))))
	def f(self):
		"""get fourier transform, looks weird"""
		return spectrum(self.x,abs(np.fft.fft(self.y))/len(self.x))
	def c(self,func,params,resolution):
		"""convolve with any function"""
		return fwhm.convolve(self.x,self.y,func,params,resolution) #x-achse berechnen (iwann mal)
	def ln(self):
		"""ln of y-axis"""
		return spectrum(self.x,np.log(self.y))
	def log(self,base):
		"""log with base n of y-axis"""
		return spectrum(self.x,np.log(self.y))/np.log(base)
	def bg(self,rois=0,order=0,func=0,params=[],fullinfo = 0):
		"""fits something (could be background) to a set of rois"""
		if rois == 0:
			rois = [[np.min(self.x),np.max(self.x)]]
		fitspectr = null()
		for roi in rois:
			fitspectr = fitspectr&self(roi[0],roi[1])
		if func == 0 and params == []:
			fit_params = np.polyfit(fitspectr.x,fitspectr.y,order)
			fit_result = np.zeros([len(self.x)])
			for i in range(len(fit_params)):
				fit_result += (self.x**i) * fit_params[len(fit_params)-(i+1):len(fit_params)-i]
			fit_spectrum = spectrum(self.x,fit_result)
			if fullinfo != 0:
				fit_spectrum.fit_p = fit_result
			return fit_spectrum
		else:
			fit_result = fwhm.fit(fitspectr.x,fitspectr.y,func,params)[0]
			fit_spectrum = spectrum(self.x, func(self.x, *fit_result))
			if fullinfo != 0:
				fit_spectrum.fit_p = fit_result
			return fit_spectrum
	def bgr(self,rois=0,order=0,func=0,params=[]):
		"""subtracts the background, see self.bg()"""
		return self - self.bg(rois,order,func,params)
	def cms(self):
		"""find cms in current data"""
		return (sum(self.x*self.y)/(sum(self.y)))
	def cms2(self):
		"""find cms in current data, nother method"""
		return np.trapz(self.y*self.x,self.x)/np.trapz(self.y,self.x)
	def measured(self,newmax=0,offset=0):
		"""scales to 0, max and then throws out a discrete distribution, lets see"""
		if newmax != 0:
			results = self.scale(newmax)
			results = spectrum(results.x,statistics.poisson(results.y))
		else:
			results = spectrum(self.x,statistics.poisson(self.y))
		if offset != 0:
			results += results.noise(offset)
		return results
	def noise(self,offset):
		"""generates a constant background, statistically distributed"""
		return spectrum(self.x,statistics.poisson2(offset,len(self.x)))
	def scale(self,newmax):
		"""scales the values to a new maximum"""
		return self/self.max()*newmax
	def base(self,neg=0):
		"""sets lowest / highest (neg=1) value to 0"""
		if neg == 1:
			return self-self.max()
		else:
			return self-self.min()
	def fit_more(self,new_x):
		new = self.rebin(new_x)
		if self.fit_flag == 1:
			new.y = np.zeros([len(new.x)])
			for i in range(len(self.fit_p)):
				new.y += (new.x**i) * self.fit_p[len(self.fit_p)-(i+1):len(self.fit_p)-i]
		elif self.fit_flag == 2:
			new.y = self.fit_func(new.x,*self.fit_p)
		elif self.fit_flag == 3:
			if self.fit_func == 0:
				new.y = fwhm.gauss(new.x, self.fit_p[0], self.fit_p[1], self.fit_p[2])
			else:
				new.y = fwhm.gauss2(new.x, self.fit_p[0], self.fit_p[1], self.fit_p[2], self.fit_p[3])
		return new
	def redo_x_axis(self,start=0.,end=1.):
		new = self.sort()
		new.x = new.x * (end-start) / (new.x[-1::] - new.x[0])
		new.x = new.x - (new.x[0] - start)
		return new
	def clean(self,value_max):
		new_x = self.x
		new_y = self.y
		while np.max(new_y) > value_max:
			new_x = np.delete(new_x,new_y.argmax())
			new_y = np.delete(new_y,new_y.argmax())
		return spectrum(new_x,new_y)
	def clean2(self,value_max):
		"""remove values by looking at derivative, just quickly done"""
		todo = self()
		while np.max(todo.d().y) > value_max:
			dev_y = todo.d().y
			todo.x = np.delete(todo.x,dev_y.argmax()+1)
			todo.y = np.delete(todo.y,dev_y.argmax()+1)
		return spectrum(todo.x,todo.y)
	def energ(self, slope, ynull):
		"""scale x-axis to eV"""
		new_x = self.x * slope + ynull
		new_y = self.y
		return spectrum(new_x, new_y)

def rsff(specfile,x_axis = 0, y_axis = 1):
	"""reads spectrum from (x y) file"""
	if os.path.exists(specfile) == False:
		print("file not found")
		return
	spectfile = open(specfile,"r")
	X_A = []
	Y_A = []
	while True:
		a = spectfile.readline()
		if a == "":
			break
		else:
			try:
				b = a.split()
				X_A.append(float(b[x_axis]))
				Y_A.append(float(b[y_axis]))
			except:
				break
	print("read spectrum with", len(X_A), "data points")
	toreturn = spectrum(X_A,Y_A)
	toreturn.source = specfile
	return toreturn

def null(start = 0, end = 0, tot = 0):
	"""returns an empty spectrum"""
	return spectrum(np.linspace(start,end,tot),np.zeros([tot]))

def Sy(y):
	"""gives spectrum with x axis = range(len(x))"""
	return spectrum([],y)

def funct(start, end, tot, func, params):
	return spectrum(np.linspace(start,end,tot), func(np.linspace(start,end,tot),*params))

#for lazy ppl:
def x(start, end, tot):
	"""gives an array, supposed to be used as an x-axis"""
	return np.linspace(start,end,tot+1)

def y(inp):
	"""gives spectrum with x axis = range(len(x))"""
	return spectrum([],inp)

def xy(inp1,inp2):
	"""gives spectrum just like spectrum class"""
	return spectrum(inp1,inp2)

"""
some plotting options for nicer graphs
def ng():
	plt.axis((7090,7120,.3,0.9))
	plt.xlabel("Energy [eV]",fontsize=12)
	plt.ylabel("Counts (arb. units)",fontsize=10)
	plt.tick_params("x",labelsize = 12)
	plt.tick_params("y",labelsize = 10)
	plt.title("Fe K-beta (VtC)",fontsize=20)
	spect.plt.legend()
	spect.plt.grid()
"""
