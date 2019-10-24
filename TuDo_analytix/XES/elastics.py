import numpy as np
import matplotlib.pyplot as plt 
from skimage.feature import peak_local_max
import scipy.ndimage as ndimage
import glob
import os
import fabio
from spect import *
from functions import *
from scipy.optimize import curve_fit
plt.ion()

class XES_elastics:
	def __init__(self,
				images = [], #list of paths to pictures. Sorted!!!
				energies = [], #e.g. [10500,10550,10600,10650,10700]
				reflex_fac = 1.5, #for 660 to 440
				crys_number=4
				):
		self.mask=[0]
		self.images= dict()
		if any(images):
			self.load_images(images)
		self.energies=energies
		self.reflex_fac=reflex_fac
		self.params=np.zeros((crys_number,2)) # [[m], [b]] for y=m*x+b
		self.params[:,0]=1
		self.crys_number=crys_number
		#self.rois

	def load_images(self,
					images):
		if not any(self.images):
			self.mask=np.zeros(fabio.open(images[0]).data.shape)
			self.mask[:]=1
		for file in images:
			self.images[file]=fabio.open(file)

	def calc_params(self,
					energies=False,
					reflex_fac= False,
					mode="auto"
					):
		xpos=np.linspace(0,self.mask.shape[1],self.mask.shape[1])
		el_params=np.zeros((self.crys_number,2))
		if not energies:
			ener = self.energies
		else:
			ener=energies
		if not reflex_fac:
			fac = self.reflex_fac
		else:
			fac=reflex_fac
		if not any(self.images.keys()):
			print("----- No images to calculate energies -----")
			return None
		elif len(self.images.keys()) != len(ener):
			print("----- number of images and energies no the same -----")
			return None
		else:
			crystals=np.zeros((self.crys_number,self.mask.shape[1]))
			peaks=np.zeros((self.crys_number,len(self.images)))
			imagecount=0
			for image in self.images:
				data=self.images[image].data
				data=data*self.mask
				if mode == "auto":
					rois=self.peakfinder(image=self.images[image],num_rois=self.crys_number)
				else:
					print("choose a mode to calculate: auto, rois tbe")
				roiscount=0
				for roi in rois:
					crystal=spectrum(xpos,np.array(np.sum(data[roi[0]-5:roi[0]+5],axis = 0)))
					peaks[roiscount,imagecount]=crystal(roi[1]-5,roi[1]+5).gauss(fullinfo=1).fit_p[1]
					roiscount+=1
				imagecount+=1
			for r in range(0,len(peaks)):
				popt=curve_fit(lin,peaks[r],ener)
				el_params[r]=popt[0]/fac
			self.params=el_params

	def peakfinder(self,
					image,
					num_rois=4,
					distance = 20): #finds elastics automatically by local
								#maxima
		im=image.data
		peaks = peak_local_max(im, min_distance=distance, num_peaks=num_rois)
		peaks = sorted(peaks ,key=lambda x: x[0])
		return peaks

	def add_mask(self, mask):
		#masks an area to avoid braggs etc. mask has the form [ x1, x2, y1, y2]
		self.mask[mask[0]:mask[1],mask[2]:mask[3]]=0