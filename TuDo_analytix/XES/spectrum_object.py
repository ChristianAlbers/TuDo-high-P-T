import numpy as np
import matplotlib.pyplot as plt 
from skimage.feature import peak_local_max
import scipy.ndimage as ndimage
import glob
import os
import fabio
from spect import *
from .elastics import *

'''
ToDo:
- allow dynamic shift of background (see calc_spectrum)
- function to set a mask (gui?)
- dictionary with reference spectra and a function to fit them
- dictionary with all spectra for the single images
- ...
'''
class XES_object:
	def __init__(self,
				images = [], #list of paths to pictures
				elastic = XES_elastics(), #XES_elastics object
				crys_number=4 #number of your crystals
				):
		self.mask=[0]
		self.det_image=0
		self.det_images = dict()
		if any(images):
			self.load_images(images)
			self.sum_images()
		self.spectrum = spectrum()
		self.spect_by_image= spectrum()
		self.bgr_by_image= spectrum()
		self.params = elastic.params
		self.crys_number=crys_number
		self.elastic = elastic
		self.rois = np.zeros((crys_number,2)).astype(int)

	def load_images(self,
					images):
		if not any(self.det_images):
			self.mask=np.zeros(fabio.open(images[0]).data.shape)
			self.mask[:]=1
		for file in images:
			self.det_images[file]=fabio.open(file)

	def sum_images(self):
		if self.det_image == 0 and any(self.det_images):
			self.det_image = np.zeros(self.det_images[list(self.det_images.keys())[0]].data.shape)
			for image in self.det_images:
				self.det_image += self.det_images[image].data
		else:
			print("----- No images to sum up -----")

	def auto_rois(self, distance = 10, roiwidth = 5):
		im=self.det_image
		peaks = peak_local_max(im, min_distance=distance,
								num_peaks=self.crys_number)
		peaks = sorted(peaks ,key=lambda x: x[0])
		rois=np.zeros((self.crys_number,2)).astype(int)
		for i in range(0,len(peaks)):
			rois[i]=[peaks[i][0]-roiwidth,peaks[i][0]+roiwidth]
		self.rois = rois

	def calc_spectrum(self, image, roi, params,
						roiwidth = 5, bgrwidth=3, bgrfit=False):
		data = image.data
		xpos=np.linspace(0,self.mask.shape[1],self.mask.shape[1])
		x = xpos*params[0] + params[1]
		if params[0] < 0:
			crystal=spectrum(x,
							np.sum(data[roi[0]-roiwidth:roi[1]+roiwidth],axis = 0)
							).reverse()
			bgr=spectrum(x,
						np.sum(data[roi[0]-roiwidth-bgrwidth:roi[0]-roiwidth],axis = 0)
						+np.sum(data[roi[1]+roiwidth:roi[1]+roiwidth+bgrwidth],axis = 0)
						).reverse()
		else:
			crystal=spectrum(x,np.sum(data[roi[0]-roiwidth:roi[1]+roiwidth],axis = 0)
							)
			bgr=spectrum(x,
						np.sum(data[roi[0]-roiwidth-bgrwidth:roi[0]-roiwidth],axis = 0)
						+np.sum(data[roi[1]+roiwidth:roi[1]+roiwidth+bgrwidth],axis = 0)
						)
		if bgrwidth > 0:
			bgr = bgr * (crystal(crystal.x[5],crystal.x[55]).i() / bgr(bgr.x[5],bgr.x[55]).i())
			#needs to be dynamic range!
		if bgrfit:
			bgr=bgr.bg([(bgr.x[2],bgr.x[-2])],order=5)
		return crystal, bgr

	def calc_spectrum_all(self, roiwidth=5, bgrwidth=3,
							bgrfit=False, show=False):
		self.spectrum=spectrum()
		crystals=dict()
		background=dict()
		for i in range(0,self.crys_number):
			crystals[i] = spectrum()
			background[i] = spectrum()

		for image in self.det_images:
			r=0
			for roi in self.rois:
				crys, bgr = self.calc_spectrum(self.det_images[image],
												roi, self.params[r],
												roiwidth,bgrwidth,bgrfit)
				crystals[r] += crys
				background[r] += bgr
				r+=1

		for crystal in crystals:
			self.spectrum += crystals[crystal]
			if show:
				crystals[crystal].p(1)
				background[crystal].p()
		if show:
			self.spectrum.p(1)

	def check_all_images(self, roiwidth=5, bgrwidth=3):
		crys=spectrum()
		bgr=spectrum()
		for image in self.det_images:
			fig, ax = plt.subplots(2,1)
			ax[0].imshow(self.det_images[image].data)
			fig.suptitle(image)
			r=0
			for roi in self.rois:
				crys, bgr = self.calc_spectrum(self.det_images[image],
												roi, self.params[r],
												roiwidth,bgrwidth)
				ax[0].axhspan(roi[0],roi[1],color="red",alpha=0.25)
				ax[1].plot(crys.x,crys.y,label=roi)
				r+=1
			plt.legend()

	def show_rois(self):
		plt.figure()
		plt.imshow(self.det_image)
		for roi in self.rois:
			plt.axhspan(roi[0],roi[1],color="red",alpha=0.25)
