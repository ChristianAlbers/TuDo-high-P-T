import numpy as np
import matplotlib.pyplot as plt 
from skimage.feature import peak_local_max
import scipy.ndimage as ndimage
import glob
import os
import fabio
from spect import *
from .elastics import *
from matplotlib.widgets import Button
from scipy.optimize import curve_fit
from lmfit import Model
'''
ToDo:
- function to set a mask (gui?)
- more dynamic reference fit (+ linear function for example) 
- ...
'''
class XES_object:
	def __init__(self,
				images = [], #list of paths to pictures
				elastic = XES_elastics(), #XES_elastics object
				crys_number=4 #number of your crystals
				):
		self.mask=[0] #masks an area on detectorimage (used for braggs etc)
		self.det_image=0 #summed detector image
		self.det_images = dict() #list of all detector images
		if any(images):
			self.load_images(images)
			self.sum_images()
		self.spectrum = spectrum() #spectrum of all summed images
		self.spectrum_temp = spectrum() #temporarily spectrum for fits
		self.spect_by_image= dict() #single spectra for all images
		self.bgr_by_image= dict() #background for all images
		self.params = elastic.params #parameters for calibration of detectorimages
		self.crys_number=crys_number #number of crystals used
		self.elastic = elastic #elastic object for calibration
		self.rois = np.zeros((crys_number,2),dtype=int)
		self.references = dict() #reference spectra
		self.references_temp = dict() #reference spectra only for fitting
		self.fit_refs=0 #marks the references used for fitting

	def load_images(self,
					images):
		if not any(self.det_images):
			self.mask=np.zeros(fabio.open(images[0]).data.shape)
			self.mask[:]=1
		for file in images:
			file=os.path.normpath(file)
			self.det_images[file]=fabio.open(file)

	def sum_images(self):
		if any(self.det_images):
			self.det_image = np.zeros(self.det_images[list(self.det_images.keys())[0]].data.shape)
			for image in self.det_images:
				self.det_image += self.det_images[image].data
		else:
			print("----- No images to sum up -----")
	'''
	def auto_rois(self, distance = 10, roiwidth = 5):
		#set automatic rois to the whole detector image
		im=self.det_image
		peaks = peak_local_max(im, min_distance=distance,
								num_peaks=self.crys_number)
		peaks = sorted(peaks ,key=lambda x: x[0])
		rois=np.zeros((self.crys_number,2)).astype(int)
		for i in range(0,len(peaks)):
			rois[i]=[peaks[i][0]-roiwidth,peaks[i][0]+roiwidth]
		self.rois = rois
	'''
	def auto_rois(self, roiwidth = 5, show=False):
		img=self.det_image
		out_rois=np.zeros((self.crys_number,2), dtype=int)	#out_rois[0][:]=[y0,y1] etc.
		y_median=np.percentile(img,q=98,axis=1)
		if show:
			plt.figure(999)
			plt.imshow(img,vmax=np.percentile(img,q=99.9))

		broaden=0.01
		y_spec=spectrum(range(0,img.shape[0]),y_median)
		y_peaks, y_properties = find_peaks(y_spec.b(broaden).y)
		while len(y_peaks) > self.crys_number:
			broaden+=0.01
			y_peaks, y_properties = find_peaks(y_spec.b(broaden).y)

		if show:
			plt.figure()
			y_spec.p()
			print(y_peaks)
		for i, peak in enumerate(y_peaks):
			out_rois[i,0]=int(peak-roiwidth) 	#y0
			out_rois[i,1]=int(peak+roiwidth)	#y1
			if show:
				plt.figure(999)
				plt.axhspan(peak-roiwidth,peak+roiwidth,color="red",alpha=0.25)
		self.rois = out_rois

	def shift_rois(self,shift,show=False):
		#shifts rois, 
		if len(shift) != len(self.rois):
			print("----- length of shift and rois don't match! -----")
			return None
		for i in range(0,len(shift)):
			self.rois[i,0] += shift[i]
			self.rois[i,1] += shift[i]
		if show:
			plt.close(999)
			plt.figure(999)
			plt.imshow(self.det_image)
			for roi in self.rois:
				plt.axhspan(roi[0],roi[1],color="red",alpha=0.25)

	def delete_roi(self,roinumber):
		self.rois=np.delete(self.rois,roinumber,0)
		self.params=np.delete(self.params,roinumber,0)
		self.crys_number-=1

	'''
	def calc_spectrum(self, image, roi, params, bgrwidth=3,
							bgr_range=(7120,7130),bgrfit=False):
		#calculate spectrum only for one image and roi
		#and return crystal spectrum and background 
		data = image.data
		xpos=np.linspace(0,self.mask.shape[1]-1,self.mask.shape[1])
		x = xpos*params[0] + params[1]
		if params[0] < 0:
			crystal=spectrum(x,
							np.sum(data[roi[0]:roi[1]],axis = 0)
							).reverse()
			bgr=spectrum(x,
						np.sum(data[roi[0]-bgrwidth:roi[0]],axis = 0)
						+np.sum(data[roi[1]:roi[1]+bgrwidth],axis = 0)
						).reverse()
		else:
			crystal=spectrum(x,np.sum(data[roi[0]:roi[1]],axis = 0)
							)
			bgr=spectrum(x,
						np.sum(data[roi[0]-bgrwidth:roi[0]],axis = 0)
						+np.sum(data[roi[1]:roi[1]+bgrwidth],axis = 0)
						)

		if bgrwidth > 0:
			bgr = bgr * (crystal(bgr_range[0],bgr_range[1]).i() / bgr(bgr_range[0],bgr_range[1]).i())

		if bgrfit:
			bgr=bgr.bg([(7010,7130)],order=11)
		return crystal, bgr
	def calc_spectrum_all(self, bgrwidth=3, bgr_range=(7115,7120),bgrfit=False, show=False):
		#calculates spectrum from all images
		#and interpolates it with 0.1eV stepsize to total values
		self.spectrum=spectrum()
		crystals=dict()
		background=dict()
		for i in range(0,self.crys_number):
			crystals[i] = spectrum()
			background[i] = spectrum()
		for image in self.det_images:
			self.spect_by_image[image]=spectrum()
			self.bgr_by_image[image]=spectrum()
			r=0
			for roi in self.rois:
				crys, bgr = self.calc_spectrum(self.det_images[image],
												roi, self.params[r],
												bgrwidth,bgr_range,bgrfit)
				crystals[r] += crys
				background[r] += bgr
				self.spect_by_image[image]+=crys
				self.bgr_by_image[image]+=bgr
				r+=1
		if show:
			plt.figure()
		for crystal in crystals:
			self.spectrum += crystals[crystal]
			self.spectrum -= background[crystal]
			if show:
				crystals[crystal].p(newlabel="crystal: " + str(crystal))
				background[crystal].p()
				((crystals[crystal]-background[crystal])/(crystals[crystal]-background[crystal])(7000,7100).max()).p(newlabel="crystal: " + str(crystal))
				
			plt.legend()
		x_low=np.round(self.spectrum.x[0])
		x_high=np.round(self.spectrum.x[-1])

		self.spectrum = spectrum(np.linspace(x_low,x_high,(x_high-x_low)*10+1),
								np.interp(
								np.linspace(x_low,x_high,(x_high-x_low)*10+1),
								self.spectrum(x_low,x_high).x,
								self.spectrum(x_low,x_high).y))

		if show:
			self.spectrum.p(1)
	'''
	def calc_spectrum_all(self, bgrwidth=3, bgr_range=(0,50),bgrfit=10, show=False):
		data = self.det_image
		xpos=np.linspace(0,self.mask.shape[1]-1,self.mask.shape[1])
		params=self.params
		if type(bgrfit) != int:
			bgrfit=10
			print("please give bgrfit = int, took bgrfit = 10")
		#x = xpos*params[0] + params[1]
		fig=plt.gcf().number
		for file in self.det_images:
			file=os.path.normpath(file)
			img=fabio.open(file).data
			#print("hier")
			#print(file)
			spect=spectrum()
			crystals=np.zeros((4,self.mask.shape[1]))
			bgrs=np.zeros((4,self.mask.shape[1]))
			for r, roi in enumerate(self.rois):
				crystals[r]=np.array(np.sum(img[roi[0]:roi[1]], axis = 0))
				bgrs[r]=np.array(np.sum(img[roi[0]-bgrwidth:roi[0]], axis = 0))+np.array(np.sum(img[roi[1]+bgrwidth:roi[1]], axis = 0))
				bgrs[r]=bgrs[r]*(1/np.sum(bgrs[r][bgr_range[0]:bgr_range[1]])*np.sum(crystals[r][bgr_range[0]:bgr_range[1]]))
				if bgrfit:
					#bgrs[r]=spectrum(xpos,bgrs[r]).pfit(bgrfit).y
					bgrs[r]=spectrum(xpos,bgrs[r]).bg([(xpos[10],xpos[-10])],bgrfit).y

				if params[r,0] < 0:
					crystal=spectrum(xpos*params[r,0]+params[r,1],crystals[r]-bgrs[r],bgrs[r]).reverse()
				else:
					crystal=spectrum(xpos*params[r,0]+params[r,1],crystals[r]-bgrs[r],bgrs[r])
					plt.legend()
				spect+=crystal
			self.spect_by_image[file]=spect
		crystals=np.zeros((4,self.mask.shape[1]))
		bgrs=np.zeros((4,self.mask.shape[1]))
		spect=spectrum()
		for r, roi in enumerate(self.rois):
			crystals[r]=np.array(np.sum(data[roi[0]:roi[1]], axis = 0))
			bgrs[r]=np.array(np.sum(data[roi[0]-bgrwidth:roi[0]], axis = 0))+np.array(np.sum(data[roi[1]+bgrwidth:roi[1]], axis = 0))
			bgrs[r]=bgrs[r]*(1/np.sum(bgrs[r][bgr_range[0]:bgr_range[1]])*np.sum(crystals[r][bgr_range[0]:bgr_range[1]]))
			if bgrfit:
					#bgrs[r]=spectrum(xpos,bgrs[r]).pfit(bgrfit).y
					bgrs[r]=spectrum(xpos,bgrs[r]).bg([(xpos[10],xpos[-10])],bgrfit).y
			if params[r,0] < 0:
				crystal=spectrum(xpos*params[r,0]+params[r,1],crystals[r]-bgrs[r],bgrs[r]).reverse()
			else:
				crystal=spectrum(xpos*params[r,0]+params[r,1],crystals[r]-bgrs[r],bgrs[r])
			if show:
				plt.figure(fig+1)
				plt.plot(xpos,crystals[r],label=r)
				plt.plot(xpos,bgrs[r])
				plt.legend()

				plt.figure(fig+2)
				(crystal/crystal.max()).p(newlabel=r)
				plt.legend()
			spect+=crystal
		#self.spectrum=spectrum(np.linspace(7020,7120,1001),np.interp(np.linspace(7020,7120,1001),spect(7020,7120).x,spect(7020,7120).n().y))
		self.spectrum=spect


	def check_all_images(self, roiwidth=5, bgrwidth=3):
		#plot all detector images with rois
		#and generated spectra for each crystal
		crys=spectrum()
		bgr=spectrum()
		for image in self.det_images:
			fig, ax = plt.subplots(2,1)
			ax[0].imshow(self.det_images[image].data)
			fig.suptitle(image)
			r=0
			for roi in self.rois:
				crys, bgr = self.calc_spectrum(self.det_images[image],roi, 
												self.params[r],bgrwidth)
				ax[0].axhspan(roi[0],roi[1],color="red",alpha=0.25)
				ax[1].plot(crys.x,crys.y,label=roi)
				r+=1
			plt.legend()

	def show_rois(self):
		plt.figure()
		plt.imshow(self.det_image)
		for roi in self.rois:
			plt.axhspan(roi[0],roi[1],color="red",alpha=0.25)

	def add_reference(self, spect, name):
		#adds reference spectrum to reference dictionary
		#and interpolates to the same values as the spectrum
		ref=spect(self.spectrum.x[0],self.spectrum.x[-1]).n()
		x_low=np.round(self.spectrum.x[0])
		x_high=np.round(self.spectrum.x[-1])
		ref=spectrum(np.linspace(x_low,x_high,(x_high-x_low)*10+1),
					np.interp(np.linspace(x_low,x_high,(x_high-x_low)*10+1),
					ref(x_low,x_high).x,
					ref(x_low,x_high).n().y))
		self.references[name]=ref

	def show_references(self):
		#plot all references
		plt.figure()
		for ref in self.references:
			self.references[ref].p(newlabel=ref)
		plt.legend()

	def reference_function(self, x, *params):
		#returns function with reference spectra
		#selected by fit_refs ([1,0,1,1 etc])
		func=0
		for i in range(0,len(self.references_temp)):
			key = list(self.references_temp.keys())[i]
			func += params[i] * self.fit_refs[i] * self.references_temp[key].y
		func+=params[-1]
		return func

	def reference_function2(self, x, *params):
		#returns function with reference spectra
		#selected by fit_refs ([1,0,1,1 etc])
		#m*x + b
		func=0
		for i in range(0,len(self.references_temp)):
			key = list(self.references_temp.keys())[i]
			func += params[i] * self.fit_refs[i] * self.references_temp[key].y
		func+=params[-2]*x + params[-1]
		return func

	def reference_function_shift(self, x, *params):
		#returns function with reference spectra
		#selected by fit_refs ([1,0,1,1 etc])
		#m*x + b
		func=0
		for i in range(0,len(self.references_temp)):
			key = list(self.references_temp.keys())[i]
			spec=self.references_temp[key]
			#spec=spec << params[-1]
			func += params[i] * self.fit_refs[i]  * (spec << params[-1]).y
		func+=params[-3]*x + params[-2]
		return func

	def fit_references(self, fit_refs="all",fit_funcs=[], fit_range="all", show=False):
		#fits spectrum to reference spectra selected by fit_refs ([1,0,1,1 etc]) + a constant
		#fit_range selects roi. If "all" the whole spectrum is taken
		if fit_refs=="all":
			self.fit_refs=np.ones(len(self.references))
		else:
			self.fit_refs=fit_refs
		if fit_range=="all":
			self.references_temp=self.references
			fit_range=(np.min(self.spectrum.x),np.max(self.spectrum.x))
		else:
			self.references_temp=dict()
			for key in self.references:
				self.references_temp[key]=self.references[key](fit_range[0],fit_range[1])

		self.spectrum_temp= self.spectrum(fit_range[0],fit_range[1])
		fparams=np.ones(len(self.references)+1)/np.sum(self.fit_refs)
		fparams[-1]=0
		print(fparams)
		limits=np.ones((2,len(self.references)+1))
		limits[0,:]=0
		limits[0,-1]=-1
		print(limits)
		popt, pcov = curve_fit(self.reference_function,
							self.spectrum_temp.x,
							self.spectrum_temp.y,
							p0=fparams, bounds=limits)

		popt[:-1]=popt[:-1]*self.fit_refs
		print("----- Fit results -----")
		for i in range(0,len(popt[:-1])):
			print(list(self.references.keys())[i], str(np.round(popt[i]/np.sum(popt[:-1])*100)) + "%")
		print("average difference square in fit area:")
		diff_fit_sq = np.sum((self.spectrum_temp.y
					- self.reference_function(self.spectrum_temp.x,*popt))**2)/len(self.spectrum_temp.y)
		print(diff_fit_sq)
		self.spectrum_temp=self.spectrum
		for key in self.references:
			self.references_temp[key]=self.references[key]
		print("average difference square in total area:")
		diff_tot_sq = np.sum((self.spectrum_temp.y
					- self.reference_function(self.spectrum_temp.x,*popt))**2)/len(self.spectrum_temp.y)
		print(diff_tot_sq)
		if show:
			plt.figure()
			plt.axvspan(fit_range[0],fit_range[1],color="gray", alpha=0.25,
						label="Fit area: " + '{:0.1e}'.format(diff_fit_sq))
			plt.plot(self.spectrum.x,self.spectrum.y,label="Spectrum")
			plt.plot(self.spectrum.x,self.reference_function(self.spectrum.x,*popt),
					label="Fit: " + '{:0.1e}'.format(diff_tot_sq))
			plt.legend()
		print(popt)
		return spectrum(self.spectrum.x,self.reference_function(self.spectrum.x,*popt)), popt, diff_fit_sq, diff_tot_sq

	def fit_references2(self, fit_refs="all",fit_funcs=[], fit_range="all", show=False):
		#fits spectrum to reference spectra selected by fit_refs ([1,0,1,1 etc]) + a constant
		#fit_range selects roi. If "all" the whole spectrum is taken
		if fit_refs=="all":
			self.fit_refs=np.ones(len(self.references))
		else:
			self.fit_refs=fit_refs
		if fit_range=="all":
			self.references_temp=self.references
			fit_range=(np.min(self.spectrum.x),np.max(self.spectrum.x))
		else:
			self.references_temp=dict()
			for key in self.references:
				self.references_temp[key]=self.references[key](fit_range[0],fit_range[1])

		self.spectrum_temp= self.spectrum(fit_range[0],fit_range[1])
		fparams=np.ones(len(self.references)+2)/np.sum(self.fit_refs)
		fparams[-2]=0
		fparams[-1]=0
		print(fparams)
		limits=np.ones((2,len(self.references)+2))
		limits[0,:]=0
		limits[0,-1]=-1
		limits[0,-2]=-1
		print(limits)
		popt, pcov = curve_fit(self.reference_function,
							self.spectrum_temp.x,
							self.spectrum_temp.y,
							p0=fparams, bounds=limits)

		popt[:-2]=popt[:-2]*self.fit_refs
		print("----- Fit results -----")
		for i in range(0,len(popt[:-2])):
			print(list(self.references.keys())[i], str(np.round(popt[i]/np.sum(popt[:-2])*100)) + "%")
		print("average difference square in fit area:")
		diff_fit_sq = np.sum((self.spectrum_temp.y
					- self.reference_function(self.spectrum_temp.x,*popt))**2)/len(self.spectrum_temp.y)
		print(diff_fit_sq)
		self.spectrum_temp=self.spectrum
		for key in self.references:
			self.references_temp[key]=self.references[key]
		print("average difference square in total area:")
		diff_tot_sq = np.sum((self.spectrum_temp.y
					- self.reference_function(self.spectrum_temp.x,*popt))**2)/len(self.spectrum_temp.y)
		print(diff_tot_sq)
		if show:
			plt.figure()
			plt.axvspan(fit_range[0],fit_range[1],color="gray", alpha=0.25,
						label="Fit area: " + '{:0.1e}'.format(diff_fit_sq))
			plt.plot(self.spectrum.x,self.spectrum.y,label="Spectrum")
			plt.plot(self.spectrum.x,self.reference_function(self.spectrum.x,*popt),
					label="Fit: " + '{:0.1e}'.format(diff_tot_sq))
			plt.legend()
		print(popt)
		return spectrum(self.spectrum.x,self.reference_function(self.spectrum.x,*popt)), popt, diff_fit_sq, diff_tot_sq

	def fit_references_shift(self, fit_refs="all",fit_funcs=[], fit_range="all", show=False):
		#fits spectrum to reference spectra selected by fit_refs ([1,0,1,1 etc]) + a constant
		#fit_range selects roi. If "all" the whole spectrum is taken
		if fit_refs=="all":
			self.fit_refs=np.ones(len(self.references))
		else:
			self.fit_refs=fit_refs
		if fit_range=="all":
			self.references_temp=self.references
			fit_range=(np.min(self.spectrum.x),np.max(self.spectrum.x))
		else:
			self.references_temp=dict()
			for key in self.references:
				self.references_temp[key]=self.references[key](fit_range[0],fit_range[1])

		self.spectrum_temp= self.spectrum(fit_range[0],fit_range[1])
		fparams=np.ones(len(self.references)+3)/np.sum(self.fit_refs)
		fparams[-1]=0
		fparams[-2]=0
		fparams[-3]=0
		print(fparams)
		limits=np.ones((2,len(self.references)+3))
		limits[0,:]=0
		limits[0,-1]=-1
		limits[1,-1]=+1
		limits[0,-2]=-1
		limits[0,-3]=-1
		print(limits)
		popt, pcov = curve_fit(self.reference_function_shift,
							self.spectrum_temp.x,
							self.spectrum_temp.y,
							p0=fparams, bounds=limits)

		#popt[:-3]=popt[:-3]*self.fit_refs
		print("----- Fit results -----")
		for i in range(0,len(popt[:-3])):
			print(list(self.references.keys())[i], str(np.round(popt[i]/np.sum(popt[:-3])*100)) + "%")
		print("average difference square in fit area:")
		diff_fit_sq = np.sum((self.spectrum_temp.y
					- self.reference_function_shift(self.spectrum_temp.x,*popt))**2)/len(self.spectrum_temp.y)
		print(diff_fit_sq)
		self.spectrum_temp=self.spectrum
		for key in self.references:
			self.references_temp[key]=self.references[key]
		print("average difference square in total area:")
		diff_tot_sq = np.sum((self.spectrum_temp.y
					- self.reference_function_shift(self.spectrum_temp.x,*popt))**2)/len(self.spectrum_temp.y)
		print(diff_tot_sq)
		if show:
			plt.figure()
			plt.axvspan(fit_range[0],fit_range[1],color="gray", alpha=0.25,
						label="Fit area: " + '{:0.1e}'.format(diff_fit_sq))
			plt.plot(self.spectrum.x,self.spectrum.y,label="Spectrum")
			plt.plot(self.spectrum.x,self.reference_function_shift(self.spectrum.x,*popt),
					label="Fit: " + '{:0.1e}'.format(diff_tot_sq))
			plt.legend()
		print(popt)
		return spectrum(self.spectrum.x,self.reference_function_shift(self.spectrum.x,*popt)), popt, diff_fit_sq, diff_tot_sq

	def add_mask(self, mask):
		#masks an area to avoid braggs etc. mask has the form [ x1, x2, y1, y2]
		self.mask[mask[0]:mask[1],mask[2]:mask[3]]=0
		for image in self.det_images:
			self.det_images[image].data=self.det_images[image].data*self.mask
		self.det_image=self.det_image*self.mask

	'''
	def zoom_rois(self):
		plt.ioff()
		class Index:
			rois=[]
			i=0
			def take_roi(self,event):
				titlestring = "ROI nr." + str(self.i)
				self.rois.append([ax.get_ylim()[0],ax.get_ylim()[1]])
				print(self.rois)
				ax.relim()
				ax.autoscale()
				plt.draw()
				self.i+=1
			def save_rois(self,event):
				plt.close("all")
				return self.rois

		fig, ax = plt.subplots()
		plt.subplots_adjust(bottom=0.2)
		image=self.det_image
		figure_obj = plt.imshow(image)
		callback=Index()
		take_button = Button(plt.axes([0.2,0.05,0.15,0.075]), "set roi")
		take_button.on_clicked(callback.take_roi)
		save_button = Button(plt.axes([0.6,0.05,0.15,0.075]), "save rois")
		save_button.on_clicked(callback.save_rois)
		print(callback.rois)
		#self.crys_number =self.rois.shape
		plt.show()
		return take_button, save_button
		'''



