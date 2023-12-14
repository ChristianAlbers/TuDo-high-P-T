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
from scipy.signal import find_peaks

class XES_elastics:
	def __init__(self,
				images = [], #list of paths to pictures. Sorted!!!
				energies = [], #e.g. [10500,10550,10600,10650,10700]
				reflex_fac = 1.5, #for 660 to 440
				crys_number=4,
				):
		self.mask=[0]
		self.det_image=0
		self.images= dict()
		if any(images):
			self.load_images(images)

		self.energies=energies
		self.reflex_fac=reflex_fac
		self.params=np.zeros((crys_number,2)) # [[m], [b]] for y=m*x+b
		self.params[:,0]=1
		self.crys_number=crys_number
		self.rois=0

	def load_images(self,
		images):
		if not any(self.images):
			self.mask=np.zeros(fabio.open(images[0]).data.shape)
			self.mask[:]=1
		for file in images:
			self.images[file]=fabio.open(os.path.normpath(file)).data
		self.image = np.zeros(self.images[list(self.images.keys())[0]].data.shape)
		for image in self.images:
			self.image+=self.images[image]



	def calc_params(self,
		energies=False,
		reflex_fac= False,
		mode="auto",
		roiradius=5,
		show=False,
		custom_rois=None
		):
		xpos=np.linspace(0,self.mask.shape[1]-1,self.mask.shape[1])
		el_params=np.zeros((self.crys_number,2))
		box_rois=np.zeros((len(self.images),self.crys_number,4), dtype=int)
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
			for i, image in enumerate(self.images):
				data=self.images[image]
				data=data*self.mask
				if show:
					if i==0:
						show_img=data
					else:
						show_img+=data
				if mode == "auto":
					rois=self.auto_rois(img=self.images[image],roiradius=roiradius,roi_number=self.crys_number)
					box_rois[i]=rois
					#take automatic rois from roifinder
				elif mode == "lines":
					#search for rois in given lines [[line0_y0,line0_y1],[line1_y0,line1_y1],...]
					#box_rois[i]=rois
					for j, line in enumerate(custom_rois):
						y_spec=spectrum(range(0,data.shape[1]),np.sum(data[line-roiradius:line+roiradius],axis=0))
						y_spec_find=y_spec-y_spec.b(0.1)
						h=1
						y_peaks, y_properties = find_peaks(y_spec_find.y,height=h)
						while len(y_peaks) > 1:
							h+=1
							y_peaks, y_properties = find_peaks(y_spec_find.y,height=h)
						y_peak=int(y_peaks[0])
						box_rois[i][j]=[int(line-roiradius),int(line+roiradius),int(y_peak-roiradius),int(y_peak+roiradius)]
						rois=box_rois[i]
				else:
					if type(custom_rois) in [list, np.array, np.ndarray]:
						box_rois[i]=custom_rois[i]
						rois=custom_rois[i]
					else:
						print("choose mode auto or give lines as [[line0_y0,line0_y1],[line1_y0,line1_y1],...] or  give custom_rois as [[y0,y1,x0,x1],[y0,y1,x0,x1],...]")
						return None
				roiscount=0
				for roi in rois:
					crystal=spectrum(xpos,np.array(np.sum(data[roi[0]:roi[1]],axis = 0)))
					peaks[roiscount,imagecount]=crystal(roi[2],roi[3]).gauss(fullinfo=1).fit_p[1]
					if show:
						plt.figure(roiscount)
						crystal.p()
						crystal(roi[2],roi[3]).gauss().p()
						plt.figure(999)
						plt.plot([roi[2],roi[3]],[roi[0],roi[0]],color="red")
						plt.plot([roi[2],roi[3]],[roi[1],roi[1]],color="red")
						plt.plot([roi[2],roi[2]],[roi[0],roi[1]],color="red")
						plt.plot([roi[3],roi[3]],[roi[0],roi[1]],color="red")

						#crystal(roi[1]-10,roi[1]+10).gauss(fullinfo=1).p()
					roiscount+=1
				imagecount+=1
				if show:
					plt.figure(999)
					plt.imshow(show_img,vmax=np.percentile(show_img,q=99.9))
			for r in range(0,len(peaks)):
				popt=curve_fit(lin,peaks[r],ener)
				el_params[r]=popt[0]/fac
			if show:
				fig = plt.figure(99)
				for r in range(0,len(rois)):
					plt.plot(peaks[r],ener/fac,marker="o",linestyle="")
				fig.gca().set_prop_cycle(None)
			if show:
				plt.figure(99)
				for r in range(0,len(rois)):
					plt.plot(xpos,xpos*el_params[r][0]+el_params[r][1])
			self.rois=box_rois
			self.params=el_params
	'''
	old
	def peakfinder(self,
					image,
					num_rois=4,
					distance = 20): #finds elastics automatically by local
								#maxima
		im=image.data
		peaks = peak_local_max(im, min_distance=distance, num_peaks=num_rois)
		peaks = sorted(peaks ,key=lambda x: x[0])
		return peaks
	'''
	def auto_rois(self, img, roi_number,roiradius, show=False):
		out_rois=np.zeros((roi_number,4), dtype=int)	#out_rois[0][:]=[y0,y1,x0,x1] etc.
		y_median=np.percentile(img,q=98,axis=1)
		if show:
			plt.figure(999)
			plt.imshow(img,vmax=np.percentile(img,q=99.9))
		broaden=0.01
		y_spec=spectrum(range(0,img.shape[0]),y_median)
		y_peaks, y_properties = find_peaks(y_spec.b(broaden).y)
		while len(y_peaks) > roi_number:
			broaden+=0.01
			y_peaks, y_properties = find_peaks(y_spec.b(broaden).y)

		if show:
			plt.figure()
			y_spec.p()
			print(y_peaks)
		for i, peak in enumerate(y_peaks):

			#y_popt=y_spec(peak-5,peak+5).gauss(fullinfo=1).fit_p
			out_rois[i,0]=int(peak-roiradius) 	#y0
			out_rois[i,1]=int(peak+roiradius)	#y1

			line_img=img[int(peak-roiradius):int(peak+roiradius),:]
			x_median=np.percentile(line_img,q=90,axis=0)
			x_spec=spectrum(range(0,img.shape[1]),x_median)
			x_spec_find=x_spec-x_spec.b(0.1)
			if show:
				x_spec_find.p()
			h=1
			x_peaks, x_properties = find_peaks(x_spec_find.y,height=h)
			while len(x_peaks) > 1:
				h+=1
				x_peaks, x_properties = find_peaks(x_spec_find.y,height=h)

			if show:
				plt.plot([x_peaks[0],x_peaks[0]],[0,h],color="k")

			out_rois[i,2]=int(x_peaks[0]-roiradius) 	#x0
			out_rois[i,3]=int(x_peaks[0]+roiradius)	#x1

		if show:
			plt.figure(999)
			#plt.imshow(img,vmax=np.percentile(img,q=99.9))
			for roi in out_rois:
				plt.plot([roi[2],roi[3]],[roi[0],roi[0]],color="red")
				plt.plot([roi[2],roi[3]],[roi[1],roi[1]],color="red")
				plt.plot([roi[2],roi[2]],[roi[0],roi[1]],color="red")
				plt.plot([roi[3],roi[3]],[roi[0],roi[1]],color="red")
		return out_rois

	def add_mask(self, mask):
		#masks an area to avoid braggs etc. mask has the form [ x1, x2, y1, y2]
		self.mask[mask[0]:mask[1],mask[2]:mask[3]]=0

	def show_images(self):
		for im in self.images:
			plt.figure()
			rois=self.peakfinder(image=self.images[im],num_rois=self.crys_number)
			for roi in rois:
				plt.axhspan(roi[0]-5,roi[0]+5,color="red",alpha=0.25)
			plt.imshow(self.images[im],vmax=self.images[im].max()/5)