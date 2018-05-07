#import sys
from spect import *
import matplotlib.pyplot as plt
import numpy as np
import glob
from PIL import Image
import collections

data = dict()
pilatusimages = dict()
sumxtal1 = dict()
sumxtal2 = dict()
sumxtal3 = dict()
sumxtal4 = dict()
bgrxtal1ll = dict()
bgrxtal2ll = dict()
bgrxtal3ll = dict()
bgrxtal4ll = dict()
bgrxtal1ul = dict()
bgrxtal2ul = dict()
bgrxtal3ul = dict()
bgrxtal4ul = dict()
bgrxtal1 = dict()
bgrxtal2 = dict()
bgrxtal3 = dict()
bgrxtal4 = dict()
sumxtal1bgr = dict()
sumxtal2bgr = dict()
sumxtal3bgr = dict()
sumxtal4bgr = dict()

prefix = "image_000"	#adapt
path2spectra = "C:\2018-04-12-XES-DELTA\run06_18_xes_vh" #adapt

#set rois
#fe-foil rois
#rois = [15, 38, 43, 67, 86, 121, 146, 174] # adapt: lL xtal1, uL xtal1, lL xtal2, ... later via sys.argv[...]
#everything else
rois = [22, 36, 48, 66, 93, 110, 153, 171]
#set bgr rois height
fact = 3
bgr_roi_ll = [rois[0]-fact,rois[1]-fact,rois[2]-fact,rois[3]-fact,rois[4]-fact,rois[5]-fact,rois[6]-fact,rois[7]-fact]
bgr_roi_ul = [rois[0]+fact,rois[1]+fact,rois[2]+fact,rois[3]+fact,rois[4]+fact,rois[5]+fact,rois[6]+fact,rois[7]+fact]


for file in glob.glob(prefix + "*.tif")[:int(62):int(1)]:	#step size later sys.argv[1] or range from a to b [a:b]
	T = file.split("_")[1]
	pilatusimages[T] = Image.open(file)
	data[T] = np.array(pilatusimages[T])	
	data[T][46][18] = 3	#correct dead pixel to bgr value, only BL9 DELTA Pilatus 100K s specific

	# plt.imshow(data[T], interpolation = 'nearest')
	# print T
	# raw_input("Press Enter to continue...")		# show picture for setting the rois
	# plt.show()										# not needed later on
	# plt.close()
	
	#xtal summed signals normalized to pixel height, pixel length is equal for all
	sumxtal1[T] = np.sum(data[T][:][rois[6]:rois[7]], axis = 0)
	sumxtal1[T] = np.array(sumxtal1[T])/(len(data[T][:][rois[6]:rois[7]]))
	sumxtal2[T] = np.sum(data[T][:][rois[4]:rois[5]], axis = 0)
	sumxtal2[T] = np.array(sumxtal2[T])/(len(data[T][:][rois[4]:rois[5]]))
	sumxtal3[T] = np.sum(data[T][:][rois[2]:rois[3]], axis = 0)
	sumxtal3[T] = np.array(sumxtal3[T])/(len(data[T][:][rois[2]:rois[3]]))
	sumxtal4[T] = np.sum(data[T][:][rois[0]:rois[1]], axis = 0)
	sumxtal4[T] = np.array(sumxtal4[T])/(len(data[T][:][rois[0]:rois[1]]))
	#xtal lower limit (ll) bgrs
	bgrxtal1ll[T] = np.sum(data[T][:][bgr_roi_ll[6]:rois[6]], axis = 0)
	bgrxtal2ll[T] = np.sum(data[T][:][bgr_roi_ll[4]:rois[4]], axis = 0)
	bgrxtal3ll[T] = np.sum(data[T][:][bgr_roi_ll[2]:rois[2]], axis = 0)
	bgrxtal4ll[T] = np.sum(data[T][:][bgr_roi_ll[0]:rois[0]], axis = 0)
	#xtal upper limit (ul) bgrs
	bgrxtal1ul[T] = np.sum(data[T][:][rois[7]:bgr_roi_ul[7]], axis = 0)
	bgrxtal2ul[T] = np.sum(data[T][:][rois[5]:bgr_roi_ul[5]], axis = 0)
	bgrxtal3ul[T] = np.sum(data[T][:][rois[3]:bgr_roi_ul[3]], axis = 0)
	bgrxtal4ul[T] = np.sum(data[T][:][rois[1]:bgr_roi_ul[1]], axis = 0)
	#xtal summed bgrs normalized to pixel height
	bgrxtal1[T] = (bgrxtal1ll[T]/(len(data[T][:][bgr_roi_ll[6]:rois[6]])) + bgrxtal1ul[T]/len(data[T][:][rois[7]:bgr_roi_ul[7]]))/2
	bgrxtal2[T] = (bgrxtal2ll[T]/(len(data[T][:][bgr_roi_ll[4]:rois[4]])) + bgrxtal2ul[T]/len(data[T][:][rois[5]:bgr_roi_ul[5]]))/2 
	bgrxtal3[T] = (bgrxtal3ll[T]/(len(data[T][:][bgr_roi_ll[2]:rois[2]])) + bgrxtal3ul[T]/len(data[T][:][rois[3]:bgr_roi_ul[3]]))/2 
	bgrxtal4[T] = (bgrxtal4ll[T]/(len(data[T][:][bgr_roi_ll[6]:rois[6]])) + bgrxtal4ul[T]/len(data[T][:][rois[1]:bgr_roi_ul[1]]))/2 
	#subst bgrs from signals
	sumxtal1bgr[T] = sumxtal1[T]-bgrxtal1[T]
	sumxtal2bgr[T] = sumxtal2[T]-bgrxtal2[T]
	sumxtal3bgr[T] = sumxtal3[T]-bgrxtal3[T]
	sumxtal4bgr[T] = sumxtal4[T]-bgrxtal4[T]

	#detector pixel length, energy x-axis
	xpos = np.linspace(1, 487, 487)

	#from here Georgs spect script
	specxtal1 = spectrum(xpos, sumxtal1bgr[T])
	specxtal2 = spectrum(xpos, sumxtal2bgr[T])
	specxtal3 = spectrum(xpos, sumxtal3bgr[T])
	specxtal4 = spectrum(xpos, sumxtal4bgr[T])
	
		
	# save all xtal spectra
	# specxtal1.save(file + ".XTAL1", force = True)
	# specxtal2.save(file + ".XTAL2", force = True)
	# specxtal3.save(file + ".XTAL3", force = True)
	# specxtal4.save(file + ".XTAL4", force = True)

	# sum up all four xtals
	# a = specxtal1+specxtal2+specxtal3+specxtal4
	# a.save(file +".ALLXTALS", force=True) # save shifted + summed spectra for further analysis (energy-scale, cutting, bgr, normalisation, fitting...)
	# b = a.xenerg()	# scale x-axis from channel to eV with params from regression
	# b.save(file +".ALLXTALS.ESCALE", force = True)
