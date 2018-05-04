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

prefix = "image_000"	#adapt
path2spectra = "C:\2018-04-12-XES-DELTA\run06_18_xes_vh" #adapt

#set rois
rois = [18, 38, 43, 67, 86, 121, 146, 174] # adapt: lL xtal1, uL xtal1, lL xtal2, ... later via sys.argv[...]
#set bgr rois size
fact = 3
bgr_roi_ll = [rois[0]-fact,rois[1]-fact,rois[2]-fact,rois[3]-fact,rois[4]-fact,rois[5]-fact,rois[6]-fact,rois[7]-fact]
bgr_roi_ul = [rois[0]+fact,rois[1]+fact,rois[2]+fact,rois[3]+fact,rois[4]+fact,rois[5]+fact,rois[6]+fact,rois[7]+fact]
#insert the xtalshift from line 49-57
xtalshift = [45, 21, 28]

# Elastic Scans
path2spectra = "C:\2018-04-12-XES-DELTA\run06_18_xes_vh"
samplelist = [\
"image_00045.tif.ALLXTALS",\
"image_00048.tif.ALLXTALS",\
"image_00049.tif.ALLXTALS",\
"image_00050.tif.ALLXTALS"]
##########

for file in glob.glob(prefix + "*.tif")[::int(1)]:	#step size(sys.argv[1])::int(sys.argv[1])
	T = file.split("_")[1]
	pilatusimages[T] = Image.open(file)
	data[T] = np.array(pilatusimages[T])	
	data[T][46][18] = 3	#correct dead pixel to bgr value, only BL9 DELTA Pilatus 100K s specific

	# plt.imshow(data[T], interpolation = 'nearest')	
	# plt.show()
	# plt.close()
	
	#xtal summed signals
	sumxtal1[T] = np.sum(data[T][:][rois[6]:rois[7]], axis = 0)
	sumxtal1[T] = np.array(sumxtal1[T])/(len(data[T][:][rois[6]:rois[7]]))
	sumxtal2[T] = np.sum(data[T][:][rois[4]:rois[5]], axis = 0)
	sumxtal2[T] = np.array(sumxtal2[T])/(len(data[T][:][rois[4]:rois[5]]))
	sumxtal3[T] = np.sum(data[T][:][rois[2]:rois[3]], axis = 0)
	sumxtal3[T] = np.array(sumxtal3[T])/(len(data[T][:][rois[2]:rois[3]]))
	sumxtal4[T] = np.sum(data[T][:][rois[0]:rois[1]], axis = 0)
	sumxtal4[T] = np.array(sumxtal4[T])/(len(data[T][:][rois[0]:rois[1]]))
	#xtal ll bgrs
	bgrxtal1ll[T] = np.sum(data[T][:][bgr_roi_ll[6]:rois[6]], axis = 0)
	bgrxtal2ll[T] = np.sum(data[T][:][bgr_roi_ll[4]:rois[4]], axis = 0)
	bgrxtal3ll[T] = np.sum(data[T][:][bgr_roi_ll[2]:rois[2]], axis = 0)
	bgrxtal4ll[T] = np.sum(data[T][:][bgr_roi_ll[0]:rois[0]], axis = 0)
	#xtal ul bgrs
	bgrxtal1ul[T] = np.sum(data[T][:][rois[7]:bgr_roi_ul[7]], axis = 0)
	bgrxtal2ul[T] = np.sum(data[T][:][rois[5]:bgr_roi_ul[5]], axis = 0)
	bgrxtal3ul[T] = np.sum(data[T][:][rois[3]:bgr_roi_ul[3]], axis = 0)
	bgrxtal4ul[T] = np.sum(data[T][:][rois[1]:bgr_roi_ul[1]], axis = 0)
	#xtal summed bgrs
	bgrxtal1[T] = (bgrxtal1ll[T]/(len(data[T][:][bgr_roi_ll[6]:rois[6]])) + bgrxtal1ul[T]/len(data[T][:][rois[7]:bgr_roi_ul[7]]))/2
	bgrxtal2[T] = (bgrxtal2ll[T]/(len(data[T][:][bgr_roi_ll[4]:rois[4]])) + bgrxtal2ul[T]/len(data[T][:][rois[5]:bgr_roi_ul[5]]))/2 
	bgrxtal3[T] = (bgrxtal3ll[T]/(len(data[T][:][bgr_roi_ll[2]:rois[2]])) + bgrxtal3ul[T]/len(data[T][:][rois[3]:bgr_roi_ul[3]]))/2 
	bgrxtal4[T] = (bgrxtal4ll[T]/(len(data[T][:][bgr_roi_ll[6]:rois[6]])) + bgrxtal4ul[T]/len(data[T][:][rois[1]:bgr_roi_ul[1]]))/2 


	xpos = np.linspace(1, 487, 487)
	#plt.plot(xpos + xtalshift[0], sumxtal1[T], label = "xtal 1 scann: " +T)
	#plt.plot(xpos + xtalshift[1], sumxtal2[T], label = "xtal 2 scann: " +T)
	#plt.plot(xpos + xtalshift[2], sumxtal3[T], label = "xtal 3 scann: " +T)
	#plt.plot(xpos, sumxtal4[T], label = "xtal 4 scann: " +T)
	#plt.show()
	

	# find maximum of xtal and calc diff to xtal 4 (once done with Fe-foil Scann 62)
	# xtalmax1= np.argmax(sumxtal1[T])
	# print(xtalmax1)
	# xtalmax2= np.argmax(sumxtal2[T])
	# print(xtalmax2)
	# xtalmax3= np.argmax(sumxtal3[T])
	# print(xtalmax3)
	# xtalmax4= np.argmax(sumxtal4[T])
	# print(xtalmax4)

	#from here Georgs sepct script
	specxtal1 = spectrum(xpos, sumxtal1[T])
	specxtal2 = spectrum(xpos, sumxtal2[T])
	specxtal3 = spectrum(xpos, sumxtal3[T])
	specxtal4 = spectrum(xpos, sumxtal4[T])
	
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
