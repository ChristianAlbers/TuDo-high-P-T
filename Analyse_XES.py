from IPython import get_ipython
get_ipython().magic('reset -sf')
from spect import *
import matplotlib.pyplot as plt
import numpy as np
import glob
from PIL import Image
import collections
import elastic
from scipy.interpolate import interp1d
plt.close('all')

plt. ion()
Scan = dict()
spectsges=dict()



#--------------------------------------------------------------------------------------------------------------------------------------------
#Config starts here


#set the path 
path2spectra = r"data/alignment_" #adapt
#set rois
rois = [(60, 66), (94, 100), (123, 129), (148, 154)]
#set bgr rois
bgrwidth = 3
bgroffset = 10
#Set elesticscan
elasticnumber = 674
energystart = 10560
energystep =15
#set detectorlength
xpos = np.linspace(1, 487, 487)
#set your energyrange
erange=(7000,7150)
#select only some crystalls? Try show_crystalls(number) first to decide
selectcrystalls=[1, 2]
#select your Scannumbers [[spectra1],[spectra2]]
scannumbers = [[458,470],[676,677,678]]
#Wanna save? Type 1
save=1



#Config ends here
#--------------------------------------------------------------------------------------------------------------------------------------------


#fit your elastic und save Parameter
elasticfit = dict()
params=elastic.load(path2spectra, elasticnumber, energystart, energystep, rois)
for x in range(0,len(scannumbers)):
	spects=dict()
	for number in scannumbers[x]:
		data = np.zeros((195, 487))
		crystall = np.zeros((len(rois),) + (487,))
		bckround = np.zeros(np.shape(crystall))
		crystall_WBCK  = np.zeros(np.shape(crystall))
		spect=spectrum()
		T = "Scan" +("%05d" % number)
		for file in glob.glob(path2spectra +  ("%05d" % number) + "/pilatus_100k/alignment_" + "*.tif")[::int(1)]: 
			pilatusimages = Image.open(file)
			data = np.array(pilatusimages)	
			#data[46][18] = 3	#correct dead pixel to bgr value, only BL9 DELTA Pilatus 100K s specific
			for (x1,x2),roiscount in zip(rois,range(0,len(rois))):
				crystall[roiscount] += np.array(np.sum(data[x1:x2], axis = 0))/(len(data[x1:x2]))
				bckround[roiscount] += 0.5*(np.sum(data[x2+bgroffset:x2+bgroffset+bgrwidth], axis = 0)/len(data[x2+bgroffset:x2+bgroffset+bgrwidth])+np.sum(data[x1-bgroffset-bgrwidth:x1-bgroffset], axis = 0)/len(data[x1-bgroffset-bgrwidth:x1-bgroffset]))
				crystall_WBCK[roiscount] += (crystall[roiscount] - bckround[roiscount])
				#plt.plot(crystall_WBCK[roiscount])			
		Scan[T] = crystall_WBCK
		try:
			selectcrystalls
		except:
			for roiscount in (range(0,len(rois))):
				spect=spect.__and__(spectrum((np.flipud(xpos)*params[roiscount][0]+params[roiscount][1])/1.5,np.flipud(Scan[T][roiscount,:])).n())
		else:
			for roiscount in (selectcrystalls[:]):
				spect=spect.__and__(spectrum((np.flipud(xpos)*params[roiscount][0]+params[roiscount][1])/1.5,np.flipud(Scan[T][roiscount,:])).n())
		spect=spect.__call__(erange[0],erange[1]).n()
		spects[T]=spect
	#plot all spectra in one plot
	plt.figure(x+1)
	spectsges[x]=spectrum()
	for scan in spects:
		crystallinfo=""
		try:
			selectcrystalls
		except:
			crystallinfo+=str(np.arange(0,len(rois)))
		else:
			crystallinfo+=str(selectcrystalls[:])
		spectsges[x]=spectsges[x].__and__(spects[scan])
		spectsges[x]=spectsges[x].n()
		spects[scan].p(newlabel=scan)
	#Save all spectra to seperate files
		if save==1:
			spects[scan].save(scan + ".spect",force=1)
			infofile = open(scan + ".info","w")
			infofile.write("crystalls: " + str(crystallinfo) + "\n"
			"energy range: " + str(erange))
			infofile.close()
	spectsges[x].p(newlabel="gesamt")
	plt.legend()


#compare several spectra ([spectrum1,spectrum2],convolute=?)
def compare_spectra(numbers,convolute=0.00001):
	plt.close(9)
	plt.figure(9)
	for scan in (numbers):
		spectsges[scan].b(convolute).n().p(newlabel=scannumbers[scan])
	plt.legend()


def show_crystalls(number):
#shows all crystalls of one Scan with the average Pixel-FWHM of elastic Scans
	plt.close(99)
	T = "Scan" +("%05d" % number)
	plt.figure(99)
	plt.title(number)
	for roiscount in (range(0,len(rois))):
		spectrum((np.flipud(xpos)*params[roiscount][0]+params[roiscount][1])/1.5,np.flipud(Scan[T][roiscount,:])).n().p(newlabel="Crystall: " + str(roiscount) + "; FWHM=" + str(format(params[roiscount][2], ".2f")))
	plt.legend()
#give numbers in [number1, number2]
def show_spectra(numbers):
	plt.close(999)
	plt.figure(999)
	for number in (numbers):
		T = "Scan" +("%05d" % number)
		spects[T].p(newlabel=T)
	plt.legend()

