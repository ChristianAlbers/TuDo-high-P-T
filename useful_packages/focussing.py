import numpy as np
import matplotlib.pyplot as plt
import TuDo_analytix.XES as XES
import glob
import fabio
from spect import *
import time
from scipy.signal import find_peaks

plt.ion()
plt.close("all")

path="D:/data/P01/2022_06/raw/"
roi_number=3



def auto_rois(spec, roi_number=roi_number, show=False):
	spec=spec
	broaden=0.01
	#y_spec=spectrum(range(0,img.shape[0]),y_median)
	y_spec=spectrum(range(0,len(spec)),spec)
	y_peaks, y_properties = find_peaks(y_spec.b(broaden).y)
	while len(y_peaks) > roi_number:
		broaden+=0.01
		y_peaks, y_properties = find_peaks(y_spec.b(broaden).y)
	return y_peaks


temp_fig=plt.figure(1)
ax2 = plt.subplot2grid((2, 1), (0, 0))
ax1 = plt.subplot2grid((2, 1), (1, 0))

while plt.fignum_exists(1):
	folders=sorted(glob.glob(path + "*/*"))
	try:
		folder=folders[-1]
		#folder = "D:/data/P01/2022_06/raw/eh2_11012742_03103/pilatus/"
		img=fabio.open(glob.glob(folder + "/*")[0]).data
		for file in glob.glob(folder + "/*")[1:]:
			img+=fabio.open(file).data
	except IndexError:
		print("No folder found")
	except PermissionError:
		print("Permission denied, try again")
	else:
		ax1.clear()
		ax2.clear()
		ax1 = plt.subplot2grid((2, 1), (0, 0))
		ax2 = plt.subplot2grid((2, 1), (1, 0))
		ax1.imshow(img,cmap="jet")
		spec=np.sum(img[5:-5,50:-50],axis=1)
		spect=spectrum(range(0,len(spec)),spec)
		ax2.plot(spect.y)
		rois=auto_rois(spect.y)
		for roi in rois:
			temp_spec=spect(roi-10,roi+10)
			diff=temp_spec.y[0]
			temp_spec-=diff
			try:
				gauss=temp_spec.gauss(fullinfo=1)
				(gauss+diff).p(newlabel=str(gauss.fit_p[0].astype(int)) + ", " + str(np.round(gauss.fit_p[2],2)))
				plt.legend()
			except:
				None
			
	plt.pause(1)
plt.close("all")
