import numpy as np
import matplotlib.pyplot as plt
import TuDo_analytix_dev.XES as XES
import glob
import fabio
from spect import *
import time


plt.ion()
plt.close("all")

#####config#####
path="D:/data/P01/2022_06/raw/"
prefix="eh2_11012742_"
detector="pilatus"
snr_offset=5000

#####elastic
elasticnumber=3108
#elasticnumber=int(input("elasticnumber: "))
el_path=sorted(glob.glob(path + prefix + "%05d"%int(elasticnumber) + "/"+ detector + "/" + "*.tif"))
elastic=XES.XES_elastics(el_path, np.linspace(10500,10680,7),crys_number=3)
elastic.calc_params(mode="auto",roiradius=4,show=False)
print(elastic.params)


plt.close("all")

temp_fig=plt.figure(1)
temp_fig.set_size_inches(12,6)
ax1 = plt.subplot2grid((4, 2), (0, 0))
ax2 = plt.subplot2grid((4, 2), (1, 0), rowspan=3)
ax3 = plt.subplot2grid((4, 2), (0, 1), rowspan=4, colspan=2)
left, bottom, width, height = [0.5, 0.25, 0.45, 0.4]

lastfolder = ""
oldfiles = 0
while plt.fignum_exists(1):
	spect_path=[]
	folders=sorted(glob.glob(path + prefix + "*/"))
	folder = folders[-1]
	if folder != lastfolder:
		oldfiles = 0
		spect_path=[]
		print(folder)
		lastfolder = folder
	for file in list(glob.glob(folder + "/" + detector + "/" + "*.tif")):
		spect_path.extend([file])
	if len(spect_path) > oldfiles:
		oldfiles = len(spect_path)
		print(len(spect_path))
		pass
	else:
		try:
			time.sleep(1)
		except:
			break
		plt.pause(1)
		continue
#	if len(spect_path) == 0:
#		time.sleep(1)
#		continue
	spect=XES.XES_object(spect_path, elastic,crys_number=3)
	spect.auto_rois(roiwidth=3)
	spect.calc_spectrum_all(bgrwidth=5,bgrfit=0,bgrsmooth=0.105,bgr_skip=5,bgr_range=(120,130),show=False)
	ax1.clear()
	ax2.clear()
	ax3.clear()
	ax1 = plt.subplot2grid((3, 2), (0, 0))
	ax2 = plt.subplot2grid((3, 2), (1, 0), rowspan=2)
	ax3 = plt.subplot2grid((3, 2), (0, 1), rowspan=3, colspan=2)

	ax1.imshow(spect.det_image)

	ax2_inset = ax2.inset_axes([left, bottom, width, height])
	ax2_inset.set_xlim([7075,7115])
	labelcolors=[]
	for r, roi in enumerate(spect.rois):
		ax1.axhspan(roi[0],roi[1],color="C"+str(r),alpha=0.2)
		crystal=spectrum(range(0,487)*spect.params[r][0]+spect.params[r][1],np.sum(spect.det_image[roi[0]:roi[1],:],axis=0))
		background=spectrum(range(0,487)*spect.params[r][0]+spect.params[r][1],np.sum(spect.det_image[roi[1]+1:roi[1]+1+(roi[1]-roi[0]),:],axis=0))
		max_pos=np.argmax(crystal.y)
		snr=int((crystal.y[max_pos]-background.y[max_pos])/np.sqrt(background.y[max_pos]))
		ax2.plot(crystal(7020,7120).x,crystal(7020,7120).y,label="SNR: " +str(snr))
		ax2_inset.plot(crystal(7070,7120).x,crystal(7070,7120).y)
		ax2_inset.set_yticks([])
		if snr <snr_offset:
			labelcolors.append("r")
		else:
			labelcolors.append("g")
	legend=ax2.legend(labelcolor=labelcolors)
	print(labelcolors)
	ax3.plot(spect.spectrum.x,spect.spectrum.y)
	ax3.set_xlim([7020,7120])
	ax3_inset = ax3.inset_axes([left, bottom, width, height])
	ax3_inset.set_xlim([7075,7115])
	ax3_inset.plot(spect.spectrum(7070,7120).x,spect.spectrum(7070,7120).y)
	ax3_inset.set_yticks([])
	plt.suptitle(folder)
	plt.pause(1)
plt.ioff()
plt.close("all")
