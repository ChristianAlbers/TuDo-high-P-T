import numpy as np
import matplotlib.pyplot as plt
import TuDo_analytix.XES as XES
from skimage.feature import peak_local_max
import glob
import fabio
from spect import *
import time


def auto_rois(img, distance = 10, roiwidth = 2):
	#set automagic rois to the whole detector image
	peaks = peak_local_max(img, min_distance=distance,
							num_peaks=4)
	peaks = sorted(peaks ,key=lambda x: x[0])
	line_rois=np.zeros((4,1)).astype(int)
	rois=np.zeros((4,2)).astype(int)
	for i in range(0,len(peaks)):
		line_rois[i]=int(peaks[i][0])
	for l, roi in enumerate(line_rois):
		line=np.sum(img[roi[0]-1:roi[0]+1+1,:],axis=0)
		rois[l]=[roi[0],np.argmax(line)]
		#plt.plot(line)
	return rois

plt.ion()
plt.close("all")

path="D:/data/P01/2022_06/raw/"
prefix="eh2_11012742_"
detector="pilatus"
runs=[list(range(2851,2855)),list(range(2858,2862))]              #,list(range(2868,2872)),list(range(2875,2879)),list(range(2885,2889)),list(range(2892,2896))]



old_file_list=[]
temp_fig=plt.figure(1)
temp_fig.set_size_inches(6,6)
#ax1=temp_fig.add_subplot(212)
#ax2=temp_fig.add_subplot(411)
ax1 = plt.subplot2grid((2, 1), (0, 0))
left, bottom, width, height = [0.5, 0.575, 0.2, 0.15]
ax1_inset = temp_fig.add_axes([left, bottom, width, height])
ax2 = plt.subplot2grid((2, 1), (1, 0))
left, bottom, width, height = [0.5, 0.15, 0.2, 0.15]
ax2_inset = temp_fig.add_axes([left, bottom, width, height])
lastfolder = ""
oldfiles = 0
while plt.fignum_exists(1):
	fio_list=[]
	file_list=[]
	for scans in runs:
		for scan in scans:
			for file in sorted(glob.glob(path + prefix + "%05d"%int(scan) + "/" +  detector + "/*.tif")):
				file_list.append(file)
			#try:
			if os.path.isfile(path + prefix + "%05d"%int(scan) + ".fio"):
				fio_list.append(path + prefix + "%05d"%int(scan) + ".fio")
			#except:
			#	None
	if len(old_file_list) < len(file_list):
		old_file_list=file_list
		time.sleep(1)
		pass
	else:
		time.sleep(1)
	try:
		roi_img=fabio.open(file_list[100]).data
	except:
		roi_img=fabio.open(file_list[-1]).data
	rois=auto_rois(roi_img[:,:])
	#rois=[[119,310]]
	print(rois)
	energy=[]
	izero=[]
	herfd=spectrum()
	ax2.clear()
	ax2_inset.clear()
	for r, scans in enumerate(runs):
		try:
			herfd_temp=[]
			energy=[]
			izero=[]
			fios=fio_list[4*r:4*(r+1)]
			for fio in fios:
				try:
					e=np.genfromtxt(fio,skip_header=123)
					if np.isnan(e[-1,0]):
						e=e[:-1]
				except:
					e=np.genfromtxt(fio,skip_header=123,skip_footer=1)
				try:
					energy.extend(np.round(e[:,3],1))
					izero.extend(e[:,5])
				except:
					time.sleep(5)
					raise ValueError
					#continue
		except ValueError:
			continue
		j=0
		file_check=False
		for i,scan in enumerate(scans):
			for file in sorted(glob.glob(path + prefix + "%05d"%int(scan) + "/" + detector + "/" + "*.tif")):
				file_check=True
				img=fabio.open(file).data
				int_alpha=0
				for roi in rois:
					int_alpha+=	np.sum(img[roi[0]-2:roi[0]+2+1:,roi[1]-2:roi[1]+2+1])
				try:
					herfd_temp.append(int_alpha/izero[j])
				except:
					break
				j+=1
				#herfd_temp.append(int_alpha)
		if file_check:
			herfd_temp=spectrum(energy[:j],herfd_temp)
			ax2.plot(herfd_temp.x,herfd_temp.y)
			try:
				ax2_inset.plot(herfd_temp(7106,7116).x,herfd_temp(7106,7116).y/herfd_temp(7106,7116).y[-1])
			except:
				None
			ax2_inset.set_xlim([7106,7116])
			if r ==0:
				herfd=herfd_temp
			else:
				herfd=(herfd+herfd_temp)&(herfd(herfd_temp.x[-1]+1e-5,7300))
	ax1.clear()
	ax1_inset.clear()
	ax1.plot(herfd.x,herfd.y)
	ax1.set_xlim([7090,7300])
	#ax1_inset = temp_fig.add_axes([left, bottom, width, height])
	ax1_inset.plot(herfd(7106,7116).x,herfd(7106,7116).y)
	ax1_inset.set_xlim([7106,7116])
	ax2.set_xlim([7090,7300])
	plt.pause(1)
	#print("I'm still alive, running Scann No: " + str(scan))
plt.ioff()
plt.close("all")
