import numpy as np
import matplotlib.pyplot as plt

from bokeh.io import curdoc, show
from bokeh.layouts import column, row
#from bokeh.models import ColumnDataSource, Slider, TextInput, Button, HoverTool, CheckboxGroup,RangeSlider,LassoSelectTool
from bokeh.models import *
from bokeh.plotting import figure
from bokeh.models.widgets import FileInput, MultiSelect, RadioButtonGroup
from bokeh.plotting import figure, output_file, show
from bokeh.palettes import Dark2, Category10
import TuDo_analytix_dev.XES as XES
from spect import *
from functions import *

import glob
import fabio
import os
from pandas import DataFrame as df

import itertools  
import time
colors = Category10[10]+Dark2[8]


det_img=np.zeros((195,487),dtype=int)
elastic_images=dict()
#colors=["blue","red","green","orange","black","orange"]

#############get info from config#############
image_type="*.tif*"
sample_prefix="C:/Beamtimes/P01_Juli19/raw/eh2_11006060_%05d/pilatus_100k/"
elastic_prefix="C:/Beamtimes/P01_Juli19/raw/eh2_11006060_%05d/pilatus_100k/"%int(1281)
roi_num = 4
det_size=(487,195)
#############create global variables#############
elastic=XES.XES_elastics()
spectra=dict()
exported_spectra=dict()

#############functions for elastics#############
def add_elastic():
	global elastic
	energies_read=str(elastic_energies.value).split(",")
	energies=[]
	for energy in energies_read:
		try:
			energies.append(int(energy))
		except:
			None
	elastic_files=sorted(glob.glob(str(elastic_path.value) %int(elastic_number.value) + image_type))
	if len(energies)==0:
		print("no energies, type energies like E1, E2, E3")
		return None
	elif len(energies) != len(elastic_files):
		print("not the same length of energies and images!")
		print(len(energies))
		print(len(elastic_files))
		return None

	elastic.load_images(elastic_files)
	elastic.energies=energies
	for i, file in enumerate(elastic_files):
		#elastic_list.value.append(file)
		elastic_list.options.append((file,str(energies[i])))
	#elastic=XES.XES_elastics(elastic_files, energies)

def del_elastic():
	global elastic
	active_option=elastic_list.value
	for i, file in enumerate(elastic.images):
		if file == active_option[0]:
			del elastic.images[file]
			elastic.energies.pop(i)
			elastic_list.options.pop(i)
			return

def load_rois():
	global elastic
	searchfile = open("C:/Beamtimes/P01_Juli19/raw/eh2_11006060_%05d/pilatus_100k/"%int(1281) + "rois.txt", "r")
	radius=0
	box_rois=np.zeros((len(list(elastic.images)),roi_num,4),dtype=int)

	for line in searchfile:
		if "roi_radius" in line:
			radius = int(line.split(" ")[1])
			print(radius)
		if "rois_y" in line:
			rois_y = line.split(" ")[1:]
			rois_y = [int(roi) for roi in rois_y]
			for i, roi in enumerate(rois_y):
				box_rois[:,i,0]=roi-radius
				box_rois[:,i,1]=roi+radius
	searchfile = open("C:/Beamtimes/P01_Juli19/raw/eh2_11006060_%05d/pilatus_100k/"%int(1281) + "rois.txt", "r")
	i=0
	for line in searchfile:
		if "rois_x" in line:
			rois_x = line.split(" ")[1:]
			rois_x = [int(roi) for roi in rois_x]
			for j, roi in enumerate(rois_x):
				box_rois[j,i,2]=roi-radius
				box_rois[j,i,3]=roi+radius
			i+=1
	if not searchfile.closed:
		searchfile.close()
	elastic.rois=box_rois
	pic=np.zeros((fabio.open(list(elastic.images)[0]).data.shape))
	for image in elastic.images:
		pic+=fabio.open(image).data
	img=pic
	source_elastic_plot1.data["counts"]=[img]
	plot1.image(image="counts",source=source_elastic_plot1,x=0,y=0,dw=det_size[0],dh=det_size[1],palette="Viridis256")
	elastic.calc_params(mode="custom",custom_rois=box_rois)
	for rois in elastic.rois:
		for roi in rois:
			plot1.rect(x=roi[2]-radius,y=roi[1]-radius,width=radius*2,height=radius*2,fill_alpha=0,line_color="red")

def set_auto_rois():
	global elastic
	new_df=df(data=dict(
						x=[],
						y=[],
						color=[],
						labels=[]
						))
	elastic.calc_params(mode="auto",show=False)
	pic=np.zeros((fabio.open(list(elastic.images)[0]).data.shape))
	for image in elastic.images:
		pic+=fabio.open(image).data
	img=pic
	source_elastic_plot1.data["counts"]=[img]
	elastic_plot1.image(image="counts",source=source_elastic_plot1,x=0,y=0,dw=det_size[0],dh=det_size[1],palette="Viridis256")

	for rois in elastic.rois:
		for roi in rois:
			elastic_plot1.rect(x=roi[2]+5,y=roi[1]-5,width=10,height=10,fill_alpha=0,line_color="red")
			#elastic_plot1.line(x=[roi[2],roi[3]],y=[roi[1],roi[1]],color="red")
			#elastic_plot1.line(x=[roi[2],roi[3]],y=[roi[0],roi[0]],color="red")
			#elastic_plot1.line(x=[roi[2],roi[2]],y=[roi[0],roi[1]],color="red")
			#elastic_plot1.line(x=[roi[3],roi[3]],y=[roi[0],roi[1]],color="red")

	x=[]
	y=[]
	labels=[]
	for i in range(0,elastic.rois.shape[1]):
		x.append(np.linspace(0,det_size[0]-1,det_size[0]))
		y.append(np.sum(img[elastic.rois[0,i,0]:elastic.rois[0,i,1]],axis=0))
		labels.append("rois: " + str(elastic.rois[0,i,0]) + ":" + str(elastic.rois[0,i,1]))

	new_df["x"]=x
	new_df["y"]=y
	new_df["color"]=colors[:len(y)]
	new_df["labels"]=labels
	source_elastic_plot2.data=new_df
	elastic_plot2.multi_line(xs="x", ys="y", source=source_elastic_plot2,color="color",legend_field="labels")
	elastic_list.value=[]

def create_elastic():
	energies_read=str(elastic_energies.value).split(",")
	print(energies_read)
	print(elastic_path.value %int(elastic_number.value) + "*.tif") # expand to all images or in config
	energies=[]
	for energy in energies_read:
		try:
			energies.append(int(energy))
		except:
			None
	el_files=sorted(glob.glob(str(elastic_path.value) %int(elastic_number.value) + "*.tif"))
	if len(energies)==0:
		print("no energies")
		return None
	elif len(energies) != len(el_files):
		print("not the same length!")
		return None


	#print(len(glob.glob(elastic_path.value %int(elastic_number.value) + "*.tif")))
	global elastic
	print(elastic)
	elastic=XES.XES_elastics(el_files, [1,2])
	#print(elastic)
	return 999

def selected_energies(attr, old, new):
	c=0
	if len(elastic_list.value)==0:
		return
	for active_option in new:
		for i, option in enumerate(elastic_list.options):
			if active_option == option[0]:
				if c == 0:
					img =fabio.open(elastic_list.options[i][0]).data
				else:
					img+=fabio.open(elastic_list.options[i][0]).data
				c+=1
	source_elastic_plot1.data["counts"]=[img]
	elastic_plot1.image(image="counts",source=source_elastic_plot1,x=0,y=0,dw=det_size[0],dh=det_size[1],palette="Viridis256")
	y=source_elastic_plot2.data["y"]
	for i in range(0,elastic.rois.shape[1]):
		y[i]=list(np.sum(img[elastic.rois[0,i,0]:elastic.rois[0,i,1]],axis=0))
		#plot1.line(x=[0,487],y=[elastic.rois[0,i,0],elastic.rois[0,i,0]],color=colors[i])
		#plot1.line(x=[0,487],y=[elastic.rois[0,i,1],elastic.rois[0,i,1]],color=colors[i])
	source_elastic_plot2.data["y"]=y
	elastic_plot2.multi_line(xs="x", ys="y", source=source_elastic_plot2,color="color")

def update_elastic_slider(attr, old, new):
	return
#############functions for spectrum#############
def add_spectrum():
	global spectra
	global elastic
	bgr_range=(int(spect_bgrrange_slider.value[0]),int(spect_bgrrange_slider.value[1]))
	scanfiles=sorted(glob.glob(str(spect_path.value) %int(spect_number.value) + image_type))
	spect=XES.XES_object(scanfiles, elastic)
	#for file in scanfiles:
		#scan_list.options.append((file,"".join(file.split("_")[-2:])))
	if type(elastic.rois)==int:
		return
	spect.auto_rois(roiwidth=5,show=False)
	#for i in range(0,elastic.rois.shape[1]):
	#	spect.rois[i,0]=elastic.rois[0,i,0]
	#	spect.rois[i,1]=elastic.rois[0,i,1]
	spect.calc_spectrum_all(bgrwidth=spect_roiwidth_slider.value,
							bgr_range=bgr_range,
							bgrfit=spect_bgrpoly_slider.value,
							show=False)
	spect_list.options.append((spect_number.value,spect_number.value))
	spectra[int(spect_number.value)]=spect

def del_spectrum():
	global spectra
	global elastic
	if len (spect_list.value)==0:
		return
	for spect in spect_list.value:
		del spectra[int(spect)]
		spect_list.options = [option for option in spect_list.options if not spect == option[0]]
		scan_list.options=[]

def del_scan():
	global spectra
	global elastic
	if len(scan_list.value)==0:
		return
	for scan in scan_list.value:
		#print(scan)
		#print(scan_list.options)
		for option in scan_list.options:
			if option[0] == scan:
				scan_number=int(option[1].split("_")[0])
				spect = spectra[int(scan_number)]
				del spect.det_images[scan]
				scan_list.options = [option for option in scan_list.options if not scan == option[0]]
		spect.sum_images()
		spectra[int(scan_number)]=spect
	spect.calc_spectrum_all(bgrwidth=3,bgrfit=True,show=False)

def selected_spectrum(attr, old, new):
	bgr_check = False
	norm = False
	bgr_range=(int(spect_bgrrange_slider.value[0]),int(spect_bgrrange_slider.value[1]))
	if 0 in spect_checks.active:
		bgr_check = True
	if 1 in spect_checks.active:
		norm = True
	new_df_2=df(data=dict(
							x=[],
							y=[],
							bgr=[],
							color=[],
							labels=[]
						))
	new_df_3=df(data=dict(
							x=[],
							y=[],
							bgr=[],
						))
	if len(spect_list.value) ==0:
		return
	if len(spect_list.value) == 1:
		c=0
		scan_list.options=[]
		scan_list.value=[]
		xpos = np.linspace(0,det_size[0]-1,det_size[0])
		x = []
		y = []
		bgr = []
		labels=[]
		spec=spectrum()
		for active_option in spect_list.value:
			for i, option in enumerate(spect_list.options):
				if active_option == option[0]:
					spect=spectra[int(active_option)]
					img = np.array(spect.det_image)
					spec+=spect.spectrum()
					for j in range(0,spect.rois.shape[0]):
						x.append(xpos*spect.params[j][0]+spect.params[j][1])
						y.append(np.sum(np.array(spect.det_image)[spect.rois[j,0]:spect.rois[j,1]],axis=0))
						bgr.append(np.sum(spect.det_image[spect.rois[j,0]-spect_roiwidth_slider.value:spect.rois[j,0]], axis = 0)
								+np.sum(spect.det_image[spect.rois[j,1]+spect_roiwidth_slider.value:spect.rois[j,1]], axis = 0))
						scale_value=np.sum(y[j][bgr_range[0]:bgr_range[1]])/np.sum(bgr[j][bgr_range[0]:bgr_range[1]])
						#bgr[j]=bgr[j]/np.sum(bgr[j][bgr_range[0]:bgr_range[1]])*np.sum(y[j][bgr_range[0]:bgr_range[1]])
						#bgr[j]=bgr[j]/np.sum(bgr[y][bgr_range[0]:bgr_range[1]])*np.sum(y[j][bgr_range[0]:bgr_range[1]])
						bgr[j] *= scale_value
						#spec.background*=scale_value
						if spect_bgrpoly_slider.value:
							bgr[j]=spectrum(x[j],bgr[j]).pfit(spect_bgrpoly_slider.value).y
						if not bgr_check:
							y[j]-=bgr[j]
							bgr[j]=np.zeros(len(y[j]))
						if norm:
							bgr[j]/=np.max(y[j])
							y[j]/=np.max(y[j])
						#labels.append("rois " + str(spect.rois[j,0]) + ":" + str(spect.rois[j,1]) + "max: " + str(np.max(y[j])))								
						labels.append("rois " + str(spect.rois[j,0]) + ":" + str(spect.rois[j,1]))
					#else:
					#	img+= np.array(spect.det_image)
					#	spec+=spect.spectrum
					#	for j in range(0,spect.rois.shape[0]):
					#		y[j]+=np.sum(np.array(spect.det_image)[spect.rois[j,0]:spect.rois[j,1]],axis=0)
					for file in spect.det_images:
						name=file.split("_")[-2] + "_" + file.split("_")[-1]
						scan_list.options.append((file,name))

	else:
		c=0
		scan_list.options=[]
		scan_list.value=[]
		xpos = np.linspace(0,det_size[0]-1,det_size[0])
		x = []
		y = []
		bgr = []
		labels=[]
		spec=spectrum()
		for active_option in spect_list.value:
			for i, option in enumerate(spect_list.options):
				if active_option == option[0]:
					spect=spectra[int(active_option)]
					if c == 0:
						img = np.array(spect.det_image)
					else:
						img+= np.array(spect.det_image)
					spec+=spectrum(spect.spectrum.x,spect.spectrum.y)
					x.append(spect.spectrum.x)
					if norm:
						y.append(list(spect.spectrum.n().y))
						bgr.append(list(spect.spectrum.n().background))
					else:
						y.append(list(spect.spectrum.y))
						bgr.append(list(spect.spectrum.background))
					if bgr_check:
						y[c]=[y[c][index]+bgr[c][index] for index in range(0,len(y[c]))]
					else:
						bgr[c]=np.zeros(len(y[c]))
					labels.append(str(active_option))
					for file in spect.det_images:
						name=file.split("_")[-2] + "_" + file.split("_")[-1]
						scan_list.options.append((file,name))

					c+=1
	
	if norm:
		spec=spec.n()

	source_spect_plot1.data["counts"]=[img]
	spect_plot1.image(image="counts",source=source_spect_plot1,x=0,y=0,dw=det_size[0],dh=det_size[1],palette="Viridis256")

	new_df_2["x"]=x
	new_df_2["y"]=y
	new_df_2["bgr"]=bgr
	new_df_2["color"]=colors[:len(y)]
	new_df_2["labels"]=labels
	source_spect_plot2.data=new_df_2
	#spect_plot2.multi_line(xs="x",ys="y",source=source_spect_plot2,color="color",legend_field="labels")
	
	spect_plot2.legend.location = "top_left"
	spect_plot2.y_range.update()
	new_df_3["x"]=spec.x
	if bgr_check:
		new_df_3["y"]=spec.y+spec.background
		new_df_3["bgr"]=spec.background
	else:
		new_df_3["y"]=spec.y
		new_df_3["bgr"]=list(np.zeros(len(spec.y)))
	source_spect_plot3.data=new_df_3
	spect_plot3.y_range.update()
	if bgr_check:
		spect_plot2_bgr .visible=True
		spect_plot3_bgr .visible=True
	else:
		spect_plot2_bgr .visible=False
		spect_plot3_bgr .visible=False

def selected_scan(attr, old, new):
	bgr_check = False
	norm = False
	bgr_range=(int(spect_bgrrange_slider.value[0]),int(spect_bgrrange_slider.value[1]))
	if 0 in spect_checks.active:
		bgr_check = True
	if 1 in spect_checks.active:
		norm = True
	new_df_2=df(data=dict(
							x=[],
							y=[],
							bgr=[],
							color=[],
							labels=[]
						))
	new_df_3=df(data=dict(
							x=[],
							y=[],
							bgr=[],
						))
	new_title_4 = spect_plot2.title
	new_title_5 = spect_plot3.title
	if len(scan_list.value) ==0:
		return
	spect_list.value=[]
	if len(scan_list.value) == 1:
		c=0
		xpos = np.linspace(0,det_size[0]-1,det_size[0])
		x = []
		y = []
		bgr = []
		labels=[]
		#for j in range(0,len(y)):
		#	y[j]=list(np.zeros(len(y[j])))
		spec=spectrum()
		for active_option in scan_list.value:
			active_option=os.path.normpath(active_option)
			for i, option in enumerate(scan_list.options):
				if active_option == os.path.normpath(option[0]):
					spect=spectra[int(option[1].split("_")[0])]
					img = np.array(spect.det_images[active_option].data)
					spec=spect.spect_by_image[active_option]()
					for j in range(0,spect.rois.shape[0]):
						x.append(xpos*spect.params[j][0]+spect.params[j][1])
						y.append(np.sum(np.array(spect.det_images[active_option].data)[spect.rois[j,0]:spect.rois[j,1]],axis=0))
						bgr.append(np.sum(np.array(spect.det_images[active_option].data)[spect.rois[j,0]-spect_roiwidth_slider.value:spect.rois[j,0]], axis = 0)
								+np.sum(np.array(spect.det_images[active_option].data)[spect.rois[j,1]+spect_roiwidth_slider.value:spect.rois[j,1]], axis = 0))
						scale_value=np.sum(y[j][bgr_range[0]:bgr_range[1]])/np.sum(bgr[j][bgr_range[0]:bgr_range[1]])
						#bgr[j]=bgr[j]/np.sum(bgr[j][bgr_range[0]:bgr_range[1]])*np.sum(y[j][bgr_range[0]:bgr_range[1]])
						#bgr[j]=bgr[j]/np.sum(bgr[y][bgr_range[0]:bgr_range[1]])*np.sum(y[j][bgr_range[0]:bgr_range[1]])
						bgr[j] = np.multiply(bgr[j],scale_value)
						labels.append("rois " + str(spect.rois[j,0]) + ":" + str(spect.rois[j,1]))
						if spect_bgrpoly_slider.value:
							bgr[j]=spectrum(x[j],bgr[j]).pfit(spect_bgrpoly_slider.value).y
						if not bgr_check:
							y[j]=np.subtract(y[j],bgr[j])
							bgr[j]=np.zeros(len(y[j]))
						if norm:
							bgr[j]/=np.max(y[j])
							y[j]/=np.max(y[j])
	else:
		c=0
		xpos = np.linspace(0,det_size[0]-1,det_size[0])
		x = []
		y = []
		bgr = []
		labels=[]
		for j in range(0,len(y)):
			y[j]=list(np.zeros(len(y[j])))
		spec=spectrum()
		for active_option in scan_list.value:
			active_option=os.path.normpath(active_option)
			for i, option in enumerate(scan_list.options):
				if active_option == os.path.normpath(option[0]):
					spect=spectra[int(option[1].split("_")[0])]
					if c == 0:
						img = np.array(spect.det_image)
					else:
						img+= np.array(spect.det_image)
					spec+=spect.spect_by_image[active_option]()
					x.append(spec.x)
					if bgr_check:
						if norm:
							y.append(list(spect.spect_by_image[active_option].n().y))
							bgr.append(list(spect.spect_by_image[active_option].n().background))
						else:
							y.append(list(spect.spect_by_image[active_option].y))
							bgr.append(list(spect.spect_by_image[active_option].background))
						y[c]=[y[c][index]+bgr[c][index] for index in range(0,len(y[c]))]
					else:
						if norm:
							y.append(list(spect.spect_by_image[active_option].n().y))
						else:
							y.append(list(spect.spect_by_image[active_option].y))
						bgr.append(np.zeros(len(y[c])))
					labels.append(str(option[1]))
					c+=1		
	if norm:
		spec=spec.n()
	source_spect_plot1.data["counts"]=[np.array(img)]

	new_df_2["x"]=x
	new_df_2["y"]=y
	new_df_2["bgr"]=bgr
	new_df_2["color"]=colors[:len(y)]
	new_df_2["labels"]=labels
	source_spect_plot2.data=new_df_2
	#spect_plot1.image(image="counts",source=source_spect_plot1,x=0,y=0,dw=det_size[0],dh=det_size[1],palette="Viridis256")
	#spect_plot2.multi_line(xs="x",ys="y",source=source_spect_plot2,color="color",legend_field="labels")
	spect_plot2.legend.location = "top_left"
	spect_plot2.y_range.update()
	new_df_3["x"]=spec.x
	if bgr_check:
		new_df_3["y"]=spec.y+spec.background
	else:
		new_df_3["y"]=spec.y
	new_df_3["bgr"]=spec.background
	source_spect_plot3.data=new_df_3
	#spect_plot3.line(x="x",y="y",source=source_spect_plot3)
	spect_plot3.y_range.update()
	if bgr_check:
		spect_plot2_bgr.visible=True
		spect_plot3_bgr.visible=True
	else:
		spect_plot2_bgr.visible=False
		spect_plot3_bgr.visible=False

def update_roiwidth(attr, old, new):
	global spectra
	for active_option in spect_list.value:
		active_option=int(active_option)
		spectra[active_option].auto_rois(roiwidth=spect_roiwidth_slider.value,show=False)
		spectra[active_option].calc_spectrum_all(bgrwidth=3,bgrfit=spect_bgrpoly_slider.value,show=False)
	selected_spectrum(0,0,0)
	return

def update_bgrfit(attr, old, new):
	global spectra
	for active_option in spect_list.value:
		active_option=int(active_option)
		spectra[active_option].auto_rois(roiwidth=spect_roiwidth_slider.value,show=False)
		spectra[active_option].calc_spectrum_all(bgrwidth=3,bgrfit=spect_bgrpoly_slider.value,show=False)
	selected_spectrum(0,0,0)
	return

def update_spect_slider(attr, old, new):
	global spectra
	bgr_range=(int(spect_bgrrange_slider.value[0]),int(spect_bgrrange_slider.value[1]))
	
	for active_option in spect_list.value:
		active_option=int(active_option)
		spectra[active_option].auto_rois(roiwidth=spect_roiwidth_slider.value,show=False)
		spectra[active_option].calc_spectrum_all(	bgrwidth=3,
													bgr_range=bgr_range,
													bgrfit=spect_bgrpoly_slider.value,
													show=False)
	if len(spect_list.value) > 0:
		selected_spectrum(0,0,0)
	else:
		selected_scan(0,0,0)
	return

def update_checkbox(new):
	global spectra
	bgr_range=(int(spect_bgrrange_slider.value[0]),int(spect_bgrrange_slider.value[1]))
	for active_option in spect_list.value:
		active_option=int(active_option)
		spectra[active_option].auto_rois(roiwidth=spect_roiwidth_slider.value,show=False)
		spectra[active_option].calc_spectrum_all(	bgrwidth=3,
													bgr_range=bgr_range,
													bgrfit=spect_bgrpoly_slider.value,
													show=False)
	selected_spectrum(0,0,0)
	selected_scan(0,0,0)
	return

def switch():
	if switch_button.active==0:
		curdoc().roots[0].children[1]=elastic_layout
	elif switch_button.active==1:
		curdoc().roots[0].children[1]=spect_layout
	else:
		print("Gedulden du dich musst junger Padawan")
	curdoc()
	curdoc().roots[0].update()

switch_button = RadioButtonGroup(labels=['elastic', 'spectrum','exported'], active=0)
switch_button.on_change('active', lambda attr, old, new: switch())

##########set up elastic col1##########
elastic_path= TextInput(title="Path of elastics",value=sample_prefix)
elastic_number= TextInput(title="elasticnumber",value="1281")
elastic_add= Button(label="add")
elastic_add.on_click(add_elastic)
elastic_energies= TextInput(title="energies", placeholder="E0,E1, ...",value="1,2,3,4,5,6,7,8,9,10")
elastic_list = MultiSelect(size=6,title="files")
elastic_list.on_change('value',selected_energies)
elastic_del= Button(label="delete")
elastic_del.on_click(del_elastic)
elastic_auto_rois= Button(label="auto rois")
elastic_auto_rois.on_click(set_auto_rois)
elastic_load_rois= Button(label="load rois")
elastic_load_rois.on_click(load_rois)

##########set up spectrum col1##########
spect_path = TextInput(title="Path of spectrum",value=sample_prefix)
spect_number= TextInput(title="Scannumber",value="1287")
spect_add = Button(label="add spectrum")
spect_add.on_click(add_spectrum)
spect_list = MultiSelect(size=6,title="spectra")
spect_list.on_change('value',selected_spectrum)
spect_del = Button(label="delete spectrum")
spect_del.on_click(del_spectrum)
#spect_label= TextInput(title="Name")
scan_list = MultiSelect(size=6,title="scan files")
scan_list.on_change('value',selected_scan)
scan_del = Button(label="delete scan")
scan_del.on_click(del_scan)

##########set up plots elastic col2##########
hover = HoverTool(tooltips=[("(x,y)", "($x, $y)"),
    						("count", "@counts")])

source_elastic_plot1=ColumnDataSource(data=dict	(
											#counts=[img_1]
											counts=[]
										))

elastic_plot1 = figure(plot_width=600, plot_height=300,lod_threshold = 10)
#elastic_plot1.tools=[PanTool(),BoxZoomTool(),WheelZoomTool(),HoverTool(),ResetTool()]
elastic_plot1.add_tools(hover)
elastic_plot1.image(image="counts",source=source_elastic_plot1,x=0,y=0,dw=det_size[0],dh=det_size[1],palette="Viridis256")

source_elastic_plot2=ColumnDataSource(data=dict	(
											#x=range(0,487),
											#y=list()np.zeros((487,4))
											x=[],
											y=[],
											color=[]
										))
elastic_plot2 = figure(plot_width=600, plot_height=300,lod_threshold = 10)

##########set up spectrum col2##########
source_spect_plot1=ColumnDataSource(data=dict	(
											counts=[]
										))
spect_plot1 = figure(plot_width=600, plot_height=300,lod_threshold = 10)
spect_plot1.add_tools(hover)
spect_plot1.image(image="counts",source=source_spect_plot1,x=0,y=0,dw=det_size[0],dh=det_size[1],palette="Viridis256")

source_spect_plot2=ColumnDataSource(data=dict(
												x=[],
												y=[],
												bgr=[],
												color=[],
												labels=[]
												))
spect_plot2 = figure(plot_width=600, plot_height=300,lod_threshold = 10, output_backend="webgl")
spect_plot2_lines = spect_plot2.multi_line(xs="x",ys="y",source=source_spect_plot2,color="color",legend_field="labels")
spect_plot2_bgr = spect_plot2.multi_line(xs="x",ys="bgr",source=source_spect_plot2,color="color")
#spect_plot2 = figure(plot_width=600, plot_height=300, output_backend="webgl")

source_spect_plot3=ColumnDataSource(data=dict	(
											#x=range(0,487),
											#y=list()np.zeros((487,4))
											x=[],
											y=[],
											bgr=[],
											color=[]
										))
spect_plot3 = figure(plot_width=600, plot_height=300,lod_threshold = 10)
spect_plot3_line = spect_plot3.line(x="x",y="y",source=source_spect_plot3)
spect_plot3_bgr = spect_plot3.line(x="x",y="bgr",source=source_spect_plot3)

##########set up elastic col3##########
elastic_roiradius_slider = Slider(	start=1,
									end = 10,
									value = 5,
									step = 1 ,
									callback_policy="mouseup")
elastic_roiradius_slider.on_change('value_throttled',update_elastic_slider)

##########set up spectrum col3##########
spect_roiwidth_slider = Slider(	start=1,
								end = 10,
								value = 3,
								step = 1 ,
								title="Roiwidth",
								callback_policy="mouseup")
spect_roiwidth_slider.on_change("value_throttled",update_spect_slider)
spect_checks = CheckboxGroup(
							labels =["show background", "normalize"],
							active=[0,0])
spect_checks.on_click(update_checkbox)
spect_bgrrange_slider = RangeSlider(start=0,
									end = det_size[0],
									value = (0,int(det_size[0]/10)),
									step = 1 ,
									title="background shift (pixel)",
									callback_policy="mouseup")
spect_bgrrange_slider.on_change("value_throttled", update_spect_slider)
spect_bgrpoly_slider = Slider(	start=0,
								end = 20,
								value = 10,
								step = 1,
								title="background polynom (0 = no fit)",
								callback_policy="mouseup")
spect_bgrpoly_slider.on_change("value_throttled", update_spect_slider)

########################################
'''
Button load config
Slider BGR width
Button show Rois
Save Spectrum (alles oder nur x-y-werte?)
Something with references
plot title check ???
merge scans button or save in list (if len(selected)>1: merge else: just save x,y)
'''
########################################

'''
test_button= Button(label="test")
test_button.on_click(test)


elastic_col1 = column(	
				elastic_path,
				elastic_number,
				elastic_energies,
				elastic_add,
				elastic_list,
				elastic_del,
				elastic_auto_rois,
				elastic_load_rois
				)
elastic_col2 = column(plot1,plot2)
'''



elastic_col1 = column(	
					elastic_path,
					elastic_number,
					elastic_energies,
					elastic_add,
					elastic_list,
					elastic_del,
					elastic_auto_rois,
					elastic_load_rois
					)
elastic_col2 = column(
					elastic_plot1,
					elastic_plot2
					)
#elastic_col3 = column(
#					elastic_roiradius_slider
#					)


spect_col1 = column(	
					spect_path,
					spect_number,
					spect_add,
					spect_list,
					spect_del,
					scan_list,
					scan_del
					)
spect_col2 = column(
					spect_plot1,
					spect_plot2,
					spect_plot3
					)
spect_col3 = column(
					spect_roiwidth_slider,
					spect_checks,
					spect_bgrrange_slider,
					spect_bgrpoly_slider
					)

#elastic_layout = row(
#					elastic_col1,
#					elastic_col2,
#					elastic_col3
#					)
elastic_layout = row(
					elastic_col1,
					elastic_col2
					)
spect_layout = row(	
					spect_col1,
					spect_col2,
					spect_col3
					)
#roi_num= TextInput(title="number of crystals",value="4")
roots=curdoc().roots
#col3 = column(roi_num) 
#curdoc().add_root(column(row(elastic_col1,elastic_col2),row(spectrum_col1,spectrum_col2)))

curdoc().add_root(column(switch_button,elastic_layout))
curdoc().title = "XES_elastics"


roots=curdoc().roots

add_elastic()
set_auto_rois()
add_spectrum()
