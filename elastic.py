#import sys
from spect import *
import matplotlib.pyplot as plt
import numpy as np
import glob
from PIL import Image
import collections
from scipy.optimize import curve_fit

data = dict()
pilatusimages = dict()
sumxtal1 = dict()
sumxtal2 = dict()
sumxtal3 = dict()
crystall = dict()
elastic = dict()
#elasticmax1 = dict()
params=dict()
t=dict()

def load(path, elasticnumber, energystart, energystep, rois):
    elasticmax = dict()
    for roiscount in (range(0,len(rois))):
        elasticmax[roiscount]=[]
    xpos = np.linspace(1, 487, 487)
    data = np.zeros((195, 487))
    crystall = np.zeros((len(rois),) + (487,))
    bckround = np.zeros(np.shape(crystall))
    crystall_WBCK  = np.zeros(np.shape(crystall))
    spect=spectrum()
    T = "Scan" +("%05d" % elasticnumber)
    i=0 
    for file in glob.glob(path +  ("%05d" % elasticnumber) + "/pilatus_100k/alignment_" + "*.tif")[::int(1)]: 
        pilatusimages = Image.open(file)
        data = np.array(pilatusimages)

        #data[46][18] = 3   #correct dead pixel to bgr value, only BL9 DELTA Pilatus 100K s specific
        for (x1,x2),roiscount in zip(rois,range(0,len(rois))):
            crystall[roiscount] = np.array(np.sum(data[x1:x2], axis = 0))/(len(data[x1:x2]))
            p=spectrum(xpos[np.argmax(crystall[roiscount])-10:np.argmax(crystall[roiscount])+10:1], crystall[roiscount][np.argmax(crystall[roiscount])-10:np.argmax(crystall[roiscount])+10:1]).gauss(fullinfo=1)
            elasticmax[roiscount].append([p.cms(),energystart+i*energystep,np.abs(p.fit_p[2]*2*np.sqrt(2*np.log(2)))])
        i+=1
        print(i)
    for c in (range(0,roiscount+1)):
        print(np.array(elasticmax[c]))
        params[c]=fwhm.fit2(np.array(elasticmax[c])[:,0],np.array(elasticmax[c])[:,1],lin)
        t[c]=sum(np.array(elasticmax[c])[:,2])/len(np.array(elasticmax[c])[:,2])
        #params[c][0].append(p)
    #print("Summe")
    #print(t)
        #params[c].append()
    print(params)
    print(t)
    out=[]
    for key in params:
        out.append([params[key][0][0], params[key][0][1],t[key]])
    return out
def lin(x,a,b):
    return a*x+b
