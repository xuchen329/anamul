#!/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.optimize as scop
import scipy.stats as scist
import os,sys
import ROOT

FolderPre = 'Matrix-'
NamePre = '/spes-M'
CHName = ['A','B','C','D']

golden_mean = (math.sqrt(5)-1.0)/2.0
fig_width = 8
fig_height = fig_width*golden_mean
fig_size = [fig_width,fig_height]
params = {'font.family':'serif',
	  'figure.figsize': fig_size}
plt.rcParams.update(params)
rc('text',usetex=True)

def DoFitting(folderpath,matrixno):
    filefullname = folderpath
    filefullname+=FolderPre
    filefullname+=str(matrixno)
    filefullname+=NamePre
    filefullname+=str(matrixno)
    for i in range(4):
        for j in range(4):
            casefilename = filefullname+'-'+CHName[i]+str(j+1)
            cmd = './bin/getfit-hn '+casefilename+' 25'
            os.system(cmd)

gltemp = 0.
def DCRvsVolt(folderdir):
    if folderdir[-2]=='A' or folderdir[-2]=='B':
	fomt = 'b.'
    else:
	fomt = 'r.'
    filename = folderdir+'/spes.log'
    if os.path.isfile(filename):
        volt,temp,DCR,errDCR = np.genfromtxt(filename,usecols=(0,1,9,10),unpack=True)
        volt = -1.*volt
	DCR = DCR/1e6
	errDCR = errDCR/1e6
        plt.errorbar(volt[DCR.nonzero()],DCR[DCR.nonzero()],yerr=errDCR[DCR.nonzero()],label=folderdir,fmt=fomt)
	global gltemp
	gltemp+=np.mean(temp)

def DoDCRPlot(folderpath,matrixno):
    filefullname = folderpath
    filefullname+=FolderPre
    filefullname+=str(matrixno)
    filefullname+=NamePre
    filefullname+=str(matrixno)
    for i in range(4):
        for j in range(4):
            casefilename = filefullname+'-'+CHName[i]+str(j+1)
            DCRvsVolt(casefilename)
    plt.grid(True)
    global gltemp
    plt.annotate("Temp: {0:.1f}$^\circ$C".format((gltemp)/16.),xy=(0.1,0.6),xycoords='axes fraction',fontsize=16)
    plt.xlabel('Bias Voltage [V]',fontsize=16)
    plt.ylabel('DCR [Mcps]',fontsize=16)
    plt.xlim(65.95,68.25)
    plt.ylim(0.1,3)
    plt.legend()
    plt.show()

if __name__=='__main__':
    argc = len(sys.argv)
    folderpath = sys.argv[1]
    matrixno = sys.argv[2]
    if folderpath[-1]!='/':
        folderpath+='/'
    DoDCRPlot(folderpath,matrixno) 
