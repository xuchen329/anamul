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

def DCRvsVolt(folderdir):
    filename = folderdir+'/spes.log'
    if os.path.isfile(filename):
        volt,temp,DCR,errDCR = np.genfromtxt(filename,usecols=(0,1,9,10),unpack=True)
        volt = -1.*volt
        plt.errorbar(volt[DCR.nonzero()],DCR[DCR.nonzero()],yerr=errDCR[DCR.nonzero()],fmt='.')

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
   # plt.annotate("Temp: {0:.1f}$^\circ$C".format(np.mean(temp)),xy=(0.1,0.6),xycoords='axes fraction',fontsize=16)
    plt.xlabel('Bias Voltage [V]',fontsize=16)
    plt.ylabel('DCR [cps]',fontsize=16)
    plt.xlim(65.95,68.25)
    plt.ylim(1e5,3e6)
    plt.show()

if __name__=='__main__':
    folderpath = sys.argv[1]
    matrixno = sys.argv[2]
    if folderpath[-1]!='/':
        folderpath+='/'
    DoFitting(folderpath,matrixno)
