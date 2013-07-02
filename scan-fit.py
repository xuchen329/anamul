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

if __name__=='__main__':
    folderpath = sys.argv[1]
    matrixno = sys.argv[2]
    if folderpath[-1]!='/':
        folderpath+='/'
    DoFitting(folderpath,matrixno)
