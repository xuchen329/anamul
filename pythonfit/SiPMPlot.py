#!/usr/bin/env python
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.optimize as scop
import scipy.stats as scist
import os,sys
import ROOT

golden_mean = (math.sqrt(5)-1.0)/2.0
fig_width = 8 # width in inches
fig_height = fig_width*golden_mean
fig_size = [fig_width, fig_height]
params = { 'font.family':'serif',
	'figure.figsize':fig_size}
plt.rcParams.update(params)
rc('text',usetex=True)

def linearfunc(x,cnt,slp):
    return cnt+slp*x
    
def GainvsVolt(folderpath):
    filename = folderpath+'/pyspes.log'
    fftplt=1;
    if 'pyspes.log' in filename:
        if os.path.isfile(filename):
            gain,errgain,volt,temp = np.genfromtxt(filename,usecols=(6,7,12,13),unpack=True,dtype='float')
            volt = -1.*volt
            mygr = ROOT.TGraphErrors(len(volt),volt.flatten(),gain.flatten(),np.zeros(len(volt),dtype=float),errgain.flatten())
            mygr.Draw("AP")
            for i in range(len(errgain)):
                if errgain[i]<1:
                    errgain[i] = 1
            popt,perr = scop.curve_fit(linearfunc,volt,gain,sigma=errgain)
            dof = len(gain)-1-2
            chi,pval = scist.chisquare(gain,linearfunc(volt,popt[0],popt[1]),dof)
            plt.errorbar(volt,gain,yerr=errgain,fmt='.')
            vbd = -1.*popt[0]/popt[1]
            errvbd = np.sqrt((perr[0][0]/(popt[0]**2)+perr[1][1]/(popt[1]**2))*(vbd**2))
            normg = popt[1]*20e-15/50./1.6e-19
            errnormg = np.sqrt(perr[1][1])*20e-15/50./1.6e-19
            plt.plot(volt,linearfunc(volt,popt[0],popt[1]),'r--')
            plt.grid(True)
            plt.xlabel('Bias Voltage [V]',fontsize=16)
            plt.ylabel('Gain [adu.]',fontsize=16)
            plt.annotate("DSF\n$U_{bd}$ : "+"{0:.2f} $\pm$ {1:.2f} V\nGain: {2:.2e} $e_0$/V\nTemp: {3:.1f}$^\circ$C\n$\chi^2$/DOF : {4:.1f}/{5}".format(vbd,errvbd,normg,np.mean(temp),chi,dof),xy=(0.6,0.3),xycoords='axes fraction',fontsize=16)
            plt.xlim(np.min(volt)-0.1,np.max(volt)+0.1)
            plt.show()
        else:
            print 'pyspes.log not exists, using FFT fit data ...'
            fftplt=1
    if (fftplt):
        filename = folderpath+'/spes.log'
        if os.path.isfile(filename):
            voltage,temperature,GainFFT,errGainFFT = np.genfromtxt(filename,usecols=(0,1,7,8),unpack=True)
            voltage = -1.*voltage
            popt,perr = scop.curve_fit(linearfunc,voltage,GainFFT,sigma=errGainFFT)
            dof = len(GainFFT)-1-2
            chi,pval = scist.chisquare(GainFFT,linearfunc(voltage,popt[0],popt[1]),dof)
            plt.errorbar(voltage,GainFFT,yerr=errGainFFT,fmt='.')
            vbd = -1.*popt[0]/popt[1]
            errvbd = np.sqrt((perr[0][0]/(popt[0]**2)+perr[1][1]/(popt[1]**2))*(vbd**2))
            normg = popt[1]*20e-15/50./1.6e-19
            errnormg = np.sqrt(perr[1][1])*20e-15/50./1.6e-19
            plt.plot(voltage,linearfunc(voltage,popt[0],popt[1]),'r--')
            plt.grid(True)
            plt.xlabel('Bias Voltage [V]',fontsize=16)
            plt.ylabel('Gain [adu.]',fontsize=16)
            plt.annotate("FFT\n$U_{bd}$ : "+"{0:.2f} $\pm$ {1:.2f} V\nGain: {2:.2e} $e_0$/V\nTemp: {3:.1f}$^\circ$C\n$\chi^2$/DOF : {4:.1f}/{5}".format(vbd,errvbd,normg,np.mean(temperature),chi,dof),xy=(0.6,0.3),xycoords='axes fraction',fontsize=16)
            plt.xlim(voltage[0]-0.1,voltage[-1]+0.1)
            plt.show()
        else:
            print "no log file found!"

def XTvsVolt(folderdir):
    filename = folderdir+'/pyspes.log'
    if os.path.isfile(filename):
        xt,errxt,volt,temp = np.genfromtxt(filename,usecols=(10,11,12,13),unpack=True)
        volt = -1.*volt
        xt = xt*100
        errxt = errxt*100
        plt.errorbar(volt,xt,yerr=errxt,fmt='r.')
        plt.grid(True)
        plt.annotate("Temp: {0:.1f}$^\circ$C".format(np.mean(temp)),xy=(0.6,0.3),xycoords='axes fraction',fontsize=16)
        plt.xlabel('Bias Voltage [V]',fontsize=16)
        plt.ylabel('XT prop [\%]',fontsize=16)
        plt.xlim(np.min(volt)-0.1,np.max(volt)+0.1)
    filename2 = folderdir+'/spes.log'
    if os.path.isfile(filename):
        voltage,temperature,crosstalk,errcrosstalk = np.genfromtxt(filename2,usecols=(0,1,11,12),unpack=True)
        voltage = -1.*voltage
        crosstalk = 100.*crosstalk
        errcrosstalk = 100.*errcrosstalk
        plt.errorbar(voltage,crosstalk,yerr=errcrosstalk,fmt='b.')
    plt.legend(('Detector response fit','DCR measurement'))
    plt.show()


def DCRvsVolt(folderdir):
    filename = folderdir+'/spes.log'
    if os.path.isfile(filename):
        volt,temp,DCR,errDCR = np.genfromtxt(filename,usecols=(0,1,9,10),unpack=True)
        volt = -1.*volt
        plt.errorbar(volt[DCR.nonzero()],DCR[DCR.nonzero()],yerr=errDCR[DCR.nonzero()],fmt='r.')
        plt.grid(True)
        plt.annotate("Temp: {0:.1f}$^\circ$C".format(np.mean(temp)),xy=(0.1,0.6),xycoords='axes fraction',fontsize=16)
        plt.xlabel('Bias Voltage [V]',fontsize=16)
        plt.ylabel('DCR [cps]',fontsize=16)
        plt.xlim(np.min(volt)-0.1,np.max(volt)+0.1)
        plt.ylim(1e5,3e6)
        #print volt-69.57
        plt.show()

if __name__=='__main__':
    folderpath = sys.argv[1]
    
    #GainvsVolt(folderpath)
    #XTvsVolt(folderpath)
    DCRvsVolt(folderpath)
