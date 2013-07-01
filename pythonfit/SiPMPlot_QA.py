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
    filename = folderpath+'/spes.log'
    if os.path.isfile(filename):
        voltage,temperature,GainMG,errGainMG,GainFFT,GainFFTerr = np.genfromtxt(filename,usecols=(0,1,4,6,7,8),unpack=True)
        if voltage[0]<0.:
            voltage = -1.*voltage
        voltageerr = np.linspace(np.std(temperature)*0.06,np.std(temperature)*0.06,len(voltage))
        
#plot mgfit
        mygrMG = ROOT.TGraphErrors(len(voltage[GainMG.nonzero()]),voltage[GainMG.nonzero()].flatten(),GainMG[GainMG.nonzero()].flatten(),voltageerr[GainMG.nonzero()].flatten(),errGainMG[GainMG.nonzero()].flatten())
        c1 = ROOT.TCanvas("MG_fit","MG_Fit")
        c1.cd()
        mygrMG.Draw("AP")
        myfitMG = ROOT.TF1("myfitMG","pol1",np.min(voltage),np.max(voltage))
        mygrMG.Fit(myfitMG,"QR")
        MGpar0 = myfitMG.GetParameter(0)
        MGpar1 = myfitMG.GetParameter(1)
        MGerr0 = myfitMG.GetParError(0)
        MGerr1 = myfitMG.GetParError(1)
        chiMG = myfitMG.GetChisquare()
        dofMG = myfitMG.GetNDF()
        plt.errorbar(voltage[GainMG.nonzero()],GainMG[GainMG.nonzero()],xerr=voltageerr[GainMG.nonzero()],yerr=errGainMG[GainMG.nonzero()],fmt='b.')
        vbdMG = -1.*MGpar0/MGpar1
        errvbdMG = np.sqrt(((MGerr0**2)/(MGpar0**2)+(MGerr1**2)/(MGpar1**2))*(vbdMG**2))
        normgMG = MGpar1*25e-15/50./1.6e-19
        errnormgMG = MGerr1*25e-15/50./1.6e-19
        plt.plot(voltage,linearfunc(voltage,MGpar0,MGpar1),'b--')

#plot fft
#mygrFFT = ROOT.TGraphErrors(len(voltage),voltage.flatten(),GainFFT.flatten(),np.zeros(len(voltage),dtype=float),GainFFTerr.flatten())
        mygrFFT = ROOT.TGraph(len(voltage),voltage.flatten(),GainFFT.flatten())
        c2 = ROOT.TCanvas("FFT_fit","FFT_Fit")
        c2.cd()
        mygrFFT.Draw("AP")
        myfitFFT = ROOT.TF1("myfitFFT","pol1",np.min(voltage),np.max(voltage))
        #myfitFFT = ROOT.TF1("myfitFFT16","pol1",71.3,71.8)
        mygrFFT.Fit(myfitFFT,"QR")
        FFTpar0 = myfitFFT.GetParameter(0)
        FFTpar1 = myfitFFT.GetParameter(1)
        FFTerr0 = myfitFFT.GetParError(0)
        FFTerr1 = myfitFFT.GetParError(1)
        FFTchi = myfitFFT.GetChisquare()
        FFTdof = myfitFFT.GetNDF()
        plt.errorbar(voltage,GainFFT,xerr=voltageerr,yerr=GainFFTerr,fmt='r.')
        FFT_vbd = -1.*FFTpar0/FFTpar1
        FFT_errvbd = np.sqrt(((FFTerr0**2)/(FFTpar0**2)+(FFTerr1**2)/(FFTpar1**2))*(FFT_vbd**2))
        FFT_normg = FFTpar1*25e-15/50./1.6e-19
        FFT_errnormg = FFTerr1*25e-15/50./1.6e-19
        plt.plot(voltage,linearfunc(voltage,FFTpar0,FFTpar1),'r--')

#global cosmetic setting
        plt.grid(True)
        plt.xlabel('Bias Voltage [V]',fontsize=16)
        plt.ylabel('Gain [adu.]',fontsize=16)
        
        plt.annotate("MGF\n$U_{bd}$ : "+"{0:.2f} $\pm$ {1:.2f} V\nGain: {2:.1f} $\pm$ {6:.1f} k$e_0$/V\nTemp: {3:.1f}$^\circ$C\n$\chi^2$/DOF : {4:.1f}/{5}".format(vbdMG,errvbdMG,normgMG/1000.,np.mean(temperature),chiMG,dofMG,errnormgMG/1000.),xy=(0.45,0.7),xycoords='axes fraction',fontsize=14)
        plt.annotate("FFT\n$U_{bd}$ : "+"{0:.2f} $\pm$ {1:.2f} V\nGain: {2:.1f} $\pm$ {6:.1f} k$e_0$/V\nTemp: {3:.1f}$^\circ$C\n$\chi^2$/DOF : {4:.1f}/{5}".format(FFT_vbd,FFT_errvbd,FFT_normg/1000.,np.mean(temperature),FFTchi,FFTdof,FFT_errnormg/1000.),xy=(0.05,0.7),xycoords='axes fraction',fontsize=14)
         
        plt.xlim(np.min(voltage)-0.1,np.max(voltage)+0.1)
    plt.show()

if __name__=='__main__':
    folderpath = sys.argv[1]
    
    GainvsVolt(folderpath)