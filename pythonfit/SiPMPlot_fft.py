#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.optimize as scop
import scipy.stats as scist
import os,sys
import ROOT

rc('text',usetex=True)
rc('font', family='serif')
def linearfunc(x,cnt,slp):
    return cnt+slp*x

def GainvsVolt(folderpath):
    filename = folderpath+'/spes-fft.log'
    if os.path.isfile(filename):
        voltage,temperature,GainFFT16,GainFFT16err,GainFFT,GainFFTerr = np.genfromtxt(filename,usecols=(0,1,2,3,4,5),unpack=True)
        if voltage[0]<0.:
            voltage = -1.*voltage
        mygrFFT16 = ROOT.TGraphErrors(len(voltage),voltage.flatten(),GainFFT16.flatten(),np.zeros(len(voltage),dtype=float),GainFFT16err.flatten())
        c1 = ROOT.TCanvas("FFT64_fit","FFT64_Fit")
        c1.cd()
        mygrFFT16.Draw("AP")
        myfitFFT16 = ROOT.TF1("myfitFFT16","pol1",np.min(voltage),np.max(voltage))
        #myfitFFT16 = ROOT.TF1("myfitFFT16","pol1",71.3,71.8)
        mygrFFT16.Fit("myfitFFT16","QR")
        FFT16par0 = myfitFFT16.GetParameter(0)
        FFT16par1 = myfitFFT16.GetParameter(1)
        FFT16err0 = myfitFFT16.GetParError(0)
        FFT16err1 = myfitFFT16.GetParError(1)
        FFT16chi = myfitFFT16.GetChisquare()
        FFT16dof = myfitFFT16.GetNDF()
        plt.errorbar(voltage,GainFFT16,yerr=GainFFT16err,fmt='b.')
        FFT16_vbd = -1.*FFT16par0/FFT16par1
        FFT16_errvbd = np.sqrt(((FFT16err0**2)/(FFT16par0**2)+(FFT16err1**2)/(FFT16par1**2))*(FFT16_vbd**2))
        FFT16_normg = FFT16par1*25e-15/50./1.6e-19
        FFT16_errnormg = FFT16err1*25e-15/50./1.6e-19
        plt.plot(voltage,linearfunc(voltage,FFT16par0,FFT16par1),'b--')
        voltageerr = np.linspace(np.std(temperature)*0.056,np.std(temperature)*0.056,len(voltage))
        mygrFFT = ROOT.TGraphErrors(len(voltage),voltage.flatten(),GainFFT.flatten(),np.zeros(len(voltage),dtype=float),GainFFTerr.flatten())
        mygrFFT = ROOT.TGraphErrors(len(voltage),voltage.flatten(),GainFFT.flatten(),voltageerr.flatten(),GainFFTerr.flatten())
        c2 = ROOT.TCanvas("FFT_fit","FFT_Fit")
        c2.cd()
        mygrFFT.Draw("AP")
        myfitFFT = ROOT.TF1("myfitFFT16","pol1",np.min(voltage),np.max(voltage))
        #myfitFFT = ROOT.TF1("myfitFFT16","pol1",71.3,71.8)
        mygrFFT.Fit("myfitFFT16","QR")
        FFTpar0 = myfitFFT.GetParameter(0)
        FFTpar1 = myfitFFT.GetParameter(1)
        FFTerr0 = myfitFFT.GetParError(0)
        FFTerr1 = myfitFFT.GetParError(1)
        FFTchi = myfitFFT.GetChisquare()
        FFTdof = myfitFFT.GetNDF()
        plt.errorbar(voltage,GainFFT,yerr=GainFFTerr,fmt='r.')
        FFT_vbd = -1.*FFTpar0/FFTpar1
        FFT_errvbd = np.sqrt(((FFTerr0**2)/(FFTpar0**2)+(FFTerr1**2)/(FFTpar1**2))*(FFT_vbd**2))
        FFT_normg = FFTpar1*25e-15/50./1.6e-19
        FFT_errnormg = FFTerr1*25e-15/50./1.6e-19
        plt.plot(voltage,linearfunc(voltage,FFTpar0,FFTpar1),'r--')
        plt.grid(True)
        plt.xlabel('Bias Voltage [V]',fontsize=16)
        plt.ylabel('Gain [adu.]',fontsize=16)
        plt.annotate("FFText\n$U_{bd}$ : "+"{0:.2f} $\pm$ {1:.2f} V\nGain: {2:.1f} $\pm$ {6:.1f} k$e_0$/V\nTemp: {3:.1f}$^\circ$C\n$\chi^2$/DOF : {4:.1f}/{5}".format(FFT16_vbd,FFT16_errvbd,FFT16_normg/1000.,np.mean(temperature),FFT16chi,FFT16dof,FFT16_errnormg/1000.),xy=(0.35,0.75),xycoords='axes fraction',fontsize=14)
        plt.annotate("FFT\n$U_{bd}$ : "+"{0:.2f} $\pm$ {1:.2f} V\nGain: {2:.1f} $\pm$ {6:.1f} k$e_0$/V\nTemp: {3:.1f}$^\circ$C\n$\chi^2$/DOF : {4:.1f}/{5}".format(FFT_vbd,FFT_errvbd,FFT_normg/1000.,np.mean(temperature),FFTchi,FFTdof,FFT_errnormg/1000.),xy=(0.05,0.75),xycoords='axes fraction',fontsize=14)
        plt.xlim(np.min(voltage)-0.1,np.max(voltage)+0.1)
    plt.show()

if __name__=='__main__':
    folderpath = sys.argv[1]
    
    GainvsVolt(folderpath)
