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
    filename = folderpath+'/pyspes.log'
    if os.path.isfile(filename):
        gain,errgain,volt,temp = np.genfromtxt(filename,usecols=(6,7,12,13),unpack=True,dtype='float')
        volt = -1.*volt
        mygr = ROOT.TGraphErrors(len(volt),volt.flatten(),gain.flatten(),np.zeros(len(volt),dtype=float),errgain.flatten())
        c1 = ROOT.TCanvas("DR_fit","DR_fit")
        c1.cd()
        mygr.Draw("AP")
        myfit = ROOT.TF1("myfit","pol1",np.min(volt),np.max(volt))
        mygr.Fit("myfit","QR")
        DRpar0 = myfit.GetParameter(0)
        DRpar1 = myfit.GetParameter(1)
        DRerr0 = myfit.GetParError(0)
        DRerr1 = myfit.GetParError(1)
        chiDR = myfit.GetChisquare()
        dofDR = myfit.GetNDF()
        plt.errorbar(volt,gain,yerr=errgain,fmt='c.')
        vbd = -1.*DRpar0/DRpar1
        errvbd = np.sqrt(((DRerr0**2)/(DRpar0**2)+(DRerr1**2)/(DRpar1**2))*(vbd**2))
        normg = DRpar1*20e-15/50./1.6e-19
        errnormg = DRerr1*20e-15/50./1.6e-19
        plt.plot(volt,linearfunc(volt,DRpar0,DRpar1),'r-')
        plt.grid(True)
        plt.xlabel('Bias Voltage [V]',fontsize=16)
        plt.ylabel('Gain [adu.]',fontsize=16)
        plt.annotate("DRF\n$U_{bd}$ : "+"{0:.2f} $\pm$ {1:.2f} V\nGain: {2:.2e} $e_0$/V\nTemp: {3:.1f}$^\circ$C\n$\chi^2$/DOF : {4:.1f}/{5}".format(vbd,errvbd,normg,np.mean(temp),chiDR,dofDR),xy=(0.05,0.75),xycoords='axes fraction',fontsize=14)
        plt.xlim(np.min(volt)-0.1,np.max(volt)+0.1)
    filename = folderpath+'/spes.log'
    if os.path.isfile(filename):
        voltage,temperature,GainMG,errGainMG,GainFFT,errGainFFT = np.genfromtxt(filename,usecols=(0,1,4,6,7,8),unpack=True)
        voltage = -1.*voltage
        mygrMG = ROOT.TGraphErrors(len(voltage),voltage.flatten(),GainMG.flatten(),np.zeros(len(volt),dtype=float),errGainMG.flatten())
        c2 = ROOT.TCanvas("MG_fit","MG_Fit")
        c2.cd()
        mygrMG.Draw("AP")
        myfitMG = ROOT.TF1("myfitMG","pol1",np.min(voltage),np.max(voltage))
        mygrMG.Fit("myfitMG","QR")
        MGpar0 = myfitMG.GetParameter(0)
        MGpar1 = myfitMG.GetParameter(1)
        MGerr0 = myfitMG.GetParError(0)
        MGerr1 = myfitMG.GetParError(1)
        chiMG = myfitMG.GetChisquare()
        dofMG = myfitMG.GetNDF()
        mygrFT = ROOT.TGraphErrors(len(voltage),voltage.flatten(),GainFFT.flatten(),np.zeros(len(volt),dtype=float),errGainFFT.flatten())
        c3 = ROOT.TCanvas("FFT_fit","FFT_Fit")
        c3.cd()
        mygrFT.Draw("AP")
        myfitFT = ROOT.TF1("myfitFT","pol1",np.min(voltage),np.max(voltage))
        mygrFT.Fit("myfitFT","QR")
        FTpar0 = myfitFT.GetParameter(0)
        FTpar1 = myfitFT.GetParameter(1)
        FTerr0 = myfitFT.GetParError(0)
        FTerr1 = myfitFT.GetParError(1)
        chiFT = myfitMG.GetChisquare()
        dofFT = myfitMG.GetNDF()
        plt.errorbar(voltage,GainMG,yerr=errGainMG,fmt='b.')
        plt.errorbar(voltage,GainFFT,yerr=errGainFFT,fmt='g.')
        vbdMG = -1.*MGpar0/MGpar1
        errvbdMG = np.sqrt(((MGerr0**2)/(MGpar0**2)+(MGerr1**2)/(MGpar1**2))*(vbdMG**2))
        vbdFT = -1.*FTpar0/FTpar1
        errvbdFT = np.sqrt(((FTerr0**2)/(FTpar0**2)+(FTerr1**2)/(FTpar1**2))*(vbdFT**2))
        normgMG = MGpar1*20e-15/50./1.6e-19
        errnormgMG = MGerr1*20e-15/50./1.6e-19
        normgFT = FTpar1*20e-15/50./1.6e-19
        errnormgFT = FTerr1*20e-15/50./1.6e-19
        plt.plot(voltage,linearfunc(voltage,MGpar0,MGpar1),'r--')
        plt.plot(voltage,linearfunc(voltage,FTpar0,FTpar1),'r-.')
        plt.grid(True)
        plt.xlabel('Bias Voltage [V]',fontsize=16)
        plt.ylabel('Gain [adu.]',fontsize=16)
        plt.annotate("MGF\n$U_{bd}$ : "+"{0:.2f} $\pm$ {1:.2f} V\nGain: {2:.2e} $e_0$/V\nTemp: {3:.1f}$^\circ$C\n$\chi^2$/DOF : {4:.1f}/{5}".format(vbdMG,errvbdMG,normgMG,np.mean(temp),chiMG,dofMG),xy=(0.35,0.75),xycoords='axes fraction',fontsize=14)
        plt.annotate("FFTF\n$U_{bd}$ : "+"{0:.2f} $\pm$ {1:.2f} V\nGain: {2:.2e} $e_0$/V\nTemp: {3:.1f}$^\circ$C\n$\chi^2$/DOF : {4:.1f}/{5}".format(vbdFT,errvbdFT,normgFT,np.mean(temp),chiFT,dofFT),xy=(0.65,0.75),xycoords='axes fraction',fontsize=14)
        plt.xlim(np.min(volt)-0.1,np.max(volt)+0.1)
    plt.show()

if __name__=='__main__':
    folderpath = sys.argv[1]
    
    GainvsVolt(folderpath)
