#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import sys,os
from functionconv import *
from matplotlib.backends.backend_pdf import PdfPages
import scipy.stats as scist

"""
def GetParFromLog(voltage,fname='spes.log'):
    allpar = np.genfromtxt(fname)
    for entry in allpar:
        if numpy.abs((voltage-entry[0]))<0.0001:
            return entry
    return [0]
"""
def GetParFromLog(rootfile,logfile='spes.log'):
    allpar = np.genfromtxt(logfile,dtype=(float,float,float,float,float,float,float,float,float,float,float,float,float,'S20',))
    for entry in allpar:
        if rootfile in entry:
            return entry
    print rootfile+" not found!! return ..."
    return [0]

def GetHistogramData(fname):
    print "open file "+fname+" ..."
    chain = ROOT.TChain("data")
    chain.Add(fname)
    chain.SetBranchStatus("*",0)
    chain.SetBranchStatus("voltage",1)
    chain.GetEntry(0)
    thisvoltage = chain.voltage
    print "measurement voltage: {0:.1f}".format(thisvoltage)
    chain.SetBranchStatus("*",0)
    chain.SetBranchStatus("qdcch",1)
    nevents = chain.GetEntries()
    #GetEntryMth = chain.GetEntry #Get the method pointer
    hist = ROOT.TH1I("spes","spes",4096,1,4096)
    #FillMth = hist.Fill
    for i in range(nevents):
        chain.GetEntry(i)
        hist.Fill(chain.qdcch)
    ret=[]
    for i in xrange(1,4097):
        ret.append(hist.GetBinContent(i))
    thishistmean = hist.GetMean()
    return ret,thisvoltage,thishistmean,nevents

def WriteParToFile(filename,par,perr):
    headerlb = '#pedestal noise mu gain sigma crosstalk volt temp\n'
    savefmt = ['{:.1f}','{:.1f}','{:.2f}','{:.1f}','{:.1f}','{:.3f}','{:.2f}','{:.1f}']
    if not os.path.isfile(filename):
        f = open(filename,'w')
        f.write(headerlb)
    else:
        f = open(filename,'a')
    for i in range(len(perr[1:])):
        idx = i+1
        f.write((savefmt[i]).format(par[idx]))
        f.write('\t')
        f.write((savefmt[i]).format(perr[idx]))
        f.write('\t')
    f.write(savefmt[-2].format(par[-2]))
    f.write('\t')
    f.write(savefmt[-1].format(par[-1]))
    f.write('\n')
    f.close()


def DoFitReturnPlot(fitpars,filepath,pdffile=None):
    filename = filepath+fitpars[-1]
    datay,volt,tmphistmean,tmpnentry = GetHistogramData(filename)
    tmppar = fitpars
    if tmppar[2]==0:
        return
    #print tmppar
    tmpmu = (tmphistmean-tmppar[2])/tmppar[4] #estimate mu: (mean of hist - pedestal)/gain
    guess = []
    guess.append(tmpnentry)#Nentry
    guess.append(tmppar[2])#pedestal
    guess.append(tmppar[3])#noise
    guess.append(tmpmu)#mu
    guess.append(tmppar[4])#gainMG
    guess.append(tmppar[5])#pixelnoise
    guess.append(0.05)#xtalk
    x = np.linspace(1,4096,4096)
    par,fitresult,perr = SPESfit(x,datay,guess)
    rightrange = numpy.max(numpy.array(datay).nonzero())
    heightrange = numpy.max(numpy.array(datay))
    dof = len(datay)-1-len(guess)
    chi,pval = scist.chisquare(datay,fitresult)
    #schi,pval = scist.chisquare(datay,fitresult,dof) #when scipy>1.0 enable this
    myplot = plt.plot(x,datay,'b.')
    plt.plot(x,fitresult,'r--')
    plt.legend(('Measurement data','Fit'))
    plt.grid(True)
    plt.xlabel('QDC channel [adu]')
    plt.ylabel('# of Entries')
    plt.annotate('$\chi ^2$/DOF: {6:.1f}/{7}\nGain: {0:.1f} $\pm$ {8:.1f}\nPixel noise: {1:.1f} $\pm$ {9:.1f}\nCrosstalk prop: {2:.2f} $\pm$ {10:.2f}\n$\mu$: {3:.2f} $\pm$ {11:.2f}\n0 pe.: {4:.1f} $\pm$ {12:.1f}\nnoise: {5:.1f} $\pm$ {13:.1f}'.format(par[4],par[5],par[6],par[3],par[1],par[2],chi,dof,perr[4],perr[5],perr[6],perr[3],perr[1],perr[2]),xy=(0.6,0.5),xycoords='axes fraction',bbox=dict(boxstyle='round',fc='0.8'))
    plt.annotate('Bias voltage: {0:.2f}V\nTemperature: {1:.1f}$^\circ$C'.format(tmppar[0],tmppar[1]),xy=(0.6,0.2),xycoords='axes fraction',bbox=dict(boxstyle='round',fc='0.8'))
    plt.xlim(0,rightrange)
    plt.ylim(0,heightrange)
    par=np.append(par,tmppar[0])
    par=np.append(par,tmppar[1])
    WriteParToFile(filepath+'pyspes.log',par,perr)
    if pdffile is not None:
        plt.savefig(pdffile,format='pdf')
        plt.clf()
    else:
        plt.show()




if __name__=='__main__':
    dirpath = sys.argv[1]
    tmppath, tmpfolder, files = os.walk(dirpath).next()
    pp = PdfPages(tmppath+'/out.pdf')
    for fname in files:
        if '.root' in fname:
            retval = GetParFromLog(fname,tmppath+'/spes.log')
            if retval[0] is not 0:
                DoFitReturnPlot(retval,tmppath+'/',pp)
    pp.close()
          
          
            

"""
        for i in range(15):
            fname='spes-'+str(i)+'.root'
            if os.path.isfile(fname):
                DoFitReturnPlot(fname,pp)
        pp.close()
    elif '.root' in filename:
        DoFitReturnPlot(filename)
        plt.show()
""" 
        
    
"""
    datay,volt = GetHistogramData("spes-1.root")
    print volt
    x = np.linspace(1,4096,4096)
    plt.plot(x,datay)
    plt.show()
"""
