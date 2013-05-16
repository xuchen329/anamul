import numpy
import scipy.interpolate as intertool 
import matplotlib.pyplot as plt
import scipy.optimize as sop
import scipy.stats as scist
#from ROOT import TF1,TCanvas,TFile,TAxis

def PixelGausF(x,*par):
    gain,sigma = par
    return numpy.fft.fft((1./(numpy.sqrt(2*numpy.pi)*sigma))*numpy.exp(-0.5*(x-gain)**2/sigma**2))

def DetectorXTF(x,*par):
    gain,sigma,xt = par
    return PixelGausF(x,gain,sigma)*(1-xt)/(1-xt*PixelGausF(x,gain,sigma))

def DetLightResponseF(x,*par):
    mu,gain,sigma,xt = par
    return numpy.exp(mu*(DetectorXTF(x,gain,sigma,xt)-1))

def PedGausF(x,*par):
    ped,sigma = par
    return numpy.fft.fft((1./(numpy.sqrt(2*numpy.pi)*sigma))*numpy.exp(-0.5*(x-ped)**2/sigma**2))

def TotalResponse(x,*par):
    nentry,ped,noise,mu,gain,sigma,xt = par
    return nentry*numpy.fft.ifft(PedGausF(x,ped,noise)*DetLightResponseF(x,mu,gain,sigma,xt))

def TotalResponseFFT(x,*par):
    ped,noise,mu,gain,sigma,xt = par
    return (PedGausF(x,ped,noise)*DetLightResponseF(x,mu,gain,sigma,xt))

def TotalResponseCall(x,*par):
    nentry,ped,noise,mu,gain,sigma,xt = par
    return (nentry*numpy.fft.ifft(PedGausF(x,ped,noise)*DetLightResponseF(x,mu,gain,sigma,xt))).astype(numpy.float64)
"""

f.SetParameters(300.,30.,1.5,200.,14.,0.2,50000)
c = TCanvas()
f.Draw()
"""
"""
myfile = TFile("spes.root","read")
myhist = myfile.Get("spes_4")
#myhist.Draw()
xaxis = myhist.GetXaxis()
lowbond = myhist.FindFirstBinAbove()
highbond = myhist.FindLastBinAbove()
xaxis.SetRangeUser(lowbond,highbond)
f = TF1('myfunc',DetectorResponseFunction(),600,1000,7)
f.SetParameter(0,6.50908e+02)
f.SetParameter(1,1.19566e+01)
f.SetParameter(6,80000)
myhist.Fit(f,'R')
mypar = f.GetParameters()
print mypar
"""

#input qdc channel range, entries of each bin, list of guessed parameters.
#parameters: nevents,pedestal,noise,mu,gain,sigma,crosstalk  
def SPESfit(qdcch,dataset,pguess):
    TotalResponseFit = lambda p, x: p[0]*numpy.fft.ifft(PedGausF(x,p[1],p[2])*DetLightResponseF(x,p[3],p[4],p[5],p[6]))
    errTotalResponseFit = lambda p,x,y,erry:((1.0/erry)*(TotalResponseFit(p,x)-y)).astype(numpy.float64)
    dataset_err = numpy.sqrt(dataset)
    for i in range(len(dataset_err)):
        if dataset_err[i]<1:
            dataset_err[i]=1.0
    par, cov_x, info, messg, success = sop.leastsq(errTotalResponseFit,pguess[:],args=(qdcch,dataset,dataset_err),full_output=1)
    if success<1 or success>4 :
        print "Fit did not converge ...: "+messg
        return 0
    else:
        s_sq = ((errTotalResponseFit(par,qdcch,dataset,dataset_err).astype(numpy.float64))**2).sum()/(len(dataset)-len(pguess))
        cov_x = cov_x*s_sq
        par_err=()
        for i in range(7):
            par_err+=(numpy.sqrt(cov_x[i][i]),)
        return par,TotalResponseFit(par,qdcch),par_err

if __name__=='__main__':
    datay = numpy.genfromtxt('tmp.data')
    rightrange = numpy.max(datay.nonzero())
    pedestal = 189.16
    noise = 15.37
    mu = 4
    gain = 145.258
    sigma = 16.02
    crosstalk = 0.2
    Nentry = 200000
    Parameter = [Nentry,pedestal,noise,mu,gain,sigma,crosstalk]
    x = numpy.linspace(1,4096,4096)
    par,fitresult,perr = SPESfit(x, datay,Parameter)
    dof = len(datay)-1-len(Parameter)
    chi,pval = scist.chisquare(datay,fitresult,dof)
    
    plt.plot(x,datay,'b.')
    plt.plot(x,fitresult,'r--')
    plt.legend(('Measurement data','Fit'))
    plt.grid(True)
    plt.xlabel('QDC channel [adu]')
    plt.ylabel('# of Entries')
    plt.annotate('$\chi ^2$/DOF: {6:.1f}/{7}\nGain: {0:.1f} $\pm$ {8:.1f}\nPixel noise: {1:.1f} $\pm$ {9:.1f}\nCrosstalk prop: {2:.2f} $\pm$ {10:.2f}\n$\mu$: {3:.2f} $\pm$ {11:.2f}\n0 pe.: {4:.1f} $\pm$ {12:.1f}\nnoise: {5:.1f} $\pm$ {13:.1f}'.format(par[4],par[5],par[6],par[3],par[1],par[2],chi,dof,perr[4],perr[5],perr[6],perr[3],perr[1],perr[2]),xy=(0.6,0.5),xycoords='axes fraction',bbox=dict(boxstyle='round',fc='0.8'))
    plt.xlim(0,rightrange)
    plt.savefig('spes.pdf',format='pdf')
    plt.show()
    



#curve = TotalResponseFit(Parameter[:],x)
#print type(curve)
#curve = TotalResponse(x,*Parameter)
#cur_mag = numpy.sqrt(curve.real**2+curve.imag**2)
#print cur_mag
#plt.plot(x,cur_mag)
#plt.plot(x,curve)
#plt.plot(x,datay)
#curve_fit(TotalResponseFit,x,datay)
"""
p1,success = sop.leastsq(errTotalResponseFit,Parameter[:], args=(x,datay))
plt.plot(x,datay,'r.')
plt.plot(x,TotalResponseFit(p1,x),'--')
plt.legend(('measurement','fit'))
plt.grid(True)
plt.xlim(0,1500)
plt.xlabel('QDC channel [adu]')
plt.ylabel('# of Entries')
plt.annotate('Gain: {0:.1f}\nSigma: {1:.1f}\n$\mu$: {2:.2f}\nPedestal: {3:.1f}\nNoise: {4:.1f}\nCrosstalk: {5:.1%}'.format(p1[4],p1[5],p1[3],p1[1],p1[2],p1[6]),xy=(1000,500),xycoords='data',bbox=dict(boxstyle="round", fc="0.8"))
"""

