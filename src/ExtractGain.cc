#include <iostream>
#include <TStyle.h>
#include <TRint.h>
#include <cmath>
#include "ExtractGain.h"

using namespace std;

const Float_t SEARCHSIG  = 15;
const Float_t BACKGROUND = 0.3;
const Float_t FITRANGE   = 0.2;
const Float_t START = 100;

TCanvas* NewCanvas(TString cname,Int_t x,Int_t y)
{
    TCanvas * can = new TCanvas(cname,cname,x,y);
    gStyle->SetOptTitle(0);
    return can;
}

Int_t FitGainLog(TH1I* hist, Int_t npks, Float_t *mean, Float_t *sigma)
{//return length of mean array
    //draw histogram
    /*
    TString cname = hist->GetName();
    TString signame = cname+"_err";
    TCanvas *can = NewCanvas(cname,800,600);
    //can->Divide(2,1);
    //can->cd(1);
    */
    /*
    hist->Draw();
    hist->GetYaxis()->SetTitle("# of Entries");
    hist->GetXaxis()->SetTitle("QDC channel [adu]");
    hist->GetXaxis()->SetRangeUser(0,hist->FindLastBinAbove());
    */
    TH1F* histtmp = new TH1F("tmp","tmp",4096,1,4096);
    for(int i=START;i<=4096;i++){
	Float_t bincontent = hist->GetBinContent(i);
	if(bincontent>0) histtmp->SetBinContent(i,TMath::Log10(bincontent));
    }
/*    NewCanvas("mytmp",800,600);
      histtmp->Draw();*/
    
    //find peaks
    TSpectrum *spec = new TSpectrum(npks);
    const Int_t cnt = spec->Search(histtmp,SEARCHSIG,"",BACKGROUND);
    Float_t *px = spec->GetPositionX();
    Int_t* idx = new Int_t[20];
    delete histtmp;
   
    TMath::Sort(cnt,px,idx,kFALSE);
    Float_t gn =0.0;
    if(cnt>1) gn=px[idx[1]]-px[idx[0]]; //more than 1 peak found
    else return 0;
    //start fits
    TF1* gfit = NULL;
    TF1* gfittmp = NULL;
    Int_t NPks = 0;
    for(int i=0;i<cnt;i++){
	TString str = "g";
	str+=i;
	if(NPks==0){ //fitting 0pe peak
	    gfittmp = new TF1(str,"gaus",0.8*px[idx[i]],1.2*px[idx[i]]);
	    hist->Fit(gfittmp,"QR");
	    Float_t tmpmean = gfittmp->GetParameter(1);
	    Float_t tmpsigm = gfittmp->GetParameter(2);
	    if(tmpmean<px[idx[i]]*0.8 || tmpmean>px[idx[i]]*1.2){
		//	rint.Run(kTRUE);
		delete gfittmp;
		return 0;
	    }
	    delete gfittmp;
	    Float_t fitleftrange = tmpmean-3*tmpsigm;
	    Float_t fitrightrange = tmpmean+3*tmpsigm;
	    gfit = new TF1(str,"gaus",fitleftrange,fitrightrange);
	    hist->Fit(gfit,"QR");
	    mean[i]      = gfit->GetParameter(1);
	    sigma[i]     = gfit->GetParameter(2);
	    //normfctor[i] = GetHistNormFactor(hist,mean[i],sigma[i]);
	    std::cout<<"Pedestal: "<<mean[i]<<std::endl;
	    std::cout<<"Noise: "<<sigma[i]<<std::endl;
	    NPks++;
	}
	else {
	    Float_t fitleftrange = px[idx[i]]-FITRANGE*gn;
	    Float_t fitrightrange = px[idx[i]]+FITRANGE*gn;
	    gfit = new TF1(str,"gaus",fitleftrange,fitrightrange);
	    hist->Fit(gfit,"QR+");
	    Float_t tmpmean = gfit->GetParameter(1);
	    if(tmpmean<fitleftrange || tmpmean>fitrightrange) continue;
	    mean[NPks] = tmpmean;
	    sigma[NPks] = gfit->GetParameter(2);
	    NPks++;
	    //normfctor[i] = GetHistNormFactor(hist,mean[i],sigma[i]);
	}
    }

    Float_t *gain = new Float_t[20];
    for(int i=0;i<NPks-1;i++){
	gain[i] = mean[i+1]-mean[i];
    }
    
    //std::cout<<"Mean of gain: "<<TMath::Mean(NPks-1,gain)<<std::endl;
    //std::cout<<"Error:        "<<TMath::RMS(NPks-1,gain)<<std::endl;

    delete[] px;
    delete[] idx;
    delete[] gain;
    return NPks;
}

Float_t* GainFromFFT(TH1I* hist){
  TH1* hm=NULL;
  TVirtualFFT::SetTransform(0);
  hm = hist->FFT(hm,"MAG");      //root FFT gives same result as fftw3 package
  TString str = hist->GetName();
  str+="_fft";
  TCanvas* FFTtmpCanvas = NewCanvas(str,800,600);
  hm->Draw();
  Float_t *fitresult = new Float_t[2];
  FindPeak(hm,fitresult);
  Float_t norm = (fitresult[0]*hist->GetBinWidth(1));
  Float_t normerr = (fitresult[0]+fitresult[1])*hist->GetBinWidth(1);
  Float_t gain = hist->GetNbinsX()/norm;
  Float_t errgain = hist->GetNbinsX()/normerr;
  Float_t err = TMath::Abs(errgain-gain);
  delete fitresult;
  Float_t *retv = new Float_t[2];
  retv[0] = gain;
  retv[1] = err; 
  // delete FFTtmpCanvas;
  return retv;
}

Float_t* GainFromFFTNoPed(TH1I* hist,Float_t kped, Float_t kgain){
    TH1I* histnew = new TH1I("hist-no-ped","spes no ped",4096,1,4096);
    Int_t binstart = floor(kped+0.5*kgain);
    for(int i=binstart;i<4096;i++){
	histnew->SetBinContent(i,hist->GetBinContent(i));
    }
    TH1* hm=NULL;
    TVirtualFFT::SetTransform(0);
    hm = histnew->FFT(hm,"MAG");
    TString str = histnew->GetName();
    str+="_fft_noped";
    TCanvas* FFTtmpCanvas = NewCanvas(str,800,600);
    hm->Draw();
    Float_t *fitresult = new Float_t[2];
    FindPeak(hm,fitresult);
    Float_t norm = (fitresult[0]*histnew->GetBinWidth(1));
    Float_t normerr = (fitresult[0]+fitresult[1])*histnew->GetBinWidth(1);
    Float_t gain = histnew->GetNbinsX()/norm;
    Float_t errgain = histnew->GetNbinsX()/normerr;
    Float_t err = TMath::Abs(errgain-gain);
    delete fitresult;
    Float_t *retv = new Float_t[2];
    retv[0] = gain;
    retv[1] = err; 
    //delete FFTtmpCanvas;
    return retv;
}

void FindPeak(TH1* hist,Float_t *res){
  int startbin = 3;
  int endbin = 0.1*hist->GetNbinsX();

  for(int i=startbin;i<endbin;i++) {
     if(    hist->GetBinContent(i)<hist->GetBinContent(i+1)
         && hist->GetBinContent(i)<hist->GetBinContent(i+2) 
	 && hist->GetBinContent(i)<hist->GetBinContent(i+3) ) {
        startbin=i;
	break;
     }
  }
  hist->SetAxisRange(startbin,endbin);
  int maxbin = hist->GetMaximumBin();
  
  TF1* f1 = new TF1("f1","gaus",maxbin-4,maxbin+4);
  f1->SetParLimits(1,maxbin-2,maxbin+2);
  hist->Fit(f1,"RQ");
  res[0] = f1->GetParameter(1);
  res[1] = f1->GetParError(1);
}

Float_t* GetDCR(TH1I* hist,Float_t pedestal,Float_t gain,Float_t effgate){
    Float_t pedestalcenter = pedestal;
    if(pedestal<0){
	Float_t fitcenter = hist->GetBinCenter(hist->GetMaximumBin());
	TF1* pedfit = new TF1("pedfit","gaus",0.8*fitcenter,1.2*fitcenter);
	hist->Fit(pedfit,"QR");
	fitcenter = pedfit->GetParameter(1);
	Float_t fitsigma = pedfit->GetParameter(2);
	delete pedfit;
	pedfit = new TF1("pedfit","gaus",fitcenter-3.*fitsigma,fitcenter+3.*fitsigma);
	hist->Fit(pedfit,"QR");
	pedestalcenter = pedfit->GetParameter(1);
    }
    Int_t nbinsx = hist->GetNbinsX();
    Double_t AllNevents = hist->Integral(0,nbinsx);
    Double_t DcrNevents = hist->Integral(pedestalcenter+0.5*gain,nbinsx);
    Double_t XtNevents = hist->Integral(pedestalcenter+1.5*gain,nbinsx);
    //cout<<"All Events: "<<AllNevents<<" DcrNevents: "<<DcrNevents<<endl;
    Double_t prop = DcrNevents/AllNevents;
    Double_t errprop = TMath::Sqrt(DcrNevents)/AllNevents;
    Double_t dcr = (1./effgate)*TMath::Log(1./(1.-prop));
    Double_t errdcr = (1./effgate)*TMath::Log(1./(1.-(prop+errprop)))-dcr;
    Double_t xtprop = XtNevents/DcrNevents;
    Double_t errxt = TMath::Sqrt((1/XtNevents+1/DcrNevents)*(xtprop*xtprop));
    cout<<"DCR: "<<dcr<<" +/-: "<<errdcr<<endl;
    Float_t *ret = new Float_t[4];
    ret[0] = dcr;
    ret[1] = errdcr;
    ret[2] = xtprop;
    ret[3] = errxt;
    return ret;
}

Int_t FitGainWithKnowledge(TH1I* hist,Int_t npks, Float_t *mean, Float_t *sigma, Float_t kped, Float_t kgn)
{//return length of mean array

    //draw histogram
    /*
    TString cname = hist->GetName();
    TString signame = cname+"_err";
    TCanvas *can = NewCanvas(cname,800,600);

    hist->Draw();
    hist->GetYaxis()->SetTitle("# of Entries");
    hist->GetXaxis()->SetTitle("QDC channel [adu]");
    hist->GetXaxis()->SetRangeUser(0,hist->FindLastBinAbove());
    */
    //start fits
    TF1* gfit= NULL;
    TF1* gfittmp= NULL;
    Int_t pkcnt = 0;

    //fit pedestal
    TString str = "g";
    str+=pkcnt;
    gfittmp = new TF1(str,"gaus",0.9*kped,1.1*kped);
    hist->Fit(gfittmp,"QR");
    Float_t tmpmean = gfittmp->GetParameter(1);
    Float_t tmpsigm = gfittmp->GetParameter(2);
    if(tmpmean<kped*0.8 || tmpmean>kped*1.2){
	return 0;
    }
    delete gfittmp;
    Float_t fitleftrange = tmpmean-3*tmpsigm;
    Float_t fitrightrange = tmpmean+3*tmpsigm;
    gfit = new TF1(str,"gaus",fitleftrange,fitrightrange);
    hist->Fit(gfit,"QR");
    mean[0] = gfit->GetParameter(1);
    sigma[0]= gfit->GetParameter(2);
    std::cout<<"Pedestal: "<<mean[0]<<std::endl;
    std::cout<<"Noise: "<<sigma[0]<<std::endl;
    pkcnt++;
    
    //fit peaks
    for(int i=0;i<npks;i++){
	TString str = "g";
	str+=pkcnt;
	Float_t fitlrange = mean[0]+pkcnt*kgn-FITRANGE*kgn;
	Float_t fitrrange = mean[0]+pkcnt*kgn+FITRANGE*kgn;
	gfit = new TF1(str,"gaus",fitlrange,fitrrange);
	hist->Fit(gfit,"QR+");
	Float_t tmpval1 = gfit->GetParameter(1);
	if (tmpval1<fitlrange || tmpval1>fitrrange) break;
	mean[pkcnt] = tmpval1;
	sigma[pkcnt] = gfit->GetParameter(2);
	pkcnt++;
    }

    Float_t *gain = new Float_t[20];
    for(int i=0;i<pkcnt-1;i++){
	gain[i] = mean[i+1]-mean[i];
    }
    
    std::cout<<"Mean of gain: "<<TMath::Mean(pkcnt-1,gain)<<std::endl;
    std::cout<<"Error:        "<<TMath::RMS(pkcnt-1,gain)<<std::endl;
    
   
    delete[] gain;
    return pkcnt;
}
