#include <iostream>
#include <TStyle.h>
#include <TRint.h>
#include <cmath>
#include "ExtractGain.h"

using namespace std;

const Float_t SEARCHSIG  = 15;
const Float_t BACKGROUND = 0.4;
const Float_t FITRANGE   = 0.25;
const Float_t START = 50;

TCanvas* NewCanvas(TString cname,Int_t x,Int_t y)
{
    TCanvas * can = new TCanvas(cname,cname,x,y);
    gStyle->SetOptTitle(0);
    return can;
}

Int_t FitGainLog(TH1I* hist, Int_t npks, Float_t *mean, Float_t *sigma,Float_t *errmean,Float_t *errsigma,Int_t pedes)
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
    if (errmean==0) errmean = new Float_t[20];
    if (errsigma==0) errsigma = new Float_t[20];
    
    TH1F* histtmp = new TH1F("tmp","tmp",4096,-0.5,4095.5);
    for(int i=START;i<=4096;i++){
	Float_t bincontent = hist->GetBinContent(i);
	if(bincontent>0) histtmp->SetBinContent(i,TMath::Log10(bincontent));
    }
/*    NewCanvas("mytmp",800,600);
      histtmp->Draw();*/
    
    //find peaks
    TSpectrum *spec = new TSpectrum(npks);
    Int_t nccnt = spec->Search(histtmp,SEARCHSIG,"",BACKGROUND);
    Float_t *px = spec->GetPositionX();
    if(pedes>0){
	for(int i=0;i<nccnt;i++){
	    if(px[i]<pedes-100){
		px[i] = 4096;
		nccnt-=1;
	    }
	}
    }
    const Int_t cnt = nccnt;
    Int_t* idx = new Int_t[20];
    delete histtmp;
   
    TMath::Sort(cnt,px,idx,kFALSE);
    Float_t gn =0.0;
    if(cnt>1) gn=px[idx[1]]-px[idx[0]]; //more than 1 peak found
    //else return 0;
    else { //only fit for pedestal
//!!! TF1* fitonce = new TF1("once","gaus",px[idx[0]]-8,px[idx[0]]+8);
	TF1* fitonce = new TF1("once","gaus",px[idx[0]]-0.5*gn,px[idx[0]]+0.5*gn);
      hist->Fit(fitonce,"QR");
      Float_t tmpmean = fitonce->GetParameter(1);
      Float_t tmpsigm = fitonce->GetParameter(2);
      Float_t tmperrmean = fitonce->GetParError(1);
      Float_t tmperrsigm = fitonce->GetParError(2);
      mean[0]    = tmpmean;
      sigma[0]   = tmpsigm;
      errmean[0] = tmperrmean;
      errsigma[0] = tmperrsigm;
      return 1;
    }
    //start fits
    TF1* gfit = NULL;
    TF1* gfittmp = NULL;
    Int_t NPks = 0;
    for(int i=0;i<cnt;i++){
	TString str = "g";
	str+=i;
	if(NPks==0){ //fitting 0pe peak
//!!!	    gfittmp = new TF1(str,"gaus",0.8*px[idx[i]],1.2*px[idx[i]]);
	    gfittmp = new TF1(str,"gaus",px[idx[i]]-0.5*gn,px[idx[i]]+0.5*gn);
	    hist->Fit(gfittmp,"QR");
	    Float_t tmpmean = gfittmp->GetParameter(1);
	    Float_t tmpsigm = gfittmp->GetParameter(2);
//!!!	    if(tmpmean<px[idx[i]]*0.8 || tmpmean>px[idx[i]]*1.2){
	    if(tmpmean<px[idx[i]]-0.5*gn || tmpmean>px[idx[i]]+0.5*gn){
		//	rint.Run(kTRUE);
		delete gfittmp;
		return 0;
	    }
	    delete gfittmp;
	    Float_t fitleftrange = tmpmean-1.5*tmpsigm;
	    Float_t fitrightrange = tmpmean+1.5*tmpsigm;
	    gfit = new TF1(str,"gaus",fitleftrange,fitrightrange);
	    hist->Fit(gfit,"QR");
	    mean[i]      = gfit->GetParameter(1);
	    sigma[i]     = gfit->GetParameter(2);
	    errmean[0]   = gfit->GetParError(1);
	    errsigma[0]  = gfit->GetParError(2);
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
	    errmean[NPks]   = gfit->GetParError(1);
	    errsigma[NPks]  = gfit->GetParError(2);
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

Int_t FitGainLogHighNoise(TH1I* hist, Int_t npks, Float_t *mean, Float_t *sigma,Float_t *errmean,Float_t *errsigma)
{//return length of mean array
   
    if (errmean==0) errmean = new Float_t[20];
    if (errsigma==0) errsigma = new Float_t[20];
    
    TH1F* histtmp = new TH1F("tmp","tmp",4096,1,4096);
    for(int i=START;i<=4096;i++){
	Float_t bincontent = hist->GetBinContent(i);
	if(bincontent>0) histtmp->SetBinContent(i,TMath::Log10(bincontent));
    }
    
    //find peaks
    TSpectrum *spec = new TSpectrum(npks);
    const Int_t cnt = spec->Search(histtmp,SEARCHSIG,"",BACKGROUND);
    Float_t *px = spec->GetPositionX();
    Int_t* idx = new Int_t[20];
   
    TMath::Sort(cnt,px,idx,kFALSE);
    Float_t *valeyleft = new Float_t[cnt];
    Float_t *valeyright = new Float_t[cnt];
    cout<<"Starting New Functions ..."<<endl;
    FindValey(histtmp,cnt,px,idx,valeyleft,valeyright);
    for(int pci=0;pci<cnt;pci++){
	cout<<px[idx[pci]]<<"   "<<valeyleft[pci]<<"   "<<valeyright[pci]<<endl;
    }
    delete histtmp;
    Float_t gn =0.0;
    if(cnt>1) {
	gn=px[idx[1]]-px[idx[0]]; //more than 1 peak found
    }
    else { //only fit for pedestal
//!!! TF1* fitonce = new TF1("once","gaus",px[idx[0]]-8,px[idx[0]]+8);
	TF1* fitonce = new TF1("once","gaus",px[idx[0]]-0.5*gn,px[idx[0]]+0.5*gn);
	hist->Fit(fitonce,"QR");
	Float_t tmpmean = fitonce->GetParameter(1);
	Float_t tmpsigm = fitonce->GetParameter(2);
	Float_t tmperrmean = fitonce->GetParError(1);
	Float_t tmperrsigm = fitonce->GetParError(2);
	mean[0]    = tmpmean;
	sigma[0]   = tmpsigm;
	errmean[0] = tmperrmean;
	errsigma[0] = tmperrsigm;
	return 1;
    }
    //start fits
    TF1* gfit = NULL;
    Int_t NPks = 0;
    for(int i=0;i<cnt;i++){
	TString str = "g";
	str+=i;
	if(NPks==0){ //fitting 0pe peak
	    gfit = new TF1(str,"gaus",valeyleft[0],valeyright[0]);
	    hist->Fit(gfit,"QR");
	    Float_t tmpmean = gfit->GetParameter(1);
	    Float_t tmpsigm = gfit->GetParameter(2);
	    if(tmpmean<valeyleft[0] || tmpmean>valeyright[0]){
		delete gfit;
		return 0;
	    }
	    else{
		mean[i]      = gfit->GetParameter(1);
		sigma[i]     = gfit->GetParameter(2);
		errmean[0]   = gfit->GetParError(1);
		errsigma[0]  = gfit->GetParError(2);
		//normfctor[i] = GetHistNormFactor(hist,mean[i],sigma[i]);
		std::cout<<"Pedestal: "<<mean[i]<<std::endl;
		std::cout<<"Noise: "<<sigma[i]<<std::endl;
		NPks++;
	    }
	}
	else {
	    gfit = new TF1(str,"gaus",valeyleft[i],valeyright[i]);
	    hist->Fit(gfit,"QR+");
	    Float_t tmpmean = gfit->GetParameter(1);
	    if(tmpmean<valeyleft[i] || tmpmean>valeyright[i]) continue;
	    mean[NPks] = tmpmean;
	    sigma[NPks] = gfit->GetParameter(2);
	    errmean[NPks]   = gfit->GetParError(1);
	    errsigma[NPks]  = gfit->GetParError(2);
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
    TH1I* histext = NULL;
    Int_t nbinsx = hist->GetNbinsX();
    Float_t lowedge = hist->GetBinLowEdge(1);
    histext = new TH1I("forfft","forfft",nbinsx,lowedge,nbinsx+0.5);
    for(int i=1;i<=nbinsx;i++){
	histext->SetBinContent(i,hist->GetBinContent(i));
    }
    TH1* hm=NULL;
    TVirtualFFT::SetTransform(0);
    hm = histext->FFT(hm,"MAG");      //root FFT gives same result as fftw3 package
    TString str = hist->GetName();
    str+="_fft";
    TCanvas* FFTtmpCanvas = NewCanvas(str,800,600);
    hm->Draw();
    Int_t *fitresult = new Int_t[2];
    FindPeak(hm,fitresult);
    Float_t norm = (floor(fitresult[0]+0.5))*histext->GetBinWidth(1);
    if(fitresult[1]<1.5) fitresult[1]=1.5;
    Float_t normerr = (floor(fitresult[0]+fitresult[1]))*histext->GetBinWidth(1);
    Float_t gain = (histext->GetNbinsX())/norm;
    Float_t errgain = (histext->GetNbinsX())/normerr;
    Float_t err = TMath::Abs(errgain-gain);
    delete fitresult;
    Float_t *retv = new Float_t[2];
    retv[0] = gain;
    retv[1] = err;
    delete histext;
  // delete FFTtmpCanvas;
    return retv;
}

Float_t* GainFromFFTShifted(TH1I* hist, Float_t ped){
    Float_t tmpmax = 0;
    Int_t pedbin = hist->FindBin(ped);
    Int_t shcenter = 0;
    for(int i=pedbin-35;i<pedbin+10;i++){
	if((hist->GetBinContent(i))>tmpmax){
	    tmpmax = hist->GetBinContent(i);
	    shcenter = i;
	}
    }
    TH1I* histext = NULL;
    Int_t nbinsx = hist->GetNbinsX();
    Float_t lowedge = hist->GetBinLowEdge(1);
    histext = new TH1I("forfft","forfft",64*nbinsx,lowedge,64*nbinsx+0.5);
    Int_t cpidx = shcenter;
    for(int i=1;i<=(nbinsx-shcenter);i++){
	histext->SetBinContent(i,hist->GetBinContent(cpidx));
	cpidx++;
    }
    TH1* hm=NULL;
    TVirtualFFT::SetTransform(0);
    hm = histext->FFT(hm,"MAG");      //root FFT gives same result as fftw3 package
    TString str = hist->GetName();
    str+="_fft_shifted";
    TCanvas* FFTtmpCanvas = NewCanvas(str,800,600);
    hm->Draw();
    Int_t *fitresult = new Int_t[2];
    FindPeak(hm,fitresult);
    Float_t norm = (fitresult[0])*histext->GetBinWidth(1);
    Float_t normerr = (fitresult[0]+fitresult[1])*histext->GetBinWidth(1);
    Float_t gain = (histext->GetNbinsX())/norm;
    Float_t errgain = (histext->GetNbinsX())/normerr;
    Float_t err = TMath::Abs(errgain-gain);
    delete fitresult;
    Float_t *retv = new Float_t[2];
    retv[0] = gain;
    retv[1] = err;
    // delete histext;
    // delete FFTtmpCanvas;
    return retv;
}

Float_t* GainFromFFTNoPed(TH1I* hist,Float_t kped, Float_t kgain){
    TH1I* histnew = new TH1I("hist-no-ped","spes no ped",4096,-0.5,4095.5);
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
    Int_t *fitresult = new Int_t[2];
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

void FindPeak(TH1* hist,Int_t *res){
  int startbin = 3;
  int endbin = 0.05*hist->GetNbinsX();

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
  int tmpcnt=1;
  for(int i=maxbin-10;i<maxbin+10;i++){
      if(hist->GetBinContent(i)==hist->GetBinContent(maxbin)) tmpcnt++;
  }
  /* 
  TF1* f1 = new TF1("f1","gaus",maxbin-20,maxbin+20);
  f1->SetParLimits(1,maxbin-2,maxbin+2);
  hist->Fit(f1,"RQ");*/
  res[0] = maxbin;
  res[1] = tmpcnt;
//  res[0] = f1->GetParameter(1);
//  res[1] = f1->GetParError(1);
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

Float_t* GetDCRHighNoise(TH1I* hist,Float_t pedestal,Float_t gain,Float_t effgate){
    Float_t fitcenter = pedestal;
    if(pedestal<0){
	fitcenter = hist->GetBinCenter(hist->GetMaximumBin());
    }
    TF1* pedfit = new TF1("pedfit","gaus",fitcenter-0.2*gain,fitcenter+0.2*gain);
    hist->Fit(pedfit,"QR");
    fitcenter = pedfit->GetParameter(1);
    Float_t fitsigma = pedfit->GetParameter(2);

    Float_t integralstart = hist->GetBinLowEdge(hist->FindFirstBinAbove());
    if(integralstart<(fitcenter-3.*fitsigma)) integralstart = fitcenter-3.*fitsigma;
    Int_t nbinsx = hist->GetNbinsX();

    Double_t AllNevents = hist->Integral(hist->FindFirstBinAbove(),nbinsx);
    Double_t NoDcrNevents = pedfit->Integral(integralstart,fitcenter+3.*fitsigma)/hist->GetBinWidth(1);
    
    //Double_t XtNevents = hist->Integral(pedestalcenter+1.5*gain,nbinsx);
    //cout<<"All Events: "<<AllNevents<<" DcrNevents: "<<DcrNevents<<endl;
    Double_t prop = NoDcrNevents/AllNevents;
    Double_t errprop = TMath::Sqrt(AllNevents-NoDcrNevents)/AllNevents;
    Double_t dcr = (1./effgate)*TMath::Log(1./prop);
    Double_t errdcr = (1./effgate)*TMath::Log(1./(prop-errprop))-dcr;
    //Double_t xtprop = XtNevents/DcrNevents;
    //Double_t errxt = TMath::Sqrt((1/XtNevents+1/DcrNevents)*(xtprop*xtprop));
    cout<<"DCR: "<<dcr<<" +/-: "<<errdcr<<endl;
    Float_t *ret = new Float_t[4];
    ret[0] = dcr;
    ret[1] = errdcr;
    //ret[2] = xtprop;
    //ret[3] = errxt;
    ret[2] = 0.;
    ret[3] = 0.;
    //delete pedfit;
    return ret;
}

Int_t FitGainWithKnowledge(TH1I* hist,Int_t npks, Float_t *mean, Float_t *sigma, Float_t kped, Float_t kgn, Float_t *errmean, Float_t *errsigma)
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
    if (errmean==0) errmean = new Float_t[20];
    if (errsigma==0) errsigma = new Float_t[20];
    //start fits
    TF1* gfit= NULL;
    TF1* gfittmp= NULL;
    Int_t pkcnt = 0;

    //fit pedestal
    TString str = "g";
    str+=pkcnt;
//    gfittmp = new TF1(str,"gaus",0.9*kped,1.1*kped);
    gfittmp = new TF1(str,"gaus",kped-50,kped+10);
    hist->Fit(gfittmp,"QR");
    Float_t tmpmean = gfittmp->GetParameter(1);
    Float_t tmpsigm = gfittmp->GetParameter(2);
    if(tmpmean<kped-50 || tmpmean>kped+10){
	return 0;
    }
//    delete gfittmp;
    /*
    Float_t fitleftrange = tmpmean-1.*tmpsigm;
    Float_t fitrightrange = tmpmean+1.*tmpsigm;
    gfit = new TF1(str,"gaus",fitleftrange,fitrightrange);
    hist->Fit(gfit,"QR");*/
    mean[0] = gfittmp->GetParameter(1);
    sigma[0]= gfittmp->GetParameter(2);
    errmean[0] = gfittmp->GetParError(1);
    errsigma[0]= gfittmp->GetParError(2);
    std::cout<<"Pedestal: "<<mean[0]<<std::endl;
    std::cout<<"Noise: "<<sigma[0]<<std::endl;
    pkcnt++;
    
    //fit peaks
    for(int i=0;i<npks;i++){
	TString str = "g";
	str+=pkcnt;
//	Float_t fitlrange = mean[0]+pkcnt*kgn-FITRANGE*kgn;
//	Float_t fitrrange = mean[0]+pkcnt*kgn+FITRANGE*kgn;
	Float_t fitlrange = mean[0]+pkcnt*kgn-40;
	Float_t fitrrange = mean[0]+pkcnt*kgn+40;
	gfit = new TF1(str,"gaus",fitlrange,fitrrange);
	hist->Fit(gfit,"QR+");
	Float_t tmpval1 = gfit->GetParameter(1);
	if (tmpval1<fitlrange || tmpval1>fitrrange) break;
	mean[pkcnt] = tmpval1;
	sigma[pkcnt] = gfit->GetParameter(2);
	errmean[pkcnt] = gfit->GetParError(1);
	errsigma[pkcnt] = gfit->GetParError(2);
	pkcnt++;
    }

    Float_t *gain = new Float_t[20];
    for(int i=0;i<pkcnt-1;i++){
	gain[i] = mean[i+1]-mean[i];
    }
    
    //std::cout<<"Mean of gain (Rough): "<<TMath::Mean(pkcnt-1,gain)<<std::endl;
    //std::cout<<"Error:        "<<TMath::RMS(pkcnt-1,gain)<<std::endl;
    
   
    delete[] gain;
    return pkcnt;
}

void FindValey(TH1F* histor, Int_t cnt, Float_t *px, Int_t* idx, Float_t *left, Float_t *right)
{
    //TRint rint("App",0x0,0);
    TH1F* hist = (TH1F*)histor->Rebin(4);
    Int_t rcnt = 0;
    for(int i=0;i<cnt;i++){
	if(px[idx[i]]<100) continue;
	Int_t binpx = hist->FindBin(px[idx[i]]);
	left[rcnt] = px[idx[i]]-1;
	right[rcnt] = px[idx[i]]+1;
	for(int bin=0;bin<5;bin++){
	    Int_t qualityleft = 0;
	    if(hist->GetBinContent(binpx-bin-1)<hist->GetBinContent(binpx-bin) &&
	       hist->GetBinContent(binpx-bin-2)<hist->GetBinContent(binpx-bin-1) &&
	       hist->GetBinContent(binpx-bin-3)<hist->GetBinContent(binpx-bin-2) )
		left[rcnt] = hist->GetBinCenter(binpx-bin-3);
	}
	for(int bin=0;bin<5;bin++){
	    if(hist->GetBinContent(binpx+bin+1)<hist->GetBinContent(binpx+bin) &&
	       hist->GetBinContent(binpx+bin+2)<hist->GetBinContent(binpx+bin+1) &&
	       hist->GetBinContent(binpx+bin+3)<hist->GetBinContent(binpx+bin+2))
		right[rcnt] = hist->GetBinCenter(binpx+bin+3);
	}
	if(right[rcnt]-left[rcnt]<5) continue;
	else rcnt++;
    }
    //rint.Run(kTRUE);
}
	
