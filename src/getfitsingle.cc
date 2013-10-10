#include <iomanip>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TRint.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TPad.h>
#include "DaqMul.h"
#include "ExtractGain.h"

using namespace std;

//do multiple gaussian(MG) fit and FFT on spes histograms to calculate SiPM gain in qdcch
//calculate dcr from spes-<i>-dcr.root file
//log results in <dir>/spes.log
//magic number in GetDCR -->effective gate
Float_t dcreffgate = 100e-9;
void DivPad(Double_t l=0.1, Double_t r=0.1);

Int_t getfit(TString inputfilename, Int_t npk=5, Float_t kped=-1, Float_t kgn=-1){
    TString fullfilename = inputfilename;
    TFile *fin = new TFile(fullfilename,"read");
    TString fname= "";
    for(int i=fullfilename.Last('/')+1;i<fullfilename.Length();i++){
	fname+=fullfilename[i];
    }
    if(fin->IsZombie()) {
	fin->Close();
	return 1;
    }
    TTree *tree = (TTree*)fin->Get("data");
    DaqMul *spes = new DaqMul(tree);

    Float_t *mean = new Float_t[20];
    Float_t *sigma = new Float_t[20];
    Float_t *errmean = new Float_t[20];
    Float_t *errsigma = new Float_t[20];
    
    TH1I* hist = new TH1I("spes","spes in qdc channel;Charge [QDC channel];Entries",4096,0.5,4096.5);
    spes->GetHistogram(hist);
    Float_t *cond = new Float_t[3];
    spes->GetCondition(cond); //Temperature,errTemp,Voltage

    TRint rint("App",0x0,0);

//AutoCorrelation fit
    Float_t GainAuto = 0;
    GainAuto = GainFromAutoCor(hist);
    cout<<"Gain From Auto-correlation: "<<GainAuto<<endl;
    
//Multi gaussian fit
    TCanvas *can_ori = new TCanvas("spes","spes",800,500);
    TH1I *histmg = (TH1I*)hist->Clone("spes_mg");
    int retcnt = 0;
    if(kped<0){
	cout<<"Fit with peak finder"<<endl;
	retcnt = FitGainLog(histmg,npk,mean,sigma,errmean,errsigma);
    }
    else{
	Float_t preknowgain = (kgn<0)?GainAuto:kgn;
	retcnt = FitGainWithKnowledge(histmg,npk,mean,sigma,kped,preknowgain,errmean,errsigma);
    }
    const Int_t cnt =retcnt;
    Float_t LFgain,errLFgain,el_noise,pix_noise;
    TGraphErrors* gr = NULL;
    if(cnt>1){
	cout<<"Mean[0]: "<<mean[0]<<endl;
	for(int i=1;i<cnt;i++){
	    cout<<"Mean["<<i<<"]: +"<<mean[i]-mean[i-1]<<endl;
	}
    ///
	Float_t *pkorder = new Float_t[20];
	for(int i=0;i<cnt;i++){
	    pkorder[i] = i;
	}
	TCanvas* canGain = new TCanvas("grgain","grgain",800,500);
	gr = new TGraphErrors(cnt,pkorder,mean,0,errmean);
	gr->Draw("ap");
	TF1* tmpfit = new TF1("tmpfit","pol1",0-0.1,cnt-0.9);
	gr->Fit(tmpfit,"QR");
	//tmpfit->GetParameter(0);
	LFgain = tmpfit->GetParameter(1);
	errLFgain = tmpfit->GetParError(1);
	cout<<"Gain from MG: "<<LFgain<<" +/- "<<errLFgain<<endl;
	
	Float_t *sigmasq = new Float_t[20];
	Float_t *errsigsq = new Float_t[20];
	for(int i=0;i<cnt;i++){
	    sigmasq[i] = sigma[i]*sigma[i];
	    errsigsq[i] = TMath::Sqrt(2)*sigma[i]*errsigma[i];
	}
	TCanvas* canSigma = new TCanvas("grsigma","grsigma",800,500);
	TGraphErrors *grsigma = new TGraphErrors(cnt,pkorder,sigmasq,0,errsigsq);
	grsigma->Draw("ap");
	TF1* tmpfitsigma = new TF1("tmpfitsigma","pol1",0-0.1,cnt-0.9);
	grsigma->Fit(tmpfitsigma,"QR");
	el_noise = TMath::Sqrt(tmpfitsigma->GetParameter(0));
	pix_noise = TMath::Sqrt(tmpfitsigma->GetParameter(1));
	cout<<"Electronic noise: "<<el_noise<<" Pixel noise: "<<pix_noise<<endl;
    }
    else{
	cout<<fname<<" MG fit failed"<<endl;
    }
    ///
    TCanvas *canHistmg = new TCanvas("SPES","SPES",800,768);
    DivPad();
    TPad* subpadhist = (TPad*)canHistmg->cd(2);
    subpadhist->SetTickx(0);
    hist->SetLineColor(kRed);
    hist->SetFillStyle(3005);
    hist->SetFillColor(kRed);
    hist->Rebin(4);
    hist->Draw();
    Float_t yend = hist->GetBinContent(hist->GetMaximumBin())+200.;
    Float_t axisstart = int(mean[0]-LFgain*0.5);
    Float_t axisend = int(mean[cnt-1]+LFgain*0.5)+1;
    TAxis* oraxisX = hist->GetXaxis();
    oraxisX->SetTitle("Charge [a.u.]");
    oraxisX->SetRangeUser(axisstart,axisend);
    oraxisX->SetLabelFont(42);
    oraxisX->SetTitleFont(42);
    oraxisX->CenterTitle(1);
    oraxisX->SetTitleSize(0.05);
    oraxisX->SetTitleOffset(0.9);
    TAxis* oraxisY = hist->GetYaxis();
    oraxisY->SetRangeUser(0,yend);
    oraxisY->SetTitle("Entries");
    oraxisY->SetLabelFont(42);
    oraxisY->SetTitleFont(42);
    oraxisY->CenterTitle(1);
    oraxisY->SetTitleSize(0.05);
    //oraxisY->SetTitleOffset(0.67);
    Float_t peaxisstart = (axisstart-mean[0])/LFgain;
    Float_t peaxisend   = (axisend-mean[0])/LFgain;
    /* TGaxis *topax = new TGaxis(axisstart,yend,axisend,yend,peaxisstart,peaxisend,11,"-");
    topax->SetLabelSize(0);
    topax->SetTitleSize(0);
    topax->Draw();
    */
    TPad* subpadfit = (TPad*)canHistmg->cd(1);
    TH2F* histfr = new TH2F("frame",";N_{pe};#mu/",
			    oraxisX->GetNbins(),peaxisstart,peaxisend,
			    int(axisend-axisstart),axisstart,axisend);
    histfr->Draw();
    gr->SetMarkerStyle(23);
    gr->SetMarkerSize(1.3);
    gr->Draw("p");
    TAxis* fitaxisX = histfr->GetXaxis();
    fitaxisX->SetRangeUser(peaxisstart,peaxisend);
    fitaxisX->SetNdivisions(11);
    TAxis* fitaxisY = histfr->GetYaxis();
    fitaxisY->SetLabelSize(0);
    fitaxisY->SetTickLength(0);
    TGaxis *fitfrleftax = new TGaxis(fitaxisX->GetXmin(),fitaxisY->GetXmin(),fitaxisX->GetXmin(),fitaxisY->GetXmax(),fitaxisY->GetXmin(),fitaxisY->GetXmax(),5,"-");
    fitfrleftax->SetLabelFont(42);
    fitfrleftax->SetTitleFont(42);
    fitfrleftax->SetLabelSize(0.04/0.3*0.7);
    fitfrleftax->SetTitleSize(0.05*(0.7/0.3));
    fitfrleftax->SetTitleOffset(1./0.7*0.3);
    fitfrleftax->SetTitle("Charge [a.u.]");
    fitfrleftax->CenterTitle();
    fitfrleftax->Draw();
    TGaxis *fitfrtopax = new TGaxis(fitaxisX->GetXmin(),fitaxisY->GetXmax(),fitaxisX->GetXmax(),fitaxisY->GetXmax(),fitaxisX->GetXmin(),fitaxisX->GetXmax(),cnt,"-");
    fitfrtopax->SetLabelFont(42);
    fitfrtopax->SetTitleFont(42);
    fitfrtopax->SetLabelSize(0.04/0.3*0.7);
    fitfrtopax->SetLabelOffset(0.04);
    fitfrtopax->SetTitleSize(0.045*(0.7/0.3));
    fitfrtopax->SetTitleOffset(1);
    fitfrtopax->SetTitle("Number of Pixels");
    fitfrtopax->CenterTitle();
    fitfrtopax->Draw();

//FFT fit <----this is history now
    /*Float_t *GainFFT = NULL;
      GainFFT = GainFromFFT(hist);
      GainFFT = GainFromAutoCor(hist);
      cout<<"Gain from FFT: "<<GainFFT[0]<<" +/- "<<GainFFT[1]<<endl;

      Float_t *GainFFTnoPed = NULL;
      GainFFTnoPed = new Float_t[2];
      Float_t fftshiftped = 0.;
      if(cnt>1) fftshiftped = mean[0];
      else fftshiftped = hist->GetBinCenter(hist->GetMaximumBin());
      GainFFTnoPed = GainFromFFTShifted(hist,fftshiftped);  ///FFT fit after shift 0pe to 0
      cout<<"Gain from FFT Shifted: "<<GainFFTnoPed[0]<<" +/- "<<GainFFTnoPed[1]<<endl;
    */
    
//DCR
    //remove .root from name
    TString dcrfname = fullfilename.Remove(fullfilename.Sizeof()-6);
    dcrfname+="-dcr.root";
    TFile *fin2 = new TFile(dcrfname,"read");
    Float_t *DCRret = NULL;
    Bool_t havedcr=0;
    if(!fin2->IsZombie()){
	havedcr=1;
	TTree *tree2 = (TTree*)fin2->Get("data");
	DaqMul *spes2 = new DaqMul(tree2);
	TH1I* histdcr = new TH1I("spes-dcr","dcr spes in qdc channel",4096,-0.5,4095.5);
	TCanvas *candcr = new TCanvas("DCR","DCR",800,500);
	spes2->GetHistogram(histdcr);
	DCRret = GetDCR(histdcr,mean[0],LFgain,dcreffgate); //magic number effective gate
	cout<<"DCR : "<<DCRret[0]<<" +/- "<<DCRret[1]<<endl;
	cout<<"XT  : "<<DCRret[2]<<" +/- "<<DCRret[3]<<endl;
	candcr->SetLogy();
	histdcr->Draw();
	TAxis *dcrxax = histdcr->GetXaxis();
	dcrxax->SetRange(histdcr->FindFirstBinAbove(),histdcr->FindLastBinAbove());
	
    }

    if(cnt>1){
	//Float_t pxnoise = TMath::Sqrt(sigma[1]*sigma[1]-sigma[0]*sigma[0]);
	cout<<"<<<<<<<<<< entry be replaced in spes.log"<<endl;
	cout<<setprecision(2)<<std::fixed;
	cout<<cond[2]<<"\t";                       //voltage
	cout<<setprecision(1)<<std::fixed;
	cout<<cond[0]<<"\t"                       //temperature
	    <<mean[0]<<"\t"                       //pedestal
	    <<el_noise<<"\t"                      //noise
	    <<LFgain<<"\t"     //GainMG
	    <<pix_noise<<"\t"                       //pixel noise
	    <<errLFgain<<"\t"      //GainMG error
	    <<GainAuto<<"\t"                    //GainAuto
	    <<"1"<<"\t";                    //GainAuto error
	if(havedcr){
	    cout<<setprecision(0)<<std::fixed;
	    cout<<DCRret[0]<<"\t"   //DCR
		<<DCRret[1]<<"\t";   //DCR error
	    cout<<setprecision(3)<<std::fixed;
	    cout<<DCRret[2]<<"\t"   //XT
		<<DCRret[3]<<"\t";   //XT error
	}
	else{
	    cout<<setprecision(0)<<std::fixed;
	    cout<<"0\t0\t";
	    cout<<setprecision(3)<<std::fixed;
	    cout<<"0\t0\t";
	}
	cout<<fname.Data()<<endl;
	cout<<">>>>>>>>>>"<<endl;
    }
    else{
	cout<<"MG fit failed, try to give predefined parameters!"<<endl;
    }
    rint.Run(kTRUE);
    return 0;
}

Int_t main(int argc, char** argv){
    if(argc<2){
	cout<<"usage: ./getfitsingle <pathtofile> [ped] [gain] [npks]"<<endl;
	cout<<"[ped],[gain]: predefine 0ped peak and gain for fit"<<endl;
	cout<<"[npks]: tells how many peaks are expected"<<endl;
	exit(0);
    }
    TString fullfilename = argv[1];
    Float_t kped = -1;
    Float_t kgn = -1;
    Float_t npks = 5;
    if(argc>2){
	kped = atof(argv[2]);
    }
    if(argc>3){
	kgn = atof(argv[3]);
    }
    if(argc>4){
	npks = atoi(argv[4]);
    }
    TString tmpname = fullfilename;
    TString dcrfname = tmpname.Remove(fullfilename.Sizeof()-6);
    dcrfname+="-dcr.root";
    TFile *f = new TFile(dcrfname,"read");
    if(!(f->IsZombie())){
	cout<<"DCR file found, effective gate: "<<endl;
	cin>>dcreffgate;
    }
    f->Close();
    delete f;
    getfit(fullfilename,npks,kped,kgn);
}

void DivPad(Double_t l, Double_t r){

  TString name = gPad->GetName();
  TPad* pad = NULL;

  TString subname = name;
  subname+="_";
  subname+=1;
  pad = new TPad(subname,subname,0,0.7,1,1,0);
  //pad->SetGrid();
  pad->SetLeftMargin(0.1);
  pad->SetRightMargin(0.05);
  pad->SetBottomMargin(0.01);
  pad->SetTopMargin(0.3);
  pad->SetNumber(1);
  pad->Draw();

  subname = name;
  subname+="_";
  subname+=2;
  pad = new TPad(subname,subname,0,0,1,0.7,0);
  pad->SetLeftMargin(0.1);
  pad->SetRightMargin(0.05);
  pad->SetTopMargin(0.);
  pad->SetBottomMargin(0.1);
  pad->SetNumber(2);
  pad->Draw();
}
