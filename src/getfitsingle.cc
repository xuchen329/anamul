#include <iostream>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TRint.h>
#include "DaqMul.h"
#include "ExtractGain.h"

using namespace std;

//do multiple gaussian(MG) fit and FFT on spes histograms to calculate SiPM gain in qdcch
//calculate dcr from spes-<i>-dcr.root file
//log results in <dir>/spes.log
//magic number in GetDCR -->effective gate

Int_t getfit(TString dir, TString fname,Int_t npk=5, Float_t kped=-1, Float_t kgn=-1){ 
    TString filefullname = dir+fname;
    TFile *fin = new TFile(filefullname,"read");
    if(fin->IsZombie()) {
	fin->Close();
	return 1;
    }
    TTree *tree = (TTree*)fin->Get("data");
    DaqMul *spes = new DaqMul(tree);

    Float_t *mean = new Float_t[20];
    Float_t *sigma = new Float_t[20];
    
    TH1I* hist = new TH1I("spes","spes in qdc channel",4096,1,4096);
    spes->GetHistogram(hist);
    Float_t *cond = new Float_t[3];
    spes->GetCondition(cond);
    TRint rint("App",0x0,0);
    TCanvas *can_ori = new TCanvas("spes","spes",800,600);

//Multi gaussian fit    
    TH1I *histmg = (TH1I*)hist->Clone("spes_mg");
    int retcnt = 0;
    if(kped<0){
	cout<<"Fit with peak finder"<<endl;
	retcnt = FitGainLog(histmg,npk,mean,sigma);
    }
    else{
	retcnt = FitGainWithKnowledge(histmg,npk,mean,sigma,kped,kgn);
    }
    const Int_t cnt =retcnt;
    histmg->Draw();
    histmg->GetXaxis()->SetRangeUser(0,histmg->FindLastBinAbove());
    Float_t *GainMG = NULL;
    if(cnt>1){ 
	GainMG = new Float_t [cnt-1];
	cout<<"Mean[0]: "<<mean[0]<<endl;
	for(int i=1;i<cnt;i++){
	    cout<<"Mean["<<i<<"]: +"<<mean[i]-mean[i-1]<<endl;
	    GainMG[i-1] = mean[i]-mean[i-1];
	}
	cout<<"Gain from MG :"<<TMath::Mean(cnt-1,GainMG)<<" +/- "<<TMath::RMS(cnt-1,GainMG)<<endl;
    }
    else{
	cout<<fname<<" MG fit failed"<<endl;
    }

//FFT fit
    Float_t *GainFFT = NULL;
    GainFFT = GainFromFFT(hist);
    cout<<"Gain from FFT: "<<GainFFT[0]<<" +/- "<<GainFFT[1]<<endl;

    Float_t *GainFFTnoPed = NULL;
    GainFFTnoPed = GainFromFFTNoPed(hist,kped,kgn);
    cout<<"Gain from FFT no ped: "<<GainFFTnoPed[0]<<" +/- "<<GainFFTnoPed[1]<<endl;

//DCR
    //remove .root from name
    TString dcrfname = filefullname.Remove(filefullname.Sizeof()-6);
    dcrfname+="-dcr.root";
    TFile *fin2 = new TFile(dcrfname,"read");
    Float_t *DCRret = NULL;
    if(!fin2->IsZombie()){
	TTree *tree2 = (TTree*)fin2->Get("data");
	DaqMul *spes2 = new DaqMul(tree2);
	TH1I* histdcr = new TH1I("spes-dcr","dcr spes in qdc channel",4096,1,4096);
	spes2->GetHistogram(histdcr);
	DCRret = GetDCR(histdcr,mean[0],(mean[1]-mean[0]),100e-9); //magic number effective gate
	cout<<"DCR : "<<DCRret[0]<<" +/- "<<DCRret[1]<<endl;
	cout<<"XT  : "<<DCRret[2]<<" +/- "<<DCRret[3]<<endl;
	TCanvas *candcr = new TCanvas("DCR","DCR",800,600);
	histdcr->Draw();
    }
    rint.Run(kTRUE);
    return 0;
}

Int_t main(int argc, char** argv){
    if(argc<3 || argc==4){
	cout<<"usage: ./getfitsingle <dir> <filename> [ped] [gain] [npks]"<<endl;
	exit(0);
    }
    TString dirname = argv[1];
    if(dirname[dirname.Length()-1]!='/'){
	dirname.Append('/');
    }
    TString filename = argv[2];
    Float_t kped = -1;
    Float_t kgn = -1;
    Float_t npks = 5;
    if(argc>3){
	kped = atof(argv[3]);
	kgn = atof(argv[4]);
	if(argc>5){
	    npks = atoi(argv[5]);
	}
    }
    getfit(dirname,filename,npks,kped,kgn);
}
