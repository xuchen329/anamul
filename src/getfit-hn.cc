#include <iomanip>
#include <iostream>
#include <fstream>
#include <math.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include "DaqMul.h"
#include "ExtractGain.h"

using namespace std;

//do multiple gaussian(MG) fit and FFT on spes histograms, get gain in qdcch
//calculate dcr from spes-<i>-dcr.root file
//log results in <dir>/spes.log
//!!! magic number in GetDCR, effective gate
//!!! magic number in FFT, pedestal position

//high noise data version
const Float_t dcreffgate = 100e-9;

Int_t fitandlog(TString dir, TString fname,Bool_t log, Int_t pedes){
    Bool_t notskip = 1;
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
    Float_t *errmean = new Float_t[20];
    Float_t *errsigma = new Float_t[20];
    
    TH1I* hist = new TH1I("spes","spes in qdc channel",4096,-0.5,4095.5);
//HIGH NOISE VERSION, cut histogram from 50
    spes->GetHistogram(hist,50);  
    Float_t *cond = new Float_t[3];
    spes->GetCondition(cond);  //Temperature,errTemp,Voltage

//Multi gaussian fit
    const int cnt = FitGainLog(hist,7,mean,sigma,errmean,errsigma,pedes);
    Bool_t mgfitfail = 0;
    Float_t LFgain,errLFgain,el_noise,pix_noise;
    if(cnt>1){ 
	///
	Float_t *pkorder = new Float_t[20];
	for(int i=0;i<cnt;i++){
	    pkorder[i] = i;
	}
	
	TGraphErrors *gr = new TGraphErrors(cnt,pkorder,mean,0,errmean);
	TF1* tmpfit = new TF1("tmpfit","pol1",0-0.1,cnt-0.9);
	gr->Fit(tmpfit,"QR");
	LFgain = tmpfit->GetParameter(1);
	errLFgain = tmpfit->GetParError(1);
	cout<<"Gain from MG: "<<LFgain<<" +/- "<<errLFgain<<endl;
	
	Float_t *sigmasq = new Float_t[20];
	Float_t *errsigsq = new Float_t[20];
	for(int i=0;i<cnt;i++){
	    sigmasq[i] = sigma[i]*sigma[i];
	    errsigsq[i] = TMath::Sqrt(2)*sigma[i]*errsigma[i];
	}
	TGraphErrors *grsigma = new TGraphErrors(cnt,pkorder,sigmasq,0,errsigsq);
	TF1* tmpfitsigma = new TF1("tmpfitsigma","pol1",0-0.1,cnt-0.9);
	grsigma->Fit(tmpfitsigma,"QR");
	el_noise = TMath::Sqrt(tmpfitsigma->GetParameter(0));
	pix_noise = TMath::Sqrt(tmpfitsigma->GetParameter(1));
	cout<<"Electronic noise: "<<el_noise<<" Pixel noise: "<<pix_noise<<endl;
	///
    }
    else{
	cout<<fname<<" MG fit failed"<<endl;
	mgfitfail=1;
    }

//FFT fit
    Float_t *GainFFT = NULL;
    GainFFT = GainFromFFTShifted(hist,pedes);
    cout<<"Gain from FFT: "<<GainFFT[0]<<" +/- "<<GainFFT[1]<<endl;
    if(GainFFT[1] > 10) notskip=0;

//DCR
    TString dcrfname = filefullname.Remove(filefullname.Sizeof()-6);//remove .root from name
    dcrfname+="-dcr.root";
    TFile *fin2 = new TFile(dcrfname,"read");
    Float_t *DCRret = NULL;
    if(!fin2->IsZombie() && (!mgfitfail)){
	TTree *tree2 = (TTree*)fin2->Get("data");
	DaqMul *spes2 = new DaqMul(tree2);
	TH1I* histdcr = new TH1I("spes-dcr","dcr spes in qdc channel",4096,-0.5,4095.5);
	spes2->GetHistogram(histdcr,50);
	DCRret = GetDCR(histdcr,mean[0],LFgain,dcreffgate);//magic number effective gate
	delete histdcr;
    }

//record everthing in log file
    if(log && notskip){
	TString logfilename = dir;
	logfilename+="spes.log";
	std::ofstream fout(logfilename.Data(),std::ofstream::app);
	if(fout.tellp()<5){
	    fout<<"#";
	    fout<<"Vop[V]\tT[C]\t0pe\tnoise\tGainMG\tPex\terrMG\tGainFFT\terrFFT\tDCR\terrDCR\tXT\terrXT"<<endl;
	}
	fout<<setprecision(2)<<std::fixed;
	fout<<cond[2]<<"\t";  //voltage
	fout<<setprecision(1)<<std::fixed;
	fout<<cond[0]<<"\t"; //temperature
	if(cnt>1){
	    fout<<mean[0]<<"\t"    //pedestal
		<<el_noise<<"\t";  //noise
	    fout<<LFgain<<"\t"     //GainMG 
		<<pix_noise<<"\t"  //pixel noise
		<<errLFgain<<"\t"; //errGainMG
	}
	else{ //MG failed
	    fout<<"0\t0\t0\t0\t0\t";
	}
	fout<<GainFFT[0]<<"\t"   //GainFFT
	    <<GainFFT[1]<<"\t";  //errGainFFT
	if(DCRret==NULL){        //noDCR
	    fout<<setprecision(0)<<std::fixed;
	    fout<<"0\t0\t";
	    fout<<setprecision(3)<<std::fixed;
	    fout<<"0\t0\t";
	}
	else {
	    fout<<setprecision(0)<<std::fixed;
	    fout<<DCRret[0]<<"\t"   //DCR
		<<DCRret[1]<<"\t";  //errDCR
	    fout<<setprecision(3)<<std::fixed;
	    fout<<DCRret[2]<<"\t"  //XT
		<<DCRret[3]<<"\t"; //errXT
	}
	fout<<fname.Data()<<"\n";
	fout.close();
    }
    return 0;
}

Int_t GetPedestal(TString dirname){
    TString filefullname = dirname;
    filefullname+="pedestal.root";
    TFile* fin = new TFile(filefullname,"read");
    if(fin->IsZombie()) {
	fin->Close();
	return 0.;
    }
    TTree *tree = (TTree*)fin->Get("data");
    DaqMul *spes = new DaqMul(tree);
    TH1I* hist = new TH1I("spes","spes in qdc channel",4096,-0.5,4095.5);
//HIGH NOISE VERSION, cut histogram from 10
    spes->GetHistogram(hist,10);
    Int_t pedes = hist->GetBinCenter(hist->GetMaximumBin());
    cout<<"Pedestal Estimation: "<<pedes<<endl;
    delete hist;
    delete spes;
    return pedes;
}

Int_t main(int argc, char** argv){
    if(argc<3){
	cout<<"usage: ./getfit <dir> <nfile>"<<endl;
	exit(0);
    }
    TString dirname = argv[1];
    if(dirname[dirname.Length()-1]!='/'){
	dirname.Append('/');
    }
    Int_t nfile = atoi(argv[2]);

    //initialize log file
    TString logfilename = dirname;
    logfilename+="spes.log";
    std::ofstream fout(logfilename,std::ofstream::out);
    fout.close();

    Int_t pedes = GetPedestal(dirname);
    for(int i=0;i<nfile;i++){
	TString name = "spes-";
	name+=i;
	name+=".root";
	Int_t ret = fitandlog(dirname,name,1,pedes);
    }
}
