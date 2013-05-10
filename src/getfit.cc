#include <iostream>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include "DaqMul.h"
#include "ExtractGain.h"

using namespace std;

//do multiple gaussian(MG) fit and FFT on spes histograms to calculate SiPM gain in qdcch
//calculate dcr from spes-<i>-dcr.root file
//log results in <dir>/spes.log
//!!! magic number in GetDCR, effective gate

Int_t fitandlog(TString dir, TString fname,Bool_t log){ 
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
    

    //Multi gaussian fit
    const int cnt = FitGainLog(hist,5,mean,sigma);
    Float_t *GainMG = NULL;
    Bool_t mgfitfail = 0;
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
	mgfitfail=1;
    }
    
    //FFT fit
    Float_t *GainFFT = NULL;
    GainFFT = GainFromFFT(hist);
    cout<<"Gain from FFT: "<<GainFFT[0]<<" +/- "<<GainFFT[1]<<endl;

    //DCR
    TString dcrfname = filefullname.Remove(filefullname.Sizeof()-6); //remove .root from name
    dcrfname+="-dcr.root";
    TFile *fin2 = new TFile(dcrfname,"read");
    Float_t *DCRret = NULL;
    if(!fin2->IsZombie() && (!mgfitfail)){
	TTree *tree2 = (TTree*)fin2->Get("data");
	DaqMul *spes2 = new DaqMul(tree2);
	TH1I* histdcr = new TH1I("spes-dcr","dcr spes in qdc channel",4096,1,4096);
	spes2->GetHistogram(histdcr);
	DCRret = GetDCR(histdcr,mean[0],(mean[1]-mean[0]),100e-9);  //magic number effective gate
	delete histdcr;
    }
    
    
    if(log){
	TString logfilename = dir;
	logfilename+="spes.log";
	std::ofstream fout(logfilename.Data(),std::ofstream::app);
	if(fout.tellp()<5){
	    fout<<"#";
	    fout<<"voltage temperature pedestal noise GainMG Pxlnoise errGainMG GainFFT errGainFFT DCR errDCR XT errXT"<<endl;
	}
	fout<<cond[2]<<"\t"<<cond[0]<<"\t"; //!! voltage temperature
	if(cnt>1){
	    fout<<mean[0]<<"\t"<<sigma[0]<<"\t"; //!! pedestal noise
	    fout<<TMath::Mean(cnt-1,GainMG)<<"\t"<<sigma[1]<<"\t"<<TMath::RMS(cnt-1,GainMG)<<"\t"; //!! GainMG pixelnoise errGainMG
	}
	else{ //MG failed
	    fout<<"0\t0\t0\t0\t0\t";
	}
	fout<<GainFFT[0]<<"\t"<<GainFFT[1];//!! GainFFT errGainFFT
	if(DCRret==NULL) fout<<"\t0\t0\t0\t0\t"<<fname.Data()<<"\n";
	else {
	    fout<<"\t";
	    fout<<DCRret[0]<<"\t"<<DCRret[1]<<"\t";//!!DCR errDCR
	    fout<<DCRret[2]<<"\t"<<DCRret[3]<<"\t"<<fname.Data()<<"\n";//!!xt errxt
	}
	fout.close();
    }
    return 0;
    //hist->Draw();
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
    
    for(int i=0;i<nfile;i++){
	TString name = "spes-";
	name+=i;
	name+=".root";
	Int_t ret = fitandlog(dirname,name,1);
    }
}
