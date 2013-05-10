#include <iomanip>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include "DaqMul.h"
#include "ExtractGain.h"

using namespace std;

//do multiple gaussian(MG) fit and FFT on spes histograms, get gain in qdcch
//calculate dcr from spes-<i>-dcr.root file
//log results in <dir>/spes.log
//!!! magic number in GetDCR, effective gate
const Float_t dcreffgate = 100e-9;

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
    spes->GetCondition(cond);  //Temperature,errTemp,Voltage

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
	cout<<"Gain from MG :"<<TMath::Mean(cnt-1,GainMG)
	    <<" +/- "<<TMath::RMS(cnt-1,GainMG)<<endl;
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
    TString dcrfname = filefullname.Remove(filefullname.Sizeof()-6);//remove .root from name
    dcrfname+="-dcr.root";
    TFile *fin2 = new TFile(dcrfname,"read");
    Float_t *DCRret = NULL;
    if(!fin2->IsZombie() && (!mgfitfail)){
	TTree *tree2 = (TTree*)fin2->Get("data");
	DaqMul *spes2 = new DaqMul(tree2);
	TH1I* histdcr = new TH1I("spes-dcr","dcr spes in qdc channel",4096,1,4096);
	spes2->GetHistogram(histdcr);
	DCRret = GetDCR(histdcr,mean[0],(mean[1]-mean[0]),dcreffgate);//magic number effective gate
	delete histdcr;
    }

//record everthing in log file
    if(log){
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
	    fout<<mean[0]<<"\t"   //pedestal 
		<<sigma[0]<<"\t"; //noise
	    Float_t pxnoise = TMath::Sqrt(sigma[1]*sigma[1]-sigma[0]*sigma[0]);
	    fout<<TMath::Mean(cnt-1,GainMG)<<"\t" //GainMG 
		<<pxnoise<<"\t"                   //pixel noise
		<<TMath::RMS(cnt-1,GainMG)<<"\t"; //errGainMG
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
