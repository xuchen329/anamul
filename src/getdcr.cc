#include <iostream>
#include <fstream>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1I.h>
#include "DaqMul.h"
#include "ExtractGain.h"

using namespace std;

//<Gain> = fitpar0 + <voltage>*fitpar1
Int_t fitandlog(TString dir, TString fname, Float_t fitpar0, Float_t fitpar1,Float_t effgate){ 
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
    
    TH1I* histdcr = new TH1I("spes-dcr","spes in qdc channel",4096,1,4096);
    spes->GetHistogram(histdcr);
    Float_t *cond = new Float_t[3];
    spes->GetCondition(cond);
    //calculate gain
    Float_t gain = fitpar0+cond[2]*fitpar1;

    Float_t *DCRret=NULL;
    DCRret = GetDCR(histdcr,-1.,gain,effgate);    
    TString logfilename = dir;
    logfilename+="spesdcr.log";
    std::ofstream fout(logfilename.Data(),std::ofstream::app);
    fout<<cond[2]<<"\t"<<cond[0]<<"\t"; //!! voltage temperature
    if(DCRret!=NULL){
      fout<<DCRret[0]<<"\t"<<DCRret[1]<<"\t";//!!DCR errDCR
      fout<<DCRret[2]<<"\t"<<DCRret[3]<<"\n";//!!xt errxt
    }
    else{
      fout<<"\t0\t0\t0\t0\t"<<"\n";
    }
    fout.close();
    return 0;
}

Int_t main(int argc, char** argv){
    if(argc<6){
	cout<<"usage: ./getfit <dir> <nfile> <fitpar0> <fitpar1> <effgate>"<<endl;
	exit(0);
    }
    TString dirname = argv[1];
    if(dirname[dirname.Length()-1]!='/'){
	dirname.Append('/');
    }
    Int_t nfile = atoi(argv[2]);
    Float_t fitpar0 = atof(argv[3]);
    Float_t fitpar1 = atof(argv[4]);
    Float_t effgate = atof(argv[5]);
    //initialize log file
    TString logfilename = dirname;
    logfilename+="spesdcr.log";
    std::ofstream fout(logfilename.Data(),std::ofstream::out);
    fout.close();
    
    for(int i=0;i<nfile;i++){
	TString name = "spes-";
	name+=i;
	name+="-dcr.root";
	fitandlog(dirname,name,fitpar0,fitpar1,effgate);
    }
}
