#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;

int main(int argc,char** argv){
    if(argc<2){
	cout<<"Flat a root file to ascii"<<endl;
	cout<<"usage: ./macroFlatRooFile <rootfile> <outputfile>"<<endl;
	exit(0);
    }
    TString ROOfname = argv[1];
    TString Ascfname = argv[2];
    TFile *fin = new TFile(ROOfname,"read");
    if(fin->IsZombie()){
	cout<<"File not exist!"<<endl;
	exit(0);
    }
//    TTree* tree = (TTree*)fin->Get("data");
    TTree* tree = (TTree*)fin->Get("DAQ");

    UInt_t event;
    UInt_t nsamples;
    UInt_t channel;
    Float_t V[10000];
    Float_t t[10000];
    tree->SetBranchAddress("event",&event);
    tree->SetBranchAddress("nsamples",&nsamples);
    tree->SetBranchAddress("channel",&channel);
    tree->SetBranchAddress("V",&V);
    tree->SetBranchAddress("t",&t);

    ofstream fout(Ascfname.Data(),ofstream::out);
    const Long64_t nentries = tree->GetEntries();
    tree->GetEntry(1);
    const Int_t loopnpt = nsamples;
    cout<<channel<<endl;
//    fout<<qdcch<<"\t"<<temperature<<"\t"<<voltage<<endl;
/*    for(Long64_t ientry=1;ientry<nentries;ientry++){
	tree->GetEntry(ientry);
	fout<<qdcch<<"\t"<<temperature<<endl;
	}*/
    for(Int_t ipt=0;ipt<loopnpt;ipt++){
	tree->GetEntry(1);
	if(channel==1){
	    fout<<t[ipt]<<"\t";
	    fout<<V[ipt]<<"\t";
	}
	else{
	    cout<<"Wrong order!!"<<endl;
	    return 0;
	}
	for(int ipulse=0;ipulse<nentries;ipulse++){
	    tree->GetEntry(ipulse);
	    if(channel==1) continue;
	    if(channel==0){
		fout<<V[ipt]<<"\t";
	    }
	}
	fout<<"\n";
	if(!(ipt%100)) cout<<100.*(1.*ipt/(1.*loopnpt))<<"% finished\r"<<flush;
    }
    fout.close();
    return 0;
}
