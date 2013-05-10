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
    TTree* tree = (TTree*)fin->Get("data");

    Int_t qdcch;
    Float_t voltage;
    Float_t temperature;
    tree->SetBranchAddress("qdcch",&qdcch);
    tree->SetBranchAddress("voltage",&voltage);
    tree->SetBranchAddress("temperature",&temperature);

    ofstream fout(Ascfname.Data(),ofstream::out);
    const Long64_t nentries = tree->GetEntries();
    tree->GetEntry(0);
    fout<<qdcch<<"\t"<<temperature<<"\t"<<voltage<<endl;
    for(Long64_t ientry=1;ientry<nentries;ientry++){
	tree->GetEntry(ientry);
	fout<<qdcch<<"\t"<<temperature<<endl;
    }
    fout.close();
    return 0;
}
