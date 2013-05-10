#include <stdio.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>

int main(int argc, char** argv){
    TFile *fout = new TFile("sipmsim.root");
    TTree* tree = new TTree("data","data");
    
    Int_t qdcch;
    Int_t charge[32];
    Int_t id = -1;
    Float_t gate_width;
    Float_t voltage=0;
    Float_t temperature=0;
    Float_t systime=0;
    Float_t sysdate=0;
    tree->Branch("qdcch",&qdcch,"qdcch/I");
    for(int i=0;i<32;i++){
	tree->Branch(Form("qdcch%d",i),&charge[i],Form("qdcch%d/I",i));
    }
    tree->Branch("temperature",&temperature,"temperature/F");
    tree->Branch("time",&systime,"time/F");
    tree->Branch("date",&sysdate,"date/F");
    tree->Branch("voltage",&voltage,"voltage/F");
    tree->Branch("gate_width",&gate_width,"gate_width/F");
    tree->Branch("id",&id,"id/I");

    
}

