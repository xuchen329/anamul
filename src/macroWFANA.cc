#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TRint.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TCanvas.h>

using namespace std;

int main(int argc,char** argv){
    if(argc<2){
	cout<<"ana a pulse waveform"<<endl;
	cout<<"usage: ./macroWFANA <rootfile>"<<endl;
	exit(0);
    }
    TString ROOfname = argv[1];
    TString Ascfname = "no";
    if(argc>2) Ascfname = argv[2];
    TFile *fin = new TFile(ROOfname,"read");
    if(fin->IsZombie()){
	cout<<"File not exist!"<<endl;
	exit(0);
    }
    ofstream fout(Ascfname.Data(),ofstream::out);

    
//set up tree
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
    const Long64_t nentries = tree->GetEntries();
    TRint rint("App",0x0,0);

//looping file    
    Int_t idx = 0;
    /*
    tree->GetEntry(0);
    Float_t Time_Start = t[0];
    Float_t Time_End   = t[nsamples];
    Float_t Time_Nbins = int((Time_End-Time_Start)/(t[1]-t[0]));
    TH2F* histfr = new TH2F("fr","fr",Time_Nbins,Time_Start,Time_End,1e5,1e-6,1e-1);
    histfr->Draw();
    can->SetLogy(1);*/
    TCanvas* cangate = new TCanvas("gate","gate",800,600);
    tree->GetEntry(1);
    for(Int_t pt=0;pt<nsamples;pt++){
	t[pt]=t[pt]*1e9;
    }
    TGraph* grgate = new TGraph(nsamples,t,V);
    grgate->Draw("al");
    if(Ascfname!="no"){
	for(Int_t ipt=0;ipt<nsamples;ipt++){
	    fout<<t[ipt]<<"\t";
	}
	fout<<"\n";
	for(Int_t ipt=0;ipt<nsamples;ipt++){
	    fout<<V[ipt]<<"\t";
	}
	fout<<"\n";
    }
    
    TCanvas* can  = new TCanvas("fr","fr",800,600);
    TMultiGraph *mulg = new TMultiGraph();
    TGraph* gr = NULL;
    for(Long64_t ientry = 0;ientry<nentries;ientry++){
	tree->GetEntry(ientry);
	if(channel==1) continue;
	Int_t skipflg = 0;
	for(Int_t pt=0;pt<nsamples;pt++){
	    if (V[pt]>0) V[pt] = 0;
	    else V[pt] = -1.*V[pt];
	    t[pt] = t[pt]*1e9;
	    //selections
	    /* if(t[pt]>25.8 && t[pt]<26.2){
		if(V[pt]>0.028 || V[pt]<0.016){
		    skipflg=1;
		    break;
		}
	    }
	    if(t[pt]>39.8 && t[pt]<40.2){
		if(V[pt]>0.012 || V[pt]<0.0035){
		    skipflg=1;
		    break;
		}
	    }
	    if(t[pt]>60. && t[pt]<100.){
		if(V[pt]>0.005){
		    skipflg=1;
		    break;
		}
	    }*/
	    
	    if(t[pt]<48){
		if(V[pt]>0.003){
		    skipflg=1;
		    break;
		}
	    }
	   
	    /*
	    if(t[pt]>52 && t[pt]<54){
		if(V[pt]>0.0255 || V[pt]<0.009) {
		    skipflg=1;
		    break;
		}
	    }
	    if(t[pt]>69 && t[pt]<71){
		if(V[pt]>0.008) {
		    skipflg=1;
		    break;
		}
	    }
	    if(t[pt]>80 && t[pt]<120){
		if(V[pt]>0.0056) {
		    skipflg=1;
		    break;
		}
		}*/
	}
	if(skipflg) continue;
	if(Ascfname!="no"){
	    for(Int_t ipt=0;ipt<nsamples;ipt++){
		fout<<V[ipt]<<"\t";
	    }
	    fout<<"\n";
	}		
	gr = new TGraph(nsamples,t,V);
	gr->SetLineColor(50+idx);
	mulg->Add((TGraph*)gr->Clone(Form("wf_%d",idx)));
	idx++;
	if(idx==25) break;
    }
    fout.close();
    mulg->Draw("al");
    mulg->GetYaxis()->SetRangeUser(1e-5,1e-1);
    can->SetLogy(1);
    rint.Run(kTRUE);
    return 0;
}
