#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <TROOT.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
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
    TH1F* hist = new TH1F("test","test",100,-0.005,0.005);
    TMultiGraph *mulg = new TMultiGraph();
    TGraph* gr = NULL;
    Float_t *tpV = new Float_t[10000];
    Float_t *tpT = new Float_t[10000];
    for(Long64_t ientry = 0;ientry<nentries;ientry++){
	tree->GetEntry(ientry);
	if(channel==1) continue;
	Int_t skipflg = 0;
	Int_t cross = 0;
	Double_t sum = 0;
	for(Int_t pt=0;pt<nsamples;pt++){
	    // if (V[pt]>0) V[pt] = 0;
	    V[pt] = -1*V[pt];
	    t[pt] = t[pt]*1e9;
	    //selections DCR
	    if(t[pt]>23.1 && t[pt]<23.2){
		if(V[pt]<0.005){
		    skipflg=1;
		    break;
		}
	    }
	    if(t[pt]>24.2 && t[pt]<24.3){
		if(V[pt]>0.035){
		    skipflg=1;
		    break;
		}
	    }
	}
	if(skipflg) continue;
	else{
	    //sum=sum*(t[1]-t[0]);
	    hist->Fill(sum/nsamples);
	}
//	if(Ascfname!="no"){
//	    for(Int_t ipt=0;ipt<nsamples;ipt++){
//		fout<<V[ipt]<<"\t";
//	    }
//	    fout<<"\n";
//	}
	//gr = new TGraph(nsamples,t,V);
	//gr->SetLineColor(45+idx);
	//mulg->Add((TGraph*)gr->Clone(Form("wf_%d",idx)));
//	if(idx==1) break;
    }
    /*
    idx==0;
    for(Long64_t ientry = 0;ientry<nentries;ientry++){
	tree->GetEntry(ientry);
	if(channel==1) continue;
	Int_t skipflg = 0;
	Int_t cross = 0;
	for(Int_t pt=0;pt<nsamples;pt++){
	    // if (V[pt]>0) V[pt] = 0;
	    V[pt] = -1.*V[pt];
	    t[pt] = t[pt]*1e9;
	    //selections DCR
	    
	    if(t[pt]>23.1 && t[pt]<23.2){
		if(V[pt]<0.005){
		    skipflg=1;
		    break;
		}
	    }
	    if(t[pt]>24.2 && t[pt]<24.3){
		if(V[pt]>0.035){
		    skipflg=1;
		    break;
		}
	    }
	}
	if(skipflg) continue;
	for(Int_t pt=0;pt<nsamples;pt++){
	    V[pt] = V[pt]-tpV[pt];
	    if(t[pt]>23 && t[pt]<28){
		if(V[pt]>0.015){
		    skipflg=1;
		    break;
		}
	    }
	}
	if(skipflg) continue;
	idx++;
	if(Ascfname!="no"){
	    for(Int_t ipt=0;ipt<nsamples;ipt++){
		fout<<V[ipt]<<"\t";
	    }
	    fout<<"\n";
	}
	gr = new TGraph(nsamples,t,V);
	gr->SetLineColor(kRed);
	mulg->Add((TGraph*)gr->Clone(Form("wf_%d",idx)));
	//if(idx==55) break;
	}*/
    fout.close();
    cout<<idx<<endl;
//   mulg->Draw("al");
    hist->Draw();
    Float_t fitcenter = hist->GetBinCenter(hist->GetMaximumBin());
    Float_t binwidth = hist->GetBinWidth(1);
    TF1* func= new TF1("fit","gaus",fitcenter-10*binwidth,fitcenter+4*binwidth);
    hist->Fit(func,"IRQ");
    cout<<func->Integral(-0.005,0.005)/binwidth<<endl;
//    mulg->GetYaxis()->SetRangeUser(1e-5,1e-1);
    //can->SetLogy(1);
    rint.Run(kTRUE);
    return 0;
}
