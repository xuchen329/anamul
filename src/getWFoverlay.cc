#include <iomanip>
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
#include <TString.h>

//#include "DAQ_ADU.h"

using namespace std;

//analyse waveform of one pixel firing for AF prop

int main(int argc,char** argv){
    if(argc<2){
	cout<<"analyse waveform of one pixel firing for AF prop"<<endl;
	cout<<"usage: ./getWF <rootfile> [dcr]"<<endl;
	exit(0);
    }
    TString ROOfname = argv[1];
    TFile *fin = new TFile(ROOfname,"read");
    if(fin->IsZombie()){
	cout<<"File not exist!"<<endl;
	exit(0);
    }
    
    const Int_t wfchannel = 0; //waveform channel

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

    Float_t tpV[10000] = {};
    Float_t tpT[10000] = {};
    Float_t sumV[10000] = {};
    Float_t sumVpos[10000] = {};
    Float_t sumT[10000] = {};
    Int_t length = 0;
    Float_t counter= 0;

    TMultiGraph *mg = new TMultiGraph();
    TGraph* tmpgr = NULL;
    
    for(Long64_t ientry = 0;ientry<nentries;ientry++){
	tree->GetEntry(ientry);
	if(channel==1) continue;
	Float_t delt = (t[1]-t[0])*1e9;
	Int_t delpt = int(2./delt);
	Int_t skipflg = 0;
	Float_t BL = 0;
	Float_t BLcounter =0;
	Float_t localmin =0;
	Float_t minpos =0;
	for(Int_t pt=0;pt<nsamples;pt++){
	    //tpV[pt] = V[pt];
	    tpT[pt] = t[pt]*1e9;
	    //select base line
	    if(tpT[pt]>0 && tpT[pt]<20){
		if(V[pt]<-0.002){
		    skipflg = 1;
		    break;
		}
		BL+=V[pt];
		BLcounter++;
	    }
	    if(tpT[pt]>=23.2 && tpT[pt]<=23.3){
		if(V[pt]>-0.005){
		    skipflg = 1;
		    break;
		}
	    }
	    if(tpT[pt]>=24 && tpT[pt]<=24.5) {
		if(V[pt]<localmin){
		    localmin = V[pt];
		    minpos = tpT[pt];
		}
	    }
	    if(localmin<-0.035){
		skipflg =1;
		break;
	    }
	    //select one pixel
	    if(tpT[pt]>30 && tpT[pt]<60.){
		Float_t p0 = (localmin*60.-minpos*(BL/BLcounter))/(60.-minpos);
		Float_t p1 = (BL/BLcounter-localmin)/(60.-minpos);
		Float_t thr = tpT[pt]*p1+p0;
		if (thr>0) break;
		if(V[pt]<thr){
		    //    cout<<V[pt]<<"  "<<thr<<endl;
		    skipflg = 1;
		    break;
		}
	    }
	    if(tpT[pt]>60) break;
	    
	}
	if (skipflg) continue;
	else{
	    tmpgr = new TGraph(nsamples,t,V);
	    mg->Add((TGraph*)tmpgr->Clone(Form("gr_%d",ientry)));
	    cout<<counter<<"\r"<<flush;
	    counter+=1.;
	    for(Int_t pt=0;pt<nsamples;pt++){
		sumV[pt]+=-1*V[pt];
	    }
	    delete tmpgr;
	}
    }
    tree->GetEntry(0);
    length = nsamples;
    for(int ipt=0;ipt<length;ipt++){
	sumV[ipt] = sumV[ipt]/counter;
	if(sumV[ipt]>0) sumVpos[ipt] = sumV[ipt];
	sumT[ipt] = t[ipt]*1e9;
    }
    TCanvas *can = new TCanvas("lin","lin",800,600);
    TGraph *gr = new TGraph(length,sumT,sumV);
    gr->Draw("ap");
    TCanvas *can2 = new TCanvas("log","log",800,600);
    //TGraph *grlog = new TGraph(length,sumT,sumVpos);
    mg->Draw("ap");
    cout<<endl;
    cout<<counter<<endl;
    cout<<"Time: ";
    for(int i=0;i<length;i++){
	cout<<sumT[i]<<"\t";
    }
    cout<<endl;
    cout<<"Amp: ";
    for(int i=0;i<length;i++){
	cout<<sumV[i]<<"\t";
    }
    cout<<endl;
    //grlog->Draw("ap");
    rint.Run(kTRUE);
    return 0;
}
