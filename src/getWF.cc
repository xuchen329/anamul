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

#include "DAQ_ADU.h"

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
    Float_t dcr = 0.;
    if(argc>2){
	dcr = atof(argv[2]);
    }
    
    const Int_t wfchannel = 0; //waveform channel

    //set up tree
    TTree* tree = (TTree*)fin->Get("DAQ_ADU");
    DAQ_ADU* wftree = new DAQ_ADU(tree);
    
    TRint rint("App",0x0,0);
    
//Find amplitude and trigger position of 1 pixel fired
    TH1F* hFilterTime = new TH1F("hFilterTime","hist of time",3000,0.5,300.5);
    TH1F* hFilterAmpl = new TH1F("hFilterAmpl","hist of ampl",601,-500.5,100.5);

    wftree->FindPixelFiringDist(hFilterTime,hFilterAmpl);
    Float_t FireTime = hFilterTime->GetBinCenter(hFilterTime->GetMaximumBin())-2.; //where wf start to rise
    Float_t OnePxAmp = hFilterAmpl->GetBinCenter(hFilterAmpl->GetMaximumBin());

    TF1* ampfit = new TF1("ampfit","gaus",0.5*OnePxAmp,1.2*OnePxAmp);
    hFilterAmpl->Fit(ampfit,"RQI");
    OnePxAmp = ampfit->GetParameter(1);
    Float_t SigmaPx = ampfit->GetParameter(2);
//    TGraph* gr = wftree->GetOnePixelFired(FireTime,OnePxAmp,SigmaPx);
    TH1F* hist = wftree->IntegralPixelFiring(FireTime,OnePxAmp,SigmaPx,dcr);
    ofstream fout("addspes.log",ofstream::app);
    if(fout.tellp()<5){
	fout<<"#";
	fout<<"DCR\tXT\tAP\tfname"<<endl;
    }
    fout<<setprecision(0)<<std::fixed;
    fout<<dcr<<"\t";
    fout<<setprecision(3)<<std::fixed;
    fout<<wftree->GetXTprop()<<"\t";
    fout<<wftree->GetAPprop()<<"\t";
    TString ofname(ROOfname(ROOfname.Last('/'),ROOfname.Sizeof()));
    fout<<ofname<<endl;
//    hFilterTime->Draw();
//    TCanvas* can = new TCanvas("wv","wv",800,600);
//    hFilterAmpl->Draw();
/*



    
    Float_t Time[10000];
    Float_t Ampl[10000];

    TMultiGraph *mg = new TMultiGraph();
    TGraph* gr = NULL;
    TH1F* hist = new TH1F("histofamp","histofamp",301,-250.5,50.5);
    TH1F* hist2 = new TH1F("histofampint","histofampint",1100,-1,10);
    TH1F* histT = new TH1F("histoftime","histoftime",150,0.5,150.5);
    for(Long64_t ientry = 0;ientry<nentries;ientry++){
	tree->GetEntry(ientry);

	Float_t minimum = 0;
	Float_t minimumt=0;
	Float_t left    = 4;
	Float_t right   = 6;
	Bool_t skipflag = 0;
	
	for(Int_t pt=0;pt<nsamples;pt++){
	    Ampl[pt] = (V[pt]-Voffset[channel])/256.;
	    Time[pt] = t[pt]*1e9;
	    if(Time[pt]< 3){
		if(Ampl[pt]<-15){
		    skipflag = 1;
		    break;
		}
	    }
	    if(Time[pt]>left && Time[pt]<right){
		if(Ampl[pt]<minimum){
		    minimum = Ampl[pt];
		    minimumt = Time[pt];
		}
	    }
	    //if(Time[pt] > right) break;
	    
	}
	hist->Fill(minimum);
	histT->Fill(minimumt);
//	if(minimum<-94.43 || minimum>-72.62) continue;
//	if(skipflag) continue;

///	Float_t sum = 0;
//	for(Int_t pt=0;pt<nsamples;pt++){
//	    sum+=Ampl[pt];
//	}
//	hist2->Fill(-1.*sum/3.72939e+04);
/*	if(-1.*sum<10000){
	    gr = new TGraph(nsamples,Time,Ampl);
	    mg->Add((TGraph*)gr->Clone(Form("wf_%d",ientry)));
	    }*/
//	hist->Fill(minimum);
//	cout<<minimum<<endl;
//	gr = new TGraph(nsamples,Time,Ampl);
//	break;
//	mg->Add((TGraph*)gr->Clone(Form("wf_%d",ientry)));
//}
//    TCanvas* can = new TCanvas("wv","wv",800,600);
//    mg->Draw("al");
//    TCanvas* can2 = new TCanvas("hist","hist",800,600);
    hist->Draw();
//    Float_t integral = 0;
//    for(int i=0;i<hist2->GetNbinsX();i++){
//	integral+=hist2->GetBinCenter(i)*hist2->GetBinContent(i);
//}
//    cout<<integral<<endl
//    gr->Draw("al");
    rint.Run(kTRUE);
    return 0;
}
