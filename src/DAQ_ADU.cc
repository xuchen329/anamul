#define DAQ_ADU_cxx
#include <iostream>
#include "DAQ_ADU.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>

using namespace std;

void DAQ_ADU::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L DAQ_ADU.C
//      Root > DAQ_ADU t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}

void DAQ_ADU::FindPixelFiringDist(TH1F* histoftime,TH1F* histofampl){
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	fChain->GetEntry(jentry);
	if (channel!=0) continue;
	Float_t minimum = 0;
	Float_t mintime = 0;
	
	for(Int_t pt=0;pt<nsamples;pt++){
	    Float_t a = (V[pt]-Voffset[channel])/256.;
	    if(a<minimum){
		minimum= a;
		mintime= t[pt]*1e9;
	    }
	}
	//cout<<"MinA: "<<minimum<<" MinT: "<<mintime<<"        \r"<<flush;
	histoftime->Fill(mintime);
	histofampl->Fill(minimum);
    }
}

Int_t DAQ_ADU::IsOnePixelFired(Float_t tturn,Float_t onepx,Float_t onepxsig){
    for(Int_t pt=0;pt<nsamples;pt++){
	Float_t time = t[pt]*1e9;
	Float_t amplitude = (V[pt]-Voffset[channel])/256.;
	if(time<=tturn){
	    if(TMath::Abs(amplitude)>0.5*TMath::Abs(onepx)) return -1;
	}
	if(time>tturn && time<tturn+7){
	    if(TMath::Abs(amplitude)>(1.5*TMath::Abs(onepx))) return 0;
	}
    }
    return 1;
}

Bool_t DAQ_ADU::IsStrictOnePixelFired(Float_t tturn,Float_t onepx,Float_t onepxsig){
    for(Int_t pt=0;pt<nsamples;pt++){
	Float_t time = t[pt]*1e9;
	Float_t amplitude = (V[pt]-Voffset[channel])/256.;
	if(time<=tturn){
	    if(TMath::Abs(amplitude)>TMath::Abs(0.5*onepx)) return 0;
	}
	if(time>tturn && time<tturn+7){
	    if(TMath::Abs(amplitude)>(TMath::Abs(onepx)+3.*onepxsig)) return 0;
	}
	if(time>tturn+25){
	    if(TMath::Abs(amplitude)>(TMath::Abs(onepx)*0.5)) return 0;
	}
    }
    return 1;
}

TGraph* DAQ_ADU::GetOnePixelFired(Float_t tturn,Float_t onepx,Float_t onepxsig){
    if (fChain == 0) return 0;

    Long64_t nentries = fChain->GetEntries();

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	fChain->GetEntry(jentry);
	if (channel!=0) continue;
	if (IsOnePixelFired(tturn,onepx,onepxsig)){
	    Float_t *time = new Float_t[nsamples];
	    Float_t *amp = new Float_t[nsamples];
	    for(Int_t pt=0;pt<nsamples;pt++){
		time[pt] = t[pt]*1e9;
		amp[pt] = (V[pt]-Voffset[channel])/256.;
	    }
	    TGraph* gr = new TGraph(nsamples,time,amp);
	    return (TGraph*)gr->Clone("onewv");
	}
    }
    return 0;
}

TH1F* DAQ_ADU::IntegralPixelFiring(Float_t tturn,Float_t onepx,Float_t onepxsig,Float_t dcr){
    if (fChain == 0) return 0;
    Int_t counter = 0;
    Int_t notvalid=0;
    TH1F* rethist = new TH1F("hist","hist",75000,-75000.5,-0.5);
    Long64_t nentries = fChain->GetEntries();
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	fChain->GetEntry(jentry);
	if (channel!=0) continue;
	Int_t retv = IsOnePixelFired(tturn,onepx,onepxsig);
	if (retv>0){
	    counter++;
	    Float_t timeunit = (t[1]-t[0])*1e9;
	    Float_t sumamp = 0;
	    for(Int_t pt=0;pt<nsamples;pt++){
		sumamp+=(V[pt]-Voffset[channel])/256.;
	    }
	    rethist->Fill(sumamp*timeunit);
	}
	if(retv<0) notvalid++;
    }
    xtprop = (float(nentries-notvalid-counter)/float(nentries-notvalid));
    cout<<" XT prop: "<<100.*xtprop<<"%"<<endl;

    TH1F* hist= (TH1F*)rethist->Rebin(15,"hnew");
    Float_t fitcenter = hist->GetBinCenter(hist->GetMaximumBin());
    TF1* func = new TF1("func","gaus",1.3*fitcenter,0.5*fitcenter);
    hist->Fit(func,"RQI");
    Float_t onepix = func->GetParameter(1);
    integratedpixelfired = 0;
    for(int i=0;i<75000;i++){
	integratedpixelfired+=(rethist->GetBinCenter(i)/onepix)*rethist->GetBinContent(i);
    }

    Float_t timewindow = 150.-tturn-7.;
    validevents = float(counter);
    cout<<"Total N pixel fired: "<<integratedpixelfired<<endl;
    cout<<"Total events: "<<counter<<endl;
    cout<<"window: "<<timewindow<<endl;
    approp = integratedpixelfired/(validevents*(1+timewindow*1e-9*dcr*(1+xtprop)))-1.;
    cout<<"AP prop: "<<approp<<endl;
    return (TH1F*)rethist->Clone("integraled");
}
