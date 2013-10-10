#define DaqMul_cxx
#include <iostream>
#include "DaqMul.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1I.h>


using namespace std;

void DaqMul::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L DaqMul.C
//      Root > DaqMul t
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
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      //nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

   }
}


void DaqMul::GetHistogram(TH1I *hist, Int_t cut){
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
    
    //TH1I* hist = new TH1I("QDC","QDC spectrum",4096,1,4096);
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	//nb = fChain->GetEntry(jentry);   nbytes += nb;
	// if (Cut(ientry) < 0) continue;
	b_qdcch->GetEntry(jentry);
//	hist->Fill(qdcch);
	if(qdcch>=cut) hist->Fill(qdcch);	
    }
}

void DaqMul::GetCondition(Float_t*rett){
    if (fChain == 0) return;   
    Long64_t nentries = fChain->GetEntriesFast();
    
    TH1F* histtemp = new TH1F("temp","temp",510,0,50);
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
	if (ientry < 0) break;
	b_temperature->GetEntry(jentry);
	if(temperature>0) histtemp->Fill(temperature);
    }
    LoadTree(0);
    b_voltage->GetEntry(0);
//    Float_t *rett = new Float_t[3];   
    rett[0] = histtemp->GetMean();
    rett[1] = histtemp->GetRMS();
    rett[2] = voltage;
    cout<<"Temperature :"<<rett[0]<<" +/- "<<rett[1]<<endl;
    cout<<"Voltage     :"<<rett[2]<<endl;
    delete histtemp;
    return;
}
