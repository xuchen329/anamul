//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Aug 30 11:10:34 2013 by ROOT version 5.34/06
// from TTree DAQ_ADU/--No Title--
// found on file: ap-0.root
//////////////////////////////////////////////////////////

#ifndef DAQ_ADU_h
#define DAQ_ADU_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class DAQ_ADU {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          event;
   UInt_t          channel;
   UInt_t          nsamples;
   Float_t         V[10000];   //[nsamples]
   Float_t         t[10000];   //[nsamples]
   Float_t         Vgain[4];
   Float_t         Voffset[4];

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_nsamples;   //!
   TBranch        *b_V;   //!
   TBranch        *b_t;   //!
   TBranch        *b_Vgain;   //!
   TBranch        *b_Voffset;   //!

   DAQ_ADU(TTree *tree=0);
   virtual ~DAQ_ADU();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

    //user functions
    void FindPixelFiringDist(TH1F* histoftime,TH1F* histofampl);
    Int_t IsOnePixelFired(Float_t tturn,Float_t onepx,Float_t onepxsig);
    Bool_t IsStrictOnePixelFired(Float_t tturn,Float_t onepx,Float_t onepxsig);
    TGraph* GetOnePixelFired(Float_t tturn,Float_t onepx,Float_t onepxsig);
    TH1F* IntegralPixelFiring(Float_t tturn,Float_t onepx,Float_t onepxsig,Float_t dcr=0.);
    Float_t xtprop;
    Float_t approp;
    Float_t integratedpixelfired;
    Float_t validevents;

    inline Float_t GetXTprop(){return xtprop;}
    inline Float_t GetAPprop(){return approp;}
    inline Float_t GetIntegralPixel(){return integratedpixelfired;}
    inline Float_t GetNevents(){return validevents;}
};

#endif

#ifdef DAQ_ADU_cxx
DAQ_ADU::DAQ_ADU(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ap-0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ap-0.root");
      }
      f->GetObject("DAQ_ADU",tree);

   }
   Init(tree);
}

DAQ_ADU::~DAQ_ADU()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DAQ_ADU::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DAQ_ADU::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DAQ_ADU::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("nsamples", &nsamples, &b_nsamples);
   fChain->SetBranchAddress("V", V, &b_V);
   fChain->SetBranchAddress("t", t, &b_t);
   fChain->SetBranchAddress("Vgain", Vgain, &b_Vgain);
   fChain->SetBranchAddress("Voffset", Voffset, &b_Voffset);
   Notify();
}

Bool_t DAQ_ADU::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DAQ_ADU::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DAQ_ADU::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DAQ_ADU_cxx
