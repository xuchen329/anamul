//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 17 16:07:44 2013 by ROOT version 5.99/01
// from TTree data/data
// found on file: outtest.root
//////////////////////////////////////////////////////////

#ifndef DaqMul_h
#define DaqMul_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1I.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class DaqMul {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           qdcch;
    /*
   Int_t           qdcch0;
   Int_t           qdcch1;
   Int_t           qdcch2;
   Int_t           qdcch3;
   Int_t           qdcch4;
   Int_t           qdcch5;
   Int_t           qdcch6;
   Int_t           qdcch7;
   Int_t           qdcch8;
   Int_t           qdcch9;
   Int_t           qdcch10;
   Int_t           qdcch11;
   Int_t           qdcch12;
   Int_t           qdcch13;
   Int_t           qdcch14;
   Int_t           qdcch15;
   Int_t           qdcch16;
   Int_t           qdcch17;
   Int_t           qdcch18;
   Int_t           qdcch19;
   Int_t           qdcch20;
   Int_t           qdcch21;
   Int_t           qdcch22;
   Int_t           qdcch23;
   Int_t           qdcch24;
   Int_t           qdcch25;
   Int_t           qdcch26;
   Int_t           qdcch27;
   Int_t           qdcch28;
   Int_t           qdcch29;
   Int_t           qdcch30;
   Int_t           qdcch31;*/
   Float_t         temperature;
   Float_t         time;
   Float_t         date;
   Float_t         voltage;
   Float_t         gate_width;
   Int_t           id;

   // List of branches
   TBranch        *b_qdcch;   //!
    /*  TBranch        *b_qdcch0;   //!
   TBranch        *b_qdcch1;   //!
   TBranch        *b_qdcch2;   //!
   TBranch        *b_qdcch3;   //!
   TBranch        *b_qdcch4;   //!
   TBranch        *b_qdcch5;   //!
   TBranch        *b_qdcch6;   //!
   TBranch        *b_qdcch7;   //!
   TBranch        *b_qdcch8;   //!
   TBranch        *b_qdcch9;   //!
   TBranch        *b_qdcch10;   //!
   TBranch        *b_qdcch11;   //!
   TBranch        *b_qdcch12;   //!
   TBranch        *b_qdcch13;   //!
   TBranch        *b_qdcch14;   //!
   TBranch        *b_qdcch15;   //!
   TBranch        *b_qdcch16;   //!
   TBranch        *b_qdcch17;   //!
   TBranch        *b_qdcch18;   //!
   TBranch        *b_qdcch19;   //!
   TBranch        *b_qdcch20;   //!
   TBranch        *b_qdcch21;   //!
   TBranch        *b_qdcch22;   //!
   TBranch        *b_qdcch23;   //!
   TBranch        *b_qdcch24;   //!
   TBranch        *b_qdcch25;   //!
   TBranch        *b_qdcch26;   //!
   TBranch        *b_qdcch27;   //!
   TBranch        *b_qdcch28;   //!
   TBranch        *b_qdcch29;   //!
   TBranch        *b_qdcch30;   //!
   TBranch        *b_qdcch31;   //!*/
   TBranch        *b_temperature;   //!
   TBranch        *b_time;   //!
   TBranch        *b_date;   //!
   TBranch        *b_voltage;   //!
   TBranch        *b_gate_width;   //!
   TBranch        *b_id;   //!

    DaqMul(TTree *tree=0);
    virtual ~DaqMul();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    virtual void     GetHistogram(TH1I* hist);
    virtual void     GetCondition(Float_t *rett);
};

#endif

#ifdef DaqMul_cxx
DaqMul::DaqMul(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("outtest.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("outtest.root");
      }
      f->GetObject("data",tree);

   }
   Init(tree);
}

DaqMul::~DaqMul()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DaqMul::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DaqMul::LoadTree(Long64_t entry)
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

void DaqMul::Init(TTree *tree)
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

   fChain->SetBranchAddress("qdcch", &qdcch, &b_qdcch);
/*   fChain->SetBranchAddress("qdcch0", &qdcch0, &b_qdcch0);
   fChain->SetBranchAddress("qdcch1", &qdcch1, &b_qdcch1);
   fChain->SetBranchAddress("qdcch2", &qdcch2, &b_qdcch2);
   fChain->SetBranchAddress("qdcch3", &qdcch3, &b_qdcch3);
   fChain->SetBranchAddress("qdcch4", &qdcch4, &b_qdcch4);
   fChain->SetBranchAddress("qdcch5", &qdcch5, &b_qdcch5);
   fChain->SetBranchAddress("qdcch6", &qdcch6, &b_qdcch6);
   fChain->SetBranchAddress("qdcch7", &qdcch7, &b_qdcch7);
   fChain->SetBranchAddress("qdcch8", &qdcch8, &b_qdcch8);
   fChain->SetBranchAddress("qdcch9", &qdcch9, &b_qdcch9);
   fChain->SetBranchAddress("qdcch10", &qdcch10, &b_qdcch10);
   fChain->SetBranchAddress("qdcch11", &qdcch11, &b_qdcch11);
   fChain->SetBranchAddress("qdcch12", &qdcch12, &b_qdcch12);
   fChain->SetBranchAddress("qdcch13", &qdcch13, &b_qdcch13);
   fChain->SetBranchAddress("qdcch14", &qdcch14, &b_qdcch14);
   fChain->SetBranchAddress("qdcch15", &qdcch15, &b_qdcch15);
   fChain->SetBranchAddress("qdcch16", &qdcch16, &b_qdcch16);
   fChain->SetBranchAddress("qdcch17", &qdcch17, &b_qdcch17);
   fChain->SetBranchAddress("qdcch18", &qdcch18, &b_qdcch18);
   fChain->SetBranchAddress("qdcch19", &qdcch19, &b_qdcch19);
   fChain->SetBranchAddress("qdcch20", &qdcch20, &b_qdcch20);
   fChain->SetBranchAddress("qdcch21", &qdcch21, &b_qdcch21);
   fChain->SetBranchAddress("qdcch22", &qdcch22, &b_qdcch22);
   fChain->SetBranchAddress("qdcch23", &qdcch23, &b_qdcch23);
   fChain->SetBranchAddress("qdcch24", &qdcch24, &b_qdcch24);
   fChain->SetBranchAddress("qdcch25", &qdcch25, &b_qdcch25);
   fChain->SetBranchAddress("qdcch26", &qdcch26, &b_qdcch26);
   fChain->SetBranchAddress("qdcch27", &qdcch27, &b_qdcch27);
   fChain->SetBranchAddress("qdcch28", &qdcch28, &b_qdcch28);
   fChain->SetBranchAddress("qdcch29", &qdcch29, &b_qdcch29);
   fChain->SetBranchAddress("qdcch30", &qdcch30, &b_qdcch30);
   fChain->SetBranchAddress("qdcch31", &qdcch31, &b_qdcch31);*/
   fChain->SetBranchAddress("temperature", &temperature, &b_temperature);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("date", &date, &b_date);
   fChain->SetBranchAddress("voltage", &voltage, &b_voltage);
   fChain->SetBranchAddress("gate_width", &gate_width, &b_gate_width);
   fChain->SetBranchAddress("id", &id, &b_id);
   Notify();
}

Bool_t DaqMul::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DaqMul::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DaqMul::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DaqMul_cxx
