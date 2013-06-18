#include <TH1I.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TSpectrum.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TVirtualFFT.h>

Int_t    FitGainLog(TH1I* hist, Int_t npks, Float_t *mean, Float_t *sigma,Float_t *errmean=0,Float_t *errsigma=0);
Float_t* GainFromFFT(TH1I* hist);
void     FindPeak(TH1* hist, Float_t *res);
TCanvas* NewCanvas(TString cname,Int_t x,Int_t y);
Float_t* GetDCR(TH1I* hist, Float_t pedestal, Float_t gain,Float_t effgate);
Int_t    FitGainWithKnowledge(TH1I* hist,Int_t npks, Float_t *mean, Float_t *sigma, Float_t kped, Float_t kgn,Float_t *errmean=0,Float_t *errsigma=0);
Float_t* GainFromFFTNoPed(TH1I* hist,Float_t kped, Float_t kgain);
