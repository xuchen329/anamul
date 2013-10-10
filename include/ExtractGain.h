#include <TH1I.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TSpectrum.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <TVirtualFFT.h>

Int_t    FitGainLog(TH1I* hist, Int_t npks, Float_t *mean, Float_t *sigma,Float_t *errmean=0,Float_t *errsigma=0,Int_t pedes=0);
Int_t    FitGainLogHighNoise(TH1I* hist, Int_t npks, Float_t *mean, Float_t *sigma,Float_t *errmean=0,Float_t *errsigma=0);
Float_t* GainFromFFT(TH1I* hist);
Float_t* GainFromFFTShifted(TH1I* hist,Float_t ped);
void     FindPeak(TH1* hist, Int_t *res);
TCanvas* NewCanvas(TString cname,Int_t x,Int_t y);
Float_t* GetDCR(TH1I* hist, Float_t pedestal, Float_t gain,Float_t effgate);
Float_t* GetDCRHighNoise(TH1I* hist, Float_t pedestal, Float_t gain,Float_t effgate);
Int_t    FitGainWithKnowledge(TH1I* hist,Int_t npks, Float_t *mean, Float_t *sigma, Float_t kped, Float_t kgn,Float_t *errmean=0,Float_t *errsigma=0);
Float_t* GainFromFFTNoPed(TH1I* hist,Float_t kped, Float_t kgain);
void FindValey(TH1F* hist, Int_t cnt, Float_t *px, Int_t* idx, Float_t *left, Float_t *right);

Float_t FindLocalMaximumPos(TH1* hist, Int_t startbin=0, Int_t endbin=0);

Float_t GainFromAutoCor(TH1I* hist);
Float_t FindHarmonyPeak(TH1F* hist);
