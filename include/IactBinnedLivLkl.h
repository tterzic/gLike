//////////////////////////////////////////////////////////////////////
// Binned LIV full likelihood
//////////////////////////////////////////////////////////////////////

#ifndef IACTBINNEDLIVLKL
#define IACTBINNEDLIVLKL

#include "Iact1dUnbinnedLkl.h"
#include "JointLkl.h"
#include "PoissonLkl.h"

class IactBinnedLivLkl : public Iact1dUnbinnedLkl, public JointLkl
{
 public:
  
  // constructors
  IactBinnedLivLkl(TString inputString="");
  
  // destructor
  virtual ~IactBinnedLivLkl();

  // Read input dN/dE files and related functions
  Int_t ResetdNdESignal();
  Int_t SetdNdESignalFunction(TString function,Float_t p0=0,Float_t p1=0,Float_t p2=0,Float_t p3=0,Float_t p4=0,Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);
  Int_t AdddNdESignalFunction(TString function,Float_t p0=0,Float_t p1=0,Float_t p2=0,Float_t p3=0,Float_t p4=0,Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);

  virtual TH2F*   GetHdNdEpOn(Bool_t isDifferential=kTRUE,Int_t nbinsE=0,Int_t nbinsT=0)  const;
  virtual TH2F*   GetHdNdEpOff(Bool_t isDifferential=kTRUE,Int_t nbinsE=0,Int_t nbinsT=0) const;

  // print data in the overview
  virtual void PrintOverview(Int_t level=0)  {Lkl::PrintOverview(level);}
  virtual void PrintData(Int_t level=0);

  // prepare material
  Int_t          BuildAndBinOnOffHistos();
  Int_t          ConfigureJointLkl();

 protected:

          Int_t    InterpretInputString(TString inputString);
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();

 private:  
  UInt_t fNEnergyBins;           // Number of enery bins
  UInt_t fMinEnergyBinContent;   // Minimum allowed number of energy bins in On and Off samples (if !=0 allow rebinning)
  UInt_t fNRemovedEnergyBins;    // Number of removed energy bins by rebinning
  UInt_t fNTimeBins;             // Number of time bins
  UInt_t fMinTimeBinContent;     // Minimum allowed number of time bins in On and Off samples (if !=0 allow rebinning)
  UInt_t fNRemovedTimeBins;      // Number of removed time bins by rebinning
  Bool_t fTauEDepFluct;     // if kTRUE (default kFALSE) tau is allowed to fluctuate (according to fTau and fDTau independently in each bin)

  TH2I* fHNOn;                   //-> histogram for On  events vs E' and t
  TH2I* fHNOff;                  //-> histogram for Off events vs E' and t

  Float_t* fOnSampleTime;        //-> array of measured time for On events 
  Float_t* fOffSampleTime;       //-> array of measured time for Off events

  Float_t  fTmin;           // [s] Minimum measured time of considered events
  Float_t  fTmax;           // [s] Maximum measured time of considered events
 
  Int_t    fNFineLEBins;       // number of fine bins for internal histos in energy
  Double_t fFineLEMin;       // minimum log(energy[GeV]) for internal histos
  Double_t fFineLEMax;       // maximum log(energy[GeV]) for internal histos 
  Int_t    fNFineTBins;      // number of fine bins for internal histos in time
  Double_t fFineTMin;        // minimum time[s] for internal histos
  Double_t fFineTMax;        // maximum time[s] for internal histos

  TH2F*    fHdNdESignal ;    //-> dN/dE vs E vs t histogram for signal events
  TH2F*    fHdNdEpSignal;    //-> dN/dE' vs E' vs t histogram for signal events (normalized)
  TH2F*    fHdNdEpSignalOff; //-> dN/dE' vs E' vs t histogram for signal events in the off region (normalized)

  ClassDef(IactBinnedLivLkl,1) // Full binned Likelihood for LIV
};

#endif
