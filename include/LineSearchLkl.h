//////////////////////////////////////////////////////////////////////
// Likelihood for line search
//////////////////////////////////////////////////////////////////////

#ifndef LINESEARCHLKL
#define LINESEARCHLKL

#include "TCanvas.h"

#include "Iact1dUnbinnedLkl.h"

class LineSearchLkl : public Iact1dUnbinnedLkl
{
 public:
  
  // constructors
  LineSearchLkl(TString inputString="");
  
  // destructor
  virtual ~LineSearchLkl();

  // Plots
  TCanvas* PlotHistosAndData();
  TCanvas* PlotforWindow();
  inline TF1*    GetFdNdEpBkg()         {return fFdNdEpBkg;}
  // getters
  inline Int_t GetEventsInEnergyWindow()             const {return fEventsInEnergyWindow;}

 protected:
          Int_t    InterpretInputString(TString inputString);
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();
  inline Int_t     SetfBkg(TF1* FdNdEpBkg) {fFdNdEpBkg   = FdNdEpBkg;  return 0;}
  Int_t ComputeBkgModelFromOnHisto();

 private:

  Int_t fEventsInEnergyWindow;
  TF1*    fFdNdEpBkg   ;    //-> dN/dE'dt vs E' for background model (normalized)
  ClassDef(LineSearchLkl,1) // Line Search Likelihood
};

#endif
