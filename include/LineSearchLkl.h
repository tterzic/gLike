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
  inline TF1*    GetFdNdEpBkg()         {return fFdNdEpBkg;}

 protected:
          Int_t    InterpretInputString(TString inputString);
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();
  inline Int_t     SetfBkg(TF1* FdNdEpBkg) {fFdNdEpBkg   = FdNdEpBkg;  return 0;}
  Int_t ComputeBkgModelFromOnHisto();

 private:
  TF1*    fFdNdEpBkg   ;    //-> dN/dE'dt vs E' for background model (normalized)

  ClassDef(LineSearchLkl,1) // Line Search Likelihood
};

#endif
