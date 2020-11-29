//////////////////////////////////////////////////////////////////////
// Unbinned LIV likelihood
//////////////////////////////////////////////////////////////////////

#ifndef IACTUNBINNEDLIVLKL
#define IACTUNBINNEDLIVLKL

#include "Lkl.h"
#include "Iact1dUnbinnedLkl.h"

class IactUnbinnedLivLkl : public Iact1dUnbinnedLkl//, public virtual Lkl
{
 public:

  // constructors
  IactUnbinnedLivLkl(TString inputString="");
  
  // destructor
  virtual ~IactUnbinnedLivLkl();

  // Read input dN/dE files and related functions
  Int_t ResetdNdESignal();
  Int_t SetdNdESignalFunction(TString function,Float_t p0=0,Float_t p1=0,Float_t p2=0,Float_t p3=0,Float_t p4=0,Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);
  Int_t AdddNdESignalFunction(TString function,Float_t p0=0,Float_t p1=0,Float_t p2=0,Float_t p3=0,Float_t p4=0,Float_t p5=0,Float_t p6=0,Float_t p7=0,Float_t p8=0,Float_t p9=0);

  virtual TH2D*   GetHdNdEpOn(Bool_t isDifferential=kTRUE,Int_t nbinsE=0,Int_t nbinsT=0)  const;
  virtual TH2D*   GetHdNdEpOff(Bool_t isDifferential=kTRUE,Int_t nbinsE=0,Int_t nbinsT=0) const;

  // print data in the overview
  virtual void PrintOverview(Int_t level=0)  {Lkl::PrintOverview(level);}
  virtual void PrintData(Int_t level=0) {Lkl::PrintData(level);};
  virtual void SetUnitsOfG(Double_t unit) {Lkl::SetUnitsOfG(unit);};

  // getters
  inline       Double_t  GetTmin()             const {return fTMin;}
  inline       Double_t  GetTmax()             const {return fTMax;}
  inline const Double_t* GetOnSampleEnergy()     const {return fOnSampleEnergy;}
  inline const Double_t* GetOnSampleTime()     const {return fOnSampleTime;}
  inline const Double_t* GetOffSampleTime()    const {return fOffSampleTime;}

  inline const TH2D*    GetHdNdEpSignal()     const {return fHdNdEpSignal;}
  inline const TH2D*    GetHdNdEpBkg()        const {return fHdNdEpBkg;}
  Double_t GetdNdEpSignalIntegral()    {CheckHistograms(kFALSE); if(!fHdNdEpSignal) return 0; return fHdNdEpSignal->GetBinContent(0);}

  Int_t        NormalizedNdEHisto(TH2D* histo);

  //virtual Int_t SimulateDataSamples(UInt_t seed=0,Float_t meanG=0);
  //virtual Int_t GetRealBkgAndGoffHistos(TRandom3* rdm,TH2F*& hdNdEpBkg,TH2F*& hdNdEpSignalOff);

 protected:

          Int_t    InterpretInputString(TString inputString);
  virtual void     SetFunctionAndPars(Double_t ginit=0);
  virtual Int_t    MakeChecks();
  virtual void     SetMinuitLink();
          Int_t    CheckHistograms(Bool_t checkdNdEpBkg=kTRUE);

 private:  
  Double_t* fOnSampleEnergy;     //-> array of measured energy for On events 
  Double_t* fOnSampleTime;       //-> array of measured time for On events 
  Double_t* fOffSampleTime;      //-> array of measured time for Off events

  Double_t  fTMin;               // [s] Minimum measured time of considered events
  Double_t  fTMax;               // [s] Maximum measured time of considered events
 
  Int_t    fNFineLEBins;         // number of fine bins for internal histos in energy
  Double_t fFineLEMin;           // minimum log(energy[GeV]) for internal histos
  Double_t fFineLEMax;           // maximum log(energy[GeV]) for internal histos 
  Int_t    fNFineTBins;          // number of fine bins for internal histos in time
  Double_t fFineTMin;            // minimum time[s] for internal histos
  Double_t fFineTMax;            // maximum time[s] for internal histos

  TH2D*    fHdNdESignal;         //-> dN/dE vs E vs t histogram for signal events
  TH2D*    fHdNdEpSignal;        //-> dN/dE' vs E' vs t histogram for signal events (normalized)
  TH2D*    fHdNdEpSignalOff;     //-> dN/dE' vs E' vs t histogram for signal events in the off region (normalized)
  TH2D*    fHdNdEBkg;            //-> dN/dE' vs E' vs t histogram for bkbkgts (normalized)
  TH2D*    fHdNdEpBkg;           //-> dN/dE' vs E' vs t histogram for bkbkgts (normalized)

  ClassDef(IactUnbinnedLivLkl,1) // Unbinned Likelihood for LIV
};

#endif
