/* ======================================================================== *\
!
!   Author(s): Daniel Kerszberg         01/2019 <mailto:dkerszbegr@ifae.es>
! 
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//
// IactBinnedLivLkl
//
//////////////////////////////////////////////////////////////////////////////

// include c++ needed classes

// include Root needed classes
#include "TMath.h"
#include "TFile.h"
#include "TPRegexp.h"

// include gLike needed classes
#include "IactBinnedLivLkl.h"
#include "IactEventListIrf.h"
#include "PoissonLkl.h"

ClassImp(IactBinnedLivLkl);

using namespace std;

// class name and title
static const TString  gName            = "IactBinnedLivLkl";
static const TString  gTitle           = "Iact Binned Likelihood for LIV";

static const Int_t    gNEnergyBins     = 100;        // default number of energy bins
static const Int_t    gNTimeBins       = 100;        // default number of time bins

// List of free parameters.

static const Int_t    gNPars           = 1;                      // Number of free+nuisance parameters
static const Char_t*  gParName[gNPars] = {"eta"};                // Name of parameters

// -2logL function for minuit
void binnedLivLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
static TMinuit* minuit = NULL;

//////////////////////////////////////////////////////////////////////////////
//
// String constructor
//
IactBinnedLivLkl::IactBinnedLivLkl(TString inputString) :
  Lkl(gNPars,inputString,gName,gTitle), Iact1dUnbinnedLkl(inputString), JointLkl(inputString),
  fTauEDepFluct(kFALSE)
{
  if(InterpretInputString(inputString))
    cout << "IactBinnedLivLkl::IactBinnedLivLkl Warning: there were problems interpreting the input string" << endl;      
}

//////////////////////////////////////////////////////////////////
//
// Interpret the input configuration string
//
// Return 0 in case checks are Ok, 1 otherwise
//
Int_t IactBinnedLivLkl::InterpretInputString(TString inputString)
{
  // Prepeare to interpret inputString
  TPMERegexp re("\\s+");
  TPMERegexp fldre("=");
  TString path          = "";
  TString inputfileName = " (No file has been specified) ";

// split the inputString into the different fields, and check the option and values specified in each of them
  UInt_t nfields = re.Split(inputString);
  for(UInt_t ifield=0;ifield<nfields;ifield++)
    {
      fldre.Split(re[ifield]);
      TString optname = fldre[0];
      if(optname.CompareTo("path",TString::kIgnoreCase)==0)
        path=fldre[1];
      else if(optname.CompareTo("inputfile",TString::kIgnoreCase)==0)
        inputfileName=fldre[1];
    }

  // open and read input files with data and IRFs
  TFile* ifile = new TFile(path+(path==""?"":"/")+inputfileName,"READ");
  IactEventListIrf* dataSet = (IactEventListIrf*) ifile->Get("IactEventListIrf");
  if(!dataSet)
    {
      cout << "Iact1dUnbinnedLkl::InterpretInputString Warning: no IactEventListIrf object in file " << inputfileName << endl;
    }
  else
    {
      // extract data
      Float_t eventOnT,eventOffT;
      dataSet->SetOnBranchAddress("t",&eventOnT);
      dataSet->SetOffBranchAddress("t",&eventOffT);

      fOnSampleTime  = new Float_t[GetNon()];
      fOffSampleTime = new Float_t[GetNoff()];

      for(Int_t i=0;i<GetNon();i++)
        {
          dataSet->GetOnEntry(i);
          fOnSampleTime[i] = eventOnT;
        }
      for(Int_t i=0;i<GetNoff();i++)
        {
          dataSet->GetOffEntry(i);
          fOffSampleTime[i] = eventOffT;
        }

      fTmin       = fOnSampleTime[0];
      fTmax       = fOnSampleTime[GetNon()-1];
    }

  return 0;
}

////////////////////////////////////////////////////////////////
// 
// Destructor
//
IactBinnedLivLkl::~IactBinnedLivLkl()
{
  
}

////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of TMinuit subtleties)
//
void IactBinnedLivLkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN
}

////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the parameters (free+nuisance)
//
void IactBinnedLivLkl::SetFunctionAndPars(Double_t ginit)
{
  fMinuit->SetFCN(binnedLivLkl);  
  fMinuit->SetName(Form("%s_Minuit",GetName()));

  // set and pars initial and step values
  Double_t pStart[gNPars] = {ginit};
  Double_t pDelta[gNPars] = {1};    // Precision of parameters during minimization

  // initialize the free (and nuisance) parameters
  SetParameters(gParName,pStart,pDelta);
}		
		

////////////////////////////////////////////////////////////////
//
// Check that all elements needed for the fit are present, otherwise
// try to create/compute them
//
// Return 0 in case checks are Ok, 1 otherwise
//
Int_t IactBinnedLivLkl::MakeChecks()
{
  if(IsChecked()) return 0;

  // add your checks here and try to mend whatever needs to be mended
  
  SetChecked();
  return 0;
}

////////////////////////////////////////////////////////////////
//
// Take the list of On and Off events and make histograms out of 
// them with the specified binning.
// The minimum number of entries per bin in any of the histograms
// is fMinBinContent (default 10), which can be changed through
// SetMinBinContent.
// If fMinBinContent>0 (default 0),
// the binning is changed so that the condition on
// the minimum number of entries is fulfilled.
//
// Return 0 in case of success, 1 otherwise
//
Int_t IactBinnedLivLkl::BuildAndBinOnOffHistos()
{
  if(!fHNOn || !fHNOff)
    {
      // Get the E' vs t distribution for On and Off events
      TH2F* provHNOn  = GetHdNdEpOn(kFALSE,fNEnergyBins,fNTimeBins);
      TH2F* provHNOff = GetHdNdEpOff(kFALSE,fNEnergyBins,fNTimeBins);

      if(!provHNOn || !provHNOff)
        {
          cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Warning: problems creating On and/or Off histograms, On = " << provHNOn << ", Off = " << provHNOff << "." << endl;
          if(provHNOn)  delete provHNOn;
          if(provHNOff) delete provHNOff;
          return 1;
        }

      Bool_t done=kFALSE;
      fNRemovedEnergyBins=0;
      fNRemovedTimeBins=0;
/*      // Rebin if necessary and requested
      if(fMinBinContent>0)
        {
          UInt_t nnewbins;
          Double_t* newbin = new Double_t[fNBins+1];
          GetRebinning(provHNOn,provHNOff,fMinBinContent,nnewbins,newbin);

          if(nnewbins<fNBins)
            {
              fNRemovedBins = fNBins-nnewbins;
              cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Message: fHNON/fHNOff original number of bins = " << fNBins << ", rebinned to " << nnewbins << " bins to keep minimum of " << fMinBinContent << " events per bin" << endl;

              // Get the rebinned On/Off histograms
              TH1F* hRebinOn  =  (TH1F*)  provHNOn->Rebin(nnewbins,"hRebinOn", newbin);
              TH1F* hRebinOff =  (TH1F*) provHNOff->Rebin(nnewbins,"hRebinOff",newbin);
              hRebinOn->SetDirectory(0);
              hRebinOff->SetDirectory(0);

              // replace the fHNOn and fHNOff histograms by the rebinned ones
              fHNOn  = new TH1I("fHNOn", "E' distribution of On events", nnewbins,newbin);
              fHNOff = new TH1I("fHNOff","E' distribution of Off events",nnewbins,newbin);
              fHNOn->SetDirectory(0);
              fHNOff->SetDirectory(0);
              for(Int_t ibin=0;ibin<nnewbins;ibin++)
                {
                  fHNOn->SetBinContent(ibin+1,Int_t(hRebinOn->GetBinContent(ibin+1)));
                  fHNOff->SetBinContent(ibin+1,Int_t(hRebinOff->GetBinContent(ibin+1)));
                }

              delete hRebinOn;
              delete hRebinOff;
              done = kTRUE;
            }

          delete [] newbin;
        }
*/
      if(!done)
        {
          fHNOn  = new TH2I("fHNOn", "E' vs t distribution of On events",fNTimeBins,provHNOn->GetXaxis()->GetXmin(),provHNOn->GetXaxis()->GetXmax(),fNEnergyBins,provHNOn->GetYaxis()->GetXmin(),provHNOn->GetYaxis()->GetXmax());
          fHNOff = new TH2I("fHNOff","E' vs t distribution of Off events",fNTimeBins,provHNOff->GetXaxis()->GetXmin(),provHNOff->GetXaxis()->GetXmax(),fNEnergyBins,provHNOff->GetYaxis()->GetXmin(),provHNOff->GetYaxis()->GetXmax());
          fHNOn->SetDirectory(0);
          fHNOff->SetDirectory(0);

          for(Int_t ibin=0;ibin<fNTimeBins;ibin++)
            {
              for(Int_t jbin=0;jbin<fNEnergyBins;jbin++)
                {
                  fHNOn->SetBinContent(ibin+1,jbin+1,Int_t(provHNOn->GetBinContent(ibin+1,jbin+1)));
                  fHNOff->SetBinContent(ibin+1,jbin+1,Int_t(provHNOff->GetBinContent(ibin+1,jbin+1)));
                }
            }
        }

      if(provHNOn)  delete provHNOn;
      if(provHNOff) delete provHNOff;
    }

  // check that everything went fine before leaving
  if(!fHNOn || !fHNOff)
    {
      cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Warning: problems rebinning On and/or Off histograms" << endl;
      return 1;
    }

/*  // check that every bin has more than the minimum required number of events (fMinBinContent)
  Bool_t binsAreOk = kTRUE;
  for(Int_t ibin=0;ibin<fNBins-fNRemovedBins;ibin++)
    {
      Int_t onevts  = fHNOn->GetBinContent(ibin+1);
      Int_t offevts = fHNOff->GetBinContent(ibin+1);
      if(onevts<fMinBinContent)
        {
          binsAreOk = kFALSE;
          cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Warning: number of On events in bin #" << ibin+1
               << " (" << onevts << ") is below the minimum requested (" << fMinBinContent << ")" << endl;
        }
      if(offevts<fMinBinContent)
        {
          binsAreOk = kFALSE;
          cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Warning: number of Off events in bin #" << ibin+1
               << " (" << offevts << ") is below the minimum requested (" << fMinBinContent << ")" << endl;
        }
    }


  if(!binsAreOk) return 1;*/

  return 0;
}

////////////////////////////////////////////////////////////////
//
// Feed the JointLkl with one PoissonLkl per E' bin 
//
// Return 0 in case of success, 1 otherwise
//
Int_t IactBinnedLivLkl::ConfigureJointLkl()
{
  if(!fHNOn || !fHNOff)
    {
      cout << "Iact1dBinnedLkl::ConfigureJointLkl (" << GetName() << ") Warning: On and/or Off histograms do no exist" << endl;
      return 1;
    }

  ClearSampleList(); // delete samples in the list

  // fill the Poisson Lkl classes for every bin in the list
  Double_t totalw = 0;
  for(Int_t ibin=0;ibin<fNTimeBins-fNRemovedTimeBins;ibin++)
    {
      for(Int_t jbin=0;jbin<fNEnergyBins-fNRemovedEnergyBins;jbin++)
        {
          PoissonLkl* pLkl = new PoissonLkl(fHNOn->GetBinContent(ibin+1,jbin+1),fHNOff->GetBinContent(ibin+1,jbin+1),GetTau(),(fTauEDepFluct? GetDTau() : 0),0,Form("%s_BinT_%d_BinE_%d",GetName(),ibin,jbin));
          Double_t lemin  = fHNOn->GetBinLowEdge(ibin+1);
          Double_t lemax  = fHNOn->GetBinLowEdge(ibin+1)+fHNOn->GetBinWidth(ibin+1);
    
          Double_t weight = 0.;//IntegrateLogE(GetHdNdEpSignal(),lemin,lemax);
          pLkl->SetUnitsOfG(weight>0? weight : 0);
    
          // is there signal leakeage in the Off region?
/*          if(GetHdNdEpSignalOff())
            {
              Double_t gfractioninoff = IntegrateLogE(GetHdNdEpSignalOff(),lemin,lemax)*GetdNdEpSignalOffIntegral()/(weight*GetdNdEpSignalIntegral());
              pLkl->SetGFractionInOff(gfractioninoff);
            }
    
          // are there gamma-ray foreground events in the signal region?
          if(GetHdNdEpFrg())
            {
              Double_t nFrgEvts = IntegrateLogE(GetHdNdEpFrg(),lemin,lemax)*GetdNdEpFrgIntegral();
              pLkl->SetFrgNEvents(nFrgEvts);
            }
    
          // is b a fixed (as opposed to nuisance) parameter? mostly for tests
          if(fKnownBackground) pLkl->SetKnownBackground();
  */  
    
          AddSample(pLkl);
    
          totalw+=weight;
        }
    }
  SetOwner(); // so that all created pLkl are deleted with the JointLkl object

  // check that the sum of all weights is approximately equal to 1
  if(TMath::Abs(totalw-1)>0.01)
    {
      cout << "Iact1dBinnedLkl::ConfigureJointLkl (" << GetName() << ") Warning: total sum of weights should be 1! (It's " << totalw << "), too many bins maybe??" << endl;
      return 1;
    }

  return 0;
}

//////////////////////////////////////////////////////////////////
//
// Produce the E' vs t distribution of On events and return 
// the corresponding histogram.
// If fHNOn does not exist, then it creates it. Otherwise, return 
// a copy of fHNOn. If isDifferential=kTRUE (default is kTRUE), 
// return the dN/dE' distribution (i.e. number of entries in the 
// bin divided by the size of the bin in E' unit). Otherwise the 
// returned histogram is the number of events per bin of E'.
// The returned histogram will not be deleted by the destructor
// of IactBinnedLivLkl.
//
TH2F* IactBinnedLivLkl::GetHdNdEpOn(Bool_t isDifferential,Int_t nbinsE,Int_t nbinsT) const
{
  if(!fHNOn || (nbinsE>0 && nbinsE!=fNEnergyBins-fNRemovedEnergyBins) || (nbinsT>0 && nbinsT!=fNTimeBins-fNRemovedTimeBins))
    {
      // we need a positive number of bins
      if(nbinsE<=0) nbinsE = gNEnergyBins;
      if(nbinsT<=0) nbinsT = gNTimeBins;
    
      // create histo
      TH2F* h = new TH2F("dNdEpOn","dN/dE' vs t for On events",nbinsT,fTmin,fTmax,nbinsE,TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));
      h->SetDirectory(0);
      h->SetXTitle("t [s]");
      h->SetYTitle("log_{10}(E' [GeV])");
      h->SetZTitle("dN/dE' [GeV^{-1}]");
   
      const Float_t *onSample = GetOnSample(); 
      // fill histo
      for(Int_t i=0;i<GetNon();i++)
        h->Fill(fOnSampleTime[i],onSample[i]);
    
      // divide by bin width
      if(h->GetEntries()>0)
        if(isDifferential)
          for(Int_t ibin=0;ibin<nbinsT;ibin++)
            {
              for(Int_t jbin=0;jbin<nbinsE;jbin++)
                {
                  Double_t leminbin = h->GetYaxis()->GetBinLowEdge(jbin+1);
                  Double_t lemaxbin = leminbin+h->GetYaxis()->GetBinWidth(jbin+1);
                  Double_t deltaE   = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin);
                  h->SetBinError(ibin+1,jbin+1,h->GetBinError(ibin+1,jbin+1)/deltaE);
                  h->SetBinContent(ibin+1,jbin+1,h->GetBinContent(ibin+1,jbin+1)/deltaE);
                }
            }
      return h;
    }

  // check if bin width is constant
  Bool_t binWidthIsConstant = kTRUE;
  for(Int_t ibin=0;ibin<fNEnergyBins-fNRemovedEnergyBins-1;ibin++)
    if(TMath::Abs(Float_t(fHNOn->GetYaxis()->GetBinWidth(ibin+1))-Float_t(fHNOn->GetYaxis()->GetBinWidth(ibin+2)))>1e-9)
      binWidthIsConstant = kFALSE;

  for(Int_t ibin=0;ibin<fNTimeBins-fNRemovedTimeBins-1;ibin++)
    if(TMath::Abs(Float_t(fHNOn->GetXaxis()->GetBinWidth(ibin+1))-Float_t(fHNOn->GetXaxis()->GetBinWidth(ibin+2)))>1e-9)
      binWidthIsConstant = kFALSE;

  TH2F* h;
  if(binWidthIsConstant)
    h = new TH2F("dNdEpOn","dN/dE' vs t for On events",fNTimeBins-fNRemovedTimeBins,fHNOn->GetXaxis()->GetXmin(),fHNOn->GetXaxis()->GetXmax(),fNEnergyBins-fNRemovedEnergyBins,fHNOn->GetYaxis()->GetXmin(),fHNOn->GetYaxis()->GetXmax());
  else
    h = new TH2F("dNdEpOn","dN/dE' vs t for On events",fNTimeBins-fNRemovedTimeBins,fHNOn->GetXaxis()->GetXbins()->GetArray(),fNEnergyBins-fNRemovedEnergyBins,fHNOn->GetYaxis()->GetXbins()->GetArray());

  h->SetDirectory(0);
  h->SetXTitle("t [s]");
  h->SetYTitle("log_{10}(E' [GeV])");
  h->SetZTitle("dN/dE' [GeV^{-1}]");
  for(Int_t ibin=0;ibin<fNTimeBins-fNRemovedTimeBins;ibin++)
    {
      for(Int_t jbin=0;jbin<fNEnergyBins-fNRemovedEnergyBins;jbin++)
        {
          h->SetBinContent(ibin+1,jbin+1,Float_t(fHNOn->GetBinContent(ibin+1,jbin+1)));
          h->SetBinError(ibin+1,jbin+1,TMath::Sqrt(fHNOn->GetBinContent(ibin+1,jbin+1)));
        }
    }

   // divide by bin width
  if(isDifferential)
    for(Int_t ibin=0;ibin<fNTimeBins-fNRemovedTimeBins;ibin++)
      {
        for(Int_t jbin=0;jbin<fNEnergyBins-fNRemovedEnergyBins;jbin++)
          {
            Double_t leminbin = h->GetBinLowEdge(jbin+1);
            Double_t lemaxbin = leminbin+h->GetBinWidth(jbin+1);
            Double_t deltaE   = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin);
            h->SetBinError(ibin+1,jbin+1,h->GetBinError(ibin+1,jbin+1)/deltaE);
            h->SetBinContent(ibin+1,jbin+1,h->GetBinContent(ibin+1,jbin+1)/deltaE);
          }
      }

  return h;
}

//////////////////////////////////////////////////////////////////
//
// Produce the E' vs t distribution of Off events and return 
// the corresponding histogram.
// If fHNOff does not exist, then it creates it. Otherwise, return 
// a copy of fHNOff. If isDifferential=kTRUE (default is kTRUE), 
// return the dN/dE' distribution (i.e. number of entries in the 
// bin divided by the size of the bin in E' unit). Otherwise the 
// returned histogram is the number of events per bin of E'.
// The returned histogram will not be deleted by the destructor
// of IactBinnedLivLkl.
//
TH2F* IactBinnedLivLkl::GetHdNdEpOff(Bool_t isDifferential,Int_t nbinsE,Int_t nbinsT) const
{
  if(!fHNOff || (nbinsE>0 && nbinsE!=fNEnergyBins-fNRemovedEnergyBins) || (nbinsT>0 && nbinsT!=fNTimeBins-fNRemovedTimeBins))
    {
      // we need a positive number of bins
      if(nbinsE<=0) nbinsE = gNEnergyBins;
      if(nbinsT<=0) nbinsT = gNTimeBins;
    
      // create histo
      TH2F* h = new TH2F("dNdEpOff","dN/dE' vs t for Off events",nbinsT,fTmin,fTmax,nbinsE,TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));
      h->SetDirectory(0);
      h->SetXTitle("t [s]");
      h->SetYTitle("log_{10}(E' [GeV])");
      h->SetZTitle("dN/dE' [GeV^{-1}]");
    
      const Float_t *offSample = GetOffSample(); 
      // fill histo
      for(Int_t i=0;i<GetNoff();i++)
        h->Fill(fOffSampleTime[i],offSample[i]);
    
      // divide by bin width
      if(h->GetEntries()>0)
        if(isDifferential)
          for(Int_t ibin=0;ibin<nbinsT;ibin++)
            {
              for(Int_t jbin=0;jbin<nbinsE;jbin++)
                {
                  Double_t leminbin = h->GetYaxis()->GetBinLowEdge(jbin+1);
                  Double_t lemaxbin = leminbin+h->GetYaxis()->GetBinWidth(jbin+1);
                  Double_t deltaE   = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin);
                  h->SetBinError(ibin+1,jbin+1,h->GetBinError(ibin+1,jbin+1)/deltaE);
                  h->SetBinContent(ibin+1,jbin+1,h->GetBinContent(ibin+1,jbin+1)/deltaE);
                }
            }
      return h;
    }

  // check if bin width is constant
  Bool_t binWidthIsConstant = kTRUE;
  for(Int_t ibin=0;ibin<fNEnergyBins-fNRemovedEnergyBins-1;ibin++)
    if(TMath::Abs(Float_t(fHNOff->GetYaxis()->GetBinWidth(ibin+1))-Float_t(fHNOff->GetYaxis()->GetBinWidth(ibin+2)))>1e-9)
      binWidthIsConstant = kFALSE;

  for(Int_t ibin=0;ibin<fNTimeBins-fNRemovedTimeBins-1;ibin++)
    if(TMath::Abs(Float_t(fHNOff->GetXaxis()->GetBinWidth(ibin+1))-Float_t(fHNOff->GetXaxis()->GetBinWidth(ibin+2)))>1e-9)
      binWidthIsConstant = kFALSE;

  TH2F* h;
  if(binWidthIsConstant)
    h = new TH2F("dNdEpOff","dN/dE' vs t for Off events",fNTimeBins-fNRemovedTimeBins,fHNOff->GetXaxis()->GetXmin(),fHNOff->GetXaxis()->GetXmax(),fNEnergyBins-fNRemovedEnergyBins,fHNOff->GetYaxis()->GetXmin(),fHNOff->GetYaxis()->GetXmax());
  else
    h = new TH2F("dNdEpOff","dN/dE' vs t for Off events",fNTimeBins-fNRemovedTimeBins,fHNOff->GetXaxis()->GetXbins()->GetArray(),fNEnergyBins-fNRemovedEnergyBins,fHNOff->GetYaxis()->GetXbins()->GetArray());

  h->SetDirectory(0);
  h->SetXTitle("t [s]");
  h->SetYTitle("log_{10}(E' [GeV])");
  h->SetZTitle("dN/dE' [GeV^{-1}]");
  for(Int_t ibin=0;ibin<fNTimeBins-fNRemovedTimeBins;ibin++)
    {
      for(Int_t jbin=0;jbin<fNEnergyBins-fNRemovedEnergyBins;jbin++)
        {
          h->SetBinContent(ibin+1,jbin+1,Float_t(fHNOff->GetBinContent(ibin+1,jbin+1)));
          h->SetBinError(ibin+1,jbin+1,TMath::Sqrt(fHNOff->GetBinContent(ibin+1,jbin+1)));
        }
    }

   // divide by bin width
  if(isDifferential)
    for(Int_t ibin=0;ibin<fNTimeBins-fNRemovedTimeBins;ibin++)
      {
        for(Int_t jbin=0;jbin<fNEnergyBins-fNRemovedEnergyBins;jbin++)
          {
            Double_t leminbin = h->GetBinLowEdge(jbin+1);
            Double_t lemaxbin = leminbin+h->GetBinWidth(jbin+1);
            Double_t deltaE   = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin);
            h->SetBinError(ibin+1,jbin+1,h->GetBinError(ibin+1,jbin+1)/deltaE);
            h->SetBinContent(ibin+1,jbin+1,h->GetBinContent(ibin+1,jbin+1)/deltaE);
          }
      }

  return h;
}

//////////////////////////////////////////////////////////////////
// 
// Recreate a fresh new version of fHdNdESignal
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t IactBinnedLivLkl::ResetdNdESignal()
{
  // Delete existing fHdNdESignal and create empty one
  if(fHdNdESignal)
    delete fHdNdESignal;

  // Create histo
  fHdNdESignal = new TH2F("fHdNdESignal","dN/dE vs t for signal events",fNFineTBins,fFineTMin,fFineTMax,fNFineLEBins,fFineLEMin,fFineLEMax);
  fHdNdESignal->SetDirectory(0);
  fHdNdESignal->SetXTitle("t [s]");
  fHdNdESignal->SetYTitle("log_{10}(E [GeV])");
  fHdNdESignal->SetZTitle("dN/dE [GeV^{-1}]");

  // Delete existing fHdNdEpSignal and fHdNdEpSignalOff 
  if(fHdNdEpSignal)
    {
      delete fHdNdEpSignal;
      fHdNdEpSignal=NULL;
    }

  if(fHdNdEpSignalOff)
  {
    delete fHdNdEpSignalOff;
    fHdNdEpSignalOff=NULL;
  }
  SetChecked(kFALSE);

  // exit
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Create and save a dN/dE histogram for signal according to 
// given <function> and parameters <p0>, <p1>, ...
// (see Iact1dUnbinnedLkl::AdddNdESignalFunction for details on available
// functions)
//
// If fdNdESignal exists it will be replaced
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t IactBinnedLivLkl::SetdNdESignalFunction(TString function,Float_t p0,Float_t p1,Float_t p2,Float_t p3,Float_t p4,Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{
  // Delete existing fHdNdESignal and create empty one
  ResetdNdESignal();

  // call to add the function
  AdddNdESignalFunction(function,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9);

  // exit
  return 0;
}

//////////////////////////////////////////////////////////////////
// 
// Add to a previously existing dN/dE histogram for signal according to 
// given <function> and parameters <p0>, <p1>, ...
//
// If <function>=="powerlawLIV"
// <p0> = "delay" eta (linear: GeV/s ; quadraitc: GeV^2/s)
// <p1> = dependency of the delay to the energy (1: linear; 2: quadratic)
// <p2> = E0, the energy scale
// <p3> = [N(t)], the array of normalizations (depending of time)
// <p4> = [alpha(t)], the array of slopes (depending of time)
//
// If <function>=="exponentialDecayLIV"
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t IactBinnedLivLkl::AdddNdESignalFunction(TString function,Float_t p0,Float_t p1,Float_t p2,Float_t p3,Float_t p4,Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{
  return 1;
}

////////////////////////////////////////////////////////////////////////
//
// Likelihood function (-2logL) 
// To be minimized by TMinuit
// Free parameters:
// par[0] = g 
//
// this is a trivial funcion that just returs 0, so probably not a good idea
// to try to minimize it. You must replace this but your likelihood function
//
void binnedLivLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;

  // get internal object, histos, values, etc
  Double_t g       = par[0];
  f = -2*TMath::Log(1);
}		
