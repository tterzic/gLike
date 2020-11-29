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

#include <iomanip>

ClassImp(IactBinnedLivLkl);

using namespace std;

// class name and title
static const TString  gName            = "IactBinnedLivLkl";
static const TString  gTitle           = "Iact Binned Likelihood for LIV";

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
  fTauEDepFluct(kFALSE), fMinEnergyBinContent(gDefMinBinContent), fMinTimeBinContent(gDefMinBinContent), fNEnergyBins(gDefNBins/*gNEnergyBins*/), fNTimeBins(gDefNBins/*gNTimeBins*/), fNRemovedEnergyBins(0), fNRemovedTimeBins(0), fTmin(gTmin), fTmax(gTmax), fHNOn(NULL), fHNOff(NULL), fOnSampleTime(NULL), fOffSampleTime(NULL)
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
      cout << "Iact1dUnbinnedLkl::InterpretInputString Warning: THERE is IactEventListIrf object in file " << inputfileName << endl;

  TF1 *f2 = new TF1("f2", "[0]*TMath::Exp(-0.5*((x-[1])/[2])**2.)/([2]*TMath::Sqrt(2*TMath::Pi()))", 0., 100.);
  f2->SetParameters(1.,50.,10.);

      // extract data

      // extract info from file 
      fTmin = 0.;//dataSet->GetEpmin();
      fTmax = 10.;//dataSet->GetEpmax();

      Double_t eventOnT,eventOffT;
      dataSet->SetOnBranchAddressD("t",&eventOnT);
      dataSet->SetOffBranchAddressD("t",&eventOffT);

      fOnSampleTime  = new Double_t[GetNon()];
      fOffSampleTime = new Double_t[GetNoff()];

      for(Int_t i=0;i<GetNon();i++)
        {
          dataSet->GetOnEntry(i);
	  if(i==0) fTmin = (eventOnT-58497.)*86400.;
	  if(i==(GetNon()-1)) fTmax = (eventOnT-58497.)*86400.;
	  cout << setprecision(20) << " on " << i << " t = " << eventOnT << "in days or " << eventOnT*86400. << " in sec" << endl;
          fOnSampleTime[i] = (eventOnT-58497.)*24.*60.*60.;
	  if(fOnSampleTime[i]==-1) fOnSampleTime[i] = (i+1)*(90./GetNon());
	  //if(fOnSampleTime[i]==-1) fOnSampleTime[i] = f2->GetRandom();
        }
      for(Int_t i=0;i<GetNoff();i++)
        {
          dataSet->GetOffEntry(i);
	  cout << setprecision(20) << " off " << i << " t = " << eventOffT << endl;
          fOffSampleTime[i] = (eventOffT-58497.)*24.*60.*60.;
	  if(fOffSampleTime[i]==-1) fOffSampleTime[i] = (i+1.5)*(90./GetNon());
	  //if(fOffSampleTime[i]==-1) fOffSampleTime[i] = f2->GetRandom();
        }

      fTmin       = fOnSampleTime[0];
      fTmax       = fOnSampleTime[GetNon()-1];
    }

      cout << "fTmin = " << fTmin << endl;
      cout << "fTmax = " << fTmax << endl;
      cout << "fEmin = " << GetEmin() << endl;
      cout << "fEmax = " << GetEmax() << endl;
      BuildAndBinOnOffHistos();
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

  // Check the IactBinnedLivLkl specific part
  /////////////////////////////////////////
  // get the dN/dE' histograms for On and Off and check binning
  if(BuildAndBinOnOffHistos())
    {
      cout << "IactBinnedLivLkl::MakeChecks (" << GetName() << ") Warning: problems building On and/or Off histos!" << endl;
      return 1;
    }

cout << "Coucou" << endl;

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
      TH2F* provHNOn  = GetHdNdEpOn(kFALSE,10,10);
      TH2F* provHNOff = GetHdNdEpOff(kFALSE,10,10);
      //TH2F* provHNOn  = GetHdNdEpOn(kFALSE,fNEnergyBins,fNTimeBins);
      //TH2F* provHNOff = GetHdNdEpOff(kFALSE,fNEnergyBins,fNTimeBins);

      if(!provHNOn || !provHNOff)
        {
          cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Warning: problems creating On and/or Off histograms, On = " << provHNOn << ", Off = " << provHNOff << "." << endl;
          if(provHNOn)  delete provHNOn;
          if(provHNOff) delete provHNOff;
          return 1;
        }

provHNOn->SaveAs("./provHNOn.root");
provHNOff->SaveAs("./provHNOff.root");

      Bool_t done=kFALSE;
      fNRemovedEnergyBins=0;
      fNRemovedTimeBins=0;
      // Rebin if necessary and requested
      if(fMinEnergyBinContent>0)
        {
          UInt_t nnewbins;
          Double_t* newbin = new Double_t[fNEnergyBins+1];
          cout << "Rebinning: " << fNEnergyBins << " --> " << nnewbins << endl;
          GetRebinning(provHNOn,provHNOff,fMinEnergyBinContent,nnewbins,newbin);
          for (int lol=0; lol<nnewbins+1; lol++) cout << "final rebinning = " << lol << " bin = " << newbin[lol] << endl;

          cout << "Rebinning: " << fNEnergyBins << " --> " << nnewbins << endl;

          if(nnewbins<fNEnergyBins)
            {
              fNRemovedEnergyBins = fNEnergyBins-nnewbins;
              cout << "Iact1dBinnedLkl::BuildAndBinOnOffHistos (" << GetName() << ") Message: fHNON/fHNOff original number of energy bins = " << fNEnergyBins << ", rebinned to " << nnewbins << " bins to keep minimum of " << fMinEnergyBinContent << " events per energy bin" << endl;

              // Get the rebinned On/Off histograms
              TH2F* hRebinOn = new TH2F("hRebinOn", "E' vs t distribution of On events",fNTimeBins,provHNOn->GetXaxis()->GetXmin(),provHNOn->GetXaxis()->GetXmax(),nnewbins,newbin); //  =  (TH1F*)  provHNOn->Rebin(nnewbins,"hRebinOn", newbin);
              TH2F* hRebinOff = new TH2F("hRebinOff", "E' vs t distribution of Off events",fNTimeBins,provHNOff->GetXaxis()->GetXmin(),provHNOff->GetXaxis()->GetXmax(),nnewbins,newbin); // =  (TH1F*) provHNOff->Rebin(nnewbins,"hRebinOff",newbin);
              for(Int_t ibin=0;ibin<fNTimeBins;ibin++)
                {
                  for(Int_t jbin=0;jbin<fNEnergyBins;jbin++)
                    {
                      hRebinOn->Fill(provHNOn->GetXaxis()->GetBinCenter(ibin+1),provHNOn->GetYaxis()->GetBinCenter(jbin+1),provHNOn->GetBinContent(ibin+1,jbin+1));
                      hRebinOff->Fill(provHNOff->GetXaxis()->GetBinCenter(ibin+1),provHNOff->GetYaxis()->GetBinCenter(jbin+1),provHNOff->GetBinContent(ibin+1,jbin+1));
                    }
                }
              hRebinOn->SetDirectory(0);
              hRebinOff->SetDirectory(0);

              // replace the fHNOn and fHNOff histograms by the rebinned ones
              fHNOn  = new TH2I("fHNOn", "E' vs t distribution of On events",fNTimeBins,provHNOn->GetXaxis()->GetXmin(),provHNOn->GetXaxis()->GetXmax(),nnewbins,newbin);
              fHNOff = new TH2I("fHNOff","E' vs t distribution of Off events",fNTimeBins,provHNOff->GetXaxis()->GetXmin(),provHNOff->GetXaxis()->GetXmax(),nnewbins,newbin);
              fHNOn->SetDirectory(0);
              fHNOff->SetDirectory(0);

              for(Int_t ibin=0;ibin<fNTimeBins;ibin++)
                {
                  for(Int_t jbin=0;jbin<nnewbins;jbin++)
                    {
                      fHNOn->SetBinContent(ibin+1,jbin+1,Int_t(hRebinOn->GetBinContent(ibin+1,jbin+1)));
                      fHNOff->SetBinContent(ibin+1,jbin+1,Int_t(hRebinOff->GetBinContent(ibin+1,jbin+1)));
                    }
                }

              delete hRebinOn;
              delete hRebinOff;
              done = kTRUE;
            }

          delete [] newbin;
        }

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
fHNOn->SaveAs("./fHNOn.root");
fHNOff->SaveAs("./fHNOff.root");


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
          //Double_t lemin  = fHNOn->GetBinLowEdge(ibin+1);
          //Double_t lemax  = fHNOn->GetBinLowEdge(ibin+1)+fHNOn->GetBinWidth(ibin+1);
   
	  Double_t sum = 0.;
          for(Int_t k = 0; k < fNTimeBins-fNRemovedTimeBins; k++) {
            sum += fHdNdESignal->GetBinContent(k+1,jbin+1);
	  }

          Double_t weight = 1./(fHdNdESignal->GetBinContent(ibin+1,jbin+1)*fHdNdESignal->GetXaxis()->GetBinWidth(ibin+1)*fHdNdESignal->GetYaxis()->GetBinWidth(jbin+1)/sum);
          //Double_t weight = fHdNdESignal->GetBinContent(ibin+1,jbin+1)/fHdNdESignal->GetEntries();
          //Double_t weight = 1./(fNTimeBins*fNEnergyBins);//IntegrateLogE(GetHdNdEpSignal(),lemin,lemax);
          pLkl->SetUnitsOfG(weight>0? weight : 0);
    
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
      if(nbinsE<=0) nbinsE = gDefNBins; //gNEnergyBins;
      if(nbinsT<=0) nbinsT = gDefNBins; //gNTimeBins;

cout << "TEST1" << endl;
cout << nbinsT << " " << fTmin << " " << fTmax << " " << nbinsE << " " << GetEmin() << " " << GetEmax() << endl;

      // create histo
      //TH2F* h = new TH2F("dNdEpOn","dN/dE' vs t for On events",nbinsT,fTmin,fTmax,nbinsE,TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));
      TH2F* h = new TH2F("dNdEpOn","dN/dE' vs t for On events",nbinsT,fTmin,fTmax,nbinsE,TMath::Log10(150.),TMath::Log10(2500.));
      h->SetDirectory(0);
      h->SetXTitle("t [s]");
      h->SetYTitle("log_{10}(E' [GeV])");
      h->SetZTitle("dN/dE' [GeV^{-1}]");
   
      const Float_t *onSample = Iact1dUnbinnedLkl::GetOnSample(); 
      // fill histo
      for(Int_t i=0;i<GetNon();i++)
        {
		cout << "onSample[i] = " << onSample[i] << endl;
          h->Fill(fOnSampleTime[i],onSample[i]);
	}
    
cout << "TEST123" << endl;
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
cout << "TEST1234" << endl;

          /*for(Int_t ibin=0;ibin<nbinsT;ibin++)
            {
              for(Int_t jbin=0;jbin<nbinsE;jbin++)
                {
                  
		}
	    }*/

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
      if(nbinsE<=0) nbinsE = gDefNBins; //gNEnergyBins;
      if(nbinsT<=0) nbinsT = gDefNBins; //gNTimeBins;
    
cout << "TEST2" << endl;

      // create histo
      //TH2F* h = new TH2F("dNdEpOff","dN/dE' vs t for Off events",nbinsT,fTmin,fTmax,nbinsE,TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));
      TH2F* h = new TH2F("dNdEpOff","dN/dE' vs t for Off events",nbinsT,fTmin,fTmax,nbinsE,TMath::Log10(150.),TMath::Log10(2500.));
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
cout << "TEST234" << endl;
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

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 
// Simulate list of On and Off events
// according to the total pdf described by fHdNdEpBkg and fHdNdEpSignal
// and the observation time in fObsTime and tau taking randomly from a
// gaussian of mean fTau and mean fDTau
// seed    (default 0) = seed for random generator
// meanG   (default 0) = mean G value
//
// IF meanG<0, do not simulate independent ON events, use OFF sample
// also as ON
//
// Return 0 in case of success
//        1 otherwise
//
Int_t IactBinnedLivLkl::SimulateDataSamples(UInt_t seed,Float_t meanG)
{
  if(meanG<0) meanG=0;
  
  // Sanity checks

  // compute weights for different pdf components
  TRandom3* rdm     = new TRandom3(seed);
  TRandom*  saverdm = gRandom;
  gRandom = rdm;

  TH2F* realHdNdEpBkg       = NULL;
  TH2F* realHdNdEpSignalOff = NULL;

  if(GetRealBkgAndGoffHistos(rdm,realHdNdEpBkg,realHdNdEpSignalOff)) return 1;

  //Float_t meanB    = GetdNdEpBkgIntegral()*fObsTime;  
  Float_t meanB    = 200;  
  Float_t meanBoff = realHdNdEpBkg->GetBinContent(0)*GetObsTime()*GetTau();
  //Float_t meanF    = GetdNdEpFrgIntegral()*fObsTime;
  //Float_t meanGoff = ((fHdNdEpSignal && realHdNdEpSignalOff)? meanG*realHdNdEpSignalOff->GetBinContent(0)/GetdNdEpSignalIntegral() : 0); 
  //Float_t meanNon  = meanB+meanF+meanG;
  //Float_t meanNoff = meanBoff+meanGoff;
  Float_t meanNon  = meanB+meanG;
  Float_t meanNoff = meanBoff;

  // setup histogram to build pdfs
  TH2F* hBkgBinIntegrated     = new TH2F("hBkgBinIntegrated",     "Histogram for Off background event generation",          fNFineLEBins,fFineLEMin,fFineLEMax,fNFineTBins,fFineTMin,fFineTMax);
  TH2F* hRealBkgBinIntegrated = new TH2F("hRealBkgBinIntegrated", "Histogram for  On background event generation",          fNFineLEBins,fFineLEMin,fFineLEMax,fNFineTBins,fFineTMin,fFineTMax);
  TH2F* hFrgBinIntegrated     = new TH2F("hFrgBinIntegrated",     "Histogram for foreground event generation",              fNFineLEBins,fFineLEMin,fFineLEMax,fNFineTBins,fFineTMin,fFineTMax);
  TH2F* hSigBinIntegrated     = new TH2F("hSigBinIntegrated",     "Histogram for signal event generation",                  fNFineLEBins,fFineLEMin,fFineLEMax,fNFineTBins,fFineTMin,fFineTMax);
  TH2F* hSigOffBinIntegrated  = new TH2F("hSigOffBinIntegrated",  "Histogram for signal event generation in the Off region",fNFineLEBins,fFineLEMin,fFineLEMax,fNFineTBins,fFineTMin,fFineTMax);
  TH2F* hOnBinIntegrated      = new TH2F("hOnBinIntegrated",      "Histogram for  On event generation",                     fNFineLEBins,fFineLEMin,fFineLEMax,fNFineTBins,fFineTMin,fFineTMax);
  TH2F* hOffBinIntegrated     = new TH2F("hOffBinIntegrated",     "Histogram for Off event generation",                     fNFineLEBins,fFineLEMin,fFineLEMax,fNFineTBins,fFineTMin,fFineTMax);
  
  hBkgBinIntegrated->SetDirectory(0);
  hRealBkgBinIntegrated->SetDirectory(0);
  hFrgBinIntegrated->SetDirectory(0);
  hSigBinIntegrated->SetDirectory(0);
  hSigOffBinIntegrated->SetDirectory(0);
  hOnBinIntegrated->SetDirectory(0);
  hOffBinIntegrated->SetDirectory(0);

  // build pdf
  /*for(Int_t ibin = 0;ibin<fNFineBins;ibin++)
    {
      Double_t leminbin = hBkgBinIntegrated->GetBinLowEdge(ibin+1);
      Double_t lemaxbin = leminbin+hBkgBinIntegrated->GetBinWidth(ibin+1);
      if(TMath::Power(10,leminbin) < fEpmax && TMath::Power(10,lemaxbin) > fEpmin)
	{
	  Float_t deltaE   = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin); 
	  hBkgBinIntegrated->SetBinContent(ibin+1,fHdNdEpBkg->GetBinContent(ibin+1)*deltaE);
	  hRealBkgBinIntegrated->SetBinContent(ibin+1,realHdNdEpBkg->GetBinContent(ibin+1)*deltaE);
	  if(fHdNdEpSignal)
	    hSigBinIntegrated->SetBinContent(ibin+1,fHdNdEpSignal->GetBinContent(ibin+1)*deltaE);
	  else
	    hSigBinIntegrated->SetBinContent(ibin+1,0);
	  if(fHdNdEpSignalOff)
	    hSigOffBinIntegrated->SetBinContent(ibin+1,realHdNdEpSignalOff->GetBinContent(ibin+1)*deltaE);
	  else
	    hSigOffBinIntegrated->SetBinContent(ibin+1,0);
	  if(fHdNdEpFrg)
	    hFrgBinIntegrated->SetBinContent(ibin+1,fHdNdEpFrg->GetBinContent(ibin+1)*deltaE);
	  else
	    hFrgBinIntegrated->SetBinContent(ibin+1,0);
	}
      else
	{
	  hBkgBinIntegrated->SetBinContent(ibin+1,0);
	  hRealBkgBinIntegrated->SetBinContent(ibin+1,0);
	  hSigBinIntegrated->SetBinContent(ibin+1,0);
	  hSigOffBinIntegrated->SetBinContent(ibin+1,0);
	  hFrgBinIntegrated->SetBinContent(ibin+1,0);
	}
    }*/


}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Histograms hdNdEpBkg and hdNdEpSignalOff will contain same entries as fHdNdEpBkg and fHdNdEpSignalOff,
// respectively, but with normalization fluctuating according to the uncertainty in tau
// Note that fHdNdEpBkg is the expected distribution of background events in the On region
// and fHdNdEpSignalOff the expected distribution of signal events in the total Off region
// (total meaning that if tau=3 the effective area to consder is that of the three subregions)
//
Int_t IactBinnedLivLkl::GetRealBkgAndGoffHistos(TRandom3* rdm,TH2F*& hdNdEpBkg,TH2F*& hdNdEpSignalOff) 
{
  // create new histos with contents of the existing ones
  //if(fHdNdEpBkg) hdNdEpBkg = new TH1F(*fHdNdEpBkg);
  if(fHdNdEpSignal) hdNdEpBkg = new TH2F(*fHdNdEpSignal);
  else
    {
      cout << "IactBinnedLivLkl::GetRealBkgAndGoffHistos Warning: fHdNdEpBkg histo does not exist" << endl;
      return 1;
    }
  if(fHdNdEpSignalOff) hdNdEpSignalOff = new TH2F(*fHdNdEpSignalOff);

  // if no uncertainty in tau, that's all we need to do
  if(GetDTau()<=0) return 0;

  // chose the true tau for this simulated sample according to the pdf
  //GetTrueTau() = rdm->Gaus(GetTau(),GetDTau());
  
  if(GetTrueTau()>0)
    {
      hdNdEpBkg->SetBinContent(0,hdNdEpBkg->GetBinContent(0)*GetTrueTau()/GetTau());
      if(hdNdEpSignalOff)
	hdNdEpSignalOff->SetBinContent(0,hdNdEpSignalOff->GetBinContent(0)*GetTrueTau()/GetTau());
    }
  else
    {
      cout << "IactBinnedLivLkl::GetRealBkgAndGoffHistos Error: negative or null tau value, we cannot work like that!" << endl;
      return 1;
    }
  
  return 0;
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
// If <function>=="gaussLIV"
// <p0> = "delay" eta (linear: GeV/s ; quadraitc: GeV^2/s)
// <p1> = dependency of the delay to the energy (1: linear; 2: quadratic)
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

  if(function=="constantLIV")
    {
	    cout << "constantLIV" << endl;
      Float_t Delay       = p0;
      Float_t EDependency = p1;
      //Float_t Emin        = p2;
      //Float_t Emax        = p3;
      Float_t Power       = p2;
      //Float_t Mean        = p5;
      //Float_t StdDev      = p6;
      //Float_t Tmin        = p7;
      //Float_t Tmax        = p8;
      Int_t nevents = 100000;
      TF1 *f1 = new TF1("f1", "[0]*x**(-[1])", GetEmin()/1000., GetEmax()/1000.);
      f1->SetParameters(1.,Power);
      f1->SetNpx(10000);
      Double_t E2 = 0.;
   
      fHdNdESignal = new TH2F("fHdNdESignal","fHdNdESignal",fNTimeBins,(1)*(90./nevents),(nevents)*(90./nevents),fNEnergyBins,GetEmin()/1000.,GetEmax()/1000.);
      for (int i=0; i<nevents; i++)
        {
          E2 = f1->GetRandom();
          fHdNdESignal->Fill((i+1)*(90./nevents),E2);
        }
      TCanvas *c1 = new TCanvas("c1","c1",900,300);
      c1->Divide(3,3);
      c1->cd(1);
      f1->Draw();
      c1->cd(3);
      fHNOn->Draw("COLZ");
      c1->cd(4);
      fHNOff->Draw("COLZ");
      c1->cd(5);
      fHdNdESignal->Draw("COLZ");
      c1->cd(6);
      c1->cd(7);
      c1->cd(8);
      c1->cd(9);
    }

  if(function=="gaussLIV")
    {
	    cout << "gaussLIV" << endl;

  Float_t Delay       = p0;
  Float_t EDependency = p1;
  Float_t Emin        = p2;
  Float_t Emax        = p3;
  Float_t Power       = p4;
  Float_t Mean        = p5;
  Float_t StdDev      = p6;
  Float_t Tmin        = p7;
  Float_t Tmax        = p8;

  TF1 *f1 = new TF1("f1", "[0]*x**(-[1])", Emin, Emax);
  f1->SetParameters(1.,Power);
  f1->SetNpx(10000);
  //f1->SetRange(1.,10000.);

  TF1 *f2 = new TF1("f2", "[0]*TMath::Exp(-0.5*((x-[1])/[2])**2.)/([2]*TMath::Sqrt(2*TMath::Pi()))", Tmin, Tmax);
  f2->SetParameters(1.,Mean,StdDev);

  Double_t E = 0.;
  Double_t E2 = 0.;
  Double_t E3 = 0.;
  Double_t t = 0.;
  Double_t t_shifted = 0.;
  Double_t t_shifted2 = 0.;
  Int_t nevents = 100000;
  TH1D *h1 = new TH1D("h1","Random Distribution from GetRandom()",20,Emin,Emax);
  TH1D *h2 = new TH1D("h2","Random Distribution from GetRandom()",50,Tmin,Tmax);
  TH1D *h3 = new TH1D("h3","Random Distribution from GetRandom() + LIV shift",50,Tmin,Tmax);
  TH2D *h4 = new TH2D("h4","E vs Random Distribution from GetRandom()",6,Tmin,Tmax,6,Emin,Emax);
  TH2D *h5 = new TH2D("h5","E vs Random Distribution from GetRandom() + LIV shift",6,Tmin,Tmax,6,Emin,Emax);
  //TH2D *h4 = new TH2D("h4","E vs Random Distribution from GetRandom()",50,Tmin,Tmax,20,Emin,Emax);
  //TH2D *h5 = new TH2D("h5","E vs Random Distribution from GetRandom() + LIV shift",50,Tmin,Tmax,20,Emin,Emax);
  fHNOn = new TH2I("fHNOn","fHNOn",6,Tmin,Tmax,6,Emin,Emax);
  fHNOff = new TH2I("fHNOff","fHNOff",6,Tmin,Tmax,6,Emin,Emax);
  fHdNdESignal = new TH2F("fHdNdESignal","fHdNdESignal",6,Tmin,Tmax,6,Emin,Emax);
  fNFineTBins=6, fFineTMin=Tmin,fFineTMax=Tmax,fNFineLEBins=6,fFineLEMin=Emin,fFineLEMax=Emax;
  fNTimeBins=6,fNEnergyBins=6;
  fNRemovedTimeBins=0,fNRemovedEnergyBins=0;
  //fUnitsOfG = (Tmax-Tmin)*(Emax-Emin);
  //SetUnitsOfG((Tmax-Tmin)*(Emax-Emin));
  //ResetdNdESignal();
  for (int i=0; i<nevents; i++)
    {
      E = f1->GetRandom();
      E2 = f1->GetRandom();
      E3 = f1->GetRandom();
      t = f2->GetRandom();
      //t_shifted = f2->GetRandom();
      t_shifted2 = f2->GetRandom() + Delay*TMath::Power(E,EDependency);
      t_shifted = t + Delay*TMath::Power(E,EDependency);
      h1->Fill(E);
      h2->Fill(t);
      h3->Fill(t_shifted);
      h4->Fill(t,E);
      h5->Fill(t_shifted,E);
      //fHdNdESignal->Fill(t,E);
      fHdNdESignal->Fill(t_shifted,E2);
      //fHdNdEpSignal->Fill(t_shifted2,E3);
      fHNOn->Fill(t_shifted,E2);
      fHNOff->Fill(t,E);
   }

  TCanvas *c1 = new TCanvas("c1","c1",900,300);
  c1->Divide(3,3);
  c1->cd(1);
  f1->Draw();
  c1->cd(2);
  f2->Draw();
  c1->cd(3);
  fHNOn->Draw("COLZ");
  c1->cd(4);
  h1->Draw();
  c1->cd(5);
  h2->Draw();
  c1->cd(6);
  h3->Draw();
  c1->cd(7);
  fHdNdESignal->Draw("COLZ");
  c1->cd(8);
  h4->Draw("COLZ");
  c1->cd(9);
  h5->Draw("COLZ");

  Float_t log10Emin = TMath::Log10(Emin);
  Float_t log10Emax = TMath::Log10(Emax);
    }
  return 0;
}

//////////////////////////////////////////////////////////////////
// Fill newbin with optimal binning which contains at least minnevts
// in all bins of hOn and hOff, and fill inewbin with the number of
// bins
void IactBinnedLivLkl::GetRebinning(TH2F* hOn,TH2F* hOff,UInt_t minnevts,UInt_t& inewbin,Double_t* newbin)
{
  Int_t nibinsT=hOn->GetNbinsX();
  Int_t nibinsE=hOn->GetNbinsY();

  inewbin = 0;
  newbin[0] = hOn->GetYaxis()->GetBinLowEdge(1);
  Int_t inewbin2 = 0;
  Double_t* newbin2 = new Double_t[fNEnergyBins+1];
  newbin2[0] = hOn->GetYaxis()->GetBinLowEdge(1);
  Int_t inewbin3 = 0;
  Double_t* newbin3 = new Double_t[fNEnergyBins+1];
  newbin3[0] = hOn->GetYaxis()->GetBinLowEdge(1);

  cout << "ca commence num bin x = " << nibinsT << " num bin y = " << nibinsE << endl;

  for(Int_t ibinT=0;ibinT<nibinsT;ibinT++)
    {
  cout << "on change de bin de T: " << ibinT << endl;
      Int_t tmpinewbin = 0;
      Double_t* tmpnewbin = new Double_t[fNEnergyBins+1];
      tmpnewbin[0] = hOn->GetYaxis()->GetBinLowEdge(1);
      for(Int_t ibinE=0;ibinE<nibinsE;ibinE++,tmpinewbin++)
        {
          Int_t non  = hOn->GetBinContent(ibinT+1,ibinE+1);
          Int_t noff = hOff->GetBinContent(ibinT+1,ibinE+1);
	  cout << "non = " << non << " noff = " << noff << endl;
          while((non<minnevts || noff<minnevts) && ibinE<nibinsE-1)
            {
              ibinE++;
              non  += hOn->GetBinContent(ibinT+1,ibinE+1);
              noff += hOff->GetBinContent(ibinT+1,ibinE+1);
	      cout << " WHILE ACTIVATED non = " << non << " noff = " << noff << endl;
            }

          tmpnewbin[tmpinewbin+1] = hOn->GetYaxis()->GetBinLowEdge(ibinE+1)+hOn->GetYaxis()->GetBinWidth(ibinE+1);

          if((non<minnevts || noff<minnevts) && tmpinewbin>1) // last bin does not comply with minimal statistics condition
            {
              tmpnewbin[tmpinewbin] = tmpnewbin[tmpinewbin+1];
              tmpinewbin--;
            }
        }
      if(ibinT==0)
        {
          inewbin3 = tmpinewbin;
          //newbin3 = tmpnewbin;
          //inewbin2 = 0;
          //inewbin3 = 0;
          inewbin2 = tmpinewbin;
          //newbin2 = tmpnewbin;
          for (int lol=0; lol<tmpinewbin+1; lol++)
            {
               //inewbin3++;
               //inewbin2++;
	       newbin2[lol] = tmpnewbin[lol];
	       newbin3[lol] = tmpnewbin[lol];
               cout << "tessssst lol = " << lol << " bin = " << newbin3[lol] << endl;
            }
        }
      else
        {
          inewbin3 = 0;
          delete [] newbin3;
          newbin3 = new Double_t[fNEnergyBins+1];
          newbin3[0] = hOn->GetYaxis()->GetBinLowEdge(1);
          Int_t i=1, j=1;
          for (int k=1; k<inewbin2+1 && i<inewbin2+1 && j<tmpinewbin+1; k++)
            {
              if(TMath::Abs(newbin2[i]-tmpnewbin[j]) < 0.001)
                {
		     cout << "je suis en 1 et k = " << k << " et i = " << i << " et newbin2 = " << newbin2[i] << " et newbin3 = " << newbin3[k] << " et j = " << j << " et tmpnewbin = " << tmpnewbin[j] << endl;
                  if(k>0 && TMath::Abs(newbin2[i]-newbin3[k-1]) > 0.001)
		    {
		      newbin3[k] = newbin2[i];
		      inewbin3++;
		      i++;
		      j++;
		    }
                }
	      else if(newbin2[i] > tmpnewbin[j])
                {
		     //cout << "je suis en 2 et k = " << k << endl;
		     cout << "je suis en 2 et k = " << k << " et i = " << i << " et newbin2 = " << newbin2[i] << " et newbin3 = " << newbin3[k] << " et j = " << j << " et tmpnewbin = " << tmpnewbin[j] << endl;
		  //newbin[k] = newbin2[i];
		  j++;
		  k--;
                }
	      else if(newbin2[i] < tmpnewbin[j])
                {
		     //cout << "je suis en 3 et k = " << k << endl;
		     cout << "je suis en 3 et k = " << k << " et i = " << i << " et newbin2 = " << newbin2[i] << " et newbin3 = " << newbin3[k] << " et j = " << j << " et tmpnewbin = " << tmpnewbin[j] << endl;
		  //newbin[k] = tmpnewbin[j];
		  i++;
		  k--;
                }
            }
        }
      cout << "res pour bin T = " << ibinT << " tmpinewbin = " << tmpinewbin << endl;
      inewbin2 = inewbin3;
      //newbin2 = newbin3;
      delete [] newbin2;
      newbin2 = new Double_t[fNEnergyBins+1];
      for (int lol=0; lol<inewbin3+1; lol++)
        {
	  //inewbin2++;
	  newbin2[lol]=newbin3[lol];
          cout << "lol = " << lol << " bin = " << newbin3[lol] << " et bin 2 = " << newbin2[lol] << endl;
        }
      delete [] tmpnewbin;
    }
    inewbin = inewbin3;
    //newbin = newbin3;
    for (int lol=0; lol<inewbin3+1; lol++)
      {
        newbin[lol] = newbin3[lol];
	//inewbin++;
      }
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
  IactBinnedLivLkl* mylkl = dynamic_cast<IactBinnedLivLkl*>(minuit->GetObjectFit());
  Double_t g             = par[0];
  TObjArrayIter* iter    = (TObjArrayIter*) mylkl->GetSampleArray()->MakeIterator();
  PoissonLkl* sample;

  // -2 log-likelihood
  f = 0;

  while((sample=dynamic_cast<PoissonLkl*>(iter->Next())))
    {
      Double_t   w_i = sample->GetUnitsOfG();
      Double_t   g_i = g*w_i;

      if(w_i>0)
        f+=sample->MinimizeLkl(g_i,kTRUE,kFALSE);
    }
}		
