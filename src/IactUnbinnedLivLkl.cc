/* ======================================================================== *\
!
!   Author(s): Daniel Kerszberg         10/2019 <mailto:dkerszbegr@ifae.es>
! 
\* ======================================================================== */

//////////////////////////////////////////////////////////////////////////////
//
// IactUnbinnedLivLkl
//
//////////////////////////////////////////////////////////////////////////////

// include c++ needed classes

// include Root needed classes
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPRegexp.h"

// include gLike needed classes
#include "IactUnbinnedLivLkl.h"
#include "IactEventListIrf.h"

#include <iomanip>

ClassImp(IactUnbinnedLivLkl);

using namespace std;

Double_t E_ebl[50] = {2.999999999999999889e-02,
3.454190000000000038e-02, 
3.977130000000000248e-02, 
4.579249999999999987e-02, 
5.272530000000000272e-02, 
6.070770000000000333e-02, 
6.989860000000000517e-02, 
8.048089999999999411e-02, 
9.266530000000000600e-02, 
1.066939999999999972e-01, 
1.228469999999999979e-01, 
1.414459999999999884e-01, 
1.628600000000000048e-01, 
1.875169999999999892e-01, 
2.159059999999999868e-01, 
2.485930000000000084e-01, 
2.862290000000000112e-01, 
3.295620000000000216e-01, 
3.794569999999999887e-01, 
4.369049999999999878e-01, 
5.030499999999999972e-01, 
5.792089999999999739e-01, 
6.668990000000000196e-01, 
7.678639999999999910e-01, 
8.841160000000000130e-01, 
1.017970000000000041e+00, 
1.172080000000000011e+00, 
1.349529999999999896e+00, 
1.553840000000000110e+00, 
1.789090000000000069e+00, 
2.059950000000000170e+00, 
2.371809999999999974e+00, 
2.730890000000000040e+00, 
3.144340000000000135e+00, 
3.620379999999999932e+00, 
4.168490000000000251e+00, 
4.799579999999999735e+00, 
5.526209999999999845e+00, 
6.362849999999999895e+00, 
7.326159999999999783e+00, 
8.435309999999999420e+00, 
9.712369999999999948e+00, 
1.118280000000000030e+01, 
1.287579999999999991e+01, 
1.482510000000000083e+01, 
1.706960000000000122e+01, 
1.965390000000000015e+01, 
2.262940000000000040e+01, 
2.605529999999999902e+01, 
3.000000000000000000e+01}; 

Double_t tau_ebl[50] = {2.709469999999999847e-03,
 6.039420000000000205e-03,
 1.148869999999999926e-02,
 1.982659999999999978e-02,
 3.105509999999999868e-02,
 4.613769999999999677e-02,
 6.490089999999999748e-02,
 8.931840000000000612e-02,
 1.207069999999999949e-01,
 1.635959999999999914e-01,
 2.227780000000000038e-01,
 3.027799999999999936e-01,
 4.113240000000000229e-01,
 5.510509999999999575e-01,
 7.336340000000000083e-01,
 9.559499999999999664e-01,
 1.233570000000000055e+00,
 1.556470000000000020e+00,
 1.933249999999999913e+00,
 2.358589999999999964e+00,
 2.817969999999999864e+00,
 3.304910000000000014e+00,
 3.795739999999999892e+00,
 4.285989999999999966e+00,
 4.743430000000000035e+00,
 5.174170000000000158e+00,
 5.553410000000000402e+00,
 5.898010000000000197e+00,
 6.219140000000000335e+00,
 6.547839999999999883e+00,
 6.915230000000000210e+00,
 7.335709999999999731e+00,
 7.839819999999999567e+00,
 8.407019999999999271e+00,
 9.089290000000000092e+00,
 9.866509999999999891e+00,
 1.084919999999999973e+01,
 1.207450000000000045e+01,
 1.375270000000000081e+01,
 1.614150000000000063e+01,
 1.972660000000000124e+01,
 2.530010000000000048e+01,
 3.358080000000000354e+01,
 4.639229999999999876e+01,
 6.437099999999999511e+01,
 9.087210000000000321e+01,
 1.262480000000000047e+02,
 1.740180000000000007e+02,
 2.347439999999999998e+02,
 3.097239999999999895e+02};

TGraph *gr_ebl = new TGraph(50,E_ebl,tau_ebl);

// class name and title
static const TString  gName            = "IactUnbinnedLivLkl";
static const TString  gTitle           = "Iact Unbinned Likelihood for LIV";

// List of free parameters.
static const Int_t    gNPars           = 1;                      // Number of free+nuisance parameters
static const Char_t*  gParName[gNPars] = {"eta"};                // Name of parameters
static const Int_t    gNBins           = 100;                    // default number of histograms for dN/dE plots

static const Int_t    gNFineLEBins       = 1000;                   // default number of fine bins for internal histos
static const Double_t gFineLEMin       = TMath::Log10(10);       // default minimum log(energy[GeV]) for internal histos
static const Double_t gFineLEMax       = TMath::Log10(1000000);   // default maximum log(energy[GeV]) for internal histos
static const Int_t    gNFineTBins       = 1000;                   // default number of fine bins for internal histos
static const Float_t  gFineTMin  = 1e01;                   // [s] default value of minimum arrival time
static const Float_t  gFineTMax  = 1e03;                   // [s] default value of maximum arrival time
static const Double_t gCenterBin       = 0.5;                    // decide which value represents bin in histogram (= 0 for lower bin edge, 0.5 for the middle, 1 for the right edge)

// static functions (for internal processing of input data)
static Int_t  SmearHistogram(TH2D* sp,TH2D* smsp,TGraph* grreso,TGraph* grbias);
static Int_t  SmearHistogram(TH2D* sp,TH2D* smsp,TH2F* mm);
static Int_t copyBinByBin(TH2D* ih,TH2D* oh,Double_t scale=0,Bool_t isDiff=kTRUE);

// -2logL function for minuit
void unbinnedLivLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// static minuit object for passing class info into the -2logL function
static TMinuit* minuit = NULL;

//////////////////////////////////////////////////////////////////////////////
//
// String constructor
//
IactUnbinnedLivLkl::IactUnbinnedLivLkl(TString inputString) :
  Lkl(gNPars,inputString,gName,gTitle), Iact1dUnbinnedLkl(inputString),
  fNFineLEBins(gNFineLEBins), fFineLEMin(gFineLEMin), fFineLEMax(gFineLEMax),
  fNFineTBins(gNFineTBins), fFineTMin(gFineTMin), fFineTMax(gFineTMax),
  fOnSampleEnergy(NULL), fOnSampleTime(NULL), fOffSampleTime(NULL),
  fHdNdESignal(NULL), fHdNdEpSignal(NULL), fHdNdEpBkg(NULL)
{
  if(InterpretInputString(inputString))
    cout << "IactUnbinnedLivLkl::IactUnbinnedLivLkl Warning: there were problems interpreting the input string" << endl;      
}

//////////////////////////////////////////////////////////////////
//
// Interpret the input configuration string
//
// Return 0 in case checks are Ok, 1 otherwise
//
Int_t IactUnbinnedLivLkl::InterpretInputString(TString inputString)
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
      cout << "IactUnbinnedLivLkl::InterpretInputString Warning: no IactEventListIrf object in file " << inputfileName << endl;
    }
  else
    {
      cout << "IactUnbinnedLivLkl::InterpretInputString Warning: THERE is IactEventListIrf object in file " << inputfileName << endl;

  TF1 *f2 = new TF1("f2", "[0]*TMath::Exp(-0.5*((x-[1])/[2])**2.)/([2]*TMath::Sqrt(2*TMath::Pi()))", 0., 100.);
  f2->SetParameters(1.,50.,10.);

      // extract data

      // extract info from file 
      fTMin = 0.;//dataSet->GetEpmin();
      fTMax = 10.;//dataSet->GetEpmax();

      Double_t eventOnT,eventOffT;
      dataSet->SetOnBranchAddress("t",&eventOnT);
      dataSet->SetOffBranchAddress("t",&eventOffT);

      fOnSampleTime  = new Double_t[GetNon()];
      fOffSampleTime = new Double_t[GetNoff()];

      for(Int_t i=0;i<GetNon();i++)
        {
          dataSet->GetOnEntry(i);
	  if(i==0) fTMin = (eventOnT-58497.)*86400.;
	  if(i==(GetNon()-1)) fTMax = (eventOnT-58497.)*86400.;
	  cout << setprecision(20) << " on " << i << " t = " << eventOnT << "in days or " << eventOnT*86400. << " in sec" << endl;
          fOnSampleTime[i] = (eventOnT-58497.)*24.*60.*60. - fTMin + 62.1;
	  if(fOnSampleTime[i]==-1) fOnSampleTime[i] = (i+1)*(90./GetNon());
	  //if(fOnSampleTime[i]==-1) fOnSampleTime[i] = f2->GetRandom();
        }
      for(Int_t i=0;i<GetNoff();i++)
        {
          dataSet->GetOffEntry(i);
	  cout << setprecision(20) << " off " << i << " t = " << eventOffT << endl;
          fOffSampleTime[i] = (eventOffT-58497.)*24.*60.*60. - fTMin + 62.1;
	  if(fOffSampleTime[i]==-1) fOffSampleTime[i] = (i+1.5)*(90./GetNon());
	  //if(fOffSampleTime[i]==-1) fOffSampleTime[i] = f2->GetRandom();
        }

      fTMin       = fOnSampleTime[0];
      fTMax       = fOnSampleTime[GetNon()-1];
    }

      cout << "fTmin = " << fTMin << endl;
      cout << "fTmax = " << fTMax << endl;
      cout << "fEmin = " << GetEmin() << endl;
      cout << "fEmax = " << GetEmax() << endl;
      //BuildAndBinOnOffHistos();
  return 0;
}

////////////////////////////////////////////////////////////////
// 
// Destructor
//
IactUnbinnedLivLkl::~IactUnbinnedLivLkl()
{
  if(fOnSampleTime)        delete [] fOnSampleTime;
  if(fOffSampleTime)       delete [] fOffSampleTime;
}

////////////////////////////////////////////////////////////////
//
// Link fMinuit to the static global variable minuit (otherwise 
// it won't work because of TMinuit subtleties)
//
void IactUnbinnedLivLkl::SetMinuitLink()
{
  minuit = fMinuit; // copy the pointer so that it can be used by FCN
}

////////////////////////////////////////////////////////////////
//
// Pass the -2logL function to fMinuit,
// initialize the parameters (free+nuisance)
//
void IactUnbinnedLivLkl::SetFunctionAndPars(Double_t ginit)
{
  fMinuit->SetFCN(unbinnedLivLkl);  
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
Int_t IactUnbinnedLivLkl::MakeChecks()
{
  if(IsChecked()) return 0;

  // add your checks here and try to mend whatever needs to be mended

  // Check the IactUnbinnedLivLkl specific part
  /////////////////////////////////////////
  // get the dN/dE' histograms for On and Off and check binning
  /*if(BuildAndBinOnOffHistos())
    {
      cout << "IactUnbinnedLivLkl::MakeChecks (" << GetName() << ") Warning: problems building On and/or Off histos!" << endl;
      return 1;
    }

cout << "Coucou" << endl;*/

  SetChecked();
  return 0;
}

////////////////////////////////////////////////////////////////
//
// Check that the needed histograms are present and that we are ready 
// for calling minimization 
//
// If needed, convolute dNdESignal*Aeff with Eres and Ebias
// to get dNdEpSignal
//
// If needed, convolute dNdESignal*AeffOff with Eres and Ebias
// to get dNdEpSignalOff
//
// Return 0 in case of success, 1 otherwise
//
Int_t IactUnbinnedLivLkl::CheckHistograms(Bool_t checkdNdEpBkg)
{
  // if fHdNdEpSignal is missing, try to construct it from fHdNdESignal, fHAeff fGEreso and fGEbias
  if(!fHdNdEpSignal && (fHdNdESignal && GetHAeff() && ((GetGEreso() && GetGEbias()) || GetMigMatrix())))
    {
      if(GetMigMatrix())
        cout << "IactUnbinnedLivLkl::CheckHistograms Message: will create fHdNdEpSignal from fHdNdESignal, fHAeff & fMigMatrix... " << flush;
      else
        cout << "Iact1dUnbinnedLkl::CheckHistograms Message: will create fHdNdEpSignal from fHdNdESignal, fHAeff, fGEreso & fGEbias... " << flush;

      // multiply dNdESignal times Aeff
      TH2D* hdNdESignalAeff = new TH2D("hdNdESignalAeff","Product of dN/dE for Signal and Aeff",fNFineTBins,fFineTMin,fFineTMax,fNFineLEBins,fFineLEMin,fFineLEMax);
      hdNdESignalAeff->SetDirectory(0);
      //hdNdESignalAeff->Multiply(GetHAeff(),fHdNdESignal);
      TH1D* haeff = (TH1D*)GetHAeff()->Clone();
      for(Int_t ibin=1;ibin<=fNFineTBins;ibin++)
        {
          for(Int_t jbin=1;jbin<=fNFineLEBins;jbin++)
            {
              Double_t lowedge = fHdNdESignal->GetYaxis()->GetBinLowEdge(jbin);
              Int_t bin = haeff->FindBin(lowedge,0,0);
              hdNdESignalAeff->SetBinContent(ibin,jbin,fHdNdESignal->GetBinContent(ibin,jbin)*GetHAeff()->GetBinContent(bin));
            }
        }
      //hdNdESignalAeff->SaveAs("./bkg_Accept.root");

      // create fHdNdEpSignal   
      fHdNdEpSignal         = new TH2D("fHdNdEpSignal","dN/dE' for Signal",fNFineTBins,fFineTMin,fFineTMax,fNFineLEBins,fFineLEMin,fFineLEMax);
      fHdNdEpSignal->SetDirectory(0);

      // smear hdNdESignalAeff
      TH2F* MigMatrix = (TH2F*)GetMigMatrix()->Clone();
      if(GetMigMatrix())
        {
          if(SmearHistogram(hdNdESignalAeff,fHdNdEpSignal,MigMatrix))
          //if(SmearHistogram(hdNdESignalAeff,fHdNdEpSignal,GetGEreso(),GetGEbias()))
            return 1;
        }
      /*else
        {
          if(SmearHistogram(hdNdESignalAeff,fHdNdEpSignal,fGEreso,fGEbias))
            return 1;
        }*/

      cout << "Done! " << endl;
      // clean
      delete hdNdESignalAeff;
    }

  if(!fHdNdEpBkg && (fHdNdEBkg && GetHAeff() && ((GetGEreso() && GetGEbias()) || GetMigMatrix())))
    {
      if(GetMigMatrix())
        cout << "IactUnbinnedLivLkl::CheckHistograms Message: will create fHdNdEpBkg from fHdNdEBkg, fHAeff & fMigMatrix... " << flush;
      else
        cout << "Iact1dUnbinnedLkl::CheckHistograms Message: will create fHdNdEpBkg from fHdNdEBkg, fHAeff, fGEreso & fGEbias... " << flush;

      // multiply dNdEBkg times Aeff
      TH2D* hdNdEBkgAeff = new TH2D("hdNdEBkgAeff","Product of dN/dE for Bkg and Aeff",fNFineTBins,fFineTMin,fFineTMax,fNFineLEBins,fFineLEMin,fFineLEMax);
      hdNdEBkgAeff->SetDirectory(0);
      TH1D* haeff = (TH1D*)GetHAeff()->Clone();
      for(Int_t ibin=1;ibin<=fNFineTBins;ibin++)
        {
          for(Int_t jbin=1;jbin<=fNFineLEBins;jbin++)
            {
              Double_t lowedge = fHdNdEBkg->GetYaxis()->GetBinLowEdge(jbin);
              Int_t bin = haeff->FindBin(lowedge,0,0);
              hdNdEBkgAeff->SetBinContent(ibin,jbin,fHdNdEBkg->GetBinContent(ibin,jbin)/*GetHAeff()->GetBinContent(bin)*/);
            }
        }
      //hdNdEBkgAeff->SaveAs("./bkggg_Accept.root");

      // create fHdNdEpBkg   
      fHdNdEpBkg         = new TH2D("fHdNdEpBkg","dN/dE' for Bkg",fNFineTBins,fFineTMin,fFineTMax,fNFineLEBins,fFineLEMin,fFineLEMax);
      fHdNdEpBkg->SetDirectory(0);

      // smear hdNdEBkgAeff
      TH2F* MigMatrix = (TH2F*)GetMigMatrix()->Clone();
      if(GetMigMatrix())
        {
          if(SmearHistogram(hdNdEBkgAeff,fHdNdEpBkg,MigMatrix))
          //if(SmearHistogram(hdNdESignalAeff,fHdNdEpSignal,GetGEreso(),GetGEbias()))
            return 1;
        }
      /*else
        {
          if(SmearHistogram(hdNdESignalAeff,fHdNdEpSignal,fGEreso,fGEbias))
            return 1;
        }*/

      cout << "Done! " << endl;
      // clean
      delete hdNdEBkgAeff;
    }

  // normalize unnormalized histos
  NormalizedNdEHisto(fHdNdEpSignal);
  //fHdNdEpSignal->SaveAs("./bkg_Accept_Smear.root");
  NormalizedNdEHisto(fHdNdEpBkg);
  //fHdNdEpBkg->SaveAs("./bkggggg_Accept_Smear.root");

  if(!fHdNdEpSignal)
    cout << "IactUnbinnedLivLkl::CheckHistograms Warning: fHdNdEpSignal histogram missing!!" << endl;
  if(!fOnSampleTime)
    cout << "IactUnbinnedLivLkl::CheckHistograms Warning: fOnSampleTime histogram missing!!" << endl;
  if(!fOffSampleTime)
    cout << "IactUnbinnedLivLkl::CheckHistograms Warning: fOffSampleTime histogram missing!!" << endl;

  return 1;
}

////////////////////////////////////////////////////////////////
// smear a spectral shape <sp> (in true energy) using 
// the energy dispersion function from migration matrix <mm>
// and put the result in histogram <smsp> (in measured energy)
// <mm> is computed as N_ij/(N_j*DeltaE_j)
// with N_ij number of events passing all analysis cuts, and
// with true energy in DeltaE_j and recontructed energy in DeltaE_i;
// N_j total number of events passing all analysis cuts with true energy in DeltaE_j;
// and DeltaE_j the size [GeV] of the DeltaE_j energy bin
// 
Int_t SmearHistogram(TH2D* sp,TH2D* smsp,TH2F* mm)
{
  // checks
  if(!sp || !smsp || !mm)
    {
      cout << "SmearHistogram Warning: missing histos" << endl;
      return 1;
    }

  Int_t nbinste  = sp->GetNbinsY();     // number of bins in input histo
  TH2D* provsmsp = new TH2D("provsmsp","Provisonal smeared histo",sp->GetNbinsX(),sp->GetXaxis()->GetXmin(),sp->GetXaxis()->GetXmax(),mm->GetXaxis()->GetNbins(),mm->GetXaxis()->GetXmin(),mm->GetXaxis()->GetXmax());
  provsmsp->SetDirectory(0);

  // do the convolution of sp with mm and store result in smsp
  for(Int_t ibin=0;ibin<sp->GetXaxis()->GetNbins();ibin++)
    {
      for(Int_t ibinte=0;ibinte<nbinste;ibinte++)
        {
          Double_t et = sp->GetYaxis()->GetBinCenter(ibinte+1); // log of true energy

          // bin size in true energy
          Double_t minbinte = TMath::Power(10,sp->GetYaxis()->GetBinLowEdge(ibinte+1));
          Double_t maxbinte = TMath::Power(10,sp->GetYaxis()->GetBinLowEdge(ibinte+1)+sp->GetYaxis()->GetBinWidth(ibinte+1));
          Double_t det      = maxbinte-minbinte;  // bin size

          // do the convolution of sp with the matrix row j=ibinte and fill
          // a histogram with same binning as matrix column i=ibinme
          for(Int_t iprovbinme=0;iprovbinme<mm->GetXaxis()->GetNbins();iprovbinme++)
            {
              Double_t em       = provsmsp->GetYaxis()->GetBinCenter(iprovbinme+1);
              Int_t    gidbin   = mm->FindBin(em,et);
              Double_t smfactor = mm->GetBinContent(gidbin);
              provsmsp->SetBinContent(ibin,iprovbinme+1,provsmsp->GetBinContent(ibin,iprovbinme+1)+sp->GetBinContent(ibin,ibinte+1)*det*smfactor);
            }
        }
    }

  // transfer provsmsp to smsp (with different binning, in general...)
  if(copyBinByBin(provsmsp,smsp))
    return 1;

  delete provsmsp;
  return 0;
}

//////////////////////////////////////////////////////////////////
//
// Copy the contents of ih into oh bin by bin
// the x axis of ih and oh must be in log-scale and oh in even binning (not necessary for ih)
// Coarse bins in ih are copied into finner bins in oh, but
// the shape (i.e. the discontinuities between bin limits) is preserved 
// the units of the x-axis of ih are those of oh times 10^scale
// if isDiff=kTRUE the input histogram is differential (default)
// if isDiff=kFALSE the input histogram is bin-integrated
// the output histogram is ALWAYS differential
Int_t copyBinByBin(TH2D* ih,TH2D* oh,Double_t scale,Bool_t isDiff)
{
  // input histogram binning
  Double_t imine   = ih->GetYaxis()->GetXmin()+scale; // minimum log(E) in input histo
  Double_t imaxe   = ih->GetYaxis()->GetXmax()+scale; // maximum log(E) in input histo
  Int_t    inbinse = ih->GetNbinsY();

  // output histogram binning
  Double_t omine   = oh->GetYaxis()->GetXmin(); // minimum log(E) in output histo
  Double_t omaxe   = oh->GetYaxis()->GetXmax(); // maximum log(E) in output histo
  Int_t    onbinse = oh->GetNbinsY();
  Double_t ode     = (omaxe-omine)/onbinse;

  // copy values
  for(Int_t jbin=0;jbin<oh->GetNbinsX();jbin++)
    {
  for(Int_t ibin=0;ibin<onbinse;ibin++)
    {
      Double_t etest = omine+ode*(ibin+gCenterBin);

      // copy bin by bin, zero outside limits
      if(etest<imine || etest>imaxe)
        oh->SetBinContent(jbin,ibin+1,0);
      else
        {
          // corresponding bin in ih histo
          Int_t etestbin;
          for(etestbin=0;etestbin<inbinse;etestbin++)
            if(ih->GetYaxis()->GetBinLowEdge(etestbin+1)+scale<etest && ih->GetYaxis()->GetBinLowEdge(etestbin+1)+ih->GetYaxis()->GetBinWidth(etestbin+1)+scale>etest)
              break;
          if(etestbin>=inbinse)
            {
              cout << "copyBinByBin Warning: the two histograms must have very different energy ranges" << endl;
              return 1;
            }
          Float_t dE = 1;
          if(!isDiff)
            {
              Double_t leminbin = ih->GetYaxis()->GetBinLowEdge(etestbin+1)+scale;
              Double_t lemaxbin = ih->GetYaxis()->GetBinLowEdge(etestbin+1)+ih->GetYaxis()->GetBinWidth(etestbin+1)+scale;
              dE                = TMath::Power(10,lemaxbin)-TMath::Power(10,leminbin);

            }
          oh->SetBinContent(jbin,ibin+1,ih->GetBinContent(jbin,etestbin+1)/dE);
        }
    }
    }
  return 0;
}

////////////////////////////////////////////////////////////////
//
// Normalize dN/dE' histos for signal and background and dN/dE for
// signal. 
// Save integral in bin 0.
// If bin 0 contains a non-zero value, do not normalize
//
// Return 0 in case of success
//        1 otherwise
//
Int_t IactUnbinnedLivLkl::NormalizedNdEHisto(TH2D* histo)
{
  // basic check
  if(!histo) return 1;

  // do not normalize if it's already normalized
  if(histo->GetBinContent(0)>0) return 0;

  // normalize and keep normalization in bin 0
  //Double_t intSignal = IntegrateLogE(histo,TMath::Log10(fEpmin),TMath::Log10(fEpmax));
  Double_t intSignal = histo->Integral();
  histo->Scale(1./intSignal);
  histo->SetBinContent(0,intSignal);

  return 0;
}

////////////////////////////////////////////////////////////////
//
// Feed the JointLkl with one PoissonLkl per E' bin 
//
// Return 0 in case of success, 1 otherwise
//
/*Int_t IactUnbinnedLivLkl::ConfigureJointLkl()
{

  if(!fHNOn || !fHNOff)
    {
      cout << "IactUnbinnedLivLkl::ConfigureJointLkl (" << GetName() << ") Warning: On and/or Off histograms do no exist" << endl;
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
      cout << "IactUnbinnedLivLkl::ConfigureJointLkl (" << GetName() << ") Warning: total sum of weights should be 1! (It's " << totalw << "), too many bins maybe??" << endl;
      return 1;
    }

  return 0;
}*/

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
// of IactUnbinnedLivLkl.
//
TH2D* IactUnbinnedLivLkl::GetHdNdEpOn(Bool_t isDifferential,Int_t nbinsE,Int_t nbinsT) const
{
      // we need a positive number of bins
      if(nbinsE<=0) nbinsE = gNBins; //gNEnergyBins;
      if(nbinsT<=0) nbinsT = gNBins; //gNTimeBins;

cout << "TEST1" << endl;
cout << nbinsT << " " << fTMin << " " << fTMax << " " << nbinsE << " " << GetEmin() << " " << GetEmax() << endl;

      // create histo
      TH2D* h = new TH2D("dNdEpOn","dN/dE' vs t for On events",nbinsT,0.,1500.,nbinsE,TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));
      //TH2D* h = new TH2D("dNdEpOn","dN/dE' vs t for On events",nbinsT,fTMin,fTMax,nbinsE,TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));
      //TH2F* h = new TH2F("dNdEpOn","dN/dE' vs t for On events",nbinsT,fTmin,fTmax,nbinsE,TMath::Log10(150.),TMath::Log10(2500.));
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
// of IactUnbinnedLivLkl.
//
TH2D* IactUnbinnedLivLkl::GetHdNdEpOff(Bool_t isDifferential,Int_t nbinsE,Int_t nbinsT) const
{
      // we need a positive number of bins
      if(nbinsE<=0) nbinsE = gNBins; //gNEnergyBins;
      if(nbinsT<=0) nbinsT = gNBins; //gNTimeBins;
    
cout << "TEST2" << endl;

      // create histo
      TH2D* h = new TH2D("dNdEpOff","dN/dE' vs t for Off events",nbinsT,0.,1500.,nbinsE,TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));
      //TH2D* h = new TH2D("dNdEpOff","dN/dE' vs t for Off events",nbinsT,fTMin,fTMax,nbinsE,TMath::Log10(GetEmin()),TMath::Log10(GetEmax()));
      //TH2F* h = new TH2F("dNdEpOff","dN/dE' vs t for Off events",nbinsT,fTmin,fTmax,nbinsE,TMath::Log10(150.),TMath::Log10(2500.));
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
/*Int_t IactUnbinnedLivLkl::SimulateDataSamples(UInt_t seed,Float_t meanG)
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

  // build pdf*/
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
/*
}*/

///////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Histograms hdNdEpBkg and hdNdEpSignalOff will contain same entries as fHdNdEpBkg and fHdNdEpSignalOff,
// respectively, but with normalization fluctuating according to the uncertainty in tau
// Note that fHdNdEpBkg is the expected distribution of background events in the On region
// and fHdNdEpSignalOff the expected distribution of signal events in the total Off region
// (total meaning that if tau=3 the effective area to consder is that of the three subregions)
//
/*Int_t IactUnbinnedLivLkl::GetRealBkgAndGoffHistos(TRandom3* rdm,TH2F*& hdNdEpBkg,TH2F*& hdNdEpSignalOff) 
{
  // create new histos with contents of the existing ones
  //if(fHdNdEpBkg) hdNdEpBkg = new TH1F(*fHdNdEpBkg);
  if(fHdNdEpSignal) hdNdEpBkg = new TH2F(*fHdNdEpSignal);
  else
    {
      cout << "IactUnbinnedLivLkl::GetRealBkgAndGoffHistos Warning: fHdNdEpBkg histo does not exist" << endl;
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
      cout << "IactUnbinnedLivLkl::GetRealBkgAndGoffHistos Error: negative or null tau value, we cannot work like that!" << endl;
      return 1;
    }
  
  return 0;
}*/

//////////////////////////////////////////////////////////////////
// 
// Recreate a fresh new version of fHdNdESignal
//
// Return 0 in case of success
//        1 if file is not found
//
Int_t IactUnbinnedLivLkl::ResetdNdESignal()
{
  // Delete existing fHdNdESignal and create empty one
  if(fHdNdESignal)
    delete fHdNdESignal;

  // Create histo
  fHdNdESignal = new TH2D("fHdNdESignal","dN/dE vs t for signal events",fNFineTBins,fFineTMin,fFineTMax,fNFineLEBins,fFineLEMin,fFineLEMax);
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

  /*if(fHdNdEpSignalOff)
  {
    delete fHdNdEpSignalOff;
    fHdNdEpSignalOff=NULL;
  }*/
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
Int_t IactUnbinnedLivLkl::SetdNdESignalFunction(TString function,Float_t p0,Float_t p1,Float_t p2,Float_t p3,Float_t p4,Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{
  // Delete existing fHdNdESignal and create empty one
  ResetdNdESignal();

  // call to add the function
  AdddNdESignalFunction(function,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9);

  // exit
  return 0;
}

Bool_t firstTime=kTRUE;
Double_t integralFirstTime=0.;

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
Int_t IactUnbinnedLivLkl::AdddNdESignalFunction(TString function,Float_t p0,Float_t p1,Float_t p2,Float_t p3,Float_t p4,Float_t p5,Float_t p6,Float_t p7,Float_t p8,Float_t p9)
{

  fFineTMin=p2;//-0.017*p4*p1 ;
  fFineTMax=p3;//-0.017*p4*p0;
  cout << "tmin = " << fFineTMin << " and tmax = " << fFineTMax << endl;
  fFineLEMin=TMath::Log10(p0);
  fFineLEMax=TMath::Log10(p1);

  //cout << "n_t = " << fNFineTBins << " tmin = " << fFineTMin << "tmax = " << fFineTMax << "n_e = " << fNFineLEBins << " emin = " << fFineLEMin << " emax = " << fFineLEMax << endl;

  //fHdNdESignal = new TH2D("fHdNdESignal","dN/dE vs t for signal events",fNFineTBins,0.,1500.,fNFineLEBins,fFineLEMin,fFineLEMax);
  // good !!! fHdNdESignal = new TH2D("fHdNdESignal","dN/dE vs t for signal events",fNFineTBins,fFineTMin,fFineTMax,fNFineLEBins,fFineLEMin,fFineLEMax);
  //fHdNdESignal = new TH2D("fHdNdESignal","dN/dE vs t for signal events",fNFineTBins,1.,1600.,fNFineLEBins,fFineLEMin,fFineLEMax);
  fHdNdESignal = new TH2D("fHdNdESignal","dN/dE vs t for signal events",fNFineTBins,fFineTMin,fFineTMax,fNFineLEBins,fFineLEMin,fFineLEMax);

  //fHdNdEpBkg = new TH2D("fHdNdEpBkg","dN/dE vs t for signal events",fNFineTBins,1.,1600.,fNFineLEBins,fFineLEMin,fFineLEMax);
  fHdNdEBkg = new TH2D("fHdNdEpBkg","dN/dE vs t for signal events",fNFineTBins,fFineTMin,fFineTMax,fNFineLEBins,fFineLEMin,fFineLEMax);
  /*const Float_t*      offSample       = GetOffSample();
  const Double_t*     offSampleTime   = GetOffSampleTime();
  UInt_t              Noff            = GetNoff();
  for(int i=0; i < Noff; i++)
    {
      fHdNdEpBkg->SetBinContent(fHdNdEpBkg->FindBin(offSampleTime[i],offSample[i]),1.);
      //fHdNdEpBkg->SetBinContent(ibin+1,jbin+1,1.);
    }*/

  Double_t Emin=p0;
  Double_t Emax=p1;
  Double_t Tmin=p2;
  Double_t Tmax=p3;
  Double_t eta =p4;
  Double_t alpha =p5;
  Double_t beta =p6;
  Double_t T_0 =p7;
  Double_t T_1 =p8;

  Double_t log10Emin = TMath::Log10(Emin);
  Double_t log10Emax = TMath::Log10(Emax);

      Int_t jbinmin = fNFineLEBins*(log10Emin-fFineLEMin)/(fFineLEMax-fFineLEMin);
      Int_t jbinmax = fNFineLEBins*(log10Emax-fFineLEMin)/(fFineLEMax-fFineLEMin);
      Double_t realEmin;//  = TMath::Power(10,fHdNdESignal->GetYaxis()->GetBinLowEdge(jbinmin+1));
      Double_t realEmax;//  = TMath::Power(10,fHdNdESignal->GetYaxis()->GetBinLowEdge(jbinmax+1)+fHdNdESignal->GetYaxis()->GetBinWidth(jbinmax+1));
      Double_t dE;//        = realEmax-realEmin;
      Double_t E;//        = realEmax-realEmin;

      Int_t ibinmin = fNFineTBins*(Tmin-fFineTMin)/(fFineTMax-fFineTMin);
      Int_t ibinmax = fNFineTBins*(Tmax-fFineTMin)/(fFineTMax-fFineTMin);
      Double_t realTmin;//  = fHdNdESignal->GetXaxis()->GetBinLowEdge(ibinmin+1);
      Double_t realTmax;//  = fHdNdESignal->GetXaxis()->GetBinLowEdge(ibinmax+1)+fHdNdESignal->GetYaxis()->GetBinWidth(ibinmax+1);
      Double_t dt;//        = realTmax-realTmin;
      Double_t t;//        = realTmax-realTmin;

      Double_t constant = TMath::Power(T_1,7.3-1.3*TMath::Log(T_1))*TMath::Power(T_1,beta);
      //cout << "ibinmin = " << ibinmin << " and ibinmax = " << ibinmax << endl;

      for(Int_t ibin=ibinmin;ibin<=ibinmax;ibin++)
        {
              realTmin = fHdNdESignal->GetXaxis()->GetBinLowEdge(ibin+1);
	      realTmax = fHdNdESignal->GetXaxis()->GetBinLowEdge(ibin+1)+fHdNdESignal->GetXaxis()->GetBinWidth(ibin+1);
	      dt = realTmax-realTmin;
	      t = (realTmax+realTmin)/2.;
	      //cout << "t = " << t << endl;
	      Double_t constant_1 = TMath::Power(realTmin,7.3-1.3*TMath::Log(realTmin))*TMath::Power(realTmin,beta);
	      Double_t constant_2 = TMath::Power(realTmax,7.3-1.3*TMath::Log(realTmax))*TMath::Power(realTmax,beta);
          for(Int_t jbin=jbinmin;jbin<=jbinmax;jbin++)
            {
              realEmin = TMath::Power(10,fHdNdESignal->GetYaxis()->GetBinLowEdge(jbin+1));
	      realEmax = TMath::Power(10,fHdNdESignal->GetYaxis()->GetBinLowEdge(jbin+1)+fHdNdESignal->GetYaxis()->GetBinWidth(jbin+1));
	      dE = realEmax-realEmin;
	      E = (realEmax+realEmin)/2.;
	      //Double_t dE_model = (TMath::Power(realEmax,-5.43) - TMath::Power(realEmin,-5.43))/-5.43;
	      Double_t dE_model = ((TMath::Power(realEmax,-2.2) - TMath::Power(realEmin,-2.2))/-2.)*TMath::Exp(gr_ebl->Eval(E/1000.));
	      //cout << "E = " << E << " and absorption = " << gr_ebl->Eval(E/1000.) << endl;
              fHdNdEBkg->SetBinContent(ibin+1,jbin+1,TMath::Power(E,-1.51));
              //fHdNdEBkg->SetBinContent(ibin+1,jbin+1,t*E);
              //fHdNdEBkg->SetBinContent(ibin+1,jbin+1,dt*dE_model);

	      //if(t-(2.5e-5)*eta*E*E < T_1)
	      if(t-0.017*eta*E < T_1)
                {
                  //if(t-(2.5e-5)*eta*E*E < T_0)
                  if(t-0.017*eta*E < T_0)
		    {
		      //cout << " case 1 t = " << t << " and E = " << E << endl;
	              fHdNdESignal->SetBinContent(ibin+1,jbin+1,/*fHdNdESignal->GetBinContent(ibin+1,jbin+1)*/+0.);
		    }
		  else
	            {
			    if(TMath::Power(t,7.3-1.3*TMath::Log(t)) <0.) cout << "part 1 = " << TMath::Power(E,-alpha) << " part 2 = " << TMath::Power(t,7.3-1.3*TMath::Log(t)) << endl;
			    //fHdNdESignal->SetBinContent(ibin+1,jbin+1,/*fHdNdESignal->GetBinContent(ibin+1,jbin+1)*/ + TMath::Power(E,-alpha)*TMath::Power(t,7.3-1.3*TMath::Log(t)) );
			    // good!! fHdNdESignal->SetBinContent(ibin+1,jbin+1,/*fHdNdESignal->GetBinContent(ibin+1,jbin+1)*/ + TMath::Power(E,-alpha)*TMath::Power(t,7.3-1.3*TMath::Log(t)) );
			    //fHdNdESignal->SetBinContent(ibin+1,jbin+1,/*fHdNdESignal->GetBinContent(ibin+1,jbin+1)*/ + (TMath::Power(realEmax,-alpha+1)/(-alpha+1) + TMath::Power(realEmin,-alpha+1)/(-alpha+1))*(realTmax-realTmin)*(constant_2+constant_1)/2.);
		            fHdNdESignal->SetBinContent(ibin+1,jbin+1,/*fHdNdESignal->GetBinContent(ibin+1,jbin+1)*/ + 1.*dE_model/*(1./dE)*(1./dt)*/);
			    //cout << "dE=" << dE << " and dT=" << dt << " and dE*dT=" << (1./dE)*(1./dt) << endl;
		    }
                }
	      else
               {
		      //cout << " case 3 t = " << t << " and E = " << E << " and val = " << (TMath::Power(realEmax,-alpha+1)/(-alpha+1) + TMath::Power(realEmin,-alpha+1)/(-alpha+1)) << " val2 = " << constant << " val3 = " << (TMath::Power(realTmax,-beta+1)/(-beta+1) + TMath::Power(realTmin,-beta+1)/(-beta+1)) << endl;
		 //if(realTmin < 0.1) realTmin=1.;
		 if(realTmin < 0.1) fHdNdESignal->SetBinContent(ibin+1,jbin+1,/*fHdNdESignal->GetBinContent(ibin+1,jbin+1)*/ + 0.);
	 	 //else fHdNdESignal->SetBinContent(ibin+1,jbin+1,/*fHdNdESignal->GetBinContent(ibin+1,jbin+1)*/ + (TMath::Power(realEmax,-alpha+1)/(-alpha+1) + TMath::Power(realEmin,-alpha+1)/(-alpha+1))*constant*(TMath::Power(realTmax,-beta+1)/(-beta+1) + TMath::Power(realTmin,-beta+1)/(-beta+1)) );
	 	 // good!! fHdNdESignal->SetBinContent(ibin+1,jbin+1,/*fHdNdESignal->GetBinContent(ibin+1,jbin+1)*/ + (TMath::Power(E,-alpha))*constant*(TMath::Power(t,-beta)) );
		 fHdNdESignal->SetBinContent(ibin+1,jbin+1,/*fHdNdESignal->GetBinContent(ibin+1,jbin+1)*/ + 1.*dE_model/*(1./dE)*(1./dt)*/);
			    //cout << "dE=" << dE << " and dT=" << dt << " and dE*dT=" << (1./dE)*(1./dt) << endl;
               }
            }
        }

cout << "int = " << fHdNdESignal->Integral() << endl;
if (fHdNdESignal->Integral()>0.) fHdNdESignal->Scale(1./fHdNdESignal->Integral());
//if (fHdNdESignal->Integral()>0. && !firstTime) {fHdNdESignal->Scale(1./integralFirstTime); cout << " NOT FIIIIIRST" << endl;}
//else if (fHdNdESignal->Integral()>0.) {fHdNdESignal->Scale(1./fHdNdESignal->Integral()); firstTime=kFALSE; integralFirstTime=fHdNdESignal->Integral(); cout << "first time!" << endl;}
//fHdNdESignal->SaveAs("bkg10.root");
cout << "int = " << fHdNdESignal->Integral() << endl;
if (fHdNdEBkg->Integral()>0.) fHdNdEBkg->Scale(1./fHdNdEBkg->Integral());

if(p9>0.1)
{
Double_t Erdm, Trdm;
fOnSampleEnergy  = new Double_t[GetNon()];
for (int test=0;test<GetNon();test++)
  {
    fHdNdESignal->GetRandom2(Trdm,Erdm);
    fOnSampleEnergy[test]=Erdm;
    fOnSampleTime[test]=Trdm/*+0.017*eta*Erdm*/;
    cout << "i = " << test << " and Trdm = " << Trdm << " and Erdm = " << Erdm << endl;
  }
}
//Iact1dUnbinnedLkl::GetHdNdEpBkg()->SaveAs("bkg_th1.root");

CheckHistograms(kTRUE);

/*  if(function=="constantLIV")
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
   
      fHdNdESignal = new TH2D();//"fHdNdESignal","fHdNdESignal",fNTimeBins,(1)*(90./nevents),(nevents)*(90./nevents),fNEnergyBins,GetEmin()/1000.,GetEmax()/1000.);
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
      //fHNOn->Draw("COLZ");
      c1->cd(4);
      //fHNOff->Draw("COLZ");
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
  //fHNOn = new TH2I("fHNOn","fHNOn",6,Tmin,Tmax,6,Emin,Emax);
  //fHNOff = new TH2I("fHNOff","fHNOff",6,Tmin,Tmax,6,Emin,Emax);
  fHdNdESignal = new TH2D("fHdNdESignal","fHdNdESignal",6,Tmin,Tmax,6,Emin,Emax);
  //fNFineTBins=6, fFineTMin=Tmin,fFineTMax=Tmax,fNFineLEBins=6,fFineLEMin=Emin,fFineLEMax=Emax;
  //fNTimeBins=6,fNEnergyBins=6;
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
      //fHNOn->Fill(t_shifted,E2);
      //fHNOff->Fill(t,E);
   }

  TCanvas *c1 = new TCanvas("c1","c1",900,300);
  c1->Divide(3,3);
  c1->cd(1);
  f1->Draw();
  c1->cd(2);
  f2->Draw();
  c1->cd(3);
  //fHNOn->Draw("COLZ");
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
    }*/
  return 0;
}

////////////////////////////////////////////////////////////////////////
//
// Likelihood function (-2logL) 
// To be minimized by TMinuit
// Free parameters:
// par[0] = eta 
//
void unbinnedLivLkl(Int_t &fpar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // to avoid warnings only
  gin[0]*=1;
  fpar*=1;
  iflag*=1;

Double_t stepLog = TMath::Exp((TMath::Log(50000) - TMath::Log(5.))/(30.));
  Double_t boundaries[30+1];
  boundaries[0] = 5.;
  for (int i = 1 ; i < 30 ; i++){
    boundaries[i] = boundaries[i-1]*stepLog;
  }
  boundaries[30] = 50000.;

  for (int i=0; i < 30; i++) {
	  Double_t energ = (TMath::Log10(boundaries[i])+TMath::Log10(boundaries[i+1]))/2.;
//ST0302
//cout << (7.43831e-01)+(-2.74945e-01)*energ+(3.82371e-02)*energ*energ "   ";
//ST0303
//cout << (1.16033e+00)+(-5.49942e-01)*energ+(6.77807e-02)*energ*energ "   " ;
//ST0306
//cout << (1.21740e+00)+(-6.29997e-01)*energ+(8.74824e-02)*energ*energ "   " ;
//ST0307
//cout << (8.15351e-01 )+(-3.59342e-01)*energ+(4.20904e-02)*energ*energ "   " ;
//ST0310
//cout << (1.68368e+00)+(-7.95249e-01)*energ+(9.82831e-02)*energ*energ << "   " ;
//ST0311
cout << /*energ << "   " <<*/ (7.43831e-01)+(-2.74945e-01)*energ+(3.82371e-02)*energ*energ << ",   " ;
  }


  Double_t x[101], y[101];
  //Double_t eta_inject = 1.;

  IactUnbinnedLivLkl* mylkl           = dynamic_cast<IactUnbinnedLivLkl*>(minuit->GetObjectFit());
  Double_t old_lkl[mylkl->GetNon()];
  Double_t new_lkl[mylkl->GetNon()];

  const Float_t*      onSample        = mylkl->GetOnSample();
  //const Double_t*      onSample        = mylkl->GetOnSampleEnergy();
  const Float_t*      offSample       = mylkl->GetOffSample();
  const Double_t*     onSampleTime    = mylkl->GetOnSampleTime();
  const Double_t*     offSampleTime   = mylkl->GetOffSampleTime();
  UInt_t              Non             = mylkl->GetNon();
  UInt_t              Noff            = mylkl->GetNoff();
  //Float_t             tau             = mylkl->GetTau();
  Float_t             tau             = 3.;
  Float_t             dTau            = mylkl->GetDTau();

  Double_t realEmin=2000.;
  Double_t realEmax=0.;
  for(int i=0; i < Non; i++)
    {
      if(onSample[i]>realEmax) realEmax = onSample[i];
      if(onSample[i]<realEmin) realEmin = onSample[i];
    }

    cout << "Emin = " << realEmin << " Emax = " << realEmax << endl;

  //mylkl->SetdNdESignalFunction("",250.,2000.,0.,1200.,eta_inject,2.4,1.5,1.,30.,0.111);
  mylkl->SetdNdESignalFunction("",300.,TMath::Power(10.,realEmax),onSampleTime[0],onSampleTime[Non-2],0.,2.4,1.5,1.,30.,0.); // skipped < 300 GeV --> Non-2
  //mylkl->SetdNdESignalFunction("",TMath::Power(10.,realEmin),TMath::Power(10.,realEmax),onSampleTime[0],onSampleTime[Non-2],0.,2.4,1.5,1.,30.,0.); // skipped < 300 GeV --> Non-2
  // for all events mylkl->SetdNdESignalFunction("",TMath::Power(10.,realEmin),TMath::Power(10.,realEmax),onSampleTime[0],onSampleTime[Non-1],0.,2.4,1.5,1.,30.,0.);
  //mylkl->SetdNdESignalFunction("",250.,2000.,62.,1225.,0.,2.4,1.5,1.,30.,0.);	
  //for(Double_t eta=-2.; eta<2.; eta+=0.1)
  //for(Double_t eta=-2.; eta<2.0; eta+=0.1)
  for(Int_t eta=0; eta<101; eta++)
  {
  // get internal object, histos, values, etc
  cout << "par0 = " << eta/10.-4. << endl;
  //cout << "par0 = " << par[0] << endl;

  mylkl->SetdNdESignalFunction("",300.,TMath::Power(10.,realEmax),onSampleTime[0],onSampleTime[Non-2],eta/10.-4.,2.4,1.5,1.,30.,0.); // events > 300 GeV --> Non-2
  //mylkl->SetdNdESignalFunction("",TMath::Power(10.,realEmin),TMath::Power(10.,realEmax),onSampleTime[0],onSampleTime[Non-2],eta/10.-4.,2.4,1.5,1.,30.,0.); // events > 300 GeV --> Non-2
  // all events mylkl->SetdNdESignalFunction("",TMath::Power(10.,realEmin),TMath::Power(10.,realEmax),onSampleTime[0],onSampleTime[Non-1],eta/10.-4.,2.4,1.5,1.,30.,0.);
  //mylkl->SetdNdESignalFunction("",225.,2000.,62.,1225.,eta/10.-4.,2.4,1.5,1.,30.,0.);
  //mylkl->SetdNdESignalFunction("",200.,2000.,50.,1500.,par[0],2.4,1.5,0.,30.);

  const TH2D*         hdNdEpSignal    = mylkl->GetHdNdEpSignal();
  //hdNdEpSignal->SaveAs("./hdNdEpSignal.root");
  //const TH1F*         hdNdEpBkg       = mylkl->GetHdNdEpBkg();
  //hdNdEpBkg->SaveAs("bkg_th1.root");
  const TH2D*         hdNdEpBkg       = mylkl->GetHdNdEpBkg();

  /*const Int_t      nbinsT           = hdNdEpBkg->GetNbinsX();
  const Double_t   tmin            = hdNdEpBkg->GetXaxis()->GetXmin();
  const Double_t   tmax            = hdNdEpBkg->GetXaxis()->GetXmax();
  const Int_t      nbins          = hdNdEpBkg->GetNbinsY();
  const Double_t   xmin            = hdNdEpBkg->GetYaxis()->GetXmin();
  const Double_t   xmax            = hdNdEpBkg->GetYaxis()->GetXmax();*/
  const Int_t      nbinsT           = hdNdEpSignal->GetNbinsX();
  const Double_t   tmin            = hdNdEpSignal->GetXaxis()->GetXmin();
  const Double_t   tmax            = hdNdEpSignal->GetXaxis()->GetXmax();
  const Int_t      nbins          = hdNdEpSignal->GetNbinsY();
  const Double_t   xmin            = hdNdEpSignal->GetYaxis()->GetXmin();
  const Double_t   xmax            = hdNdEpSignal->GetYaxis()->GetXmax();

  // Estimated number of background events in signal and background regions
  Double_t g       = mylkl->GetdNdEpSignalIntegral();///1000000.;//par[0]; //GetG();
  g       = Non; //g*Non;
  //Double_t g       = hdNdEpSignal->GetBinContent(0);//par[0]; //GetG();
  //Double_t b       = (Non + Noff - (1.+tau)*g + TMath::Sqrt(TMath::Power(Non + Noff - (1.+tau)*g,2) + 4.*(1.+tau)*Noff*g))/(2.*(1.+tau));   //par[1]; 
  //Double_t b       = Noff;
  Double_t b       = Noff/tau;
  cout << "g = " << g << " and b = " << b << endl;
  //Double_t tauest  = par[2];
  //Double_t boff    = b*tauest;
  Double_t boff    = b*tau;
  //Double_t fnorm   = g+b+boff;
  Double_t fnorm   = 726;

  // sum signal and background contributions and normalize resulting pdf (On + Off)
  TH2D* hdNdEpOn  = new TH2D("hdNdEpOn", "On  event rate vs E' vs t",nbinsT,tmin,tmax,nbins,xmin,xmax);
  hdNdEpOn->Reset();
  //hdNdEpOn->Add(hdNdEpSignal,hdNdEpSignal,g/(g+b),g/(g+b));
  //hdNdEpOn->Add(hdNdEpSignal,hdNdEpSignal,1./2.,1./2.);
  //hdNdEpOn->Add(hdNdEpSignal,(Non-b)/Non);
  //hdNdEpOn->Add(hdNdEpSignal,hdNdEpBkg,(Non-b),b);
  hdNdEpOn->Add(hdNdEpSignal,hdNdEpBkg,(726-119/3.),119/3.);

  // normalize
  if(fnorm>0)
    hdNdEpOn->Scale(1./(fnorm));
  else
    mylkl->NormalizedNdEHisto(hdNdEpOn);

  //hdNdEpOn->SaveAs("./hdNdEpOn.root");
  TH2D* hdNdEpOff = new TH2D("hdNdEpOff","Off event rate vs E' vs t",nbinsT,tmin,tmax,nbins,xmin,xmax);
  hdNdEpOff->Reset();
  hdNdEpOff->Add(hdNdEpBkg,b/Non);

  /*// normalize
  if(fnorm>0)
    hdNdEpOff->Scale(1./fnorm);
  else
    mylkl->NormalizedNdEHisto(hdNdEpOff);*/

  // -2 log-likelihood
  f = 0;

  Int_t skipped = 0;
  // On events
  for(ULong_t ievent=0; ievent<Non; ievent++)
    {
	    if (onSample[ievent] < 2.477) {skipped++; cout << "SKIPPED " << std::setprecision(6) << TMath::Power(10.,onSample[ievent]) << " " << onSampleTime[ievent] << endl; continue;}
	    //if (onSample[ievent] > 3.) continue;
      cout << std::setprecision(6) << TMath::Power(10.,onSample[ievent]) << " " << onSampleTime[ievent] << endl; 
      Float_t val = hdNdEpOn->GetBinContent(hdNdEpOn->FindBin(onSampleTime[ievent]/*-0.017*eta*TMath::Power(10.,onSample[ievent])*/,onSample[ievent]));// + hdNdEpOff->GetBinContent(hdNdEpOff->FindBin(onSampleTime[ievent]/*-0.017*eta*TMath::Power(10.,onSample[ievent])*/,onSample[ievent]));
      //Float_t val = hdNdEpOn->GetBinContent(hdNdEpOn->FindBin(onSampleTime[ievent]-0.000025*eta*TMath::Power(10.,onSample[ievent])*TMath::Power(10.,onSample[ievent]),onSample[ievent]));// + hdNdEpOff->GetBinContent(hdNdEpOff->FindBin(onSampleTime[ievent]/*-0.017*eta*TMath::Power(10.,onSample[ievent])*/,onSample[ievent]));
      //cout << "lkl val = " << val << endl;
      if(val>0)
	{
          f += -2*TMath::Log(val);
	  if(eta==0) old_lkl[ievent]=-2*TMath::Log(val);
	  else {
            new_lkl[ievent]=-2*TMath::Log(val);
	    if(TMath::Abs(old_lkl[ievent]-new_lkl[ievent])<0.001) old_lkl[ievent]=-2*TMath::Log(val);
	    else {
	      cout << "Why lower? old val = " << old_lkl[ievent] << " and new = " << new_lkl[ievent] << " for E = " << onSample[ievent] << " and T = " << onSampleTime[ievent] << endl;
	      old_lkl[ievent]=-2*TMath::Log(val);
	    }
	  }
	  //cout << "Why THIS? i = " << ievent << " bin = " << hdNdEpOn->FindBin(onSample[ievent],onSampleTime[ievent]) << " E = " << onSample[ievent] << " and T = " << onSampleTime[ievent] << " and log = " << val << endl;
	}
      else
	{
	  cout << "Why 0? i = " << ievent << " bin = " << hdNdEpOn->FindBin(onSampleTime[ievent],onSample[ievent]) << " E = " << onSample[ievent] << " and T = " << onSampleTime[ievent] << " and log = " << val << endl;
        f += 100.;
        //f += 5000.;
        //f += 1e99;
	}
    }

cout << "skipped = " << skipped << endl;

  // Off events
  /*for(ULong_t ievent=0; ievent<Noff; ievent++)
    {
      Float_t val = hdNdEpOff->GetBinContent(hdNdEpOff->FindBin(offSampleTime[ievent],offSample[ievent]));
      if(val>0)
        f += -2*TMath::Log(val);
      else
	{
	  cout << "Why 0? i = " << ievent << " bin = " << hdNdEpOff->FindBin(offSampleTime[ievent],offSample[ievent]) << " E = " << offSample[ievent] << " and T = " << offSampleTime[ievent] << " and log = " << val << endl;
        f += 100.;
        //f += 0;
        //f += 1e99;
	}
    }*/

  // nuisance tau
    f += -2*TMath::Log(TMath::Gaus(1.47,1.51,0.04,kTRUE));

  // nuisance tau
  //if(dTau>0)
    //f+=-2*TMath::Log(TMath::Gaus(tauest, tau, dTau, kTRUE));

  // tot Nevts and nuisance Noff
  /*if(g+b>0)
    f += -2*TMath::Log(TMath::Poisson(Non,g+b));
  else
    f += 0;
    //f += 1e99;

  if(boff>0)
    f += -2*TMath::Log(TMath::Poisson(Noff,boff));
  else
    f += 0;
    //f += 1e99;*/

  cout << "-2loglkl = " << f << endl;

  delete hdNdEpOn;
  delete hdNdEpOff;
  x[eta]=eta/10.-4.;
  y[eta]=f;
  }
  Double_t min=y[0];
  Double_t eta_min=99.;
  for(int i=0;i<101;i++) {if (y[i]<min) {min=y[i]; eta_min=x[i];} cout << "i= " << i << " x = " << x[i] << " y = " << y[i] << endl;}
  for(int i=0;i<101;i++) {y[i]=y[i]-min;}

TLatex latex;

   TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
   TGraph *gr = new TGraph(101,x,y);
   gr->SetTitle(";#eta_{1};-2log(#lambda)");

/*Double_t actual_value = eta_min;
Double_t testing= gr->Eval(actual_value-1)-gr->Eval(eta_min)-2.71;
while(TMath::Abs(testing)>0.1)
{
//	cout << "actual value = " << actual_value << " and testing = " << testing;
if(testing<0.) actual_value+=-0.01;
else actual_value+=+0.01;
//	cout << " AFTER actual value = " << actual_value << " and testing = " << testing << endl;;
testing= gr->Eval(actual_value)-gr->Eval(eta_min)-2.71;
}

cout << "lower limit = " << actual_value << endl;*/

/*Double_t actual_value_max = eta_min;
Double_t testing_max= gr->Eval(actual_value_max+1)-gr->Eval(eta_min)-2.71;
while(TMath::Abs(testing_max)>0.1)
{
	cout << "actual value = " << actual_value_max << " and testing = " << testing_max;
if(testing_max<0.) actual_value_max+=+0.01;
else actual_value_max+=-0.01;
	cout << " AFTER actual value = " << actual_value_max << " and testing = " << testing_max << endl;
testing_max= gr->Eval(actual_value_max)-gr->Eval(eta_min)-2.71;
}

cout << "upper limit = " << actual_value_max << endl;*/

   //c1->SetLogy();
   //gr->GetHistogram()->SetMaximum(min+16.);
   //gr->GetHistogram()->SetMinimum(min-1.);
   gr->GetHistogram()->SetMaximum(16.);
   gr->GetHistogram()->SetMinimum(-1.);
   gr->Draw("AC*");
   //latex.DrawLatex(-1.,min+90.,Form("#eta_{inject} = %.1f",eta_inject));
   //latex.DrawLatex(-1.,min+80.,Form("#eta_{rec} = %.1f",eta_min));
   //latex.DrawLatex(-1.,min+80.,Form("#eta_{rec} = %.1f^{+%.1f}_{-%.1f}",eta_min,actual_value_max-eta_min,eta_min-actual_value));
   //latex.DrawLatex(-1.,min+80.,Form("#eta_{rec} = %.1f^{+%.1f}",eta_min,actual_value_max-eta_min));
   //c1->SaveAs("./plot_lkl.pdf");
}
