// Macro testLIVSearches.C
// Author: D. Kerszberg
// Date: March 2019
// For beginners, to understand the basic usage of the IactBinnedLivLkl class
// IMPORTANT NOTE: for some still-to-be-understood "feature", macros
// containing Lkl-based objects MUST be run in compiled mode, i.e.
// run this macro with:
// .x testLIVSearches.C+

#include <iostream>
//#include "Iact1dUnbinnedLkl.h"
//#include "Iact1dBinnedLkl.h"
#include "IactBinnedLivLkl.h"
#include "JointLkl.h"
#include "TCanvas.h"

void testLIVSearches(Bool_t data, Float_t delay)
{

if(data)
{

  cout << "MOCK DATA TEST" << endl;

  // input data
  const Double_t logJ         = 19.;   // [GeV^2 cm^-5] log_10 of J-factor of the assumed DM source
  const Double_t DlogJ        = 0;     // [GeV^2 cm^-5] statistical error in log_10 of J-factor of the assumed DM source
  //const Double_t mass         = 1000.; // [GeV] mass of the DM particle
  //const TString  dNdEFileName = TString(Form("./DM/dNdE/Cirelli/dNdESignal_bb_%.1fmass.root",mass)); // dN/dE input file  
  const TString  inputFile1   = "./data/genericIact_dataIRF_01.root";  // input file with event list and their associated IRF
  const TString  inputFile2   = "./data/genericIact_dataIRF_02.root";  // input file with event list and their associated IRF
  const TString  inputFile3   = "./data/Output_glikeInputs-Bin1-Koji.root";
  //const TString  inputFile4   = "/home/dkerszberg/Softs/gLikeLiv/data/GRB190114C/Output_glikeInputs_W0and180_wobblepos0.root";
  const TString  inputFile4   = "/home/dkerszberg/Softs/gLikeLiv/data/GRB190114C/Output_glikeInputs_W0and180_wobblepos1.root";
  //const TString  inputFile4   = "/home/dkerszberg/Softs/gLikeLiv/data/Mrk421/Output_glikeInputs_W0and180_wobblepos0.root";
  //const TString  inputFile4   = "/home/dkerszberg/Softs/gLikeLiv/data/Mrk421/Output_glikeInputs_W0and180_wobblepos1.root";
  //const TString  inputFile4   = "/home/dkerszberg/Softs/gLikeLiv/data/Mrk421/Output_glikeInputs_W90and270_wobblepos0.root";
  //const TString  inputFile4   = "/home/dkerszberg/Softs/gLikeLiv/data/Mrk421/Output_glikeInputs_W90and270_wobblepos1.root";
  const Double_t errorDef     = 4;

  // create and configure an Iact1dUnbinnedLkl object for 1D unbinned likelihood analysis
  //Iact1dUnbinnedLkl* unbn = new Iact1dUnbinnedLkl(Form("logJ=%.2f DlogJ=0 inputfile=%s",logJ,inputFile1.Data()));
  //unbn->ReaddNdESignal(dNdEFileName);
  //unbn->SetDMAnnihilationUnitsForG(mass);    // set units for DM annihilation <sv>
  IactBinnedLivLkl* binLIV1 = new IactBinnedLivLkl(Form("logJ=%.2f DlogJ=0 inputfile=%s",logJ,inputFile4.Data()));
  binLIV1->AdddNdESignalFunction("constantLIV",delay,1.,1.);
  cout << "QUOI" << endl;
  binLIV1->ConfigureJointLkl();
  binLIV1->ComputeLklVsG();
  binLIV1->GetLklVsG()->Draw(); // plot the -2logL vs g curve
  //binLIV1->MakeChecks();
  binLIV1->PrintData();
  binLIV1->PrintOverview();
  //IactBinnedLivLkl* binLIV2 = new IactBinnedLivLkl(Form("logJ=%.2f DlogJ=0 inputfile=%s",logJ,inputFile2.Data()));
  //binLIV2->PrintData();
  //binLIV2->PrintOverview();
  //JointLkl* jointLkl = new JointLkl(Form("DlogJ=%.2f",DlogJ));
  //jointLkl->SetErrorDef(errorDef);
  //jointLkl->AddSample(binLIV1);
  //jointLkl->AddSample(binLIV2);
  //jointLkl->PrintData();
  //jointLkl->ConfigureJointLkl();
  //jointLkl->ComputeLklVsG();
  //jointLkl->GetLklVsG()->Draw(); // plot the -2logL vs g curve

  //binLIV1->ConfigureJointLkl();
  //binLIV1->ComputeLklVsG();

  //return;

  // create and configure an Iact1dBinnedLkl object for 1D unbinned likelihood analysis
  //Iact1dBinnedLkl* bn = new Iact1dBinnedLkl(Form("logJ=%.2f DlogJ=0 inputfile=%s",logJ,inputFile2.Data()));
  //bn->ReaddNdESignal(dNdEFileName);
  //bn->SetDMAnnihilationUnitsForG(mass);    // set units for DM annihilation <sv>

  // create and fill a JointLkl object for the combined analysis of both datasets
  //JointLkl* jointLkl = new JointLkl(Form("DlogJ=%.2f",DlogJ));
  //jointLkl->SetErrorDef(errorDef);
  //jointLkl->AddSample(unbn);
  //jointLkl->AddSample(bn);

  //IactBinnedLivLkl* truc = new IactBinnedLivLkl();
  //truc->SetUnitsOfG(1./(500.*25.*TMath::Log(1.5)));
  //truc->AdddNdESignalFunction("",delay,1.,1.,1.5,1.,100.,5.,90.,115.);
  //truc->BuildAndBinOnOffHistos();
  //truc->ConfigureJointLkl();
  //truc->ComputeLklVsG();


  double res, error;
  res=binLIV1->GetParVal(0);
  error = binLIV1->GetParErr(0);
  cout << "res = " << res << " error = " << error << endl;

      int order=1;
      double H_0 = 0.673/(9.777752e9*365.25*24*60*60);
      double z=0.42;
      //Double_t kl(Double_t z) {
      TF1 *h = new TF1("h", "pow((1+x),[0])/sqrt(0.6825+0.3175*pow(1+x,3))", 0, 6);
      h->SetParameter(0,order);
      //TF1 *h = new TF1("h", "pow((1+x),order)/sqrt(0.7+0.3*pow(1+x,3))", 0, 6);
      Double_t intz = h->Integral(0, z);
      cout << "intz = " << intz << endl;
      cout << "res = " << res << " error = " << error << endl;
      double EQG_linear = ((1.+order)/(2.*H_0))*(intz/TMath::Abs(res));
      cout << scientific << "EQG_linear = " << EQG_linear << endl;
  // print how JointLkl looks like
  //jointLkl->PrintData();

  // call minimization and print/plot results
  //jointLkl->ComputeLklVsG();
  //jointLkl->PrintOverview();     // print the details from the fit
  //TCanvas* c2 = new TCanvas("c2","",700,500);
  //jointLkl->GetLklVsG()->Draw(); // plot the -2logL vs g curve
  //truc->GetLklVsG()->Draw(); // plot the -2logL vs g curve
  ////truc->SetUnitsOfG(0.5*25.);
  //cout << truc->GetUnitsOfG() << endl;
  //cout << truc->GetdNdEpSignalIntegral() << endl;
}
  else
    {
      cout << "MC TEST" << endl;
      TH1D* resSimu = new TH1D("resSimu","resSimu",50,0.,1.);
      double res, error;
      for (int i = 0; i<1; i++)
        {
          IactBinnedLivLkl* truc = new IactBinnedLivLkl();
          truc->AdddNdESignalFunction("gaussLIV",delay,1.,1.,1.5,1.,100.,5.,90.,115.);
          truc->BuildAndBinOnOffHistos();
          truc->ConfigureJointLkl();
          truc->ComputeLklVsG();
          cout << "g for lkl 0 = " << truc->GetGForLkl(0.,kTRUE) << endl;
          cout << "par val 0 = " << truc->GetParVal(0) << endl;
          cout << "par err 0 = " << truc->GetParErr(0) << endl;
          resSimu->Fill(truc->GetParVal(0));
	  res=truc->GetParVal(0);
	  error = truc->GetParErr(0);
	  delete truc;
        }

      TCanvas* c2 = new TCanvas("c2","",700,500);
      resSimu->Draw();
      //truc->GetLklVsG()->Draw(); // plot the -2logL vs g curve
      //cout << "Units of G = " << truc->GetUnitsOfG() << endl;
      //cout << "Signal integral = " << truc->GetdNdEpSignalIntegral() << endl;

      int order=1;
      double H_0 = 0.673/(9.777752e9*365.25*24*60*60);
      double z=0.42;
      //Double_t kl(Double_t z) {
      TF1 *h = new TF1("h", "pow((1+x),[0])/sqrt(0.6825+0.3175*pow(1+x,3))", 0, 6);
      h->SetParameter(0,order);
      //TF1 *h = new TF1("h", "pow((1+x),order)/sqrt(0.7+0.3*pow(1+x,3))", 0, 6);
      Double_t intz = h->Integral(0, z);
      cout << "intz = " << intz << endl;
      cout << "res = " << res << " error = " << error << endl;
      //res += 2.*error;
      cout << "res = " << res << endl;
      delete h;
      //return res;
      //}

      double EQG_linear = ((1.+order)/(2.*H_0))*(intz/TMath::Abs(res));
      cout << "EQG_linear = " << EQG_linear << endl;
    }
}
