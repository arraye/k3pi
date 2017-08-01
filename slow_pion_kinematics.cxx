/********************************************************/
/***********  Slow Pion Kinematic Distributions *********/
/********************************************************/

/* D0_M Mass fit*/

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

// ROOT
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TAxis.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCut.h"
#include "TLine.h"
#include "TFrame.h"
#include "TCanvas.h"


// RooFit
/*
#include "RooPlot.h"
#include "RooHist.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooTruthModel.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooAbsData.h"
#include "RooBinning.h"
#include "RooCategory.h"
#include "RooCBShape.h"
#include "RooPolynomial.h"
#include "TFitResult.h"
#include "RooExponential.h"
#include "RooProdPdf.h"
#include "RooChebychev.h"
*/

//C++
#include <vector>
#include <sstream>
#include <iomanip> //for use of setprecision
#include <iostream>
#include <fstream>

//For sleep function
#include <stdlib.h>

#include <libgen.h>

using namespace RooFit;

int slow_pion_kinematics(int entries=0) 
{


//std::string rootfiles[4][2] = {{"md_11","data/run1/charm_magdown_2011.root"},{"mu_11","data/run1/charm_magup_2011.root"},{"md_12","data/run1/DVntuple119.root"},{"mu_12","data/run1/charm_magup_2012.root"}}

std::string rootfiles[1][2]={"full_run_2011","/nfs/lhcb/malexander01/charm/4pi/data/full-run-I/full_2011_data.root"};

for (int filen=0; filen<1; filen++) {

  //Open ntuple
  //TString filepath = "/nfs/lhcb/malexander01/charm/4pi/data/full-run-I/full_2011_data.root";
  TString filepath = rootfiles[filen][1];
  TFile * file = TFile::Open(filepath);



  //Define tree
  TTree* DecayTree = dynamic_cast<TTree*>(file->Get("Dst2010ToD0ToKpipipipiTuple/DecayTree"));

  //Might not be necessary - can take subset for debugging if slow
  if (entries > 0) {
    cout<<"Starting analysis with "<<entries<<" entries."<<endl;
    DecayTree->SetEntries(entries);
  }
  else {   cout<<"Starting analysis with all entries."<<endl; }

  //Obtain only necessary branches
  DecayTree->SetBranchStatus("*",0);
  DecayTree->SetBranchStatus("piplus1_P",1);
  DecayTree->SetBranchStatus("piplus1_PT",1);
  DecayTree->SetBranchStatus("piplus1_ID",1);
  DecayTree->SetBranchStatus("piplus1_PX",1);
  DecayTree->SetBranchStatus("piplus1_PY",1);
  DecayTree->SetBranchStatus("piplus1_PZ",1); 

  //Objects for holding each branch object
  Double_t px,py,pz; //Momenta of slow pion
  Int_t qs;  //Sign of slow pion

  //Point branch to objects - allows for reference to
  DecayTree->SetBranchAddress("piplus1_ID",&qs);
  DecayTree->SetBranchAddress("piplus1_PX",&px);
  DecayTree->SetBranchAddress("piplus1_PY",&py);
  DecayTree->SetBranchAddress("piplus1_PZ",&pz);


  //Plot histogram px vs py (for developing)
  Int_t nentries = (Int_t)DecayTree->GetEntries();
  TH2F *hpzpy = new TH2F("hpzpy","py vs pz",100,0,20000,100,-2000,2000);
  TH1F *k_dist = new TH1F("k_dist","distribution of variable k",100,-0.001,0.01);

  TH1F *theta_x_dist = new TH1F("theta_x_dist","distribution of variable theta_x",100,-0.4,0.4);

  TH2F *kinematicspos = new TH2F("kinematicspos","qs_theta_x vs k",100,-0.31,0.31,100,0,1);
  TH2F *kinematicsneg = new TH2F("kinematicsneg","qs_theta_x vs k",100,-0.31,0.31,100,0,1);
  //TH2F *kinematicsratio = new TH2F("kinematicsneg","qs_theta_x vs k",100,-0.31,0.31,100,0,1);

  TH3F *kinematics_ratio_3d = new TH3F("kinematics_ratio","qs_theta_x vs k",100,-0.31,0.31,100,0,1,100,0,2);
  

  double k;
  double theta_x;
  double qs_theta_x;
  int sign;

  //Debug!
  for(Int_t i=0;i<3;i++){
    DecayTree->GetEntry(i);
    cout<<"px: "<<px<<" py: "<<py<<" pz: "<<pz<<endl;
    }

  //Loop through DecayTree and fill histograms
  for (Int_t i=0;i<nentries;i++) {
     k=0;
     theta_x = 0;
  
     DecayTree->GetEntry(i);
     //hpxpy->Fill(piplus1_PX,piplus1_PY);
     //hpzpy->Fill(pz,py);

     //k (or C for curvature) = 1/sqrt(px^2+pz^2)
     //pz and px in MeV, convert to GeV (for compliance with D0->hh binned)
     //Add to k_dist histogram for each event
     k = (1/(sqrt(px/1000*px/1000 + pz/1000*pz/1000)));
     k_dist->Fill(k);

     //theta_x = arc
     theta_x = TMath::ATan(px/pz);
     theta_x_dist->Fill(theta_x);
  
     if (qs>0){sign=+1;} else{sign=-1;}
     //qs_theta_x = sign * theta_x;
     qs_theta_x = theta_x;
     //let's look at only +ve soft pion distr. first
     if (sign==-1){
     //cout<<"debug: sign is neg? sign = "<<sign<<endl;
     kinematicsneg->Fill(qs_theta_x,k);}
     else{
     kinematicspos->Fill(qs_theta_x,k);}

     //cout<<"sign: "<<sign<<" theta_x: "<<theta_x<<" qs_theta_x: "<<qs_theta_x<<endl;
}

TCanvas* kinposcanvas = new TCanvas("kinposcanvas","Positive soft pion kinematics canvas",900,600);
kinposcanvas->cd();
kinposcanvas->SetLogz();
kinematicspos->Draw("colz");

kinposcanvas->SaveAs((rootfiles[filen][0] + "_kinposcanvas.pdf").c_str());

TCanvas* kinnegcanvas = new TCanvas("kinnegcanvas","Positive soft pion kinematics canvas",900,600);
kinnegcanvas->cd();
kinematicsneg->Draw("colz");
kinnegcanvas->SaveAs((rootfiles[filen][0] + "_kinnegcanvas.pdf").c_str());


}
}
