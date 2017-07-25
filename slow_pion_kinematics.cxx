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

using namespace RooFit;

int slow_pion_kinematics(int entries=0) 
{
  //Open ntuple
  TFile * file = TFile::Open("/nfs/lhcb/malexander01/charm/4pi/data/full-run-I/full_2011_data.root");

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
  Double_t px,py,pz;
  Int_t qs;

  //Point branch to objects - allows for reference to
  DecayTree->SetBranchAddress("piplus1_ID",&qs);
  DecayTree->SetBranchAddress("piplus1_PX",&px);
  DecayTree->SetBranchAddress("piplus1_PY",&py);
  DecayTree->SetBranchAddress("piplus1_PZ",&pz);


  //Plot histogram px vs py (for developing)
  Int_t nentries = (Int_t)DecayTree->GetEntries();
  TH2F *hpxpy = new TH2F("hpxpy","py vs px",100,-3000,3000,100,-3000,3000);
  TH1F *k_dist = new TH1F("k_dist","distribution of variable k",100,-0.001,0.02);

  TH1F *theta_x_dist = new TH1F("theta_x_dist","distribution of variable theta_x",100,-0.4,0.4);

  TH2F *kinematics = new TH2F("kinematics","qs_theta_x vs k",100,-200,200,100,0,0.02);

  double k;
  double theta_x;
  double qs_theta_x;

  //Loop through DecayTree and fill histograms
  for (Int_t i=0;i<nentries;i++) {
     k=0;
     theta_x = 0;
  
     DecayTree->GetEntry(i);
     //hpxpy->Fill(piplus1_PX,piplus1_PY);
     hpxpy->Fill(px,py);
     //k = (piplus1_PX*piplus1_PX + piplus1_PY*piplus1_PY);
     k = 1/sqrt((px*px + py*py));
     k_dist->Fill(k);
     //k_dist->Draw();
     theta_x = TMath::ATan(px/pz);
     theta_x_dist->Fill(theta_x);
  
     qs_theta_x = qs * theta_x;
     kinematics->Fill(qs_theta_x,k);
}

kinematics->Draw();
  
  




}
