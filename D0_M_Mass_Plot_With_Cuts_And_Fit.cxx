#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
using namespace RooFit;

void D0_M_Mass_Plot_With_Cuts_And_Fit()
{

//Open file
TFile * file = TFile::Open("/nfs/lhcb/malexander01/charm/4pi/data/full-run-I/full_2011_data.root");

//Define tree
TTree* DecayTree = dynamic_cast<TTree*>(file->Get("Dst2010ToD0ToKpipipipiTuple/DecayTree"));

//Create a subset to cut down on processing time
DecayTree->SetEntries(10000);

//Declare observables (Var to plot & vars to cut on)
RooRealVar D0_M("D0_M","",1790,1940);
RooRealVar piplus1_P("piplus1_P","",0,120000);
RooRealVar piplus1_PT("piplus1_PT","",0,5000);
RooRealVar D0_L0HadronDecision_TOS("D0_L0HadronDecision_TOS","",0,1);
RooRealVar lab0_L0Global_TIS("lab0_L0Global_TIS","",0,1);
RooRealVar lab0_Hlt2CharmHadD02HHHHDecision_TOS("lab0_Hlt2CharmHadD02HHHHDecision_TOS","",0,1);
RooRealVar D0_Hlt1TrackAllL0Decision_TOS("D0_Hlt1TrackAllL0Decision_TOS","",0,1);

//Create dataset loading D0_M & Various Cut variables
RooDataSet * data1 = new RooDataSet("data1", "data1", DecayTree, RooArgSet(D0_M,piplus1_P, piplus1_PT, D0_L0HadronDecision_TOS, lab0_L0Global_TIS, lab0_Hlt2CharmHadD02HHHHDecision_TOS, D0_Hlt1TrackAllL0Decision_TOS));

cout << "Number of events in data set:" << endl;
data1->Print() ;

//Make cuts on data
TCut piplus1_P_Cut = "piplus1_P > 3000"; 
TCut piplus1_PT_Cut = "piplus1_PT > 350";
TCut Trigger_Selection_1 = "(D0_L0HadronDecision_TOS == 1) || (lab0_L0Global_TIS == 1)";
TCut Trigger_Selection_2 = "lab0_Hlt2CharmHadD02HHHHDecision_TOS == 1";
TCut Trigger_Selection_3 = "D0_Hlt1TrackAllL0Decision_TOS == 1";

RooDataSet* ReducedDataSet = (RooDataSet*)data1->reduce(piplus1_P_Cut);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(piplus1_PT_Cut);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(Trigger_Selection_1);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(Trigger_Selection_2);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(Trigger_Selection_3);


//Fitting D0 Mass
//Define shape parameters
RooRealVar peak_gaussian_mean("peak_gaussian_mean","mean of gaussian for signal peak",1865,1805,1925, "MeV/c^{2}");
RooRealVar peak_gaussian_sigma("peak_gaussian_sigma","width of gaussian for signal peak",2.5, 0.0, 20., "MeV/c^{2}");

//Define Gaussian PDF
RooGaussian peak_gaussian("peak_gaussian", "gaussian for signal peak", D0_M, peak_gaussian_mean, peak_gaussian_sigma);

//Define yield parameter
RooRealVar peak_yield("peak_yield", "yield of signal peak", 1000,0,1000000);


//Fitting linear background
//Define shape parameter
RooRealVar background_c1("background_c1","slope of background", 0, -10, 10, "MeV/c^{2}");

//Define linear PDF
RooPolynomial background_poly("background_poly", "linear function for background", D0_M, RooArgList(background_c1));

//Define yield parameter
RooRealVar background_yield("background_yield", "yield of background", 100, 0, 10000);



//Total PDF
RooArgList shapes;
RooArgList yields;
shapes.add(background_poly);
yields.add(background_yield);
shapes.add(peak_gaussian);
yields.add(peak_yield);

RooAddPdf totalPdf("totalPdf","sum of signal and background PDFs",shapes,yields);

totalPdf.fitTo(*ReducedDataSet, Extended() );


//RooPlot

TCanvas* can2 = new TCanvas("can2","Mass fits data",800,600);
Plot = D0_M.frame(100);
Plot->SetTitle("Plot");
Plot->GetXaxis()->SetTitle("m(D^{0}) (MeV/c^{2})");
ReducedDataSet->plotOn(Plot);
//Display fit parameters
totalPdf.paramOn(Plot, Format("NELU",AutoPrecision(2)),Layout(0.1,0.4,0.9));
//Plot background
totalPdf.plotOn(Plot, Components(background_poly), LineColor(kGreen));
//Plot total PDF
totalPdf.plotOn(Plot, LineColor(kRed));

Plot->Draw();

}



