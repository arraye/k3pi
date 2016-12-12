/* D0_M Mass fit*/

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

void Laura_and_Laurens_combined()
{

//****************************************************
//*************       Set up data        *************
//****************************************************

//Open file
TFile * file = TFile::Open("/nfs/lhcb/malexander01/charm/4pi/data/full-run-I/full_2011_data.root");

//Define tree
TTree* DecayTree = dynamic_cast<TTree*>(file->Get("Dst2010ToD0ToKpipipipiTuple/DecayTree"));

//Create a subset to cut down on processing time
DecayTree->SetEntries(100000);

//Declare observables (Var to plot & vars to cut on)
RooRealVar D0_M("D0_M","",1790,1940);
RooRealVar piplus1_P("piplus1_P","",0,120000);
RooRealVar piplus1_PT("piplus1_PT","",0,5000);
RooRealVar D0_L0HadronDecision_TOS("D0_L0HadronDecision_TOS","",0,1);
RooRealVar lab0_L0Global_TIS("lab0_L0Global_TIS","",0,1);
RooRealVar lab0_Hlt2CharmHadD02HHHHDecision_TOS("lab0_Hlt2CharmHadD02HHHHDecision_TOS","",0,1);
RooRealVar D0_Hlt1TrackAllL0Decision_TOS("D0_Hlt1TrackAllL0Decision_TOS","",0,1);
RooRealVar D0_ID("D0_ID","",-500,500);

//Create dataset loading D0_M & Various Cut variables
RooDataSet * FullDataSet = new RooDataSet("data1", "data1", DecayTree, RooArgSet(D0_M,piplus1_P, piplus1_PT, D0_L0HadronDecision_TOS, lab0_L0Global_TIS, lab0_Hlt2CharmHadD02HHHHDecision_TOS, D0_Hlt1TrackAllL0Decision_TOS,D0_ID));

cout << "Number of events in data set:" << endl;
FullDataSet->Print() ;

//Make cuts on data
TCut piplus1_P_Cut = "piplus1_P > 3000"; 
TCut piplus1_PT_Cut = "piplus1_PT > 350";
TCut Trigger_Selection_1 = "(D0_L0HadronDecision_TOS == 1) || (lab0_L0Global_TIS == 1)";
TCut Trigger_Selection_2 = "lab0_Hlt2CharmHadD02HHHHDecision_TOS == 1";
TCut Trigger_Selection_3 = "D0_Hlt1TrackAllL0Decision_TOS == 1";
TCut D0_Data = "D0_ID>0";
TCut D0_bar_Data = "D0_ID<0";

RooDataSet* ReducedDataSet = (RooDataSet*)FullDataSet->reduce(piplus1_P_Cut);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(piplus1_PT_Cut);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(Trigger_Selection_1);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(Trigger_Selection_2);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(Trigger_Selection_3);



//****************************************************
//*************        PDF Shapes        *************
//****************************************************

// Yields
RooRealVar sig_D0_M_Yield("sigD0_MYield","sig yield D0_M", 0.5*ReducedDataSet->numEntries(),0,2*ReducedDataSet->numEntries());
RooRealVar bkg_D0_M_Yield("bkgD0_MYield","bkg yield D0_M", 0.7*ReducedDataSet->numEntries(),0,2*ReducedDataSet->numEntries());

// Gaussian signal PDF
RooRealVar sigmean("sigmean","D^{0} mass",1864.84,1820,1910);
RooRealVar sigwidth("sigwidth","D^{0} width",2.5,0.,10.);
RooGaussian signal("signal","Signal PDF",D0_M,sigmean,sigwidth);

// CB signal PDF
RooRealVar CB_mean("CB_mean","D^{0} mass",1865,1820,1910);
RooRealVar CB_width("CB_width","D^{0} width",2.5,0.,30.);
RooRealVar CB_alpha("CB_alpha","CB alpha",2,0,10);
RooRealVar CB_n("CB_n","CB n",4,0.01,10);
RooCBShape CBall("CBall","Crystall Ball shape",D0_M,CB_mean,CB_width,CB_alpha,CB_n);

// Combined signal PDF (gauss + CB)
RooRealVar sigfrac("sigfrac","signal fraction",0.7,0,1);
RooAddPdf sigD0_MPdf("sigD0_MPdf","D0_M sig with yields",RooArgList(signal,CBall),RooArgList(sigfrac));

//Background PDF
RooRealVar bkgD0_M("bkgD0_M","slope of background", 0, 0, 10, "MeV/c^{2}");
RooPolynomial bkgD0_MPdf("background_poly", "poly function for background", D0_M, RooArgList(bkgD0_M));

//Total PDF
RooAddPdf totalPdf("totalPdf","D0_M total",RooArgList(sigD0_MPdf,bkgD0_MPdf),RooArgList(sig_D0_M_Yield,bkg_D0_M_Yield));
totalPdf.fitTo(*ReducedDataSet, Extended() );


//****************************************************
//*************       Canvas Plot        *************
//****************************************************

//Create canvas
TCanvas* D0_M_Canvas = new TCanvas("can2","Mass fits data",800,600);

// Define legend + pads
TLegend *leg = new TLegend(0.63,0.58,0.90,0.90);
TLine *l=new TLine(1820.,0.0,1910.,0.0); //horizontal line for a pull plot
TPad *pad1 = new TPad("pad1","pad1",0,0.33,1,1);
TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.33);
pad1->SetBottomMargin(-0.005);
pad1->SetBorderMode(0);
pad2->SetTopMargin(-0.02);
pad2->SetBottomMargin(0.1);
pad2->SetBorderMode(0);
pad1->Draw();
pad2->Draw();

//Create D0_M Plot
D0_M_Frame = D0_M.frame(100);
D0_M_Frame->SetTitle("Plot");
D0_M_Frame->GetXaxis()->SetTitle("m(D^{0}) (MeV/c^{2})");
ReducedDataSet->plotOn(D0_M_Frame);
totalPdf.paramOn(D0_M_Frame, Format("NELU",AutoPrecision(2)),Layout(0.1,0.4,0.9)); //Display fit parameters
totalPdf.plotOn(D0_M_Frame, Components(bkgD0_MPdf), LineColor(kGreen),LineStyle(kDashed)); //Plot background
totalPdf.plotOn(D0_M_Frame, Components(sigD0_MPdf), LineColor(kBlue),LineStyle(kDashed)); //Plot signal
totalPdf.plotOn(D0_M_Frame, LineColor(kRed)); //Plot total PDF on frame D0_M_Frame

//Plot D0_M_Frame on pad1
pad1->cd(); //Makes pad1 the active pad
D0_M_Frame->Draw();
//Set up legend
leg->SetFillColor(kWhite);
leg->SetTextSize(0.04);
leg->SetTextColor(kBlack);
leg->SetHeader("D^{0} #rightarrow K^{+}#pi^{-}#pi^{+}#pi^{-}");
leg->AddEntry("h_data", "m(D^{0}) data", "lep");
leg->AddEntry("totalPdf_Norm[D0_M]","Total Fitting function PDF", "l");
leg->AddEntry("totalPdf_Norm[D0_M]_Comp[bkgD0_MPdf]","Background PDF","l");
leg->AddEntry("totalPdf_Norm[D0_M]_Comp[sigD0_MPdf]","Signal PDF","l");
leg->Draw();

//Create residuals plot and plot on pad2
RooHist* pull_hist=D0_M_Frame->pullHist(); //residHist but normalised wrt error
RooPlot* pull_frame = D0_M.frame(Title("Pull Distribution"));
pull_frame->addPlotable(pull_hist,"P");
pad2->cd();
pull_frame->Draw();
l->SetLineColor(kRed);
l->Draw("same");

D0_M_Canvas->cd();
}
