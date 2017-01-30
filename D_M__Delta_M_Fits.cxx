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
#include "TMath.h"
using namespace RooFit;

void D_M__Delta_M_Fits()
{

//****************************************************
//*************       Set up data        *************
//****************************************************

//Open file
TFile * file = TFile::Open("/nfs/lhcb/malexander01/charm/4pi/data/full-run-I/full_2011_data.root");

//Define tree
TTree* DecayTree = dynamic_cast<TTree*>(file->Get("Dst2010ToD0ToKpipipipiTuple/DecayTree"));

//Create a subset to cut down on processing time for development
DecayTree->SetEntries(100000);

//SetBranchStatus - Avoids root processing timeout on entire fileset
DecayTree->SetBranchStatus("*",0);
DecayTree->SetBranchStatus("D0_M",1);
DecayTree->SetBranchStatus("D0_L0HadronDecision_TOS",1);
DecayTree->SetBranchStatus("lab0_L0Global_TIS",1);
DecayTree->SetBranchStatus("D0_Hlt1TrackAllL0Decision_TOS",1);
DecayTree->SetBranchStatus("lab0_Hlt2CharmHadD02HHHHDecision_TOS",1);
DecayTree->SetBranchStatus("piplus1_P",1);
DecayTree->SetBranchStatus("piplus1_PT",1);
DecayTree->SetBranchStatus("lab0_M",1); //D*(2010);
DecayTree->SetBranchStatus("D0_TAU",1);
DecayTree->SetBranchStatus("D0_ID",1);


//***************Declare observables (Var to plot & vars to cut on)***************
//D0_M
RooRealVar D0_M("D0_M","",1790,1940);
double D0_M_Min(1820),D0_M_Max(1910); //Maximum and minimum values of D0_M to consider

//Delta_M Min & Max
double Delta_M_Min(139.57), Delta_M_Max(155);

//D*(2010)
RooRealVar lab0_M("lab0_M","lab0_M",D0_M_Min+Delta_M_Min,D0_M_Max+Delta_M_Max);

//Cuts and ID variables
RooRealVar piplus1_P("piplus1_P","",0,120000);
RooRealVar piplus1_PT("piplus1_PT","",0,5000);
RooRealVar D0_L0HadronDecision_TOS("D0_L0HadronDecision_TOS","",0,1);
RooRealVar lab0_L0Global_TIS("lab0_L0Global_TIS","",0,1);
RooRealVar lab0_Hlt2CharmHadD02HHHHDecision_TOS("lab0_Hlt2CharmHadD02HHHHDecision_TOS","",0,1);
RooRealVar D0_Hlt1TrackAllL0Decision_TOS("D0_Hlt1TrackAllL0Decision_TOS","",0,1);
RooRealVar D0_ID("D0_ID","",-500,500);
RooRealVar D0_TAU("D0_TAU","D0_TAU",-500,500);


//***************Create dataset loading D0_M & Various Cut variables***************
RooDataSet * FullDataSet = new RooDataSet("FullDataSet", "FullDataSet", DecayTree, RooArgSet(D0_M,piplus1_P, piplus1_PT, D0_L0HadronDecision_TOS, lab0_L0Global_TIS, lab0_Hlt2CharmHadD02HHHHDecision_TOS, D0_Hlt1TrackAllL0Decision_TOS,D0_ID,lab0_M));

cout << "Number of events in data set:" << endl;
FullDataSet->Print() ;

//Make cuts on data and create 3 datasets; Full, D0 only and D0_Bar only
TCut piplus1_P_Cut = "piplus1_P > 3000"; 
TCut piplus1_PT_Cut = "piplus1_PT > 350";
TCut Trigger_Selection_1 = "(D0_L0HadronDecision_TOS == 1) || (lab0_L0Global_TIS == 1)";
TCut Trigger_Selection_2 = "lab0_Hlt2CharmHadD02HHHHDecision_TOS == 1";
TCut Trigger_Selection_3 = "D0_Hlt1TrackAllL0Decision_TOS == 1";
TCut Lab0_M_D0_M = "((lab0_M - D0_M) > 139.57) && ((lab0_M - D0_M) < 155)";
TCut D0_Data = "D0_ID>0";
TCut D0_bar_Data = "D0_ID<0";

RooDataSet* ReducedDataSet = (RooDataSet*)FullDataSet->reduce(piplus1_P_Cut);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(piplus1_PT_Cut);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(Trigger_Selection_1);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(Trigger_Selection_2);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(Trigger_Selection_3);
RooDataSet* ReducedDataSet = (RooDataSet*)ReducedDataSet->reduce(Lab0_M_D0_M);

RooDataSet* D0_DataSet = (RooDataSet*)ReducedDataSet->reduce(D0_Data);
RooDataSet* D0_Bar_DataSet = (RooDataSet*)ReducedDataSet->reduce(D0_bar_Data);



//****************************************************
//*************        PDF Shapes        *************
//****************************************************

//************D0 Mass PDF************

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
RooRealVar sigfrac_D0_M("sigfrac_D0_M","signal fraction",0.7,0,1);
RooAddPdf sig_D0_M_PDF("sig_D0_M_PDF","D0_M sig with yields",RooArgList(signal,CBall),RooArgList(sigfrac_D0_M));

//Background PDF
RooRealVar bkgD0_M("bkgD0_M","slope of background", 0, 0, 10, "MeV/c^{2}");
RooPolynomial bkg_D0_M_PDF("background_poly", "poly function for background", D0_M, RooArgList(bkgD0_M));

//Total D0_M PDF
RooAddPdf Total_D0_M_PDF("Total_D0_M_PDF","D0_M total",RooArgList(sig_D0_M_PDF,bkg_D0_M_PDF),RooArgList(sig_D0_M_Yield,bkg_D0_M_Yield));
Total_D0_M_PDF.fitTo(*ReducedDataSet, Extended() );


//************Delta mass PDF************

//Create Delta_M Variables
RooFormulaVar Total_Delta_M_f("Total_Delta_M_f","#Deltam","@0 - @1", RooArgList(lab0_M, D0_M));
RooFormulaVar D0_Delta_M_f("D0_Delta_M_f","#Deltam","@0 - @1", RooArgList(lab0_M, D0_M));
RooFormulaVar D0_bar_Delta_M_f("D0_bar_Delta_M_f","#Deltam","@0 - @1", RooArgList(lab0_M, D0_M));
RooRealVar* Total_Delta_M = (RooRealVar*)ReducedDataSet->addColumn(Total_Delta_M_f);
RooRealVar* D0_Delta_M = (RooRealVar*)D0_DataSet->addColumn(D0_Delta_M_f);
RooRealVar* D0_bar_Delta_M = (RooRealVar*)D0_Bar_DataSet->addColumn(D0_bar_Delta_M_f);
Total_Delta_M->setRange(Delta_M_Min, Delta_M_Max);  
D0_Delta_M->setRange(Delta_M_Min,Delta_M_Max);
D0_bar_Delta_M->setRange(Delta_M_Min, Delta_M_Max);

//Time integrated fit parameters
//For Johnson plots
  RooRealVar par0("par0", "p1", 2.07022e+01, 0, 100);
  RooRealVar par1("par1", "a", 3.60598e-01, 0, 0.5);
  RooRealVar par2("par2","gamma", -1.24512e-01, -5, 5);
  RooRealVar par3("par3","delta", 1.60301, 0, 5);
  RooRealVar par4("par4","mu",1.45245e+02, 139, 160);
  RooRealVar par5("par5","sigma",5.31443e-01, 0, 5);
  RooRealVar rrv_Delta_M_Min("rrv_Delta_M_Min","rrv_Delta_M_Min",139.57);

//Fitting Delta_M Signal and Background PDFs
  RooGenericPdf bkg_Delta_M_PDF("bkg_Delta_M_PDF","bkg_Delta_M_PDF","((Total_Delta_M_f-rrv_Delta_M_Min)+par0*((Total_Delta_M_f-rrv_Delta_M_Min)^2))^par1",RooArgList(*Total_Delta_M,rrv_Delta_M_Min,par0,par1));
  RooGenericPdf sig_Delta_M_PDF("sig_Delta_M_PDF","sig_Delta_M_PDF","(exp(-0.5*(par2+par3*TMath::ASinH((Total_Delta_M_f-par4)/par5))^2))/((1+((Total_Delta_M_f-par4)/par5)^2))^(.5)",RooArgList(*Total_Delta_M,par2,par3,par4,par5));

  RooAddPdf Delta_M_PDF("Delta_M_PDF","deltam sig+bkg",RooArgList(bkg_Delta_M_PDF,sig_Delta_M_PDF),RooArgList(bkg_D0_M_Yield,sig_D0_M_Yield));
  RooFitResult* Time_Integrated_Delta_M_Fit_Result = Delta_M_PDF.fitTo(*ReducedDataSet, Save() ) ;



//****************************************************
//*************       Canvas Plot        *************
//****************************************************

//************Plot D0 Mass Plot************

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
Total_D0_M_PDF.paramOn(D0_M_Frame, Format("NELU",AutoPrecision(2)),Layout(0.1,0.4,0.9)); //Display fit parameters
Total_D0_M_PDF.plotOn(D0_M_Frame, Components(bkg_D0_M_PDF), LineColor(kGreen),LineStyle(kDashed)); //Plot background
Total_D0_M_PDF.plotOn(D0_M_Frame, Components(sig_D0_M_PDF), LineColor(kBlue),LineStyle(kDashed)); //Plot signal
Total_D0_M_PDF.plotOn(D0_M_Frame, LineColor(kRed)); //Plot total PDF on frame D0_M_Frame

//Plot D0_M_Frame on pad1
pad1->cd(); //Makes pad1 the active pad
D0_M_Frame->Draw();
//Set up legend
leg->SetFillColor(kWhite);
leg->SetTextSize(0.04);
leg->SetTextColor(kBlack);
leg->SetHeader("D^{0} #rightarrow K^{+}#pi^{-}#pi^{+}#pi^{-}");
leg->AddEntry("h_data", "m(D^{0}) data", "lep");
leg->AddEntry("Total_D0_M_PDF_Norm[D0_M]","Total Fitting function PDF", "l");
leg->AddEntry("Total_D0_M_PDF_Norm[D0_M]_Comp[bkg_D0_M_PDF]","Background PDF","l");
leg->AddEntry("Total_D0_M_PDF_Norm[D0_M]_Comp[sig_D0_M_PDF]","Signal PDF","l");
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


//************Plot Delta Mass Plot************


  TCanvas canvas_time_int ;
  RooPlot* Delta_M_Only_Plot = Total_Delta_M->frame() ;
  ReducedDataSet.plotOn(Delta_M_Only_Plot) ;
  Delta_M_PDF.plotOn(Delta_M_Only_Plot) ;
  Delta_M_Only_Plot->Draw() ;
  canvas_time_int.SaveAs("LAUREN_DELTA_M.pdf") ;


}



