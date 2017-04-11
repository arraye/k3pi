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
#include "TH1F.h"

#include <vector>
#include <sstream>

using namespace RooFit;

void D0_To_K3pi_A_Gamma()
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
DecayTree->SetBranchStatus("D0_TAU",1); // D0 lifetime
DecayTree->SetBranchStatus("D0_ID",1);


//***************Declare observables (Var to plot & vars to cut on)***************
//D0_M
double D0_M_Min(1820),D0_M_Max(1910);  //Maximum and minimum values of D0_M to consider
RooRealVar D0_M("D0_M","",D0_M_Min,D0_M_Max);

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
RooRealVar D0_ID("D0_ID","D0_ID",-500,500);
RooRealVar D0_TAU("D0_TAU","",-500,500);


//***************Create dataset loading D0_M & Various Cut variables***************

  RooArgSet data_ArgSet(D0_M,lab0_M,D0_TAU,D0_L0HadronDecision_TOS,lab0_L0Global_TIS,lab0_Hlt2CharmHadD02HHHHDecision_TOS,D0_Hlt1TrackAllL0Decision_TOS,piplus1_P,piplus1_PT);
  data_ArgSet.add(D0_ID);

RooDataSet FullDataSet("FullDataSet","FullDataSet",data_ArgSet,Import (*DecayTree));

cout << "Number of events in data set:" << endl;
FullDataSet.Print() ;

//Make cuts on data and create 3 datasets; Full, D0 only and D0_Bar only
TCut piplus1_P_Cut = "piplus1_P > 3000"; 
TCut piplus1_PT_Cut = "piplus1_PT > 350";
TCut Trigger_Selection_1 = "(D0_L0HadronDecision_TOS == 1) || (lab0_L0Global_TIS == 1)";
TCut Trigger_Selection_2 = "lab0_Hlt2CharmHadD02HHHHDecision_TOS == 1";
TCut Trigger_Selection_3 = "D0_Hlt1TrackAllL0Decision_TOS == 1";
TCut Lab0_M_D0_M = "((lab0_M - D0_M) > 139.57) && ((lab0_M - D0_M) < 155)";
TCut D0_Data = "D0_ID>0";
TCut D0_bar_Data = "D0_ID<0";

RooDataSet* ReducedDataSet = (RooDataSet*)FullDataSet.reduce(piplus1_P_Cut);
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
RooRealVar CB_n("CB_n","CB1 n",4,0.01,10);
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
  Delta_M_PDF.fitTo(*ReducedDataSet, Extended());

par0.setConstant(kTRUE);
par1.setConstant(kTRUE);
par2.setConstant(kTRUE);
par3.setConstant(kTRUE);
par4.setConstant(kTRUE);
par5.setConstant(kTRUE);



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
D0_M_Canvas->SaveAs("Plots/LAUREN_D0_M.pdf");



//************Plot Delta Mass Plot************


  TCanvas canvas_time_int ;
  RooPlot* Delta_M_Only_Plot = Total_Delta_M->frame() ;
  ReducedDataSet->plotOn(Delta_M_Only_Plot) ;
  Delta_M_PDF.plotOn(Delta_M_Only_Plot) ;
  Delta_M_Only_Plot->Draw() ;
  canvas_time_int.SaveAs("Plots/LAUREN_DELTA_M.pdf") ;




//****************************************************
//*************       Bin Edges          *************
//****************************************************

  double bins=20; // # of bins
  //double bins=10; //TEMPORARY!!!!!!

  double frequency=100; // Fineness of binning
  double totalCand=ReducedDataSet->sumEntries();
  double leftCand=totalCand; //Candidates left initialised to total number
  double entryPerBin=totalCand/bins; 
  double dTau_max=6.; // Maximum lifetime considered: 6 ps
  double dTau=dTau_max/bins/frequency; //Possibly redundant (for for loop)?
  double nEntries;

  std::vector<double>bin_edges;
  bin_edges.push_back(0.);
  double current_position=0.;

  std::stringstream le_s; // left edge
  le_s<<current_position;
  std::stringstream re_s; // right edge

  int iterBin=0;

//For percentage display
float count_max= bins*frequency+1;
int count_display=0;
float progress_percentage=0;

//Create array holding nbins+1 containing low-edge of nbins first bins and the upper-edge of last bin
//(root.cern.ch/root/HowtoHistogram.html

cout<<"Creating bin_edges vector..."<<endl;
  for (int iterTau=0;iterTau<(bins*frequency+1);iterTau++)  // so for loop will iterate up to a maximum value of dTau_max
    {
//For percentage display
    progress_percentage=(count_display/count_max)*100;
    cout<<"\r"<<"Current progress: "<<fixed<<setprecision(2)<<progress_percentage<<"%";
    count_display++;
 

    current_position=iterTau*dTau;
    re_s.str(std::string());re_s.clear();  //Clear previous value of right edge
    re_s<<current_position;

    nEntries=ReducedDataSet->sumEntries(TString(le_s.str())+"<D0_TAU*1000.0&&D0_TAU*1000.0<"+TString(re_s.str()));  // Calculate how many entries will be in bin (effectively cutting on range D0_TAUs)
   
 //cout<<TString(le_s.str())+"<D0_TAU*1000.&&D0_TAU*1000.<"+TString(re_s.str())<<" # of entries= "<<nEntries<<endl;

    //std::cout<<"nEntries= "<<nEntries<<" and entryPerBin= "<<entryPerBin<<endl;

    if (nEntries >= entryPerBin) 
	    {
	      bin_edges.push_back(current_position);
	      iterBin++;
	      leftCand-=nEntries;
	      //std::cout<<"Bin #"<<iterBin<<"; leftCand= "<<leftCand<<endl;
	      le_s.str(std::string());le_s.clear();
	      le_s<<current_position;
	    }
   
    }
    bin_edges.push_back(dTau_max);  //add dTau_max to end of bin_edges vector
    std::cout<<"size of the bin_edges= "<<bin_edges.size()<<endl; // bin_edges.size(): number of elements in vector (expect 101)



//****************************************************
//*************          Yields          *************
//****************************************************

//Vectors for yield parameters
std::vector<double>D0_SigYield_vector;
std::vector<double>D0_SigYield_error;
std::vector<double>D0_Bar_SigYield_vector;
std::vector<double>D0_Bar_SigYield_error;

std::vector<double>D0_BkgYield_vector;
std::vector<double>D0_BkgYield_error;
std::vector<double>D0_Bar_BkgYield_vector;
std::vector<double>D0_Bar_BkgYield_error;

//Vector for looping over D0/D0bar data
std::vector<RooDataSet*>Flavour;

std::cout<<"D0_dataset_size = "<<D0_DataSet->sumEntries()<<endl;
std::cout<<"D0__Bar_dataset_size = "<<D0_Bar_DataSet->sumEntries()<<endl;



//****************************************************
//************* Bin Centres & data in bins fits *************
//****************************************************


//For D0-TAU plot
double no_of_d0_tau_bins = bin_edges.size()
double entries_per_bin[no_of_d0_tau_bins]
double d0_tau_bin_centres_for_TGraph[no_of_d0_tau_bins];


//Iterate over D0 and D0_bar
int i_flavour_iterate=0;
for(i_flavour_iterate=0; i_flavour_iterate < 2; i_flavour_iterate++)
  {

    RooDataSet* current_dataset;
    current_dataset = i_flavour_iterate==0 ? D0_DataSet : D0_Bar_DataSet;
    string str_flavour = i_flavour_iterate == 0 ? string("D0_") : string("D0_Bar_");

    std::vector<RooDataSet*>dataInBin;
    // Generate a vector saving locations of bin centres
    std::vector<double>bin_centres;


   // Find the bin centres
    std::vector<double>::iterator bin_iterate;    
    int i_bin_iterate=0;
    // Find the bin centres (why?) 
    for(bin_iterate=bin_edges.begin(); i_bin_iterate <= bin_edges.size(); bin_iterate++)
      {
       std::stringstream betL_s;    //What are betL_s and betR_s? temp bin edges?
       betL_s<<bin_edges[i_bin_iterate];
       std::stringstream betR_s;
       betR_s<<bin_edges[i_bin_iterate+1];


       //bin_centre = centre of left and right bin_edge
       bin_centres.push_back((bin_edges[i_bin_iterate+1]-bin_edges[i_bin_iterate]/2)+bin_edges[i_bin_iterate]);

       dataInBin.push_back((RooDataSet*) (current_dataset.reduce(RooFit::Cut(TString(betL_s.str())+"< D0_TAU*1000.&&D0_TAU*1000.<="+TString(betR_s.str())))));

       RooDataSet* tmpDataSet=(RooDataSet*) (dataInBin[i_bin_iterate]); // Create a dataset from the data that conforms to above cut (D0_TAU cut) I think for test purposes??!!!

       //Test
       cout<<"test: "<<TString(betL_s.str())+"< D0_TAU*1000.&&D0_TAU*1000.<="+TString(betR_s.str())<<"; bin size = "<<tmpDataSet->sumEntries()<<endl;

       entries_per_bin[i_bin_iterate].push_back(*tmpDataSet->sumEntries());
       d0_tau_bin_centres_for_TGraph[i_bin_iterate] = bin_centres[i_bin_iterate];
      
       i_bin_iterate++;

     } //end find bin centres


//Plot a D0_TAU plot.
const unsigned int number_of_bin_centres = bin_centres.size();
TGraph *D0_TAU_Plot = new TGraph(number_of_bin_centres,d0_tau_bin_centres_for_TGraph,entries_per_bin);





      // Fit data for each bin and save plots
      std::vector<double>::iterator fit_iterate;
      int i_fit_iterate=0;
      for(fit_iterate=bin_edges.begin(); i_fit_iterate<bin_edges.size(); fit_iterate++, i_fit_iterate++)
          {
                std::cout<<"test 1"<<endl;

		RooDataSet* tmpDataset=(RooDataSet*) (dataInBin[i_fit_iterate]);
                RooRealVar* Total_Delta_M = (RooRealVar*)tmpDataset->addColumn(Total_Delta_M_f);
                Total_Delta_M->setRange(Delta_M_Min, Delta_M_Max); 
                
                //Set values... why? Why not use earlier values? Resetting for new fit I guess
                sig_D0_M_Yield.setVal(0.7 * tmpDataset->numEntries());  
                bkg_D0_M_Yield.setVal(0.1 * tmpDataset->numEntries());

                RooFitResult* D0_M_Fit_Result = Total_D0_M_PDF.fitTo(*tmpDataset, Save());
                RooFitResult* Delta_M_Fit_Result = Delta_M_PDF.fitTo(*tmpDataset, Save()); 

                std::cout<<"test 2"<<endl;
                std::stringstream isstr;
                isstr << i_fit_iterate;
                string isstr_s = isstr.str();


                //Plot D0(bar)_M for each decay time bin ++++++++++++++++++and save
                TCanvas canvas_after_bins;
                RooPlot* mass_plot = D0_M.frame();
                tmpDataset->plotOn(mass_plot);
                Total_D0_M_PDF->plotOn(mass_plot);
                mass_plot->Draw();
               canvas_after_bins.SaveAs(("Plots/" + str_flavour + "Fit_" + isstr_s + ".pdf").c_str()) ; //c_str = stringcopy


 //Not working?
                //Plot Delta_M for each decay time bin, and save
                TCanvas canvas_after_bins;
                RooPlot* delta_m_plot = Total_Delta_M->frame();
                tmpDataset->plotOn(delta_m_plot);
                Delta_M_PDF->plotOn(delta_m_plot);
                delta_m_plot->Draw();
                canvas_after_bins.SaveAs(("Plots/" + str_flavour + "Delta_M_Fit_" + isstr_s + ".pdf").c_str()) ; //c_str = stringcopy



                //Save yields of fits in yield vectors, for use in A_Gamma determination
                RooRealVar* sig_D0_M_Yield_PTR=&sig_D0_M_Yield;
                RooRealVar* bkg_D0_M_Yield_PTR=&bkg_D0_M_Yield;

                if (i_flavour_iterate==0)
			{
			D0_SigYield_vector.push_back(sig_D0_M_Yield_PTR->getVal());
			D0_SigYield_error.push_back(sig_D0_M_Yield_PTR->getError());                
			D0_BkgYield_vector.push_back(bkg_D0_M_Yield_PTR->getVal());
			D0_BkgYield_error.push_back(bkg_D0_M_Yield_PTR->getError());      
			}
		else
		        {  
		        D0_Bar_SigYield_vector.push_back(sig_D0_M_Yield_PTR->getVal());
		        D0_Bar_SigYield_error.push_back(sig_D0_M_Yield_PTR->getError());                
		        D0_Bar_BkgYield_vector.push_back(bkg_D0_M_Yield_PTR->getVal());
		        D0_Bar_BkgYield_error.push_back(bkg_D0_M_Yield_PTR->getError());    
		        }

                } // end interate fit for each bin 

std::cout<<"end of bin fit iterate"<<endl;

} //end iterate over D0 and D0_bar
std::cout<<"end of d0_d0bar iterate"<<endl;



//****************************************************
//*********** Filling & plotting A_Gamma *************
//****************************************************


    //std::vector<double>A_gamma;
    const unsigned int number_of_bins = bin_centres.size();
    double A_Gamma_Value;
    double A_Gamma_For_TGraph[number_of_bins];
    double bin_centres_for_TGraph[number_of_bins];

    double denominator;
    double numerator;


    for (int A_Gamma_iterate=0; A_Gamma_iterate <= number_of_bins; A_Gamma_iterate++)
      {
      numerator = D0_SigYield_vector[A_Gamma_iterate]- D0_Bar_SigYield_vector[A_Gamma_iterate];
      denominator = D0_SigYield_vector[A_Gamma_iterate]+ D0_Bar_SigYield_vector[A_Gamma_iterate];

      if (numerator == 0 || denominator == 0)
      	{ A_Gamma_Value = 0; }
      else
     	{ A_Gamma_Value = numerator/denominator; }

        A_Gamma_For_TGraph[A_Gamma_iterate] = A_Gamma_Value;
        bin_centres_for_TGraph[A_Gamma_iterate] = bin_centres[A_Gamma_iterate];

	std::cout<<"A_Gamma_For_TGraph["<<A_Gamma_iterate<<"] ="<<A_Gamma_For_TGraph[A_Gamma_iterate]<<endl;
	std::cout<<"bin_centres_for_TGraph["<<A_Gamma_iterate<<"] ="<<bin_centres_for_TGraph[A_Gamma_iterate]<<endl;

      }      

	A_Gamma_Canvas = new TCanvas("A_Gamma_Canvas","Simple A_Gamma plot",200,10,700,500);
	A_Gamma_Canvas->cd();
	A_Gamma_Canvas->Update();
	A_Gamma_Canvas->SetGrid();
	A_Gamma_Canvas->GetFrame()->SetFillColor(21);
	A_Gamma_Canvas->GetFrame()->SetBorderSize(12);
	    
	TGraph *A_Gamma_Graph = new TGraph(number_of_bins,bin_centres_for_TGraph,A_Gamma_For_TGraph);

	A_Gamma_Canvas->SetTitle("A_{#Gamma} vs #tau_{D^{0}}");
	A_Gamma_Graph->SetMarkerColor(4);
	A_Gamma_Graph->SetMarkerStyle(21);
	A_Gamma_Graph->GetXaxis()->SetLimits(0.,7.);
	A_Gamma_Graph->GetXaxis()->SetTitle("#tau_{D^{0}} #times 1000 [ps]");
	A_Gamma_Graph->GetXaxis()->CenterTitle();
	A_Gamma_Graph->GetYaxis()->SetTitle("A_{#Gamma}");
	A_Gamma_Graph->GetYaxis()->CenterTitle();
	A_Gamma_Graph->Draw("ALP");
	A_Gamma_Canvas->Update();

	A_Gamma_Canvas->SaveAs("Plots/A_Gamma.pdf");




} 
//end D_M__Delta_M_Fits()














          	  


