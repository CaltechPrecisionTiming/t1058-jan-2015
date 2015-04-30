#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#endif

#include <iostream>
#include <fstream> 
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TVirtualFitter.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include <math.h> 



void MakePlotAlternativeFormat(string filename, string plotname, double scalefactor, double fitmin, double fitmax) {
  // Get the tree


  TFile *inputfile = new TFile(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("tree");

  // get the variables from the ntuple
  float t1 = 0;
  float t2 = 0;
  float t3 = 0;
  float t4 = 0;
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;

  tree->SetBranchAddress("ch1Time",&t1);
  tree->SetBranchAddress("ch2Time",&t2);
  tree->SetBranchAddress("ch3Time",&t3);
  tree->SetBranchAddress("ch4Time",&t4);
  tree->SetBranchAddress("ch1Amp",&ch1Amp);
  tree->SetBranchAddress("ch2Amp",&ch2Amp);
  tree->SetBranchAddress("ch3Amp",&ch3Amp);
  tree->SetBranchAddress("ch4Amp",&ch4Amp);

  //create histograms
  TH1F *dt;
  TH1F *histAmplitude;
  histAmplitude = new TH1F("histAmplitude","; Amplitude [V];Number of Events",50,0,0.5);
  dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 500, -1,1);
  
  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();
  std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
      tree->GetEntry(iEntry);    
      
      //require cherenkov and signal in the front MCP
      if (ch1Amp > 50 ) {
	histAmplitude->Fill(ch2Amp/1000);
	if (ch2Amp > 50 
	    && ch1Amp < 490 && ch2Amp < 490
	    ) {
	  dt->Fill(t2 - t1);
	}
      }
  }

  TCanvas * c = 0;


  //Energy plot
  c = new TCanvas("c","c",600,600);
  
  histAmplitude->SetAxisRange(0.0,0.5,"X");
  histAmplitude->SetTitle("");
  histAmplitude->GetXaxis()->SetTitle("Pulse Height [V]");
  histAmplitude->GetYaxis()->SetTitle("Number of Events");
  histAmplitude->GetYaxis()->SetTitleOffset(1.3);
  histAmplitude->SetMaximum(1.2*histAmplitude->GetMaximum());
  histAmplitude->Draw();
  histAmplitude->SetStats(0);
  histAmplitude->Fit("gaus","","",fitmin,fitmax);
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.040);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  //tex->DrawLatex(0.55, 0.80, Form("#sigma/#mu = %.0f %s",100*fitter->GetParameter(2)/fitter->GetParameter(1),"%"));
  tex->DrawLatex(0.55, 0.85, Form("Mean = %.2f #pm %.2f %s",fitter->GetParameter(1),TMath::Max(0.01,fitter->GetParError(1)),"V"));
  tex->DrawLatex(0.15, 0.92, Form("Attenuation Factor = %.3f",scalefactor));
  
  c->SaveAs( Form("%s_amplitude.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_amplitude.pdf", plotname.c_str()) );


  //time resolution plot
  c = new TCanvas("c","c",600,600);
  
  dt->SetAxisRange(0.1,0.4,"X");
  dt->SetTitle("");
  dt->GetXaxis()->SetTitle("#Delta t [ns]");
  dt->GetXaxis()->SetTitleSize(0.045);
  dt->GetXaxis()->SetLabelSize(0.040);
  dt->GetXaxis()->SetLabelOffset(0.012);
  dt->GetXaxis()->SetTitleOffset(1.0);
  dt->GetYaxis()->SetTitle("Number of Events");
  dt->GetYaxis()->SetTitleOffset(1.04);
  dt->GetYaxis()->SetTitleSize(0.045);
  dt->GetYaxis()->SetLabelSize(0.040);
  dt->SetMaximum(1.2*dt->GetMaximum());
  dt->Draw();
  dt->SetStats(0);
  dt->Fit("gaus","","",0.2,0.3);
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.55, 0.80, Form("#sigma = %.1f #pm %.1f ps",1000*fitter->GetParameter(2),1000*fitter->GetParError(2)));
  //tex->DrawLatex(0.55, 0.85, Form("Mean = %.2f #pm %.2f %s",fitter->GetParameter(1),fitter->GetParError(1),"V"));
  tex->DrawLatex(0.55, 0.85, Form("Mean = %.2f ns",fitter->GetParameter(1)));
  
  c->SaveAs( Form("%s_dt.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_dt.pdf", plotname.c_str()) );


  
}



void ProtonTimingAnalysis() {


  //*************************************
  // Make plot from heejong analysis ntuple
  //*************************************
  MakePlotAlternativeFormat("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_3.ana_heejongCode.root","Proton_DeltaT_120GeV", 1.0, 0.1, 0.3);




}
