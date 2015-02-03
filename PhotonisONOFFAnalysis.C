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
#include "TEfficiency.h"

void MakeAmplitudePlot(string filename, string plotname, double scalefactor, double fitmin, double fitmax) {
  // Get the tree


  TFile *inputfile = new TFile(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("tree");

  // get the variables from the ntuple
  float t1gausroot = 0;
  float t2gausroot = 0;
  float t3gausroot = 0;
  float t4gausroot = 0;
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;
  float ch1Int = 0;
  float ch2Int = 0;
  float ch3Int = 0;
  float ch4Int = 0;
  unsigned int ch1QualityBit = 0;
  unsigned int ch2QualityBit = 0;
  unsigned int ch3QualityBit = 0;
  unsigned int ch4QualityBit = 0;

  tree->SetBranchAddress("t1gausroot",&t1gausroot);
  tree->SetBranchAddress("t2gausroot",&t2gausroot);
  tree->SetBranchAddress("t3gausroot",&t3gausroot);
  tree->SetBranchAddress("t4gausroot",&t4gausroot);
  tree->SetBranchAddress("ch1Amp",&ch1Amp);
  tree->SetBranchAddress("ch2Amp",&ch2Amp);
  tree->SetBranchAddress("ch3Amp",&ch3Amp);
  tree->SetBranchAddress("ch4Amp",&ch4Amp);
  tree->SetBranchAddress("ch1Int",&ch1Int);
  tree->SetBranchAddress("ch2Int",&ch2Int);
  tree->SetBranchAddress("ch3Int",&ch3Int);
  tree->SetBranchAddress("ch4Int",&ch4Int);
  tree->SetBranchAddress("ch1QualityBit",&ch1QualityBit);
  tree->SetBranchAddress("ch2QualityBit",&ch2QualityBit);
  tree->SetBranchAddress("ch3QualityBit",&ch3QualityBit);
  tree->SetBranchAddress("ch4QualityBit",&ch4QualityBit);

  //create histograms
  TH1F *dt;
  TH1F *histAmplitude;
  histAmplitude = new TH1F("histAmplitude","; Amplitude [V];Number of Events",50,0,0.5);
  dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 1000, -10,10);
  
  double Denominator = 0;
  double Numerator = 0;

  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();
  std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
      tree->GetEntry(iEntry);    
      
      //require cherenkov and signal in the front MCP
      if (ch1Amp > 0.05 && ch4Amp > 0.1) {
	histAmplitude->Fill(ch3Amp);
	if (ch3Amp > 0.02 
	    && ch1Amp < 0.49 && ch3Amp < 0.49
	    //&& ch1QualityBit == 0 && ch2QualityBit == 0
	    ) {
	  dt->Fill(t3gausroot - t1gausroot);
	}

	if (ch2Amp > 0.02) {
	  Denominator++;
	  if (ch3Amp > 0.02) Numerator++;
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
  histAmplitude->GetYaxis()->SetTitleOffset(1.5);
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
  tex->DrawLatex(0.55, 0.80, Form("#sigma/#mu = %.1f %s",100*fitter->GetParameter(2)/fitter->GetParameter(1),"%"));
  tex->DrawLatex(0.55, 0.85, Form("Mean = %.2f #pm %.2f %s",fitter->GetParameter(1),fitter->GetParError(1),"V"));
  tex->DrawLatex(0.15, 0.92, Form("Attenuation Factor = %.3f",scalefactor));
  
  c->SaveAs( Form("%s_amplitude.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_amplitude.pdf", plotname.c_str()) );


  //time resolution plot
  c = new TCanvas("c","c",600,600);
  
  dt->SetAxisRange(2.5,4,"X");
  dt->SetTitle("");
  dt->GetXaxis()->SetTitle("#Delta t [ns]");
  dt->GetYaxis()->SetTitle("Number of Events");
  dt->GetYaxis()->SetTitleOffset(1.5);
  dt->SetMaximum(1.2*dt->GetMaximum());
  dt->Draw();
  dt->SetStats(0);
  dt->Fit("gaus","","",2.5,4);
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.040);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.55, 0.80, Form("#sigma = %.1f #pm %.1f ps",1000*fitter->GetParameter(2),1000*fitter->GetParError(2)));
  tex->DrawLatex(0.55, 0.85, Form("Mean = %.2f #pm %.2f %s",fitter->GetParameter(1),fitter->GetParError(1),"V"));
  
  c->SaveAs( Form("%s_dt.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_dt.pdf", plotname.c_str()) );

  double effErrlow = Numerator / Denominator  - TEfficiency::ClopperPearson((UInt_t)Denominator, (UInt_t)Numerator, 0.68269, kFALSE);
  double effErrhigh = TEfficiency::ClopperPearson((UInt_t)Denominator, (UInt_t)Numerator, 0.68269, kTRUE) - Numerator / Denominator;
  cout << "Photonis Detection Efficiency (Amplitude > 20mV) : " << Numerator / Denominator 
       << " + " << effErrhigh << " - " << effErrlow
       << "\n";

  
}


Double_t TOFResolutionFunction(Double_t *x,Double_t *par)
{
  Double_t arg = 0;
//if (par[2] != 0) arg = (x[0] - par[1])/par[2];
  arg = x[0];
  Double_t fitval = par[0]*(1.0 / sqrt(arg)) + par[1];
  return fitval;
}


void MakeAmplitudeVsShowerDepthGraph() {

  //use beam energy for xaxis
  float x_photocathodeON[5] = { 2.5, 4.5, 6.5, 8.5, 10.5 };
  float xerr_photocathodeON[5] = { 0.25, 0.25, 0.25, 0.25, 0.25 };
  float y_photocathodeON[5] = {  0.12, 0.19, 0.19, 0.21, 0.14 };
  float yerr_photocathodeON[5] = { 0.01, 0.01, 0.01, 0.03, 0.01 };

  float eff_photocathodeON[4] = {  0.985, 0.973, 0.990, 0.995, 0.944 };
  float efferrl_photocathodeON[4] = { 0.009, 0.011, 0.008, 0.007, 0.032 };
  float efferrh_photocathodeON[4] = { 0.006, 0.008, 0.005, 0.003, 0.022 };

  float x_photocathodeOFF[5] = { 2.5, 4.5, 6.5, 8.5, 10.5 };
  float xerr_photocathodeOFF[5] = { 0.25, 0.25, 0.25, 0.25, 0.25 };
  float y_photocathodeOFF[5] = {  0.13, 0.32, 0.31, 0.24, 0.24 };
  float yerr_photocathodeOFF[5] = { 0.01, 0.03, 0.01, 0.03, 0.01 };

  float eff_photocathodeOFF[4] = {  0.991, 1.0, 1.0, 1.0, 0.997 };
  float efferrl_photocathodeOFF[4] = { 0.007, 0.004, 0.005, 0.008, 0.008 };
  float efferrh_photocathodeOFF[4] = { 0.004, 0.0, 0.0, 0.0, 0.003 };


  TGraphErrors *graphON = new TGraphErrors(6,x_photocathodeON,y_photocathodeON,xerr_photocathodeON,yerr_photocathodeON);
  graphLead->SetLineWidth(2);
  TGraphErrors *graphOFF = new TGraphErrors(6,x_photocathodeOFF,y_photocathodeOFF,xerr_photocathodeOFF,yerr_photocathodeOFF);
  graphTungsten->SetLineWidth(2);
  graphTungsten->SetLineColor(kBlue);

  TCanvas * c = new TCanvas("c","c",800,600);
  graphLead->Draw("AP");
  graphTungsten->Draw("Psame");

  TLegend *legend = new TLegend (0.65,0.7,0.8,0.85);
  legend->SetTextSize(0.03);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(graphLead, "Lead Absorber","LP");
  legend->AddEntry(graphTungsten, "Tungsten Absorber","LP");
  legend->Draw();

  graphLead->SetTitle("");
  graphLead->GetXaxis()->SetTitle("Absorber Thickness [X_{0}]");
  graphLead->GetXaxis()->SetTitleOffset(1.2);
  graphLead->GetYaxis()->SetTitle("Mean Amplitude [V]");
  graphLead->GetYaxis()->SetTitleOffset(1.25);
  graphLead->GetYaxis()->SetRangeUser(0,0.6);

  c->SaveAs( "AmplitudeVsLeadThickness.gif" );
  c->SaveAs( "AmplitudeVsLeadThickness.pdf" );

}


void MakeTimeResolutionVsShowerDepthGraph() {

 
  float x_photocathodeON[5] = { 2.5, 4.5, 6.5, 8.5, 10.5 };
  float xerr_photocathodeON[5] = { 0.25 , 0.25, 0.25, 0.25, 0.25 };
  float y_photocathodeON[5] = { 59 , 62 , 56 , 54, 60 };
  float yerr_photocathodeON[5] = { 3 , 3, 2 , 2 , 3 };

   float x_photocathodeOFF[5] = { 2.5, 4.5, 6.5, 8.5, 10.5 };
  float xerr_photocathodeOFF[5] = { 0.25 , 0.25, 0.25, 0.25, 0.25 };
  float y_photocathodeOFF[5] = { 60 , 49 , 48 , 47 , 74 };
  float yerr_photocathodeOFF[5] = { 2 , 2, 2 , 2 , 4 };



  TGraphErrors *graphON = new TGraphErrors(6,x_photocathodeON,y_photocathodeON,xerr_photocathodeON,yerr_photocathodeON);
  graphON->SetLineWidth(2);
  TGraphErrors *graphOFF = new TGraphErrors(6,x_photocathodeOFF,y_photocathodeOFF,xerr_photocathodeOFF,yerr_photocathodeOFF);
  graphOFF->SetLineWidth(2);
  graphOFF->SetLineColor(kBlue);

  TCanvas * c = new TCanvas("c","c",800,600);
  graphON->Draw("AP");
  graphOFF->Draw("Psame");

  TLegend *legend = new TLegend (0.50,0.7,0.8,0.85);
  legend->SetTextSize(0.04);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(graphOFF, "Photocathode OFF","LP");
  legend->AddEntry(graphON, "Photocathode ON","LP");
  legend->Draw();

  graphON->SetTitle("");
  graphON->GetXaxis()->SetTitle("Tungsten Absorber Thickness [X_{0}]");
  graphON->GetXaxis()->SetTitleOffset(1.2);
  graphON->GetYaxis()->SetTitle("Time Resolution [ps]");
  graphON->GetYaxis()->SetTitleOffset(1.25);
  graphON->GetYaxis()->SetRangeUser(0,120);

  c->SaveAs( "PhotonisTimeResolutionVsAbsorberThickness.gif" );
  c->SaveAs( "PhotonisTimeResolutionVsAbsorberThickness.pdf" );

}






void PhotonisONOFFAnalysis() {

  //*************************************
  // Photocathode ON Vs Shower Depth
  //*************************************
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_61.root.ana_lowpassfilter.root","PhotonisPhotocathodeON_Electron_NoAbsorber_8GeV", 1.0, 0.03, 0.3);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_60.root.ana_lowpassfilter.root","PhotonisPhotocathodeON_Electron_2X0Tungsten_8GeV", 1.0, 0.05, 0.4);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_62.root.ana_lowpassfilter.root","PhotonisPhotocathodeON_Electron_4X0Tungsten_8GeV", 1.0, 0.05, 0.4);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_63.root.ana_lowpassfilter.root","PhotonisPhotocathodeON_Electron_6X0Tungsten_8GeV", 1.0, 0.05, 0.4);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_64.root.ana_lowpassfilter.root","PhotonisPhotocathodeON_Electron_8X0Tungsten_8GeV", 1.0, 0.05, 0.4);

  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_52.root.ana_lowpassfilter.root","PhotonisPhotocathodeOFF_Electron_NoAbsorber_8GeV", 1.0, 0.03, 0.3);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_51.root.ana_lowpassfilter.root","PhotonisPhotocathodeOFF_Electron_2X0Tungsten_8GeV", 1.0, 0.01, 0.5);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_53.root.ana_lowpassfilter.root","PhotonisPhotocathodeOFF_Electron_4X0Tungsten_8GeV", 1.0, 0.01, 0.5);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_54.root.ana_lowpassfilter.root","PhotonisPhotocathodeOFF_Electron_6X0Tungsten_8GeV", 1.0, 0.01, 0.5);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_56And57.root.ana_lowpassfilter.root","PhotonisPhotocathodeOFF_Electron_8X0Tungsten_8GeV", 1.0, 0.01, 0.5);


  // MakeAmplitudeVsShowerDepthGraph();
  MakeTimeResolutionVsShowerDepthGraph();


  //
  //*************************************
  //Results for Y11 Fibers from May Testbeam
  //*************************************
//   MakeTimeResolutionPlot("cpt_may_run_131.ana.root","TOF_ShashlikY11Fiber_Electron_4GeV",4,false);
//   MakeTimeResolutionPlot("cpt_may_run_132.ana.root","TOF_ShashlikY11Fiber_Electron_8GeV",8,false);
//   MakeTimeResolutionPlot("cpt_may_run_133.ana.root","TOF_ShashlikY11Fiber_Electron_16GeV",16,false);
//   MakeTimeResolutionVsEnergyPlot_Y11();




}
