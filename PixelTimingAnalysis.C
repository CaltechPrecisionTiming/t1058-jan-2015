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


void MakePlot(string filename, string plotname, int refChannel, int detChannel, double scalefactor, double fitmin, double fitmax, double timefitmin, double timefitmax) {
  // Get the tree


  TFile *inputfile = new TFile(filename.c_str(),"READ");
  TTree *tree = (TTree*)inputfile->Get("tree");

  // get the variables from the ntuple
  float t1gausroot = 0;
  float t2gausroot = 0;
  float t3gausroot = 0;
  float t4gausroot = 0;
  float t5gausroot = 0;
  float t6gausroot = 0;
  float t7gausroot = 0;
  float t8gausroot = 0;
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;
  float ch5Amp = 0;
  float ch6Amp = 0;
  float ch7Amp = 0;
  float ch8Amp = 0;
  float ch1Int = 0;
  float ch2Int = 0;
  float ch3Int = 0;
  float ch4Int = 0;
  float ch5Int = 0;
  float ch6Int = 0;
  float ch7Int = 0;
  float ch8Int = 0;
  unsigned int ch1QualityBit = 0;
  unsigned int ch2QualityBit = 0;
  unsigned int ch3QualityBit = 0;
  unsigned int ch4QualityBit = 0;
  unsigned int ch5QualityBit = 0;
  unsigned int ch6QualityBit = 0;
  unsigned int ch7QualityBit = 0;
  unsigned int ch8QualityBit = 0;

  tree->SetBranchAddress("t1gausroot",&t1gausroot);
  tree->SetBranchAddress("t2gausroot",&t2gausroot);
  tree->SetBranchAddress("t3gausroot",&t3gausroot);
  tree->SetBranchAddress("t4gausroot",&t4gausroot);
  tree->SetBranchAddress("t5gausroot",&t5gausroot);
  tree->SetBranchAddress("t6gausroot",&t6gausroot);
  tree->SetBranchAddress("t7gausroot",&t7gausroot);
  tree->SetBranchAddress("t8gausroot",&t8gausroot);
  tree->SetBranchAddress("ch1Amp",&ch1Amp);
  tree->SetBranchAddress("ch2Amp",&ch2Amp);
  tree->SetBranchAddress("ch3Amp",&ch3Amp);
  tree->SetBranchAddress("ch4Amp",&ch4Amp);
  tree->SetBranchAddress("ch5Amp",&ch5Amp);
  tree->SetBranchAddress("ch6Amp",&ch6Amp);
  tree->SetBranchAddress("ch7Amp",&ch7Amp);
  tree->SetBranchAddress("ch8Amp",&ch8Amp);
  tree->SetBranchAddress("ch1Int",&ch1Int);
  tree->SetBranchAddress("ch2Int",&ch2Int);
  tree->SetBranchAddress("ch3Int",&ch3Int);
  tree->SetBranchAddress("ch4Int",&ch4Int);
  tree->SetBranchAddress("ch5Int",&ch5Int);
  tree->SetBranchAddress("ch6Int",&ch6Int);
  tree->SetBranchAddress("ch7Int",&ch7Int);
  tree->SetBranchAddress("ch8Int",&ch8Int);
  tree->SetBranchAddress("ch1QualityBit",&ch1QualityBit);
  tree->SetBranchAddress("ch2QualityBit",&ch2QualityBit);
  tree->SetBranchAddress("ch3QualityBit",&ch3QualityBit);
  tree->SetBranchAddress("ch4QualityBit",&ch4QualityBit);
  tree->SetBranchAddress("ch5QualityBit",&ch5QualityBit);
  tree->SetBranchAddress("ch6QualityBit",&ch6QualityBit);
  tree->SetBranchAddress("ch7QualityBit",&ch7QualityBit);
  tree->SetBranchAddress("ch8QualityBit",&ch8QualityBit);

  //create histograms
  TH1F *dt;
  TH1F *histAmplitude;
  histAmplitude = new TH1F("histAmplitude","; Amplitude [V];Number of Events",50,0,0.5);
  dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 2000, -10,10);
  
  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();
  std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
      tree->GetEntry(iEntry);    
      
      //require cherenkov and signal in the front MCP
      double refAmp = 0;
      double detAmp = 0;
      if (refChannel == 1) refAmp = ch1Amp;
	  else if(refChannel == 2) refAmp = ch2Amp;
	  else if(refChannel == 3) refAmp = ch3Amp;
	  else if(refChannel == 4) refAmp = ch4Amp;
	  else if(refChannel == 5) refAmp = ch5Amp;
	  else if(refChannel == 6) refAmp = ch6Amp;
	  else if(refChannel == 7) refAmp = ch7Amp;
	  else if(refChannel == 8) refAmp = ch8Amp;
      if (detChannel == 1) detAmp = ch1Amp;
	  else if(detChannel == 2) detAmp = ch2Amp;
	  else if(detChannel == 3) detAmp = ch3Amp;
	  else if(detChannel == 4) detAmp = ch4Amp;
	  else if(detChannel == 5) detAmp = ch5Amp;
	  else if(detChannel == 6) detAmp = ch6Amp;
	  else if(detChannel == 7) detAmp = ch7Amp;
	  else if(detChannel == 8) detAmp = ch8Amp;

      double deltaTAcrossDRSBoard = 0;
      if (refChannel >= 5 && detChannel < 5) deltaTAcrossDRSBoard = t5gausroot - t1gausroot;
      if (refChannel < 5 && detChannel >= 5) deltaTAcrossDRSBoard = t1gausroot - t5gausroot;

      if (ch1Amp > 0.05 && ch5Amp > 0.05 && ch4Amp > 0.1 ) {
	histAmplitude->Fill(detAmp);
	if (refAmp > 0.02 && detAmp > 0.02 
	    && refAmp < 0.49 && detAmp < 0.49
	    //&& ch1QualityBit == 0 && ch2QualityBit == 0
	    ) {
	  double t0 = t1gausroot;
	  double t1 = t1gausroot;
	  if (refChannel == 1) t0 = t1gausroot;
	  else if(refChannel == 2) t0 = t2gausroot;
	  else if(refChannel == 3) t0 = t3gausroot;
	  else if(refChannel == 4) t0 = t4gausroot;
	  else if(refChannel == 5) t0 = t5gausroot;
	  else if(refChannel == 6) t0 = t6gausroot;
	  else if(refChannel == 7) t0 = t7gausroot;
	  else if(refChannel == 8) t0 = t8gausroot;
	  if (detChannel == 1) t1 = t1gausroot;
	  else if(detChannel == 2) t1 = t2gausroot;
	  else if(detChannel == 3) t1 = t3gausroot;
	  else if(detChannel == 4) t1 = t4gausroot;
	  else if(detChannel == 5) t1 = t5gausroot;
	  else if(detChannel == 6) t1 = t6gausroot;
	  else if(detChannel == 7) t1 = t7gausroot;
	  else if(detChannel == 8) t1 = t8gausroot;

	  dt->Fill(t1 - t0 + deltaTAcrossDRSBoard);
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
  tex->DrawLatex(0.45, 0.80, Form("#sigma/#mu = %.1f %s",100*fitter->GetParameter(2)/fitter->GetParameter(1),"%"));
  tex->DrawLatex(0.45, 0.85, Form("Mean = %.1f #pm %.1f %s",1000*fitter->GetParameter(1),fitter->GetParError(1),"mV"));
  tex->DrawLatex(0.15, 0.92, Form("Attenuation Factor = %.3f",scalefactor));
  
  c->SaveAs( Form("%s_amplitude.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_amplitude.pdf", plotname.c_str()) );

  return;

  //time resolution plot
  c = new TCanvas("c","c",600,600);
  
  dt->SetAxisRange(timefitmin,timefitmax,"X");
  dt->SetTitle("");
  dt->GetXaxis()->SetTitle("#Delta t [ns]");
  dt->GetYaxis()->SetTitle("Number of Events");
  dt->GetYaxis()->SetTitleOffset(1.5);
  dt->SetMaximum(1.2*dt->GetMaximum());
  dt->Draw();
  dt->SetStats(0);
  dt->Fit("gaus","","",timefitmin,timefitmax);
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
  float x_lead[4] = { 2, 4, 6, 12 };
  float xerr_lead[4] = { 0.25,0.25, 0.25, 0.25 };
  float y_lead[4] = {  0.27, 0.41, 0.32, 0.06 };
  float yerr_lead[4] = { 0.01, 0.01, 0.01, 0.01 };

  float x_tungsten[4] = { 2, 4, 6, 8 };
  float xerr_tungsten[4] = { 0.25,0.25, 0.25, 0.25 };
  float y_tungsten[4] = {  0.27, 0.44, 0.35, 0.24 };
  float yerr_tungsten[4] = { 0.01, 0.01, 0.01, 0.01 };


  TGraphErrors *graphLead = new TGraphErrors(4,x_lead,y_lead,xerr_lead,yerr_lead);
  graphLead->SetLineWidth(2);
  TGraphErrors *graphTungsten = new TGraphErrors(4,x_tungsten,y_tungsten,xerr_tungsten,yerr_tungsten);
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



void PlotTimeResolutionVsPixel() {

  //*************************************************************
  //Run 31 Beam at center of 4 pixels , 4X0 tungsten
  //*************************************************************

  TH2F *TimeResolutionVsPixelGraphRun31 = new TH2F("TimeResolutionVsPixelGraphRun31",";Detector Channel; Reference Channel", 5, -0.5, 4.5, 5, -0.5, 4.5);
  TimeResolutionVsPixelGraphRun31->GetXaxis()->SetBinLabel(1,"Ref");
  TimeResolutionVsPixelGraphRun31->GetXaxis()->SetBinLabel(2,"Pixel 44");
  TimeResolutionVsPixelGraphRun31->GetXaxis()->SetBinLabel(3,"Pixel 54");
  TimeResolutionVsPixelGraphRun31->GetXaxis()->SetBinLabel(4,"Pixel 43");
  TimeResolutionVsPixelGraphRun31->GetXaxis()->SetBinLabel(5,"Pixel 53");
  TimeResolutionVsPixelGraphRun31->GetYaxis()->SetBinLabel(5,"Ref");
  TimeResolutionVsPixelGraphRun31->GetYaxis()->SetBinLabel(4,"Pixel 44");
  TimeResolutionVsPixelGraphRun31->GetYaxis()->SetBinLabel(3,"Pixel 54");
  TimeResolutionVsPixelGraphRun31->GetYaxis()->SetBinLabel(2,"Pixel 43");
  TimeResolutionVsPixelGraphRun31->GetYaxis()->SetBinLabel(1,"Pixel 53");


 
  TimeResolutionVsPixelGraphRun31->SetStats(0);
  TimeResolutionVsPixelGraphRun31->GetYaxis()->SetTitleOffset(2.0);

  TimeResolutionVsPixelGraphRun31->SetBinContent(2,5,29); TimeResolutionVsPixelGraphRun31->SetBinError(2,5,1);
  TimeResolutionVsPixelGraphRun31->SetBinContent(3,5,26); TimeResolutionVsPixelGraphRun31->SetBinError(3,5,1);
  TimeResolutionVsPixelGraphRun31->SetBinContent(4,5,33); TimeResolutionVsPixelGraphRun31->SetBinError(4,5,2);
  TimeResolutionVsPixelGraphRun31->SetBinContent(5,5,36); TimeResolutionVsPixelGraphRun31->SetBinError(5,5,2);
  TimeResolutionVsPixelGraphRun31->SetBinContent(3,4,27); TimeResolutionVsPixelGraphRun31->SetBinError(3,4,1);
  TimeResolutionVsPixelGraphRun31->SetBinContent(4,4,45); TimeResolutionVsPixelGraphRun31->SetBinError(4,4,2);
  TimeResolutionVsPixelGraphRun31->SetBinContent(5,4,46); TimeResolutionVsPixelGraphRun31->SetBinError(5,4,2);
  TimeResolutionVsPixelGraphRun31->SetBinContent(4,3,38); TimeResolutionVsPixelGraphRun31->SetBinError(4,3,2);
  TimeResolutionVsPixelGraphRun31->SetBinContent(5,3,41); TimeResolutionVsPixelGraphRun31->SetBinError(5,3,2);
  TimeResolutionVsPixelGraphRun31->SetBinContent(5,2,34); TimeResolutionVsPixelGraphRun31->SetBinError(5,2,2);

  TimeResolutionVsPixelGraphRun31->GetZaxis()->SetRangeUser(20,50);

  TimeResolutionVsPixelGraphRun31->SetMarkerSize(2.0);

  TCanvas *c = 0;
  c = new TCanvas("c","c",800,600);
  c->SetLeftMargin(0.15);

  TimeResolutionVsPixelGraphRun31->Draw("colztexte");
  TimeResolutionVsPixelGraphRun31->GetYaxis()->SetTitleOffset(1.70);
  TimeResolutionVsPixelGraphRun31->GetYaxis()->SetTitleSize(0.045);
  TimeResolutionVsPixelGraphRun31->GetYaxis()->SetLabelSize(0.060);
  TimeResolutionVsPixelGraphRun31->GetXaxis()->SetTitleOffset(1.02);
  TimeResolutionVsPixelGraphRun31->GetXaxis()->SetTitleSize(0.045);
  TimeResolutionVsPixelGraphRun31->GetXaxis()->SetLabelSize(0.060);
  TimeResolutionVsPixelGraphRun31->GetZaxis()->SetLabelSize(0.050);


  c->SaveAs( "TimeResolutionVsPixelGraph_Run31.gif");
  c->SaveAs( "TimeResolutionVsPixelGraph_Run31.pdf");

  return;

  //*************************************************************
  //Run 32,33 Beam at center of pixel 54 , 4X0 tungsten
  //*************************************************************
  TH2F *TimeResolutionVsPixelGraphRun32And33 = new TH2F("TimeResolutionVsPixelGraphRun32And33",";Detector Channel; Reference Channel", 5, -0.5, 4.5, 5, -0.5, 4.5);

  TimeResolutionVsPixelGraphRun32And33->GetXaxis()->SetBinLabel(1,"Ref");
  TimeResolutionVsPixelGraphRun32And33->GetXaxis()->SetBinLabel(2,"Pixel 44");
  TimeResolutionVsPixelGraphRun32And33->GetXaxis()->SetBinLabel(3,"Pixel 54");
  TimeResolutionVsPixelGraphRun32And33->GetXaxis()->SetBinLabel(4,"Pixel 43");
  TimeResolutionVsPixelGraphRun32And33->GetXaxis()->SetBinLabel(5,"Pixel 53");
  TimeResolutionVsPixelGraphRun32And33->GetYaxis()->SetBinLabel(5,"Ref");
  TimeResolutionVsPixelGraphRun32And33->GetYaxis()->SetBinLabel(4,"Pixel 44");
  TimeResolutionVsPixelGraphRun32And33->GetYaxis()->SetBinLabel(3,"Pixel 54");
  TimeResolutionVsPixelGraphRun32And33->GetYaxis()->SetBinLabel(2,"Pixel 43");
  TimeResolutionVsPixelGraphRun32And33->GetYaxis()->SetBinLabel(1,"Pixel 53");
  TimeResolutionVsPixelGraphRun32And33->SetStats(0);
  TimeResolutionVsPixelGraphRun32And33->GetYaxis()->SetTitleOffset(2.0);

  TimeResolutionVsPixelGraphRun32And33->SetBinContent(2,5,67); TimeResolutionVsPixelGraphRun32And33->SetBinError(2,5,3);
  TimeResolutionVsPixelGraphRun32And33->SetBinContent(3,5,37); TimeResolutionVsPixelGraphRun32And33->SetBinError(3,5,1);
  TimeResolutionVsPixelGraphRun32And33->SetBinContent(4,5,90); TimeResolutionVsPixelGraphRun32And33->SetBinError(4,5,3);
  TimeResolutionVsPixelGraphRun32And33->SetBinContent(5,5,72); TimeResolutionVsPixelGraphRun32And33->SetBinError(5,5,5);
  TimeResolutionVsPixelGraphRun32And33->SetBinContent(3,4,90); TimeResolutionVsPixelGraphRun32And33->SetBinError(3,4,4);
  TimeResolutionVsPixelGraphRun32And33->SetBinContent(4,4,99); TimeResolutionVsPixelGraphRun32And33->SetBinError(4,4,5);
  TimeResolutionVsPixelGraphRun32And33->SetBinContent(5,4,93); TimeResolutionVsPixelGraphRun32And33->SetBinError(5,4,4);
  TimeResolutionVsPixelGraphRun32And33->SetBinContent(4,3,96); TimeResolutionVsPixelGraphRun32And33->SetBinError(4,3,4);
  TimeResolutionVsPixelGraphRun32And33->SetBinContent(5,3,80); TimeResolutionVsPixelGraphRun32And33->SetBinError(5,3,5);
  TimeResolutionVsPixelGraphRun32And33->SetBinContent(5,2,91); TimeResolutionVsPixelGraphRun32And33->SetBinError(5,2,4);

  TimeResolutionVsPixelGraphRun32And33->GetZaxis()->SetRangeUser(0,100);

  TimeResolutionVsPixelGraphRun32And33->SetMarkerSize(2.0);

  TCanvas *c = 0;
  c = new TCanvas("c","c",800,600);
  c->SetLeftMargin(0.15);

  TimeResolutionVsPixelGraphRun32And33->Draw("colztexte");

  c->SaveAs( "TimeResolutionVsPixelGraph_Run32And33.gif");
  c->SaveAs( "TimeResolutionVsPixelGraph_Run32And33.pdf");


  //*************************************************************
  //Run 34 Beam at center of pixel 54 , 2X0 tungsten
  //*************************************************************
  TH2F *TimeResolutionVsPixelGraphRun34 = new TH2F("TimeResolutionVsPixelGraphRun34",";Detector Channel; Reference Channel", 5, -0.5, 4.5, 5, -0.5, 4.5);

  TimeResolutionVsPixelGraphRun34->GetXaxis()->SetBinLabel(1,"Ref");
  TimeResolutionVsPixelGraphRun34->GetXaxis()->SetBinLabel(2,"Pixel 44");
  TimeResolutionVsPixelGraphRun34->GetXaxis()->SetBinLabel(3,"Pixel 54");
  TimeResolutionVsPixelGraphRun34->GetXaxis()->SetBinLabel(4,"Pixel 43");
  TimeResolutionVsPixelGraphRun34->GetXaxis()->SetBinLabel(5,"Pixel 53");
  TimeResolutionVsPixelGraphRun34->GetYaxis()->SetBinLabel(5,"Ref");
  TimeResolutionVsPixelGraphRun34->GetYaxis()->SetBinLabel(4,"Pixel 44");
  TimeResolutionVsPixelGraphRun34->GetYaxis()->SetBinLabel(3,"Pixel 54");
  TimeResolutionVsPixelGraphRun34->GetYaxis()->SetBinLabel(2,"Pixel 43");
  TimeResolutionVsPixelGraphRun34->GetYaxis()->SetBinLabel(1,"Pixel 53");
  TimeResolutionVsPixelGraphRun34->SetStats(0);
  TimeResolutionVsPixelGraphRun34->GetYaxis()->SetTitleOffset(2.0);

  TimeResolutionVsPixelGraphRun34->SetBinContent(2,5,92); TimeResolutionVsPixelGraphRun34->SetBinError(2,5,7);
  TimeResolutionVsPixelGraphRun34->SetBinContent(3,5,54); TimeResolutionVsPixelGraphRun34->SetBinError(3,5,2);
  TimeResolutionVsPixelGraphRun34->SetBinContent(4,5,97); TimeResolutionVsPixelGraphRun34->SetBinError(4,5,5);
  TimeResolutionVsPixelGraphRun34->SetBinContent(5,5,78); TimeResolutionVsPixelGraphRun34->SetBinError(5,5,23);
  TimeResolutionVsPixelGraphRun34->SetBinContent(3,4,118); TimeResolutionVsPixelGraphRun34->SetBinError(3,4,7);
  TimeResolutionVsPixelGraphRun34->SetBinContent(4,4,124); TimeResolutionVsPixelGraphRun34->SetBinError(4,4,9);
  TimeResolutionVsPixelGraphRun34->SetBinContent(5,4,133); TimeResolutionVsPixelGraphRun34->SetBinError(5,4,10);
  TimeResolutionVsPixelGraphRun34->SetBinContent(4,3,125); TimeResolutionVsPixelGraphRun34->SetBinError(4,3,7);
  TimeResolutionVsPixelGraphRun34->SetBinContent(5,3,125); TimeResolutionVsPixelGraphRun34->SetBinError(5,3,13);
  TimeResolutionVsPixelGraphRun34->SetBinContent(5,2,118); TimeResolutionVsPixelGraphRun34->SetBinError(5,2,8);

  TimeResolutionVsPixelGraphRun34->GetZaxis()->SetRangeUser(0,135);

  TimeResolutionVsPixelGraphRun34->SetMarkerSize(2.0);

  TCanvas *c = 0;
  c = new TCanvas("c","c",800,600);
  c->SetLeftMargin(0.15);

  TimeResolutionVsPixelGraphRun34->Draw("colztexte");

  c->SaveAs( "TimeResolutionVsPixelGraph_Run34.gif");
  c->SaveAs( "TimeResolutionVsPixelGraph_Run34.pdf");




}


void PlotAmplitudeVsPixel() {

  //*************************************************************
  //Run 31 Beam at center of 4 pixels , 4X0 tungsten
  //*************************************************************
  TH2F *AmplitudeVsPixelGraphRun31 = new TH2F("AmplitudeVsPixelGraphRun31",";X Axis [mm]; Y Axis [mm]", 2, -0.5, 1.5, 2, -0.5, 1.5);
  AmplitudeVsPixelGraphRun31->GetXaxis()->SetBinLabel(1,"");
  AmplitudeVsPixelGraphRun31->GetXaxis()->SetBinLabel(2,"");
  AmplitudeVsPixelGraphRun31->GetYaxis()->SetBinLabel(1,"");
  AmplitudeVsPixelGraphRun31->GetYaxis()->SetBinLabel(2,"");
  AmplitudeVsPixelGraphRun31->SetStats(0);
  AmplitudeVsPixelGraphRun31->GetYaxis()->SetTitle("Y Axis");
  AmplitudeVsPixelGraphRun31->GetYaxis()->SetTitleOffset(1.0);
  AmplitudeVsPixelGraphRun31->GetXaxis()->SetTitle("X Axis");
  AmplitudeVsPixelGraphRun31->GetZaxis()->SetTitleOffset(1.25);
  AmplitudeVsPixelGraphRun31->GetZaxis()->SetTitle("Amplitude [mV]");

  AmplitudeVsPixelGraphRun31->SetBinContent(1,1,101); AmplitudeVsPixelGraphRun31->SetBinError(1,1,10);
  AmplitudeVsPixelGraphRun31->SetBinContent(1,2,78); AmplitudeVsPixelGraphRun31->SetBinError(1,2,10);
  AmplitudeVsPixelGraphRun31->SetBinContent(2,1,125); AmplitudeVsPixelGraphRun31->SetBinError(2,1,10);
  AmplitudeVsPixelGraphRun31->SetBinContent(2,2,108); AmplitudeVsPixelGraphRun31->SetBinError(2,2,10);

  AmplitudeVsPixelGraphRun31->GetZaxis()->SetRangeUser(0,150);

  AmplitudeVsPixelGraphRun31->SetMarkerSize(2.0);

  TCanvas *c = 0;
  c = new TCanvas("c","c",800,600);
  c->SetRightMargin(0.15);

  c->Range(-10,-1,10,1);
  // TGaxis *xaxis = new TGaxis(-4.5,-0.2,5.5,-0.2,-6,6,100,"");

  AmplitudeVsPixelGraphRun31->Draw("colztextee");
  AmplitudeVsPixelGraphRun31->GetYaxis()->SetTitleOffset(0.6);
  AmplitudeVsPixelGraphRun31->GetYaxis()->SetTitleSize(0.055);
  AmplitudeVsPixelGraphRun31->GetYaxis()->SetLabelSize(0.060);
  AmplitudeVsPixelGraphRun31->GetXaxis()->SetTitleOffset(0.6);
  AmplitudeVsPixelGraphRun31->GetXaxis()->SetTitleSize(0.055);
  AmplitudeVsPixelGraphRun31->GetXaxis()->SetLabelSize(0.060);
  AmplitudeVsPixelGraphRun31->GetZaxis()->SetLabelSize(0.045);
  AmplitudeVsPixelGraphRun31->GetZaxis()->SetTitleSize(0.045);
  AmplitudeVsPixelGraphRun31->GetZaxis()->SetTitleOffset(1.05);

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.25, 0.80, "Pixel 44");
  tex->DrawLatex(0.65, 0.80, "Pixel 54");
  tex->DrawLatex(0.25, 0.40, "Pixel 43");
  tex->DrawLatex(0.65, 0.40, "Pixel 53");

  //axis->SetLabelColor(kBlack);
  // xaxis->Draw();

  c->SaveAs( "AmplitudeVsPixelGraph_Run31.gif");
  c->SaveAs( "AmplitudeVsPixelGraph_Run31.pdf");

  return;

  //*************************************************************
  //Run 32,33 Beam at center of pixel 54 , 4X0 tungsten
  //*************************************************************

  TH2F *AmplitudeVsPixelGraphRun32And33 = new TH2F("AmplitudeVsPixelGraphRun32And33",";X Axis [mm]; Y Axis [mm]", 2, -0.5, 1.5, 2, -0.5, 1.5);
  AmplitudeVsPixelGraphRun32And33->GetXaxis()->SetBinLabel(1,"");
  AmplitudeVsPixelGraphRun32And33->GetXaxis()->SetBinLabel(2,"");
  AmplitudeVsPixelGraphRun32And33->GetYaxis()->SetBinLabel(1,"");
  AmplitudeVsPixelGraphRun32And33->GetYaxis()->SetBinLabel(2,"");
  AmplitudeVsPixelGraphRun32And33->SetStats(0);
  AmplitudeVsPixelGraphRun32And33->GetYaxis()->SetTitle("Y Axis");
  AmplitudeVsPixelGraphRun32And33->GetYaxis()->SetTitleOffset(1.0);
  AmplitudeVsPixelGraphRun32And33->GetXaxis()->SetTitle("X Axis");
  AmplitudeVsPixelGraphRun32And33->GetZaxis()->SetTitleOffset(1.25);
  AmplitudeVsPixelGraphRun32And33->GetZaxis()->SetTitle("Amplitude [mV}");

  AmplitudeVsPixelGraphRun32And33->SetBinContent(1,1,38); AmplitudeVsPixelGraphRun32And33->SetBinError(1,1,10);
  AmplitudeVsPixelGraphRun32And33->SetBinContent(1,2,35); AmplitudeVsPixelGraphRun32And33->SetBinError(1,2,10);
  AmplitudeVsPixelGraphRun32And33->SetBinContent(2,1,57); AmplitudeVsPixelGraphRun32And33->SetBinError(2,1,10);
  AmplitudeVsPixelGraphRun32And33->SetBinContent(2,2,202); AmplitudeVsPixelGraphRun32And33->SetBinError(2,2,10);

  AmplitudeVsPixelGraphRun32And33->GetZaxis()->SetRangeUser(0,250);

  AmplitudeVsPixelGraphRun32And33->SetMarkerSize(2.0);

  TCanvas *c = 0;
  c = new TCanvas("c","c",800,600);
  c->SetRightMargin(0.15);

  c->Range(-10,-1,10,1);
  // TGaxis *xaxis = new TGaxis(-4.5,-0.2,5.5,-0.2,-6,6,100,"");

  AmplitudeVsPixelGraphRun32And33->Draw("colztextee");

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.25, 0.80, "Pixel 44");
  tex->DrawLatex(0.65, 0.80, "Pixel 54");
  tex->DrawLatex(0.25, 0.40, "Pixel 43");
  tex->DrawLatex(0.65, 0.40, "Pixel 53");

  //axis->SetLabelColor(kBlack);
  // xaxis->Draw();

  c->SaveAs( "AmplitudeVsPixelGraph_Run32And33.gif");
  c->SaveAs( "AmplitudeVsPixelGraph_Run32And33.pdf");


  //*************************************************************
  //Run 34 Beam at center of pixel 54 , 2X0 tungsten
  //*************************************************************

  TH2F *AmplitudeVsPixelGraphRun34 = new TH2F("AmplitudeVsPixelGraphRun34",";X Axis [mm]; Y Axis [mm]", 2, -0.5, 1.5, 2, -0.5, 1.5);
  AmplitudeVsPixelGraphRun34->GetXaxis()->SetBinLabel(1,"");
  AmplitudeVsPixelGraphRun34->GetXaxis()->SetBinLabel(2,"");
  AmplitudeVsPixelGraphRun34->GetYaxis()->SetBinLabel(1,"");
  AmplitudeVsPixelGraphRun34->GetYaxis()->SetBinLabel(2,"");
  AmplitudeVsPixelGraphRun34->SetStats(0);
  AmplitudeVsPixelGraphRun34->GetYaxis()->SetTitle("Y Axis");
  AmplitudeVsPixelGraphRun34->GetYaxis()->SetTitleOffset(1.0);
  AmplitudeVsPixelGraphRun34->GetXaxis()->SetTitle("X Axis");
  AmplitudeVsPixelGraphRun34->GetZaxis()->SetTitleOffset(1.25);
  AmplitudeVsPixelGraphRun34->GetZaxis()->SetTitle("Amplitude [mV}");

  AmplitudeVsPixelGraphRun34->SetBinContent(1,2,37); AmplitudeVsPixelGraphRun34->SetBinError(1,2,10);
  AmplitudeVsPixelGraphRun34->SetBinContent(2,1,252); AmplitudeVsPixelGraphRun34->SetBinError(2,1,10);
  AmplitudeVsPixelGraphRun34->SetBinContent(1,1,36); AmplitudeVsPixelGraphRun34->SetBinError(1,1,10);
  AmplitudeVsPixelGraphRun34->SetBinContent(2,2,43); AmplitudeVsPixelGraphRun34->SetBinError(2,2,10);

  AmplitudeVsPixelGraphRun34->GetZaxis()->SetRangeUser(0,250);

  AmplitudeVsPixelGraphRun34->SetMarkerSize(2.0);

  TCanvas *c = 0;
  c = new TCanvas("c","c",800,600);
  c->SetRightMargin(0.15);

  c->Range(-10,-1,10,1);
  // TGaxis *xaxis = new TGaxis(-4.5,-0.2,5.5,-0.2,-6,6,100,"");

  AmplitudeVsPixelGraphRun34->Draw("colztextee");

  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.050);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.25, 0.80, "Pixel 44");
  tex->DrawLatex(0.65, 0.80, "Pixel 54");
  tex->DrawLatex(0.25, 0.40, "Pixel 43");
  tex->DrawLatex(0.65, 0.40, "Pixel 53");

  //axis->SetLabelColor(kBlack);
  // xaxis->Draw();

  c->SaveAs( "AmplitudeVsPixelGraph_Run34.gif");
  c->SaveAs( "AmplitudeVsPixelGraph_Run34.pdf");



}





void PixelTimingAnalysis() {

  //*************************************
  // Run 31 : Beam at Center of 4 pixels , 4X0 tungsten
  //*************************************
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_31.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 1, 2, 1.0, 0.02, 0.2, 1.7 , 2.2);
  // MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_31.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 1, 3, 1.0, 0.02, 0.2, 1.7 , 2.2);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_31.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 5, 6, 1.0, 0.02, 0.2, 1.9 , 2.4);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_31.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 5, 7, 1.0, 0.02, 0.2, 1.9 , 2.4);
  
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_31.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 2, 3, 1.0, 0.02, 0.2, -0.2 , 0.2);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_31.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 2, 6, 1.0, 0.02, 0.2, -0.1 , 0.5);
  // MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_31.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 2, 7, 1.0, 0.02, 0.2, -0.1 , 0.5);
  // MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_31.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 3, 6, 1.0, 0.02, 0.2, 0.1 , 0.5);
  // MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_31.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 3, 7, 1.0, 0.02, 0.2, 0.0 , 0.4);
  // MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_31.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 6, 7, 1.0, 0.02, 0.2, -0.3 , 0.1);


  //*************************************
  // Run 32 : Beam at Center of pixel 54 , 4X0 tungsten
  //*************************************
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_32And33.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 1, 2, 1.0, 0.00, 0.2, 1.7 , 2.2);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_32And33.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 1, 3, 1.0, 0.02, 0.4, 1.7 , 2.2);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_32And33.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 5, 6, 1.0, 0.0, 0.2, 1.5 , 2.4);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_32And33.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 5, 7, 1.0, 0.02, 0.2, 1.5 , 2.4);
  
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_32And33.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 2, 3, 1.0, 0.02, 0.2, -0.3 , 0.3);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_32And33.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 2, 6, 1.0, 0.02, 0.2, -0.1 , 0.5);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_32And33.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 2, 7, 1.0, 0.02, 0.2, -0.1 , 0.5);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_32And33.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 3, 6, 1.0, 0.02, 0.2, -0.2 , 0.5);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_32And33.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 3, 7, 1.0, 0.02, 0.2, -0.2 , 0.4);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_32And33.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 6, 7, 1.0, 0.02, 0.2, -0.3 , 0.4);

  //*************************************
  // Run 34 : Beam at Center of pixel 54 , 2X0 tungsten
  //*************************************
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_34.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 1, 2, 1.0, 0.00, 0.2, 1.5 , 2.2);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_34.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 1, 3, 1.0, 0.02, 0.4, 1.5 , 2.2);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_34.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 5, 6, 1.0, 0.0, 0.2, 1.5 , 2.4);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_34.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 5, 7, 1.0, 0.02, 0.2, 1.5 , 2.4);
  
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_34.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 2, 3, 1.0, 0.02, 0.2, -0.5 , 0.5);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_34.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 2, 6, 1.0, 0.02, 0.2, -0.5 , 0.5);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_34.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 2, 7, 1.0, 0.02, 0.2, -0.5 , 0.5);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_34.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 3, 6, 1.0, 0.02, 0.2, -0.5 , 0.6);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_34.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 3, 7, 1.0, 0.02, 0.2, -0.5 , 0.6);
  //MakePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_34.root.ana_lowpassfilter.root","Electron_4X0Tungsten_8GeV_PhotonisPixel_RefToPixel44", 6, 7, 1.0, 0.02, 0.2, -0.5 , 0.6);

  PlotAmplitudeVsPixel();
  //PlotTimeResolutionVsPixel();






  //
  //*************************************
  //Results for Y11 Fibers from May Testbeam
  //*************************************
//   MakeTimeResolutionPlot("cpt_may_run_131.ana.root","TOF_ShashlikY11Fiber_Electron_4GeV",4,false);
//   MakeTimeResolutionPlot("cpt_may_run_132.ana.root","TOF_ShashlikY11Fiber_Electron_8GeV",8,false);
//   MakeTimeResolutionPlot("cpt_may_run_133.ana.root","TOF_ShashlikY11Fiber_Electron_16GeV",16,false);
//   MakeTimeResolutionVsEnergyPlot_Y11();




}
