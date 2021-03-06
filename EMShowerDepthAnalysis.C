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
  dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 400, -1,1);
  
  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();
  std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
      tree->GetEntry(iEntry);    
      
      //require cherenkov and signal in the front MCP
      if (ch1Amp > 0.02 && ch4Amp > 0.1) {
	histAmplitude->Fill(ch2Amp);
	if (ch2Amp > 0.02 
	    && ch1Amp < 0.49 && ch2Amp < 0.49
	    //&& ch1QualityBit == 0 && ch2QualityBit == 0
	    ) {
	  dt->Fill(t2gausroot - t1gausroot);
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
  //tex->DrawLatex(0.15, 0.92, Form("Attenuation Factor = %.3f",scalefactor));
  
  c->SaveAs( Form("%s_amplitude.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_amplitude.pdf", plotname.c_str()) );


  //time resolution plot
  c = new TCanvas("c","c",600,600);
  
  dt->SetAxisRange(0.1,0.6,"X");
  dt->SetTitle("");
  dt->GetXaxis()->SetTitle("#Delta t [ns]");
  dt->GetYaxis()->SetTitle("Number of Events");
  dt->GetYaxis()->SetTitleOffset(1.3);
  dt->SetMaximum(1.2*dt->GetMaximum());
  dt->Draw();
  dt->SetStats(0);
  dt->Fit("gaus","","",0.15,0.45);
  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.040);
  tex->SetTextFont(42);
  tex->SetTextColor(kBlack);
  tex->DrawLatex(0.55, 0.80, Form("#sigma = %.1f #pm %.1f ps",1000*fitter->GetParameter(2),1000*fitter->GetParError(2)));
  //tex->DrawLatex(0.55, 0.85, Form("Mean = %.2f #pm %.2f %s",fitter->GetParameter(1),fitter->GetParError(1),"V"));
  tex->DrawLatex(0.55, 0.85, Form("Mean = %.1f ps",fitter->GetParameter(1)));
  
  c->SaveAs( Form("%s_dt.gif", plotname.c_str()) );
  c->SaveAs( Form("%s_dt.pdf", plotname.c_str()) );


  
}



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
  dt = new TH1F("dt","; #Delta t [ns]; Number of Events", 400, -1,1);
  
  //read all entries and fill the histograms
  Long64_t nentries = tree->GetEntries();
  std::cout<<"Number of events in Physics Sample: "<<nentries<<std::endl;  
  for (Long64_t iEntry=0;iEntry<nentries;iEntry++) {
      tree->GetEntry(iEntry);    
      
      //require cherenkov and signal in the front MCP
      if (ch1Amp > 20 && ch4Amp > 100) {
	histAmplitude->Fill(ch2Amp/1000);
	if (ch2Amp > 20 
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
  dt->GetXaxis()->SetTitleOffset(1.00);
  dt->GetYaxis()->SetTitle("Number of Events");
  dt->GetYaxis()->SetTitleOffset(1.02);
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
  tex->DrawLatex(0.50, 0.80, Form("#sigma = %.1f #pm %.1f ps",1000*fitter->GetParameter(2),1000*fitter->GetParError(2)));
  //tex->DrawLatex(0.55, 0.85, Form("Mean = %.2f #pm %.2f %s",fitter->GetParameter(1),fitter->GetParError(1),"V"));
  tex->DrawLatex(0.50, 0.85, Form("Mean = %.2f ns",fitter->GetParameter(1)));
  
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
  float x_lead_SetOne[4] = { 2.5, 4.5, 6.5, 12.5 };
  float xerr_lead_SetOne[4] = { 0.25,0.25, 0.25, 0.25 };
  float y_lead_SetOne[4] = {  0.27, 0.41, 0.32, 0.06 }; //old numbers  
  float yerr_lead_SetOne[4] = { 0.01, 0.01, 0.01, 0.01 };

  float x_tungsten_SetOne[4] = { 2.5, 4.5, 6.5, 8.5 };
  float xerr_tungsten_SetOne[4] = { 0.25,0.25, 0.25, 0.25 };
  float y_tungsten_SetOne[4] = {  0.27, 0.44, 0.35, 0.24 };
  float yerr_tungsten_SetOne[4] = { 0.01, 0.01, 0.01, 0.01 };

  float x_lead_SetTwo[3] = { 6.5, 8.5, 10.5 };
  float xerr_lead_SetTwo[3] = { 0.25,0.25, 0.25 };
  //float y_lead_SetTwo[3] = {  0.32, 0.17, 0.11  }; //with low pass filter
  float y_lead_SetTwo[3] = {  0.474, 0.316, 0.16  }; //no low pass filter  
  float yerr_lead_SetTwo[3] = { 0.01, 0.01, 0.01 };

  float x_tungsten_SetTwo[3] = { 4.5, 10.5, 12.5 };
  float xerr_tungsten_SetTwo[3] = { 0.25, 0.25, 0.25 };
  //float y_tungsten_SetTwo[3] = {  0.35, 0.15, 0.06}; //with low pass filter
  float y_tungsten_SetTwo[3] = {  0.515, 0.22, 0.09}; //no low pass filter
  float yerr_tungsten_SetTwo[3] = { 0.01, 0.01, 0.01 };

  //By looking at the amplitude of CH1, we can get the average amplitude for
  //the point with 0 absorber. it is about 0.026 +- 0.005.

  //we normalized SetTwo to SetOne in 2 different situations.
  //For 2X0 tungsten, we get a scaling factor of 0.44/0.515 to go from Set2 to Set1
  //For 4X0 lead, we got a scaling factor of 0.32/0.474 to go from Set2 to Set1
  
  //If we use the average scaling factor from both runs, we get the following results
  //scaling factor = 0.765
  float x_lead_combined[7] = { 0.1, 2.4, 4.5, 6.5, 8.4, 10.5, 12.5 };
  float xerr_lead_combined[7] = { 0.25, 0.25,0.25, 0.25, 0.25, 0.25, 0.25 };
  float y_lead_combined[7] = {  0.026, 0.27, 0.41, 0.34, 0.24, 0.12, 0.06  };
  float yerr_lead_combined[7] = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };
  float x_tungsten_combined[7] = { 0, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5 };
  float xerr_tungsten_combined[7] = { 0.25, 0.25, 0.25, 0.25, 0.25 , 0.25, 0.25};
  float y_tungsten_combined[7] = {  0.026, 0.27, 0.42, 0.35, 0.24, 0.17, 0.07 };
  float yerr_tungsten_combined[7] = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01 };


  TGraphErrors *graphLead = new TGraphErrors(7,x_lead_combined,y_lead_combined,xerr_lead_combined,yerr_lead_combined);
  graphLead->SetLineWidth(3);
  TGraphErrors *graphTungsten = new TGraphErrors(7,x_tungsten_combined,y_tungsten_combined,xerr_tungsten_combined,yerr_tungsten_combined);
  graphTungsten->SetLineWidth(3);
  graphTungsten->SetLineColor(kBlue);

  TCanvas * c = new TCanvas("c","c",800,600);
  graphLead->Draw("AP");
  graphTungsten->Draw("Psame");

  TLegend *legend = new TLegend (0.50,0.7,0.8,0.85);
  legend->SetTextSize(0.05);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(graphLead, "Lead Absorber","LP");
  legend->AddEntry(graphTungsten, "Tungsten Absorber","LP");
  legend->Draw();

  graphLead->SetTitle("");
  graphLead->GetXaxis()->SetTitle("Absorber Thickness [X_{0}]");
  graphLead->GetXaxis()->SetTitleSize(0.045);
  graphLead->GetXaxis()->SetLabelSize(0.045);
  graphLead->GetXaxis()->SetTitleOffset(1.0);
  graphLead->GetYaxis()->SetTitle("Mean Amplitude [V]");
  graphLead->GetYaxis()->SetTitleOffset(1.02);
  graphLead->GetYaxis()->SetTitleSize(0.045);
  graphLead->GetYaxis()->SetLabelSize(0.045);
  graphLead->GetYaxis()->SetRangeUser(0,0.6);

  c->SaveAs( "AmplitudeVsLeadThickness.gif" );
  c->SaveAs( "AmplitudeVsLeadThickness.pdf" );

}


void MakeTimeResolutionVsShowerDepthGraph() {

  //use beam energy for xaxis
  // float x_lead[4] = { 2, 4, 6, 12 };
  // float xerr_lead[4] = { 0.25,0.25, 0.25, 0.25 };
  //With old reconstruction
  // float y_lead[4] = {  18.9, 18.5, 19.3, 20.2 };
  // float yerr_lead[4] = { 0.7, 0.7, 0.7, 0.8 };
  //With improved reconstruction (spline + low pass filter)
  // float y_lead[4] = {  15.3, 15.9, 16.7 , 19.0  };
  // float yerr_lead[4] = { 0.5, 0.6, 0.6, 0.7  };
  //with low pass filter, no spline
  // float y_lead[4] = {  15.2, 14.9, 16.6 , 17.5  };
  // float yerr_lead[4] = { 0.5, 0.6, 0.6, 0.7  };
  // //with Heejong code
  // float y_lead[4] = { 13.8 , 13.7 , 15.9 , 17.3   };
  // float yerr_lead[4] = { 0.5, 0.6, 0.6, 0.7  };

  // //with Heejong code + some modifications
  // //For proton beam At CDF 75% amplitude we get 10.2ps +- 0.4ps resolution
  // float y_lead[4] = { 11.7 , 12.9  , 14.0 , 16.8    };
  // float yerr_lead[4] = { 0.4, 0.4 , 0.5 , 0.6  };


  // float x_tungsten[4] = { 2, 4, 6, 8 };
  // float xerr_tungsten[4] = { 0.25,0.25, 0.25, 0.25 };
  //With old reconstruction
  // float y_tungsten[4] = {  18.9, 19.5, 18.1, 18.3 };
  // float yerr_tungsten[4] = { 0.7, 0.7, 0.6, 0.5 };
  //With improved reconstruction (spline + low pass filter)
  // float y_tungsten[4] = {  15.3 , 15.9, 15.7 , 14.9  };
  // float yerr_tungsten[4] = { 0.5 , 0.6 , 0.5 , 0.5 };
  //With low pass filter
  // float y_tungsten[4] = {  15.3 , 15.1, 15.5 , 15.1  };
  // float yerr_tungsten[4] = { 0.5 , 0.6 , 0.5 , 0.5 };
  // //With Heejong code
  // float y_tungsten[4] = { 13.8 , 15.2 , 14.1 , 13.7   };
  // float yerr_tungsten[4] = { 0.5 , 0.6 , 0.6 , 0.6 };

  // //With Heejong code + some modifications
  // float y_tungsten[4] = { 11.7 , 14.2 , 14.3 , 12.1 };
  // float yerr_tungsten[4] = { 0.4 , 0.6 , 0.6 , 0.5 };

  //With Heejong code + some modifications with additional points
  float x_tungsten_SetOne[4] = { 2.5, 4.5, 6.5, 8.5 };
  float xerr_tungsten_SetOne[4] = { 0.25 , 0.25, 0.25, 0.25 };
  float y_tungsten_SetOne[4] = { 11.7 , 14.2 , 14.3 , 12.1 };
  float yerr_tungsten_SetOne[4] = { 0.4 , 0.6 , 0.6 , 0.5 };
  float x_tungsten_SetTwo[3] = { 4.5, 10.5, 12.5 };
  float xerr_tungsten_SetTwo[3] = { 0.25 , 0.25, 0.25 };
  float y_tungsten_SetTwo[3] = { 12.1 , 13.3, 13.6 };
  float yerr_tungsten_SetTwo[3] = { 0.5 , 0.5 , 0.5 };
  float x_tungsten_combined[6] = { 2.5, 4.5, 6.5, 8.5, 10.5, 12.5 };
  float xerr_tungsten_combined[6] = { 0.25 , 0.25, 0.25, 0.25, 0.25, 0.25 };
  float y_tungsten_combined[6] = { 11.7 , 13.2 , 14.3 , 12.1, 13.3, 13.6 };
  float yerr_tungsten_combined[6] = { 0.4 , 0.4 , 0.6 , 0.5, 0.5, 0.5 };

  float x_lead_SetOne[4] = { 2.5, 4.5, 6.5, 12.5 };
  float xerr_lead_SetOne[4] = { 0.25,0.25, 0.25, 0.25 };
  float y_lead_SetOne[4] = { 11.7 , 12.9  , 14.0 , 16.8 };
  float yerr_lead_SetOne[4] = { 0.4, 0.4 , 0.5 , 0.6  };
  float x_lead_SetTwo[3] = { 6.5, 8.5, 10.5 };
  float xerr_lead_SetTwo[3] = { 0.25, 0.25, 0.25 };
  float y_lead_SetTwo[3] = { 12.2, 12.9, 12.4 };
  float yerr_lead_SetTwo[3] = { 0.4, 0.5 , 0.4  };
  float x_lead_combined[6] = { 2.5, 4.5, 6.5, 8.5, 10.5, 12.5 };
  float xerr_lead_combined[6] = { 0.25,0.25, 0.25, 0.25, 0.25, 0.25 };
  float y_lead_combined[6] = { 11.7 , 12.9  , 13.1 , 12.9 , 12.4 , 16.8 };
  float yerr_lead_combined[6] = { 0.4, 0.4 , 0.3 , 0.5, 0.4, 0.6  };


  TGraphErrors *graphLead = new TGraphErrors(6,x_lead_combined,y_lead_combined,xerr_lead_combined,yerr_lead_combined);
  graphLead->SetLineWidth(3);
  TGraphErrors *graphTungsten = new TGraphErrors(6,x_tungsten_combined,y_tungsten_combined,xerr_tungsten_combined,yerr_tungsten_combined);
  graphTungsten->SetLineWidth(3);
  graphTungsten->SetLineColor(kBlue);

  TCanvas * c = new TCanvas("c","c",800,600);
  graphLead->Draw("AP");
  graphTungsten->Draw("Psame");

  TLegend *legend = new TLegend (0.30,0.7,0.8,0.85);
  legend->SetTextSize(0.05);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(graphLead, "Lead Absorber","LP");
  legend->AddEntry(graphTungsten, "Tungsten Absorber","LP");
  legend->Draw();

  graphLead->SetTitle("");
  graphLead->GetXaxis()->SetTitle("Absorber Thickness [X_{0}]");
  graphLead->GetXaxis()->SetTitleSize(0.045);
  graphLead->GetXaxis()->SetLabelSize(0.045);
  graphLead->GetXaxis()->SetTitleOffset(1.0);
  graphLead->GetYaxis()->SetTitle("Time Resolution [ps]");
  graphLead->GetYaxis()->SetTitleOffset(1.02);
  graphLead->GetYaxis()->SetTitleSize(0.045);
  graphLead->GetYaxis()->SetLabelSize(0.045);
  graphLead->GetYaxis()->SetRangeUser(5,20);

  c->SaveAs( "TimeResolutionVsLeadThickness.gif" );
  c->SaveAs( "TimeResolutionVsLeadThickness.pdf" );

}






void EMShowerDepthAnalysis() {

  //*************************************
  // Amplitude Vs Shower Depth
  //*************************************
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_4.ana.root","Amplitude_Electron_NoAbsorber_8GeV", 1.0, 0.1, 0.45);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_5.ana.root","Amplitude_Electron_2X0Lead_8GeV", 3.162, 0.05, 0.2);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_6.ana.root","Amplitude_Electron_4X0Load_8GeV",3.162, 0.02,0.45 );
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_7.ana.root","Amplitude_Electron_10X0Lead_8GeV", 1.0, 0.00, 0.1);

  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_8.ana.root","Amplitude_Electron_2X0Tungsten_8GeV", 3.162, 0.05, 0.2);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_9.ana.root","Amplitude_Electron_4X0Tungsten_8GeV",3.162, 0.02,0.45 );
  // MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_10.ana.root","Amplitude_Electron_6X0Tungsten_8GeV", 1.0, 0.00, 0.4);

  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_68.ana.root","Amplitude_Electron_4X0Lead_8GeV", 3.162, 0.05, 0.3);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_69.ana.root","Amplitude_Electron_6X0Lead_8GeV", 3.162, 0.02, 0.2);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_70.ana.root","Amplitude_Electron_6X0Lead_8GeV", 1.0, 0.02, 0.4);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_71.ana.root","Amplitude_Electron_8X0Lead_8GeV", 1.0, 0.02, 0.4);

  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_67.ana.root","Amplitude_Electron_2X0Tungsten_8GeV", 3.162, 0.05, 0.3);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_65.ana.root","Amplitude_Electron_8X0Tungsten_8GeV", 1.0, 0.02, 0.4);
  //MakeAmplitudePlot("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_66.ana.root","Amplitude_Electron_10X0Tungsten_8GeV", 1.0, 0.02, 0.4);


  //*************************************
  // Make plot from heejong analysis ntuple
  //*************************************
  MakePlotAlternativeFormat("/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/t1058_jan2015_run_8.ana_heejongCode.root","Amplitude_Electron_NoAbsorber_8GeV", 1.0, 0.1, 0.3);



  //*************************************
  // Final Graphs
  //*************************************
  //MakeAmplitudeVsShowerDepthGraph();
  //MakeTimeResolutionVsShowerDepthGraph();





}
