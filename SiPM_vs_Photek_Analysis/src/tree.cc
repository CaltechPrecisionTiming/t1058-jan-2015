#define tree_cxx
#include "tree.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TLatex.h>
#include <iostream>

const float lumi = 5;
//Axis
const float axisTitleSize = 0.06;
const float axisTitleOffset = .8;

const float axisTitleSizeRatioX   = 0.18;
const float axisLabelSizeRatioX   = 0.12;
const float axisTitleOffsetRatioX = 0.94;

const float axisTitleSizeRatioY   = 0.15;
const float axisLabelSizeRatioY   = 0.108;
const float axisTitleOffsetRatioY = 0.32;

//Margins
const float leftMargin   = 0.12;
const float rightMargin  = 0.05;
const float topMargin    = 0.07;
const float bottomMargin = 0.12;

//CMS STANDARD
TString CMSText = "CMS";
TString extraText   = "Preliminary";
//TString lumiText = "2.32 fb^{-1} (13 TeV)";
TString lumiText = "12.92 fb^{-1} (13 TeV)";

void tree::Loop(int amp, int Att, float ND)
{
 
//   In a ROOT session, you can do:
//      root> .L tree.C
//      root> tree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  int times = amp-Att;
  //V0 = Vx10^(-times/20)
  TH1F* ch1_Int = new TH1F("ch1Int", "ch1Int", 1000, 0, 10);
  TH1F* deltaT12 = new TH1F("deltaT12", "deltaT12", 100, -1.6, -1);
   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "n: " << nentries << std::endl;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      ch1_Int->Fill(ch1Int);
      //ch1_Int->Fill(ch1Int*TMath::Power(10,-times*1.0/20.1));
      deltaT12->Fill(t1gausroot-t2gausroot);
      // if (Cut(ientry) < 0) continue;
   }

   //---------------------
   //FileName
   //---------------------
   char* fname = "ND0p5";
   
   TCanvas* c = new TCanvas( "c", "c", 2119, 33, 800, 700 );
   c->SetHighLightColor(2);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetLeftMargin( leftMargin );
   c->SetRightMargin( 1.6*rightMargin );
   c->SetTopMargin( topMargin );
   c->SetBottomMargin( bottomMargin );
   c->SetFrameBorderMode(0);
   c->SetFrameBorderMode(0);
   
   double ch1AmpMean = ch1_Int->GetMean();
   double ch1AmpRMS = ch1_Int->GetRMS();
   ch1_Int->GetXaxis()->SetRangeUser(ch1AmpMean-5.*ch1AmpRMS,ch1AmpMean+5.*ch1AmpRMS);
   TF1* f1Int = new TF1("f1Int","gaus(0)",ch1AmpMean-2.*ch1AmpRMS,ch1AmpMean+2.*ch1AmpRMS);
   ch1_Int->Fit(f1Int,"REL");
   ch1_Int->SetTitle("");
   ch1_Int->SetStats(0);
   ch1_Int->SetYTitle("events/xx pC");
   ch1_Int->SetXTitle("charge [pC]");
   ch1_Int->SetMarkerColor(kBlack);
   ch1_Int->SetLineColor(kBlack);
   ch1_Int->SetMarkerStyle(20);
   ch1_Int->SetMarkerSize(1.2);
   ch1_Int->Draw("E");
   
   TLatex latex;
   latex.SetNDC();
   latex.SetTextAngle(0);
   latex.SetTextColor(kBlack);    
   latex.SetTextFont(42);
   latex.SetTextAlign(31); 
   latex.SetTextSize(0.04);    
   latex.DrawLatex(0.89, 0.85, Form("Q = %.2f #pm %.2f [pc]", f1Int->GetParameter(1), f1Int->GetParameter(2)));
   
   c->SaveAs(Form("chargeCh1_%s.pdf", fname));
   c->SaveAs(Form("chargeCh1_%s.png", fname));
   c->SaveAs(Form("chargeCh1_%s.C", fname));
   //---------------------
   //------deltaT---------
   //---------------------
   double deltaT12Mean = deltaT12->GetMean();
   double deltaT12RMS  = deltaT12->GetRMS();
   deltaT12->GetXaxis()->SetRangeUser(deltaT12Mean-5.*deltaT12RMS,deltaT12Mean+5.*deltaT12RMS);
   TF1* f1DeltaT = new TF1("deltaT","gaus(0)",deltaT12Mean-2.*deltaT12RMS,deltaT12Mean+2.*deltaT12RMS);
   deltaT12->Fit(f1DeltaT,"REL");
   deltaT12->SetTitle("");
   deltaT12->SetStats(0);
   deltaT12->SetYTitle("events/xx ns");
   deltaT12->SetXTitle("t_{1} - t_{0} [ns]");
   deltaT12->SetMarkerColor(kBlack);
   deltaT12->SetLineColor(kBlack);
   deltaT12->SetMarkerStyle(20);
   deltaT12->SetMarkerSize(1.2);
   deltaT12->Draw("E");
   latex.DrawLatex(0.89, 0.85, Form("#sigma_{t} = %.2f #pm %.2f [ps]", f1DeltaT->GetParameter(2)*1000., f1DeltaT->GetParError(2)*1000.));
   c->SaveAs(Form("deltaT12_%s.pdf", fname));
   c->SaveAs(Form("deltaT12_%s.png", fname));
   c->SaveAs(Form("deltaT12_%s.C", fname));
   
   
   TFile* fout = new TFile("myAnalysisFile.root", "recreate");
   f1Int->Write();
   ch1_Int->Write();
   f1DeltaT->Write();
   deltaT12->Write();
   fout->Close();
}
