#define tree_cxx
#include "tree.hh"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>

void tree::Loop()
{

  TCanvas* c = new TCanvas( "c", "c", 800, 600 );
  c->SetHighLightColor(2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F* h_time = new TH1F("h_time", "h_time", 300, 25., 27.);
   
   TH1F* h_amp  = new TH1F("h_amp", "h_amp", 1000, 0., .5);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if( !(ch1Amp > 0.02 && ch1Amp < 0.48 && ch3Amp > 0.01 && ch3Amp < 0.48) ) continue;
      h_time->Fill(ch3THM-ch1THM);
      h_amp->Fill(ch3Amp);
   }

   TF1* gaus1 = new TF1( "gaus1", "gaus(0)", h_time->GetMean()-2.*h_time->GetRMS() , h_time->GetMean()+2.*h_time->GetRMS() );
   h_time->Fit( gaus1, "RQM");
   h_time->GetXaxis()->SetRangeUser(h_time->GetMean()-5.*h_time->GetRMS() , h_time->GetMean()+5.*h_time->GetRMS());
   h_time->Draw();
   c->SaveAs("h_time.pdf");

   
   TF1* gaus2 = new TF1( "gaus2", "gaus(0)", h_amp->GetMean()-2.*h_amp->GetRMS() , h_amp->GetMean()+2.*h_amp->GetRMS() );
   h_amp->Fit( gaus2, "RQM");
   h_amp->GetXaxis()->SetRangeUser(h_amp->GetMean()-5.*h_amp->GetRMS() , h_amp->GetMean()+5.*h_amp->GetRMS());
   h_amp->Draw();
   c->SaveAs("h_amp.pdf");
   
   
}

std::pair<float,float> tree::GetAmp( TString pos )
{

  TCanvas* c = new TCanvas( "c", "c", 800, 600 );
  c->SetHighLightColor(2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  
   if (fChain == 0) return std::make_pair(-666, -666);

   Long64_t nentries = fChain->GetEntriesFast();

   
   TH1F* h_amp  = new TH1F("h_amp_" + pos, "h_amp_" + pos, 1000, 0., .5);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if( !(ch1Amp > 0.02 && ch1Amp < 0.48 && ch3Amp > 0.01 && ch3Amp < 0.48) ) continue;
      h_amp->Fill(ch3Amp);
   }

      
   TF1* gaus1 = new TF1( "gaus1", "gaus(0)", h_amp->GetMean()-2.*h_amp->GetRMS() , h_amp->GetMean()+2.*h_amp->GetRMS() );
   h_amp->Fit( gaus1, "RQM");
   h_amp->GetXaxis()->SetRangeUser(h_amp->GetMean()-5.*h_amp->GetRMS() , h_amp->GetMean()+5.*h_amp->GetRMS());
   h_amp->Draw();
   c->SaveAs("h_amp_" + pos + ".pdf");
   
   return std::make_pair( gaus1->GetParameter(1), gaus1->GetParameter(2) );
}

std::pair<float,float> tree::GetTimeResolution(TString pos)
{

  TCanvas* c = new TCanvas( "c", "c", 800, 600 );
  c->SetHighLightColor(2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  
  if (fChain == 0) return std::make_pair(-666, -666);

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F* h_time = new TH1F("h_time_" + pos, "h_time_" + pos, 300, 25., 27.);
   //TH1F* h_time = new TH1F("h_time_" + pos, "h_time_" + pos, 300, 27., 29.);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if( !(ch1Amp > 0.02 && ch1Amp < 0.48 && ch3Amp > 0.01 && ch3Amp < 0.48) ) continue;
      h_time->Fill(ch3THM-ch1THM);
      //h_time->Fill(t3gausroot-ch1THM);
   }

   TF1* gaus1 = new TF1( "gaus1", "gaus(0)", h_time->GetMean()-2.*h_time->GetRMS() , h_time->GetMean()+2.*h_time->GetRMS() );
   h_time->Fit( gaus1, "RQM");
   h_time->GetXaxis()->SetRangeUser(h_time->GetMean()-5.*h_time->GetRMS() , h_time->GetMean()+5.*h_time->GetRMS());
   h_time->Draw();
   c->SaveAs("h_time_" + pos + ".pdf");

   return std::make_pair( gaus1->GetParameter(2), gaus1->GetParError(2) );
   
}


std::pair<float,float> tree::GetMeanTime(TString pos)
{

  TCanvas* c = new TCanvas( "c", "c", 800, 600 );
  c->SetHighLightColor(2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  
  if (fChain == 0) return std::make_pair(-666, -666);

   Long64_t nentries = fChain->GetEntriesFast();

   TH1F* h_time = new TH1F("h_time_" + pos, "h_time_" + pos, 300, 25., 27.);
   //TH1F* h_time = new TH1F("h_time_" + pos, "h_time_" + pos, 300, 27., 29.);
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if( !(ch1Amp > 0.02 && ch1Amp < 0.48 && ch3Amp > 0.01 && ch3Amp < 0.48) ) continue;
      h_time->Fill(ch3THM-ch1THM);
      //h_time->Fill(t3gausroot-ch1THM);
   }

   TF1* gaus1 = new TF1( "gaus1", "gaus(0)", h_time->GetMean()-2.*h_time->GetRMS() , h_time->GetMean()+2.*h_time->GetRMS() );
   h_time->Fit( gaus1, "RQM");
   h_time->GetXaxis()->SetRangeUser(h_time->GetMean()-5.*h_time->GetRMS() , h_time->GetMean()+5.*h_time->GetRMS());
   h_time->Draw();
   c->SaveAs("h_time_" + pos + ".pdf");

   return std::make_pair( gaus1->GetParameter(1), gaus1->GetParameter(2) );
   
}
