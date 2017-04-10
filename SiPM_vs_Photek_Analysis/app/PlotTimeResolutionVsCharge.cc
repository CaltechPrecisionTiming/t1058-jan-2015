#include <iostream>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>

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

int main ( int argc, char** argv )
{
  int n = 7;
  float y[]    = {9.77, 12.84, 19.82, 27.98, 41.24,63.64, 84.84};
  float erry[] = {0.06, 0.06, 0.09, 0.14, 0.20,0.31, 0.48};
  float x[]    = {5.50, 3.06, 1.32, 0.69, 0.32, 0.15, 0.07};
  float errx[] = {0.27, 0.23, 0.16, 0.12, 0.08, 0.06, 0.03};

  TGraphErrors* g = new TGraphErrors(n,x,y,errx,erry);

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
   TF1* f1Res = new TF1("f1Res","[0]/x + [1]/sqrt(x) + [2]",x[0],x[n-1]);

   g->SetTitle(";charge [pC];time resolution [ps]");
   //g->GetYaxis()->SetTitle("time resolution [ps]");
   //g->GetXaxis()->SetTitle("charge [pC]");
   
   g->Draw("AP");
   g->SetMarkerStyle(20);
   g->SetMarkerColor(kBlack);
   
   g->Fit(f1Res,"RE");
   
   TLatex latex;
   latex.SetNDC();
   latex.SetTextAngle(0);
   latex.SetTextColor(kBlack);    
   latex.SetTextFont(42);
   latex.SetTextAlign(31); 
   latex.SetTextSize(0.03);    
   latex.DrawLatex(0.87, 0.85, Form("#sigma_{t} = #frac{%.2f #pm %.2f [ps]}{Q (pC)} + #frac{%.2f #pm %.2f [ps]}{#sqrt{Q (pC)}} + %.2f #pm %.2f [ps]", f1Res->GetParameter(0), f1Res->GetParError(0), f1Res->GetParameter(1), f1Res->GetParError(1), f1Res->GetParameter(2), f1Res->GetParError(2) ));
   c->SaveAs("deltaT_vs_Charge_SiPM_1x1mm.pdf");
   c->SaveAs("deltaT_vs_Charge_SiPM_1x1mm.png");
   c->SaveAs("deltaT_vs_Charge_SiPM_1x1mm.C");
  return 0;
}
