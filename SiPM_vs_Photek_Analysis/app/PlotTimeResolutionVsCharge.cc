#include <iostream>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TMath.h>
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
  int n = 14;
/*
  //  ND:         0p5,  0p8,   1p1,   1p4,   1p6,  1p8_01, 1p8_02, 1p8_03, 04, 05, 06 , 07, 08, 09, 10, 11, 12  
  //  x: charge y: deltaT  
  float y[]    = {31.92, 36.35, 42.32, 60.91, 87.24, 51.80, 63.08, 62.32, 55.92, 54.85, 56.14, 53.31, 54.37, 54.04, 57.48, 53.79  };
  float erry[] = { 0.11,  0.12,  0.15,  0.06,  0.09,  1.22,  0.46,  0.51,  0.70,  0.59,  0.64,  0.94,  1.04,  0.83,  1.60,  3.25  };
  float yTHM[]    = {18.87, 21.72, 30.89, 43.04, 52.34, 64.25, 61.94, 58.13, 48.13, 46.47, 44.92, 42.59, 41.22, 41.60, 41.01, 39.72, 39.17};
  float erryTHM[] = { 0.06,  0.07,  0.10,  0.05,  0.09,  0.92,  0.37,  0.35,  0.37,  0.39,  0.35, 0.52,  0.44,  0.49, 0.62,  0.68,  0.83};
  float x[] =    {1.05, 0.89, 0.39, 0.16, 0.09, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12};
  float errx[] = {0.04, 0.06, 0.09, 0.06, 0.04, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
*/
  //  ND:         1p8_01, 1p8_02, 1p8_03, 04, 05, 06 , 07, 08, 09,  10, 11, 12,  0p8, 0p5  
  //  x: charge y: deltaT  
  float yTHM[]    = {64.25, 61.94, 58.13, 48.13, 46.47, 44.92, 42.59, 41.22, 41.60, 41.01, 39.72, 39.17, 21.72, 18.87};
  float erryTHM[] = { 0.92,  0.37,  0.35,  0.37,  0.39,  0.35,  0.52,  0.44,  0.49,  0.62,  0.68,  0.83,  0.07, 0.06};
  float x[] =    {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.89, 1.05};
  float errx[] = {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.06, 0.04};
/*
  int n = 6;
  //  ND:         0p5,  0p8,   1p1,   1p4,   1p6,  1p8
  //  x: charge y: deltaT  
  float y[]    = {31.92, 36.35, 42.32, 55.67, 71.80, 80.88};
  float erry[] = { 0.11,  0.12,  0.15,  0.11,  0.21,  0.12};
  float yTHM[]    = {18.87, 21.72, 30.89, 43.04, 52.34, 60.21};
  float erryTHM[] = { 0.06,  0.07,  0.10,  0.05,  0.09,  0.09};
  float x[]    = {1.05, 0.89, 0.39, 0.16, 0.09, 0.06};
  float errx[] = {0.04, 0.06, 0.09, 0.06, 0.04, 0.03};
  
  int n = 5;
  //  ND:         0p5,  0p8,   1p1,   1p4,   1p6
  //  x: charge y: deltaT  
  float y[]    = {31.92, 36.35, 42.32, 60.91, 87.24};
  float erry[] = { 0.11,  0.12,  0.15,  0.06,  0.09};
  float yTHM[]    = {18.87, 21.72, 30.89, 44.31, 63.24};
  float erryTHM[] = { 0.06,  0.07,  0.10,  0.05,  0.06};
  float x[]    = {1.05, 0.88, 0.40, 0.17, 0.10};
  float errx[] = {0.04, 0.07, 0.10, 0.06, 0.05};

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
   TF1* f1Res = new TF1("f1Res","[0]/x + [1]/sqrt(x) + [2]",x[0],1.1);

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
*/ 
   TGraphErrors* gr = new TGraphErrors(n,x,yTHM,errx,erryTHM);

   TCanvas* cr = new TCanvas( "cr", "cr", 2119, 33, 800, 700 );
   cr->SetHighLightColor(2);
   cr->SetFillColor(0);
   cr->SetBorderMode(0);
   cr->SetBorderSize(2);
   cr->SetLeftMargin( leftMargin );
   cr->SetRightMargin( 1.6*rightMargin );
   cr->SetTopMargin( topMargin );
   cr->SetBottomMargin( bottomMargin );
   cr->SetFrameBorderMode(0);
   cr->SetFrameBorderMode(0);
   //gr->GetXaxis()->SetRangeUser(0,1.2);
   TF1* f1ResTHM = new TF1("f1ResTHM","[0]/x + [1]/sqrt(x) + [2]",x[0],x[n-1]);

   gr->SetTitle(";charge [pC];time resolution [ps]");
   //gr->GetYaxis()->SetTitle("time resolution [ps]");
   //gr->GetXaxis()->SetTitle("charge [pC]");
   
   gr->Draw("AP");
   gr->SetMarkerStyle(20);
   gr->SetMarkerColor(kBlack);
   
   gr->Fit(f1ResTHM,"RE");


   TLatex latexr;
   latexr.SetNDC();
   latexr.SetTextAngle(0);
   latexr.SetTextColor(kBlack);    
   latexr.SetTextFont(42);
   latexr.SetTextAlign(31); 
   latexr.SetTextSize(0.03);    
   latexr.DrawLatex(0.87, 0.85, Form("#sigma_{t} = #frac{%.2f #pm %.2f [ps]}{Q (pC)} + #frac{%.2f #pm %.2f [ps]}{#sqrt{Q (pC)}} + %.2f #pm %.2f [ps]", f1ResTHM->GetParameter(0), f1ResTHM->GetParError(0), f1ResTHM->GetParameter(1), f1ResTHM->GetParError(1), f1ResTHM->GetParameter(2), f1ResTHM->GetParError(2) ));
   cr->SaveAs("deltaTTHM_vs_Charge_SiPM_1x1mm.pdf");
   cr->SaveAs("deltaTTHM_vs_Charge_SiPM_1x1mm.png");
   cr->SaveAs("deltaTTHM_vs_Charge_SiPM_1x1mm.C");

  return 0;
}
