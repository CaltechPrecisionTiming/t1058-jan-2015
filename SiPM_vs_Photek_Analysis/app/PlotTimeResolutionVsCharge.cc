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
  //int n = 14;
/*
  //  ND:         0p5,  0p8,   1p1,   1p4,   1p6,  1p8_01, 1p8_02, 1p8_03, 04, 05, 06 , 07, 08, 09, 10, 11, 12  
  //  x: charge y: deltaT  
  float y[]    = {31.92, 36.35, 42.32, 60.91, 87.24, 51.80, 63.08, 62.32, 55.92, 54.85, 56.14, 53.31, 54.37, 54.04, 57.48, 53.79  };
  float erry[] = { 0.11,  0.12,  0.15,  0.06,  0.09,  1.22,  0.46,  0.51,  0.70,  0.59,  0.64,  0.94,  1.04,  0.83,  1.60,  3.25  };
  float yTHM[]    = {18.87, 21.72, 30.89, 43.04, 52.34, 64.25, 61.94, 58.13, 48.13, 46.47, 44.92, 42.59, 41.22, 41.60, 41.01, 39.72, 39.17};
  float erryTHM[] = { 0.06,  0.07,  0.10,  0.05,  0.09,  0.92,  0.37,  0.35,  0.37,  0.39,  0.35, 0.52,  0.44,  0.49, 0.62,  0.68,  0.83};
  float x[] =    {1.05, 0.89, 0.39, 0.16, 0.09, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12};
  float errx[] = {0.04, 0.06, 0.09, 0.06, 0.04, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};

  //  ND:         1p8_01, 1p8_02, 1p8_03, 04, 05, 06 , 07, 08, 09,  10, 11, 12,  0p8, 0p5  
  //  x: charge y: deltaT  
  float yTHM[]    = {64.25, 61.94, 58.13, 48.13, 46.47, 44.92, 42.59, 41.22, 41.60, 41.01, 39.72, 39.17, 21.72, 18.87};
  float erryTHM[] = { 0.92,  0.37,  0.35,  0.37,  0.39,  0.35,  0.52,  0.44,  0.49,  0.62,  0.68,  0.83,  0.07, 0.06};
  float x[] =    {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.89, 1.05};
  float errx[] = {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.06, 0.04};

  int n = 6;
  //  ND:         0p5,  0p8,   1p1,   1p4,   1p6,  1p8
  //  x: charge y: deltaT  
  float y[]    = {31.92, 36.35, 42.32, 55.67, 71.80, 80.88};
  float erry[] = { 0.11,  0.12,  0.15,  0.11,  0.21,  0.12};
  float yTHM[]    = {18.87, 21.72, 30.89, 43.04, 52.34, 60.21};
  float erryTHM[] = { 0.06,  0.07,  0.10,  0.05,  0.09,  0.09};
  float x[]    = {1.05, 0.89, 0.39, 0.16, 0.09, 0.06};
  float errx[] = {0.04, 0.06, 0.09, 0.06, 0.04, 0.03};

  int n = 9;
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
  
  int n = 9;
  float y[]    = {7.11, 8.15, 9.41, 11.81, 15.87, 19.49, 48.06, 64.65, 85.98};
  float erry[] = {0.05, 0.06, 0.07, 0.08, 0.11, 0.14, 0.37, 0.52, 1.47};
  float x[]    = {4.384, 2.742, 1.808, 0.9831, 0.5937, 0.3862, 0.1077, 0.06313, 0.04589};
  float errx[] = {0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001};

//Jiajing's measurements, no amp, Marcel's algorithm
  int n2 = 7;
  float y2[]    = {11.99, 14.21, 21.17, 29.68, 43.99, 63.85, 82.15};
  float erry2[] = {0.05, 0.05, 0.09, 0.1, 0.15, 0.28, 0.62};
  float x2[]    = {6.484, 3.653, 1.606, 0.854, 0.414, 0.1931, 0.0973};
  float errx2[] = {0.001, 0.001, 0.001, 0.001, 0.0004, 0.0005, 0.0004};

  /* int n2 = 6;
  float y2[]    = {22.94, 30.05, 31.43, 44.37, 65.92, 95.0}; //22.36
  float erry2[] = {0.05, 0.08, 0.11, 0.15, 0.25, 0.5}; //0.05
  float x2[]    = {12.654, 6.956, 6.233, 2.991, 1.081, 0.1946}; //38.12
  float errx2[] = {0.005, 0.004, 0.006, 0.004, 0.004, 0.003}; //0.01
  for (int i = 0; i < n2; i++){x2[i] = x2[i] / 63.095;}
  */
//BTL Laser Study
  /*
  //373nm, constant 70V
  int n = 7;
  float y[] = {31.58, 31.43, 49.27, 52.71, 60.48, 98.29, 317.06}; //Time resolution
  //float x[] = {0, 1.4, 40, 50, 70, 85, 90.5}; //Tune
  float x[] = {0.1032, 0.07676, 0.05022, 0.04476, 0.03658, 0.02042, 0.01056}; //Amplitude

  //373nm, constant 0 tune
  int n2 = 8;
  float y2[] = {215.9, 72.98, 53.49, 43.39, 31.58, 34.01, 32.94, 33.28}; //Time Resolution
  // float x2[] = {66, 67, 68, 69, 70, 71, 72, 73}; //Voltage
  float x2[] = {0.007656, 0.01656, 0.0284, 0.0423, 0.1032, 0.0793, 0.1016, 0.1232}; //Amplitude
  
  //405nm run 1, constant 0 tune
  int n = 7;
  float y[] = {33.36, 33.25, 33.18, 31.97, 29.44, 31.91, 31.06}; //Time resolution
  //float x[] = {66, 67, 68, 69, 70, 71, 72, 73}; //Voltage
  float x[] = {0.06743, 0.1338, 0.2048, 0.278, 0.904, 1.091, 1.258}; //Amplitude

  //405nm run 1, constant 70V
  int n2 = 6;
  float y2[] = {29.44, 32.03, 33.61, 33.57, 46.22, 59.25}; //Time resolution
  //float x2[] = {0, 40, 50, 70, 73.5, 74}; //Tune
  float x2[] = {0.904, 0.2512, 0.2256, 0.09234, 0.01473, 0.009577}; //Amplitude

  */
  //405nm run 2, constant 0 tune
  int n = 7;
  float y[] = {34.71, 34.62, 38.64, 33.98, 35.1, 34.18, 33.4}; //Time resolution
  //float x[] = {66, 67, 68, 69, 70, 71, 72}; //Voltage
  float x[] = {0.0773, 0.1715, 0.2838, 0.4088, 0.5214, 0.6678, 0.8151}; //Amplitude

  //405nm run 2, constant 70V
  int n2 = 6;
  float y2[] = {35.1, 32.2, 31.23, 40.08, 44.37, 54.7}; //Time resolution
  //float x2[] = {0, 40, 50, 70, 73.5, 74}; //Tune
  float x2[] = {0.5214, 0.305, 0.2613, 0.09602, 0.02159, 0.01331};
  

  /*
  for(int i = 0; i < n; i++)
  {
    x[i] = x[i] / (.0973/63.095);
    errx[i] = errx[i] / (.0973/63.095);
  }

  for(int i = 0; i < n2; i++)
  {
    x2[i] = x2[i] / (0.0118);
    errx2[i] = errx2[i] / (0.0118);
  }
  */

  // TGraphErrors* g = new TGraphErrors(n,x,y,errx,erry);
  // TGraphErrors* g2 = new TGraphErrors(n2,x2,y2,errx2,erry2);
   TGraphErrors* g = new TGraphErrors(n,x,y);
   TGraphErrors* g2 = new TGraphErrors(n2,x2,y2);

   TCanvas* c = new TCanvas( "c", "c", 2119, 33, 800, 700 );
   // c->SetLogx();
   
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
   // TF1* f1Res = new TF1("f1Res","[0]/x + [1]/sqrt(x) + [2]",x[0],x[n-1]);
   // TF1* f1Res2 = new TF1("f1Res2","[0]/x + [1]/sqrt(x) + [2]",x2[0],x2[n2-1]);
   //f1Res2->SetParLimits(0,0.0,1e12);
   //f1Res2->SetParLimits(1,0.0,1e12);
   //f1Res2->SetParLimits(2,0.0,1e12);
   //f1Res->SetLineColor(1);

   // g->SetTitle(";Number of Photons;Time Resolution [ps]");
   //g->GetYaxis()->SetTitle("time resolution [ps]");
   //g->GetXaxis()->SetTitle("charge [pC]");
   g->GetYaxis()->SetRange(30, 55);

   g2->Draw("AP");
   g2->SetMarkerStyle(20);
   g2->SetMarkerColor(kBlack);
   g->Draw("P");
   g->SetMarkerStyle(20);
   g->SetMarkerColor(kRed);   

   // g->Fit(f1Res,"E");
   //g2->Fit(f1Res2, "E");
   
   TLatex latex;
   latex.SetNDC();
   latex.SetTextAngle(0);
   latex.SetTextColor(kBlack);    
   latex.SetTextFont(42);
   latex.SetTextAlign(31); 
   latex.SetTextSize(0.02);    
   //latex.DrawLatex(0.90, 0.85, Form("3x3mm SiPM: #sigma_{t} = #frac{%.2f #pm %.2f [ps]}{Q (pC)} + #frac{%.2f #pm %.2f [ps]}{#sqrt{Q (pC)}} + %.2f #pm %.2f [ps]", f1Res->GetParameter(0), f1Res->GetParError(0), f1Res->GetParameter(1), f1Res->GetParError(1), f1Res->GetParameter(2), f1Res->GetParError(2) ));
   //latex.SetTextColor(kRed);
    //latex.DrawLatex(0.90, 0.7, Form("1x1mm SiPM: #sigma_{t} = #frac{%.2f #pm %.2f [ps]}{Q (pC)} + #frac{%.2f #pm %.2f [ps]}{#sqrt{Q (pC)}} + %.2f #pm %.2f [ps]", f1Res2->GetParameter(0), f1Res2->GetParError(0), f1Res2->GetParameter(1), f1Res2->GetParError(1), f1Res2->GetParameter(2), f1Res2->GetParError(2) ));
    // latex.DrawLatex(0.90, 0.85, Form("3x3mm SiPM: #sigma_{t} = #frac{%.2f #pm %.2f [ps]}{N} + #frac{%.2f #pm %.2f [ps]}{#sqrt{N}} + %.2f #pm %.2f [ps]", f1Res->GetParameter(0), f1Res->GetParError(0), f1Res->GetParameter(1), f1Res->GetParError(1), f1Res->GetParameter(2), f1Res->GetParError(2) ));
   // latex.SetTextColor(kRed);
   // latex.DrawLatex(0.90, 0.7, Form("1x1mm SiPM: #sigma_{t} = #frac{%.2f #pm %.2f [ps]}{N} + #frac{%.2f #pm %.2f [ps]}{#sqrt{N}} + %.2f #pm %.2f [ps]", f1Res2->GetParameter(0), f1Res2->GetParError(0), f1Res2->GetParameter(1), f1Res2->GetParError(1), f1Res2->GetParameter(2), f1Res2->GetParError(2) ));
    c->SaveAs("405nm_amplitudes2.pdf");
   // c->SaveAs("deltaT_vs_NPhotons_SiPM_3vs1logx.png");
   // c->SaveAs("deltaT_vs_NPhotons_SiPM_3vs1logx.C");
    // c->SaveAs("deltaT_vs_Charge_SiPM_3vs1.pdf");
    //c->SaveAs("deltaT_vs_Charge_SiPM_3vs1.png");
   //c->SaveAs("deltaT_vs_Charge_SiPM_3vs1.C");
   
   /*
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
   */
  return 0;
}

