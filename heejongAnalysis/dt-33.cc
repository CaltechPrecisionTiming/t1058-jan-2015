#include <fstream>
#include <iomanip>
#include <string>
#include <iostream>


#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TPostScript.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <TLine.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TMapFile.h>
#include <TPaveStats.h>


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>  
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <time.h>


#include "TWaveForm2.h"

TFile* file;
TTree* tree;

TCanvas*      c1, *c2, *c3, *c4;
TPad*         pad1;
TStyle*       style;
TPostScript*  psf1;

int graphic_init();

int 
main(int argc, char **argv){

  double tnorm = 200.0;

  // read DRS4 time calibration data (dV)
  if(0){
    FILE* fp1;
    char stitle[200];
    sprintf( stitle, "/home/heejong/tmp/sine-01-%s.txt", argv[1]);

    fp1 = fopen( stitle, "r");
    printf("time calibration : %s\n", stitle);

    double amp_diff[1024], tcal[1024];  
    double amp_sum = 0;
    for( int j = 0; j < 1024; j++){  
      int dummy;
      fscanf( fp1, "%d  %lf\n", &dummy, &amp_diff[j]);       
      amp_sum += amp_diff[j];
      if(0) printf(" %8.4f\n", amp_diff[j]);
    }
  
    fclose(fp1);

    // normalize to (dT)
    for( int j = 0; j < 1024; j++){         
      tcal[j] = amp_diff[j]/amp_sum*tnorm;
      if(0) printf("tcal   %5d  =>   %7.3f\n", j, tcal[j]);
    }
  }


  // variables
  int event, tc1;
  float b1_t[1024], b1_c[4][1024];
  double m_b1_t[1024], m_b1_c[4][1024];

  double th_led[48];
  double th_cfd[48];
  double th_e[48], e_led[48], e_cfd[48];
  //  double e_led2[48], e_cfd2[48];


  // array for leading edge threshold, cfd fraction 
  for( int i = 0 ; i < 48; i++){
    th_led[i] = i*3.0+1.;
    th_cfd[i] = 0.02*i + 0.01;
    th_e[i] = 0;
  }


  // Histogram declaration
  graphic_init();
 
  c1 = new TCanvas("c1", "c1", 700, 500);
  c1->Divide(1, 1);

  char hid[100];  sprintf( hid, "h2");
  char htitle[100];  sprintf( htitle,"Event");
  
  TH2F* h2 = new TH2F( hid, htitle, 1024, 0, 1024, 200, -500, 4200);

  h2->SetBit(TH1::kNoStats);  
  h2->GetXaxis()->SetTitle("Time");
  h2->GetYaxis()->SetTitle("Amplitude");

  TH1F* h_amp[4], *h_amp2[4], *h_tc[4], *h_rms[4], *h_tc_time[1024];

  TH1F* h_tc_diff = new TH1F("h_tc_diff", "h_tc_diff", 100, 0, 100);

  for( int i = 0; i < 4; i++){
    sprintf( hid, "h_amp_%d", i);
    h_amp[i] = new TH1F( hid, hid, 300, 0, 600.); 

    sprintf( hid, "h_amp2_%d", i);
    h_amp2[i] = new TH1F( hid, hid, 300, 0, 600.); 

    sprintf( hid, "h_tc_%d", i);
    h_tc[i] = new TH1F( hid, hid, 1024, -0.5, 1023.5); 

    sprintf( hid, "h_rms_%d", i);
    h_rms[i] = new TH1F( hid, hid, 500, 0., 6.); 
  }


  for( int i = 0; i < 1024; i++){
    sprintf( hid, "h_tc_time_%d", i);
    sprintf( htitle, "TC  %4d", i);
    h_tc_time[i] = new TH1F( hid, htitle, 400, 0.19, 0.23); 
  }

  TH1F* h_ped[4];
  TH1F* h_ped2[4];
  for( int j = 0; j < 4; j++){
    sprintf( hid, "h_ped_c%d", j+1);
    h_ped[j] = new TH1F( hid, hid, 200, -20., 20.);

    sprintf( hid, "h_ped2_c%d", j+1);
    h_ped2[j] = new TH1F( hid, hid, 200, -20., 20.);
  }

  TH1F* h_led[48];
  for( int i = 0; i < 48; i++){
    sprintf( hid, "h_led_%d", i);
    sprintf( htitle, "LED %3.1f", i*3.0+1.);
    h_led[i] = new TH1F( hid, htitle, 3000, 1.5, 4.5);
  }

  TH1F* h_cfd[48];
  for( int i = 0; i < 48; i++){
    sprintf( hid, "h_cfd_%d", i);

    sprintf( htitle, "CFD %6.3f", 0.02*i+0.01);
    h_cfd[i] = new TH1F( hid, htitle, 1500, 0, 5);
  }


  TH1F* h_rise[4];
  for( int i = 0; i < 4; i++){
    sprintf( hid, "h_rise_%d", i);
    sprintf( htitle, "Rise Time %d", i);
    h_rise[i] = new TH1F( hid, htitle, 200, 0.5, 1.5);
  }


  //output ntuple file
  TFile *fout = new TFile("output.root","RECREATE");
  TTree* treeOut = new TTree("tree","tree");
  
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;
  float ch1Time = 0;
  float ch2Time = 0;
  float ch3Time = 0;
  float ch4Time = 0;
  treeOut->Branch("ch1Amp",&ch1Amp,"ch1Amp/F");
  treeOut->Branch("ch2Amp",&ch2Amp,"ch2Amp/F");
  treeOut->Branch("ch3Amp",&ch3Amp,"ch3Amp/F");
  treeOut->Branch("ch4Amp",&ch4Amp,"ch4Amp/F");
  treeOut->Branch("ch1Time",&ch1Time,"ch1Time/F");
  treeOut->Branch("ch2Time",&ch2Time,"ch2Time/F");
  treeOut->Branch("ch3Time",&ch3Time,"ch3Time/F");
  treeOut->Branch("ch4Time",&ch4Time,"ch4Time/F");


  // loop over root files
    
  char title[200];  
  //  sprintf( title, "/home/heejong/tb2015/data/%s-%s-%s.root", argv[1], argv[2], argv[3]);
  sprintf( title, "/uscms_data/d2/sxie/releases/CMSSW_7_2_3/src/Timing/t1058-jan-2015/%s.root", argv[1]);

  file = new TFile( title);
  printf("%s data processed.\n", title);

  tree = (TTree*)file->Get("pulse");

  TBranch* t_event = tree->GetBranch("event");
  TBranch* t_tc1   = tree->GetBranch("tc1");
  TBranch* t_b1_t  = tree->GetBranch("b1_t");
  TBranch* t_b1_c  = tree->GetBranch("b1_c");

  t_event->SetAddress( &event);
  t_tc1->SetAddress( &tc1);
  t_b1_t->SetAddress( b1_t);
  t_b1_c->SetAddress( b1_c);


  // loop over events
  for( int loop = 0; loop < tree->GetEntries(); loop++){

    tree->GetEntry(loop);

    if(loop%1000 == 0) printf(" event = %d \n", loop);

    for(int i=0; i<4; i++)
      for(int j=0; j<1024; j++)
	m_b1_c[i][j] = (double)b1_c[i][j];

    // use own calibrated time for DRS4
    //    m_b1_t[0] = 0;
    for(int j=0; j<1024; j++)
      m_b1_t[j] = (double)b1_t[j];
    //      m_b1_t[j] = (tcal[(j-1+tc1)%1024]+m_b1_t[j-1]);


    // 1. pedestal correction
    double ped_bc[4], rms[4];
      
    for( int i = 0; i<4; i++)
      ped_bc[i] = 0;
      
    // Pedestal subtraction
    for( int i = 0; i<4; i++)
      if (i<2) 
	for( int k = 400; k < 450; k++)
	  h_ped[i]->Fill( m_b1_c[i][k]);
      else if (i==2)
	for( int k = 800; k < 850; k++)
	  h_ped[i]->Fill( m_b1_c[i][k]);
      else 
	for( int k = 0; k < 50; k++)
	  h_ped[i]->Fill( m_b1_c[i][k]);
      
    for( int i = 0; i<4; i++){
      ped_bc[i]= h_ped[i]->GetMean();
      h_rms[i]->Fill( h_ped[i]->GetRMS());
      rms[i] =  h_ped[i]->GetRMS();
	
      // subtract pedestal 
      for( int k = 0; k < 1024; k++){
	m_b1_c[i][k] -= ped_bc[i];
	m_b1_c[i][k] *= (-1.); // invert the poloarity for negative pulse
      }	  
    }

    for( int i = 0; i < 4; i++)
      h_ped[i]->Reset();


    // 2. find positions of two peaks
    int    tc_max1 = -10;    
    int    tc_max2 = -10;    
    int    tc_max3 = -10;    
    int    tc_max4 = -10;    

    double max = -10.;
    for(int j=00; j<400; j++){
      if( m_b1_c[0][j] > max ){	    
	tc_max1 = j;
	max = m_b1_c[0][j];
      }
    }
         
    max  = -10.;
    for(int j=00; j<400; j++){
      if( m_b1_c[1][j] > max){	    
	tc_max2 = j;
	max = m_b1_c[1][j];
      }
    }  

    max  = -10.;
    for(int j=000; j<400; j++){
      if( m_b1_c[2][j] > max){	    
	tc_max3 = j;
	max = m_b1_c[2][j];
      }
    }  

    max  = -10.;
    for(int j=700; j<1024; j++){
      if( m_b1_c[3][j] > max){	    
	tc_max4 = j;
	max = m_b1_c[3][j];
      }
    }  

    
    // //For Photonis
    // if(  m_b1_c[0][tc_max1] < 50 || m_b1_c[2][tc_max3] < 100 ) continue;
    // if(  m_b1_c[0][tc_max1] > 490 || m_b1_c[2][tc_max3] > 490 ) continue;
    // if(  m_b1_c[3][tc_max4] < 100 ) continue; //cherenkov cut

    // //For photek -> photek
    if(  m_b1_c[0][tc_max1] < 10 || m_b1_c[1][tc_max2] < 10 ) continue;
    // if(  m_b1_c[0][tc_max1] > 490 || m_b1_c[1][tc_max2] > 490) continue;
    //if(  m_b1_c[3][tc_max4] < 100 ) continue; //cherenkov cut

    // //For photek -> photek protons
    // if(  m_b1_c[0][tc_max1] < 100 || m_b1_c[1][tc_max2] < 50 ) continue;
    // if(  m_b1_c[0][tc_max1] > 490 || m_b1_c[1][tc_max2] > 490) continue;


    // fill maximum amplitude peak
    h_amp[0]->Fill( m_b1_c[0][tc_max1]);
    h_amp[1]->Fill( m_b1_c[1][tc_max2]);
    h_amp[2]->Fill( m_b1_c[2][tc_max3]);
    h_amp[3]->Fill( m_b1_c[3][tc_max4]);
          

    // find the sample corresponding threshold (50mV)
    tc_max1 = -10;    
    for(int j=00; j<600; j++){
      if( m_b1_c[0][j] > 10 ){	    
	tc_max1 = j;
	break;
      }
    }
         
    tc_max2 = -10;    
    for(int j=00; j<600; j++){
      if( m_b1_c[1][j] > 10){	    
	tc_max2 = j;
	break;
      }
    }            

    tc_max3 = -10;    
    for(int j=00; j<600; j++){
      if( m_b1_c[2][j] > 10){	    
	tc_max3 = j;
	break;
      }
    }            
           
    h_tc[0]->Fill( tc_max1);
    h_tc[1]->Fill( tc_max2);
    h_tc[2]->Fill( tc1);


    h_tc_diff->Fill( tc_max2 -tc_max1);

    // if( tc_max1 < 50 || tc_max2 < 50 ) continue;
     
    // make waveforms
    char stitle[100];
    sprintf( stitle, "Event %5d", loop);
      
    TWaveForm* wf01 = new TWaveForm( 1024, m_b1_t, m_b1_c[0]);
    TWaveForm* wf02 = wf01->Resize( m_b1_t[tc_max1] - 7, m_b1_t[tc_max1] + 7);
    TWaveForm* wf03 = wf02->UpSample2(0.005); 
    TWaveForm* wf031 = wf03->LowPass2(100);
    TWaveForm* wf04 = wf031->Resize( m_b1_t[tc_max1] -  4, m_b1_t[tc_max1] + 4);

    TWaveForm* wf11 = new TWaveForm( 1024, m_b1_t, m_b1_c[1]);
    TWaveForm* wf12 = wf11->Resize( m_b1_t[tc_max2] - 7, m_b1_t[tc_max2] + 7);
    TWaveForm* wf13 = wf12->UpSample2(0.005);  
    TWaveForm* wf131 = wf13->LowPass2(100);
    TWaveForm* wf14 = wf131->Resize( m_b1_t[tc_max2] - 4, m_b1_t[tc_max2] + 4);

    // TWaveForm* wf21 = new TWaveForm( 1024, m_b1_t, m_b1_c[2]);
    // TWaveForm* wf22 = wf21->Resize( m_b1_t[tc_max3] - 7, m_b1_t[tc_max3] + 7);
    // TWaveForm* wf23 = wf22->UpSample2(0.005);  
    // TWaveForm* wf231 = wf23->LowPass2(100);
    // TWaveForm* wf24 = wf231->Resize( m_b1_t[tc_max3] - 4, m_b1_t[tc_max3] + 4);

    double *wf04_y = wf04->Get_samp_y();
    double *wf14_y = wf14->Get_samp_y();
    //double *wf24_y = wf24->Get_samp_y();

    h_rise[0]->Fill( wf03->Risetime19());
    h_rise[1]->Fill( wf13->Risetime19());
    //h_rise[2]->Fill( wf23->Risetime19());

    for( int i = 0; i < 400; i++){
      h_ped2[0]->Fill( wf04_y[i]);
      h_ped2[1]->Fill( wf14_y[i]);
       //h_ped2[2]->Fill( wf24_y[i]);
    }

    h_amp2[0]->Fill( wf04->Amplitude());
    h_amp2[1]->Fill( wf14->Amplitude());
    //h_amp2[2]->Fill( wf24->Amplitude());


    // get timing from the discriminators
    double led[4][48], cfd[4][48];
      
    double* t0_cfd = wf04->CFD( th_cfd, 48);	  
	  
    for( int k = 0; k < 48; k++){
      cfd[0][k] = t0_cfd[k];
    }

    //    double* t1_led = wf14->LED( th_led, 48);	  
    double* t1_cfd = wf14->CFD( th_cfd, 48);	  
	  
    for( int k = 0; k < 48; k++){
      cfd[1][k] = t1_cfd[k];
    }

    // // double* t2_led = wf24->LED( th_led, 48);	  
    // double* t2_cfd = wf24->CFD( th_cfd, 48);	  
	  
    // for( int k = 0; k < 48; k++){
    //   cfd[2][k] = t2_cfd[k];
    // }

      
    for( int i = 0; i < 48; i++){
      h_cfd[i]->Fill( cfd[1][i]-cfd[0][i]);
    }


    //Fill output tree
    ch1Amp = wf04->Amplitude();
    ch2Amp = wf14->Amplitude();
    ch3Amp = m_b1_c[2][tc_max3];
    //ch3Amp = wf24->Amplitude();
    ch4Amp = m_b1_c[3][tc_max4];
    ch1Time = cfd[0][37];
    ch2Time = cfd[1][37];
    ch3Time = cfd[2][37];
    ch4Time = cfd[3][37];
    treeOut->Fill();

    delete wf01;
    delete wf02;
    delete wf03;
    delete wf031;
    delete wf04;   

    delete wf11;
    delete wf12;
    delete wf13;
    delete wf131;
    delete wf14;   
  }

  file->Close();


  char pstitle[200];

  // sprintf( pstitle, "%s-%s-%s-%s.ps", argv[0], argv[1], argv[2], argv[3]);
  sprintf( pstitle, "%s-%s.ps", argv[0], argv[1]);
  psf1 = new TPostScript( pstitle, 112);

  c2 = new TCanvas("c2", "c2", 700, 500);
  c2->Divide(2, 2);

  //  gStyle->SetStatX(0.5);
    
  psf1->NewPage();
  for( int i = 0; i < 4; i++){
    c2->cd(i+1);

    h_amp[i]->Draw();

    if(0){
      //h_amp[i]->Fit("gaus", "QE", "", 200, 450);   
      
      h_amp[i]->GetXaxis()->SetTitle("Amplitude (mV)");
      
      double mean  = h_amp[i]->GetFunction("gaus")->GetParameter(1);
      double sigma = h_amp[i]->GetFunction("gaus")->GetParameter(2);
      
      printf("i = %d,   mean =  %5.1f  +- %5.1f,  resol = %5.1f\n",
	     i, mean, sigma, sigma/mean*2.35*100.);
    }
    
    c2->Update();
  }
  

  psf1->NewPage();
  for( int i = 0; i < 4; i++){
    c2->cd(i+1);

    h_amp2[i]->Draw();

    //    if( i < 2) h_amp2[i]->Fit("gaus", "QE", "", 300, 450);   
    //    if( i > 2) h_amp2[i]->Fit("gaus", "QE", "", 250, 490);    

    c2->Update();
  }


  gStyle->SetStatX(0.95);

  psf1->NewPage();
  for( int i = 0; i < 4; i++){
    c2->cd(i+1);
    h_tc[i]->Draw();
    c2->Update();
  }


  psf1->NewPage();
  for( int i = 0; i < 3; i++){
    c2->cd(i+1);
    h_rise[i]->Draw();
    c2->Update();
  }


  if(1){
    psf1->NewPage();

    c3 = new TCanvas("c3", "c3", 10, 10, 700, 500);
    c3->SetFillStyle( 4000);
    c3->Divide( 3, 2);


    int maxbin = 0;
    double tled[48], tcfd[48];
    double tled_mean[48], tcfd_mean[48];
    //  double tled2[48], tcfd2[48];
    //  double tled_mean2[48], tcfd_mean2[48];
    double mean;


    for( int i = 0; i < 8; i++){

      if(i) psf1->NewPage();

      for( int j = 0; j < 6; j++){
	cout << "here1 " << i << " " << j << "\n";

	c3->cd(j+1);
	maxbin = h_cfd[i*6+j]->GetMaximumBin();
	mean =  h_cfd[i*6+j]->GetBinCenter(maxbin);
	cout << "here2 " << i << " " << j << "\n";
	h_cfd[i*6+j]->SetAxisRange( 0.2, 0.4, "x");
	h_cfd[i*6+j]->Draw();

	cout << "here3 " << i << " " << j << "\n";

	if(1){
	  //h_cfd[i*6+j]->Fit("gaus", "QE", 0, mean-0.1, mean+0.1);     
	  //h_cfd[i*6+j]->Fit("gaus", "QE", 0, 0.35, 0.5);     
	  h_cfd[i*6+j]->Fit("gaus", "QE", 0, 0.2, 0.4);     
	  tcfd_mean[i*6+j] = h_cfd[i*6+j]->GetFunction("gaus")->GetParameter(1)*1000;
	  tcfd[i*6+j] = h_cfd[i*6+j]->GetFunction("gaus")->GetParameter(2)*1000*2.35;
	  e_cfd[i*6+j] = h_cfd[i*6+j]->GetFunction("gaus")->GetParError(2)*1000*2.35;
	  cout << "here4 " << i << " " << j << "\n";
	}
	c3->Update();
      }
  
    }
  }
  //  printf("%s :  dt_mean   =  %8.3f,   dt_fwhm =  %8.3f\n", argv[1], tcfd_mean[30], tcfd[30]);

  //treeOut->Write();
  fout->WriteTObject(treeOut);

  psf1->Close();
  fout->Close();

  return 0;
}



int
graphic_init( void){

  style = new TStyle("style", "style");
  
  style->SetLabelFont(132,"X");
  style->SetLabelFont(132,"Y");
  style->SetTitleFont(132,"X");
  style->SetTitleFont(132,"Y");
  style->SetTitleFont(132,"");
  style->SetTitleFontSize( 0.07);
  style->SetStatFont(132);
  style->GetAttDate()->SetTextFont(132);

  style->SetStatW(0.2);
  style->SetStatH(0.25);

  style->SetFuncColor(2);
  style->SetFuncWidth(2);
  style->SetLineWidth(2);
  
  style->SetOptFile(0);
  style->SetOptTitle(1);

  style->SetFrameBorderMode(0);
  style->SetCanvasBorderMode(0);
  style->SetPadBorderMode(0);
  style->SetTitleStyle(4000);
  style->SetPadColor(0);
  style->SetCanvasColor(0);

  style->SetTitleFillColor(0);

  style->SetTitleBorderSize(0);
  //  style->SetTitleX(0.3);
  //  style->SetTitleY(0.06);
  style->SetStatColor(0);
  style->SetStatBorderSize(1);

  style->SetOptStat("emri");
  // style->SetOptStat(1);
  style->SetOptFit(1);
  style->SetTitleOffset( 1.0,"Y");

  style->SetMarkerStyle(20);
  style->SetMarkerSize( 0.3);
  style->SetMarkerColor(4);

  // style->SetOptDate(21);

  //  style->SetPadGridX(1);
  //  style->SetPadGridY(1);


  style->cd();

  return 0;
}
