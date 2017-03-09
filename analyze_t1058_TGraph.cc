#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Math/Interpolator.h"
#include <assert.h>
#include "TGraph.h"
#include "TGraphErrors.h"

#define binsize 1

enum PulseQuality {
  kNegativePolarity       = 0x000001, // bit 0
  kSuddenJump             = 0x000002, // bit 1
  kFlatTop                = 0x000004, // bit 2
  kSecondPulse            = 0x000008, // bit 3
  kNoPulse                = 0x000010, // bit 4
  kLargeNegativeAmplitude = 0x000020, // bit 5
  kSaturated              = 0x000040  // bit 6
};

bool doFilter = false;

int FindMaxAbsolute( int n, float *a, bool _findMin = false );
TGraphErrors* GetTGraphFilter( float* channel, float* time, TString pulseName, bool makePlot );
float GetPulseIntegral(int peak, float *a, std::string option);
float GetBaseline(TGraphErrors * pulse, int i_low, int i_high, TString fname );
TGraphErrors GetTGraph(  float* channel, float* time, bool invert = false );
int FindRealMin( int n, float *a);
void RisingEdgeFitTime(TGraphErrors* pulse, const float index_min, float* tstamp, int event, TString fname, bool makePlot );

float LED( TH1F * pulse, double threshold, int nsamples, int splineBinFactor );
float LinearFit_Baseline(TH1F * pulse, const int index_min, const int range);
float LinearFit_Intercept(TH1F * pulse, const float base, const int index_first, const int index_last);
float GausFit_MeanTime(TGraphErrors* pulse, const int index_min, const int index_first, const int index_last);
void FitRisingEdge(TH1F* pulse, int nbinsL, int nbinsH, float &THM, float &risetime, float base);
void FitFullPulse(TH1F* pulse, float &par0, float &par1, float &par2);
TH1F* InterpolateWaveform(int nsamples, float* outputwaveform, float *inputwaveform, int splineBinFactor, std::string name);

const int Nsamples = 1024;

int main (int argc, char **argv)
{
  TFile *f;

  if (argc >= 3)
    {
      f = new TFile(argv[1]);

      std::cout << ">> Opening file " << argv[1] << " ......" << std::endl;
      if (!f->IsOpen())
	{			// terminate if the file can't be opened
	  std::cerr << "!! File open error:" << argv[1] << std::endl;
	  return 1;
	}
    }
  else
    {				// terminate if there is no input file or more than 1 input file
      std::cerr << "!! No input file" << std::endl;
      return 1;
    }


  bool includePulseshapeInOutput = false;
  if (argc >= 4) includePulseshapeInOutput = bool(atoi(argv[3]));
  
  //const int splineBinFactor = 40;
  const int splineBinFactor = 1;

  int t_[Nsamples]; for (int i = 0; i < 1024; i++) { t_[i] = i; }
  float Channel1VoltagesRaw_[Nsamples];
  float Channel2VoltagesRaw_[Nsamples];
  float Channel3VoltagesRaw_[Nsamples];
  float Channel4VoltagesRaw_[Nsamples];

  float Channel1Voltages_[Nsamples*splineBinFactor];
  float Channel2Voltages_[Nsamples*splineBinFactor];
  float Channel3Voltages_[Nsamples*splineBinFactor];
  float Channel4Voltages_[Nsamples*splineBinFactor];

  float ti1_[Nsamples], ti2_[Nsamples], ti3_[Nsamples], ti4_[Nsamples];

  bool convert2Volts = true;
  float amp[4];
  float integral[4];
  float gauspeak[4];
  float linearTime0[4];

  TTree* t1 = (TTree*)f->Get("p");      // andriy's converter
  if (t1) 
    {
      t1->SetBranchAddress("c1",Channel1VoltagesRaw_);
      t1->SetBranchAddress("c2",Channel2VoltagesRaw_);
      t1->SetBranchAddress("c3",Channel3VoltagesRaw_);
      t1->SetBranchAddress("c4",Channel4VoltagesRaw_);

      t1->SetBranchAddress("t1",ti1_);
      t1->SetBranchAddress("t2",ti2_);
      t1->SetBranchAddress("t3",ti3_);
      t1->SetBranchAddress("t4",ti4_);      
    }
  
  if (!t1) 
    {
      t1 = (TTree*)f->Get("T");   // artur's converter
      convert2Volts = false;

      t1->SetBranchAddress("c1",Channel1VoltagesRaw_);
      t1->SetBranchAddress("c2",Channel2VoltagesRaw_);
      t1->SetBranchAddress("c3",Channel3VoltagesRaw_);
      t1->SetBranchAddress("c4",Channel4VoltagesRaw_);

      t1->SetBranchAddress("t1",ti1_);
      t1->SetBranchAddress("t2",ti2_);
      t1->SetBranchAddress("t3",ti3_);
      t1->SetBranchAddress("t4",ti4_);
    }

  // Create the output file with a TTree
  TFile* fout;
  if(strncmp(argv[2], "same", 100) == 0){
    std::string fn(argv[1]);
    int pf = fn.find(".root");
    int pi = fn.rfind("/")+1;
    fn = "AnaFiles/" + fn.substr(pi, pf-pi) + "_ana.root";
    std::cout << "fname: " << fn << std::endl;
    //return 0;
    fout = new TFile(fn.c_str(),"recreate");
  }else{
    fout = new TFile(argv[2],"recreate");
  }
  TTree* treeOut = new TTree("tree","tree");
  
  unsigned int eventNumber = 0;
  float ch1Time_gausfitroot = 0;
  float ch2Time_gausfitroot = 0;
  float ch3Time_gausfitroot = 0;
  float ch4Time_gausfitroot = 0;
  
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;

  float ch1THM = 0;//Time at half the Maximum
  float ch2THM = 0;//Time at half the Maximum
  float ch3THM = 0;//Time at half the Maximum
  float ch4THM = 0;//Time at half the Maximum

  float ch1Risetime = 0;
  float ch2Risetime = 0;
  float ch3Risetime = 0;
  float ch4Risetime = 0;

  float ch1_TFF = 0.0;
  float ch2_TFF = 0.0;
  float ch3_TFF = 0.0;
  float ch4_TFF = 0.0;

  float ch1_TFF_v2 = 0.0;
  float ch2_TFF_v2 = 0.0;
  float ch3_TFF_v2 = 0.0;
  float ch4_TFF_v2 = 0.0;
  
  float ch1BL = 0.0;
  float ch2BL = 0.0;
  float ch3BL = 0.0;
  float ch4BL = 0.0;

  float ch1_AFF = 0.0;
  float ch2_AFF = 0.0;
  float ch3_AFF = 0.0;
  float ch4_AFF = 0.0;

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
  
  float ch1chisq = -1;
  float ch2chisq = -1;
  float ch3chisq = -1;
  float ch4chisq = -1;
  


  treeOut->Branch("event",&eventNumber,"event/i");
  treeOut->Branch("t1gausroot",&ch1Time_gausfitroot,"t1gausroot/F");
  treeOut->Branch("t2gausroot",&ch2Time_gausfitroot,"t2gausroot/F");
  treeOut->Branch("t3gausroot",&ch3Time_gausfitroot,"t3gausroot/F");
  treeOut->Branch("t4gausroot",&ch4Time_gausfitroot,"t4gausroot/F");
  
  treeOut->Branch("ch1Amp",&ch1Amp,"ch1Amp/F");
  treeOut->Branch("ch2Amp",&ch2Amp,"ch2Amp/F");
  treeOut->Branch("ch3Amp",&ch3Amp,"ch3Amp/F");
  treeOut->Branch("ch4Amp",&ch4Amp,"ch4Amp/F");

  treeOut->Branch("ch1THM",&ch1THM,"ch1THM/F");
  treeOut->Branch("ch2THM",&ch2THM,"ch2THM/F");
  treeOut->Branch("ch3THM",&ch3THM,"ch3THM/F");
  treeOut->Branch("ch4THM",&ch4THM,"ch4THM/F");

  treeOut->Branch("ch1Risetime",&ch1Risetime,"ch1Risetime/F");
  treeOut->Branch("ch2Risetime",&ch2Risetime,"ch2Risetime/F");
  treeOut->Branch("ch3Risetime",&ch3Risetime,"ch3Risetime/F");
  treeOut->Branch("ch4Risetime",&ch4Risetime,"ch4Risetime/F");

  treeOut->Branch("ch1BL",&ch1BL,"ch1BL/F");
  treeOut->Branch("ch2BL",&ch2BL,"ch2BL/F");
  treeOut->Branch("ch3BL",&ch3BL,"ch3BL/F");
  treeOut->Branch("ch4BL",&ch4BL,"ch4BL/F");
  
  treeOut->Branch("ch1_TFF", &ch1_TFF, "ch1_TFF/F");
  treeOut->Branch("ch2_TFF", &ch2_TFF, "ch2_TFF/F");
  treeOut->Branch("ch3_TFF", &ch3_TFF, "ch3_TFF/F");
  treeOut->Branch("ch4_TFF", &ch4_TFF, "ch4_TFF/F");
  
  treeOut->Branch("ch1_TFF_v2", &ch1_TFF_v2, "ch1_TFF_v2/F");
  treeOut->Branch("ch2_TFF_v2", &ch2_TFF_v2, "ch2_TFF_v2/F");
  treeOut->Branch("ch3_TFF_v2", &ch3_TFF_v2, "ch3_TFF_v2/F");
  treeOut->Branch("ch4_TFF_v2", &ch4_TFF_v2, "ch4_TFF_v2/F");

  
  treeOut->Branch("ch1_AFF", &ch1_AFF, "ch1_AFF/F");
  treeOut->Branch("ch2_AFF", &ch2_AFF, "ch2_AFF/F");
  treeOut->Branch("ch3_AFF", &ch3_AFF, "ch3_AFF/F");
  treeOut->Branch("ch4_AFF", &ch4_AFF, "ch4_AFF/F");
  
  treeOut->Branch("ch1QualityBit",&ch1QualityBit,"ch1QualityBit/i");
  treeOut->Branch("ch2QualityBit",&ch2QualityBit,"ch2QualityBit/i");
  treeOut->Branch("ch3QualityBit",&ch3QualityBit,"ch3QualityBit/i");
  treeOut->Branch("ch4QualityBit",&ch4QualityBit,"ch4QualityBit/i");
  
  treeOut->Branch("ch1Int",&ch1Int,"ch1Int/F");
  treeOut->Branch("ch2Int",&ch2Int,"ch2Int/F");
  treeOut->Branch("ch3Int",&ch3Int,"ch3Int/F");
  treeOut->Branch("ch4Int",&ch4Int,"ch4Int/F");
  
  treeOut->Branch("ch1chisq",&ch1chisq,"ch1chisq/F");
  treeOut->Branch("ch2chisq",&ch2chisq,"ch2chisq/F");
  treeOut->Branch("ch3chisq",&ch3chisq,"ch3chisq/F");
  treeOut->Branch("ch4chisq",&ch4chisq,"ch4chisq/F");  

  if (includePulseshapeInOutput) {
    
    treeOut->Branch("c1",Channel1VoltagesRaw_,"c1[1024]/F");
    treeOut->Branch("c2",Channel2VoltagesRaw_,"c2[1024]/F");
    treeOut->Branch("c3",Channel3VoltagesRaw_,"c3[1024]/F");
    treeOut->Branch("c4",Channel4VoltagesRaw_,"c4[1024]/F");
    treeOut->Branch("t",t_,"t[1024]/I");
    
    treeOut->Branch("ti1",ti1_,"ti1[1024]/D");
    treeOut->Branch("ti2",ti2_,"ti2[1024]/D");
    treeOut->Branch("ti3",ti3_,"ti3[1024]/D");
    treeOut->Branch("ti4",ti4_,"ti4[1024]/D");
  }

  //read all entries and fill the histograms
  Long64_t nentries = t1->GetEntries();
  //nentries = 100;

  for (Long64_t iEntry=0;iEntry<nentries;iEntry++)   
    {
      if(iEntry%1000==0) std::cout<<"Processing Event: "<<iEntry<<" out of: "<<nentries<<std::endl;
      
      t1->GetEntry(iEntry);
      eventNumber = iEntry+1;

      //Make Pulse shape Graph
      TString pulseName1 = Form("pulse_event%d_ch1", iEntry);
      TString pulseName2 = Form("pulse_event%d_ch2", iEntry);
      TString pulseName3 = Form("pulse_event%d_ch3", iEntry);
      TString pulseName4 = Form("pulse_event%d_ch4", iEntry);
      
      TGraphErrors* pulse1 = new TGraphErrors( GetTGraph( Channel1VoltagesRaw_, ti1_) );
      TGraphErrors* pulse2 = new TGraphErrors( GetTGraph( Channel2VoltagesRaw_, ti2_) );
      TGraphErrors* pulse3 = new TGraphErrors( GetTGraph( Channel3VoltagesRaw_, ti3_) );
      TGraphErrors* pulse4 = new TGraphErrors( GetTGraph( Channel4VoltagesRaw_, ti4_) );

      //estimate baseline
      float baseline1;
      float baseline2;
      float baseline3;
      float baseline4;
      
      baseline1 = GetBaseline( pulse1, 5 ,50, pulseName1);
      baseline2 = GetBaseline( pulse2, 5 ,50, pulseName2);
      baseline3 = GetBaseline( pulse3, 5 ,50, pulseName3);
      baseline4 = GetBaseline( pulse4, 5 ,50, pulseName4);

      // std::cout<<"Baseline "<<iEntry<<" "<<baseline1<<" "<<baseline2<<std::endl;
      // std::cout<<"index min "<<iEntry<<" "<<index_min1<<" "<<index_min2<<" "<<ch2Amp<<std::endl;
            
      // Correct pulse shape for baseline offset
      for(int j = 0; j < 1024; j++)
      	{
      	  Channel1VoltagesRaw_[j] = Channel1VoltagesRaw_[j] - baseline1;
      	  Channel2VoltagesRaw_[j] = Channel2VoltagesRaw_[j] - baseline2;
	  Channel3VoltagesRaw_[j] = Channel3VoltagesRaw_[j] - baseline3;
      	  Channel4VoltagesRaw_[j] = Channel4VoltagesRaw_[j] - baseline4;
      	}
      
      delete pulse1;
      delete pulse2;
      delete pulse3;
      delete pulse4;

      //----------------------------------------
      //
      //----------------------------------------

      
      //----------------------------------------------------------------------------------------
      //Create baseline corrected TGraphs (make sure you invert your pulse if they are negative)
      //----------------------------------------------------------------------------------------
      pulse1 = new TGraphErrors( GetTGraph( Channel1VoltagesRaw_, ti1_, false ) );
      pulse2 = new TGraphErrors( GetTGraph( Channel2VoltagesRaw_, ti2_, true ) );
      pulse3 = new TGraphErrors( GetTGraph( Channel3VoltagesRaw_, ti3_, false ) );
      pulse4 = new TGraphErrors( GetTGraph( Channel4VoltagesRaw_, ti4_, false ) );
      
      // if (doFilter) {
      // 	pulse1 = GetTGraphFilter( Channel1VoltagesRaw_, ti1_, pulseName1 , false);
      // 	pulse2 = GetTGraphFilter( Channel2VoltagesRaw_, ti2_, pulseName2 , false);
      // }

      //------------------------------------------------------------------
      //Getting index to maximum or minimum dependending on signal polarity
      //------------------------------------------------------------------
      //Find the absolute maximum. This is only used as a rough determination to decide if we'll use the early time samples
      //or the late time samples to do the baseline fit
      //NOTE: if your pulse is negative set _findMin (last input in the FindMaxAbsolute function) flag to <true>
      int index_min1 = FindMaxAbsolute(1024, Channel1VoltagesRaw_, false); // return index of the max
      int index_min2 = FindMaxAbsolute(1024, Channel2VoltagesRaw_, true); // return index of the max
      int index_min3 = FindMaxAbsolute(1024, Channel3VoltagesRaw_, false); // return index of the max
      int index_min4 = FindMaxAbsolute(1024, Channel4VoltagesRaw_, false); // return index of the max

      //Compute Amplitude : use units V
      Double_t tmpAmp = 0.0;
      Double_t tmpMin = 0.0;
      pulse1->GetPoint(index_min1, tmpMin, tmpAmp);
      ch1Amp = tmpAmp;
      pulse2->GetPoint(index_min2, tmpMin, tmpAmp);
      ch2Amp = tmpAmp;
      pulse3->GetPoint(index_min3, tmpMin, tmpAmp);
      ch3Amp = tmpAmp;
      pulse4->GetPoint(index_min4, tmpMin, tmpAmp);
      ch4Amp = tmpAmp;  
      
      //Get Pulse Integral
      if ( index_min1 != 0 ) ch1Int = GetPulseIntegral( index_min1 , Channel1VoltagesRaw_, "full");
      else ch1Int = 0.0;
      


      if ( index_min2 != 0 ) ch2Int = GetPulseIntegral( index_min2 , Channel2VoltagesRaw_, "full");
      else ch2Int = 0.0;
      

      if ( index_min3 != 0 ) ch3Int = GetPulseIntegral( index_min3 , Channel3VoltagesRaw_, "full");
      else ch3Int = 0.0;

      if ( index_min4 != 0 ) ch4Int = GetPulseIntegral( index_min4 , Channel4VoltagesRaw_, "full");
      else ch4Int = 0.0;
      


      //----------------
      // Gauss TimeStamp
      //----------------
      ch2Time_gausfitroot = GausFit_MeanTime( pulse2, index_min2, 3, 3);


      //---------------------
      // RisingEdge TimeStamp
      //---------------------
      float fs1[5];
      float fs2[5];
      float fs3[5];
      float fs4[5];
      
      RisingEdgeFitTime( pulse1, index_min1, fs1, iEntry, "linearFit_" + pulseName1, false);
      RisingEdgeFitTime( pulse2, index_min2, fs2, iEntry, "linearFit_" + pulseName2, false);
      RisingEdgeFitTime( pulse3, index_min3, fs3, iEntry, "linearFit_" + pulseName3, false);
      RisingEdgeFitTime( pulse4, index_min4, fs4, iEntry, "linearFit_" + pulseName4, false);
      ch1THM = fs1[3];
      ch2THM = fs2[3];
      ch3THM = fs3[3];
      ch4THM = fs4[3];
       
      //-------------------
      //for debugging the fits visually
      //--------------------

      /*
      TCanvas* c = new TCanvas("c","c",600,600);
      pulse2->GetXaxis()->SetRange(0,1024);
      pulse2->SetMarkerStyle(20);
      pulse2->Draw("AP");
      c->SaveAs("pulse1.pdf");
      */
      
      delete pulse1;
      delete pulse2;
      delete pulse3;
      delete pulse4;
      
      treeOut->Fill();
    }
  
  
  treeOut->Write();
  fout->Write();
  fout->Close();
}

TGraphErrors GetTGraph(  float* channel, float* time, bool invert )
{		
  //Setting Errors
  float errorX[1024], errorY[1024], channelFloat[1024];
  float _errorY = 0.00; //5%error on Y
  for ( int i = 0; i < 1024; i++ )
    {
      errorX[i]       = .0;
      errorY[i]       = _errorY*channel[i];
      if ( invert ) channelFloat[i] = -channel[i];
      else channelFloat[i] = channel[i];
    }
  //TGraphErrors* tg = new TGraphErrors( 1024, time, channelFloat, errorX, errorY );
  TGraphErrors tg( 1024, time, channelFloat, errorX, errorY );
  return tg;
};

int FindMaxAbsolute( int n, float *a, bool _findMin ) {
  
  if (n <= 0 || !a) return -1;
  int loc = 0;
  if ( _findMin )
    {
      float xmin = a[5];
      for  (int i = 5; i < n-10; i++) {
	// if (xmin > a[i] && a[i+1] < 0.5*a[i] && a[i] < -40. )  
	if (xmin > a[i] && a[i+1] < 0.5*a[i] && a[i] < -40. / 4096. )  
	  {
	    //std::cout << i << " " << a[i] << std::endl;
	    xmin = a[i];
	    loc = i;
	    //if ( a[i+5]>a[i] && a[i+10]>a[i+5] ) {
	    //break;
	  }
      }
    }
  else
    {
      float xmax = a[5];
      for  ( int i = 5; i < n-10; i++ )
	{
	  if ( a[i] > xmax && a[i+1] < 0.9*a[i] && a[i] > 40. / 4096. )  
	    {
	      //std::cout << i << " " << a[i] << std::endl;
	      xmax = a[i];
	      loc = i;
	      //if ( a[i+5]>a[i] && a[i+10]>a[i+5] ) {
	      //break;
	    }
	}
    }
  return loc;
}

float GetBaseline(TGraphErrors * pulse, int i_low, int i_high, TString fname )
{
  double x_low, x_high, y, dummy;
  pulse->GetPoint(i_low, x_low, y);
  pulse->GetPoint(i_high, x_high, y);
  
  TF1* flinear = new TF1("flinear","[0]", x_low, x_high );
  
  pulse->Fit("flinear","RQ","", x_low, x_high );
  
  /* std::cout << "make plot" << std::endl;
  std::cout << x_low << x_high << fname << std::endl;
  TCanvas* c = new TCanvas("canvas","canvas",800,400) ;
  pulse->GetXaxis()->SetLimits(x_low-3, x_high+3);
  pulse->SetMarkerSize(1);
  pulse->SetMarkerStyle(20);
  pulse->Draw("AP");
  c->SaveAs(fname+"LinearFit.pdf"); */
  
  float a = flinear->GetParameter(0);
  delete flinear;
  
  return a;
}


float GetPulseIntegral(int peak, float *a, std::string option) 
{
  float integral = 0.;

  if (option == "full") {
    for (int i=5; i < 1024; i++) {
      integral += a[i] * 0.2 * 1e-9 * (1.0/50.0) * 1e12; //in units of pC, for 50Ohm termination
    }
  }
  else {
    for (int i=peak-20; i < peak+25; i++) {
      integral += a[i] * 0.2 * 1e-9 * (1.0/50.0) * 1e12; //in units of pC, for 50Ohm termination
    }
  }
  return -1.0 * integral;

}

TGraphErrors* GetTGraphFilter( float* channel, float* time, TString pulseName, bool makePlot )
{
  float Gauss[1024];
  //Setting Errors
  float errorX[1024], errorY[1024], channelFloat[1024];
  float _errorY = 0.00; //5%error on Y
  
  for ( int i = 0; i < 1024; i++ )
    {
      
      errorX[i]       = .0;
      errorY[i]       = _errorY*channel[i];
      channelFloat[i] = -channel[i];
    }
  
  TF1 *fb = new TF1("fb","gaus(0)", 0.0, 204.6);
  fb->SetParameter(1, 100);
  float sigma =1.0;
  fb->SetParameter(2, sigma);
  fb->SetParameter(0, 1/(sqrt(3.1415*2.0)*sigma) );
  //eval Gaussian
  float step = 0.2;//200ps
  for ( int i = 0; i < 1024; i++ )
    {
      Gauss[i] = fb->Eval( float(i)*step );
    }
  
  float channelFloatFiltered[2048];
  for ( int i = 0; i < 2048; i++ )
    {
      float convolvedPoint = 0;
      for ( int j = 0; j <= i; j++ )
	{
	  if ( i < 1024 )
	    {
	      convolvedPoint += channelFloat[i-j]*Gauss[1023-j];
	    }
	  else
	    {
	      if ( 1023-(i-1023)-j >= 0 ) convolvedPoint += channelFloat[1023-j]*Gauss[1023-(i-1023)-j];
	    }
	}
      //if ( i < 1024 ) channelFloatFiltered[i] = convolvedPoint;
      channelFloatFiltered[i] = convolvedPoint;
    }
  
  float channelFloatFilteredFix[1024];
  for ( int i = 0; i < 1024; i++ )
    {
      channelFloatFilteredFix[i] = 0.2*channelFloatFiltered[i+523];
    }
  
  TGraphErrors* tg = new TGraphErrors( 1024, time, channelFloat, errorX, errorY );
  TGraphErrors* tg2 = new TGraphErrors( 1024, time, channelFloatFilteredFix, errorX, errorY );

  if (makePlot) {
    TCanvas* c = new TCanvas("canvas","canvas",800,400) ;         
    tg2->GetXaxis()->SetLimits(50, 70);
    tg->GetXaxis()->SetLimits(50, 70);
    //tg2->Fit("fb","","", 0.0, 204.6 );
    tg2->SetMarkerSize(0.5);
    tg->SetMarkerSize(0.5);
    tg2->SetMarkerStyle(20);
    tg->SetMarkerStyle(20);
    tg2->Draw("AP");
    tg2->SetMarkerColor(kBlue);
    tg->Draw("sameP");
    c->SaveAs(pulseName + "GausPulse.pdf");
  }
  return tg;
};

int FindRealMin( int n, float *a) {  
  if (n <= 0 || !a) return -1;
  float xmin = a[5];
  int loc = 0;
  
  float noise = 0;
  
  for ( int i = 5; i < 50; i++)
    {
      if( std::abs(a[i]) > noise ) 
	{
	  noise = std::abs(a[i]);
	}
    }

  for  (int i = 5; i < n-10; i++) {
    if (xmin > a[i] && a[i+1] < 0.5*a[i] && a[i] < -3*noise && a[i] < -50./4096. )  
      {
	//std::cout << a[i] << std::endl;
	xmin = a[i];
	loc = i;
	//if ( a[i+5]>a[i] && a[i+10]>a[i+5] ) {
	//break;
      }
  }
 
  float xmin_init = xmin;
  float xmin_new = a[5];
  int loc_new = loc;

  bool stop = false;
  while( !stop )
    {
      for ( int i = 5; i < loc_new -25; i++ )
        {
          if ( a[i] < xmin_new && 0.5*a[i] > a[i+1] && a[i] < 0.15* xmin_init )
            {
              xmin_new = a[i];
              loc_new = i;
            }
        }

      xmin_init = xmin_new;

      if( loc_new == loc ) break;
      //std::cout << "new peak @ " << loc_new << ", ymin: " << xmin_new << std::endl;
      if ( xmin_new > -2*noise || xmin_new > -40 / 4096. ) loc_new = 0;
      xmin_new = a[5];
      loc = loc_new;
    }

  //std::cout << "LOC2: " << loc << std::endl;                                                                                                                               
  /*                                                                
  while ( xmin_init != xmin_new ) {
    for (int i = 5; i < loc - 50; i++) {
      if (xmin_new > a[i] && a[i+1] < 0.5*a[i] && a[i] < xmin_init*2/3 )  {
        xmin_new = a[i];
        loc = i;
      }
    }
    xmin_init = xmin_new
    xmin_new = a[5]
  }
  */
  return loc_new;
}

float GausFit_MeanTime(TGraphErrors* pulse, const int index_min, const int index_first, const int index_last)
{
  double x_low, x_high, y;
  pulse->GetPoint(index_min-index_first, x_low, y);
  pulse->GetPoint(index_min+index_last, x_high, y);
  
  TF1* fpeak = new TF1("fpeak","gaus", x_low, x_high);
  pulse->Fit("fpeak","Q","", x_low, x_high);
  float timepeak = fpeak->GetParameter(1);
  delete fpeak;
  
  return timepeak;
}

void RisingEdgeFitTime(TGraphErrors * pulse, const float index_min, float* tstamp, int event, TString fname, bool makePlot )
{
  double x_low, x_high, y, dummy;
  double ymax;
  pulse->GetPoint(index_min, x_low, ymax);
  for ( int i = 1; i < 100; i++ )
    {
      pulse->GetPoint(index_min-i, x_low, y);
      // std::cout<<"BNNNBBB  "<<event<<" "<<ymax<<" "<<x_low<<" "<<y<<std::endl;
      if ( y < 0.15*ymax ) break;
      // if ( y < 0.2*ymax ) break;
   }

  for ( int i = 1; i < 100; i++ )
    {
      pulse->GetPoint(index_min-i, x_high, y);
      if ( y < 0.95*ymax ) break;
    }

 

  // if(x_low_1<x_low) x_low = x_low_1;
  // if(x_high_1>x_high) x_high = x_high_1;

  //----------------------------------------------------------------
  //The following is the best configuration for picsec laser trigger
  //----------------------------------------------------------------
  pulse->GetPoint(index_min-8, x_low, y);
  pulse->GetPoint(index_min-2, x_high, y);
  pulse->GetPoint(index_min, dummy, y);
  
  TF1* flinear = new TF1("flinear","[0]*x+[1]", x_low, x_high );
  float max = -9999;
  double* yy = pulse->GetY();
  
  for ( int i = 0; i < 1024; i++ )
    {
      if ( yy[i] > max ) max = yy[i];
    }
  //std::cout << "max: " << max << std::endl;

  /*if( max < 10 || index_min < 0 || index_min > 1023 )
    {
      std::cout << "DEB: skipping event--> " << event << std::endl;
      return;
    }
  */

  // if(std::abs(x_high-x_low) < 0.8) {x_low -= 0.1; x_high += 0.2;}
  
  pulse->Fit("flinear","Q","", x_low, x_high );
  double slope = flinear->GetParameter(0);
  double b     = flinear->GetParameter(1);
  
  if ( makePlot )
    {
      std::cout << "make plot" << std::endl;
      TCanvas* c = new TCanvas("canvas","canvas",800,400) ;
      pulse->GetXaxis()->SetLimits(x_low-10, x_high+10);
      pulse->SetMarkerSize(0.3);
      pulse->SetMarkerStyle(20);
      pulse->Draw("AP");
      c->SaveAs(fname+"LinearFit.pdf");
      //delete c;
    }
  tstamp[0] = (0.0*y-b)/slope;
  //std::cout << "----" << tstamp[0]  << std::endl;
  tstamp[1] = (0.15*y-b)/slope;
  tstamp[2] = (0.30*y-b)/slope;
  tstamp[3] = (0.40*y-b)/slope;//Best Configuration for picsec laser trigger
  tstamp[4] = (0.60*y-b)/slope;
  
  delete flinear;
};
