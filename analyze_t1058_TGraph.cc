#include <iostream>
#include <sys/stat.h>
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
int FindFirstMaximum( int n, float *a, bool _findMin = false );
bool DetectDoublePeak( int n, float *a, bool _findMin = false );

TGraphErrors* GetTGraphFilter( float* channel, float* time, TString pulseName, bool makePlot );
float GetPulseIntegral(int peak, float *a, float *t, std::string option, int peakmin, int peakmax);
float GetBaseline(TGraphErrors * pulse, int i_low, int i_high, TString fname );
TGraphErrors GetTGraph(  float* channel, float* time, bool invert = false );
int FindRealMin( int n, float *a);
float RisingEdgeFitTime(TGraphErrors* pulse, const float index_min, const float lowFraction, const float highFraction, float* tstamp, int event, TString fname, bool makePlot, bool trigger = false);
void FitDirectHitPlusScintillationSignal(TGraphErrors* pulse, const float index_min, float* result, int event, TString fname, bool makePlot);

float LED( TH1F * pulse, double threshold, int nsamples, int splineBinFactor );
float LinearFit_Baseline(TH1F * pulse, const int index_min, const int range);
float LinearFit_Intercept(TH1F * pulse, const float base, const int index_first, const int index_last);
float GausFit_MeanTime(TGraphErrors* pulse, const int index_min, const int index_first, const int index_last, float * fit_result, TString pulseName, bool makePlot );
void FitRisingEdge(TH1F* pulse, int nbinsL, int nbinsH, float &THM, float &risetime, float base);
void FitFullPulse(TH1F* pulse, float &par0, float &par1, float &par2);
TH1F* InterpolateWaveform(int nsamples, float* outputwaveform, float *inputwaveform, int splineBinFactor, std::string name);

const int Nsamples = 1024;

int main (int argc, char **argv)
{
  TFile *f;

  float f_xpos = 0.0;//in mm
  float f_ypos = 0.0;//in mm
  float V_bias = 0.0;//bias voltage, in V

  std::string runNum;
  std::string inputDir;
  
  if (argc >= 3)
    {
      f = new TFile(argv[1]);

      std::cout << ">> Opening file " << argv[1] << " ......" << std::endl;
      if (!f->IsOpen())
	{			// terminate if the file can't be opened
	  std::cerr << "!! File open error:" << argv[1] << std::endl;
	  return 1;
	}

    std::string fn(argv[1]);
    int pf = fn.find(".root");
    int pi = fn.rfind("/")+1;
    int di = fn.rfind("/data/");
    runNum = fn.substr(pi, pf-pi);
    if(di != std::string::npos) inputDir = fn.substr(0, di+1);
    else inputDir = fn.substr(0, pi);
    
    std::cout<<"input Dir: "<<inputDir<<std::endl;
    std::cout<<"runNum: "<<runNum<<std::endl;

    mkdir((inputDir+"pulses").c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    }
  else
    {				// terminate if there is no input file or more than 1 input file
      std::cerr << "!! No input file" << std::endl;
      return 1;
    }

  
  if (argc >= 5)
  {
	//std::string s_xpos(argv[3]);
  	//std::string s_ypos(argv[4]);
	f_xpos = float(atof(argv[3]));
	f_ypos = float(atof(argv[4]));
  }

  if (argc >= 6)
  {
	V_bias = float(atof(argv[5]));
  }


  bool includePulseshapeInOutput = false;
  if (argc >= 7) includePulseshapeInOutput = bool(atoi(argv[6]));
  
  //const int splineBinFactor = 40;
  const int splineBinFactor = 1;

  int t_[Nsamples]; for (int i = 0; i < 1024; i++) { t_[i] = i; }
  float Channel1VoltagesRaw_[Nsamples];
  float Channel2VoltagesRaw_[Nsamples];
  float Channel3VoltagesRaw_[Nsamples];
  float Channel4VoltagesRaw_[Nsamples];

  float Channel1VoltagesRawNorm_[Nsamples];
  float Channel2VoltagesRawNorm_[Nsamples];
  float Channel3VoltagesRawNorm_[Nsamples];
  float Channel4VoltagesRawNorm_[Nsamples];


  float Channel1Voltages_[Nsamples*splineBinFactor];
  float Channel2Voltages_[Nsamples*splineBinFactor];
  float Channel3Voltages_[Nsamples*splineBinFactor];
  float Channel4Voltages_[Nsamples*splineBinFactor];

  float ti1_[Nsamples], ti2_[Nsamples], ti3_[Nsamples], ti4_[Nsamples];
  float tishift1_[Nsamples], tishift2_[Nsamples], tishift3_[Nsamples], tishift4_[Nsamples];

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
    //fn = "AnaFiles/" + fn.substr(pi, pf-pi) + "_ana.root";
    if(includePulseshapeInOutput) fn = fn.substr(0, pf) + "_ana_withpulse.root";
    else fn = fn.substr(0, pf) + "_ana.root";
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

  float ch1LinearTime15 = 0.0;// rising edge 15%
  float ch2LinearTime15 = 0.0;// rising edge 15%
  float ch3LinearTime15 = 0.0;// rising edge 15%
  float ch4LinearTime15 = 0.0;// rising edge 15%

  float ch1LinearTime30 = 0.0;// rising edge 15%
  float ch2LinearTime30 = 0.0;// rising edge 15%
  float ch3LinearTime30 = 0.0;// rising edge 15%
  float ch4LinearTime30 = 0.0;// rising edge 15%

  float ch1LinearTime40 = 0.0;// rising edge 15%
  float ch2LinearTime40 = 0.0;// rising edge 15%
  float ch3LinearTime40 = 0.0;// rising edge 15%
  float ch4LinearTime40 = 0.0;// rising edge 15%

  float ch1LinearTime60 = 0.0;// rising edge 15%
  float ch2LinearTime60 = 0.0;// rising edge 15%
  float ch3LinearTime60 = 0.0;// rising edge 15%
  float ch4LinearTime60 = 0.0;// rising edge 15%

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
 
  float ch1Int_gauspeak = 0;
  float ch2Int_gauspeak = 0;
  float ch3Int_gauspeak = 0;
  float ch4Int_gauspeak = 0;
  float ch5Int_gauspeak = 0;
  float ch6Int_gauspeak = 0;
  float ch7Int_gauspeak = 0;
  float ch8Int_gauspeak = 0;
  
  unsigned int ch1QualityBit = 0;
  unsigned int ch2QualityBit = 0;
  unsigned int ch3QualityBit = 0;
  unsigned int ch4QualityBit = 0;
  
  float ch1chisq = -1;
  float ch2chisq = -1;
  float ch3chisq = -1;
  float ch4chisq = -1;
  


  treeOut->Branch("event",&eventNumber,"event/i");
  treeOut->Branch("x",&f_xpos,"x/F");
  treeOut->Branch("y",&f_ypos,"y/F");
  treeOut->Branch("V_bias",&V_bias,"V_bias/F");
  treeOut->Branch("t1gausroot",&ch1Time_gausfitroot,"t1gausroot/F");
  treeOut->Branch("t2gausroot",&ch2Time_gausfitroot,"t2gausroot/F");
  treeOut->Branch("t3gausroot",&ch3Time_gausfitroot,"t3gausroot/F");
  treeOut->Branch("t4gausroot",&ch4Time_gausfitroot,"t4gausroot/F");
  
  treeOut->Branch("ch1Amp",&ch1Amp,"ch1Amp/F");
  treeOut->Branch("ch2Amp",&ch2Amp,"ch2Amp/F");
  treeOut->Branch("ch3Amp",&ch3Amp,"ch3Amp/F");
  treeOut->Branch("ch4Amp",&ch4Amp,"ch4Amp/F");

  treeOut->Branch("ch1LinearTime15",&ch1LinearTime15,"ch1LinearTime15/F");
  treeOut->Branch("ch2LinearTime15",&ch2LinearTime15,"ch2LinearTime15/F");
  treeOut->Branch("ch3LinearTime15",&ch3LinearTime15,"ch3LinearTime15/F");
  treeOut->Branch("ch4LinearTime15",&ch4LinearTime15,"ch4LinearTime15/F");

  treeOut->Branch("ch1LinearTime30",&ch1LinearTime30,"ch1LinearTime30/F");
  treeOut->Branch("ch2LinearTime30",&ch2LinearTime30,"ch2LinearTime30/F");
  treeOut->Branch("ch3LinearTime30",&ch3LinearTime30,"ch3LinearTime30/F");
  treeOut->Branch("ch4LinearTime30",&ch4LinearTime30,"ch4LinearTime30/F");

  treeOut->Branch("ch1LinearTime40",&ch1LinearTime40,"ch1LinearTime40/F");
  treeOut->Branch("ch2LinearTime40",&ch2LinearTime40,"ch2LinearTime40/F");
  treeOut->Branch("ch3LinearTime40",&ch3LinearTime40,"ch3LinearTime40/F");
  treeOut->Branch("ch4LinearTime40",&ch4LinearTime40,"ch4LinearTime40/F");

  treeOut->Branch("ch1LinearTime60",&ch1LinearTime60,"ch1LinearTime60/F");
  treeOut->Branch("ch2LinearTime60",&ch2LinearTime60,"ch2LinearTime60/F");
  treeOut->Branch("ch3LinearTime60",&ch3LinearTime60,"ch3LinearTime60/F");
  treeOut->Branch("ch4LinearTime60",&ch4LinearTime60,"ch4LinearTime60/F");

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
 
  treeOut->Branch("ch1Int_gauspeak",&ch1Int_gauspeak,"ch1Int_gauspeak/F");
  treeOut->Branch("ch2Int_gauspeak",&ch2Int_gauspeak,"ch2Int_gauspeak/F");
  treeOut->Branch("ch3Int_gauspeak",&ch3Int_gauspeak,"ch3Int_gauspeak/F");
  treeOut->Branch("ch4Int_gauspeak",&ch4Int_gauspeak,"ch4Int_gauspeak/F");
  
  treeOut->Branch("ch1chisq",&ch1chisq,"ch1chisq/F");
  treeOut->Branch("ch2chisq",&ch2chisq,"ch2chisq/F");
  treeOut->Branch("ch3chisq",&ch3chisq,"ch3chisq/F");
  treeOut->Branch("ch4chisq",&ch4chisq,"ch4chisq/F");  

  if (includePulseshapeInOutput) {
    
    treeOut->Branch("c1",Channel1VoltagesRaw_,"c1[1024]/F");
    treeOut->Branch("c2",Channel2VoltagesRaw_,"c2[1024]/F");
    treeOut->Branch("c3",Channel3VoltagesRaw_,"c3[1024]/F");
    treeOut->Branch("c4",Channel4VoltagesRaw_,"c4[1024]/F");

    treeOut->Branch("c1norm",Channel1VoltagesRawNorm_,"c1norm[1024]/F");
    treeOut->Branch("c2norm",Channel2VoltagesRawNorm_,"c2norm[1024]/F");
    treeOut->Branch("c3norm",Channel3VoltagesRawNorm_,"c3norm[1024]/F");
    treeOut->Branch("c4norm",Channel4VoltagesRawNorm_,"c4norm[1024]/F");

    treeOut->Branch("t",t_,"t[1024]/I");
    
    treeOut->Branch("ti1",ti1_,"ti1[1024]/F");
    treeOut->Branch("ti2",ti2_,"ti2[1024]/F");
    treeOut->Branch("ti3",ti3_,"ti3[1024]/F");
    treeOut->Branch("ti4",ti4_,"ti4[1024]/F");

    treeOut->Branch("tishift1",tishift1_,"tishift1[1024]/F");
    treeOut->Branch("tishift2",tishift2_,"tishift2[1024]/F");
    treeOut->Branch("tishift3",tishift3_,"tishift3[1024]/F");
    treeOut->Branch("tishift4",tishift4_,"tishift4[1024]/F");

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
      pulse3 = new TGraphErrors( GetTGraph( Channel3VoltagesRaw_, ti3_, true ) );
      pulse4 = new TGraphErrors( GetTGraph( Channel4VoltagesRaw_, ti4_, true ) );
      
      // if (doFilter) {
      //pulse1 = GetTGraphFilter( Channel1VoltagesRaw_, ti1_, Form("myPulseFilter%d",iEntry) , true);
      // 	pulse2 = GetTGraphFilter( Channel2VoltagesRaw_, ti2_, pulseName2 , false);
      // }

      //pulse3 = GetTGraphFilter( Channel3VoltagesRaw_, ti3_, Form("myPulse3Filter%d",iEntry) , false);
      

      
      //------------------------------------------------------------------
      //Getting index to maximum or minimum dependending on signal polarity
      //------------------------------------------------------------------
      //Find the absolute maximum. This is only used as a rough determination to decide if we'll use the early time samples
      //or the late time samples to do the baseline fit
      //NOTE: if your pulse is negative set _findMin (last input in the FindMaxAbsolute function) flag to <true>
      //std::cout << "=====event " << iEntry+1 << "==========" << std::endl;	
      int index_min1 = FindMaxAbsolute(1024, Channel1VoltagesRaw_, false); // return index of the max
      int index_min2 = FindMaxAbsolute(1024, Channel2VoltagesRaw_, true); // return index of the max
      int index_min3 = FindMaxAbsolute(1024, Channel3VoltagesRaw_, true); // return index of the max
      int index_min4 = FindMaxAbsolute(1024, Channel4VoltagesRaw_, true); // return index of the max

      //std::cout << "event: " << iEntry << " --> " << ti3_[index_min3] << " " << Channel3VoltagesRaw_[index_min3] << std::endl;
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
      if ( index_min1 != 0 ) ch1Int = GetPulseIntegral( index_min1 , Channel1VoltagesRaw_, ti1_, "part", 20, 15);
      else ch1Int = 0.0;

      if ( index_min2 != 0 ) 
	{
		ch2Int = GetPulseIntegral( index_min2 , Channel2VoltagesRaw_,ti2_,  "part", 20, 150);
		ch2Int_gauspeak = 2.0 * GetPulseIntegral( index_min2 , Channel2VoltagesRaw_,ti2_,  "part", 20, 0);
	}
      else ch2Int = 0.0;
      
      if ( index_min3 != 0 ) ch3Int = GetPulseIntegral( index_min3 , Channel3VoltagesRaw_, ti3_,  "full", 6, 10);
      else ch3Int = 0.0;

      if ( index_min4 != 0 ) ch4Int = GetPulseIntegral( index_min4 , Channel4VoltagesRaw_, ti4_, "full",35, 75);
      else ch4Int = 0.0;

      //----------------
      // Gauss TimeStamp
      //----------------

      float fit_result1[3];   
      float fit_result2[3];   
 
      ch1Time_gausfitroot = GausFit_MeanTime( pulse1, index_min1, 4, 3, fit_result1, pulseName1, false);
      ch2Time_gausfitroot = GausFit_MeanTime( pulse2, index_min2, 3, 3, fit_result2, pulseName2, false);
      // ch3Time_gausfitroot = GausFit_MeanTime( pulse3, index_min3, 3, 3, pulseName3, false);
      // ch4Time_gausfitroot = GausFit_MeanTime( pulse4, index_min4, 8, 8, pulseName3, false);
      
      //ch1Int_gauspeak = sqrt(2.0*3.141592653) * fit_result1[2] * fit_result1[0] * 1e-9 * (1.0/50.0) * 1e12; //in units of pC, for 50Ohm termination
      //ch2Int_gauspeak = sqrt(2.0*3.141592653) * fit_result2[2] * fit_result2[0] * 1e-9 * (1.0/50.0) * 1e12; //in units of pC, for 50Ohm termination

      //shift the pulse by subtracting the refrence time
      
      for(int i=0;i<Nsamples;i++)
      {
	if(ch1Time_gausfitroot > 0 && ch1Time_gausfitroot < 200)
	{
		tishift1_[i] = ti1_[i] -  ch1Time_gausfitroot;
		tishift2_[i] = ti1_[i] -  ch1Time_gausfitroot;
		tishift3_[i] = ti1_[i] -  ch1Time_gausfitroot;
		tishift4_[i] = ti1_[i] -  ch1Time_gausfitroot;
	}
	if(ch1Amp>0) Channel1VoltagesRawNorm_[i] = Channel1VoltagesRaw_[i]/ch1Amp;
	if(ch2Amp>0) Channel2VoltagesRawNorm_[i] = Channel2VoltagesRaw_[i]/ch2Amp;
	if(ch3Amp>0) Channel3VoltagesRawNorm_[i] = Channel3VoltagesRaw_[i]/ch3Amp;
	if(ch4Amp>0) Channel4VoltagesRawNorm_[i] = Channel4VoltagesRaw_[i]/ch4Amp;
      }
      //---------------------
      // RisingEdge TimeStamp
      //---------------------
      float fs1[5];
      float fs2[5];
      float fs3[5];
      float fs4[5];
      

      ch1Risetime = RisingEdgeFitTime( pulse1, index_min1, 0.15, 0.95, fs1, iEntry, "linearFit_" + pulseName1, false, true);
      ch2Risetime = RisingEdgeFitTime( pulse2, index_min2, 0.1, 0.9, fs2, iEntry, "linearFit_" + pulseName2, false, false);
      //ch3Risetime = RisingEdgeFitTime( pulse3, index_min3, 0.1, 0.9, fs3, iEntry, "linearFit_" + pulseName3, false, false);
      //ch4Risetime = RisingEdgeFitTime( pulse4, index_min4, 0.1, 0.4, fs4, iEntry, "linearFit_" + pulseName4, false);


      ch1LinearTime15 = fs1[1];
      ch2LinearTime15 = fs2[1];
      ch3LinearTime15 = fs3[1];
      ch4LinearTime15 = fs4[1];

      ch1LinearTime30 = fs1[2];
      ch2LinearTime30 = fs2[2];
      ch3LinearTime30 = fs3[2];
      ch4LinearTime30 = fs4[2];

      ch1LinearTime40 = fs1[3];
      ch2LinearTime40 = fs2[3];
      ch3LinearTime40 = fs3[3];
      ch4LinearTime40 = fs4[3];

      ch1LinearTime60 = fs1[4];
      ch2LinearTime60 = fs2[4];
      ch3LinearTime60 = fs3[4];
      ch4LinearTime60 = fs4[4];

      ch1THM = fs1[2];
      ch2THM = fs2[2];
      ch3THM = fs3[3];
      ch4THM = fs4[1];
      
      //-------------------
      //for debugging the fits visually
      //--------------------

     
      
      if(iEntry+1<=20){
      TCanvas* c = new TCanvas("c","c",600,600);
      pulse1->GetXaxis()->SetRangeUser(0,200);
      pulse1->SetMarkerStyle(20);
      pulse1->Draw("AP");
      c->SaveAs(Form("%s/pulses/%s_pulse1_event%lld.pdf", inputDir.c_str(), runNum.c_str(), iEntry+1));
      c->SaveAs(Form("%s/pulses/%s_pulse1_event%lld.png", inputDir.c_str(), runNum.c_str(), iEntry+1));
      c->SaveAs(Form("%s/pulses/%s_pulse1_event%lld.C", inputDir.c_str(), runNum.c_str(), iEntry+1));

      //pulse2->GetXaxis()->SetRange(0,200);
      pulse2->SetMarkerStyle(20);
      pulse2->GetXaxis()->SetRangeUser(50,70);
      pulse2->Draw("AP");
      c->SaveAs(Form("%s/pulses/%s_pulse2_event%lld.pdf", inputDir.c_str(), runNum.c_str(), iEntry+1));
      c->SaveAs(Form("%s/pulses/%s_pulse2_event%lld.png", inputDir.c_str(), runNum.c_str(), iEntry+1));
      c->SaveAs(Form("%s/pulses/%s_pulse2_event%lld.C", inputDir.c_str(), runNum.c_str(), iEntry+1));
      
      pulse3->SetMarkerStyle(20);
      pulse3->GetXaxis()->SetRangeUser(30,90);
      // pulse3->GetYaxis()->SetRangeUser(-0.03,0.04);
      pulse3->Draw("AP");
      c->SaveAs(Form("%s/pulses/%s_pulse3_event%lld.pdf", inputDir.c_str(), runNum.c_str(), iEntry+1));
      c->SaveAs(Form("%s/pulses/%s_pulse3_event%lld.png", inputDir.c_str(), runNum.c_str(), iEntry+1));
      c->SaveAs(Form("%s/pulses/%s_pulse3_event%lld.C", inputDir.c_str(), runNum.c_str(), iEntry+1));
      
      pulse4->SetMarkerStyle(20);
      pulse4->GetXaxis()->SetRangeUser(0,200);
      pulse4->Draw("AP");
      c->SaveAs(Form("%s/pulses/%s_pulse4_event%lld.pdf", inputDir.c_str(), runNum.c_str(), iEntry+1));
      c->SaveAs(Form("%s/pulses/%s_pulse4_event%lld.png", inputDir.c_str(), runNum.c_str(), iEntry+1));
      c->SaveAs(Form("%s/pulses/%s_pulse4_event%lld.C", inputDir.c_str(), runNum.c_str(), iEntry+1));
      }
      
      // if(iEntry+1>10000) break;
      
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

bool DetectDoublePeak( int n, float *a, bool _findMin )
{
  if (n <= 0 || !a) return false;
  int loc = 0;
  int ncounts = 0;
  if ( _findMin )
    {
      float xmin = a[5];
      for  (int i = 5; i < n-10; i++) {
	if ( a[i] < xmin && a[i+1] < 0.8*a[i] && a[i] < -10. / 4096. )  
	  {
	    xmin = a[i];
	    if ( loc == i-1 ) ncounts++;
	    else if ( loc != i-1 && ncounts >= 3 ) return true;
	    else ncounts = 0;
	    //std::cout << "loc: " <<  i << " " << a[i] << " " << ncounts << std::endl;
	    loc = i;
	  }
      }
    }
  else
    {
      float xmax = a[5];
      for  ( int i = 6; i < n-10; i++ )
	{
	  if ( a[i] > xmax && a[i+1] < 0.9*a[i] && a[i] > 10. / 4096. )  
	    {
	      xmax = a[i];
	      if ( loc == i-1 ) ncounts++;
	      else if ( loc != i-1 && ncounts >= 3 ) return true;
	      else ncounts = 0;
	      loc = i;
	    }
	}
    }
  
  return false;
}

int FindMaxAbsolute( int n, float *a, bool _findMin ) {
  
  if (n <= 0 || !a) return -1;
  int loc = 0;
  if ( _findMin )
    {
      float xmin = a[5];
      for  (int i = 5; i < n-10; i++) {
	  // to 2mV cut
	  if ( a[i] < xmin  && a[i+1] > 0.98*a[i] && a[i] < -0.002 )  
	  {
	    xmin = a[i];
	    loc = i;
	  }
      }
    }
  else
    {
      float xmax = a[5];
      for  ( int i = 5; i < n-10; i++ )
	{
	  // to 2 mV cut
	  if ( a[i] > xmax && a[i+1] < 0.9*a[i] && a[i] > 0.002 )  
	    {
	      //std::cout << i << " " << a[i] << std::endl;
	      xmax = a[i];
	      loc = i;
	    }
	}
    }
  return loc;
}

int FindFirstMaximum( int n, float *a, bool _findMin ) {
  
  if (n <= 0 || !a) return -1;
  int loc = 0;
  if ( _findMin )
    {
      float xmin = a[5];
      for  (int i = 5; i < n-10; i++) {
	  // to 2mV cut
	  if ( a[i] < xmin  && a[i+1] > 0.98*a[i] && a[i] < -0.02 )  
	  {
	    //std::cout << i << std::endl;
	    return i;
	  }
      }
    }
  else
    {
      float xmax = a[5];
      for  ( int i = 5; i < n-10; i++ )
	{
	  // to 2 mV cut
	  if ( a[i] > xmax && a[i+1] < 0.9*a[i] && a[i] > 0.02 )  
	    {
	      //std::cout << i << " " << a[i] << std::endl;
              return i;
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


float GetPulseIntegral(int peak, float *a, float *t, std::string option, int peakmin, int peakmax) 
{
  float integral = 0.;

  if (option == "full") {
    for (int i=5; i < 1024; i++) {
      integral += a[i] * 0.2 * 1e-9 * (1.0/50.0) * 1e12; //in units of pC, for 50Ohm termination
    }
  }
  else {
    //for (int i=peak-4; i < peak+4; i++) {
    //}
    for (int i=peak-peakmin; i < peak+peakmax; i++) {
    //for (int i = 389; i < 397; i++){
      //it makes more sense to do the integral around the peak rather than a fixed window
      //integral += a[i] * 0.2 * 1e-9 * (1.0/50.0) * 1e12; //in units of pC, for 50Ohm termination
      //trapezoid
      //integral += 0.5*(a[i]+a[i-1]) * (t[i]-t[i-1]) * 1e-9 * (1.0/50.0) * 1e12; //in units of pC, for 50Ohm termination
      
      //Simpson's Rule for equaled space-->Cartwright correction for unequaled space, only worked for odd points
      integral += ( (t[i+2]-t[i]) / 6.0 ) * ( ( 2-(t[i+2]-t[i+1])/(t[i+1]-t[i]) )* a[i] + (t[i+2]-t[i])*(t[i+2]-t[i])/((t[i+2]-t[i+1])*(t[i+1]-t[i])) * a[i+1] + ( 2-(t[i+1]-t[i])/(t[i+2]-t[i+1]) ) * a[i+2] ) * 1e-9 * (1.0/50.0) * 1e12; //in units of pC, for 50Ohm termination
      i++;
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
  float sigma =0.8;
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
    tg2->GetXaxis()->SetRangeUser(300, 700);
    tg->GetXaxis()->SetRangeUser(300, 700);
    //tg2->Fit("fb","","", 0.0, 204.6 );
    tg2->SetMarkerSize(0.5);
    tg->SetMarkerSize(0.5);
    tg2->SetMarkerStyle(20);
    tg->SetMarkerStyle(20);
    tg2->Draw("AP");
    tg2->GetYaxis()->SetRangeUser(0, 0.2);
    tg2->SetMarkerColor(kBlue);
    tg->Draw("sameP");
    tg->GetYaxis()->SetRangeUser(0, 0.5);
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

float GausFit_MeanTime(TGraphErrors* pulse, const int index_min, const int index_first, const int index_last, float * fit_result, TString pulseName, bool makePlot)
{
  double x_low, x_high, y;
  pulse->GetPoint(index_min-index_first, x_low, y);
  pulse->GetPoint(index_min+index_last, x_high, y);
  
  double x_dummy, ymax, ymax_minus1;
  pulse->GetPoint(index_min, x_dummy, ymax);
  pulse->GetPoint(index_min-1, x_dummy, ymax_minus1);
  if (ymax_minus1 > ymax)
    {
    pulse->GetPoint(index_min-index_first-1, x_low, y);
    }

  TF1* fpeak = new TF1("fpeak","gaus", x_low-.1, x_high+.1);
  pulse->Fit("fpeak","Q","", x_low-.1, x_high+.1);
  float timepeak = fpeak->GetParameter(1);

  fit_result[0] = fpeak->GetParameter(0);
  fit_result[1] = fpeak->GetParameter(1);
  fit_result[2] = fpeak->GetParameter(2);

  if ( makePlot )
    {
      std::cout << "make plot" << std::endl;
      TCanvas* c = new TCanvas("canvas","canvas",800,400) ;
      pulse->GetXaxis()->SetLimits(x_low-10, x_high+10);
      pulse->SetMarkerSize(1.5);
      pulse->SetMarkerStyle(20);
      pulse->Draw("AP");
      c->SaveAs(pulseName+"GausFit.pdf");
      //delete c;
    }
  
  delete fpeak;
  
  return timepeak;
}

float RisingEdgeFitTime(TGraphErrors* pulse, const float index_min, const float lowFraction, const float highFraction, float* tstamp, int event, TString fname, bool makePlot, bool trigger)
{
  double x_low, x_high, xdummy, y, dummy;
  double ymax;
  pulse->GetPoint(index_min, x_low, ymax);
  //std::cout << x_low << " " << ymax << std::endl;
  for ( int i = 1; i < 500; i++ )
    {
      pulse->GetPoint(index_min-i, x_low, y);
      if ( y < lowFraction*ymax ) break;
    }

  for ( int i = 1; i < 500; i++ )
    {
      pulse->GetPoint(index_min-i, x_high, y);
      if ( y < highFraction*ymax ) break;
    }

  //----------------------------------------------------------------
  //The following is the best configuration for picsec laser trigger
  //Fit from 2 samples left of the peak location to 9 samples left
  //of peak location. That gives time resolution of 13.0ps from trigger
  //signal to Photek signal (using gaus peak)
  //----------------------------------------------------------------
  if (trigger)
    {
      pulse->GetPoint(index_min-9, x_low, y);
      pulse->GetPoint(index_min-2, x_high, y);
      pulse->GetPoint(index_min, dummy, y);
    }
 
  TF1* flinear = new TF1("flinear","[0]*x+[1]", x_low, x_high );  
  pulse->Fit(flinear,"Q","", x_low, x_high );
  double slope = flinear->GetParameter(0);
  double b     = flinear->GetParameter(1);
  
  if ( makePlot )
    {
      //std::cout << "make plot" << std::endl;
      TCanvas* c = new TCanvas("canvas","canvas",800,400) ;
      pulse->GetXaxis()->SetLimits(x_low-10, x_high+10);
      //pulse->GetXaxis()->SetLimits(x_low, x_high);
      pulse->SetMarkerSize(0.5);
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
  return (float)(((0.9*ymax-b)/slope) - ((0.1*ymax-b)/slope));
};


void FitDirectHitPlusScintillationSignal(TGraphErrors* pulse, const float index_min, float* result, int event, TString fname, bool makePlot) {

  double x_low, x_high, xdummy, y, dummy;
  double ymax;
  pulse->GetPoint(index_min, x_low, ymax);
 
  TF1* flinear = new TF1("flinear","[0]*x+[1]", x_low, x_high );  
  pulse->Fit(flinear,"Q","", x_low, x_high );
  double slope = flinear->GetParameter(0);
  double b     = flinear->GetParameter(1);
  
  if ( makePlot )
    {
      //std::cout << "make plot" << std::endl;
      TCanvas* c = new TCanvas("canvas","canvas",800,400) ;
      pulse->GetXaxis()->SetLimits(x_low-10, x_high+10);
      //pulse->GetXaxis()->SetLimits(x_low, x_high);
      pulse->SetMarkerSize(0.5);
      pulse->SetMarkerStyle(20);
      pulse->Draw("AP");
      c->SaveAs(fname+"LinearFit.pdf");
      //delete c;
    }

  result[0] = 0; //direct hit amplitude
  result[1] = 0; //scintillation amplitude

  return;
}
