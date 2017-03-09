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

float find_CFD(TH1F * pulse, int first_bin, int last_bin, int minBin);
int FindMin( int n, int splineBinFactor, float *a);
int FindMax( int n, int splineBinFactor, float *a);
int FindRisingEdge( int n, int binMax, float *a);
int FindFirstPulsePeak( int n, float *a);
unsigned int CheckPulseQuality( int binMin, int binMax, float *a, float minPulse);
float ChannelIntegral(float *a, int peak);

float LED( TH1F * pulse, double threshold, int nsamples, int splineBinFactor );
float LinearFit_Baseline(TH1F * pulse, const int index_min, const int range);
float LinearFit_Intercept(TH1F * pulse, const float base, const int index_first, const int index_last);
float GausFit_MeanTime(TH1F * pulse, const int index_first, const int index_last);
void FitRisingEdge(TH1F* pulse, int nbinsL, int nbinsH, float &THM, float &risetime, float base, bool setbins);
void FitFullPulse(TH1F* pulse, float &par0, float &par1, float &par2);
TH1F* InterpolateWaveform(int nsamples, float* outputwaveform, float *inputwaveform, int splineBinFactor, std::string name);

const int Nsamples = 1024;

int main (int argc, char **argv) {

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
  if (argc >= 4) 
      includePulseshapeInOutput = bool(atoi(argv[3]));
  
  //const int splineBinFactor = 40;
  const int splineBinFactor = 1;

  int t_[Nsamples]; 

  for (int i = 0; i < Nsamples; i++) {
      t_[i] = i;
  }

  float Channel1VoltagesRaw_[Nsamples];
  float Channel2VoltagesRaw_[Nsamples];
  float Channel3VoltagesRaw_[Nsamples];
  float Channel4VoltagesRaw_[Nsamples];
  float Channel5VoltagesRaw_[Nsamples];
  float Channel6VoltagesRaw_[Nsamples];
  float Channel7VoltagesRaw_[Nsamples];
  float Channel8VoltagesRaw_[Nsamples];
  float Channel1Voltages_[Nsamples*splineBinFactor];
  float Channel2Voltages_[Nsamples*splineBinFactor];
  float Channel3Voltages_[Nsamples*splineBinFactor];
  float Channel4Voltages_[Nsamples*splineBinFactor];
  float Channel5Voltages_[Nsamples*splineBinFactor];
  float Channel6Voltages_[Nsamples*splineBinFactor];
  float Channel7Voltages_[Nsamples*splineBinFactor];
  float Channel8Voltages_[Nsamples*splineBinFactor];

  bool convert2Volts = true;

  TTree* t1 = (TTree*)f->Get("p");      // andriy's converter
  if (t1) 
    {
      t1->SetBranchAddress("c1",Channel1VoltagesRaw_);
      t1->SetBranchAddress("c2",Channel2VoltagesRaw_);
      t1->SetBranchAddress("c3",Channel3VoltagesRaw_);
      t1->SetBranchAddress("c4",Channel4VoltagesRaw_);
      t1->SetBranchAddress("c5",Channel5VoltagesRaw_);
      t1->SetBranchAddress("c6",Channel6VoltagesRaw_);
      t1->SetBranchAddress("c7",Channel7VoltagesRaw_);
      t1->SetBranchAddress("c8",Channel8VoltagesRaw_);
    }
  
  if (!t1) 
    {
      t1 = (TTree*)f->Get("T");   // artur's converter
      convert2Volts = false;

      t1->SetBranchAddress("c1",Channel1VoltagesRaw_);
      t1->SetBranchAddress("c2",Channel2VoltagesRaw_);
      t1->SetBranchAddress("c3",Channel3VoltagesRaw_);
      t1->SetBranchAddress("c4",Channel4VoltagesRaw_);
    }

  //create two histograms
  TH1F *CH1pulseRaw   = new TH1F("CH1pulseRaw","CH1pulseRaw",Nsamples,0,Nsamples);
  TH1F *CH2pulseRaw   = new TH1F("CH2pulseRaw","CH2pulseRaw",Nsamples,0,Nsamples);
  TH1F *CH3pulseRaw   = new TH1F("CH3pulseRaw","CH3pulseRaw",Nsamples,0,Nsamples);
  TH1F *CH4pulseRaw   = new TH1F("CH4pulseRaw","CH4pulseRaw",Nsamples,0,Nsamples);
  TH1F *CH5pulseRaw   = new TH1F("CH5pulseRaw","CH5pulseRaw",Nsamples,0,Nsamples);
  TH1F *CH6pulseRaw   = new TH1F("CH6pulseRaw","CH6pulseRaw",Nsamples,0,Nsamples);
  TH1F *CH7pulseRaw   = new TH1F("CH7pulseRaw","CH7pulseRaw",Nsamples,0,Nsamples);
  TH1F *CH8pulseRaw   = new TH1F("CH8pulseRaw","CH8pulseRaw",Nsamples,0,Nsamples);
  TH1F *CH1pulse = 0;
  TH1F *CH2pulse = 0;
  TH1F *CH3pulse = 0;
  TH1F *CH4pulse = 0;
  TH1F *CH5pulse = 0;
  TH1F *CH6pulse = 0;
  TH1F *CH7pulse = 0;
  TH1F *CH8pulse = 0;
  
  TH1F *CH1Amp     = new TH1F("CH1Amp","CH1Amp",40,-0.3,-0.6);
  TH1F *CH2Amp     = new TH1F("CH2Amp","CH2Amp",40,-0.3,-0.6);

  TH1F *GausPeak_CH12_dt   = new TH1F("GausPeak_CH12_dt","GausPeak_CH12_dt; t1-t2 [ns]; Events",4000,-4,4);
  TH1F *GausPeak_CH34_dt   = new TH1F("GausPeak_CH34_dt","GausPeak_CH34_dt; t3-t4 [ns]; Events",4000,-4,4);
  TH1F *GausPeak_TOF_CH13   = new TH1F("GausPeak_TOF_CH13","GausPeak_TOF_CH13; t1-t3 [ns]; Events",4000,-4,4);
  TH1F *GausPeak_TOF_CH14   = new TH1F("GausPeak_TOF_CH14","GausPeak_TOF_CH14; t1-t4 [ns]; Events",4000,-40,4);

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
  } else{
    fout = new TFile(argv[2],"recreate");
  }
  TTree* treeOut = new TTree("tree","tree");
  
  unsigned int eventNumber = 0;
  float ch1Time_gausfitroot = 0;
  float ch2Time_gausfitroot = 0;
  float ch3Time_gausfitroot = 0;
  float ch4Time_gausfitroot = 0;
  float ch5Time_gausfitroot = 0;
  float ch6Time_gausfitroot = 0;
  float ch7Time_gausfitroot = 0;
  float ch8Time_gausfitroot = 0;
  
  float ch1Amp = 0;
  float ch2Amp = 0;
  float ch3Amp = 0;
  float ch4Amp = 0;
  float ch5Amp = 0;
  float ch6Amp = 0;
  float ch7Amp = 0;
  float ch8Amp = 0;

  float ch1THM = 0;//Time at half the Maximum
  float ch2THM = 0;//Time at half the Maximum
  float ch3THM = 0;//Time at half the Maximum
  float ch4THM = 0;//Time at half the Maximum
  float ch5THM = 0;//Time at half the Maximum
  float ch6THM = 0;//Time at half the Maximum
  float ch7THM = 0;//Time at half the Maximum
  float ch8THM = 0;//Time at half the Maximum

  float ch1Risetime = 0;
  float ch2Risetime = 0;
  float ch3Risetime = 0;
  float ch4Risetime = 0;
  float ch5Risetime = 0;
  float ch6Risetime = 0;
  float ch7Risetime = 0;
  float ch8Risetime = 0;

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
  unsigned int ch5QualityBit = 0;
  unsigned int ch6QualityBit = 0;
  unsigned int ch7QualityBit = 0;
  unsigned int ch8QualityBit = 0;
  
  float ch1chisq = -1;
  float ch2chisq = -1;
  float ch3chisq = -1;
  float ch4chisq = -1;
  


  treeOut->Branch("event",&eventNumber,"event/i");
  treeOut->Branch("t1gausroot",&ch1Time_gausfitroot,"t1gausroot/F");
  treeOut->Branch("t2gausroot",&ch2Time_gausfitroot,"t2gausroot/F");
  treeOut->Branch("t3gausroot",&ch3Time_gausfitroot,"t3gausroot/F");
  treeOut->Branch("t4gausroot",&ch4Time_gausfitroot,"t4gausroot/F");
  treeOut->Branch("t5gausroot",&ch5Time_gausfitroot,"t5gausroot/F");
  treeOut->Branch("t6gausroot",&ch6Time_gausfitroot,"t6gausroot/F");
  treeOut->Branch("t7gausroot",&ch7Time_gausfitroot,"t7gausroot/F");
  treeOut->Branch("t8gausroot",&ch8Time_gausfitroot,"t8gausroot/F");
  
  treeOut->Branch("ch1Amp",&ch1Amp,"ch1Amp/F");
  treeOut->Branch("ch2Amp",&ch2Amp,"ch2Amp/F");
  treeOut->Branch("ch3Amp",&ch3Amp,"ch3Amp/F");
  treeOut->Branch("ch4Amp",&ch4Amp,"ch4Amp/F");
  treeOut->Branch("ch5Amp",&ch5Amp,"ch5Amp/F");
  treeOut->Branch("ch6Amp",&ch6Amp,"ch6Amp/F");
  treeOut->Branch("ch7Amp",&ch7Amp,"ch7Amp/F");
  treeOut->Branch("ch8Amp",&ch8Amp,"ch8Amp/F");

  treeOut->Branch("ch1THM",&ch1THM,"ch1THM/F");
  treeOut->Branch("ch2THM",&ch2THM,"ch2THM/F");
  treeOut->Branch("ch3THM",&ch3THM,"ch3THM/F");
  treeOut->Branch("ch4THM",&ch4THM,"ch4THM/F");
  treeOut->Branch("ch5THM",&ch5THM,"ch5THM/F");
  treeOut->Branch("ch6THM",&ch6THM,"ch6THM/F");
  treeOut->Branch("ch7THM",&ch7THM,"ch7THM/F");
  treeOut->Branch("ch8THM",&ch8THM,"ch8THM/F");
  treeOut->Branch("ch1Risetime",&ch1Risetime,"ch1Risetime/F");
  treeOut->Branch("ch2Risetime",&ch2Risetime,"ch2Risetime/F");
  treeOut->Branch("ch3Risetime",&ch3Risetime,"ch3Risetime/F");
  treeOut->Branch("ch4Risetime",&ch4Risetime,"ch4Risetime/F");
  treeOut->Branch("ch5Risetime",&ch5Risetime,"ch5Risetime/F");
  treeOut->Branch("ch6Risetime",&ch6Risetime,"ch6Risetime/F");
  treeOut->Branch("ch7Risetime",&ch7Risetime,"ch7Risetime/F");
  treeOut->Branch("ch8Risetime",&ch8Risetime,"ch8Risetime/F");

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
  treeOut->Branch("ch5QualityBit",&ch5QualityBit,"ch5QualityBit/i");
  treeOut->Branch("ch6QualityBit",&ch6QualityBit,"ch6QualityBit/i");
  treeOut->Branch("ch7QualityBit",&ch7QualityBit,"ch7QualityBit/i");
  treeOut->Branch("ch8QualityBit",&ch8QualityBit,"ch8QualityBit/i");
  
  treeOut->Branch("ch1Int",&ch1Int,"ch1Int/F");
  treeOut->Branch("ch2Int",&ch2Int,"ch2Int/F");
  treeOut->Branch("ch3Int",&ch3Int,"ch3Int/F");
  treeOut->Branch("ch4Int",&ch4Int,"ch4Int/F");
  treeOut->Branch("ch5Int",&ch5Int,"ch5Int/F");
  treeOut->Branch("ch6Int",&ch6Int,"ch6Int/F");
  treeOut->Branch("ch7Int",&ch7Int,"ch7Int/F");
  treeOut->Branch("ch8Int",&ch8Int,"ch8Int/F");
  
  treeOut->Branch("ch1chisq",&ch1chisq,"ch1chisq/F");
  treeOut->Branch("ch2chisq",&ch2chisq,"ch2chisq/F");
  treeOut->Branch("ch3chisq",&ch3chisq,"ch3chisq/F");
  treeOut->Branch("ch4chisq",&ch4chisq,"ch4chisq/F");
  

  if (includePulseshapeInOutput) {

    treeOut->Branch("c1",Channel1VoltagesRaw_,"c1[1024]/F");
    treeOut->Branch("c2",Channel2VoltagesRaw_,"c2[1024]/F");
    treeOut->Branch("c3",Channel3VoltagesRaw_,"c3[1024]/F");
    treeOut->Branch("c4",Channel4VoltagesRaw_,"c4[1024]/F");
    treeOut->Branch("c5",Channel5VoltagesRaw_,"c5[1024]/F");
    treeOut->Branch("c6",Channel6VoltagesRaw_,"c6[1024]/F");
    treeOut->Branch("c7",Channel7VoltagesRaw_,"c7[1024]/F");
    treeOut->Branch("c8",Channel8VoltagesRaw_,"c8[1024]/F");
    treeOut->Branch("t",t_,"t[1024]/I");

  }

  //read all entries and fill the histograms
  Long64_t nentries = t1->GetEntries();
  
  //DEBUGGING: set num entries to change (comment out when done)
  //nentries = 2000;

  for (Long64_t iEntry=0;iEntry<nentries;iEntry++)   
    {
      if(iEntry%100==0) std::cout<<"Processing Event: "<<iEntry<<" out of: "<<nentries<<std::endl;
      
      t1->GetEntry(iEntry);
      eventNumber = iEntry+1;

      //////////////////////
      // convert to Volts
      //////////////////////
      if(convert2Volts){
	for (int ii=0;ii<Nsamples;ii++){
	  Channel1VoltagesRaw_[ii] = 0.001*Channel1VoltagesRaw_[ii];
	  Channel2VoltagesRaw_[ii] = 0.001*Channel2VoltagesRaw_[ii];
	  Channel3VoltagesRaw_[ii] = 0.001*Channel3VoltagesRaw_[ii];
	  Channel4VoltagesRaw_[ii] = 0.001*Channel4VoltagesRaw_[ii];
	  Channel5VoltagesRaw_[ii] = 0.001*Channel5VoltagesRaw_[ii];
	  Channel6VoltagesRaw_[ii] = 0.001*Channel6VoltagesRaw_[ii];
	  Channel7VoltagesRaw_[ii] = 0.001*Channel7VoltagesRaw_[ii];
	  Channel8VoltagesRaw_[ii] = 0.001*Channel8VoltagesRaw_[ii];
	}
      }
      //////////////////////
      // end convert to Volts
      //////////////////////


      // Set raw pulses
      for (int ii=0;ii<Nsamples;ii++) {
	CH1pulseRaw->SetBinContent(ii+1,Channel1VoltagesRaw_[ii]);
	CH2pulseRaw->SetBinContent(ii+1,Channel2VoltagesRaw_[ii]);
	CH3pulseRaw->SetBinContent(ii+1,Channel3VoltagesRaw_[ii]);
	CH4pulseRaw->SetBinContent(ii+1,Channel4VoltagesRaw_[ii]);
	CH5pulseRaw->SetBinContent(ii+1,Channel5VoltagesRaw_[ii]);
	CH6pulseRaw->SetBinContent(ii+1,Channel6VoltagesRaw_[ii]);
	CH7pulseRaw->SetBinContent(ii+1,Channel7VoltagesRaw_[ii]);
	CH8pulseRaw->SetBinContent(ii+1,Channel8VoltagesRaw_[ii]);	  
	CH1pulseRaw->SetBinError(ii+1, 0.001);
	CH2pulseRaw->SetBinError(ii+1, 0.001);
	CH3pulseRaw->SetBinError(ii+1, 0.001);
	CH4pulseRaw->SetBinError(ii+1, 0.001);
	CH5pulseRaw->SetBinError(ii+1, 0.001);
	CH6pulseRaw->SetBinError(ii+1, 0.001);
	CH7pulseRaw->SetBinError(ii+1, 0.001);
	CH8pulseRaw->SetBinError(ii+1, 0.001);
      }

      //do spline to add more points to the waveform
      CH1pulse = InterpolateWaveform(Nsamples, Channel1Voltages_, Channel1VoltagesRaw_, splineBinFactor, "CH1pulse");
      CH2pulse = InterpolateWaveform(Nsamples, Channel2Voltages_, Channel2VoltagesRaw_, splineBinFactor, "CH2pulse");
      CH3pulse = InterpolateWaveform(Nsamples, Channel3Voltages_, Channel3VoltagesRaw_, splineBinFactor, "CH3pulse");
      CH4pulse = InterpolateWaveform(Nsamples, Channel4Voltages_, Channel4VoltagesRaw_, splineBinFactor, "CH4pulse");
      CH5pulse = InterpolateWaveform(Nsamples, Channel5Voltages_, Channel5VoltagesRaw_, splineBinFactor, "CH5pulse");
      CH6pulse = InterpolateWaveform(Nsamples, Channel6Voltages_, Channel6VoltagesRaw_, splineBinFactor, "CH6pulse");
      CH7pulse = InterpolateWaveform(Nsamples, Channel7Voltages_, Channel7VoltagesRaw_, splineBinFactor, "CH7pulse");
      CH8pulse = InterpolateWaveform(Nsamples, Channel8Voltages_, Channel8VoltagesRaw_, splineBinFactor, "CH8pulse");
      

      // Find Min of the Channel data (Voltage)
      int index_min1 = FindMin (Nsamples, splineBinFactor, Channel1Voltages_);	// return index of the min
      int index_min2 = FindMin (Nsamples, splineBinFactor, Channel2Voltages_);	// return index of the min
      int index_min3 = FindMin (Nsamples, splineBinFactor, Channel3Voltages_);	// return index of the min
      int index_min4 = FindMin (Nsamples, splineBinFactor, Channel4Voltages_);	// return index of the min
      int index_min5 = FindMin (Nsamples, splineBinFactor, Channel5Voltages_);	// return index of the min
      int index_min6 = FindMin (Nsamples, splineBinFactor, Channel6Voltages_);	// return index of the min
      int index_min7 = FindMin (Nsamples, splineBinFactor, Channel7Voltages_);	// return index of the min
      int index_min8 = FindMin (Nsamples, splineBinFactor, Channel8Voltages_);	// return index of the min

      
      assert(CH1pulse);
      assert(CH2pulse);
      assert(CH3pulse);
      assert(CH4pulse);
      assert(CH5pulse);
      assert(CH6pulse);
      assert(CH7pulse);
      assert(CH8pulse);

      // Find Max of the Channel data (Voltage)
      int index_max1 = FindMax (Nsamples, splineBinFactor, Channel1Voltages_);	// return index of the max
      int index_max2 = FindMax (Nsamples, splineBinFactor, Channel2Voltages_);	// return index of the max
      int index_max3 = FindMax (Nsamples, splineBinFactor, Channel3Voltages_);	// return index of the max
      int index_max4 = FindMax (Nsamples, splineBinFactor, Channel4Voltages_);	// return index of the max
      int index_max5 = FindMax (Nsamples, splineBinFactor, Channel5Voltages_);	// return index of the max
      int index_max6 = FindMax (Nsamples, splineBinFactor, Channel6Voltages_);	// return index of the max
      int index_max7 = FindMax (Nsamples, splineBinFactor, Channel7Voltages_);	// return index of the max
      int index_max8 = FindMax (Nsamples, splineBinFactor, Channel8Voltages_);	// return index of the max
     
      
      //////////////////
      // Find the rising edge on CH1
      int fbin1 = FindRisingEdge(Nsamples, index_min1, Channel1Voltages_);
      int fbin2 = FindRisingEdge(Nsamples, index_min2, Channel2Voltages_);
      int fbin3 = FindRisingEdge(Nsamples, index_min3, Channel3Voltages_);
      int fbin4 = FindRisingEdge(Nsamples, index_min4, Channel4Voltages_);
      int fbin5 = FindRisingEdge(Nsamples, index_min5, Channel5Voltages_);
      int fbin6 = FindRisingEdge(Nsamples, index_min6, Channel6Voltages_);
      int fbin7 = FindRisingEdge(Nsamples, index_min7, Channel7Voltages_);
      int fbin8 = FindRisingEdge(Nsamples, index_min8, Channel8Voltages_);
      
      // Find the quality of the pulse
      ch1QualityBit = CheckPulseQuality ( index_min1, index_max1, Channel1VoltagesRaw_, 0.02);
      ch2QualityBit = CheckPulseQuality ( index_min2, index_max2, Channel2VoltagesRaw_, 0.02);
      ch3QualityBit = CheckPulseQuality ( index_min3, index_max3, Channel3VoltagesRaw_, 0.02);
      ch4QualityBit = CheckPulseQuality ( index_min4, index_max4, Channel4VoltagesRaw_, 0.02);
      ch5QualityBit = CheckPulseQuality ( index_min5, index_max5, Channel5VoltagesRaw_, 0.02);
      ch6QualityBit = CheckPulseQuality ( index_min6, index_max6, Channel6VoltagesRaw_, 0.02);
      ch7QualityBit = CheckPulseQuality ( index_min7, index_max7, Channel7VoltagesRaw_, 0.02);
      ch8QualityBit = CheckPulseQuality ( index_min8, index_max8, Channel8VoltagesRaw_, 0.02);
     
      // For the first version of Photonis data the pulse has many peaks --> find the first one
      int index_firstPulse1 = FindFirstPulsePeak(Nsamples, Channel1Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse2 = FindFirstPulsePeak(Nsamples, Channel2Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse3 = FindFirstPulsePeak(Nsamples, Channel3Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse4 = FindFirstPulsePeak(Nsamples, Channel4Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse5 = FindFirstPulsePeak(Nsamples, Channel5Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse6 = FindFirstPulsePeak(Nsamples, Channel6Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse7 = FindFirstPulsePeak(Nsamples, Channel7Voltages_); // this is useful ONLY for Photonis MCP
      int index_firstPulse8 = FindFirstPulsePeak(Nsamples, Channel8Voltages_); // this is useful ONLY for Photonis MCP
       
      //////////////////////////////////////////////
      ///////// Done with setup, start the fits
      //////////////////////////////////////////////

      // invert the pulses for easy fitting
      CH1pulse->Scale(-1.);
      CH2pulse->Scale(-1.);
      CH3pulse->Scale(-1.);
      CH4pulse->Scale(-1.);
      CH5pulse->Scale(-1.);
      CH6pulse->Scale(-1.);
      CH7pulse->Scale(-1.);
      CH8pulse->Scale(-1.);
      CH1pulseRaw->Scale(-1.);
      CH2pulseRaw->Scale(-1.);
      CH3pulseRaw->Scale(-1.);
      CH4pulseRaw->Scale(-1.);
      CH5pulseRaw->Scale(-1.);
      CH6pulseRaw->Scale(-1.);
      CH7pulseRaw->Scale(-1.);
      CH8pulseRaw->Scale(-1.);
        
      // fit the baseline
      float base1 = LinearFit_Baseline( CH1pulse, index_min1, 10 );
      float base2 = LinearFit_Baseline( CH2pulse, index_min2, 10 );
      float base3 = LinearFit_Baseline( CH3pulse, index_min3, 10 );
      float base4 = LinearFit_Baseline( CH4pulse, index_min4, 10 );
      float base5 = LinearFit_Baseline( CH5pulse, index_min5, 10 );
      float base6 = LinearFit_Baseline( CH6pulse, index_min6, 10 );
      float base7 = LinearFit_Baseline( CH7pulse, index_min7, 10 );
      float base8 = LinearFit_Baseline( CH8pulse, index_min8, 10 );
<<<<<<< HEAD

      /*base1 = 0;
      base2 = 0;
=======
    /*  base1 = 0;
      base2 = 0;        
>>>>>>> 7ace26a8e87e9b7269fa37c6080a6f279fa5cdfe
      base3 = 0;
      base4 = 0;
      base5 = 0;
      base6 = 0;
      base7 = 0;
<<<<<<< HEAD
      base8 = 0;
      */
=======
      base8 = 0; */

>>>>>>> 7ace26a8e87e9b7269fa37c6080a6f279fa5cdfe
      //std::cout << "baseline: " << base1 << " " << base2 << "\n";

      /////////////////////////
      // Find the amplitudes
      /////////////////////////
      ch1Amp = -1 * Channel1Voltages_[index_min1] - base1;
      ch2Amp = -1 * Channel2Voltages_[index_min2] - base2;
      ch3Amp = -1 * Channel3Voltages_[index_min3] - base3;
      ch4Amp = -1 * Channel4Voltages_[index_min4] - base4;
      
      //////////////////////////////////////////////////////////////////////
      
      // DEBUGGING CODE 
      /*
      std::cout << "==========" << iEntry << "==========" << std::endl;
      std::cout << index_min4 << "," << index_max4 << std::endl;
      
      float max = 0.0;
      int maxi = 0;
      for (int i =  0; i < sizeof(Channel4Voltages_); i++) {
        if (Channel4Voltages_[i] < max) {
            max = Channel4Voltages_[i];
            maxi = i;}
      }
      
      float max2 = 0.0;
      int maxi2 = 0;
      for (int i =  0; i < sizeof(Channel4VoltagesRaw_); i++) {
        if (Channel4VoltagesRaw_[i] < max) {
            max2 = Channel4VoltagesRaw_[i];
            maxi2 = i;}
      }
      
      std::cout << "Raw: " << maxi2 << "," << max2 << std::endl;
      std::cout << "Outputed: " << maxi << "," << max << "\n" << std::endl;
      */
      ////////////////////////////////////////////////////////////////////////
      
      ch5Amp = -1 * Channel5Voltages_[index_min5] - base5;
      ch6Amp = -1 * Channel6Voltages_[index_min6] - base6;
      ch7Amp = -1 * Channel7Voltages_[index_min7] - base7;
      ch8Amp = -1 * Channel8Voltages_[index_min8] - base8;
      
      ch1Int = -1 * ChannelIntegral(Channel1Voltages_, index_min1) - 7 * base1;
      ch2Int = -1 * ChannelIntegral(Channel2Voltages_, index_min2) - 7 * base2;
      ch3Int = -1 * ChannelIntegral(Channel3Voltages_, index_min3) - 7 * base3;
      ch4Int = -1 * ChannelIntegral(Channel4Voltages_, index_min4) - 7 * base4;
      ch5Int = -1 * ChannelIntegral(Channel5Voltages_, index_min5) - 7 * base5;
      ch6Int = -1 * ChannelIntegral(Channel6Voltages_, index_min6) - 7 * base6;
      ch7Int = -1 * ChannelIntegral(Channel7Voltages_, index_min7) - 7 * base7;
      ch8Int = -1 * ChannelIntegral(Channel8Voltages_, index_min8) - 7 * base8;


      ///////////////////
      // Gaussian fit
      // fit the gaussian peak
      ///////////////////
      // float timepeak1 =  GausFit_MeanTime(CH1pulse, index_min1 - 6*splineBinFactor, index_min1+6*splineBinFactor);
      // float timepeak2 =  GausFit_MeanTime(CH2pulse, index_min2 - 6*splineBinFactor, index_min2+6*splineBinFactor);
      // float timepeak3 =  GausFit_MeanTime(CH3pulse, index_min3 - 6*splineBinFactor, index_min3+6*splineBinFactor);
      // float timepeak4 =  GausFit_MeanTime(CH4pulse, index_min4 - 6*splineBinFactor, index_min4+6*splineBinFactor);
      // float timepeak5 =  GausFit_MeanTime(CH5pulse, index_min5 - 6*splineBinFactor, index_min5+6*splineBinFactor);
      // float timepeak6 =  GausFit_MeanTime(CH6pulse, index_min6 - 6*splineBinFactor, index_min6+6*splineBinFactor);
      // float timepeak7 =  GausFit_MeanTime(CH7pulse, index_min7 - 6*splineBinFactor, index_min7+6*splineBinFactor);
      // float timepeak8 =  GausFit_MeanTime(CH8pulse, index_min8 - 6*splineBinFactor, index_min8+6*splineBinFactor);
      float timepeak1 =  GausFit_MeanTime(CH1pulse, index_min1 - 4*splineBinFactor, index_min1+4*splineBinFactor);
      float timepeak2 =  GausFit_MeanTime(CH2pulse, index_min2 - 3*splineBinFactor, index_min2+4*splineBinFactor);
      float timepeak3 =  GausFit_MeanTime(CH3pulse, index_min3 - 9*splineBinFactor, index_min3+9*splineBinFactor);
      float timepeak4 =  GausFit_MeanTime(CH4pulse, index_min4 - 8*splineBinFactor, index_min4+7*splineBinFactor);
      float timepeak5 =  GausFit_MeanTime(CH5pulse, index_min5 - 4*splineBinFactor, index_min5+4*splineBinFactor);
      float timepeak6 =  GausFit_MeanTime(CH6pulse, index_min6 - 4*splineBinFactor, index_min6+4*splineBinFactor);
      float timepeak7 =  GausFit_MeanTime(CH7pulse, index_min7 - 4*splineBinFactor, index_min7+4*splineBinFactor);
      float timepeak8 =  GausFit_MeanTime(CH8pulse, index_min8 - 4*splineBinFactor, index_min8+4*splineBinFactor);
                
      ch1Time_gausfitroot = timepeak1*0.2/splineBinFactor;
      ch2Time_gausfitroot = timepeak2*0.2/splineBinFactor;
      ch3Time_gausfitroot = timepeak3*0.2/splineBinFactor;
      ch4Time_gausfitroot = timepeak4*0.2/splineBinFactor;
      ch5Time_gausfitroot = timepeak5*0.2/splineBinFactor;
      ch6Time_gausfitroot = timepeak6*0.2/splineBinFactor;
      ch7Time_gausfitroot = timepeak7*0.2/splineBinFactor;
      ch8Time_gausfitroot = timepeak8*0.2/splineBinFactor;
            
      ///////////////////
      // CFD
      ///////////////////
      // float timepeak1 =  LED(CH1pulse, 0.50 * ch1Amp, Nsamples ,splineBinFactor);
      // float timepeak2 =  LED(CH2pulse, 0.50 * ch2Amp, Nsamples ,splineBinFactor);
      // float timepeak3 =  LED(CH3pulse, 0.50 * ch3Amp, Nsamples ,splineBinFactor);
      // float timepeak4 =  LED(CH4pulse, 0.50 * ch4Amp, Nsamples ,splineBinFactor);
      // float timepeak5 =  LED(CH5pulse, 0.50 * ch5Amp, Nsamples ,splineBinFactor);
      // float timepeak6 =  LED(CH6pulse, 0.50 * ch6Amp, Nsamples ,splineBinFactor);
      // float timepeak7 =  LED(CH7pulse, 0.50 * ch7Amp, Nsamples ,splineBinFactor);
      // float timepeak8 =  LED(CH8pulse, 0.50 * ch8Amp, Nsamples ,splineBinFactor);
      
      
      // ch1Time_gausfitroot = timepeak1*0.2/splineBinFactor;
      // ch2Time_gausfitroot = timepeak2*0.2/splineBinFactor;
      // ch3Time_gausfitroot = timepeak3*0.2/splineBinFactor;
      // ch4Time_gausfitroot = timepeak4*0.2/splineBinFactor;
      // ch5Time_gausfitroot = timepeak5*0.2/splineBinFactor;
      // ch6Time_gausfitroot = timepeak6*0.2/splineBinFactor;
      // ch7Time_gausfitroot = timepeak7*0.2/splineBinFactor;
      // ch8Time_gausfitroot = timepeak8*0.2/splineBinFactor;

      //FitFullPulse(CH2pulse, ch1_TFF, ch1_AFF);//ch1 Time Full Fit
      //FitFullPulse(CH3pulse, ch2_TFF, ch2_AFF);//ch2 Time Full Fit
      //FitFullPulse(CH4pulse, ch4_TFF, ch4_AFF, ch4_TFF_v2);//ch4 Time Full Fit
      
      //std::cout << "ch4_AFF_v2: " << ch4_TFF_v2 << std::endl;
      //std::cout << "ch2_AFF: " << ch2_AFF << std::endl;
      
      //Fit Rising Edge
<<<<<<< HEAD
      //FitRisingEdge(CH1pulse, -1, 0, ch1THM, ch1Risetime, base1);
      FitRisingEdge(CH2pulse, -1, 0, ch2THM, ch2Risetime, base2);
      //      FitRisingEdge(CH3pulse, -1, 0, ch3THM, ch3Risetime, 415, 435);
      FitRisingEdge(CH4pulse, -1, 0, ch4THM, ch4Risetime, base4);
            
=======
      //std::cout << iEntry << std::endl;
      FitRisingEdge(CH1pulse, -1, 0, ch1THM, ch1Risetime, base1, false);
      FitRisingEdge(CH2pulse, -1, 0, ch2THM, ch2Risetime, base2, true);
      FitRisingEdge(CH3pulse, -1, 0, ch3THM, ch3Risetime, base3, false);
      FitRisingEdge(CH4pulse, -1, 0, ch4THM, ch4Risetime, base4, false);
      
>>>>>>> 7ace26a8e87e9b7269fa37c6080a6f279fa5cdfe
      if(ch1Amp < 0.05 || ch2Amp < 0.05){
        //continue;
      }

//       if(ch3Amp > 0.05 && iEntry > 300) {
//   	break;
//       }

      //if ( iEntry >= 2 ) break;
      //Fill the tree
      treeOut->Fill();

     
      if (iEntry < nentries-1) {      
	delete CH1pulse;
	delete CH2pulse;
	delete CH3pulse;
	delete CH4pulse;
	delete CH5pulse;
	delete CH6pulse;
	delete CH7pulse;
	delete CH8pulse;
      }


    }
  
  CH1pulse->Write();
  CH2pulse->Write();
  CH3pulse->Write();
  CH4pulse->Write();
  CH5pulse->Write();
  CH6pulse->Write();
  CH7pulse->Write();
  CH8pulse->Write();
  CH1pulseRaw->Write();
  CH2pulseRaw->Write();
  CH3pulseRaw->Write();
  CH4pulseRaw->Write();
  CH5pulseRaw->Write();
  CH6pulseRaw->Write();
  CH7pulseRaw->Write();
  CH8pulseRaw->Write();
 
  GausPeak_CH12_dt->Write();
  GausPeak_CH34_dt->Write();
  GausPeak_TOF_CH13->Write();
  GausPeak_TOF_CH14->Write();
  
  CH1Amp->Write();
  CH2Amp->Write();
  //CH3Amp->Write();
  //CH4Amp->Write();
  
  treeOut->Write();

  fout->Close();
}

////////////////////////////////////////////////////
// Remove off-base (e.g. y-shifted) pulses from Raw Voltage
////////////////////////////////////////////////////

////////////////////////////////////////////
// Do Spline to interpolate
////////////////////////////////////////////
TH1F* InterpolateWaveform(int nsamples, float *outputwaveform, float *inputwaveform, int splineBinFactor, std::string name) {
  ROOT::Math::Interpolator cspline( nsamples, ROOT::Math::Interpolation::kCSPLINE);  

  TH1F *pulse = new TH1F(name.c_str(),name.c_str(),Nsamples * splineBinFactor,0,Nsamples*splineBinFactor);
  
  // IF ACTUALLY HAVE TO DO A SPLINE
  
  if (splineBinFactor != 1) {
    
    double bins[nsamples];
    double w[nsamples];
    for (int i=0; i<nsamples; i++) { 
      bins[i] = i+1;
      w[i] = inputwaveform[i];
      //std::cout << bins[i] << " " << w[i] << "\n";
    }

    cspline.SetData( nsamples, bins, w);


    //do spline
    for (int i=0; i < nsamples*splineBinFactor; i++) {
      if (i > splineBinFactor) {
	    outputwaveform[i] =  cspline.Eval( pulse->GetXaxis()->GetBinCenter(i+1) / splineBinFactor );
      } else {
	outputwaveform[i] = 0;
      }
    }
  
    //do low pass filter
    const int n_filter = 100;
    for( int i = 0+n_filter; i < nsamples*splineBinFactor-n_filter; i++){
      double temp = 0;
      double count = 0;
      for( int j = i-n_filter; j < i+n_filter; j++){
	temp += outputwaveform[j];
	count += 1.0;
      }
    
      outputwaveform[i] = temp/count;    
    }

    for (int i=1; i<=pulse->GetXaxis()->GetNbins();i++) {
      pulse->SetBinContent(i, outputwaveform[i-1]);
      pulse->SetBinError(i,0.001);    
    }
    
    // IF SPLINE FACTOR IS  JUST 1 (EG NO SPLINING)
    
  } else {

    for (int ii=0;ii<nsamples;ii++) {	
      outputwaveform[ii] = inputwaveform[ii];
    }
    
    // GET RID OF THIS LOW PASS NONSENSE
    
    /*
    //do low pass filter
    const int n_filter = 3;
    for( int i = 0+n_filter; i < nsamples-n_filter; i++){
      double temp = 0;
      double count = 0;
      for( int j = i-n_filter; j < i+n_filter; j++){
    	temp += outputwaveform[j];
    	count += 1.0;
      }    
      outputwaveform[i] = temp/count;    
    }
    */
    
    for (int ii=0;ii<nsamples;ii++) {	
      pulse->SetBinContent(ii+1,outputwaveform[ii]);
      pulse->SetBinError(ii+1, 0.001);
    }
  }

  return pulse;
}




////////////////////////////////////////////
// find minimum of the pulse
// aa added protection against pulses with single high bin
////////////////////////////////////////////
int FindMin( int n, int splineBinFactor, float *a) {
  
  if (n <= 0 || !a) return -1;
  float xmin = a[5];
  int loc = 0;
  for  (int i = 5*splineBinFactor; i < (n-5)*splineBinFactor; i++) {
    if (xmin > a[i] && a[i+1] < 0.5*a[i])  {
      xmin = a[i];
      loc = i;
    }
  }
  
  return loc;
}

////////////////////////////////////////////
// find maximum of the pulse 
////////////////////////////////////////////
int FindMax( int n, int splineBinFactor, float *a) {

  if (n <= 0 || !a) return -1;
  float xmax = a[0];
  int loc = 0;
  for  (int i = 0; i < n*splineBinFactor; i++) {
    if (xmax < a[i]) {
      xmax = a[i];
      loc = i;
    }
  }
  return loc;

}

////////////////////////////////////////////
// find rising edge of the pulse
////////////////////////////////////////////
int FindRisingEdge( int n, int binMax, float *a) {

  if (n <= 0 || !a) return -1;
  float xmin = a[0];
  int loc = -99;
  for (int i = binMax-10; i <binMax; i++)
  { // sometimes there is noise, check that it is rising in three bins
    if ( a[i] < -0.01 && a[i+1] < a[i] && a[i+2] < a[i+1] )
      // if ( Channel1Voltages_[i+2] < 0.3*Channel1Voltages_[index_min1])
    {
      loc = i; 
      break;
    }
  }  
  return loc;  
}

////////////////////////////////////////////
// for photonis: find the first pulse 
////////////////////////////////////////////
int FindFirstPulsePeak( int n, float *a) {

  if (n <= 0 || !a) return -1;

  int loc = 0;
  for  (int i = 20; i < 1000; i++) {
    if ( a[i] < -0.015
         && a[i] <= a[i-1] && a[i] <= a[i+1]
         && a[i-1] <= a[i-2]
         && a[i-2] <= a[i-3]
      ) {

      loc = i;
      break;
    }
  }
  return loc;

}

////////////////////////////////////////////
// assign pulse quality bit
////////////////////////////////////////////
unsigned int CheckPulseQuality( int binMin, int binMax, float *a, float minPulse) {
  unsigned int answer = 0;

  //*******************************************
  //Check if there is real pulse in the event
  //*******************************************
  if (!(a[binMin] < -minPulse)) {
    answer |= kNoPulse;
    // std::cout << "No Pulse: " << answer << "\n";
  }

  //*******************************************
  //Check pulses with large opposite amplitude 
  //*******************************************
  if (a[binMax] > 0.1) {
    answer |= kLargeNegativeAmplitude;
    //cout << "kLargeNegativeAmplitude: " << answer << "\n";
  }

  //*******************************************
  //Check pulses with saturation
  //*******************************************
  if ( (a[binMin] == a[binMin+1] || a[binMin] == a[binMin-1])
    ) {
    answer |= kSaturated;
    //cout << "kSaturated: " << answer << "\n";
  }


  //************************************************************************************
  //check that no points near the peak has negative polarity
  //************************************************************************************
  for (int ii=binMin-3;ii<binMin+4;ii++)
  {
    //cout << "bin " << ii << " : " << a[ii] << "\n";
    if (a[ii] > 0) {
      answer |= kNegativePolarity;
      //cout << "Negative Polarity\n";
    }
  }


  //************************************************************************************
  //check that no points near the peak has sudden jump
  //************************************************************************************
  bool hasSuddenJump = false;
  if (!(a[binMin] < a[binMin+1] && a[binMin+1] < a[binMin+2] && a[binMin+2] < a[binMin+3]  && a[binMin+3] < a[binMin+4])) hasSuddenJump = true;
  if (!(a[binMin] < a[binMin-1] && a[binMin-1] < a[binMin-2] && a[binMin-2] < a[binMin-3]  && a[binMin-3] < a[binMin-4])) hasSuddenJump = true;

  if (hasSuddenJump == true) {
//     cout << "Sudden Jump\n";
//     cout << a[binMin] << " " 
//          <<  a[binMin+1] << " "
//          <<  a[binMin+2] << " "
//          <<  a[binMin+3] << " "
//            <<  a[binMin+4] << " "
//          << "\n";
//     cout << a[binMin] << " " 
//          <<  a[binMin-1] << " "
//          <<  a[binMin-2] << " "
//          <<  a[binMin-3] << " "
//          <<  a[binMin-4] << " "
//          << "\n";
    answer |= kSuddenJump;
    //cout << "kSuddenJump: " << answer << "\n";
  }
  
  //************************************************************************************
  //check that second point away from peak is at least 10% lower - prevents strange pulses where there is flattish (but not completely flat) top
  //************************************************************************************
  bool hasFlatTop = false;
  if (!(fabs(a[binMin] - a[binMin+2])/fabs(a[binMin]) > 0.05)) hasFlatTop = true;
  if (!(fabs(a[binMin] - a[binMin-2])/fabs(a[binMin]) > 0.05)) hasFlatTop = true;

  if (hasFlatTop == true) {
    // std::cout << "flat pulse\n";
    // std::cout << a[binMin] << " " <<binMin<< " " 
    //      <<  a[binMin+2] << " "
    //      <<  fabs(a[binMin] - a[binMin+2])/fabs(a[binMin]) << " "
    //      << "\n";
    // std::cout << a[binMin] << " " 
    //      <<  a[binMin-2] << " "
    //      <<  fabs(a[binMin] - a[binMin-2])/fabs(a[binMin]) << " "

    //      << "\n";
    answer |= kFlatTop;
    //cout << "kFlatTop: " << answer << "\n";
  }

  //************************************************************************************
  //check for presence of 2nd pulse
  //************************************************************************************
  bool hasSecondPulse = false;
  int secondPulseIndex = -99;
  float secondPulseHeight = 0;
  for  (int i = 20; i < 1000; i++) {
    if (secondPulseHeight > a[i] 
        && fabs(a[i] - a[i-1])/fabs(a[i]) < 0.5
        && fabs(a[i] - a[i+1])/fabs(a[i]) < 0.5
        && abs(binMin - i) > 20
        && fabs(fabs(a[i]) - fabs(a[binMin]))/fabs(a[binMin]) < 0.15
      ) {
      secondPulseHeight = a[i];
      secondPulseIndex = i;
      hasSecondPulse = true;
    }
  }

//     cout << "Second Pulse\n";
//     cout << secondPulseIndex << " : " << secondPulseHeight << "\n";
  if (hasSecondPulse) {
    // std::cout << "First Pulse\n";
    // std::cout << binMin << " : " << a[binMin] << "\n";
    // std::cout << "Second Pulse\n";
    // std::cout << secondPulseIndex << " : " << secondPulseHeight << "\n";
    answer |= kSecondPulse;
  }
  
  //cout << "final answer : " << answer << "\n";
  return answer; 
}

///////////////////////
// find the baseline //
///////////////////////

float LinearFit_Baseline(TH1F * pulse, const int index_min, const int range)
{
  TF1 *fBaseline = new TF1("fBaseline","pol0",10, index_min-range);
  pulse->Fit("fBaseline","Q","", 10, index_min-range);
  // TF1 *fBaseline = new TF1("fBaseline","pol0",5, 100);
  // pulse->Fit("fBaseline","Q","", 5, 100);
  float base = fBaseline->GetParameter(0);
  delete fBaseline;
  
  return base;
}

/////////////////////////////////////////////////////////////
// find the intercept of the linear fit on the rising edge //
/////////////////////////////////////////////////////////////

float LinearFit_Intercept(TH1F* pulse, const float base, const int index_first, const int index_last)
{
  TF1* fRise = new TF1("fRise","pol1", index_first, index_last);
  pulse->Fit("fRise","Q","", index_first, index_last);
  float timeIntercept = (base - fRise->GetParameter(0))/fRise->GetParameter(1);
  delete fRise;

  return timeIntercept;
}

//////////////////////////////////////
// find the mean time from gaus fit //
//////////////////////////////////////


float GausFit_MeanTime(TH1F* pulse, const int index_first, const int index_last)
{
  TF1* fpeak = new TF1("fpeak","gaus", index_first, index_last);
  pulse->Fit("fpeak","Q","", index_first, index_last);
  float timepeak = fpeak->GetParameter(1);
  delete fpeak;
  
  return timepeak;
}

float ChannelIntegral(float *a, int peak) 
{
  float integral = 0.;

  //for sharp gaussian peak type pulses
  integral  = a[peak - 3] + a[peak - 2] + a[peak - 1] + a[peak] + a[peak + 1] + a[peak + 2] + a[peak + 3];

  //for scintillation type pulses
  //for (int i= std::max(peak - 100,2); i < std::min (peak + 800, 1023); ++i) integral += a[i];

  return integral;
}

//////////////////////////////////////////////
// Golden Configuration for Rising Edge Fit //
//////////////////////////////////////////////

void FitRisingEdge(TH1F* pulse, int nbinsL, int nbinsH, float &THM, float &risetime, float base, bool setbins ){

<<<<<<< HEAD
  TF1 *fTopline = new TF1("fTopline","pol0",600, 800);
  float x_min = pulse->GetBinCenter( pulse->GetMaximumBin()-1 );
  float x_max = pulse->GetBinCenter( pulse->GetMaximumBin()+1 );
  //std::cout << " min: " << x_min << " max: " << x_max << std::endl;
  pulse->Fit("fTopline","Q","", x_min, x_max);
  float top = fTopline->GetParameter(0);
  //std::cout << "top : " << top << " base: " << base << "\n";
=======
  //Get the location of the max value bin location
  int maxBin = pulse->GetMaximumBin();
  
  //want the center, not the edges to get "ave" value between max bin locaton and next highest before it
  float maxBinCenter = pulse -> GetBinCenter(maxBin);
  float prevBinCenter = pulse -> GetBinCenter(maxBin - 1);
  
  //fit the "top" line, that marks the peak value (actually an ave of two highest vals)
  TF1 *fTopline = new TF1("fTopline","pol0", prevBinCenter, maxBinCenter);
  pulse->Fit("fTopline","Q","", prevBinCenter, maxBinCenter);
  float top = fTopline->GetParameter(0); // get the max height
  //std::cout << "top : " << top << "\n";
>>>>>>> 7ace26a8e87e9b7269fa37c6080a6f279fa5cdfe

  double max=-9999;
  for (int i=10;i<pulse->GetXaxis()->GetNbins();i++) {    
    if (pulse->GetBinContent(i) > max) {
      max = pulse->GetBinContent(i);
    }
  }
  
  //std::cout << "max: " << pulse->GetMaximum() << " " << max << " " << min << "\n";

  int bM=0;  
  for (int i=10;i<pulse->GetXaxis()->GetNbins();i++) {
    //std::cout << "bin: "<< i << " " << pulse->GetBinContent(i) << "\n";
<<<<<<< HEAD
    if (pulse->GetBinContent(i) > base + 0.8*(top-base)) {
=======
    if (pulse->GetBinContent(i) > base + 0.9*(top-base)) {
>>>>>>> 7ace26a8e87e9b7269fa37c6080a6f279fa5cdfe
      bM = i;
      break;
    }
  }
  int bL = 0;
  for (int i=10;i<pulse->GetXaxis()->GetNbins();i++) {
    if (pulse->GetBinContent(i) > base + 0.1*(top-base)) {
      bL = i;
      break;
    }
  }
  
  //If want to restrict number of bins that are fit on rising edge
  bool RESTRICTBINS = setbins;
  /*
  if (RESTRICTBINS) {
    int BINS2FIT = 2; // want bM + bL = BINS2FIT
    
    if ((bM - bL) != (BINS2FIT - 1)) {
        //std::cout << bM << " , " << bL << std::endl;
        bL = bM - (BINS2FIT - 1);
        
    }
  }
  */


  //int bM = pulse->FindFirstBinAbove(0.6*pulse->GetMaximum());
  //int bL = pulse->FindFirstBinAbove(0.1*pulse->GetMaximum());
  //std::cout << bM << " : " << bL << "\n";
  //bM = 415;
  //bL = 435;
  //TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL+nbinsL), pulse->GetBinCenter(bM+nbinsH));
  TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL+nbinsL), pulse->GetBinCenter(bL + 3 + nbinsL ));
  //pulse->Fit(f,"MWLR");
  pulse->Fit(f,"RQ");
  float m = f->GetParameter(1);
  float b = f->GetParameter(0);
  delete f;
  
  //std::cout << m << " * x + " << b << std::endl;
  //std::cout << "HTM: " << 0.2*(0.2*pulse->GetMaximum()-b)/m << std::endl;
<<<<<<< HEAD
  THM = 0.2*(base+0.5*(top-base)-b)/m;//converted to picoseconds
  //risetime = 0.2*(base+0.75*(top-base)-b)/m - 0.2*(base+0.25*(top-base)-b)/m;
  risetime = 0.2*(0.75*top-0.25*top)/m;
  
=======
  
  THM = 0.2*(base+0.5*(top-base)-b)/m;//converted to nanoseconds
//  if (THM < 65 or THM > 66.1) {
//  std::cout << THM << "\n" << std::endl;}
  
  //risetime = 0.2*(base+0.7*(top-base)-b)/m - 0.2*(base+0.3*(top-base)-b)/m;
  risetime = 0.2*(0.9*top-0.1*top)/m;
>>>>>>> 7ace26a8e87e9b7269fa37c6080a6f279fa5cdfe
}


//For SiPMs on scintillator
// void FitRisingEdge(TH1F* pulse, int nbinsL, int nbinsH, float &THM, float &risetime ){
//   int bM = pulse->FindFirstBinAbove(0.5*pulse->GetMaximum());
//   int bL = pulse->FindFirstBinAbove(0.1*pulse->GetMaximum());
//   TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL+nbinsL), pulse->GetBinCenter(bM+nbinsH));
//   //pulse->Fit(f,"MWLR");
//   pulse->Fit(f,"RQ");
//   float m = f->GetParameter(1);
//   float b = f->GetParameter(0);
//   delete f;
//   //std::cout << "HTM: " << 0.2*(0.2*pulse->GetMaximum()-b)/m << std::endl;
//   THM = 0.2*(0.0*pulse->GetMaximum()-b)/m;//converted to picoseconds
//   risetime = 0.2*(0.6*pulse->GetMaximum()-b)/m - THM;  
// }


//For Hamamatsu on scintillator
// void FitRisingEdge(TH1F* pulse, int nbinsL, int nbinsH, float &THM, float &risetime ){

//   double firstPeakAmplitude = 0;
//   int firstpeakbin = 0;
//   for (int i=1; i<1024; i++) {
//     if ( pulse->GetBinContent(i) > 0.02
// 	 && pulse->GetBinContent(i) > pulse->GetBinContent(i-1)
// 	 && pulse->GetBinContent(i-1) > pulse->GetBinContent(i-2)
// 	 && pulse->GetBinContent(i-2) > pulse->GetBinContent(i-3)
// 	 && pulse->GetBinContent(i) > pulse->GetBinContent(i+1)
// 	 && pulse->GetBinContent(i+1) > pulse->GetBinContent(i+2)
// 	 ) {
//       firstPeakAmplitude = pulse->GetBinContent(i);
//       firstpeakbin = i;
//       break;
//     }
//   }
  
// //   std::cout << "FirstPeak: " << firstpeakbin << " " << firstPeakAmplitude << "\n";


//   int bM = pulse->FindFirstBinAbove(0.7*firstPeakAmplitude);
//   int bL = pulse->FindFirstBinAbove(0.1*firstPeakAmplitude);
//   TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL+nbinsL), pulse->GetBinCenter(bM+nbinsH));
//   //pulse->Fit(f,"MWLR");
//   pulse->Fit(f,"RQ");
//   float m = f->GetParameter(1);
//   float b = f->GetParameter(0);
//   delete f;
//   //std::cout << "HTM: " << 0.2*(0.2*pulse->GetMaximum()-b)/m << std::endl;
//   THM = 0.2*(0.1*firstPeakAmplitude-b)/m;//converted to picoseconds
//   risetime = 0.2*(0.6*firstPeakAmplitude-b)/m - THM;  
// }



// float FitRisingEdge(TH1F* pulse, int nbinsL, int nbinsH){
//   int bM = pulse->FindFirstBinAbove(0.6*pulse->GetMaximum());
//   int bL = pulse->FindFirstBinAbove(0.1*pulse->GetMaximum());
//   TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL+nbinsL), pulse->GetBinCenter(bM+nbinsH));
//   //pulse->Fit(f,"MWLR");
//   pulse->Fit(f,"RQ");
//   float m = f->GetParameter(1);
//   float b = f->GetParameter(0);
//   delete f;
//   //std::cout << "HTM: " << 0.2*(0.2*pulse->GetMaximum()-b)/m << std::endl;
//   return  0.2*(0.2*pulse->GetMaximum()-b)/m;//converted to picoseconds
  
// }

// float FitRisingEdge(TH1F* pulse, float Max){
//   int bM = pulse->FindFirstBinAbove(0.8*Max);
//   int bL = pulse->FindFirstBinAbove(0.125*Max);
//   TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL), pulse->GetBinCenter(bM));
//   //pulse->Fit(f,"MWLR");
//   pulse->Fit(f,"RQ");
//   float m = f->GetParameter(1);
//   float b = f->GetParameter(0);
//   delete f;
//   //std::cout << "HTM: " << 0.2*(0.3*pulse->GetMaximum()-b)/m << std::endl;
//   return  0.2*(0.5*pulse->GetMaximum()-b)/m;//converted to picoseconds
// }

// void FitRisingEdge(TH1F* pulse, float Max, float &THM, float &t0, float baseline){
//   int bM = pulse->FindFirstBinAbove(0.8*Max);
//   int bL = pulse->FindFirstBinAbove(0.125*Max);
//   TF1* f = new TF1("f", "[0]+x*[1]", pulse->GetBinCenter(bL), pulse->GetBinCenter(bM));
//   //pulse->Fit(f,"MWLR");
//   pulse->Fit(f,"RQ");
//   float m = f->GetParameter(1);
//   float b = f->GetParameter(0);
//   THM = 0.20*(0.15*pulse->GetMaximum()-b)/m;//converted to picoseconds
//   t0 = 0.20*(baseline-b)/m;//converted to picoseconds
//   delete f;
// }

// void FitFullPulse(TH1F* pulse, float &par0, float &par1, float &par2){
//   //TF1* f = new TF1("f","[2]+[1]*0.004217/2*exp(0.004217/2*(2.0*357.0+1.0*2.06**2.0-2.0*(x-[0])))*ROOT::Math::erfc((357.0+0.004217*2.06**2.0-(x-[0]))/(1.41*2.06))",10,1000);
  
//   TF1 *f = new TF1("f","[2]+[1]*0.004568/2*exp(0.004568/2*(2.0*260.0+1.0*5.0**2.0-2.0*(x-[0])))*ROOT::Math::erfc((260.0+0.004568*5.0**2.0-(x-[0]))/(1.41*5.0))",0,1000);
  
//   f->SetParameter(0, 40.0);
//   f->SetParameter(1, 60.0);
//   f->SetParameter(0, 0.0);
//   //pulse->Fit(f,"MWLR");
//   pulse->Fit(f,"RQ");
//   par0 = f->GetParameter(0);
//   par1 = f->GetMaximum(10,1000);
//   par2 = f->GetX(10,0.3*par1);
//   delete f;
// }


float LED( TH1F * pulse, double threshold, int nsamples, int splineBinFactor ) {

  double crosstime;
  int bin=0;
     
  // first sample above thresh (in absolute value)
  while( pulse->GetBinContent(bin) < threshold && bin< nsamples*splineBinFactor)
      bin++;

    //linear interpolation
    if(bin < nsamples*splineBinFactor){
      double inf = pulse->GetBinContent(bin-1);
      double sup = pulse->GetBinContent(bin);		
      crosstime = pulse->GetXaxis()->GetBinCenter(bin-1) + ((threshold-inf)/(sup-inf))*(pulse->GetXaxis()->GetBinCenter(bin)-pulse->GetXaxis()->GetBinCenter(bin-1));
    }
    else{
      crosstime = 1000000.;
      std::cerr << "signal does not reach threshold\n";
    }

  return crosstime;

}
