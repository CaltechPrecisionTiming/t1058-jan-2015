#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <string.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include <tree.hh>

using namespace std;

int main( int argc, char** argv )
{
  ifstream file;// read file directly
  //default values
  int amp = 0;
  int Att = 0;
  //float ND = 0.;
  float left_ch1 = 0.;
  float right_ch1 = 9999.;
  std::cout << "[Usage]: ./SiPMvsPhotekNoAmp filename.root <Amp0/36dB> <Attenuator0/10/20dB> <ND0p5/0p8> <left_ch1Int> <right_ch1Int>\n" << std::endl;
  std::cout << "[Usage]: ./SiPMvsPhotekNoAmp ../../Output/20170407AmpOnOutput/1x1mm_SiPMon20dbAtt_71V_PhotekOn_3p5kV_AmpOn_LaserOn_Tune30_ND0p5_Output.root 36 20 ND0p5_peak0 0. 9999.\n" << std::endl;
  if (argc == 2)
  {
        //---------------------------------------------------------------
        //Open file and assume there is no attenuator 
        //---------------------------------------------------------------
        file.open (argv[1], ios::in | ios::binary);
        std::cout << "[INFO]: Opening file " << argv[1] << " ......" << std::endl;
        std::cout<< std::endl;
        if ( !file.is_open () )
        {
                std::cerr << "!! File open error:" << argv[1] << "; make sure the file is in the correct location" << std::endl;
                return 1;
        }
        std::cout << "[INFO]: Amplifier is  " << argv[2] << " dB ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: Attenuator is  " << argv[3] << " dB ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: There is no ND filter    ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: ch1Int left side is  " << left_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: ch1Int right side is  " << right_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;
  }
  else if ( argc == 3 )
  { 
        file.open (argv[1], ios::in | ios::binary);
        std::cout << "[INFO]: Opening file " << argv[1] << " ......" << std::endl;
        std::cout<< std::endl;
        if ( !file.is_open () )
        {
                std::cerr << "!! File open error:" << argv[1] << "; make sure the file is in the correct location" << std::endl;
                return 1;
        }
        amp = atoi( argv[2] );
        std::cout << "[INFO]: Amplifier is  " << amp << " dB ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: Attenuator is  " << Att << " dB ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: There is no ND filter    ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: ch1Int left side is  " << left_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: ch1Int right side is  " << right_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;
  }
  else if ( argc == 4 )
  {
        file.open (argv[1], ios::in | ios::binary);
        std::cout << "[INFO]: Opening file " << argv[1] << " ......" << std::endl;
        std::cout<< std::endl;
        if ( !file.is_open () )
        {
                std::cerr << "!! File open error:" << argv[1] << "; make sure the file is in the correct location" << std::endl;
                return 1;
        }
        amp = atoi( argv[2] );
        std::cout << "[INFO]: Amplifier is  " << amp << " dB ......" << std::endl;
        std::cout<< std::endl;
        Att = atoi( argv[3] );
        std::cout << "[INFO]: Attenuator is  " << Att << " dB ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: There is no ND filter    ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: ch1Int left side is  " << left_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: ch1Int right side is  " << right_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;

  }
  else if ( argc == 5 )
  {
        file.open (argv[1], ios::in | ios::binary);
        std::cout << "[INFO]: Opening file " << argv[1] << " ......" << std::endl;
        std::cout<< std::endl;
        if ( !file.is_open () )
        {
                std::cerr << "!! File open error:" << argv[1] << "; make sure the file is in the correct location" << std::endl;
                return 1;
        }
        amp = atoi( argv[2] );
        std::cout << "[INFO]: Amplifier is  " << amp << " dB ......" << std::endl;
        std::cout<< std::endl;
        Att = atoi( argv[3] );
        std::cout << "[INFO]: Attenuator is  " << Att << " dB ......" << std::endl;
        std::cout<< std::endl;
        //ND = atof( argv[4] );
        std::cout << "[INFO]: ND filter is  " << argv[4] << "  ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: ch1Int left side is  " << left_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: ch1Int right side is  " << right_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;

  }
  else if ( argc == 6 )
  {
        file.open (argv[1], ios::in | ios::binary);
        std::cout << "[INFO]: Opening file " << argv[1] << " ......" << std::endl;
        std::cout<< std::endl;
        if ( !file.is_open () )
        {
                std::cerr << "!! File open error:" << argv[1] << "; make sure the file is in the correct location" << std::endl;
                return 1;
        }
        amp = atoi( argv[2] );
        std::cout << "[INFO]: Amplifier is  " << amp << " dB ......" << std::endl;
        std::cout<< std::endl;
        Att = atoi( argv[3] );
        std::cout << "[INFO]: Attenuator is  " << Att << " dB ......" << std::endl;
        std::cout<< std::endl;
        //ND = atof( argv[4] );
        std::cout << "[INFO]: ND filter is  " << argv[4] << "  ......" << std::endl;
        std::cout<< std::endl;
        left_ch1 = (float)atof( argv[5] );
        std::cout << " left: " << left_ch1 << std::endl;
        std::cout << "[INFO]: ch1Int left side is  " << left_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;
        std::cout << "[INFO]: ch1Int right side is  " << right_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;
  }
  else if ( argc == 7 )
  {
        file.open (argv[1], ios::in | ios::binary);
        std::cout << "[INFO]: Opening file " << argv[1] << " ......" << std::endl;
        std::cout<< std::endl;
        if ( !file.is_open () )
        {
                std::cerr << "!! File open error:" << argv[1] << "; make sure the file is in the correct location" << std::endl;
                return 1;
        }
        amp = atoi( argv[2] );
        std::cout << "[INFO]: Amplifier is  " << amp << " dB ......" << std::endl;
        std::cout<< std::endl;
        Att = atoi( argv[3] );
        std::cout << "[INFO]: Attenuator is  " << Att << " dB ......" << std::endl;
        std::cout<< std::endl;
        //ND = atof( argv[4] );
        std::cout << "[INFO]: ND filter is  " << argv[4] << "  ......" << std::endl;
        std::cout<< std::endl;
        left_ch1 = (float)atof( argv[5] );
        std::cout << "[INFO]: ch1Int left side is  " << left_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;
        right_ch1 = (float)atof( argv[6] );
        std::cout << "[INFO]: ch1Int right side is  " << right_ch1 << "  ......" << std::endl;
        std::cout<< std::endl;
  }
  else
  {
          cerr << "[ERROR]!! No input file, please provide one" << endl;
          return 1;
  }

  string filename = argv[1];
  TFile *f = new TFile ((char *) filename.c_str (), "read"); 
  //TFile* f = new TFile("/Users/MJIAJING/precisiontiming/DAQ/Output/20170407AmpOnOutput/1x1mm_SiPMon20dbAtt_71V_PhotekOn_3p5kV_AmpOn_LaserOn_Tune30_ND0p5_Output.root", "read");
  TTree* myTree = (TTree*)f->Get("tree");
  std::cout << "here" << std::endl;
  tree* t = new tree( myTree );
  std::cout << "here" << std::endl;
  //std::cout << "amp: " << amp << " Att: " << Att << " left: " << left_ch1 << " right: " << right_ch1 << std::endl;
  t->Loop(amp,Att,argv[4],left_ch1,right_ch1);
  return 0;
}
