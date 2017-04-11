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
  float ND = 0.;
  std::cout << "[Usage]: ./SiPMvsPhotekNoAmp filename.root <AmpOn/Off> <Attenuator0/10/20dB> <ND0.5/0.8>\n" << std::endl;
  std::cout << "[Usage]: ./SiPMvsPhotekNoAmp ../../Output/20170407AmpOnOutput/1x1mm_SiPMon20dbAtt_71V_PhotekOn_3p5kV_AmpOn_LaserOn_Tune30_ND0p5_Output.root On 20 0.5\n" << std::endl;
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
        if( argv[2] == "AmpOn" || argv[2] == "On" ) amp = 36;
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
        if( argv[2] == "AmpOn" || argv[2] == "On" ) amp = 36;
        Att = atoi( argv[3] );

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
        if( argv[2] == "AmpOn" || argv[2] == "On" ) amp = 36;
        Att = atoi( argv[3] );
        ND = atof( argv[4] );

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
  t->Loop(amp,Att,ND);
  return 0;
}
