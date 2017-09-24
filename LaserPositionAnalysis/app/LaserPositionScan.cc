#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>//std::pair
//ROOT INCLUDES
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
//LOCAL INCLUDES
#include "tree.hh"

std::string ParseCommandLine( int argc, char* argv[], std::string opt )
{
  for (int i = 1; i < argc; i++ )
    {
      std::string tmp( argv[i] );
      if ( tmp.find( opt ) != std::string::npos )
        {
          return tmp.substr( tmp.find_last_of("=") + 1 );
        }
    }

  return "";
};

float GetPosition( std::string pos )
{

  if ( pos == "0p0" ) return .0;
  if ( pos == "1p2" ) return 1.2;
  if ( pos == "2p4" ) return 2.4;
  if ( pos == "3p6" ) return 3.6;
  if ( pos == "4p8" ) return 4.8;
  if ( pos == "6p0" ) return 6.0;
  if ( pos == "7p2" ) return 7.2;
  if ( pos == "8p4" ) return 8.4;
  if ( pos == "9p6" ) return 9.6;
  if ( pos == "10p8" ) return 10.8;
  if ( pos == "12p0" ) return 12.0;
  if ( pos == "13p2" ) return 13.2;
  if ( pos == "14p4" ) return 14.4;
  if ( pos == "15p6" ) return 15.6;
  if ( pos == "16p8" ) return 16.8;
  else return -666.;
}

int main ( int argc, char** argv )
{

  //--------------------
  //Input List
  //--------------------
  std::string inputList = ParseCommandLine( argc, argv, "-inputList=" );
  if (  inputList == "" )
    {
      std::cerr << "[ERROR]: please provide an input list using --inputList=<inputList>" << std::endl;
      return -1;
    }


  float position[20];
  float positionErr[20];
  float amp[20];
  float ampErr[20];
  float meanTime[20];
  float meanTimeErr[20];
  float timeResolution[20];
  float timeResolutionErr[20];
  
  TFile* file = NULL;
  TTree* mytree = NULL;
  
  tree* loopTree = NULL;
    
  std::ifstream ifs( inputList.c_str(), std::ifstream::in );
  std::string currentFileName;
  int npoints = 0;
  if ( ifs.is_open() )
    {
      while ( ifs.good() )
	{
	  ifs >> currentFileName;
	  if ( ifs.eof() ) break;
	  std::cout << currentFileName << std::endl;

	  int init_pos = currentFileName.find("373_laser_");
	  //int init_pos = currentFileName.find("407_laser_");
	  int last_pos = currentFileName.find("_ANA.root");
	  std::string pos = currentFileName.substr( init_pos+10, last_pos-(init_pos+10));
	  std::cout << "pos: " << pos << std::endl;
	  position[npoints]    = GetPosition( pos );
	  positionErr[npoints] = 0.6;
	  file =  new TFile( currentFileName.c_str(), "READ");
	  mytree = (TTree*)file->Get("tree");
    
	  loopTree = new tree( mytree );
	  std::pair<float,float> myAmp = loopTree->GetAmp( pos.c_str() );
	  amp[npoints] = myAmp.first;
	  ampErr[npoints] = myAmp.second;
	  std::pair<float,float> myMeanTime = loopTree->GetMeanTime( pos.c_str() );
	  meanTime[npoints]    = myMeanTime.first;
	  meanTimeErr[npoints] = myMeanTime.second;
	  std::pair<float,float> myTimeResolution = loopTree->GetTimeResolution( pos.c_str() );
	  timeResolution[npoints]    = myTimeResolution.first*1000.;
	  timeResolutionErr[npoints] = myTimeResolution.second*1000.;
	  delete loopTree;
	  npoints++;
	  
	}
    }
  else
    {
      std::cout << "[ERROR] unable to open file; quitting" << std::endl;
    }

  TGraphErrors* graphAmp = new TGraphErrors(npoints+1, position, amp, positionErr, ampErr);

  TCanvas* c = new TCanvas( "c", "c", 800, 600 );
  c->SetHighLightColor(2);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);

  graphAmp->SetTitle("");
  graphAmp->GetXaxis()->SetTitleOffset(0.88);
  graphAmp->GetXaxis()->SetTitleSize(0.05);
  graphAmp->GetYaxis()->SetTitleSize(0.05);
  graphAmp->GetYaxis()->SetTitleOffset(0.89);
  graphAmp->GetXaxis()->SetTitle("Position [mm]");
  graphAmp->GetYaxis()->SetTitle("Amplitude [V]");
  graphAmp->SetMarkerStyle(20);
  graphAmp->SetMarkerColor(kBlue);
  graphAmp->SetMarkerSize(1.3);
  graphAmp->SetLineColor(kBlue);
  graphAmp->GetYaxis()->SetRangeUser(0,0.45);
  graphAmp->Draw("AP");
  c->SaveAs("amp_vs_position.pdf");


  TGraphErrors* graphMeanTime = new TGraphErrors(npoints+1, position, meanTime, positionErr, meanTimeErr);
  graphMeanTime->SetTitle("");
  graphMeanTime->GetXaxis()->SetTitleOffset(0.88);
  graphMeanTime->GetXaxis()->SetTitleSize(0.05);
  graphMeanTime->GetYaxis()->SetTitleSize(0.05);
  graphMeanTime->GetYaxis()->SetTitleOffset(0.89);
  graphMeanTime->GetXaxis()->SetTitle("Position [mm]");
  graphMeanTime->GetYaxis()->SetTitle("Mean time [ns]");
  graphMeanTime->SetMarkerStyle(20);
  graphMeanTime->SetMarkerColor(kBlue);
  graphMeanTime->SetMarkerSize(1.3);
  graphMeanTime->SetLineColor(kBlue);
  //graphMeanTime->GetYaxis()->SetRangeUser(27.5,28.5);//407 nm
  graphMeanTime->GetYaxis()->SetRangeUser(24.5,28.5);//373 nm
  graphMeanTime->Draw("AP");
  c->SaveAs("meanTime_vs_position.pdf");
  

  TGraphErrors* graphTimeResolution = new TGraphErrors(npoints+1, position, timeResolution, positionErr, timeResolutionErr);
  graphTimeResolution->SetTitle("");
  graphTimeResolution->GetXaxis()->SetTitleOffset(0.88);
  graphTimeResolution->GetXaxis()->SetTitleSize(0.05);
  graphTimeResolution->GetYaxis()->SetTitleSize(0.05);
  graphTimeResolution->GetYaxis()->SetTitleOffset(0.89);
  graphTimeResolution->GetXaxis()->SetTitle("Position [mm]");
  graphTimeResolution->GetYaxis()->SetTitle("#sigma_{t} [ps]");
  graphTimeResolution->SetMarkerStyle(20);
  graphTimeResolution->SetMarkerColor(kBlue);
  graphTimeResolution->SetMarkerSize(1.3);
  graphTimeResolution->SetLineColor(kBlue);
  //graphTimeResolution->GetYaxis()->SetRangeUser(0,20);//407 nm
  graphTimeResolution->GetYaxis()->SetRangeUser(0,60);//373 nm
  graphTimeResolution->Draw("AP");
  c->SaveAs("timeResolution_vs_position.pdf");
  /*
    TFile* file =  new TFile("/Users/cmorgoth/Work/data/cptLab/LaserLYSOSCAN_CMORGOTH/373_laser_3p6_ANA.root");
    TTree* mytree = (TTree*)file->Get("tree");
    
    tree* loopTree = new tree( mytree );
    loopTree->Loop();
  */
  return 0;
}
