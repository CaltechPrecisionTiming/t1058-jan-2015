#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include <tree.hh>

int main( int argc, char** argv )
{
  TFile* f = new TFile("/Users/cmorgoth/Work/git/t1058-jan-2015/1x1mm_SiPMon_71V_PhotekOn_3p5kV_AmpOff_LaserOn_Tune30_ND1p8_ANA.root", "read");
  TTree* myTree = (TTree*)f->Get("tree");
  std::cout << "here" << std::endl;
  tree* t = new tree( myTree );
  std::cout << "here" << std::endl;
  t->Loop();
  return 0;
}
