//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Sep 23 11:45:37 2017 by ROOT version 6.10/02
// from TTree tree/tree
// found on file: /Users/cmorgoth/Work/data/cptLab/LaserLYSOSCAN_CMORGOTH/373_laser_3p6_ANA.root
//////////////////////////////////////////////////////////

#ifndef tree_h
#define tree_h

#include <utility>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          event;
   Float_t         t1gausroot;
   Float_t         t2gausroot;
   Float_t         t3gausroot;
   Float_t         t4gausroot;
   Float_t         ch1Amp;
   Float_t         ch2Amp;
   Float_t         ch3Amp;
   Float_t         ch4Amp;
   Float_t         ch1THM;
   Float_t         ch2THM;
   Float_t         ch3THM;
   Float_t         ch4THM;
   Float_t         ch1Risetime;
   Float_t         ch2Risetime;
   Float_t         ch3Risetime;
   Float_t         ch4Risetime;
   Float_t         ch1BL;
   Float_t         ch2BL;
   Float_t         ch3BL;
   Float_t         ch4BL;
   Float_t         ch1_TFF;
   Float_t         ch2_TFF;
   Float_t         ch3_TFF;
   Float_t         ch4_TFF;
   Float_t         ch1_TFF_v2;
   Float_t         ch2_TFF_v2;
   Float_t         ch3_TFF_v2;
   Float_t         ch4_TFF_v2;
   Float_t         ch1_AFF;
   Float_t         ch2_AFF;
   Float_t         ch3_AFF;
   Float_t         ch4_AFF;
   UInt_t          ch1QualityBit;
   UInt_t          ch2QualityBit;
   UInt_t          ch3QualityBit;
   UInt_t          ch4QualityBit;
   Float_t         ch1Int;
   Float_t         ch2Int;
   Float_t         ch3Int;
   Float_t         ch4Int;
   Float_t         ch1chisq;
   Float_t         ch2chisq;
   Float_t         ch3chisq;
   Float_t         ch4chisq;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_t1gausroot;   //!
   TBranch        *b_t2gausroot;   //!
   TBranch        *b_t3gausroot;   //!
   TBranch        *b_t4gausroot;   //!
   TBranch        *b_ch1Amp;   //!
   TBranch        *b_ch2Amp;   //!
   TBranch        *b_ch3Amp;   //!
   TBranch        *b_ch4Amp;   //!
   TBranch        *b_ch1THM;   //!
   TBranch        *b_ch2THM;   //!
   TBranch        *b_ch3THM;   //!
   TBranch        *b_ch4THM;   //!
   TBranch        *b_ch1Risetime;   //!
   TBranch        *b_ch2Risetime;   //!
   TBranch        *b_ch3Risetime;   //!
   TBranch        *b_ch4Risetime;   //!
   TBranch        *b_ch1BL;   //!
   TBranch        *b_ch2BL;   //!
   TBranch        *b_ch3BL;   //!
   TBranch        *b_ch4BL;   //!
   TBranch        *b_ch1_TFF;   //!
   TBranch        *b_ch2_TFF;   //!
   TBranch        *b_ch3_TFF;   //!
   TBranch        *b_ch4_TFF;   //!
   TBranch        *b_ch1_TFF_v2;   //!
   TBranch        *b_ch2_TFF_v2;   //!
   TBranch        *b_ch3_TFF_v2;   //!
   TBranch        *b_ch4_TFF_v2;   //!
   TBranch        *b_ch1_AFF;   //!
   TBranch        *b_ch2_AFF;   //!
   TBranch        *b_ch3_AFF;   //!
   TBranch        *b_ch4_AFF;   //!
   TBranch        *b_ch1QualityBit;   //!
   TBranch        *b_ch2QualityBit;   //!
   TBranch        *b_ch3QualityBit;   //!
   TBranch        *b_ch4QualityBit;   //!
   TBranch        *b_ch1Int;   //!
   TBranch        *b_ch2Int;   //!
   TBranch        *b_ch3Int;   //!
   TBranch        *b_ch4Int;   //!
   TBranch        *b_ch1chisq;   //!
   TBranch        *b_ch2chisq;   //!
   TBranch        *b_ch3chisq;   //!
   TBranch        *b_ch4chisq;   //!

   tree(TTree *tree=0);
   virtual ~tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
  virtual std::pair<float,float> GetTimeResolution( TString pos = "default");
  virtual std::pair<float,float> GetMeanTime( TString pos = "default");
  virtual std::pair<float,float> GetAmp(TString pos = "default");
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tree_cxx
tree::tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/cmorgoth/Work/data/cptLab/LaserLYSOSCAN_CMORGOTH/373_laser_3p6_ANA.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/Users/cmorgoth/Work/data/cptLab/LaserLYSOSCAN_CMORGOTH/373_laser_3p6_ANA.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

tree::~tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("t1gausroot", &t1gausroot, &b_t1gausroot);
   fChain->SetBranchAddress("t2gausroot", &t2gausroot, &b_t2gausroot);
   fChain->SetBranchAddress("t3gausroot", &t3gausroot, &b_t3gausroot);
   fChain->SetBranchAddress("t4gausroot", &t4gausroot, &b_t4gausroot);
   fChain->SetBranchAddress("ch1Amp", &ch1Amp, &b_ch1Amp);
   fChain->SetBranchAddress("ch2Amp", &ch2Amp, &b_ch2Amp);
   fChain->SetBranchAddress("ch3Amp", &ch3Amp, &b_ch3Amp);
   fChain->SetBranchAddress("ch4Amp", &ch4Amp, &b_ch4Amp);
   fChain->SetBranchAddress("ch1THM", &ch1THM, &b_ch1THM);
   fChain->SetBranchAddress("ch2THM", &ch2THM, &b_ch2THM);
   fChain->SetBranchAddress("ch3THM", &ch3THM, &b_ch3THM);
   fChain->SetBranchAddress("ch4THM", &ch4THM, &b_ch4THM);
   fChain->SetBranchAddress("ch1Risetime", &ch1Risetime, &b_ch1Risetime);
   fChain->SetBranchAddress("ch2Risetime", &ch2Risetime, &b_ch2Risetime);
   fChain->SetBranchAddress("ch3Risetime", &ch3Risetime, &b_ch3Risetime);
   fChain->SetBranchAddress("ch4Risetime", &ch4Risetime, &b_ch4Risetime);
   fChain->SetBranchAddress("ch1BL", &ch1BL, &b_ch1BL);
   fChain->SetBranchAddress("ch2BL", &ch2BL, &b_ch2BL);
   fChain->SetBranchAddress("ch3BL", &ch3BL, &b_ch3BL);
   fChain->SetBranchAddress("ch4BL", &ch4BL, &b_ch4BL);
   fChain->SetBranchAddress("ch1_TFF", &ch1_TFF, &b_ch1_TFF);
   fChain->SetBranchAddress("ch2_TFF", &ch2_TFF, &b_ch2_TFF);
   fChain->SetBranchAddress("ch3_TFF", &ch3_TFF, &b_ch3_TFF);
   fChain->SetBranchAddress("ch4_TFF", &ch4_TFF, &b_ch4_TFF);
   fChain->SetBranchAddress("ch1_TFF_v2", &ch1_TFF_v2, &b_ch1_TFF_v2);
   fChain->SetBranchAddress("ch2_TFF_v2", &ch2_TFF_v2, &b_ch2_TFF_v2);
   fChain->SetBranchAddress("ch3_TFF_v2", &ch3_TFF_v2, &b_ch3_TFF_v2);
   fChain->SetBranchAddress("ch4_TFF_v2", &ch4_TFF_v2, &b_ch4_TFF_v2);
   fChain->SetBranchAddress("ch1_AFF", &ch1_AFF, &b_ch1_AFF);
   fChain->SetBranchAddress("ch2_AFF", &ch2_AFF, &b_ch2_AFF);
   fChain->SetBranchAddress("ch3_AFF", &ch3_AFF, &b_ch3_AFF);
   fChain->SetBranchAddress("ch4_AFF", &ch4_AFF, &b_ch4_AFF);
   fChain->SetBranchAddress("ch1QualityBit", &ch1QualityBit, &b_ch1QualityBit);
   fChain->SetBranchAddress("ch2QualityBit", &ch2QualityBit, &b_ch2QualityBit);
   fChain->SetBranchAddress("ch3QualityBit", &ch3QualityBit, &b_ch3QualityBit);
   fChain->SetBranchAddress("ch4QualityBit", &ch4QualityBit, &b_ch4QualityBit);
   fChain->SetBranchAddress("ch1Int", &ch1Int, &b_ch1Int);
   fChain->SetBranchAddress("ch2Int", &ch2Int, &b_ch2Int);
   fChain->SetBranchAddress("ch3Int", &ch3Int, &b_ch3Int);
   fChain->SetBranchAddress("ch4Int", &ch4Int, &b_ch4Int);
   fChain->SetBranchAddress("ch1chisq", &ch1chisq, &b_ch1chisq);
   fChain->SetBranchAddress("ch2chisq", &ch2chisq, &b_ch2chisq);
   fChain->SetBranchAddress("ch3chisq", &ch3chisq, &b_ch3chisq);
   fChain->SetBranchAddress("ch4chisq", &ch4chisq, &b_ch4chisq);
   Notify();
}

Bool_t tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tree_cxx
