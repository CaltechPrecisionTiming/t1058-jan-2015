#include "TObject.h"
#include <math.h>
#include <string>
#include <TH1I.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <iostream>

#include <vector>


using namespace::std;


class TWaveForm : public TObject {
private:
  int     nb_samp; 	   //number of samples recorded
  double  t0;
  double* samp_x;
  double* samp_y;	   //arrays of samples
  

public:
  string message;	//for notification of error, exception, history of operations...
  // creator, with various types valid for array initialization   
  TWaveForm(int Nb_samp=0, double* SampX=NULL, double* SampY=NULL,string mes="");

  TWaveForm(const TWaveForm &wfm);	        // copy creator 
  TWaveForm& operator=(const TWaveForm &wfm);	// assignment
  //destruction
  ~TWaveForm();

  // getters
  int     Get_nb_samp() const;
  double  Get_t0() const;
  double* Get_samp_x() const;
  double* Get_samp_y() const;
  
  // graphics in time or frequency domain
  TGraph* Graph() ;
  TGraphErrors* GraphErrors() ;
  // void Graph( TGraph*) const;
  //  void Draw(const Option_t* option);
  

  
  // overloaded arithmetic operators
  TWaveForm& operator*=(const double &factor);  
  TWaveForm& operator/=(const double &factor);  
  TWaveForm& operator+=(double offset);
  TWaveForm& operator-=(double offset);
  
  TWaveForm& operator+=(const TWaveForm &wfm);
  TWaveForm& operator-=(const TWaveForm &wfm);
  
  // resampling	(nb: create a new TWaveForm on heap and return the pointer)
  TWaveForm* DownSample(int factor, int shift=0) const;
  TWaveForm* UpSample(int factor) const;
  TWaveForm* UpSample2(double step) const;
  TWaveForm* LowPass(int n_filter, int n_interp) const;
  TWaveForm* LowPass2(int n_filter) const;
  TWaveForm* Resize(double start, double stop) const;
  TWaveForm* Rising() const;
  TWaveForm* Rise() const;
  TWaveForm* Rise2() const;

  // amplitude (max sample, integrated pulse)
  double Amplitude() const;
  double Peaktime() const;
  double Risetime() const;
  double Risetime19() const;
  double Charge(double ohm=50) const;
  double Distance10( const TWaveForm*) const;
  double Distance( const TWaveForm*) const;
  double Distance( const TWaveForm*, TH1F* ) const;
  double Distance2( const TWaveForm*) const;
  double Distance2( const TWaveForm*, TH1F* ) const;
  double Distance3( const TWaveForm*) const;
  double Distance3( const TWaveForm*, TH1F* ) const;

  double Chisquare( const TWaveForm*, double, double, int&, double*) const;
  double Zero_cross() const;
  // double Chisquare( const TWaveForm*, double, double, int&) const;

  // time: leading-edge discriminator, constant fraction discriminator
  double* LED(double* thresh, int n_thresh=1);
  double* CFD(double* thresh, int n_thresh=1);
  double Fit1(double th1, double th2);
  

  // incrementation of a histogram (for noise, pedestal characterization)
  void FillHist(TH1I* histo) const;	
};

