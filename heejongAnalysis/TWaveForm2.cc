#include "TObject.h"
#include <math.h>
#include <string>
#include <TH1I.h>
#include <TGraph.h>
#include <iostream>
#include <algorithm>

#include "Math/Polynomial.h"
#include "Math/Interpolator.h"

#include "TWaveForm2.h"

//creation with optional initialization parameters
TWaveForm::TWaveForm(int Nsamp, double* SampX, double* SampY, string mes)
{
  nb_samp = Nsamp;
  t0 = SampX[0];

  // Dynamic allocation of sample array, and copy from input
  samp_x = new double[nb_samp];
  samp_y = new double[nb_samp];

  for( int i=0; i < nb_samp; i++){
    samp_x[i]=SampX[i];
    samp_y[i]=SampY[i];
  }

  message=mes;
}


// copy
TWaveForm::TWaveForm(const TWaveForm &wfm)
{
  nb_samp = wfm.nb_samp;
  t0 = wfm.t0;

  //Dynamic allocation of sample array, and copy from wfm
  samp_x = new double[nb_samp];
  samp_y = new double[nb_samp];

  for( int i=0; i<nb_samp; i++){
    samp_x[i]=wfm.samp_x[i];
    samp_y[i]=wfm.samp_y[i];
  }

  message=wfm.message;
}


// assignment
TWaveForm& TWaveForm::operator=(const TWaveForm &wfm)
{
  // if the size of the array is changed, delete it and reallocate	
  if( wfm.nb_samp != nb_samp){
    nb_samp = wfm.nb_samp;
    delete[] samp_x;
    delete[] samp_y;
    samp_x = new double[nb_samp];		
    samp_y = new double[nb_samp];		
  }

  //copy the values
  t0 = wfm.t0;

  for( int i=0; i<nb_samp; i++){
    samp_x[i]=wfm.samp_x[i];
    samp_y[i]=wfm.samp_y[i];
  }	

  message = wfm.message;

  return *this;
}


// deletion
TWaveForm::~TWaveForm()	
{
  delete[] samp_x; 
  delete[] samp_y; 
}


// getters
int TWaveForm::Get_nb_samp() const
{
  return nb_samp;
}

double TWaveForm::Get_t0() const
{
  return t0;
}

double* TWaveForm::Get_samp_x() const
{
  return samp_x;
}

double* TWaveForm::Get_samp_y() const
{
  return samp_y;
}

// graphics
TGraph* TWaveForm::Graph() 
// void TWaveForm::Graph( TGraph* gr) const
{
  TGraph* gr = new TGraph( nb_samp, samp_x, samp_y);
  // gr = new TGraph( nb_samp, samp_x, samp_y);

  return gr;
}

TGraphErrors* TWaveForm::GraphErrors() 
// void TWaveForm::Graph( TGraph* gr) const
{
  double et[nb_samp];
  for( int i = 0; i <  nb_samp; i++)
    et[i] = 0.;
  
  
  TGraphErrors* gr = new TGraphErrors( nb_samp, samp_x, samp_y, et, et);
  // gr = new TGraph( nb_samp, samp_x, samp_y);

  return gr;
}

//draw in time domain
/*
void TWaveForm::Draw( const Option_t* option)
{
  TGraph* gr = Graph();
  gr->Draw( option);
}
*/

//multiplication by a constant (unary)
TWaveForm& TWaveForm::operator*=(const double &factor)
{
  for( int i=0; i<nb_samp; i++)  
    samp_y[i] *= factor;

  message += "multiplied by factor\n";
  return *this;
}


//division by a constant (unary)
TWaveForm& TWaveForm::operator/=(const double &factor)
{ 
  for( int i=0; i<nb_samp; i++)  
    samp_y[i] /= factor;
  
  message+="divided by factor\n";
  return *this;
}

// incrementation by an offset
TWaveForm& TWaveForm::operator+=(double offset)
{
  for( int i=0; i<nb_samp; i++)  
    samp_y[i] += offset;

  message += "+offset applied\n";
  return *this;
}

// decrementation by an offset
TWaveForm& TWaveForm::operator-=(double offset)
{
  for( int i=0; i<nb_samp; i++)  
    samp_y[i]-=offset;

  message += "-offset applied\n";
  return *this;
}




// Incrementation/decrementation of waveforms by another waveform:
// the frequency and phase must match,
// but the lenghts and starting time can be different.
// The samples having the same time id are added

TWaveForm& TWaveForm::operator+=(const TWaveForm &wfm)
{

  if( t0 == wfm.t0){
    if( nb_samp == wfm.nb_samp){
      // simplest case, no reallocation
      for( int i=0; i< nb_samp;i++){
	samp_y[i] = samp_y[i] + wfm.samp_y[i];
	message += "addition OK, size unchanged\n";
      }
    }
    else{
      // select the smaller
      if( nb_samp > wfm.nb_samp) 
	nb_samp=wfm.nb_samp;	

      // allocate new temporary array
      double* tarray_x = new double[nb_samp];	
      double* tarray_y = new double[nb_samp];	
      for(int i=0;i<nb_samp;i++){		  
	tarray_x[i]=samp_x[i];
	tarray_y[i]=samp_y[i]+wfm.samp_y[i];
      }
      delete[] samp_x;
      delete[] samp_y;

      samp_x = tarray_x;
      samp_y = tarray_y;

      message += "addition OK, new size\n";
    }
    return *this;			
  }
  else{
   
    int n_sample = nb_samp + wfm.nb_samp;

    double* tarray_x = new double[n_sample];
    double* tarray_y = new double[n_sample];

    int i, j, index;
    i = j = index = 0;

    while( index < n_sample){
      
      if( samp_x[i] < wfm.samp_x[j]){
	tarray_x[index] = samp_x[i];
	tarray_y[index] = samp_y[i];
	i++;
      }
      else{
	tarray_x[index] = wfm.samp_x[j];
	tarray_y[index] = wfm.samp_y[j];
	j++;
      }

      index++;
    }

    delete[] samp_x;
    delete[] samp_y;
    
    samp_x = tarray_x;
    samp_y = tarray_y;

    nb_samp = n_sample;
    t0 = samp_x[0];
    
    message += "addition OK, size increased\n";

    return *this;			
  }
  
}


TWaveForm& TWaveForm::operator-=(const TWaveForm &wfm)
{
  //copy-creation
  TWaveForm wftemp = wfm;	

  //inversion
  wftemp *= -1.;	

  //add to *this
  (*this) += wftemp;	

  return *this;
}


// take sub samples every "factor" sample
// starting from the "offset"-th sample
TWaveForm* TWaveForm::DownSample(int factor, int offset) const
{
  string mesg="";	
  int n_d = (int)floor((nb_samp-1-offset)/factor)+1;

  if(n_d<0){
    n_d=0;
    mesg += "size is 0 !\n";
  }

  // sub samples
  double arr_x[n_d], arr_y[n_d];
  for( int i=0; i<n_d; i++){
    arr_x[i]=samp_x[offset+factor*i];
    arr_y[i]=samp_y[offset+factor*i];
  }	

  // output
  TWaveForm* w_out = new TWaveForm( n_d, arr_x, arr_y, mesg);	
  return w_out;
}


TWaveForm* TWaveForm::UpSample(int factor) const
{ 
 string mesg="";

 ROOT::Math::Interpolator cspline( nb_samp, ROOT::Math::Interpolation::kCSPLINE);
 cspline.SetData( nb_samp, samp_x, samp_y);

 int nb_new = nb_samp + (nb_samp-1)*(factor -1);

 double* out_x = new double[nb_new];
 double* out_y = new double[nb_new];
 // double* out_x = new double[nb_samp*factor];
 // double* out_y = new double[nb_samp*factor];

 for( int i = 0; i < nb_samp-1; i++)
   for( int j = 0; j < factor; j++){
     out_x[i*factor+j] = samp_x[i]+ ((samp_x[i+1] - samp_x[i])/factor)*j;
     out_y[i*factor+j] = cspline.Eval( out_x[i*factor+j]);
   }
 
 out_x[nb_new-1] = samp_x[nb_samp-1];
 out_y[nb_new-1] = cspline.Eval( out_x[nb_new-1]);

 TWaveForm* w_out = new TWaveForm( nb_new, out_x, out_y, mesg);
 

 delete out_x;
 delete out_y;

 return w_out;
}


TWaveForm* TWaveForm::UpSample2( double step) const
{ 
 string mesg="";

 ROOT::Math::Interpolator cspline( nb_samp, ROOT::Math::Interpolation::kCSPLINE);
 cspline.SetData( nb_samp, samp_x, samp_y);

 int nb_new = int((samp_x[nb_samp-1] - samp_x[0])/step);

 double* out_x = new double[nb_new];
 double* out_y = new double[nb_new];

 for( int i = 0; i < nb_new; i++){
   out_x[i] = samp_x[0]+ i*step;
   out_y[i] = cspline.Eval( out_x[i]);
 }
 
 TWaveForm* w_out = new TWaveForm( nb_new, out_x, out_y, mesg);
 
 delete out_x;
 delete out_y;

 return w_out;
}


TWaveForm* TWaveForm::Resize( double start, double stop) const
{
  string mesg="";

  int istart, istop;
  istart = istop = 0;

  for( int i = 0; i < nb_samp; i++){
    if( samp_x[i] < start) istart = i;
    if( samp_x[i] < stop ) istop = i+1;
  }
  
  int nb_arr = istop-istart+1;

  if( nb_arr <= 1){
    nb_arr = 1;
    double arr_x[1] = {0};
    double arr_y[1] = {0};

    TWaveForm* w_out = new TWaveForm( nb_arr, arr_x, arr_y, mesg);	
    return w_out;
  }


  double arr_x[nb_arr];
  double arr_y[nb_arr];

  for( int i=0; i < nb_arr; i++){
    arr_x[i] = samp_x[i+istart];
    arr_y[i] = samp_y[i+istart];
  }	

  // output
  TWaveForm* w_out = new TWaveForm( nb_arr, arr_x, arr_y, mesg);	
  return w_out;
}



TWaveForm* TWaveForm::Rising() const
{

  int i0 = 0;
  for( int i=0; i < nb_samp; i++)
    if(samp_y[i] > 20){
      i0 = i;
      break;
    }

  int i1 = 0;
  for( int i=i0; i < nb_samp; i++)
    if(samp_y[i] > samp_y[i+1]){
      i1 = i;
      break;
    }
  
  double* arr_x = new double[i1-i0+1];	
  double* arr_y = new double[i1-i0+1];

  for( int i=i0; i <= i1; i++){
    arr_x[i-i0] = samp_x[i];
    arr_y[i-i0] = samp_y[i];
  }	

  // output
  TWaveForm* w_out = new TWaveForm( i1-i0+1, arr_x, arr_y, "");	
  return w_out;
}





// amplitude: maximal absolute value
TWaveForm* TWaveForm::Rise() const 
{
  string mesg="";

  double threshold = 0.3*Amplitude();
  int th_id=0;
  
  // search threshold point
  for( int i  = 0; i < nb_samp; i++){
    if( samp_y[i] > threshold){
      th_id = i;
      break;
    }
  }
  
  int start_id = 0;
  for( int i  =  th_id; i > 0; i--){
    if( samp_y[i] > samp_y[i-1])
      start_id = i-1;
    else
      break;
  }

  int stop_id = 0;
  for( int i  =  th_id; i < nb_samp; i++){
    if( samp_y[i+1] > samp_y[i])
      stop_id = i+1;
    else
      break;
  }

  int nb_arr = stop_id-start_id+1;

  if( nb_arr <= 1){
    nb_arr = 1;
    double arr_x[1] = {0};
    double arr_y[1] = {0};

    TWaveForm* w_out = new TWaveForm( nb_arr, arr_x, arr_y, mesg);	
    return w_out;
  }

    

  //  printf("nb_arr  = %5d, amp = %7.3f,  stop = %5d,   start  =  %5d\n", 
  // nb_arr, Amplitude(), stop_id, start_id);

  double arr_x[nb_arr];
  double arr_y[nb_arr];

  for( int i=0; i < nb_arr; i++){
    arr_x[i] = samp_x[i+start_id];
    arr_y[i] = samp_y[i+start_id];

    //    printf("i = %5d,  %7.3f\n", i, arr_y[i]);
  }	


  

  // output
  TWaveForm* w_out = new TWaveForm( nb_arr, arr_x, arr_y, mesg);	
  return w_out;

}


// amplitude: maximal absolute value
TWaveForm* TWaveForm::Rise2() const 
{
  string mesg="";

  double threshold = 0.3*Amplitude();
  //  double threshold = 30;
  int th_id=0;
  
  // search threshold point
  for( int i  = 0; i < nb_samp; i++){
    if( samp_y[i] > threshold){
      th_id = i;
      break;
    }
  }
  
  int start_id = 0;
  for( int i  =  th_id; i > 0; i--){
    if( samp_y[i] > samp_y[i-1])
      start_id = i-1;
    else
      break;
  }

  int stop_id = 0;
  for( int i  =  th_id; i < nb_samp; i++){
    if( samp_y[i+1] > samp_y[i])
      stop_id = i+1;
    else
      break;
  }


  int stop_id2 = 0;
  for( int i  =  stop_id; i < nb_samp; i++){
    if( samp_y[i+1] < samp_y[i])
      stop_id2 = i+1;
    else
      break;
  }


  int nb_arr = stop_id2-start_id+1;

  if( nb_arr <= 1){
    nb_arr = 1;
    double arr_x[1] = {0};
    double arr_y[1] = {0};

    TWaveForm* w_out = new TWaveForm( nb_arr, arr_x, arr_y, mesg);	
    return w_out;
  }

  double arr_x[nb_arr];
  double arr_y[nb_arr];

  for( int i=0; i < nb_arr; i++){
    arr_x[i] = samp_x[i+start_id];
    arr_y[i] = samp_y[i+start_id];

    //    printf("i = %5d,  %7.3f\n", i, arr_y[i]);
  }	

  // output
  TWaveForm* w_out = new TWaveForm( nb_arr, arr_x, arr_y, mesg);	
  return w_out;

}


// low-pass filtering with factor f_filter combined with upsampling with factor f_interp
TWaveForm* TWaveForm::LowPass(int n_filter, int n_interp) const
{
  TWaveForm* wf_down = DownSample( n_filter, 0);
  TWaveForm* wf_up   = wf_down->UpSample( n_filter*n_interp);
  TWaveForm* wf_up0  = UpSample( n_interp);


  double* x_up0 = wf_up0->Get_samp_x(); 
  double* y_up0 = new double[wf_up0->Get_nb_samp()];
  double* x_up1 = new double[wf_up0->Get_nb_samp()];

  TWaveForm* out = new TWaveForm(*wf_up);	//creation with copy

  for( int i = 1; i < n_filter; i++){
    delete wf_down;			
    delete wf_up;			

    // decimation (n_filter interleaved downsampled signals)
    wf_down = DownSample( n_filter, i); 
    wf_up = wf_down->UpSample( n_filter*n_interp);

    (*out) += (*wf_up);
  }

  double *outx = out->Get_samp_x();
  double *outy = out->Get_samp_y();

  // re-arrange sampled points
  int last_i = 0;
  int count = 0;
  double part_sum = 0;
  double time_sum = 0;

  
  for( int j = 0; j < wf_up0->Get_nb_samp()-1; j++){
    int start = last_i;   
    for( int i = start ; i < out->Get_nb_samp(); i++){
    
      if( fabs( x_up0[j]- outx[i]) <  fabs( x_up0[j+1]- outx[i]) ){	
	part_sum += outy[i];
	time_sum += outx[i];
	count++;
	last_i++;
      }
      else{
	y_up0[j] = part_sum/count;
	x_up1[j] = time_sum/count;
	part_sum = 0;
	time_sum = 0;
	count = 0;
	break;
      }
    }
  
  }
  

  for( int i = last_i ; i < out->Get_nb_samp(); i++){        
    part_sum += outy[i];
    time_sum += outx[i];
    count++;      
  }

  y_up0[wf_up0->Get_nb_samp()-1] = part_sum/count;
  x_up1[wf_up0->Get_nb_samp()-1] = time_sum/count;

  TWaveForm* rout = new TWaveForm( wf_up0->Get_nb_samp(), x_up0, y_up0);
  // TWaveForm* rout = new TWaveForm( wf_up0->Get_nb_samp(), x_up1, y_up0);

  delete wf_down;
  delete wf_up;
  delete wf_up0;
  delete out;


  delete y_up0;
  delete x_up1;

  return rout;
}



// low-pass filtering with factor f_filter combined with upsampling with factor f_interp
TWaveForm* TWaveForm::LowPass2(int n_filter) const
{

  double* x0; 
  double* y0;
  
  x0 = Get_samp_x();
  y0 = Get_samp_y();

  double arr_x[nb_samp];
  double arr_y[nb_samp];
  double arr_y2[nb_samp];

  for( int i = 0; i < nb_samp; i++){
    arr_x[i] = x0[i];
    arr_y[i] = y0[i];
  }

  for( int i = 0+n_filter; i < nb_samp-n_filter; i++){
    double temp = 0;
    double count = 0;
    for( int j = i-n_filter; j < i+n_filter; j++){

      temp += arr_y[j];
      count += 1.0;
    }
    
    arr_y2[i] = temp/count;    
  }

  TWaveForm* rout = new TWaveForm( nb_samp, arr_x, arr_y2);

  return rout;
}






double TWaveForm::Amplitude() const
{
  /*
  double threshold = 30;
  int th_id=0;
  
  // search threshold point
  for( int i  = 0; i < nb_samp; i++){
    if( samp_y[i] > threshold){
      th_id = i;
      break;
    }
  }
  

  int id=0;

  // search the first peak
  for( int i = th_id; i<nb_samp; i++){
    if( samp_y[i+1] >= samp_y[i]){
      id = i+1;
    }
    else
      break;
  }


  printf("Amplitude = %7.3f(%5d),  th_id = %5d\n", samp_y[id], id, th_id);
  */

  int id=0;
  double max=-1000;

  /*
  for( int i=0; i <nb_samp - 30; i++)	
    if( max < fabs(samp_y[i])){
      max = fabs(samp_y[i]);
      id = i;
    }
  */

  for( int i=0; i <nb_samp - 30; i++)	
    if( max < samp_y[i]){
      max = samp_y[i];
      id = i;
    }


  return samp_y[id];	
}



double TWaveForm::Zero_cross() const
{

  double* x0 = new double[nb_samp];	
  double* y0 = new double[nb_samp];

  x0 = Get_samp_x();
  y0 = Get_samp_y();

  // decide the range using thresholds
  int i0, i1;
  i0 = i1 = 0;

  for( int i = 0; i < nb_samp; i++)
    if(y0[i] < 0){
      i0 = i-1;
      i1 = i;
      break;
    }
      
  double slope  = (y0[i1]-y0[i0])/(x0[i1]-x0[i0]);
  double proj_t = (slope*x0[i0] - y0[i0])/slope;

 
  return proj_t;
}



double TWaveForm::Peaktime() const
{
  /*
  double threshold = 30;
  int th_id=0;
  
  // search threshold point
  for( int i  = 0; i < nb_samp; i++){
    if( samp_y[i] > threshold){
      th_id = i;
      break;
    }
  }
  

  int id=0;

  // search the first peak
  for( int i = th_id; i<nb_samp; i++){
    if( samp_y[i+1] > samp_y[i]){
      id = i+1;
    }
    else
      break;
  }
  */

  int id = 0;
  double max = -1000;

  for( int i=0; i <nb_samp - 30; i++)	
    if( max < fabs(samp_y[i])){
      max = fabs(samp_y[i]);
      id = i;
    }
  

  return samp_x[id];	
}




// amplitude: maximal absolute value
double TWaveForm::Risetime() const 
{
  
  double threshold = 0.3*Amplitude();
  int th_id=0;
  
  // search threshold point
  for( int i  = 0; i < nb_samp; i++){
    if( samp_y[i] > threshold){
      th_id = i;
      break;
    }
  }
  
  int start_id = 0;
  for( int i  =  th_id; i > 0; i--){
    if( samp_y[i] > samp_y[i-1])
      start_id = i-1;
    else
      break;
  }

  int stop_id = 0;
  for( int i  =  th_id; i < nb_samp; i++){
    if( samp_y[i] > samp_y[i-1])
      stop_id = i;
    else
      break;
  }

  //  printf("Rise : stop = %7.3f, start = %7.3f\n", samp_x[stop_id], samp_x[start_id]);

  
  return (samp_x[stop_id]-samp_x[start_id]);	
}

double TWaveForm::Risetime19() const 
{
  
  double amp = Amplitude();
  double threshold = 0.3*amp;
  int th_id=0;
  
  // search threshold point
  for( int i  = 0; i < nb_samp; i++){
    if( samp_y[i] > threshold){
      th_id = i;
      break;
    }
  }
  
  int start_id = 0;
  for( int i  =  th_id; i > 0; i--){
    if( samp_y[i] > samp_y[i-1] && samp_y[i] > 0.1*amp)
      start_id = i-1;
    else
      break;
  }

  int stop_id = 0;
  for( int i  =  th_id; i < nb_samp; i++){
    if( samp_y[i] > samp_y[i-1] &&  samp_y[i] < 0.9*amp)
      stop_id = i;
    else
      break;
  }

  //  printf("Rise 19: stop = %7.3f, start = %7.3f\n", samp_x[stop_id], samp_x[start_id]);

  
  return (samp_x[stop_id]-samp_x[start_id]);	
}

// charge in a "ohm" load (default=50)
double TWaveForm::Charge(double ohm) const
{
  double sum = 0.;
  for( int i=0; i<nb_samp;i++)
    sum += samp_y[i];

  double dt = (samp_x[nb_samp-1]-samp_x[0])/(nb_samp-1);

  
  return (sum*dt)/ohm;
}


// leading-edge discriminator with several thresholds
double* TWaveForm::LED( double* thresh, int n_thresh)
{
  double* crosstime = new double[n_thresh];
  int id_samp=0;

  for( int i=0; i<n_thresh; i++){
    double thr = thresh[i];
    
    // first sample above thresh (in absolute value)
    while( samp_y[id_samp]<thr && id_samp<nb_samp)
      id_samp++;

    //linear interpolation
    if(id_samp < nb_samp) {
      double inf = samp_y[id_samp-1];
      double sup = samp_y[id_samp];		
      crosstime[i]=samp_x[id_samp-1] + ((thr-inf)/(sup-inf))*(samp_x[id_samp]-samp_x[id_samp-1]);
    }
    else{
      crosstime[i] = 1000000.;
      message += "signal does not reach threshold";
    }
  }
  

  return crosstime;
}

// constant fraction discriminator
// thresh is an array of relative thresholds (between 0 and 1)
double* TWaveForm::CFD(double* thresh, int n_thresh)
{	
  double amp = Amplitude();
  double abs_thresh[n_thresh];

  //calculate absolute threshold using pulse amplitude
  for( int i=0; i<n_thresh; i++){
    abs_thresh[i] = amp*thresh[i];
  }

  return LED( abs_thresh, n_thresh);
}

// constant fraction discriminator
// thresh is an array of relative thresholds (between 0 and 1)
double TWaveForm::Fit1(double th1, double th2)
{	
  TGraphErrors* g1 = GraphErrors();

  int istart, istop;
  istart = istop = 0;

  for( int i = 0; i < nb_samp; i++){
    if( samp_y[i] > th1){
      istart = i;
      break;
    }
  }

  for( int i = 0; i < nb_samp; i++){
    if( samp_y[i] > th2){
      istop = i;
      break;
    }
  }


  g1->Fit("pol1", "QE0", "", samp_x[istart-1], samp_x[istop]); 

  double aa = g1->GetFunction("pol1")->GetParameter(1);
  double bb = g1->GetFunction("pol1")->GetParameter(0);

  // double tt = (0.-bb)/aa;
  double tt = ( ((samp_x[istop]-samp_x[istart-1])/2.+samp_x[istop])-bb)/aa;
  delete g1;

  return tt;
}



//fill a given histogram with the sample values
void TWaveForm::FillHist( TH1I* histo) const
{
  for( int i=0; i<nb_samp; i++) 
    histo->Fill(samp_y[i]);
}





double TWaveForm::Distance10(const TWaveForm* wfm) const
{

  //  printf(" here 0\n");
  // double* x0 = new double[nb_samp];	
  // double* y0 = new double[nb_samp];

  double* x0; // = new double[nb_samp];	
  double* y0; // = new double[nb_samp];

  x0 = Get_samp_x();
  y0 = Get_samp_y();

  // adjust the base to 0
  if(0)
    for( int i = nb_samp-1; i >= 0; i--){
      y0[i] -= y0[0];
    }


  //  double* x1 = new double[wfm->nb_samp];	
  //  double* y1 = new double[wfm->nb_samp];

  double* x1; // = new double[wfm->nb_samp];	
  double* y1; // = new double[wfm->nb_samp];

  x1 = wfm->Get_samp_x();
  y1 = wfm->Get_samp_y();


  //  double factor = y0[nb_samp-1]/(y1[wfm->nb_samp-1]-y1[0]);

  // normalize the peak to the same value
  if(0)
    for( int i = wfm->nb_samp-1; i >= 0; i--){
      y1[i] -= y1[0];
      //      y1[i] *= factor;
    }


  //
  int    start_wfm = 1;
  int    y0_i = 1;
  double time_diff = 0;
  int    count = 0;
  //  double up_th = 600.;
  double up_th = 250.;
  double do_th = 200.;
  int count_w = 0;



  //  printf(" here 1\n");
  if( 1){
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > do_th){
	y0_i  = i;
	break;
      }
    }
    
    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > do_th){
	start_wfm  = i-1;
	break;
      }
    }
    
    //    printf(" cfd value = %7.3f\n", up_th*y0[nb_samp-1]);


    //    printf(" here 2\n");
    while( (y0_i != (nb_samp-1)) && (y0[y0_i] < up_th)  && count_w < 10000 ){
    //    while( (y0_i != (nb_samp-1)) && (y0[y0_i] < 120.)  && count_w < 10000 ){
      count_w++;
      for( int j = start_wfm; j < wfm->nb_samp; j++){
	if( y0[y0_i] < y1[j]){
	  // slop*distance + offset
	  double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	  
	  //	time_diff += fabs( proj_t - x0[y0_i]);
	  time_diff += ( proj_t - x0[y0_i]);
	  
	  y0_i++;
	  count++;
	  
	  //	  printf("proj_t = %8.3f, x0 = %8.3f,  x1 = %8.3f\n", proj_t, x0[y0_i], x1[j-1]);

	  break;
	}
	else
	  start_wfm++;
      }
      
      //    printf(" y0_i = %5d\n", y0_i);
    }
  }  


  //  printf(" here 3\n");
  //  time_diff += fabs( x0[0] - x1[0]);
  //  count++;
  //  time_diff += fabs( x0[nb_samp-1] - x1[wfm->nb_samp-1]);
  
  //  printf("count = %d\n", count);


  if( count_w >= 9999){
    printf(".... lost ... \n");
    return -1000.;
  }


  if( 0){
    int    start = 1;
    int    y1_i = 1;
    
    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > do_th){
	y1_i  = i;
	break;
      }
    }
    
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > do_th){
	start  = i-1;
	break;
      }
    }
    
    
    //  while( y1_i != (wfm->nb_samp-1)){
    while( (y1_i != (wfm->nb_samp-1)) && (y1[y1_i] < up_th*y1[wfm->nb_samp-1])  ){
      
      for( int j = start; j < nb_samp; j++){
	if( y1[y1_i] < y0[j]){
	  // slop*distance + offset
	  double proj_t = ((x0[j]-x0[j-1])/(y0[j]-y0[j-1]))*(y1[y1_i]-y0[j-1]) + x0[j-1];
	  
	  //	time_diff += fabs( proj_t - x1[y1_i]);
	  time_diff += ( proj_t - x1[y1_i]);
	  
	  y1_i++; 
	  count++;
	  break;
	}
	else
	  start++;
      }
      
      
      // printf(" nb_sampl = %5d, wfm->nb_samp  = %5d, y1_i = %5d, start = %5d,  %7.3f,  %7.3f\n", 
      //	   nb_samp, wfm->nb_samp, y1_i, start, y1[y1_i], y0[start-1]);
      
    }
  
    //  time_diff += fabs( x0[0] - x1[0]);
    //  time_diff += fabs( x0[nb_samp-1] - x1[wfm->nb_samp-1]);
    
    // printf("nsample = %5d,  y0_i = %5d,   time_diff = %7.3f\n", 
    // nb_samp, y0_i, time_diff/(nb_samp+wfm->nb_samp));
    
    
    // return  time_diff/(nb_samp);
    // return  time_diff/(nb_samp+wfm->nb_samp);  
  }

  
  return  time_diff/count;
}




double TWaveForm::Distance(const TWaveForm* wfm) const
{

  double* x0 = new double[nb_samp];	
  double* y0 = new double[nb_samp];

  x0 = Get_samp_x();
  y0 = Get_samp_y();

  // adjust the base to 0
  for( int i = nb_samp-1; i >= 0; i--){
    y0[i] -= y0[0];
  }


  double* x1 = new double[wfm->nb_samp];	
  double* y1 = new double[wfm->nb_samp];

  x1 = wfm->Get_samp_x();
  y1 = wfm->Get_samp_y();


  double factor = y0[nb_samp-1]/(y1[wfm->nb_samp-1]-y1[0]);

  // normalize the peak to the same value
  if(0)
  for( int i = wfm->nb_samp-1; i >= 0; i--){
    y1[i] -= y1[0];
    y1[i] *= factor;
  }


  //
  int    start_wfm = 1;
  int    y0_i = 1;
  double time_diff = 0;
  int    count = 0;
  double up_th = 0.3;
  double do_th = Amplitude()*0.2;
  int count_w = 0;

  if( 1){
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > do_th){
	y0_i  = i;
	break;
      }
    }
    
    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > do_th){
	start_wfm  = i-1;
	break;
      }
    }
    
    //    printf(" cfd value = %7.3f\n", up_th*y0[nb_samp-1]);


    while( (y0_i != (nb_samp-1)) && (y0[y0_i] < up_th*y0[nb_samp-1])  && count_w < 10000 ){
    //    while( (y0_i != (nb_samp-1)) && (y0[y0_i] < 120.)  && count_w < 10000 ){
      count_w++;
      for( int j = start_wfm; j < wfm->nb_samp; j++){
	if( y0[y0_i] < y1[j]){
	  // slop*distance + offset
	  double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	  
	  //	time_diff += fabs( proj_t - x0[y0_i]);
	  time_diff += ( proj_t - x0[y0_i]);
	  
	  y0_i++;
	  count++;
	  
	  break;
	}
	else
	  start_wfm++;
      }
      
      //    printf(" y0_i = %5d\n", y0_i);
    }
  }  

  //  time_diff += fabs( x0[0] - x1[0]);
  //  count++;
  //  time_diff += fabs( x0[nb_samp-1] - x1[wfm->nb_samp-1]);
  

  if( count_w >= 9999){
    printf(".... lost ... \n");
    return -1000.;
  }


  if( 0){
    int    start = 1;
    int    y1_i = 1;
    
    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > do_th){
	y1_i  = i;
	break;
      }
    }
    
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > do_th){
	start  = i-1;
	break;
      }
    }
    
    
    //  while( y1_i != (wfm->nb_samp-1)){
    while( (y1_i != (wfm->nb_samp-1)) && (y1[y1_i] < up_th*y1[wfm->nb_samp-1])  ){
      
      for( int j = start; j < nb_samp; j++){
	if( y1[y1_i] < y0[j]){
	  // slop*distance + offset
	  double proj_t = ((x0[j]-x0[j-1])/(y0[j]-y0[j-1]))*(y1[y1_i]-y0[j-1]) + x0[j-1];
	  
	  //	time_diff += fabs( proj_t - x1[y1_i]);
	  time_diff += ( proj_t - x1[y1_i]);
	  
	  y1_i++; 
	  count++;
	  break;
	}
	else
	  start++;
      }
      
      
      // printf(" nb_sampl = %5d, wfm->nb_samp  = %5d, y1_i = %5d, start = %5d,  %7.3f,  %7.3f\n", 
      //	   nb_samp, wfm->nb_samp, y1_i, start, y1[y1_i], y0[start-1]);
      
    }
  
    //  time_diff += fabs( x0[0] - x1[0]);
    //  time_diff += fabs( x0[nb_samp-1] - x1[wfm->nb_samp-1]);
    
    // printf("nsample = %5d,  y0_i = %5d,   time_diff = %7.3f\n", 
    // nb_samp, y0_i, time_diff/(nb_samp+wfm->nb_samp));
    
    
    // return  time_diff/(nb_samp);
    // return  time_diff/(nb_samp+wfm->nb_samp);  
  }



  return  time_diff/count;
}


double TWaveForm::Chisquare(const TWaveForm* wfm, double th_low, double th_high, int& count, double* pdT) const
//double TWaveForm::Chisquare(const TWaveForm* wfm, double th_low, double th_high, int& count) const
{

  double* x0 = new double[nb_samp];	
  double* y0 = new double[nb_samp];

  x0 = Get_samp_x();
  y0 = Get_samp_y();


  double* x1 = new double[wfm->nb_samp];	
  double* y1 = new double[wfm->nb_samp];

  x1 = wfm->Get_samp_x();
  y1 = wfm->Get_samp_y();


  // decide the range using thresholds
  int i0, i1;
  i0 = i1 = 0;

  for( int i = 0; i < nb_samp; i++)
    if(y0[i] > th_low){
      i0 = i;
      break;
    }

  for( int i = 0; i < nb_samp; i++)
    if(y0[i] > th_high){
      i1 = i;
      break;
    }
  
  //
  double dt = 0;
  count = 0;
  for( int i = i0; i < i1; i++){

    int j0, j1;
    j0 = j1 = 0;
    for( int j = 0; j < nb_samp; j++)
      if( y1[j] > y0[i]){
	j1 = j;
	j0 = j-1;
	break;
      }
    
    double slope = (y1[j1]-y1[j0])/(x1[j1]-x1[j0]);
    double proj_t = x1[j0] + (y0[i]-y1[j0])/slope;

    pdT[count] =  (x0[i] - proj_t);
    dt += (x0[i] - proj_t);
    count++;
  }

  return  dt/count;
}


double TWaveForm::Distance(const TWaveForm* wfm, TH1F* histo) const
{
  double* x0 = new double[nb_samp];	
  double* y0 = new double[nb_samp];

  x0 = Get_samp_x();
  y0 = Get_samp_y();

  // adjust the base to 0
  for( int i = nb_samp-1; i >= 0; i--){
    y0[i] -= y0[0];
  }


  double* x1 = new double[wfm->nb_samp];	
  double* y1 = new double[wfm->nb_samp];

  x1 = wfm->Get_samp_x();
  y1 = wfm->Get_samp_y();


  double factor = y0[nb_samp-1]/(y1[wfm->nb_samp-1]-y1[0]);

  // normalize the peak to the same value
  if(1)
  for( int i = wfm->nb_samp-1; i >= 0; i--){
    y1[i] -= y1[0];
    y1[i] *= factor;
  }

  //
  int    start_wfm = 1;
  int    y0_i = 1;
  double time_diff = 0;
  int    count = 0;
  double up_th = 1.0;
  
  for( int i = 0; i < nb_samp; i++){
    if( y0[i] > 20){
      y0_i  = i;
      break;
    }
  }

  for( int i = 0; i < wfm->nb_samp; i++){
    if( y1[i] > 20){
      start_wfm  = i-1;
      break;
    }
  }

  //  while( y0_i != (nb_samp-1)){    
  while( (y0_i != (nb_samp-1)) && (y0[y0_i] < up_th*y0[nb_samp-1])  ){

    for( int j = start_wfm; j < wfm->nb_samp; j++){

      if( y0[y0_i] < y1[j]){
	// slop*distance + offset
	double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	
	//	time_diff += fabs( proj_t - x0[y0_i]);
	time_diff += ( proj_t - x0[y0_i]);
	y0_i++; 
	count++;

	//	histo->Fill(  fabs( proj_t - x0[y0_i]) );
	histo->Fill(  ( proj_t - x0[y0_i]) );

	break;
      }
      else
	start_wfm++;
    }

    //    printf(" y0_i = %5d\n", y0_i);
  }


   
  // time_diff += fabs( x0[0] - x1[0]);
  // time_diff += fabs( x0[nb_samp-1] - x1[wfm->nb_samp-1]);
  
  if(0){
    int    start = 1;
    int    y1_i = 1;

    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > 20){
	y1_i  = i;
      break;
      }
    }
    
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > 20){
	start  = i-1;
	break;
      }
    }
    
    
    //  while( y1_i != (wfm->nb_samp-1)){
    while( (y1_i != (wfm->nb_samp-1)) && (y1[y1_i] <  up_th*y1[wfm->nb_samp-1])  ){    
      for( int j = start; j < nb_samp; j++){
	if( y1[y1_i] < y0[j]){
	  // slop*distance + offset
	  double proj_t = ((x0[j]-x0[j-1])/(y0[j]-y0[j-1]))*(y1[y1_i]-y0[j-1]) + x0[j-1];
	  
	  //	time_diff += fabs( proj_t - x1[y1_i]);
	  //	histo->Fill(  fabs( proj_t - x1[y1_i]) );
	  
	  time_diff += ( proj_t - x1[y1_i]);
	  histo->Fill(  ( proj_t - x1[y1_i]) );
	  
	  y1_i++;
	  count++;
	  break;
	}
	else
	  start++;
      }
      
      
      // printf(" nb_sampl = %5d, wfm->nb_samp  = %5d, y1_i = %5d, start = %5d,  %7.3f,  %7.3f\n", 
      // nb_samp, wfm->nb_samp, y1_i, start, y1[y1_i], y0[start-1]);
      
    }
    
    // time_diff += fabs( x0[0] - x1[0]);
    // time_diff += fabs( x0[nb_samp-1] - x1[wfm->nb_samp-1]);
    
    // printf("nsample = %5d,  y0_i = %5d,   time_diff = %7.3f\n", 
    // nb_samp, y0_i, time_diff/(nb_samp+wfm->nb_samp));
    
    
    // return  time_diff/(nb_samp);
    //  return  time_diff/(nb_samp+wfm->nb_samp);  
  }


  return  time_diff/count;
}


double TWaveForm::Distance2(const TWaveForm* wfm) const
{

  double* x0 = new double[nb_samp];	
  double* y0 = new double[nb_samp];

  x0 = Get_samp_x();
  y0 = Get_samp_y();

  // adjust the base to 0
  for( int i = nb_samp-1; i >= 0; i--){
    y0[i] -= y0[0];
  }


  double* x1 = new double[wfm->nb_samp];	
  double* y1 = new double[wfm->nb_samp];

  x1 = wfm->Get_samp_x();
  y1 = wfm->Get_samp_y();

  int max0, max1;
  double amp0, amp1;
  max0 = max1 = 0;
  amp0 = amp1 = 0;

  for( int i = 0; i < nb_samp; i++){
    if( y0[i] > amp0){
      amp0 = y0[i];
      max0 = i;
    }
  }

  for( int i = 0; i < wfm->nb_samp; i++){
    if( y1[i] > amp1){
      amp1 = y1[i];
      max1 = i;
    }
  }

  double factor = amp0/(amp1-y1[0]);

  // normalize the peak to the same value
  for( int i = wfm->nb_samp-1; i >= 0; i--){
    y1[i] -= y1[0];
    y1[i] *= factor;
  }

  //
  int    start_wfm = 1;
  int    y0_i = 1;
  double time_diff = 0;
  int    count = 0;
  double up_th = 1.0;
  //  double do_th = 0.2;
  int    count_w = 0;

  if( 1){

    // start from > 20mV
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > 20){
	y0_i  = i;
	break;
      }
    }

    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > 20){
	start_wfm  = i-1;
	break;
      }
    }

    // rising slope 
    while( (y0_i != (max0-1)) && (y0[y0_i] < up_th*amp0) && count_w < 10000){
      count_w++;

      for( int j = start_wfm; j < max1+1; j++){
	if( y0[y0_i] < y1[j]){
	  // slop*distance + offset
	  double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	
	  //	  time_diff += fabs( proj_t - x0[y0_i]);
	  time_diff += ( proj_t - x0[y0_i]);
	  y0_i++;
	  count++;

	  break;
	}
	else
	  start_wfm++;
      }

      // printf(" y0_i = %5d\n", y0_i);
    }



    y0_i = max0+1;
    start_wfm = max1+1;

    // normalization of falling part 
    double factor2 = (amp0-y0[nb_samp-1])/(amp0-y1[wfm->nb_samp-1]);
    
    if(1)
    for( int i = max1+1; i < wfm->nb_samp; i++){
      y1[i] = amp0-(amp0-y1[i])*factor2;
    }

    if( count_w >= 10000) return -1000.;
    count_w = 0;
    
    /*
    for( int i = 0; i < 3; i++){
      printf("After : i = %5d => %7.3f,  %7.3f\n", i, y0[nb_samp-1-i], y1[wfm->nb_samp-1-i]);
    }

    for( int i = 0; i < 3; i++){
      printf("After2 : i = %5d => %7.3f,  %7.3f\n", i, y0[max0+i], y1[max1+i]);
    }
    */

    // falling slope 
    while( (y0_i != (nb_samp-1)) && (y0[y0_i] > 20) && count_w < 10000){

      count_w++;

      for( int j = start_wfm; j < wfm->nb_samp; j++){
	if( y0[y0_i] > y1[j]){
	  // slop*distance + offset
	  double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	
	  //	  time_diff += fabs( proj_t - x0[y0_i]);
	  time_diff += ( proj_t - x0[y0_i]);
	  y0_i++;
	  count++;

	  break;
	}
	else
	  start_wfm++;
      }
      //    printf(" y0_i = %5d\n", y0_i);
    }

    if( count_w >= 10000) return -1000.;
  } 


  //  time_diff += fabs( x0[0] - x1[0]);
  //  time_diff += fabs( x0[max0] - x1[max1]);
  //  count++;
  

  // the other projection
  if( 0){
    int    start = 1;
    int    y1_i = 1;
    
    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > 20){
	y1_i  = i;
	break;
      }
    }
    
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > 20){
	start  = i-1;
	break;
      }
    }

  
    // rising slope 
    while( (y1_i != (max1-1)) && (y1[y1_i] < up_th*amp0)  ){
    
      for( int j = start; j < max0+1; j++){
	if( y1[y1_i] < y0[j]){
	  // slop*distance + offset
	  double proj_t = ((x0[j]-x0[j-1])/(y0[j]-y0[j-1]))*(y1[y1_i]-y0[j-1]) + x0[j-1];
	
	  //	  time_diff += fabs( proj_t - x1[y1_i]);
	  time_diff += ( proj_t - x1[y1_i]);
	  y1_i++;
	  count++;

	  break;
	}
	else
	  start++;
      }
      //    printf(" y0_i = %5d\n", y0_i);
    }

    y1_i = max1+1;
    start = max0+1;


    // falling slope 
    while( (y1_i != (wfm->nb_samp-1)) && (y1[y1_i] > 20)  ){
    
      for( int j = start; j < nb_samp; j++){
	if( y1[y1_i] > y0[j]){
	  // slop*distance + offset
	  double proj_t = ((x0[j]-x0[j-1])/(y0[j]-y0[j-1]))*(y1[y1_i]-y0[j-1]) + x0[j-1];
	
	  //	  time_diff += fabs( proj_t - x1[y1_i]);
	  time_diff += ( proj_t - x1[y1_i]);
	  y1_i++;
	  count++;

	  break;
	}
	else
	  start++;
      }
      //    printf(" y0_i = %5d\n", y0_i);
    }
    
  }



  return  time_diff/count;
}



double TWaveForm::Distance2(const TWaveForm* wfm, TH1F* histo) const
{

  double* x0 = new double[nb_samp];	
  double* y0 = new double[nb_samp];

  x0 = Get_samp_x();
  y0 = Get_samp_y();

  // adjust the base to 0
  for( int i = nb_samp-1; i >= 0; i--){
    y0[i] -= y0[0];
  }


  double* x1 = new double[wfm->nb_samp];	
  double* y1 = new double[wfm->nb_samp];

  x1 = wfm->Get_samp_x();
  y1 = wfm->Get_samp_y();

  int max0, max1;
  double amp0, amp1;
  max0 = max1 = 0;
  amp0 = amp1 = 0;

  for( int i = 0; i < nb_samp; i++){
    if( y0[i] > amp0){
      amp0 = y0[i];
      max0 = i;
    }
  }

  for( int i = 0; i < wfm->nb_samp; i++){
    if( y1[i] > amp1){
      amp1 = y1[i];
      max1 = i;
    }
  }

  double factor = amp0/(amp1-y1[0]);

  // normalize the peak to the same value
  for( int i = wfm->nb_samp-1; i >= 0; i--){
    y1[i] -= y1[0];
    y1[i] *= factor;
  }

  //
  int    start_wfm = 1;
  int    y0_i = 1;
  double time_diff = 0;
  int    count = 0;
  double up_th = 1.0;
  //  double do_th = 0.2;
  int count_w = 0;


  if( 1){

    // start from > 20mV
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > 20){
	y0_i  = i;
	break;
      }
    }

    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > 20){
	start_wfm  = i-1;
	break;
      }
    }

    // rising slope 
    while( (y0_i != (max0-1)) && (y0[y0_i] < up_th*amp0) && count_w < 10000 ){
      count_w++;

      for( int j = start_wfm; j < max1+1; j++){
	if( y0[y0_i] < y1[j]){
	  // slop*distance + offset
	  double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	
	  //	  time_diff += fabs( proj_t - x0[y0_i]);
	  //	  histo->Fill( fabs( proj_t - x0[y0_i]));

	  time_diff += ( proj_t - x0[y0_i]);
	  histo->Fill( ( proj_t - x0[y0_i]));

	  y0_i++;
	  count++;

	  break;
	}
	else
	  start_wfm++;
      }
      //    printf(" y0_i = %5d\n", y0_i);
    }

    y0_i = max0+1;
    start_wfm = max1+1;

    // normalization of falling part 
    double factor2 = (amp0-y0[nb_samp-1])/(amp0-y1[wfm->nb_samp-1]);
    
    if(1)
    for( int i = max1+1; i < wfm->nb_samp; i++){
      y1[i] = amp0-(amp0-y1[i])*factor2;
    }


    if( count_w >= 10000) return -1000.;
    count_w = 0;
    

    // falling slope 
    while( (y0_i != (nb_samp-1)) && (y0[y0_i] > 20) && count_w < 10000 ){
      count_w++;

      for( int j = start_wfm; j < wfm->nb_samp; j++){
	if( y0[y0_i] > y1[j]){
	  // slop*distance + offset
	  double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	
	  //	  time_diff += fabs( proj_t - x0[y0_i]);
	  //	  histo->Fill( fabs( proj_t - x0[y0_i]));

	  time_diff += ( proj_t - x0[y0_i]);
	  histo->Fill( ( proj_t - x0[y0_i]));

	  y0_i++;
	  count++;

	  break;
	}
	else
	  start_wfm++;
      }

      //    printf(" y0_i = %5d\n", y0_i);
    }
    
    if( count_w >= 10000) return -1000.;
  }  

  //  time_diff += fabs( x0[0] - x1[0]);
  //  count++;
  //  time_diff += fabs( x0[nb_samp-1] - x1[wfm->nb_samp-1]);
  


  // the other projection
  if( 0){
    int    start = 1;
    int    y1_i = 1;
    
    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > 20){
	y1_i  = i;
	break;
      }
    }
    
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > 20){
	start  = i-1;
	break;
      }
    }

  
    // rising slope 
    while( (y1_i != (max1-1)) && (y1[y1_i] < up_th*amp0)  ){
    
      for( int j = start; j < max0+1; j++){
	if( y1[y1_i] < y0[j]){
	  // slop*distance + offset
	  double proj_t = ((x0[j]-x0[j-1])/(y0[j]-y0[j-1]))*(y1[y1_i]-y0[j-1]) + x0[j-1];
	
	  //	  time_diff += fabs( proj_t - x1[y1_i]);
	  //	  histo->Fill( fabs( proj_t - x1[y1_i]));

	  time_diff += ( proj_t - x1[y1_i]);
	  histo->Fill( ( proj_t - x1[y1_i]));

	  y1_i++;
	  count++;

	  if(  fabs( proj_t - x1[y1_i]) > 2.0)	    
	    printf("falling : y1_i = %5d, y0_i = %5d\n", y1_i, j);

	  break;
	}
	else
	  start++;
      }
      //    printf(" y0_i = %5d\n", y0_i);
    }

    y1_i = max1+1;
    start = max0+1;


    // falling slope 
    while( (y1_i != (wfm->nb_samp-1)) && (y1[y1_i] > 20)  ){
    
      for( int j = start; j < nb_samp; j++){
	if( y1[y1_i] > y0[j]){
	  // slop*distance + offset
	  double proj_t = ((x0[j]-x0[j-1])/(y0[j]-y0[j-1]))*(y1[y1_i]-y0[j-1]) + x0[j-1];
	
	  //	  time_diff += fabs( proj_t - x1[y1_i]);
	  //	  histo->Fill( fabs( proj_t - x1[y1_i]));

	  time_diff += ( proj_t - x1[y1_i]);
	  histo->Fill( ( proj_t - x1[y1_i]));

	  if(  fabs( proj_t - x1[y1_i]) > 2.0)	    
	    printf("falling : y1_i = %5d, y0_i = %5d\n", y1_i, j);


	  y1_i++;
	  count++;

	  break;
	}
	else
	  start++;
      }
      //    printf(" y0_i = %5d\n", y0_i);
    }
    
  }



  return  time_diff/count;
}



double TWaveForm::Distance3(const TWaveForm* wfm) const
{

  double* x0 = new double[nb_samp];	
  double* y0 = new double[nb_samp];

  x0 = Get_samp_x();
  y0 = Get_samp_y();

  // adjust the base to 0
  for( int i = nb_samp-1; i >= 0; i--){
    y0[i] -= y0[0];
  }


  double* x1 = new double[wfm->nb_samp];	
  double* y1 = new double[wfm->nb_samp];

  x1 = wfm->Get_samp_x();
  y1 = wfm->Get_samp_y();

  int max0, max1;
  double amp0, amp1;
  max0 = max1 = 0;
  amp0 = amp1 = 0;

  for( int i = 0; i < nb_samp; i++){
    if( y0[i] > amp0){
      amp0 = y0[i];
      max0 = i;
    }
  }

  for( int i = 0; i < wfm->nb_samp; i++){
    if( y1[i] > amp1){
      amp1 = y1[i];
      max1 = i;
    }
  }

  double factor = amp0/(amp1-y1[0]);

  // normalize the peak to the same value
  for( int i = wfm->nb_samp-1; i >= 0; i--){
    y1[i] -= y1[0];
    y1[i] *= factor;
  }

  //
  int    start_wfm = 1;
  int    y0_i = 1;
  double time_diff = 0;
  int    count = 0;
  double up_th = 1.0;
  //  double do_th = 0.2;
  int    count_w = 0;

  vector<double> numbers;

  
  if( 1){

    // start from > 20mV
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > 20){
	y0_i  = i;
	break;
      }
    }

    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > 20){
	start_wfm  = i-1;
	break;
      }
    }

    // rising slope 
    while( (y0_i != (max0-1)) && (y0[y0_i] < up_th*amp0) && count_w < 10000){
      count_w++;

      for( int j = start_wfm; j < max1+1; j++){
	if( y0[y0_i] < y1[j]){
	  // slop*distance + offset
	  double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	
	  //	  time_diff += fabs( proj_t - x0[y0_i]);
	  time_diff += ( proj_t - x0[y0_i]);

	  double dtmp = ( proj_t - x0[y0_i]);
	  numbers.push_back( dtmp);

	  // *numbers = ( proj_t - x0[y0_i]);
	  //  numbers++;

	  y0_i++;
	  count++;

	  break;
	}
	else
	  start_wfm++;
      }

      // printf(" y0_i = %5d\n", y0_i);
    }


    y0_i = max0+1;
    start_wfm = max1+1;

    // normalization of falling part 
    double factor2 = (amp0-y0[nb_samp-1])/(amp0-y1[wfm->nb_samp-1]);
    
    if(1)
    for( int i = max1+1; i < wfm->nb_samp; i++){
      y1[i] = amp0-(amp0-y1[i])*factor2;
    }

    if( count_w >= 10000) return -1000.;
    count_w = 0;
    

    // falling slope 
    while( (y0_i != (nb_samp-1)) && (y0[y0_i] > 20) && count_w < 10000){

      count_w++;

      for( int j = start_wfm; j < wfm->nb_samp; j++){
	if( y0[y0_i] > y1[j]){
	  // slop*distance + offset
	  double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	
	  //	  time_diff += fabs( proj_t - x0[y0_i]);
	  time_diff += ( proj_t - x0[y0_i]);

	  double dtmp = ( proj_t - x0[y0_i]);
	  numbers.push_back( dtmp);

	  y0_i++;
	  count++;

	  break;
	}
	else
	  start_wfm++;
      }
      //    printf(" y0_i = %5d\n", y0_i);
    }

    if( count_w >= 10000) return -1000.;

    sort( numbers.begin(), numbers.end());
    

    int bin = 5;

    double step = (numbers[numbers.size()-1]-numbers[0])/bin;
    double nstep[bin], sum[bin];
    int nentry[bin];

    for( int i = 0; i < bin; i++){
      nentry[i] = 0;
      sum[i] = 0;
      nstep[i] = step*(i+1) + numbers[0]; 
      //      printf(" i = %5d,  nstep = %10.6f\n", i, nstep[i]);
    }

    int start_v = 0;

    for( uint i = 0; i < numbers.size(); i++){
     
      if( numbers[i] > nstep[start_v])
	start_v++;
      
      nentry[start_v]++;
      sum[start_v] += numbers[i];
  
    }

    time_diff = 0;
    int max_ent = 0;

    for( int i = 0; i < bin; i++){
      //      printf("i = %5d,   nentry = %5d,   sum = %10.6f ( %10.6f)\n", i, nentry[i], sum[i], sum[i]/nentry[i]);
      if( nentry[i] > max_ent){
	time_diff = sum[i]/nentry[i];
	max_ent = nentry[i];
      }
    }
    
  } 


  //  time_diff += fabs( x0[0] - x1[0]);
  //  time_diff += fabs( x0[max0] - x1[max1]);
  //  count++;
  
  return  time_diff;
  //  return  time_diff/count;
}



double TWaveForm::Distance3(const TWaveForm* wfm, TH1F* histo) const
{

  double* x0 = new double[nb_samp];	
  double* y0 = new double[nb_samp];

  x0 = Get_samp_x();
  y0 = Get_samp_y();

  // adjust the base to 0
  for( int i = nb_samp-1; i >= 0; i--){
    y0[i] -= y0[0];
  }


  double* x1 = new double[wfm->nb_samp];	
  double* y1 = new double[wfm->nb_samp];

  x1 = wfm->Get_samp_x();
  y1 = wfm->Get_samp_y();

  int max0, max1;
  double amp0, amp1;
  max0 = max1 = 0;
  amp0 = amp1 = 0;

  for( int i = 0; i < nb_samp; i++){
    if( y0[i] > amp0){
      amp0 = y0[i];
      max0 = i;
    }
  }

  for( int i = 0; i < wfm->nb_samp; i++){
    if( y1[i] > amp1){
      amp1 = y1[i];
      max1 = i;
    }
  }

  double factor = amp0/(amp1-y1[0]);

  // normalize the peak to the same value
  for( int i = wfm->nb_samp-1; i >= 0; i--){
    y1[i] -= y1[0];
    y1[i] *= factor;
  }

  //
  int    start_wfm = 1;
  int    y0_i = 1;
  double time_diff = 0;
  int    count = 0;
  double up_th = 1.0;
  //  double do_th = 0.2;
  int    count_w = 0;

  vector<double> numbers;

  
  if( 1){

    // start from > 20mV
    for( int i = 0; i < nb_samp; i++){
      if( y0[i] > 20){
	y0_i  = i;
	break;
      }
    }

    for( int i = 0; i < wfm->nb_samp; i++){
      if( y1[i] > 20){
	start_wfm  = i-1;
	break;
      }
    }

    // rising slope 
    while( (y0_i != (max0-1)) && (y0[y0_i] < up_th*amp0) && count_w < 10000){
      count_w++;

      for( int j = start_wfm; j < max1+1; j++){
	if( y0[y0_i] < y1[j]){
	  // slop*distance + offset
	  double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	
	  //	  time_diff += fabs( proj_t - x0[y0_i]);
	  time_diff += ( proj_t - x0[y0_i]);
	  histo->Fill( proj_t - x0[y0_i]);

	  double dtmp = ( proj_t - x0[y0_i]);
	  numbers.push_back( dtmp);

	  // *numbers = ( proj_t - x0[y0_i]);
	  //  numbers++;

	  y0_i++;
	  count++;

	  break;
	}
	else
	  start_wfm++;
      }

      // printf(" y0_i = %5d\n", y0_i);
    }


    y0_i = max0+1;
    start_wfm = max1+1;

    // normalization of falling part 
    double factor2 = (amp0-y0[nb_samp-1])/(amp0-y1[wfm->nb_samp-1]);
    
    if(1)
    for( int i = max1+1; i < wfm->nb_samp; i++){
      y1[i] = amp0-(amp0-y1[i])*factor2;
    }

    if( count_w >= 10000) return -1000.;
    count_w = 0;
    

    // falling slope 
    while( (y0_i != (nb_samp-1)) && (y0[y0_i] > 20) && count_w < 10000){

      count_w++;

      for( int j = start_wfm; j < wfm->nb_samp; j++){
	if( y0[y0_i] > y1[j]){
	  // slop*distance + offset
	  double proj_t = ((x1[j]-x1[j-1])/(y1[j]-y1[j-1]))*(y0[y0_i]-y1[j-1]) + x1[j-1];
	
	  //	  time_diff += fabs( proj_t - x0[y0_i]);
	  time_diff += ( proj_t - x0[y0_i]);
	  histo->Fill( proj_t - x0[y0_i]);

	  double dtmp = ( proj_t - x0[y0_i]);
	  numbers.push_back( dtmp);

	  y0_i++;
	  count++;

	  break;
	}
	else
	  start_wfm++;
      }
      //    printf(" y0_i = %5d\n", y0_i);
    }

    if( count_w >= 10000) return -1000.;

    sort( numbers.begin(), numbers.end());
    
    /*
    for( int i = 0; i < numbers.size(); i++){
      printf(" %5d   =>  %10.6f\n", i, numbers[i]);
    }

    printf("min = %10.6f,  max = %10.6f\n", numbers[0], numbers[numbers.size()-1]);
    */


    int bin = 5;

    double step = (numbers[numbers.size()-1]-numbers[0])/bin;

    double nstep[bin], sum[bin];
    int nentry[bin];

    for( int i = 0; i < bin; i++){
      nentry[i] = 0;
      sum[i] = 0;
      nstep[i] = step*(i+1) + numbers[0]; 
      //      printf(" i = %5d,  nstep = %10.6f\n", i, nstep[i]);
    }

    int start_v = 0;

    for( uint i = 0; i < numbers.size(); i++){
     
      if( numbers[i] > nstep[start_v])
	start_v++;
      
      nentry[start_v]++;
      sum[start_v] += numbers[i];
  
    }
    
    time_diff = 0;
    int max_ent = 0;

    for( int i = 0; i < bin; i++){
      printf("i = %5d,   nentry = %5d,   sum = %10.6f ( %10.6f)\n", i, nentry[i], sum[i], sum[i]/nentry[i]);
      if( nentry[i] > max_ent){
	time_diff = sum[i]/nentry[i];
	max_ent = nentry[i];
      }
    }

  } 





  //  time_diff += fabs( x0[0] - x1[0]);
  //  time_diff += fabs( x0[max0] - x1[max1]);
  //  count++;
  

  return  time_diff;
  //  return  time_diff/count;
}

