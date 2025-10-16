#include "mycpplib.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

// constructor
Mandel::Mandel(double re1, double re2, double im1, double im2) :
  _re1(re1), _re2(re2), _im1(im1), _im2(im2) {
  _nTrials=255;
  _nr=2000;        // pixel resolution
  _ni=2000;
  _img=0;
}



int Mandel::mandel_test(double c_re, double c_im){  
  // If a point is in the set, its magnitude will remain bounded by
  // 2.0 over iterations of z -> z^2 + C.  Stop the loop after a 
  // maximum of NTRIALS and consider this point to be in the set

  double z_re = c_re;
  double z_im = c_im;
  int counts = 1;
  while (z_re*z_re+z_im*z_im<=4.0  && counts<_nTrials) {
    counts++;
    double re=z_re;  // careful! keep old values for calculation
    double im=z_im;
    // z -> z^2
    z_re = re*re-im*im;
    z_im = 2*re*im;
    // add c to z^2
    z_re = z_re + c_re;
    z_im = z_im + c_im;
  }
  return counts;
}


// explore the Mandelbrot set
// we pass the image buffer as a 1D array, and access the pixels
// using pointer arithmatic
// eg for an array a[n][m], n=nrow, m=mcolumns
// a[0] = a[0][0]
// a[m] = a[1][0]
// a[m+1] = a[1][1]
// ...
void Mandel::calculate(){
  if (_img) delete [] _img;  // free any existing buffer
  _img = new double[_nr*_ni];
  double dx=(_re2-_re1)/_nr;
  double dy=(_im2-_im1)/_ni;
  // loop over grid starting in lower left corner
  for (int j=0; j<_ni; ++j){
    double im=_im1+j*dy;
    for (int i=0; i<_nr; ++i){
      double re=_re1+i*dx;
      _img[j*_nr+i]=mandel_test(re,im);
    }
  }
}

Mandel::~Mandel(){
  if (_img) delete [] _img;
}

#ifdef EXTERNC
extern "C" {
#endif

//compute the volume of a d-dimensional hypersphere by the 'stone-throwing' method. 
double HSVolume(int d, unsigned long long N, double r)
{
  srand48((long)time(NULL));

  unsigned long long count=0;
  for (unsigned long long i=0; i<N; i++) {

    //initialize the X-vector with pseudo-random vals
    double square_radius=0.; 
    for (int j=0; j<d; j++) {
      double x_j = drand48();
      square_radius += x_j*x_j; 
    } 

    if (square_radius < 1.) count++; 
  }

  return pow(2., d)*((double)count)/((double)N); 
}

#ifdef EXTERNC
}
#endif

