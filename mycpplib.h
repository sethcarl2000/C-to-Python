#ifndef MYCPPLIB_H
#define MYCPPLIB_H
#include <cmath>


class Mandel{
 public:
  Mandel(double re1=-2.1, double re2=0.6, double im1=-1.5, double im2=1.5);
  int mandel_test(double c_re, double c_im);
  void setResolution(int nr, int ni) {
    _nr=nr;
    _ni=ni;
  }
  void getResolution(int &nr, int &ni) {
    nr=_nr;
    ni=_ni;
  }
  void calculate();
  void setPrecision(int nTrials) {_nTrials=nTrials;}
  int getPrecision() {return _nTrials;}
  void getRange(double &re1, double &re2, double &im1, double &im2){
    re1=_re1;
    re2=_re2;
    im1=_im1;
    im2=_im2;
  }
  void setRange(double re1, double re2, double im1, double im2){
    _re1=re1;
    _re2=re2;
    _im1=im1;
    _im2=im2;
  }
  void getImg(double *img) {
    for (int i=0; i<_nr*_ni; ++i) img[i]=_img[i];
  }
  ~Mandel();
 private:
  double _re1, _re2, _im1, _im2;
  int _nTrials;
  int _nr, _ni;
  double *_img;  // image buffer
};

#ifdef EXTERNC
extern "C" {
#endif 

//return the volume of a d-dimensional hypersphere using 'N' points
double HSVolume(int d, unsigned long long N, double r); 

#ifdef EXTERNC
}
#endif

#endif // block to ensure header is only compiled once
