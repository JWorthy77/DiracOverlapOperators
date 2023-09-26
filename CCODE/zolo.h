#ifndef _ZOLO
#define _ZOLO

#include "elliptic.h"
//#include <boost/multiprecision/mpfr.hpp>
//namespace mp = boost::multiprecision;
//typedef mp::number<mp::mpfr_float_backend<100> > ftype;

class zolo {

public:
  zolo();
  zolo(ftype lmin, ftype lmax, int N);
  void setZolo(ftype lmin, ftype lmax, int N);
  void setInnerCoeffs();
  ftype calcM();
  ftype calcC();
  ftype calcLambda();
  ftype Rhalf(ftype& x);
  ftype zsgn(ftype& x);
  void getRoots();
  void getExtrema();
  void writeZolo(int Npts);
  void testZolo(int Nz);
  void evalZolo(int Np,std::ofstream& ofile);

  int Ndw; // = Nnum+Ndenom
  int Nnum; // numerator
  int Ndenom; // denominator
  ftype lmin;
  ftype lmax;
  ftype beta;
  ftype k;
  ftype kdash;
  ftype K;
  ftype Kdash;
  ftype M;
  ftype C;
  ftype lambda;
  ftype mult; // = m in notes
  std::vector<ftype> cnum;
  std::vector<ftype> cdenom;
  std::vector<ftype> roots;
  std::vector<ftype> extrema;

  bool VBS=false;
};

#endif
