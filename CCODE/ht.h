#ifndef _HT
#define _HT

#include "elliptic.h"
//#include <boost/multiprecision/mpfr.hpp>
//namespace mp = boost::multiprecision;
//typedef mp::number<mp::mpfr_float_backend<100> > ftype;

class ht {

public:
  ht();
  ht(int N);
  ht(int N, ftype r1, ftype r2);
  void setHT(int N);
  ftype htsgn(ftype& x);
  void evalHT(int Np,std::ofstream& ofile, ftype r1, ftype r2, int scale,bool REFLECT);
  void writeHT(int Npts,ftype r1, ftype r2, int scale=1,bool REFLECT=false);
  void writeFactors();

  int Ndw; // = Nnum+Ndenom
  int Nnum; // numerator
  int Ndenom; // denominator
  std::vector<ftype> num; // zeros of numerator
  std::vector<ftype> denom; // zeros of denominator
  ftype mult; // multiplier

  bool VBS=false;
};

#endif
