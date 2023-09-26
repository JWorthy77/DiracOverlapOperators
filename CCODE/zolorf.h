#ifndef _ZOLORF
#define _ZOLORF

#include "zolo.h"

class zolorf {

public:
  zolorf();
  zolorf(zolo& z);
  void setFactorZolo(zolo& z);


  ftype evalFactorZolo(ftype& r);
  void evalFactorZolo(int Np);
  void writeFactors();

  void setPartialFractionZolo();
  void setPartialFractionZolo(std::vector<ftype>& num, std::vector<ftype>& denom,
                              std::vector<ftype>& front, std::vector<ftype>& pfcoeffs);
  void expandPoly(std::vector<ftype>& zeros,std::vector<ftype>& epoly);
  void dividePoly(std::vector<ftype>& enmr, std::vector<ftype>& edenom, 
                  std::vector<ftype>& front, std::vector<ftype>& dnum);
  ftype evalpoly(ftype& x,std::vector<ftype>& coeffs);
  ftype evalfpoly(ftype& x,std::vector<ftype>& zeros);
  ftype polyDeriv(ftype& x,std::vector<ftype>& zeros);

  ftype evalPartialFractionZolo(ftype& r);
  void evalZolo(int Np,bool PF);
  void writePartFracCoeffs();

  void testDividePoly();
  void testPartFrac();

private:
// factored coefficients
  int Nnum;
  int Ndenom;
  ftype rmin;
  ftype rmax;
  ftype d; // multiplication constant
  std::vector<ftype> av; // numerator coeffs (x-av[i])
  std::vector<ftype> dv; // denominator coeffs (x-dv[i])
// partial fraction coefficients
  int Nfront;
  std::vector<ftype> front; // polynomial
  std::vector<ftype> pfcoeffs; // partial fraction coeffs
  bool VERBOSE=true; 
};

#endif
