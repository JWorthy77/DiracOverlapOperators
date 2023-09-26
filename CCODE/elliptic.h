#ifndef _ELLIPTIC
#define _ELLIPTIC

#include <boost/multiprecision/mpfr.hpp>
namespace mp = boost::multiprecision;
typedef mp::number<mp::mpfr_float_backend<100> > ftype;
static ftype err=1e-75;

class elliptic {

public:
  ftype getpi(ftype err);
  ftype calcK(ftype k, ftype err);
  ftype agm(ftype x, ftype y, ftype err);
  ftype calcam(ftype u, ftype m, ftype err); 
  ftype sn(ftype u, ftype m, ftype err); 
  ftype dn(ftype u, ftype m, ftype err); 
  ftype asnfn(ftype x, ftype k); 
  ftype arcsnN(ftype xmax, ftype k, int N); 
  ftype arcsn(ftype xmax, ftype k, ftype err);
};

#endif
