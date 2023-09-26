#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <gmp.h>
#include "math.h"

#include "ht.h"
#include "elliptic.h"
#include "zolo.h"
#include "zolorf.h"

int main() {

  bool ZOLO=false;
  if (ZOLO) { 
  elliptic e;
  zolo z;
  z.setZolo(ftype(1e-3),ftype(10),21);
  z.writeZolo(200);
//  z.testZolo(18);

  zolorf zrf(z);
//  zolorf zrf;
//  zrf.testDividePoly();
//  zrf.testPartFrac();
  zrf.evalFactorZolo(200);
  zrf.evalZolo(200,true);
  }

  bool HT=true; 
  if (HT) {
    ht h;
    std::cout << "setHT" << std::endl;
    int Nbase=17;
    int scale=1;
    int N=Nbase/scale;
    ftype xlow=1e-3;
    ftype xhigh=20;
    h.setHT(N);
    h.writeFactors();
    std::cout << "writeHT" << std::endl;
    int Nplot=200;
    bool REFLECT=true;
    h.writeHT(Nplot,xlow,xhigh,scale,REFLECT);
  }

  return 0;
}
