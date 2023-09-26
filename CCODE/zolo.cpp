#include "zolo.h"
#include "elliptic.h"
#include <fstream>

zolo::zolo() {
}


zolo::zolo(ftype lmin, ftype lmax, int N) {

  setZolo(lmin, lmax, N);

}

void zolo::setZolo(ftype lmin, ftype lmax, int N) {

  std::cout << std::setprecision(20);
  beta = lmax/lmin; if (VBS) std::cout << "beta: " << beta << std::endl;
  this->lmin = lmin; if (VBS) std::cout << "lmin: " << lmin << std::endl;
  this->lmax = lmax; if (VBS) std::cout << "lmax: " << lmax << std::endl;
  this->Ndw = N; if (VBS) std::cout << "N: " << N << std::endl;
  if (N%2==0) {
    this->Nnum = N/2-1;
    this->Ndenom = N/2;
  }
  else {
    this->Nnum = (N-1)/2;
    this->Ndenom = (N-1)/2;
  }
  k = ftype(1)/beta; if (VBS) std::cout << "k: " << k << std::endl;
  kdash = sqrt(ftype(1)-k)*sqrt(ftype(1)+k); if (VBS) std::cout << "kdash: " << kdash << std::endl;
  elliptic e;
  K = e.calcK(k,err); if (VBS) std::cout << "K: " << K << std::endl;
  Kdash = e.calcK(kdash,err); if (VBS) std::cout << "Kdash: " << Kdash << std::endl;
  cnum.resize(Nnum); // N/2
  cdenom.resize(Ndenom); // N/2
  setInnerCoeffs();
  M = calcM(); if (VBS) std::cout << "M: " << M << std::endl;
  C = calcC(); if (VBS) std::cout << "C: " << C << std::endl;
  lambda = calcLambda(); if (VBS) std::cout << "lambda: " << lambda << std::endl;
  mult=2*lambda/(1+lambda)/M; 
  getRoots();
  getExtrema();
}

void zolo::setInnerCoeffs() {

  elliptic e;
  for (int l=0;l<Nnum;l++) {
    ftype argC = ftype(2*l+2)*Kdash/ftype(Ndw); if (VBS) std::cout << l << " argC: " << argC << std::endl;
    ftype JargC = e.sn(argC,kdash,err); if (VBS) std::cout << l << " JargC: " << JargC << std::endl;
    cnum[l]=-JargC/(1-JargC)*JargC/(1+JargC);
    if(VBS) std::cout << "cnum(l): " << cnum[l] << std::endl; 
//    exit(0);
  }
  for (int l=0;l<Ndenom;l++) {
    ftype argCdash = (2*l+1)*Kdash/Ndw;
    ftype JargCdash = e.sn(argCdash,kdash,err); 
    cdenom[l] = -JargCdash/(1-JargCdash)*JargCdash/(1+JargCdash);
    if(VBS) std::cout << "cdenom(l): " << cdenom[l] << std::endl;
  }
}

ftype zolo::Rhalf(ftype& x) { 
  ftype R=1;
  for(int l=0;l<Nnum;l++) {
    R=R*(1-x/cnum[l])/(1-x/cdenom[l]);
  }
  if (Ndenom>Nnum) {
    R=R/(1-x/cdenom[Nnum]);
  }
  return R;
}

ftype zolo::calcM() {

  ftype one=1;
  return Rhalf(one);
}

ftype zolo::calcC() { 

  elliptic e;
  ftype bKD = Kdash/Ndw;
  ftype esn = e.sn(K,k,err);
  ftype edn = e.dn(K,k,err);
  ftype dsn = e.sn(bKD,kdash,err); 
  ftype ddn = e.dn(bKD,kdash,err); 
  ftype num = esn*ddn;
  ftype denom = edn*edn*dsn*dsn;
  return num/(1-denom);
}

ftype zolo::calcLambda() { 

  ftype C2 = C*C;
  return M/(C*Rhalf(C2));
}

ftype zolo::zsgn(ftype& x) {
      
  ftype arg=x/lmin*x/lmin;
  return x/lmin*mult*Rhalf(arg);
}

void zolo::getRoots() { 
  // this uses Chiu's formula for the roots, but uses the lambda calculated
  // with Kennedys formulation (ie no theta function)
  std::cout << "getRoots" << std::endl;
  roots.resize(Ndw);
  elliptic e;
  std::vector<ftype> vs(Ndw);
  ftype arg = sqrt( ftype(1.0+3.0*lambda)/pow(ftype(1+lambda),3) );
  ftype lkappa=sqrt(ftype(1+lambda)*ftype(1-lambda));
  ftype asn = e.arcsn(arg,lkappa,ftype(1e-20));
  ftype c1 = M*asn;
  ftype c2 = 2*Kdash/Ndw;  
  for(int s=0;s<Ndw;s++) {
    int so2 = (s+1)/2; // 1/2=0, 3/2=1, ...
    vs[s] = pow(ftype(-1),s)*c1 + so2*c2;
  }

  for(int s=0;s<Ndw;s++) {
    ftype snda = e.sn(vs[s],kdash,err); 
    ftype sn2=snda*snda;
    ftype arg = 1 - kdash*kdash*sn2;
    ftype omega = sqrt(arg)/lmin;
    roots[s] = 1/omega;
    if (VBS) std::cout << "roots["<< s << "]: " << roots[s] << std::endl;
  }
}

void zolo::getExtrema() { 
  std::cout << "getExtrema" << std::endl;
  extrema.resize(Ndw+1);
  elliptic e;
  ftype del = Kdash/Ndw;  
  for(int s=0;s<=Ndw;s++) {
    ftype sn=e.sn(s*del,kdash,err);
    ftype arg=1-kdash*kdash*sn*sn;
    extrema[s] = lmin/sqrt(arg);
  }
}

void zolo::writeZolo(int Npts) {
  
  std::ofstream ofile;
  std::ostringstream ss1;
  std::ostringstream ss2;
  ss1<<float(lmin);
  ss2<<float(lmax);
  ofile.open("zolo"+std::to_string(Ndw)+"R"+ss1.str()+"_"+ss2.str()+".dat"); 
  evalZolo(Npts,ofile); // no plot points
  ofile.close();
  ofile.open("zroots"+std::to_string(Ndw)+"R"+ss1.str()+"_"+ss2.str()+".dat"); 
  for(int j=0;j<Ndw;j++) ofile << roots[j] << " " << 1 << " " << 0 << std::endl;
  ofile.close();
  ofile.open("zextrema"+std::to_string(Ndw)+"R"+ss1.str()+"_"+ss2.str()+".dat"); 
  for(int j=0;j<=Ndw;j++) ofile << extrema[j] << " " << 1 << " " << 0 << std::endl;
  ofile.close();
}

void zolo::testZolo(int Nz) {
  
  ftype lmin=ftype(1)/1000 ; ftype lmax=10;
  setZolo(lmin,lmax,Nz);
  getRoots();
  std::ofstream ofile;
  ofile.open("zolo"+std::to_string(Nz)+".dat"); 
  evalZolo(200,ofile); // no plot points
  ofile.close();
  ofile.open("zroots"+std::to_string(Nz)+".dat");
  for(int j=0;j<Ndw;j++) ofile << roots[j] << " " << 1 << " " << 0 << std::endl;
  ofile.close();
}

void zolo::evalZolo(int Np,std::ofstream& ofile) {
  
  ftype xmin=lmin;
  ftype xmax=lmax;

  ftype dx=(log(xmax)-log(xmin))/(Np-1);
  ftype x=xmin;
  ftype xdash=0;
  for (int j=0;j<Np;j++){
    ftype Sest=zsgn(x);
    ofile << x << " " << Sest << " " << Sest-1 << std::endl; 
    x=x+dx;
    xdash=xdash+dx;
    x=xmin*exp(xdash);
  }
}

