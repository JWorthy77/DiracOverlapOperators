#include "ht.h"
#include <fstream>

ht::ht() {
}

ht::ht(int N) {
  setHT(N);
}

ht::ht(int N, ftype r1, ftype r2) {
  setHT(N);
  writeHT(200,r1,r2);
}

void ht::setHT(int N) {

  ftype one(1);
  ftype pi(M_PIl);

  std::cout << "ht::setHT: " << pi << std::endl;

  Ndw=N;
  if (N%2==0) { //  even
    std::cout << "N even" << std::endl;
    Nnum=N/2-1;
    Ndenom=N/2;
    std::cout << "resize" << std::endl;
    num.resize(Nnum);
    denom.resize(Ndenom);
    std::cout << "eval" << tan(ftype(0)) << std::endl;
    for(int j=1;j<=Nnum;j++){
      num[j-1]=-tan(j*pi/N)*tan(j*pi/N);
    }
    for(int j=1;j<=Ndenom;j++){
      denom[j-1]=-tan((j*one-one/2)*pi/N)*tan((j*one-one/2)*pi/N);
    }
    std::cout << "done" << std::endl;
    mult=one*N;
  }
  else { // odd
    std::cout << "N odd" << std::endl;
    Nnum=(N-1)/2;
    Ndenom=(N-1)/2;
    num.resize(Nnum);
    denom.resize(Ndenom);
    for(int j=1;j<=Nnum;j++){
      num[j-1]=-tan(j*pi/N)*tan(j*pi/N);
    }
    for(int j=1;j<=Ndenom;j++){
      denom[j-1]=-tan((j*one-one/2)*pi/N)*tan((j*one-one/2)*pi/N);
    }
    mult=one/ftype(N);
    mult=one*N;
  }
}

ftype ht::htsgn(ftype& x) {
      
  ftype val=mult*x;
  for(int j=0;j<Nnum;j++) {
    val=val*(1-x*x/num[j])/(1-x*x/denom[j]);
  }
  if (Ndenom>Nnum) {
    val=val/(1-x*x/denom[Nnum]);
  }
  return val;
}

void ht::evalHT(int Np,std::ofstream& ofile, ftype r1, ftype r2, int scale, bool REFLECT) {
  
  ftype dx=(log(r2)-log(r1))/(Np-1);
  ftype x=r1;
  ftype xdash=0;
  ftype s=ftype(scale);
  if (REFLECT) {
    for (int j=0;j<Np;j++){
      xdash=(Np-j-1)*dx;
      x=r1*exp(xdash);
      ftype sx=-s*x;
      ftype Sest=htsgn(sx);
      ofile << -x << " " << Sest << " " << Sest-1 << std::endl; 
//      xdash=(Np-j-1)*dx;
//      x=r1*exp(xdash);
  }

  }
  for (int j=0;j<Np;j++){
    ftype sx=s*x;
    ftype Sest=htsgn(sx);
    ofile << x << " " << Sest << " " << Sest-1 << std::endl; 
    x=x+dx;
    xdash=xdash+dx;
    x=r1*exp(xdash);
  }
}


void ht::writeHT(int Npts, ftype r1, ftype r2, int scale,bool REFLECT) {
  
  std::cout << std::setprecision(20);
  std::ofstream ofile;
  std::ostringstream ss1;
  std::ostringstream ss2;
  std::ostringstream ss3;
  ss1<<float(r1);
  ss2<<float(r2);
  ss3<<scale;
  ofile.open("ht"+std::to_string(Ndw)+"R"+ss1.str()+"_"+ss2.str()+"S"+ss3.str()+".dat"); 
  evalHT(Npts,ofile,r1,r2,scale,REFLECT); // no plot points
  ofile.close();
}

void ht::writeFactors() {

  std::ofstream ofile;

  ofile.open("htFactors"+std::to_string(Ndw)+".dat"); 
  ofile << "s(x) = d.x.sum_m (c(m)-x^2)/sum_n (g(n)-x^2)" << std::endl;
  ofile << "d: " << mult << std::endl;
  ofile << "c: " ; for (int j=0;j<Nnum;j++) ofile << num[j] << " " ; ofile <<std::endl;
  ofile << "d: " ; for (int j=0;j<Ndenom;j++) ofile << denom[j] << " " ; ofile <<std::endl;
  ofile.close();

  ofile.open("htDenominators"+std::to_string(Ndw)+".dat"); 
  for (int j=0;j<Ndenom;j++) ofile << denom[j] << std::endl;
  ofile.close();
}

