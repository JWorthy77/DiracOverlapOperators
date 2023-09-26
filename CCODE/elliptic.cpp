#include <iostream>
#include <iomanip>
#include <stdio.h>
//#include <gmp.h>

#include "elliptic.h"
#include <boost/math/special_functions/jacobi_elliptic.hpp>
//namespace mp = boost::multiprecision;

ftype elliptic::getpi(ftype err) {
// pi=(4.agm[1,1/2^1/2]^2)/(1-sum (2^j-1)cj^2))
// cj=cjm1^2/(4.aj)
  ftype ort2=1/sqrt(ftype(2));
  ftype c1=agm(1,ort2,err);
  ftype num=4*c1*c1;
  int n=0;
  int Nmax=10e3;
  ftype a=1 ; 
  ftype g=ort2;
  ftype sum=0;
  ftype ap1,gp1,cp1,c,pi,pip1,diff;
  pi=1;
  for(int j=1;j<Nmax;j++) {
    ap1=(a+g)/2;
    gp1=sqrt(a*g);
    cp1=(a-g)/2;
    a=ap1;
    g=gp1;   
    c=cp1;
    sum+=pow(2,j+1)*c*c;
    pip1=num/(1-sum);
    diff=abs(pip1-pi);
    pi=pip1;
    if (diff<err) break;
  }
  return pi;
}

ftype elliptic::calcK(ftype k, ftype err) {
// agm(x,y) = pi/4 . (x+y)/K[(x-y)/(x+y)]
// K(k) = pi / 2.agm[1,(1-k^2)^1/2]
  ftype pi=getpi(err);
  ftype fk=sqrt(1-k*k);
  ftype denom=2*agm(1,fk,err);
  return pi/denom;
}

ftype elliptic::agm(ftype x, ftype y, ftype err) {

  ftype a=(x+y)/2;
  ftype g=sqrt(x*y);
  ftype diff=a;
  ftype ap1;
  ftype gp1;
  int n=0;
  int Nmax=1000;
  while(diff>err && n<Nmax) {
    ap1=(a+g)/2;
    gp1=sqrt(a*g);
    diff=abs(ap1-a);
    a=ap1;
    g=gp1;   
    n++;
  }
  if (n>=Nmax) {
    std::cout << "agm didn't converge" << std::endl;
    exit(0);
  }
  std::cout << "agm converged after " << n << " steps: " << a << std::endl;
  return a;
}

ftype elliptic::calcam(ftype u, ftype m, ftype err) {

  ftype phiN;
  ftype phiNm1;
  int Nmax=30;
  int N;
  for(N=0;N<Nmax;N++) {

    ftype a=1;
    ftype b=sqrt(1.0-m);
    ftype c=sqrt(m);

    std::vector<ftype> av(N);
    std::vector<ftype> bv(N);
    std::vector<ftype> cv(N);
    ftype ap1;
    ftype bp1;
    ftype cp1;
    for(int i=0;i<N;i++) {
      ap1=(a+b)/2;
      bp1=sqrt(a*b);
      cp1=(a-b)/2;
      av[i]=ap1;
      bv[i]=bp1;
      cv[i]=cp1;
      a=ap1;
      b=bp1;
      c=cp1;
    }
    ftype phi=pow(ftype(2),N)*a*u;
//    std::cout << "phi: " << phi << std::endl;
    ftype phim1;
    ftype diff;
    for(int i=N-1;i>=0;i--) {
//      std::cout << cv[i]/av[i] << std::endl;
      phim1=(phi+asin(cv[i]/av[i]*sin(phi)))/2;
//      std::cout << "phim1: " << phim1 << std::endl;
      phi=phim1;
    }
    phiNm1=phiN;
    phiN=phi;
    diff = abs(phiN-phiNm1);
    if (N>0) std::cout << N << " am err: " << diff << " tol: " << err << std::endl;
    if (N>0) if (diff < err) break;
  }
  if (N==Nmax) {
    std::cout << "am did not converge" << std::endl;
    exit(0);
  }
  return phiN;
}

ftype elliptic::sn(ftype u, ftype m, ftype err) {

  ftype am = calcam(u,m,err);
  ftype sb = boost::math::jacobi_sn(m,u);
  ftype s = sin(am);
//  std::cout << "sb: " << sb << std::endl;
//  std::cout << "s: " << s << std::endl;
  return sb;
}

ftype elliptic::dn(ftype u, ftype m, ftype err) {

  ftype s = sn(u,m,err);
  ftype d = sqrt(ftype(1)-m*s*s);
  return d;
}

ftype elliptic::asnfn(ftype x,ftype k) {

  return pow( (1+x)*(1-x) , ftype(-0.5)) * pow( (1+k*x)*(1-k*x), ftype(-0.5));
}

ftype elliptic::arcsnN(ftype xmax, ftype k,int N) {
  // int_0^xmax 1/sqrt[ (1-t^2)(1-k^2t^2) ]
  ftype dx = xmax/N;
  ftype x = dx/2;
  ftype arcsn = 0;
  ftype sqrt3o5 = sqrt(ftype(3)/5)*dx/2;
  ftype FoN = ftype(5)/9;
  ftype EoN = ftype(8)/9;
  ftype half = ftype(1)/2;
  ftype rt6o5 = sqrt(ftype(6)/5);
  ftype s1 = sqrt(ftype(3)/7-2*rt6o5/7);
  ftype s2 = sqrt(ftype(3)/7+2*rt6o5/7);
  ftype w1=(18+sqrt(ftype(30)))/36;
  ftype w2=(18-sqrt(ftype(30)))/36;
  for(int i=1;i<=N;i++) {
//    ftype xm1 = x-sqrt3o5;
//    ftype xp1 = x+sqrt3o5;
//    ftype xm1 = x-half*dx;
//    ftype xp1 = x+half*dx;
    ftype xm1 = x-s1*dx/2;
    ftype xp1 = x+s1*dx/2;
    ftype xm2 = x-s2*dx/2;
    ftype xp2 = x+s2*dx/2;


    ftype vm2 = asnfn(xm2,k);
    ftype vm1 = asnfn(xm1,k);
    ftype vi = asnfn(x,k);
    ftype vp1 = asnfn(xp1,k);
    ftype vp2 = asnfn(xp2,k);
    arcsn = arcsn+(w2*vm2+w1*vm1+w1*vp1+w2*vp2)*dx/2;
//    arcsn = arcsn+(FoN*vm1+EoN*vi+FoN*vp1)*dx/2;
//    arcsn = arcsn+(half*vm1+half*vp1)*dx;
    x=x+dx;
  }
  return arcsn;
}

ftype elliptic::arcsn(ftype xmax, ftype k, ftype err) {

  int N=16;
  int Nmax=100000;
  ftype diff=1;
  ftype as=arcsnN(xmax,k,N);
  while(diff>err && N<Nmax) {
    N=2*N;
    ftype asp1=arcsnN(xmax,k,N);
    diff=abs(asp1-as);
    as=asp1;
  }
  std::cout << "arcsn N: " << N << " err: " << diff << std::endl;
  return as;
}



