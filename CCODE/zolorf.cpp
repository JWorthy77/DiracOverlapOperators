#include "zolorf.h"
#include "zolo.h"
#include <fstream>

zolorf::zolorf() {
}


zolorf::zolorf(zolo& z) {

  setFactorZolo(z); // makes av,dv
//  setPartialFractionZolo(av,dv,front,pfcoeffs); // makes front,pfcoeffs
  setPartialFractionZolo(); // makes front,pfcoeffs

}

void zolorf::setFactorZolo(zolo& z) {

  Nnum=z.Nnum;
  Ndenom=z.Ndenom;
  rmin=z.lmin;
  rmax=z.lmax;
  d=z.mult/z.lmin;
  av.resize(Nnum);
  dv.resize(Ndenom);
  for(int m=0;m<Ndenom;m++) d=-d*rmin*rmin*z.cdenom[m];
  for(int m=0;m<Nnum;m++) d=-d/(rmin*rmin*z.cnum[m]);
  for(int m=0;m<Nnum;m++) av[m]=z.lmin*z.lmin*z.cnum[m];
  for(int m=0;m<Ndenom;m++) dv[m]=z.lmin*z.lmin*z.cdenom[m];
  writeFactors();
}

void zolorf::writeFactors() {

  std::ofstream ofile;
  ofile.open("zoloFactorCoeffs.dat"); 
//  ofile << "s(x) = d.x.sum_m (c[n]-x^2)/sum_n (g[n]-x^2)" << std::endl;
  ofile << "s(x) = d.x.sum_m (x^2-c[n])/sum_n (x^2-g[n])" << std::endl;
  ofile << "d: " << d << std::endl;
  ofile << "c: "; for (int j=0;j<Nnum;j++) ofile << av[j] << " "; ofile << std::endl;
  ofile << "g: "; for (int j=0;j<Ndenom;j++) ofile << dv[j] << " "; ofile << std::endl;
  ofile.close();
}

ftype zolorf::evalFactorZolo(ftype& r) {

  ftype v=r*d;
//  for(int m=0;m<Nnum;m++) v=v*(av[m]-r*r);
//  for(int m=0;m<Ndenom;m++) v=v/(dv[m]-r*r);
  for(int m=0;m<Nnum;m++) v=v*(r*r-av[m]);
  for(int m=0;m<Ndenom;m++) v=v/(r*r-dv[m]);
  return v;
}

void zolorf::evalFactorZolo(int Np) {

  std::ofstream ofile;
  ofile.open("zoloFactor.dat"); 
  ftype dx=(log(rmax)-log(rmin))/(Np-1);
  ftype x=rmin;
  ftype xdash=0;
  for (int j=0;j<Np;j++){
    ftype Sest=evalFactorZolo(x);
    ofile << x << " " << Sest << " " << Sest-1 << std::endl; 
    x=x+dx;
    xdash=xdash+dx;
    x=rmin*exp(xdash);
  }
  ofile.close();
}

void zolorf::setPartialFractionZolo() {

  std::vector<ftype> num(av.size());
  for(int n=0;n<num.size();n++) num[n]=-av[n];
  std::vector<ftype> denom(dv.size());
  for(int n=0;n<denom.size();n++) denom[n]=-dv[n];
  setPartialFractionZolo(av,dv,front,pfcoeffs);
  writePartFracCoeffs();
}

void zolorf::setPartialFractionZolo(std::vector<ftype>& num, std::vector<ftype>& denom,
                                    std::vector<ftype>& front, std::vector<ftype>& pfcoeffs) {
// assumes denom has distinct real coefficients
// input is factored form (coefficients are zeros)
//  std::vector<ftype>& num=av;
//  std::vector<ftype>& denom=dv;

  int nn=num.size();
  int nd=denom.size();
  pfcoeffs.resize(nd);
  std::vector<ftype> enmr(nn+1);
  expandPoly(num,enmr);
  std::vector<ftype> edenom(nd+1);
  expandPoly(denom,edenom);
  int nf=nn-nd+1;
  if (nf>=1) {
    std::vector<ftype> pnum(nd);
    dividePoly(enmr,edenom,front,pnum);
    for(int j=0;j<nd;j++) {
      ftype z=denom[j];
      pfcoeffs[j]=evalpoly(z,pnum)/polyDeriv(z,denom);
    }
  }
  else {
    for (int j=0;j<nd;j++) {
      ftype z=denom[j];
      pfcoeffs[j]=evalfpoly(z,num)/polyDeriv(z,denom);
    }
  }
  if (VERBOSE) {
    for(int k=0;k<pfcoeffs.size();k++) std::cout << pfcoeffs[k] << std::endl;
    for(int k=0;k<pfcoeffs.size();k++) std::cout <<denom[k] << std::endl;
    for(int k=0;k<front.size();k++) std::cout << front[k] << std::endl;
  }
}

void zolorf::writePartFracCoeffs() {

  std::ofstream ofile;
  ofile.open("zoloPartFracCoeffs.dat"); 
  ofile << "s(x) = d.x.(sum_m f[m] x^m + sum_n c[n]/(g[n]-x^2)" << std::endl;
  ofile << "d: " << d << std::endl;
  ofile << "f: "; for (int j=0;j<front.size();j++) ofile << front[j] << " "; ofile << std::endl;
  ofile << "c: "; for (int j=0;j<pfcoeffs.size();j++) ofile << pfcoeffs[j] << " "; ofile << std::endl;
  ofile << "g: "; for (int j=0;j<dv.size();j++) ofile << dv[j] << " "; ofile << std::endl;
  ofile.close();

  ofile.open("zoloDenoms"+std::to_string(dv.size())+".dat"); 
  for (int j=0;j<dv.size();j++) ofile << dv[j] << std::endl;
  ofile.close();
}

void zolorf::expandPoly(std::vector<ftype>& zeros,std::vector<ftype>& epoly) {
//  expand factored poly
  int nn=zeros.size();
  epoly.resize(nn+1);
  epoly[0]=-zeros[0];
  epoly[1]=1;
  for(int j=2;j<=nn;j++) {
    epoly[j]=epoly[j-1]; // = 1
    for(int i=j-1;i>0;i--) {
      epoly[i]=epoly[i-1]-epoly[i]*zeros[j-1];
    }
    epoly[0]=-epoly[0]*zeros[j-1];
  }
}

void zolorf::dividePoly(std::vector<ftype>& enmr, std::vector<ftype>& edenom, 
                        std::vector<ftype>& front, std::vector<ftype>& dnum) {
// if the top poly has an equal or higher order than the bottom then want
// to split into a polynomial and a lower order rational function
  int nn = enmr.size();
  int nd = edenom.size();
  std::cout << "nn: " << nn << std::endl;
  std::cout << "nd: " << nd << std::endl;
  if (nn < nd) { // do nothing, just copy to output
    dnum.resize(nn);
    for(int j=0;j<nn;j++) dnum[j]=enmr[j];
//    front=0; // don't assign
    return;
  }
  else {
    if(VERBOSE) std::cout << "divide poly" << std::endl;
    std::vector<ftype> tnum(nn);
    for(int j=0;j<nn;j++) tnum[j]=enmr[j];
    int nf=1+nn-nd;
    front.resize(nf);
    for (int j=0;j<nf;j++) {
      front[nf-j-1]=tnum[nn-j-1]/edenom[nd-1];
      for (int k=0;k<nd;k++) {
        tnum[nn-j-k-1]=tnum[nn-j-k-1]-front[nf-j-1]*edenom[nd-k-1];
      }
    }
    dnum.resize(nd);
    for (int j=0;j<nd;j++) dnum[j]=tnum[j];
  }
}

ftype zolorf::evalpoly(ftype& x,std::vector<ftype>& coeffs) {
  // poly = c0+c1.x+c2.x^2+...
  int nc=coeffs.size();
  if (nc==1) return coeffs[0];
  ftype val=x*coeffs[nc-1]+coeffs[nc-2];
  for(int j=nc-3;j>=0;j--) {
    val=x*val+coeffs[j];
  }
  return val;
}

ftype zolorf::evalfpoly(ftype& x,std::vector<ftype>& zeros) {
  // poly = (x-c0)(x-c1)(x-c2)...
  int nc=zeros.size();
  ftype fpoly=1;
  for(int j=0;j<nc;j++) {
    fpoly=fpoly*(x-zeros[j]);
  }
  return fpoly;
}

ftype zolorf::polyDeriv(ftype& x,std::vector<ftype>& zeros) {
// calculate derivative of poly dercibed by zeros 

  int nc=zeros.size();
  ftype polyderiv=0;
  for(int j=0;j<nc;j++) {
    ftype st=1;
    for(int i=0;i<nc;i++) {
      if (i!=j) st=st*(x-zeros[i]);
    }
    polyderiv=polyderiv+st;
  }
  return polyderiv;
}

ftype zolorf::evalPartialFractionZolo(ftype& r) {

  Nfront=front.size();
  ftype pf=ftype(0);
  if(Nfront==1) pf=front[0];
  if(Nfront>1) {std::cout << "Not implemented - Nfront: " << Nfront << std::endl; std::exit(0);}
//  for(int n=0;n<Ndenom;n++) pf=pf+pfcoeffs[n]/(dv[n]-r*r);
  for(int n=0;n<Ndenom;n++) pf=pf+pfcoeffs[n]/(r*r-dv[n]);
  pf=d*r*pf;
//  if (Nfront==1) pf=-pf;
  return pf;
}

void zolorf::evalZolo(int Np,bool PF) {

  std::ofstream ofile;
  if(PF){
    ofile.open("zoloPartFrac.dat"); 
  }
  else {
    ofile.open("zoloFactor.dat"); 
  }
  ftype dx=(log(rmax)-log(rmin))/(Np-1);
  ftype x=rmin;
  ftype xdash=0;
  ftype Sest;
  for (int j=0;j<Np;j++){
    if(PF){
      Sest=evalPartialFractionZolo(x);
    }
    else {
      Sest=evalFactorZolo(x);
    }
    ofile << x << " " << Sest << " " << Sest-1 << std::endl; 
    x=x+dx;
    xdash=xdash+dx;
    x=rmin*exp(xdash);
  }
  ofile.close();
}

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makePartFrac(num,denom,front,pfcoeffs)
!     assumes denom has distinct real coefficients
!     passes in factored form (coefficients are zeros)
      implicit none
      real(prc),dimension(:) :: num,denom
      real(prc),dimension(:) :: front,pfcoeffs
      real(prc),allocatable,dimension(:) :: enum,edenom,pnum
      real(prc) z
      integer j,nn,nd,nf

      nn=size(num)
      nd=size(denom)
      allocate(enum(nn+1))
      call expandPoly(num,enum)
      allocate(edenom(nd+1))
      call expandPoly(denom,edenom)
      nf=nn-nd+1
      if (nf.ge.1) then
        allocate(pnum(nd))
        call dividePoly(enum,edenom,front,pnum)
        do j=1,nd
          z=denom(j)
          pfcoeffs(j)=evalpoly(z,pnum)/polyDeriv(z,denom)
        end do
        deallocate(pnum)
      else
        do j=1,nd
          z=denom(j)
          pfcoeffs(j)=evalfpoly(z,num)/polyDeriv(z,denom)
        end do
      end if
      deallocate(enum)
      deallocate(edenom)
      
      return
      end subroutine makePartFrac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
void zolorf::testPartFrac() {

  std::cout << "testPartFrac" << std::endl;

  std::vector<ftype> num1(2);
  std::vector<ftype> num2(4);
  std::vector<ftype> denom(3);
  std::vector<ftype> pfront1;
  std::vector<ftype> pfront2(2);
  std::vector<ftype> pf;

  num1[0]=4; num1[1]=5;
  denom[0]=1; denom[1]=2; denom[2]=3;
  setPartialFractionZolo(num1,denom,pfront1,pf);

  std::cout << "pf: " ;
  for(int l=0;l<pf.size();l++) std::cout << pf[l] << " ";
  std::cout << std::endl;
  std::cout << "should be 6 -6 1" << std::endl;
  std::cout << "size of pfront1: " << pfront1.size() << std::endl;
  std::cout << "should be unassigned (size 0)" << std::endl;

  num2[0]=num1[0]; num2[1]=num1[1]; num2[2]=6; num2[3]=7;
  setPartialFractionZolo(num2,denom,pfront2,pf);

  std::cout << std::endl << "pf2: " ;
  for(int l=0;l<pf.size();l++) std::cout << pf[l] << " ";
  std::cout << std::endl;
  std::cout << "should be 180 -120 12" << std::endl;
  std::cout << std::endl << "pfront2: " ;
  for(int l=0;l<pfront2.size();l++) std::cout << pfront2[l] << " ";
  std::cout << std::endl;
  std::cout << "should be -16 1" << std::endl;

  std::vector<ftype> num3(2); num3[0]=2; num3[1]=4;
  std::vector<ftype> denom3(3); denom3[0]=1; denom3[1]=3; denom3[2]=5;
  std::vector<ftype> pfront3;
  std::vector<ftype> pf3;
  setPartialFractionZolo(num3,denom3,pfront3,pf3);
  std::cout << std::endl << "pf3: " ;
  for(int l=0;l<pf3.size();l++) std::cout << pf3[l] << " ";
  std::cout << std::endl;
  std::cout << "should be 3/8 1/4 3/8" << std::endl;
  std::cout << "size of pfront3: " << pfront1.size() << " should be unassigned (size 0)" << std::endl;

}

void zolorf::testDividePoly() {
   
  std::cout << "testDividePoly" << std::endl;
  std::vector<ftype> num(5);
  num[0]=0; num[1]=1 ; num[2]=1 ; num[3]=1 ; num[4]=1; // (/0d0,1d0,1d0,1d0,1d0/)
  std::vector<ftype> denom(3);
  denom[0]=1; denom[1]=0 ; denom[2]=1; // (/1d0,0d0,1d0/)
  std::vector<ftype> front;
  std::vector<ftype> dnum;

  std::cout << "num: " ;
  for(int l=0;l<5;l++) std::cout << num[l] << " ";
  std::cout << std::endl;
  std::cout << "denom: " ;
  for(int l=0;l<3;l++) std::cout << denom[l] << " ";
  std::cout << std::endl;
  dividePoly(num,denom,front,dnum);
  std::cout << "front: " ;
  for(int l=0;l<3;l++) std::cout << front[l] << " ";
  std::cout << std::endl;
  std::cout << "should be: 0 1 1" << std::endl;
  std::cout << "dnum: " ;
  for(int l=0;l<2;l++) std::cout << dnum[l] << " ";
  std::cout << std::endl;
  std::cout << "should be: 0 0" << std::endl;
}
