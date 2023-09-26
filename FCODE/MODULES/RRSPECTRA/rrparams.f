      module rrparams
      implicit none

      integer neps,np,nev,nks
      integer SPECTRA_KERNEL
      logical KERNEL
      real*8 w1,w2,cv_fac,eps_err
c      parameter(neps=1000)
c      parameter(np=300,nev=3,nks=20)
      parameter(np=3000,nev=1,nks=20)
c      parameter(np=4000,nev=6,nks=30)
      parameter(w1=1.d-6,w2=1.d-5,cv_fac=0.1d0,eps_err=1.d-9)
c      parameter(w1=5.d-6,w2=5.d-5,cv_fac=0.1d0,eps_err=1.d-9)
c      parameter(w1=5.d-5,w2=5.d-4,cv_fac=0.1d0,eps_err=1.d-9)
c      parameter(w1=5.d-4,w2=5.d-3,cv_fac=0.1d0,eps_err=1.d-9)
c      parameter(SPECTRA_KERNEL=2) ! 1 is Wilson, 2 is Shamir
c      parameter(KERNEL=.true.) ! true if only calculate eigenvalues of kernel rather than overlap operator

      end module rrparams
