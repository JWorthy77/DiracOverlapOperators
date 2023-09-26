      module spectra
      use pacc
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine spectratask(fid,lambda,SRF)
!     calculate a bunch of spectra concurrently on different processors
      use gaugefield
      use rvmodule
      use ratfuncs
      use options
      use IOmodule
      use rvmod
      use rrparams
      use ritzmod
      use rrspectrum
#ifdef PARALLEL
      use basicparallelmod
#endif
      implicit none
      integer fid ! field id
      real*8 lambda(nev)
      type(sgnratfunc) SRF
      character(len=80) fname 
c      integer,parameter :: NHTterms=10
c      integer,parameter :: NZterms=24

#ifdef PARALLEL
      print *,"Sequential RRspectraRJN - rank:",rank,"file:",fid
#endif
      call startritz
c      call setHTcoeffs(NHTterms,SRF)
c      call setZoloCoeffs(NZterms,SRF,1d-4,20*one)
c      call setZoloCoeffs(NZterms,SRF,1d-3,20*one)
c      print *,"HT:",NHTterms

      call readThetaFileName(fid,fname)
#ifdef PARALLEL
      print *,rank,fname
#endif
      call readMPIConFile(fname,theta)
      call coef(u,theta)
      call spectrum(-MDW,SRF,lambda)
#ifdef PARALLEL
      print *,"Sequential RRspectraRJN - rank:",rank,"file:",fid,
     &                                                    "COMPLETED"
#endif
      return
      end subroutine spectratask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine RRspectraRJN()
      use gaugefield
      use rvmodule
      use ratfuncs
      use options
      use IOmodule
      use rvmod
      use rrparams
      use ritzmod
      use rrspectrum
      implicit none
      integer lseed,icf
      integer,parameter :: Ncf=100
      real*8 lambda(nev)
      type(sgnratfunc) SRF
      character(len=80) fname 
      integer,parameter :: NHTterms=40
      integer,parameter :: NZterms=20
      print *,"RRspectraRJN"

      lseed=137
      call setrn(lseed)
      call startritz
      call setHTcoeffs(NHTterms,SRF)
c      call setZoloCoeffs(NZterms,SRF,1d-4,20*one)
c      call setZoloCoeffs(NZterms,SRF,1d-3,10*one)
      print *,"HT:",NHTterms


      do icf=1,Ncf
        open(unit=11,file='NKspectraHT.dat',access='append',
     &           status='unknown',form='formatted')
        open(unit=88,file='DetailsNKspectraHT.dat',access='append',
     &           status='unknown',form='formatted')
        call readThetaFileName(500+10*icf,fname)
        call readMPIConFile(fname,theta)
        call coef(u,theta)
  
        call spectrum(-MDW,SRF,lambda)
        write(11,*) 500+10*icf,lambda
        close(11)
        close(88)
      end do

      return
      end subroutine RRspectraRJN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine estimateMaxEigenvalue(Nmax,DCASE,mass,lmax)
      use options
      use rvmodule
      use IOmodule
      use gaugefield
      use statsmod
      use basicdiracopsmod
      implicit none
      integer Nmax,DCASE
      real(prc) mass,lmax
     
      complex(prc),dimension(Nv,Ndc) :: R,TMP,DR
      integer i

      call setRVs(Nv*4,R)
      call normalise(R)
      open(unit=81,file='DetailsMaxEig.dat',access='append',
     &           status='unknown',form='formatted')
      do i=1,Nmax
        if (DCASE.eq.1) then 
          call DShamir(R,TMP,u,.false.,-MDW,czero)
          call DShamir(TMP,DR,u,.true.,-MDW,czero)
        elseif (DCASE.eq.2) then
          call DWilson(R,TMP,u,.false.,-MDW,czero)
          call DWilson(TMP,DR,u,.true.,-MDW,czero)
        endif
        R=DR
        lmax=mag(R)
        R=R/lmax
        print *,i,sqrt(lmax)
        write(81,*) i,sqrt(lmax)
      end do
      close(81)
      lmax=sqrt(lmax)

      return
      end subroutine estimateMaxEigenvalue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getThirringSpectra()
      use numbers
      use options
      use IOmodule
      use gaugefield
      use statsmod
      use ratfuncs
      use basicdiracopsmod
      use kernelspectrarange
      implicit none
      integer,parameter :: Ngf=1
      real(prc) Mmin,Mmax,dM
      real(prc) lmin(Ngf),lmax(Ngf),lminav,lmaxav,lminsd,lmaxsd
      real(prc) klow,khigh
      integer i,idx
      integer DCASE
      logical P1D    
      type(sgnratfunc) SRF
      character*80 fname
      complex(prc),dimension(Nv,Ndc) :: R,TMP,DR

      DCASE=3 ! 1=Shamir kernel, 2=Wilson kernel
              ! 3=DOL(Shamir)
              ! 4=DOL(Wilson), 41=DOL(Wilson)-MSCG, 42=VOL(Wilson)-MSCG, 43=(2+V(W)+Vdag(W))/4
              ! 5=DDW(Shamir), 6=DDW(Wilson), 7=KDDW4(Wilson), 8=KDDW4(Shamir)
              ! 10=VOL(Shamir), 11=VOLpoly(Shamir)
      P1D=.false. ! plus 1D?
      call setHTcoeffs(20,SRF)
      MTYPE=3
      MDW=1d0

c      open(unit=11,file='SRCompactDOLW12x12Bp25Ls20.dat',
      open(unit=11,file='SR12x12DWilsonBp40Ls20.dat',
     &            status='unknown',form='formatted')
      idx=0
      do i=95+5,95+Ngf*5,5
c        call readThirringConFile(i,theta)
        call readThetaFileName(i,fname)
        call readMPIConFile(fname,theta)
        call coef(u,theta)
c        R=1
c        call DWilson(R,TMP,u,.false.,-MDW,czero)
c        call DWilson(TMP,DR,u,.true.,-MDW,czero)
c        print *,minval(abs(TMP)),maxval(abs(TMP))
c        print *,minval(abs(DR)),maxval(abs(DR))
c        stop
c        call calcRayleighRitzEigs()
c        stop
        idx=idx+1
c        call calcDEigs(3,klow,khigh,2,P1D)
c        klow=1q-1
c        khigh=10q0
c        call setZoloCoeffs(4,SRF,klow,khigh)
c        print *,DCASE
c        call calcDEigs(20,lmin(idx),lmax(idx),DCASE,P1D,SRF) 
       call estimateKernelExtrema(10,1,lmax(idx))
       call estimateKernelExtrema(10,-1,lmin(idx))
        write(11,*) lmin(idx),lmax(idx),lmax(idx)/lmin(idx)
        flush(11)
        print *, lmin(idx),lmax(idx),lmax(idx)/lmin(idx)
      end do
      call calcVarReal(Ngf,lmin,lminav,lminsd)
      call calcVarReal(Ngf,lmax,lmaxav,lmaxsd)
      write(11,*) lminav,lminsd,lmaxav,lmaxsd
      print *,"lmin: ",lminav,lminsd
      print *,"lmax: ",lmaxav,lmaxsd
      close(11)
      return
      end subroutine getThirringSpectra
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcRayleighRitzEigs()
!     calculate Wilson or Shamir min eigs
      use gaugefield
      use basicdiracopsmod
      use overlapmoduledev
c      use shamirmodule
      use domainwallmod
      use rvmodule
      use options
      use axbmodule1
      implicit none
      complex(prc),dimension(Nv,Ndc) :: R,Rp,Rm,P,Pp,gP,gPp,TMP,DR
      real(prc) :: moo,moop,dmoo,da,m,beta
      integer :: e,l

!     calc max eig
      call setRVs(Nv*4,R)
      call normalise(R)
      Rm=R
      call calcRRSearch(R,P)
      do e=1,1000
        call calcRRSearch(R,gP)
        da=1d-6
        call calcRitzFunctional(R,moo)
        Rp=R+da*P
        call calcRitzFunctional(Rp,moop)
        m=sign(1d0,moop-moo)
        da=1d-1
        do l=1,75
          call calcRitzFunctional(R,moo)
          Rp=R-m*da*P
          call calcRitzFunctional(Rp,moop)
          dmoo=moop-moo
c          print *,da,dmoo,moo
          if (dmoo > 0) then
            R=Rm
            da=da/2
          else
            Rm=R
            R=Rp
          end if
        end do
        call calcRRSearch(Rp,gPp)
        beta=sum(conjg(gPp)*gPp)/sum(conjg(gP)*gP)
        Pp=gPp+beta*(P-Rp*sum(conjg(Rp)*P)/sum(conjg(Rp)*Rp))

        R=Rp
        P=Pp
        print *,moo,moo**0.5
      end do


      call DWilson(R,TMP,u,.false.,-MDW,czero)
      call DWilson(TMP,DR,u,.true.,-MDW,czero)
      print *,moo,moo**0.5,maxval(abs(DR-moo*R))
      stop

      return
      end subroutine calcRayleighRitzEigs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcRitzFunctional(R,moo)
!     calculate Ritz functional
      use options 
      use gaugefield
      use basicdiracopsmod
      implicit none
      complex(prc),dimension(Nv,Ndc),intent(in) :: R
      real(prc),intent(out) :: moo
      complex(prc),dimension(Nv,Ndc) :: TMP,DR

      call DWilson(R,TMP,u,.false.,-MDW,czero)
      call DWilson(TMP,DR,u,.true.,-MDW,czero)
      moo=sum(conjg(R)*DR)/sum(conjg(R)*R)

      return
      end subroutine calcRitzFunctional
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcRRSearch(R,P)
!     calculate Rayleight-Ritz search vector
      use options 
      use gaugefield
      use basicdiracopsmod
      implicit none
      complex(prc),dimension(Nv,Ndc),intent(in) :: R
      complex(prc),dimension(Nv,Ndc),intent(out) :: P
      real(prc) :: moo
      complex(prc),dimension(Nv,Ndc) :: DR,TMP

      call calcRitzFunctional(R,moo)
      call DWilson(R,TMP,u,.false.,-MDW,czero)
      call DWilson(TMP,DR,u,.true.,-MDW,czero)
      DR=DR-moo*R
      P=DR/sum(conjg(R)*R)

      return
      end subroutine calcRRSearch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine calcDEigs(Nmax,lmin,lmax,DCASE,P1D,SRF)
c!     calculate Wilson or Shamir min/max eigs
c      use gaugefield
c      use basicdiracopsmod
c      use overlapmoduledev
cc      use shamirmodule
c      use domainwallmod
c      use rvmodule
c      use options
c      use axbmodule1
c      implicit none
c      integer,intent(in) :: Nmax
c      real(prc),intent(out) :: lmax,lmin
c      integer,intent(in) :: DCASE
c      logical,intent(in) :: P1D
c      type(sgnratfunc),intent(in),optional :: SRF
c      complex(prc),dimension(Nv,Ndc) :: R,DR,TMP
c      complex(prc),dimension(Nv,Ndc,Ls) :: RL,DRL,TMPL
c      real ev(1)
c      integer i
c      procedure(),pointer :: Vptr => NULL()
c      procedure(),pointer :: Mptr => NULL()
c      procedure(),pointer :: Dptr => NULL()
c
cc      goto 1023
c
c!     calc max eig
c      call setRVs(Nv*4,R)
c      call setRVs(Nv*4*Ls,RL)
c      call normalise(R)
c      call normalisev(Nv*Ndc*Ls,RL)
c      do i=1,Nmax
c        if(.not.P1D)then
c          if (DCASE.eq.1) then 
c            call DShamir(R,TMP,u,.false.,-MDW,czero)
c            call DShamir(TMP,DR,u,.true.,-MDW,czero)
c          elseif (DCASE.eq.2) then
c            call DWilson(R,TMP,u,.false.,-MDW,czero)
c            call DWilson(TMP,DR,u,.true.,-MDW,czero)
c          elseif (DCASE.eq.0) then
c            call DHermWilson(R,TMP,u,.false.,-MDW,czero)
c            call DHermWilson(TMP,DR,u,.true.,-MDW,czero)
c          elseif (DCASE.eq.3) then 
c            call DOLS(R,TMP,u,.false.,baremass,SRF)
c            call DOLS(TMP,DR,u,.true.,baremass,SRF)
c          elseif (DCASE.eq.4) then 
c            call DOverlap(R,TMP,u,.false.,baremass,SRF)
c            call DOverlap(TMP,DR,u,.true.,baremass,SRF)
c          elseif (DCASE.eq.41) then
c            Vptr => VOLMpf
c            Mptr => DdagDpC
c            Dptr => DWilson 
cc            print *,'DOL 1'
c            call DOL(R,TMP,u,.false.,baremass,SRF,Vptr,Mptr,Dptr)
cc            print *,'DOL 2'
c            call DOL(TMP,DR,u,.true.,baremass,SRF,Vptr,Mptr,Dptr)
cc            print *,'DOL 3'
c          elseif (DCASE.eq.42) then
cc            Vptr => VOLMpf
c            Mptr => DdagDpC
c            Dptr => DWilson 
cc            print *,"hello"
c            call VOLMpf(R,TMP,u,.false.,-MDW,SRF,Mptr,Dptr)
cc            print *,"goodbye"
c            call VOLMpf(TMP,DR,u,.true.,-MDW,SRF,Mptr,Dptr)
c          elseif (DCASE.eq.43) then
c            Vptr => VOLMpf
c            Mptr => DdagDpC
c            Dptr => DWilson 
c            call VOLNK(R,DR,u,.false.,-MDW,SRF,Vptr,Mptr,Dptr)
c          elseif (DCASE.eq.10) then 
c            call VOLSpf(R,TMP,u,.false.,baremass,SRF)
c            call VOLSpf(TMP,DR,u,.true.,baremass,SRF)
c          elseif (DCASE.eq.11) then 
c           call VOLSpoly(R,TMP,u,.false.,baremass,SRF)
c            call VOLSpoly(TMP,DR,u,.true.,baremass,SRF)
c          elseif(DCASE.eq.7)then
c            DWkernel=2 ! Wilson
c            call KDDW4(R,TMP,u,.false.,baremass)
c            call KDDW4(TMP,DR,u,.true.,baremass)
c          elseif(DCASE.eq.8)then
c            DWkernel=1 ! Shamir
c            call KDDW4(R,TMP,u,.false.,baremass)
c            call KDDW4(TMP,DR,u,.true.,baremass)
c          endif
c          R=DR
c          lmax=mag(R)
c          R=R/lmax
c        elseif(P1D)then
c          if (DCASE.eq.5)then 
c            call DDW_Shamir(RL,TMPL,u,.false.,baremass)
c            call DDW_Shamir(TMPL,DRL,u,.true.,baremass)
c          elseif(DCASE.eq.6)then
c            call DDW_Wilson(RL,TMPL,u,.false.,baremass)
c            call DDW_Wilson(TMPL,DRL,u,.true.,baremass)
cc          elseif(DCASE.eq.7)then
cc            DWkernel=2 ! Wilson
cc            call KDDW4(RL,TMPL,u,.false.,baremass)
cc            call KDDW4(TMPL,DRL,u,.true.,baremass)
cc          elseif(DCASE.eq.8)then
cc            DWkernel=1 ! Shamir
cc            call KDDW4(RL,TMPL,u,.false.,baremass)
cc            call KDDW4(TMPL,DRL,u,.true.,baremass)
c          endif
c          RL=DRL
c          lmax=magv(Nv*Ndc*Ls,RL)
c          RL=RL/lmax
c        endif
c        print *,i,sqrt(lmax)
c      end do
c      lmax=sqrt(lmax)
c
cc      stop
c
cc1023  continue
c!     calc min eig
c      call setRVs(Nv*4,R)
c      call setRVs(Nv*4*Ls,RL)
c      call normalise(R)
c      call normalisev(Nv*Ndc*Ls,RL)
c      Mptr => DdagD
c      do i=1,Nmax
c        if(.not.P1D)then
c          if (DCASE.eq.1) then 
c            call IDW(R,TMP,u,.true.,-MDW,czero)
c            call DWilson(TMP,DR,u,.false.,-MDW,2*cone)
c            call IDW(DR,TMP,u,.false.,-MDW,czero)
c            call DWilson(TMP,DR,u,.false.,-MDW,2*cone)
c          elseif (DCASE.eq.2) then
c            Mptr => DdagD
c            call IM(R,DR,u,.false.,-MDW,czero,Mptr)
c          elseif (DCASE.eq.0) then
c            Mptr => H2pC
c            call IM(R,DR,u,.false.,-MDW,czero,Mptr)
c          elseif (DCASE.eq.3) then 
c            call IDOLS(R,TMP,u,.false.,baremass,SRF)
c            call IDOLS(TMP,DR,u,.true.,baremass,SRF)
c          elseif (DCASE.eq.4) then 
c            call IDOverlap(R,TMP,u,.false.,baremass,SRF)
c            call IDOverlap(TMP,DR,u,.true.,baremass,SRF)
c          elseif (DCASE.eq.41) then 
c            print *,"41"
c            MTYPE=3
c            DWkernel=2
c            Vptr => VOLMpf
c            Mptr => DdagDpC
c            Dptr => DWilson 
c            call IDOL(R,TMP,u,.false.,baremass,SRF,Vptr,Mptr,Dptr)
c            call IDOL(TMP,DR,u,.true.,baremass,SRF,Vptr,Mptr,Dptr)
c          elseif (DCASE.eq.43) then
c            Vptr => VOLMpf
c            Mptr => DdagDpC
c            Dptr => DWilson 
c            call IVNK(R,DR,u,.false.,-MDW,SRF,Vptr,Mptr,Dptr)
c          elseif(DCASE.eq.7)then
c            DWkernel=2 ! Wilson
c            call IKDDW4(R,TMP,u,.false.,baremass)
c            call IKDDW4(TMP,DR,u,.true.,baremass)
c          elseif(DCASE.eq.8)then
c           DWkernel=1 ! Shamir
c            call IKDDW4(R,TMP,u,.false.,baremass)
c            call IKDDW4(TMP,DR,u,.true.,baremass)
c          endif
c          R=DR
c          lmin=mag(R)
c          R=R/lmin
c        else
c          if (DCASE.eq.5) then 
c            DWkernel=1 ! Shamir
c            call IDDW(RL,TMPL,u,.false.,baremass)
c            call IDDW(TMPL,DRL,u,.true.,baremass)
c          elseif(DCASE.eq.6)then
c            DWkernel=2 ! Wilson
c            call IDDW(RL,TMPL,u,.false.,baremass)
c            call IDDW(TMPL,DRL,u,.true.,baremass)
cc          elseif(DCASE.eq.7)then
cc            DWkernel=2 ! Wilson
cc            call IKDDW4(RL,TMPL,u,.false.,baremass)
cc            call IKDDW4(TMPL,DRL,u,.true.,baremass)
cc          elseif(DCASE.eq.8)then
cc            DWkernel=1 ! Shamir
cc            call IKDDW4(RL,TMPL,u,.false.,baremass)
cc            call IKDDW4(TMPL,DRL,u,.true.,baremass)
c          endif
c          RL=DRL
c          lmin=magv(Nv*Ndc*Ls,RL)
c          RL=RL/lmin
c        endif
c        print *,i,1/sqrt(lmin)
c      end do
c     lmin=one/sqrt(lmin)
c      call calceigs('SM',1,ev,1,Nv)
c      print *,"ev:",ev
c      print *,"Spec Rad:",lmax/lmin
c
c      call DWilson(R,TMP,u,.false.,-MDW,czero)
c      call DWilson(TMP,DR,u,.true.,-MDW,czero)
c      print *,"eigerr:",maxval(abs(DR-lmin*lmin*R))
c
c      return
c      end subroutine calcDEigs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function mag(V)
      implicit none
      complex(prc) :: V(Nv*Ndc)
      real(prc) :: rV(Nv*Ndc),iV(Nv*Ndc)

      rV=real(V,prc)
      iV=dimag(V)
      mag=dot_product(rV,rV)+dot_product(iV,iV)
      mag=sqrt(mag)
      return
      end function mag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine normalise(V)
      implicit none
      complex(prc) :: V(Nv*Ndc)
      real(prc) :: nrm

      nrm=mag(V)
      V=V/nrm
      return
      end subroutine normalise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function magv(N,V)
      implicit none
      integer N
      complex(prc) :: V(N)
      real(prc) :: rV(N),iV(N)

      rV=real(V,prc)
      iV=dimag(V)
      magv=dot_product(rV,rV)+dot_product(iV,iV)
      magv=sqrt(magv)
      return
      end function magv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine normalisev(N,V)
      implicit none
      integer N
      complex(prc) :: V(N)
      real(prc) :: nrm

      nrm=magv(N,V)
      V=V/nrm
      return
      end subroutine normalisev
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef USEARPACK
      subroutine evalSpectra()
      use IOmodule
      use gaugefield
      use arpackmodule
      use arpackmodule2
      implicit none
      integer,parameter :: Neigs=5
      integer ig,Ng,kernel,Ntheta
      real evs(Neigs),evl(Neigs)

      Ng=1
      print *,'Calculate Spectra'
      call setHTcoeffs(2,ASRF)
c      open(unit=11,file='DWilsonEigs.dat',
c     &            status='unknown',form='formatted')
      do ig=1,Ng
        Ntheta=95+ig*5
        print *,'ig:',ig
c        call readThirringGaugeField(ig,theta)
        call readThirringConFile(Ntheta,theta)
        call coef(u,theta)
        kernel=3
c        call calceigs('SM',Neigs,evs,kernel,Nv)
c        print *,'small',evs
c        call calceigs('SI',Neigs,ev,kernel,Nv)
c        print *,'small imaginary part',ev
c        call calceigs('SR',Neigs,ev,kernel,Nv)
c        print *,'small real part',ev
c        call calceigs('LM',Neigs,evl,kernel,Nv)
        call calceigsOL('LM',Neigs,evl,kernel,Nv)
        print *,'large',evl
c        write(11,*) Ntheta,evs,evl
c        kernel=2
c        call calceigs('SM',Neigs,ev,kernel,Nv)
c        print *,'small',ev
c        call calceigs('LM',Neigs,ev,kernel,Nv)
c        print *,'large',ev
      enddo
c      close(11)

      return
      end subroutine evalSpectra
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module spectra
