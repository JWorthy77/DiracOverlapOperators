      module gwmodule
      use arraysizes
      use numbers
      use options
      use ratfuncs
      use gaugefield
      use rvmodule
      use gammas
      use overlapmoduledev
      use domainwallmod
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine GWerrorKDW(R,err,SRF)
!     this is the Ls not -> infty error
      implicit none
      complex(prc),dimension(4*Nv),intent(in) :: R
      complex(prc) err
      type(sgnratfunc),intent(in) :: SRF
      complex(prc),dimension(4*Nv) :: DR,RR,DRL,DRR,DR3

      print *,'GWerror'
!     G5.D
      call KDDW4(R,DRL,u,.false.,baremass)
      call mGmu(DRL,5)
      print *,'G5.D'
      RR=R
      call mGmu(RR,5)
      call KDDW4(RR,DRR,u,.false.,baremass)
      print *,'+D.G5'
      call KDDW4(DRL,DR3,u,.false.,baremass)
      print *,'-D.G5.D'
      DR=DRL+DRR-2*DR3
      err=maxval(abs(DR))
      print *,"max err:",err
      return
      end subroutine GWerrorKDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine GWerrorDOL(R,err,SRF)
      implicit none
      complex(prc),dimension(4*Nv),intent(in) :: R
      complex(prc) err
      type(sgnratfunc),intent(in) :: SRF
      complex(prc),dimension(4*Nv) :: DR,RR,DRL,DRR,DR3

!     G3.D
      if (DWkernel.eq.1) then
        print *,'GWerror Shamir Overlap'
        call DOLS(R,DRL,u,.false.,baremass,SRF)
      elseif (DWkernel.eq.2) then
        print *,'GWerror Wilson Overlap'
        call Doverlap(R,DRL,u,.false.,baremass,SRF)
      endif
      call mGmu(DRL,4)
      print *,'G3.D'
      RR=R
      call mGmu(RR,4)
      if (DWkernel.eq.1) then
        call DOLS(RR,DRR,u,.false.,baremass,SRF)
      elseif (DWkernel.eq.2) then
        call Doverlap(RR,DRR,u,.false.,baremass,SRF)
      endif
      print *,'+D.G3'
      if (DWkernel.eq.1) then
        call DOLS(DRL,DR3,u,.false.,baremass,SRF)
      elseif (DWkernel.eq.2) then
        call DOverlap(DRL,DR3,u,.false.,baremass,SRF)
      endif
      print *,'-D.G3.D'
      DR=DRL+DRR-2*DR3
      err=maxval(abs(DR))
      print *,"max err:",err
      return
      end subroutine GWerrorDOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine GWcorrectionDOL(R,cor,SRF)
      implicit none
      complex(prc),dimension(4*Nv),intent(in) :: R
      complex(prc) cor
      type(sgnratfunc),intent(in) :: SRF
      complex(prc),dimension(4*Nv) :: DR,RR,DRL,DRR,DR3

!     G3.D
      if (DWkernel.eq.1) then
        print *,'GWerror Shamir Overlap'
        call DOLS(R,DRL,u,.false.,baremass,SRF)
      elseif (DWkernel.eq.2) then
        print *,'GWerror Wilson Overlap'
        call Doverlap(R,DRL,u,.false.,baremass,SRF)
      endif
      call mGmu(DRL,4)
      print *,'G3.D'
c      RR=R
c      call mGmu(RR,4)
c      if (DWkernel.eq.1) then
c        call DOLS(RR,DRR,u,.false.,baremass,SRF)
c      elseif (DWkernel.eq.2) then
c        call Doverlap(RR,DRR,u,.false.,baremass,SRF)
c      endif
c      print *,'+D.G3'
      if (DWkernel.eq.1) then
        call DOLS(DRL,DR3,u,.false.,baremass,SRF)
      elseif (DWkernel.eq.2) then
        call DOverlap(DRL,DR3,u,.false.,baremass,SRF)
      endif
      print *,'D.G3.D'
      DR=abs(2*DR3)
      cor=maxval(abs(DR))
      print *,"max cor:",cor
      return
      end subroutine GWcorrectionDOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine GamHermError(R,err,SRF)
c      use shamirmodule
!     this |G5.D.G5-Ddag|_inf
c      implicit none
c      complex(prc),dimension(4*Nv),intent(in) :: R
c      complex(prc) err
c      type(sgnratfunc),intent(in) :: SRF
c      complex(prc),dimension(4*Nv) :: TMP,DRL,DRR,DR
c
c      print *,'Gamma Hermiticity'
c!     G3.D.G3
c      TMP=R
c      call mGmu(TMP,4)
c      call Doverlap(TMP,DRL,u,.false.,baremass,SRF)
c!      call DOLS(TMP,DRL,u,.false.,baremass,SRF)
c      call mGmu(DRL,4)
c      print *,'G3.D.G3'
c      call Doverlap(R,DRR,u,.true.,baremass,SRF)
c!      call DOLS(R,DRR,u,.true.,baremass,SRF)
c      print *,'-Ddag'
c      DR=DRL+DRR
c      err=maxval(abs(DR))
c      print *,"max err:",err
c      return
c      end subroutine GamHermError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine testGWerror()
c      use arpackmodule
c      implicit none
c      logical ZOLO
c      complex(prc),dimension(4*Nv) :: R
c      complex(prc) errH,errZ
c      type(sgnratfunc) :: ZRF,HRF
c      real(prc) lmin,lmax
c      real ev(1)
c      integer Nz,Nh
c      
c      MTYPE=1
c      call setRVs(Nv*4,R)
c      call calceigs('SM',1,ev,2,Nv)
c      print *,'ev min:',ev
c      lmin=ev(1)
c      call calceigs('LM',1,ev,2,Nv)
c      print *,'ev max:',ev
c      lmax=ev(1)
c      do Nz=5,25,5
c        call setZoloCoeffs(Nz,ZRF,lmin,lmax)
c        Nh=2*Nz
c        call setHTcoeffs(Nh,HRF)
c        call GWerror(R,errZ,ZRF)
c        call GWerror(R,errH,HRF)
c        print *,Nz,errH,errZ
c        open(unit=11,file='GWerrWilson.dat',status='unknown',
c     &                              access='append',form='formatted')
c        write(11,*) Nz,errH,errZ
c        close(11)
c      end do
c
c      end subroutine testGWerror
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ThirringGWerror()
      use IOmodule
      use gaugefield
      use spectra
      use kernelspectrarange
      implicit none
      complex(prc),dimension(4*Nv) :: R
      type(sgnratfunc) :: SRF
      complex(prc) err
      integer idx,j
      real(prc) lmin,lmax
      integer,parameter :: Nerr=8
      integer n,NHT
      character(len=80) fname 
      
      call readThetaFileName(500,fname)
      call readMPIConFile(fname,theta)
      call coef(u,theta)
      call setRVs(Nv*4,R)
      NHT=100
      call setHTcoeffs(NHT,SRF)
      DWkernel=2 ! 1=Shamir,2=Wilson
      call GWerrorDOL(R,err,SRF) ! direct formulation only
      print *,NHT,err
      open(unit=11,file='play.dat',status='unknown',
     &                              access='append',form='formatted')
      write(11,*) NHT,err
      close(11)



      return
c      OLTYPE=2 ! 1 =direct,2=indirect
      DWkernel=2 ! 1=Shamir,2=Wilson
      MTYPE=1 ! 1=standard,3=alt3, 2=transpose alt3
c      ThirringFileDir=
c     &           '/home/jude/2019/Thirring/Sunbird/converted/b_0.30/'
      call readConvertedThetaFileName(1,fname)
      call reaaConverteddThirringGaugeField(fname,theta)
      call coef(u,theta)
      call setRVs(Nv*4,R)

c      do n=10,80,10
c        call setHTcoeffs(Ls,SRF)
       call estimateKernelExtrema(10,1,lmax)
       call estimateKernelExtrema(10,-1,lmin)
c        lmin=real(0.001,prc)
c        lmax=real(9,prc)
        call setZoloCoeffs(Ls,SRF,lmin,lmax)
!        call GWerror(R,err,SRF)
        call GWerrorDOL(R,err,SRF) ! direct formulation only
        print *,Ls,err
        open(unit=11,file='play.dat',status='unknown',
     &                              access='append',form='formatted')
        write(11,*) Ls,err
        close(11)
c      end do
      return
      end subroutine ThirringGWerror
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine testGamHerm()
c      use IOmodule
c      use gaugefield
c      implicit none
c      complex(prc),dimension(4*Nv) :: R
c      type(sgnratfunc) :: SRF
c      complex(prc) err
c      integer idx,j
c      integer,parameter :: Nerr=8
c      complex(prc) errvec(Nerr)
c      integer n
c      
c      OLTYPE=1 ! 1 =direct,2=indirect
c      DWkernel=2 ! 1=Shamir,2=Wilson
c      MTYPE=3 ! 1=standard,3=alt3, 2=transpose alt3
c      ThirringFileDir=
c     &           '/home/jude/2019/Thirring/Sunbird/converted/b_0.30/'
c      call readThirringGaugeField(1,theta)
c      call coef(u,theta)
c      call setRVs(Nv*4,R)
c
c      do n=10,50,10
c        call setHTcoeffs(n,SRF)
c        call GamHermError(R,err,SRF)
c        print *,n,err
c        open(unit=11,file='GamHermErrAppend.dat',status='unknown',
c     &                              access='append',form='formatted')
c        write(11,*) n,err
c        close(11)
c      end do
c      return
c      end subroutine testGamHerm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine QuenchedThirringGWerror()
      use rvmodule
      use gaugefield
      use gaugemodule
      use kernelspectrarange
      implicit none
      complex(prc),dimension(4*Nv) :: R
      type(sgnratfunc) :: SRF,ZSRF
      complex(prc) err0,err1,err2,err3,err4
      complex(prc) errA,errB,errC,errD,errE
      integer idx,j
      real(prc) lmin,lmax
c      integer,parameter :: Nerr=8
      integer n,NHT,i,NZ
      
      GAUGETYPE=1 ! 1=compact 2=non-compact
      MDW=one ! domain wall height
      MTYPE=3 ! mass term type
      DWkernel=1 ! 1 for Shamir, 2 for Wilson
      OLTYPE=1 ! 1 is direct calculation, 2 is for DW (K)-type
      baremass=0.0d0
      gbeta=1.0


      print *,seed
      seed=0
      seed(1) = 167868904
      seed(2) = 32712
      seed(3) = seed(1)
      seed(4) = seed(3)
      print *,seed
      call random_seed(put=seed)

      print *,"GAUGETYPE:",GAUGETYPE
      print *,"MDW:",MDW
      print *,"DWkernel:",DWKernel
      print *,"gbeta:",gbeta

      if (GAUGETYPE.eq.1) then
        call makeQuenchedCosineThirringField()
      elseif (GAUGETYPE.eq.2) then
        call makeQuenchedGaussianThirringField()
      endif
      call coef(u,theta)

      call estimateKernelExtrema(10,1,lmax)
      call estimateKernelExtrema(10,-1,lmin)
      print *,lmin,lmax

      call setRVs(Nv*4,R)

      do i=1,40
        NHT=4*i
        NZ=2+2*i
        call setHTcoeffs(NHT,SRF)
        call setZoloCoeffs(NZ,ZSRF,lmin,lmax)
        baremass=0.0
        call GWerrorDOL(R,err0,SRF) ! direct formulation only
        baremass=0.001
        call GWerrorDOL(R,err1,SRF) ! direct formulation only
c        baremass=0.005
c        call GWerrorDOL(R,err2,SRF) ! direct formulation only
        baremass=0.01
        call GWerrorDOL(R,err3,SRF) ! direct formulation only
c        baremass=0.05
c        call GWerrorDOL(R,err4,SRF) ! direct formulation only

        baremass=0.0
        call GWerrorDOL(R,errA,ZSRF) ! direct formulation only
        baremass=0.001
        call GWerrorDOL(R,errB,ZSRF) ! direct formulation only
c        baremass=0.005
c        call GWerrorDOL(R,errC,ZSRF) ! direct formulation only
c        baremass=0.01
c        call GWerrorDOL(R,errD,ZSRF) ! direct formulation only
c        baremass=0.05
c        call GWerrorDOL(R,errE,ZSRF) ! direct formulation only

        print *,NHT,real(err0),real(err1),real(err2),real(err3),
     &              real(err4)
        print *,NZ,real(errA),real(errB),real(errC),real(errD),
     &              real(errE)
        open(unit=11,file='SC_GWHTgb1.dat',status='unknown',
     &                              access='append',form='formatted')
        write(11,*) NHT,real(err0),real(err1),real(err2),real(err3),
     &                  real(err4)
        close(11)
        open(unit=11,file='SC_GWZgb1.dat',status='unknown',
     &                              access='append',form='formatted')
        write(11,*) NZ,real(errA),real(errB),real(errC),real(errD),
     &                  real(errE)
        close(11)
      end do

      return
      end subroutine QuenchedThirringGWerror
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine QuenchedThirringGWcorrection()
      use rvmodule
      use gaugefield
      use gaugemodule
      use kernelspectrarange
      implicit none
      complex(prc),dimension(4*Nv) :: R
      type(sgnratfunc) :: SRF,ZSRF
      complex(prc) err0,err1,err2,err3,err4
      complex(prc) errA,errB,errC,errD,errE
      integer idx,j
      real(prc) lmin,lmax
c      integer,parameter :: Nerr=8
c      integer n,NHT,i,NZ
      
      GAUGETYPE=2 ! 1=compact 2=non-compact
      MDW=9*one/10 ! domain wall height
      MTYPE=3 ! mass term type
      DWkernel=2 ! 1 for Shamir, 2 for Wilson
      OLTYPE=1 ! 1 is direct calculation, 2 is for DW (K)-type
      baremass=0.0d0
      gbeta=0.001


      print *,seed
      seed=0
      seed(1) = 167868904
      seed(2) = 32712
      seed(3) = seed(1)
      seed(4) = seed(3)
      print *,seed
      call random_seed(put=seed)

      print *,"GAUGETYPE:",GAUGETYPE
      print *,"MDW:",MDW
      print *,"DWkernel:",DWKernel
      print *,"gbeta:",gbeta


      if (GAUGETYPE.eq.1) then
        call makeQuenchedCosineThirringField()
      elseif (GAUGETYPE.eq.2) then
        call makeQuenchedGaussianThirringField()
      endif
c      theta=zero
      call coef(u,theta)

      call estimateKernelExtrema(10,1,lmax)
      call estimateKernelExtrema(10,-1,lmin)
      print *,lmin,lmax

      call setRVs(Nv*4,R)
c      R=cone

c      do i=1,40
c        NHT=4*i
c        NZ=2+2*i
c        call setHTcoeffs(NHT,SRF)
        call setZoloCoeffs(30,ZSRF,lmin,lmax)
        baremass=0.0
        call GWcorrectionDOL(R,err0,ZSRF) ! direct formulation only

        print *,gbeta,real(err0)
        open(unit=11,file='GWcor.dat',status='unknown',
     &                              access='append',form='formatted')
        write(11,*) gbeta,real(err0)
        close(11)
c      end do

      return
      end subroutine QuenchedThirringGWcorrection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine gwerrCombo()
      use IOmodule
      use gaugefield
      use gaugemodule
      use spectra
      use kernelspectrarange
      implicit none
      complex(prc),dimension(4*Nv) :: R,DRj,DRjm2
      type(sgnratfunc) :: SRF
      complex(prc) err
      integer idx,j
      real(prc) lmin,lmax,errop

      GAUGETYPE=2 ! 1=compact 2=non-compact
      MTYPE=3
      gbeta=1.75
      baremass=0.0
      dwkernel=2 ! 1=Shamir,2=Wilson

      call makeQuenchedGaussianThirringField()
      call coef(u,theta)
      call setRVs(Nv*4,R)
      
      DRjm2=0
      open(unit=11,file='play.dat',status='unknown',
     &                              access='append',form='formatted')
      do j=10,60,2
        call setHTcoeffs(j,SRF)
        call GWerrorDOL(R,err,SRF) ! direct formulation only
        call DOLop(R,DRj,u,.false.,baremass,SRF)
        errop=maxval(abs(DRj-DRjm2))
        print *,j,real(err),errop
        write(11,*) j,real(err),errop
        DRjm2=DRj
      enddo
      close(11)

      return
      end subroutine gwerrCombo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module gwmodule

