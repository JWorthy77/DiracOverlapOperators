!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! must be included in order of dependence
#include "MODULES/modules.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program dirac
!     use all modules in main program to keep them constantly in scope
      use allmodules
!      use ifport
      implicit none
#ifdef PARALLEL
      include 'mpif.h'
#endif
      character(len=10) :: sdate,stime,szone,edate
!      character(len=10) :: etime,ezone
      integer,dimension(8) :: svals,evals
      real tstart,tend ! cpu seconds
      integer ierr
      integer ntasks
      integer numargs
      character(len=12) arg1,arg2,arg3
      integer iarg1,iarg2,iarg3
      logical RVS_SET

      RVS_SET=.false.

      print *,"Before MPI_INIT"
      goto 4321
      oc_idx=0
      outer_count=0
      ic_idx=0
      inner_count=0

#ifdef PARALLEL    
      call MPI_INIT(ierr)
#endif
      print *,"Start program Dirac"
      call date_and_time(sdate,stime,szone,svals)
      print *,"Date:",sdate,"time:",stime
      call cpu_time(tstart)

#ifdef PARALLEL
      call getMPIinfo()
#ifdef IFORTRAN
      call SLEEPQQ(1000*rank) ! pause each process by 100*rank miliseconds
                              ! to ensure the rv initialises differently
                              ! on each process
#else
      call SLEEP(rank)
#endif


#endif

!!!!!! command line arguments !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      numargs=command_argument_count()
      print *,'numargs:',numargs
      if ((numargs.eq.1).and.(QUENCHED)) then ! this is for calcOverlapRange
        call get_command_argument(1,arg1)
       read (arg1,'(I12)') iarg1
        print *,iarg1
        call initRVs(.false.,.false.,iarg1)
        RVS_SET=.true.
      endif
      if ((numargs.eq.3).and.(.not.QUENCHED)) then ! this is for calcOverlapRange
        call get_command_argument(1,arg1)
        call get_command_argument(2,arg2)
        call get_command_argument(3,arg3)
        read (arg1,'(I12)') iarg1
        read (arg2,'(I5)') iarg2
        print *,iarg1
        print *,iarg2
        call initRVs(.false.,.false.,iarg1)
        RVS_SET=.true.
      endif
      if (numargs.eq.2) then ! this is just for the axial ward calc
        call get_command_argument(1,arg1)
        call get_command_argument(2,arg2)
        print *,arg1,arg2
      endif

      if (.not.RVS_SET) then
        call initRVs(.false.,.false.,0) ! uses a time based seed initialiser
        RVS_SET=.true.
      endif

!!!!!! utility tests !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      call testrvs()
!      call testgrvs()
!      call testRatFuncs
!      call testEigs()
!      call testZolo(5)
!      call testZoloFunctions(11,1d-2,100d0)
!      return

!!!!!! spcify options !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      GAUGETYPE=2 ! 1=compact 2=non-compact
      MDW=one ! domain wall height
      MTYPE=3 ! mass term type
      DWkernel=2 ! 1 for Shamir, 2 for Wilson
      OLTYPE=1 ! 1 is direct calculation, 2 is for DW (K)-type
      baremass=0.0d0
      gbeta=0.8

      print *,"domain wall Ls:",Ls

c      call printOptions()

!      call initSU3()
!      call metropolisSU3(uSU3)
!      print *,"uSU3"
!      print *,uSU3

      call setGammas
      call setPauliMatrices
#ifdef TWODIMENSIONS
      call setIndices2d
#else
      call setIndices
#endif
      call makeGaugeField(.true.)

!!!!!! operator and field tests !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      call testBasicDiracOperators()
!      call testVOLoperators()
!      call testOverlapOperators()
!      call testDomainWall()
!!      call testGammaHermiticity()
!      call testOperatorConvergence()
!      call testOperatorEquivalence()
!      call testOperators()
!      call testFields()
!      call testCondensate()
!      call validateCondensate()
!      call testquenchedauxfields()
!!       call testweyl()
!!       call testolweyl()

!!!!! calculate locality and GW error !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      call calcLocality()
!      call calcQuenchedLocality()
!      call testLocality(5,4,1)
!      call testConvertedLocality(5,6)
!      call ThirringGWerror()
!      call gwerrCombo()
!      call QuenchedThirringGWerror()
!      call QuenchedThirringGWcorrection()

!!!! calculate eigenvalue extrema !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      call calcKernelRange() ! calculate max/min eigs of kernel only
!      call calcQuenchedKernelRange() ! calculate max/min eigs of kernel only
!      call calcOverlapRange(iarg2) ! calculate max/min eigs of overlap only
!      call calcWOLExtrema()

!      call getThirringSpectra()
!      call evalSpectra()
!      call RRspectraRJN() ! calculate minimum eigenvalues using Rajamanis code
!      call calcRayleighRitzEigs() 

#ifdef PARALLEL

!      ntasks=0 ! ntasks calculated from input
!      call parallelSequence(ntasks,2) ! condensate
!      ntasks=0 ! ntasks calculated from input
!      call parallelSequence(ntasks,3) ! condensate slightly differently organised
!      ntasks=0 ! ntasks calculated from input
!      call parallelSequence(ntasks,4) ! overlap spectra range 
!      ntasks=0 ! ntasks calculated from input
!      call parallelSequence(ntasks,5) ! quenched condensate
!      ntasks=100
!      call parallelSequence(ntasks,1) ! spectra
#endif 

!!!! calculate condensate !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      call measureCondensate()
!      call measureCondensateDomWall()


!!!! meson propagators !!!!!!!!!!!!!!!

c      call calcQuenchedMesonPropagator(10)
!      call calcQuenchedAxialWard()
!      read (arg1,'(I5)') iarg1
!      read (arg2,'(I5)') iarg2
!      call calcAxialWard(iarg1,iarg2)

!!!!! old and untidied !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c      call readThirringConFile(500,theta)
c      call getThirringSpectra()
c      call measureCondensateGN(1000)


c      call testStagGN()

c       call test()

c      call testVOL(.false.,8)
c      call testVOL(.true.,2)
c      call testConv(.true.)
c      call testOverlap(10)
c      call testOverlap(9)
c      call setRVs(Nv*4,RR)    
c      call testOLAction(12,RR)
c      call testOLGaugeSymmetry()
c      call testSPF()
c      call testPF()
c      call testHMCSPF()
c      call testHMCPF()
c      call testDSOperators()
c      call testOLShamirGaugeSymmetry()
c      call testMultigrid()

c      call evalCondensateDW(u,baremass,1)
c      call evalCondensateDW(u,baremass,2)
c      call testNewSGN(.false.)

c      call writeGaugeField('F1.dat')     
c      call evalCondensate(u,zero,pbp)
c      call makeDynamicGaugeField()     

c      call evalMesonMass()

c      call testGaugeFieldCorrelation()
c      call testLocality(5,1,1)
c      call testLocality(5,3,160)
c      MTYPE=2
c      call testConvertedLocality(5,6)

c       call extractGaugeLines(5)

c       call testGWerror()
c       call ThirringGWerror()
c       call measureCondensateGW()
c       call testGamHerm()

c      call MGLs(32)

c       call getSpectra()
c      call evalSpectra()

c      call setHTcoeffs(12,SRF)
c      call calcOLDeterminant(u,.false.,baremass,SRF)
c      call calcKDWDeterminant(u,.false.,baremass,detK)
c      call calcDWDeterminant(u,.false.,baremass,detDW)
c      MTYPE=1
c      call calcDWDeterminant(u,.false.,one,detDW1)
c      print *,detK,detDW,detDW1,detDW/detDW1
c      DWkernel=2 ! 1 for Shamir, 2 for Wilson
c      call calcOLDeterminant(u,.false.,baremass,SRF)

c      call evalSmearedCondensateDW()
c      call evalCondensateDDW(u,pbp)
c      print *,pbp
#ifdef PARALLEL
      print *,'rank',rank,'before MPI_FINALZE'
      print *,rank,"End program Dirac"
      call date_and_time(sdate,stime,szone,svals)
      print *,rank,"Date:",sdate,"time:",stime
      call cpu_time(tend)
      print *,rank,"CPU seconds:",tend-tstart,"CPU mins:",(tend-tstart)/
     &60,"CPU hrs:",(tend-tstart)/3600
      call MPI_FINALIZE()
#else
      print *,"End program Dirac"
      call date_and_time(sdate,stime,szone,svals)
      print *,"Date:",sdate,"time:",stime
      call cpu_time(tend)
      print *,"CPU seconds:",tend-tstart,"CPU mins:",(tend-tstart)/
     &60,"CPU hrs:",(tend-tstart)/3600
#endif
4321  continue
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testOperators()
!     test different Dirac operators
      use testgamhermmod
      use testDconvmod
      use testDequivmod
      implicit none

!      call testVOL(.false.,11)
!      call testGammaHermiticity()
!      call testOperatorConvergence()
      call testOperatorEquivalence()

      return
      end subroutine testOperators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine test()
      use numbers
      use rvmodule
      use ratfuncs
      use options
      use gammas
      use gaugefield
      use overlapmoduledev
      use domainwallmod
      use basicdiracopsmod
      implicit none
      complex(prc) R(Nv,4),TMP(Nv,4)
      complex(prc) DR1(Nv,4),DR2(Nv,4),DR3(Nv,4) 
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()

      MTYPE=1
      call setRVs(Nv*4,R)     
      call setHTcoeffs(Ls,SRF)
      DWkernel=2
      print *,"DAGGER=.FALSE."
      call KDDW4(R,DR3,u,.false.,baremass)
      call DOverlap(R,DR1,u,.false.,baremass,SRF)
      print *,"DOLW-KDDW4",maxval(abs(DR3-DR1))
c      DWkernel=1
c      call KDDW4(R,DR3,u,.false.,baremass)
c      call DOLS(R,DR2,u,.false.,baremass,SRF)
c      print *,"DOLS-KDDW4",maxval(abs(DR3-DR2))
c      print *,"DAGGER=.TRUE."
c      call KDDW4(R,DR3,u,.true.,baremass)
c      call DOverlap(R,DR1,u,.true.,baremass,SRF)
c      print *,"DOLW-KDDW4",maxval(abs(DR3-DR1))
c      DWkernel=1
c      call KDDW4(R,DR3,u,.true.,baremass)
c      call DOLS(R,DR2,u,.true.,baremass,SRF)
c      print *,"DOLS-KDDW4",maxval(abs(DR3-DR2))

c      TMP=R
c      call mGmu(TMP,4)
c      call KDDW4(TMP,DR1,u,.false.,baremass)
c      call mGmu(DR1,4)
c      call KDDW4(R,DR2,u,.true.,baremass)
c      print *,"G5 KDDW4 Hermitian",maxval(abs(DR1-DR2))
      
      TMP=R
      call mGmu(TMP,5)
      call DWilson(TMP,DR1,u,.false.,baremass,cone)
      call mGmu(DR1,5)
      call DWilson(R,DR2,u,.true.,baremass,cone)
      print *,"G5 Wilson Hermitian",maxval(abs(DR1-DR2))

      TMP=R
      call mGmu(TMP,5)
      call DWilson(TMP,DR1,u,.false.,baremass,czero)
      call DWilson(DR1,DR2,u,.true.,baremass,czero)
      call mGmu(DR2,5)
      call DWilson(R,TMP,u,.true.,baremass,czero)
      call DWilson(TMP,DR1,u,.false.,baremass,czero)
      print *,"G5 Wilson DdagD Hermitian",maxval(abs(DR1-DR2))
      
      
      Mptr => DdagDpC
      TMP=R
      call mGmu(TMP,5)
      call DdagDpC(TMP,DR1,u,.false.,baremass,czero)
      call mGmu(DR1,5)
      call DdagDpC(R,DR2,u,.true.,baremass,czero)
      print *,"G5 DdagDpC Hermitian",maxval(abs(DR1-DR2))

      Mptr => DdagDpC
      TMP=R
      call mGmu(TMP,5)
      call IM(TMP,DR1,u,.false.,-MDW,czero,Mptr)
      call mGmu(DR1,5)
      call IM(R,DR2,u,.true.,-MDW,czero,Mptr)
      print *,"G5 IM Hermitian",maxval(abs(DR1-DR2))


      Mptr => VOLpf
      TMP=R
      call mGmu(TMP,5)
      call Mptr(TMP,DR1,u,.false.,-MDW,SRF)
      call mGmu(DR1,5)
      call Mptr(R,DR2,u,.true.,-MDW,SRF)
      print *,"G5 VOLpf Hermitian",maxval(abs(DR1-DR2))

      Mptr => DOverlap
      TMP=R
      call mGmu(TMP,5)
      call Mptr(TMP,DR1,u,.false.,-MDW,SRF)
      call mGmu(DR1,5)
      call Mptr(R,DR2,u,.true.,-MDW,SRF)
      print *,"G5 DOverlap Hermitian",maxval(abs(DR1-DR2))


      return
      end subroutine test
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
