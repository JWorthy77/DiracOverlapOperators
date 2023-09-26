      module testcondmod
      use pacc
      use arraysizes
      use numbers
      use options
      use ratfuncs
      use gaugefield
      use gaugemodule
      use IOmodule
      use condpointtoolsmod
      use condensatemodule
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testCondensate()
!     when comparing with notes, rem MDW=0.9 for analytic solution
      implicit none
      integer dop
      complex(prc) pbpM1S,pbpM3S,pbpM1W,pbpM3W

!     free field point tests
      print *,"free field point tests"
      call makeGaugeField(.true.)
      call coef(u,theta)

      baremass=0.01
!      dop=2 ! dirac operator - 1 overlap, 2 domain wall

      ! test free field Wilson point values
      call testCond(2,2,1,.false.) ! DomWall, Wilson, M1, POINT
      call testCond(2,2,3,.false.) ! DomWall, Wilson, M3, POINT
!      call testCond(1,2,1,.false.) ! Overlap, Wilson, M1, POINT
!      call testCond(1,2,3,.false.) ! Overlap, Wilson, M3, POINT
!      call testCond(1,2,4,.false.) ! Overlap, Wilson, M4, POINT
      ! test free field Wilson noisy values
!      call testCond(2,2,1,.true.) ! DomWall, Wilson, M1, NOISY
!      call testCond(2,2,3,.true.) ! DomWall, Wilson, M3, NOISY
!      call testCond(1,2,1,.true.) ! Overlap, Wilson, M1, NOISY
!      call testCond(1,2,3,.true.) ! Overlap, Wilson, M3, NOISY
!      call testCond(1,2,4,.true.) ! Overlap, Wilson, M4, NOISY
      ! test free field Shamir point values
      call testCond(2,1,1,.false.) ! DomWall, Shamir, M1, POINT
      call testCond(2,1,3,.false.) ! DomWall, Shamir, M3, POINT
!      call testCond(1,1,1,.false.) ! Overlap, Shamir, M1, POINT
!      call testCond(1,1,3,.false.) ! Overlap, Shamir, M3, POINT
      ! test free field Shamir noisy values
!      call testCond(2,1,1,.true.) ! DomWall, Shamir, M1, NOISY
!      call testCond(2,1,3,.true.) ! DomWall, Shamir, M3, NOISY
!      call testCond(1,1,1,.true.) ! Overlap, Shamir, M1, NOISY
!      call testCond(1,1,3,.true.) ! Overlap, Shamir, M3, NOISY


c      call testCond_DomWallPoint()    
c      call testCond_OLPoint()

!     non-free field noisy tests
!      call setGRVs(3*Nv,theta)
!      theta=theta/4
!      call coef(u,theta)
!      print *,"non-free field noisy tests - theta = grv/4"
!      call testCondensateDW_OL(1,1,.true.) ! Shamir, M1, NOISY 
!      call testCondensateDW_OL(1,3,.true.) ! Shamir, M3, NOISY 
!      call testCondensateDW_OL(2,1,.true.) ! Wilson, M1, NOISY
!      call testCondensateDW_OL(2,3,.true.) ! Wilson, M3, NOISY 

      return
      end subroutine testCondensate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testCond(dop,idk,idm,NOISY)
      use statsmod
      implicit none
      integer,intent(in) :: dop,idk,idm
      logical,intent(in) :: NOISY

      print *," "
      if (dop.eq.1) then ! dirav operator is overlap
        call testCondensateOverlap(idk,idm,NOISY)
      elseif(dop.eq.2) then ! dirac operator is domain wall
        call testCondensateDomWall(idk,idm,NOISY)
      endif
      return
      end subroutine testCond
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine testCond_DomWallPoint()
c!     test condensate using 6x6 meshes with free field
c!     when comparing with notes, rem MDW=0.9 for analytic solution
c      implicit none
c      type(sgnratfunc) :: SRF
c      complex(prc) pbpM1S,pbpM3S,pbpM1W,pbpM3W
c
c      print *,"test DDW point condensate on free field"
c
c      MTYPE=1
c      dwkernel=1
c      call evalCondPoint_DDW_Shamir(u,pbpM1S,1)
c      print *,"DDW pbp M1 Shamir:",pbpM1S
c
c      MTYPE=3
c      dwkernel=1
c      call evalCondPoint_DDW_Shamir(u,pbpM3S,1)
c      print *,"DDW pbp M3 Shamir:",pbpM3S
c
c      MTYPE=1
c      dwkernel=2
c      call evalCondPoint_DDW_Wilson(u,pbpM1W,1)
c      print *,"DDW pbp M1 Wilson:",pbpM1W
c
c      MTYPE=3
c      dwkernel=2
c      call evalCondPoint_DDW_Wilson(u,pbpM3W,1)
c      print *,"DDW pbp M3 Wilson:",pbpM3W
c
c
c      print *,"DDW pbp M1 Shamir:",pbpM1S
c      print *,"DDW pbp M3 Shamir:",pbpM3S
c      print *,"DDW pbp M1 Wilson:",pbpM1W
c      print *,"DDW pbp M3 Wilson:",pbpM3W
c
c      return
c      end subroutine testCond_DomWallPoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     subroutine testCond_OLPoint()
c!     test condensate using 6x6 meshes with free field
c!     when comparing with notes, rem MDW=0.9 for analytic solution
c      implicit none
c      type(sgnratfunc) :: SRF
c      complex(prc) pbpM1S,pbpM4S,pbpM1W,pbpM4W
c
c      print *,"test OL point condensate on free field"
c
c      call setHTcoeffs(Ls,SRF)
c      MTYPE=1
c      dwkernel=1
c      call evalCondPoint_OL(u,pbpM1S,1,SRF)
c      print *,"OL pbp M1 Shamir:",pbpM1S
c
c      MTYPE=4
c      dwkernel=1
c      call evalCondPoint_OL(u,pbpM4S,1,SRF)
c      print *,"OL pbp M4 Shamir:",pbpM4S
c
c      MTYPE=1
c      dwkernel=2
c      call evalCondPoint_OL(u,pbpM1W,1,SRF)
c      print *,"OL pbp M1 Wilson:",pbpM1W
c
c      MTYPE=4
c      dwkernel=2
c      call evalCondPoint_OL(u,pbpM4W,1,SRF)
c      print *,"OL pbp M4 Wilson:",pbpM4W
c
c      print *,"OL pbp M1 Shamir:",pbpM1S
c      print *,"OL pbp M4 Shamir:",pbpM4S
c      print *,"OL pbp M1 Wilson:",pbpM1W
c      print *,"OL pbp M4 Wilson:",pbpM4W
c
c      return
c      end subroutine testCond_OLPoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testCondensateDomWall(idk,idm,NOISY)
      use statsmod
      implicit none
      integer,intent(in) :: idk,idm
      logical,intent(in) :: NOISY
      type(sgnratfunc) :: SRF
      complex(prc) pbp,pbpA,pbpB,meanA,meanB,sdA,sdB
      complex(prc) pbpCompsA(1000),pbpCompsB(1000)
      integer idn,knoise

      MTYPE=idm
      dwkernel=idk
      print *,"test Condensate (Domain Wall):"
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"Ls:",Ls
      print *,"baremass:",baremass
      print *,"GAUGETYPE",GAUGETYPE

      call setHTcoeffs(Ls,SRF)
      if (NOISY) then
        pbpA=0
        pbpB=0
        knoise=1000
        do idn=1,knoise
          if (dwkernel.eq.1) then
            call evalCondNoisy_DomWall_Shamir(u,pbp)
          elseif (dwkernel.eq.2) then
            call evalCondNoisy_DomWall_Wilson(u,pbp)
          endif
          pbpCompsB(idn)=pbp
          pbpB=pbpB+pbp
          call calcVar(idn,pbpCompsB(1:idn),meanB,sdB)
          print *,"DW pbp:",idn,meanB,sdB,sdB/sqrt(real(idn))
        enddo
        pbpB=pbpB/knoise
      elseif(.not.NOISY) then

        if (dwkernel.eq.1) then
          call evalCondPoint_DDW_Shamir(u,pbpB,1)
        elseif (dwkernel.eq.2) then
          call evalCondPoint_DDW_Wilson(u,pbpB,1)
        endif

      endif

      print *,"DomWall pbp:",pbpB

      return
      end subroutine testCondensateDomWall
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testCondensateOverlap(idk,idm,NOISY)
      use statsmod
      implicit none
      integer,intent(in) :: idk,idm
      logical,intent(in) :: NOISY
      type(sgnratfunc) :: SRF
      complex(prc) pbp,pbpA,pbpB,meanA,meanB,sdA,sdB
      complex(prc) pbpComps(1000)
      integer idn,knoise

      MTYPE=idm
      dwkernel=idk
      print *,"test Condensate (Overlap):"
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"Ls:",Ls
      print *,"baremass:",baremass
      print *,"GAUGETYPE",GAUGETYPE

      call setHTcoeffs(Ls,SRF)
      if (NOISY) then
        pbpA=0
        pbpB=0
        knoise=1000
        do idn=1,knoise
          call evalCondNoisy_OL(u,pbp,SRF)
          pbpComps(idn)=pbp
          pbpB=pbpB+pbp
          call calcVar(idn,pbpComps(1:idn),meanB,sdB)
          print *,"OL pbp:",idn,meanB,sdB,sdB/sqrt(real(idn))
        enddo
        pbpB=pbpB/knoise
      elseif(.not.NOISY) then

        call evalCondPoint_OL(u,pbpB,1,SRF)

      endif

      print *,"Overlap pbp:",pbpB

      return
      end subroutine testCondensateOverlap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testCondensateDW_OL(idk,idm,NOISY)
      use statsmod
!     test condensate using 6x6 meshes with free field
!     when comparing with notes, rem MDW=0.9 for analytic solution
      implicit none
      integer,intent(in) :: idk,idm
      logical,intent(in) :: NOISY
      type(sgnratfunc) :: SRF
      complex(prc) pbp,pbpA,pbpB,meanA,meanB,sdA,sdB
      complex(prc) pbpCompsA(1000),pbpCompsB(1000)
      integer idn

      MTYPE=idm
      dwkernel=idk
      print *,"test Condensate:"
      print *,"MTYPE:",MTYPE
      print *,"dwkernel:",dwkernel
      print *,"Ls:",Ls
      print *,"baremass:",baremass
      print *,"GAUGETYPE",GAUGETYPE

      call setHTcoeffs(Ls,SRF)
      if (NOISY) then
        pbpA=0
        pbpB=0
        do idn=1,100
          if (dwkernel.eq.1) then
            call evalCondNoisy_DomWall_Shamir(u,pbp)
          elseif (dwkernel.eq.2) then
            call evalCondNoisy_DomWall_Wilson(u,pbp)
          endif
          pbpCompsB(idn)=pbp
          pbpB=pbpB+pbp
          call evalCondNoisy_OL(u,pbp,SRF)
          pbpCompsA(idn)=pbp
          pbpA=pbpA+pbp
c          print *,"OL pbp component:",idn,pbpCompsA(1:idn)
c          print *,"DW pbp component:",idn,pbpCompsB(1:idn)
          call calcVar(idn,pbpCompsA(1:idn),meanA,sdA)
          call calcVar(idn,pbpCompsB(1:idn),meanB,sdB)
          print *,"OL pbp:",idn,meanA,sdA,sdA/sqrt(real(idn))
          print *,"DW pbp:",idn,meanB,sdB,sdB/sqrt(real(idn))
        enddo
        pbpA=pbpA/1000
        pbpB=pbpB/1000
      elseif(.not.NOISY) then

        if (dwkernel.eq.1) then
          call evalCondPoint_DDW_Shamir(u,pbpB,1)
        elseif (dwkernel.eq.2) then
          call evalCondPoint_DDW_Wilson(u,pbpB,1)
        endif
        if (MTYPE.eq.3) then
          MTYPE=4 ! since this is the equivalent version for shamir M3
        endif
        call evalCondPoint_OL(u,pbpA,1,SRF)
        if (MTYPE.eq.4) then
          MTYPE=3
        endif

      endif

      print *,"OL pbp:",pbpA
      print *,"DW pbp:",pbpB

      return
      end subroutine testCondensateDW_OL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testGaugeFieldCorrelation()
      implicit none
      integer i,t,c
      integer,parameter :: Nth=90
      complex(prc) val(Nth),vav,var,fluc(Nth)
      complex(prc) cor(20)
      complex(prc) pbp
      character(len=80) fname 
      
c      call setRVs(Nv*4,R)
      do i=1,Nth
        call readConvertedThetaFileName(i,fname)
        call readConvertedThirringGaugeField(fname,theta)
        call coef(u,theta)
c        call evalCondensateDW(u,baremass,2,pbp)
        val(i)=pbp
      enddo
      open(unit=10,file='gcor.dat',status='unknown')
      vav=zero
      do i=1,Nth
        vav=vav+val(i)
      end do
      vav=vav/Nth
      var=zero
      do i=1,Nth
        var=var+(val(i)-vav)*(val(i)-vav)
        fluc(i)=(val(i)-vav)*conjg(val(i)-vav)
      end do
      var=var/(Nth-one)

      cor=czero
      do t=1,20
        do c=1,70
          cor(t)=cor(t)+fluc(c)*fluc(c+t)
        end do
        cor(t)=cor(t)/90
      enddo
      print *,"correlation"
      print *,cor

      do i=1,Nth
        write(10,*) i,val(i),var
      end do
      close(10)
    
      return
      end subroutine testGaugeFieldCorrelation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      subroutine validateCondensate()
c      use kernelspectrarange
c      use statsmod
c!     when comparing with notes, rem MDW=0.9 for analytic solution
c!     and use 8x8 for thesis cases
c      implicit none
c      real(prc) lmin,lmax
c      integer j,jmin,jmax,shft
c      integer idn,is,LsTest
c      type(sgnratfunc) :: SRF
c      character(len=80) :: fname,fname2
c      integer,parameter :: knoise=40
c      real div
c      real :: diff1,diff3
c      complex(prc) pbpHM1S,pbpHM3S,pbpHM1W,pbpHM3W
c      complex(prc) pbpZM1S,pbpZM3S,pbpZM1W,pbpZM3W
c      complex(prc),dimension(knoise) :: pHM1Ss,pHM3Ss,pHM1Ws,pHM3Ws
c      complex(prc),dimension(knoise) :: pZM1Ss,pZM3Ss,pZM1Ws,pZM3Ws
c      complex(prc) pbpHM1Sav,pbpHM3Sav,pbpHM1Wav,pbpHM3Wav
c      complex(prc) pbpZM1Sav,pbpZM3Sav,pbpZM1Wav,pbpZM3Wav
c      complex(prc) pbpHM1Ssd,pbpHM3Ssd,pbpHM1Wsd,pbpHM3Wsd
c      complex(prc) pbpZM1Ssd,pbpZM3Ssd,pbpZM1Wsd,pbpZM3Wsd
c      logical,parameter :: SHAMIRANALYTIC=.false. ! Ls conv of OL Shamir on 8x8 free aux field M1/M3,HT/ZOLO (point)
c      logical,parameter :: WILSONANALYTIC=.false. ! Ls conv of OL Wilson on 8x8 free aux field M1/M3,HT/ZOLO (point)
c      logical,parameter :: SHAMIRNOISY=.false. ! Ls conv of OL Shamir on 8x8 random aux field M1/M3,HT/ZOLO (noisy)
c      logical,parameter :: WILSONNOISY=.true. ! Ls conv of OL Wilson on 8x8 random aux field M1/M3,HT/ZOLO (noisy)
c      logical,parameter :: DOMWALLSHAMIR=.false. ! Equivalence of point Shamir domwall and OL, free aux, rand aux, M1/M3 HT
c      logical,parameter :: DOMWALLWILSON=.false. ! Equivalence of point Shamir domwall and OL, free aux, rand aux, M1/M3 HT
c
c!     analytic with MDW=0.9
c      call makeGaugeField(.true.)
c      call coef(u,theta)
c      baremass=0.01
c      MDW=0.9  ! approximation eig(Hs)=eig(Hw)/[2+eig(Hw)] only holds for MDW=1.0
c
c      if (SHAMIRANALYTIC) then
c
c        dwkernel=1 ! Shamir
c        call getMinMaxHEigs(lmin,lmax)
c        print *,"lmin:",lmin
c        print *,"lmax:",lmax
c        fname="EigShamir8x8.dat"
c        open(unit=11,file=fname,status='unknown',form='formatted')
c        write(11,'(2E13.6)') lmin,lmax
c        close(11)
c
c        jmin=3
c        jmax=10
c        fname="FreeShamir8x8.dat"
c        open(unit=11,file=fname,status='unknown',form='formatted')
c        do j=jmin,jmax
c          call setHTcoeffs(j,SRF)
c          MTYPE=1
c          call evalCondPoint_OL(u,pbpHM1S,1,SRF)
c          MTYPE=3
c          call evalCondPoint_OL(u,pbpHM3S,1,SRF)
c
c          call setZoloCoeffs(j,SRF,lmin,lmax)
c          MTYPE=1
c          call evalCondPoint_OL(u,pbpZM1S,1,SRF)
c          MTYPE=3
c          call evalCondPoint_OL(u,pbpZM3S,1,SRF)
c          print *,"pbp M1/M3:",j,pbpHM1S,pbpHM3S
c          write(11,'(I6,4E19.10)') j,real(pbpHM1S),real(pbpHM3S),
c     &                             real(pbpZM1S),real(pbpZM3S)
c        end do
c        close(11)
c
c      endif
c
c      if (WILSONANALYTIC) then
c
c        dwkernel=2 ! Wilson
c        call getMinMaxHEigs(lmin,lmax)
c        print *,"lmin:",lmin
c       print *,"lmax:",lmax
c        fname="EigWilson8x8.dat"
c        open(unit=11,file=fname,status='unknown',form='formatted')
c        write(11,'(2E13.6)') lmin,lmax
c        close(11)
c
c        jmin=3
c        jmax=10
c        fname="FreeWilson8x8.dat"
c        open(unit=11,file=fname,status='unknown',form='formatted')
c        do j=jmin,jmax
c          call setHTcoeffs(j,SRF)
c          MTYPE=1
c          call evalCondPoint_OL(u,pbpHM1S,1,SRF)
c          MTYPE=3
c          call evalCondPoint_OL(u,pbpHM3S,1,SRF)
c
c          call setZoloCoeffs(j,SRF,lmin,lmax)
c          MTYPE=1
c          call evalCondPoint_OL(u,pbpZM1S,1,SRF)
c          MTYPE=3
c          call evalCondPoint_OL(u,pbpZM3S,1,SRF)
c          print *,"pbp M1/M3:",j,pbpHM1S,pbpHM3S
c          write(11,'(I6,4E19.10)') j,real(pbpHM1S),real(pbpHM3S),
c     &                             real(pbpZM1S),real(pbpZM3S)
c        end do
c        close(11)
c
c      endif
c
c
c!     noisy with MDW=1.0
c      if (seedsize.ne.8) then
c        print *,"rvs not properly initialised ... system has seedsize !=
c     &8"
c        print *,"seedsize:",seedsize   ! will need to check on other systems
c        stop
c      endif
c     
c      seed(1)= -668831433
c      seed(2)= -312266119
c      seed(3)=  1666396292
c      seed(4)= -109198905
c      seed(5)= -1566517286
c      seed(6)= -1589941339
c      seed(7)= -1136709282
c      seed(8)=  860140008
c
c      call RANDOM_SEED(put=seed)
c      call setGRVs(3*Nv,theta)
c!      theta=theta
c      call coef(u,theta)
c      baremass=0.01
c      MDW=1.0  ! approximation eig(Hs)=eig(Hw)/[2+eig(Hw)] only holds for MDW=1.0
c
c      if (SHAMIRNOISY) then
c
c        dwkernel=1 ! Shamir
c        call getMinMaxHEigs(lmin,lmax)
c        print *,"lmin:",lmin
c        print *,"lmax:",lmax
c        fname="EigShamirNoisy8x8.dat"
c        open(unit=11,file=fname,status='unknown',form='formatted')
c        write(11,'(2E13.6)') lmin,lmax
c        close(11)
c
c        jmin=6
c        jmax=20
c        fname="NoisyShamir8x8.dat"
c        open(unit=11,file=fname,status='unknown',form='formatted')
c        do j=jmin,jmax,2
c          do idn=1,knoise
c            call setHTcoeffs(j,SRF)
c            MTYPE=1
c            call evalCondNoisy_OL(u,pbpHM1S,SRF)
c            pHM1Ws(idn)=pbpHM1W
c            MTYPE=3
c            call evalCondNoisy_OL(u,pbpHM3S,SRF)
c            pHM3Ws(idn)=pbpHM3W
c
c            call setZoloCoeffs(j,SRF,lmin,lmax)
c            MTYPE=1
c            call evalCondNoisy_OL(u,pbpZM1S,SRF)
c            pZM1Ws(idn)=pbpZM1W
c            MTYPE=3
c            call evalCondNoisy_OL(u,pbpZM3S,SRF)
c            pZM3Ws(idn)=pbpZM3W
c            print *,"pbp M1/M3:",j,pbpHM1S,pbpHM3S
c          end do
c          call calcVar(knoise,pHM1Ss,pbpHM1Sav,pbpHM1Ssd)
c          call calcVar(knoise,pHM3Ss,pbpHM3Sav,pbpHM3Ssd)
c          call calcVar(knoise,pZM1Ss,pbpZM1Sav,pbpZM1Ssd)
c          call calcVar(knoise,pZM3Ss,pbpZM3Sav,pbpZM3Ssd)
c          div=sqrt(real(knoise))
c          write(11,'(I6,8E13.6)') j,real(pbpHM1Sav),real(pbpHM3Sav),
c     &                              real(pbpZM1Sav),real(pbpZM3Sav),
c     &                      real(pbpHM1Ssd)/div,real(pbpHM3Ssd)/div,
c     &                      real(pbpZM1Ssd)/div,real(pbpZM3Ssd)/div
c        end do
c        close(11)
c      endif
c
c
c
c      if (WILSONNOISY) then
c
c        dwkernel=2 ! Wilson
c        call getMinMaxHEigs(lmin,lmax)
c        print *,"lmin:",lmin
c        print *,"lmax:",lmax
c        fname="EigWilsonNoisy8x8.dat"
c        open(unit=11,file=fname,status='unknown',form='formatted')
c        write(11,'(2E13.6)') lmin,lmax
c        close(11)
c
c        jmin=6
c        jmax=20
c        fname="NoisyWilson8x8.dat"
c        open(unit=11,file=fname,status='unknown',form='formatted')
c        do j=jmin,jmax,2
c         fname2="Noise"//trim(itoa(j))//".dat"
c          print *,"fname2: ",fname2
c          open(unit=14,file=fname2,status='unknown',form='formatted')
c          do idn=1,knoise
c            call setHTcoeffs(j,SRF)
c            MTYPE=1
c            call evalCondNoisy_OL(u,pbpHM1W,SRF)
c            pHM1Ws(idn)=pbpHM1W
c            MTYPE=3
c            call evalCondNoisy_OL(u,pbpHM3W,SRF)
c            pHM3Ws(idn)=pbpHM3W
c
c            call setZoloCoeffs(j,SRF,lmin,lmax)
c            MTYPE=1
c            call evalCondNoisy_OL(u,pbpZM1W,SRF)
c            pZM1Ws(idn)=pbpZM1W
c            MTYPE=3
c            call evalCondNoisy_OL(u,pbpZM3W,SRF)
c            pZM3Ws(idn)=pbpZM3W
c            print *,"pbp M1/M3:",j,idn,pbpHM1W,pbpHM3W
c            print *,"pbp M1/M3:",j,idn,pbpZM1W,pbpZM3W
c            write(14,'(2I6,4E17.6)') j,idn,real(pbpHM1W),real(pbpHM3W),
c     &                                     real(pbpZM1W),real(pbpZM3W)
c          end do
c          close(14)
c          call calcVar(knoise,pHM1Ws,pbpHM1Wav,pbpHM1Wsd)
c          call calcVar(knoise,pHM3Ws,pbpHM3Wav,pbpHM3Wsd)
c          call calcVar(knoise,pZM1Ws,pbpZM1Wav,pbpZM1Wsd)
c          call calcVar(knoise,pZM3Ws,pbpZM3Wav,pbpZM3Wsd)
c          div=sqrt(real(knoise))
c          write(11,'(I6,8E13.4)') j,real(pbpHM1Wav),real(pbpHM3Wav),
c     &                              real(pbpZM1Wav),real(pbpZM3Wav),
c     &                      real(pbpHM1Wsd)/div,real(pbpHM3Wsd)/div,
c     &                      real(pbpZM1Wsd)/div,real(pbpZM3Wsd)/div
c        end do
c        close(11)
c      endif
c
c
c      LsTest=16
c      if (Ls.ne.LsTest) then
c        print *,"recompile with Ls=16"
c        stop
c      endif
c      call makeGaugeField(.true.)
c      call coef(u,theta)
c      baremass=0.01
c      MDW=1.0  
c      if (DOMWALLWILSON) then ! Equivalence of point Wilson domwall and OL, free aux, rand aux, M1/M3 HT
c          dwkernel=2 ! Wilson
c          call setHTcoeffs(LsTest,SRF)
c          MTYPE=1
c          call evalCondPoint_OL(u,pbpHM1W,1,SRF)
c          print *,"Wilson OL M1:",pbpHM1W
c          diff1=abs(real(pbpHM1W))
c          call evalCondPoint_DDW_Wilson(u,pbpHM1W,1)
c          diff1=diff1-abs(real(pbpHM1W))
c          print *,"Wilson DW M1:",pbpHM1W
c          MTYPE=3
c          call evalCondPoint_OL(u,pbpHM3W,1,SRF)
c          print *,"Wilson OL M3:",pbpHM3W
c          diff3=abs(real(pbpHM3W))
c          call evalCondPoint_DDW_Wilson(u,pbpHM3W,1)
c          diff3=diff3-abs(real(pbpHM3W))
c          print *,"M1/M3 difference:",diff1,diff3
c      end if
c
c      return
c      end subroutine validateCondensate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testcondmod
