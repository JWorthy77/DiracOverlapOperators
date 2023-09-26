      module validatecondmodule
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
      integer,parameter :: knoise=80
      character(len=80) :: fname,fname2
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine validateCondensate()
      use kernelspectrarange
      use statsmod
!     when comparing with notes, rem MDW=0.9 for analytic solution
!     and use 8x8 for thesis cases
      implicit none
      logical,parameter :: SHAMIRANALYTIC=.false. ! Ls conv of OL Shamir on 8x8 free aux field M1/M3,HT/ZOLO (point)
      logical,parameter :: WILSONANALYTIC=.false. ! Ls conv of OL Wilson on 8x8 free aux field M1/M3,HT/ZOLO (point)
      logical,parameter :: SHAMIRNOISYANALYTIC=.false. ! Ls conv of OL Shamir on 8x8 free aux field M1/M3,HT/ZOLO (noisy)
      logical,parameter :: WILSONNOISYANALYTIC=.false. ! Ls conv of OL Wilson on 8x8 free aux field M1/M3,HT/ZOLO (noisy)
      logical,parameter :: SHAMIRNOISY=.false. ! Ls conv of OL Shamir on 8x8 random aux field M1/M3,HT/ZOLO (noisy)
      logical,parameter :: WILSONNOISY=.true. ! Ls conv of OL Wilson on 8x8 random aux field M1/M3,HT/ZOLO (noisy)
      logical,parameter :: DOMWALLSHAMIR=.false. ! Equivalence of point Shamir domwall and OL, free aux, rand aux, M1/M3 HT
      logical,parameter :: DOMWALLWILSON=.false. ! Equivalence of point Shamir domwall and OL, free aux, rand aux, M1/M3 HT

!     analytic with MDW=0.9
      call makeGaugeField(.true.)
      call coef(u,theta)
      baremass=0.01
      MDW=0.9  ! approximation eig(Hs)=eig(Hw)/[2+eig(Hw)] only holds for MDW=1.0

      if (SHAMIRANALYTIC) then
        call validateOLAnalyticPoint(1)
      endif

      if (WILSONANALYTIC) then
        call validateOLAnalyticPoint(2)
      endif

      if (SHAMIRNOISYANALYTIC) then
        call validateOLNoisy(1,3,10,1)
      endif

      if (WILSONNOISYANALYTIC) then
        call validateOLNoisy(2,3,10,1) 
      endif


!     noisy with MDW=1.0
      if ((seedsize.ne.8).and.(seedsize.ne.12)) then
        print *,"rvs not properly initialised ... system has seedsize !=
     &8"
        print *,"seedsize:",seedsize   ! will need to check on other systems
        stop
      endif
     
      seed(1)= -668831433
      seed(2)= -312266119
      seed(3)=  1666396292
      seed(4)= -109198905
      seed(5)= -1566517286
      seed(6)= -1589941339
      seed(7)= -1136709282
      seed(8)=  860140008
      if (seedsize.eq.12) then
      seed(9)= -1566517286
      seed(10)= -1589941339
      seed(11)= -1136709282
      seed(12)=  860140008
      endif

      call RANDOM_SEED(put=seed)

      call setGRVs(3*Nv,theta)
      theta=theta
      call coef(u,theta)
      baremass=0.01
      MDW=1.0  ! approximation eig(Hs)=eig(Hw)/[2+eig(Hw)] only holds for MDW=1.0

      if (SHAMIRNOISY) then
        call validateOLNoisy(1,6,20,2)
      endif

      if (WILSONNOISY) then
        call validateOLNoisy(2,19,24,1)
      endif


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

      return
      end subroutine validateCondensate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine validateOLAnalyticPoint(kernel)
      use kernelspectrarange
      use statsmod
!     when comparing with notes, rem MDW=0.9 for analytic solution
!     and use 8x8 for thesis cases
      implicit none
      integer,intent(in) :: kernel
      real(prc) lmin,lmax
      integer j,jmin,jmax,shft
      integer idn,is,LsTest
      type(sgnratfunc) :: SRF
c      character(len=80) :: fname,fname2
      real :: diff1,diff3
      complex(prc) pbpHM1S,pbpHM3S,pbpHM1W,pbpHM3W
      complex(prc) pbpZM1S,pbpZM3S,pbpZM1W,pbpZM3W
      complex(prc),dimension(knoise) :: pHM1Ss,pHM3Ss,pHM1Ws,pHM3Ws
      complex(prc),dimension(knoise) :: pZM1Ss,pZM3Ss,pZM1Ws,pZM3Ws
      complex(prc) pbpHM1Sav,pbpHM3Sav,pbpHM1Wav,pbpHM3Wav
      complex(prc) pbpZM1Sav,pbpZM3Sav,pbpZM1Wav,pbpZM3Wav
      complex(prc) pbpHM1Ssd,pbpHM3Ssd,pbpHM1Wsd,pbpHM3Wsd
      complex(prc) pbpZM1Ssd,pbpZM3Ssd,pbpZM1Wsd,pbpZM3Wsd

!     analytic with MDW=0.9
      call makeGaugeField(.true.)
      call coef(u,theta)
      baremass=0.01
      MDW=0.9  ! approximation eig(Hs)=eig(Hw)/[2+eig(Hw)] only holds for MDW=1.0

      dwkernel=kernel ! 1=Shamir,2=Wilson
      call getMinMaxHEigs(lmin,lmax)
      print *,"lmin:",lmin
      print *,"lmax:",lmax
      if (dwkernel.eq.1) then
        fname="EigShamir8x8.dat"
      elseif (dwkernel.eq.2) then
        fname="EigWilson8x8.dat"
      endif
      open(unit=11,file=fname,status='unknown',form='formatted')
      write(11,'(2E13.6)') lmin,lmax
      close(11)

      jmin=3
      jmax=10
      if (dwkernel.eq.1) then
        fname="FreeShamir8x8.dat"
      elseif (dwkernel.eq.2) then
        fname="FreeWilson8x8.dat"
      endif
      open(unit=11,file=fname,status='unknown',form='formatted')
      do j=jmin,jmax
        call setHTcoeffs(j,SRF)
        MTYPE=1
        call evalCondPoint_OL(u,pbpHM1S,1,SRF)
        MTYPE=3
        call evalCondPoint_OL(u,pbpHM3S,1,SRF)

        call setZoloCoeffs(j,SRF,lmin,lmax)
        MTYPE=1
        call evalCondPoint_OL(u,pbpZM1S,1,SRF)
        MTYPE=3
        call evalCondPoint_OL(u,pbpZM3S,1,SRF)
        print *,"pbp M1/M3:",j,pbpHM1S,pbpHM3S
        write(11,'(I6,4E19.10)') j,real(pbpHM1S),real(pbpHM3S),
     &                             real(pbpZM1S),real(pbpZM3S)
      end do
      close(11)

      return
      end subroutine validateOLAnalyticPoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine validateOLNoisy(kernel,jmin,jmax,shft)
      use kernelspectrarange
      use statsmod
      implicit none
      integer,intent(in) :: kernel
      integer,intent(in) :: jmin,jmax,shft
      real(prc) lmin,lmax
      integer j,idn
      type(sgnratfunc) :: SRF
      real div
      real :: diff1,diff3
      complex(prc) pbpHM1,pbpHM3,pbpHM4
      complex(prc) pbpZM1,pbpZM3,pbpZM4
      complex(prc),dimension(knoise) :: pHM1s,pHM3s,pHM4s
      complex(prc),dimension(knoise) :: pZM1s,pZM3s,pZM4s
      complex(prc) pbpHM1av,pbpHM3av,pbpHM4av
      complex(prc) pbpZM1av,pbpZM3av,pbpZM4av
      complex(prc) pbpHM1sd,pbpHM3sd,pbpHM4sd
      complex(prc) pbpZM1sd,pbpZM3sd,pbpZM4sd


      dwkernel=kernel ! 1=Shamir,2=Wilson
      call getMinMaxHEigs(lmin,lmax)
      print *,"lmin:",lmin
      print *,"lmax:",lmax
      if (dwkernel.eq.1) then
        fname="EigShamirNoisy8x8.dat"
      elseif (dwkernel.eq.2) then
        fname="EigWilsonNoisy8x8.dat"
      endif
      open(unit=11,file=fname,status='unknown',form='formatted')
      write(11,'(2E13.6)') lmin,lmax
      close(11)

      if (dwkernel.eq.1) then
        fname="NoisyShamir8x8.dat"
      elseif (dwkernel.eq.2) then
        fname="NoisyWilson8x8.dat"
      endif
      open(unit=11,file=fname,status='unknown',form='formatted')
      do j=jmin,jmax,shft
        if (dwkernel.eq.1) then
          fname2="NoiseS"//trim(itoa(j))//".dat"
        elseif (dwkernel.eq.2) then
          fname2="NoiseW"//trim(itoa(j))//".dat"
        endif
        print *,"fname2: ",fname2
        open(unit=14,file=fname2,status='unknown',form='formatted')
        do idn=1,knoise
          call setHTcoeffs(j,SRF)
          MTYPE=1
          call evalCondNoisy_OL(u,pbpHM1,SRF)
          pHM1s(idn)=pbpHM1
          MTYPE=3
          call evalCondNoisy_OL(u,pbpHM3,SRF)
          pHM3s(idn)=pbpHM3
          MTYPE=4
          call evalCondNoisy_OL(u,pbpHM4,SRF)
          pHM4s(idn)=pbpHM4

          call setZoloCoeffs(j,SRF,lmin,lmax)
          MTYPE=1
          call evalCondNoisy_OL(u,pbpZM1,SRF)
          pZM1s(idn)=pbpZM1
          MTYPE=3
          call evalCondNoisy_OL(u,pbpZM3,SRF)
          pZM3s(idn)=pbpZM3
          MTYPE=4
          call evalCondNoisy_OL(u,pbpZM4,SRF)
          pZM4s(idn)=pbpZM4
          print *,"pbp M1/M3/M4:",j,idn,pbpHM1,pbpHM3,pbpHM4
          print *,"pbp M1/M3/M4:",j,idn,pbpZM1,pbpZM3,pbpZM4
          write(14,'(2I6,6E17.6)') j,idn,real(pbpHM1),real(pbpHM3),
     &                                   real(pbpHM4),real(pbpZM1),
     &                                   real(pbpZM4),real(pbpZM4)
        end do
        close(14)
        call calcVar(knoise,pHM1s,pbpHM1av,pbpHM1sd)
        call calcVar(knoise,pHM3s,pbpHM3av,pbpHM3sd)
        call calcVar(knoise,pHM4s,pbpHM4av,pbpHM4sd)
        call calcVar(knoise,pZM1s,pbpZM1av,pbpZM1sd)
        call calcVar(knoise,pZM3s,pbpZM3av,pbpZM3sd)
        call calcVar(knoise,pZM4s,pbpZM4av,pbpZM4sd)
        div=sqrt(real(knoise))
        write(11,'(I6,12E15.6)') j,real(pbpHM1av),real(pbpHM3av),
     &                            real(pbpHM4av),real(pbpZM1av),
     &                            real(pbpZM3av),real(pbpZM4av),
     &                    real(pbpHM1sd)/div,real(pbpHM3sd)/div,
     &                    real(pbpHM4sd)/div,real(pbpZM1sd)/div,
     &                    real(pbpZM3sd)/div,real(pbpZM4sd)/div
      end do
      close(11)
      
      return
      end subroutine validateOLNoisy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module validatecondmodule
