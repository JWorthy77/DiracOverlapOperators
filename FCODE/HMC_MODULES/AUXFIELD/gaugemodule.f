!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module gaugemodule
      use pacc
      use arraysizes
      use numbers
      use indices
      use rvmodule
      use options
      implicit none
      logical,parameter :: VB_GM=.false.     
      integer Naccepted
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makeGaugeField(GZERO)
      use gaugefield
!      use hmc2wilsonferms
!      use hmc2domwallferms
      use rhmc1domwallferm
      implicit none

      logical GZERO
      integer isw
      real(prc) dH
      real(prc) thetat(Nv,3)

      if (GZERO) then
        theta=0
        call coef(u,theta)
        return
      end if

      theta=0
      call coef(u,theta)
      call initRHMC()

!     loop over Nsweep Hybrid MC steps
      Naccepted=0
      do isw=1,Nswp
        print *,"sweep:",isw," of ",Nswp
        thetat=theta
!        call march2DW(dH,thetat)
!        call march2DomWallFerms(dH,thetat)
        call march1DomWallFerm(dH,thetat)
        call accept(dH,thetat)
      end do
      call coef(u,theta)
      print *,Naccepted,"accepted of",Nswp
      write(100,*) Naccepted,"accepted of",Nswp
      return
      end subroutine makeGaugeField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine accept(dH,thetat) ! MC step: accept new with prob=min(1,exp(H0-H1))
      use gaugefield
      implicit none
      real(prc) dH
      real(prc) thetat(Nv,3)
      real(prc) y
      real x

      y=exp(dH) ! dH = Hprev - Hproposed
      x=urv()
      if (DUPLICATE) then
        read(12,*) x
      endif
      print *,"dH:",dH
      print *,"accept if x:",x," < exp(dH)=y:",y
      if (x.lt.y) then
        print *,"ACCEPT"
        theta=thetat
        Naccepted=Naccepted+1
      else
        print *,"DECLINE"
      end if
      return
      end subroutine accept
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module gaugemodule
