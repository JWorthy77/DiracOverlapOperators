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
      use hmc2wilsonferms
      implicit none

      logical GZERO
      integer isw,Nsw
      real(prc) dH
      real(prc) thetat(Nv,3)

      if (GZERO) then
        theta=0
        call coef(u,theta)
        return
      end if

      Naccepted=0
      Nsw=10
!     loop over Nsweep Hybrid MC steps
      do isw=1,Nsw
        print *,"sweep:",isw," of ",Nsw
        thetat=theta
        print *,"march:",isw," of ",Nsw
        call march2DW(dH,thetat)
!        call march2DomWallFerms(dH,thetat)
        print *,"accept:",isw," of ",Nsw
        call accept(dH,thetat)
      end do
      call coef(u,theta)
      print *,Naccepted,"accepted of",Nsw
      write(100,*) Naccepted,"accepted of",Nsw
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

      y=exp(dH)
      x=urv()
      print *,"dH:",dH
      print *,"accept if x:",x," < exp(-dH)=y:",y
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
