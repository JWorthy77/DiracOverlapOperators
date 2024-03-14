!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module gaugemodule
      use pacc
      use arraysizes
      use numbers
      use indices
      use rvmodule
      use options
      implicit none

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

      Nsw=50
!     loop over Nsweep Hybrid MC steps
      do isw=1,Nsw
        print *,"sweep:",isw," of ",Nsw
        thetat=theta
!        call march(dH,thetat)
        call march2DW(dH,thetat)
        call accept(dH,thetat)
      end do
      call coef(u,theta)

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
      print *,"x:",x,"y:",y,"dH:",dH
      if (x.lt.y) then
        print *,"ACCEPT"
        theta=thetat
      else
        print *,"DECLINE"
      end if
      return
      end subroutine accept
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module gaugemodule
