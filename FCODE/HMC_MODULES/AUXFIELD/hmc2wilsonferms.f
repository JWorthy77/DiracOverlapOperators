!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module hmc2wilsonferms
      use pacc
      use arraysizes
      use numbers
      use indices
      use rvmodule
      use options
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine march2DW(dH,thetat)
      ! march dAdt=P       (theta is A)
      !       dPdt=-dSdA
      use gaugefield
      implicit none
      real(prc),intent(out) :: dH
      real(prc),intent(inout) :: thetat(Nv,3)
      real(prc) F(Nv,3)
      real(prc) pp(Nv,3)
      complex(prc) ps(Nv,4),ut(Nv,3)
      real(prc) etime,proby,ytest,avsteps,dt,h0,h1
      integer mu,ts,tsmax

      etime=0.5
      tsmax=10000
      dt=0.05
      avsteps=etime/dt
      proby=one/avsteps

      call coef(ut,thetat)
!     randomise starting momentum
      call setGRVs(3*Nv,pp)
!     randomise pseudo-fermion field
      call setCGRVs(4*Nv,ps)

      h0=ham2DW(thetat,ut,pp,ps)

!     half time-step for pp
      call force(thetat,ut,ps,F)
!      print *,thetat
      call coef(ut,thetat)
      pp=pp-dt*half*F
!     hybrid loop
      proby=1.0/avsteps
      do ts=1,tsmax
        thetat=thetat+dt*pp
!        print *,thetat
        call coef(ut,thetat)
        call force(thetat,ut,ps,F)
        ytest=urv()
        if (ytest.lt.proby) then
          print *,"ts:",ts
          pp=pp-dt*half*F
          goto 501
        else
          pp=pp-dt*F
        endif
      end do


501   continue
      h1=ham2DW(thetat,ut,pp,ps)
      dH=h0-h1
      return
      end subroutine march2DW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine force(thetat,ut,ps,F)
      implicit none
      real(prc),intent(in) :: thetat(Nv,3)
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4)
      real(prc),intent(out) :: F(Nv,3)
      real(prc) :: dSdA(Nv,3)

      call forceThirring3(thetat,F)
      call fermionforce(ut,ps,dSdA) 
      print *,dSdA
!      stop
      F=F+dSdA
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine forceThirring3(theta,dSdA)
      implicit none
      real(prc) theta(Nv,3),dSdA(Nv,3)

      dSdA=Nferms*gbeta*theta

      return
      end subroutine forceThirring3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforce(ut,ps,dSdA) ! for Seff=1/2phi^dag (Dw^dag.Dw)^-1 phi
      use gammas
      use WilsonDirac
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4)
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4) :: eta,nu

      call DdagD(ps,nu,ut,.false.,baremass)
      call DWilson(nu,eta,ut,.false.,baremass)
      call WilsonDerivs(dSdA,eta,nu,.false.)

      return
      end subroutine fermionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2DW(thetat,ut,pp,ps)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat,pp
      complex(prc),dimension(Nv,3),intent(in) :: ut
      complex(prc),dimension(Nv,4),intent(in) :: ps
      real(prc) hg,hp,hf

      hp=0.5*sum(pp*pp)
      hg=hamThirring(thetat)
      hf=ham2WilsonFerms(ps,ut)
      ham2DW=(hg+hp+hf)/Nv
!      ham2DW=(hg+hp)/Nv
      print *,ham2DW,hg/Nv,hp/Nv,hf/Nv
!      print *,ham2DW,hg/Nv,hp/Nv

      return
      end              
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function hamThirring(thetat)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat

      hamThirring=Nferms*gbeta*sum(thetat*thetat)

      return
      end function hamThirring
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2WilsonFerms(ps,ut)
      use WilsonDirac
      implicit none
      complex(prc),intent(in) :: ps(Nv*4)
      complex(prc),intent(in) :: ut(Nv*3)
      complex(prc) tmp(Nv*4)

      print *,'energy 2 Wilson Fermions'
      call IDdagD(ps,tmp,ut,.false.,baremass)
      ham2WilsonFerms=half*dot_product(ps,tmp)

      return
      end function ham2WilsonFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module hmc2wilsonferms
