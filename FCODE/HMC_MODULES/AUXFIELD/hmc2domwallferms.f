!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module hmc2domwallferms
      use pacc
      use arraysizes
      use numbers
      use indices
      use rvmodule
      use options
      implicit none
      logical,parameter :: VB_DWF=.true.
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine march2DomWallFerms(dH,thetat) ! march dAdt=P (theta is A)
      use gaugefield                           !       dPdt=-dSdA
      use hmc2wilsonferms
      implicit none
      real(prc),intent(out) :: dH
      real(prc),intent(inout) :: thetat(Nv,3)
      real(prc) F(Nv,3)
      real(prc) pp(Nv,3)
      complex(prc) ps(Nv,4,Ls),ut(Nv,3,Ls)
      real(prc) etime,proby,ytest,avsteps,dt,h0,h1
      integer mu,ts,tsmax

      etime=1.0
      tsmax=100
      dt=0.1
      avsteps=etime/dt
      proby=one/avsteps

      call coef(ut,thetat)
      call setGRVs(3*Nv,pp)  ! randomise starting momentum
      call setCGRVs(4*Nv*Ls,ps) ! randomise pseudo-fermion field

      h0=ham2DomWallFerms(thetat,ut,pp,ps) ! initial hamiltonian energy
      if (VB_DWF) then ; print *,"h0:",h0 ; end if
      call force2DomWallFerms(thetat,ut,ps,F)
      pp=pp-dt*half*F ! half time step before leap frog
      proby=1.0/avsteps
      do ts=1,tsmax  ! time march
        thetat=thetat+dt*pp
!        call coef(ut,thetat)
        call force2DomWallFerms(thetat,ut,ps,F)
        ytest=urv()
        if (ytest.lt.proby) then
          print *,"ts:",ts
          pp=pp-dt*half*F
          goto 501
        else
          pp=pp-dt*F
        endif
      end do
      pp=pp+half*dt*F ! correction to half step if tsmax reached

501   continue
      h1=ham2DomWallFerms(thetat,ut,pp,ps) ! final hamiltonian energy
      if (VB_DWF) then ; print *,"h1:",h1 ; end if
      dH=h0-h1
      return
      end subroutine march2DomWallFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine force2DomWallFerms(thetat,ut,ps,F)
      implicit none
      real(prc),intent(in) :: thetat(Nv,3)
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: F(Nv,3)
      real(prc) :: dSdA(Nv,3)

      call forceThirring3(thetat,F)
      call fermionforce(ut,ps,dSdA) 
      F=F+half*dSdA
      return
      end subroutine force2DomWallFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine forceThirring3(thetat,dSdA)
      implicit none
      real(prc) thetat(Nv,3),dSdA(Nv,3)

      dSdA=Nferms*gbeta*thetat

      return
      end subroutine forceThirring3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforce(ut,ps,dSdA) ! for Seff=1/2phi^dag (Dw^dag.Dw)^-1 phi
      use gammas
      use WilsonDirac
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      complex(prc),dimension(Nv,4) :: eta,nu
      
      call IDdagD(ps,nu,ut,.false.,baremass)
      call DWilson(nu,eta,ut,.false.,baremass)
      call WilsonDerivs(dSdA,eta,nu,.false.)

      return
      end subroutine fermionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2DomWallFerms(thetat,ut,pp,ps)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat,pp
      complex(prc),dimension(Nv,3),intent(in) :: ut
      complex(prc),dimension(Nv,4,Ls),intent(in) :: ps
      real(prc) hg,hp,hf

      hp=0.5*sum(pp*pp)
      hg=hamThirring(thetat)
      hf=ham2Ferms(ps,ut)
      ham2DomWallFerms=(hg+hp+hf)/Nv
      if (VB_DWF) print *,"hg:",hg/Nv
      if (VB_DWF) print *,"hp:",hp/Nv
      if (VB_DWF) print *,"hf:",hf/Nv
      if (VB_DWF) print *,"h:",ham2DomWallFerms
      write(101,*) hg/Nv,hp/Nv,hf/Nv
      return
      end function ham2DomWallFerms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function hamThirring(thetat)
      implicit none
      real(prc),dimension(Nv,3),intent(in) :: thetat

      hamThirring=Nferms*half*gbeta*sum(thetat*thetat)

      return
      end function hamThirring
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(prc) function ham2Ferms(ps,ut)
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ps(Nv*4*Ls)
      complex(prc),intent(in) :: ut(Nv*3)
      complex(prc) tmp(Nv*4*Ls)
      procedure(),pointer :: Dptr=>NULL()

      Dptr => DDW_OWilson
      Dptr => DDW_Wilson
      if (VB_DWF) print *,'energy 2 DomWall Fermions'
      call IMdagM_DWkernel(ps,tmp,ut,.false.,baremass,Dptr)
      ham2Ferms=half*dot_product(ps,tmp)

      return
      end function ham2Ferms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module hmc2domwallferms
