!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module scalargaugemodule
#ifdef STAGGERED
      use pacc
      use arraysizes
      use numbers
      use indices
      use rvmodule
      use options
      implicit none

      integer Naccept

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makeScalarBosonField(GZERO,sigma)
      implicit none

      logical GZERO
      real(prc),intent(out) :: sigma(Nv)
      integer isw,Nsw
      real(prc) dH
      real(prc) sigmat(Nv)

      if (GZERO) then
        sigma=one
        return
      end if
      Nsw=10
      Naccept=0
!     loop over Nsweep Hybrid MC steps
      do isw=1,Nsw
        print *,"sweep:",isw," of ",Nsw
        sigmat=sigma
        call marchSigma(dH,sigmat)
        call acceptSigma(dH,sigmat,sigma)
      end do
      print *,"Acceptance rate:",Naccept," of ",Nsw

      return
      end subroutine makeScalarBosonField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine marchSigma(dH,sigmat)
      use staggeredmodule
      implicit none
      real(prc) dH,H0,H1
      real(prc) sigmat(Nv)
      real(prc) dSdsi(Nv)
      real(prc) dt,hg,hp,hf
      real(prc) pp(Nv),phi(Nv,Nferms/2)
      real(prc) etime,proby,ytest,avsteps
      integer mu,ts,tsmax,i

      etime=0.5
      tsmax=10000
      dt=0.05
      avsteps=etime/dt
      proby=one/avsteps

!     set noise for dynamic part
      if(DYNAMIC.gt.0) then
        call setPhi(sigmat,phi)
        phi=phi/Nv
      end if
c      do i=1,Nferms/2
c        print *,sum(phi(:,i))/Nv
c      end do

!     set heatbath (initialise) for p 
      call setGRVs(Nv,pp)

!     calc initial hamiltonian
      call hamiltonGN(sigmat,pp,phi,H0,hg,hp,hf)

!     half time-step for pp
      call forceGN(sigmat,phi,dSdsi)
      pp=pp-dt/2*dSdsi

!     hybrid loop
      proby=1.0/avsteps
      do ts=1,tsmax
        sigmat=sigmat+dt*pp
        call forceGN(sigmat,phi,dSdsi)
        ytest=urv()
        if (ytest.lt.proby) then
          pp=pp-dt/2*dSdsi
          goto 501
        else
          pp=pp-dt*dSdsi
        endif
        call hamiltonGN(sigmat,pp,phi,H1,hg,hp,hf)
      end do

501   continue
      call hamiltonGN(sigmat,pp,phi,H1,hg,hp,hf)
      print *,hg/Nv,hp/Nv
      dH=H0-H1
      return
      end subroutine marchSigma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine hamiltonGN(sigma,pp,phi,h,hg,hp,hf)
      implicit none
      real(prc) sigma(Nv),pp(Nv),phi(Nv,Nferms/2)
      real(prc) h,hg,hp,hf

      hp=0.5*sum(pp*pp)
      call SG_GN(sigma,hg)
      h=hg+hp

      if (DYNAMIC.gt.0) then
        call SF_GN(sigma,phi,hf)
        h=h+hf
      endif

      return
      end subroutine hamiltonGN            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SG_GN(sigma,hg)
      implicit none
      real(prc) sigma(Nv)
      real(prc) hg

      hg=Nferms*gbeta*sum(sigma*sigma)/4

      return
      end subroutine SG_GN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setPhi(sigma,phi)
!     phi=M^T.M^{-1}.grv
      use staggeredmodule
      use axbmoduleReal
      implicit none
      real(prc) sigma(Nv),phi(Nv,Nferms/2),tmp(Nv)
      procedure(),pointer :: Mptr => NULL()
      integer i

!!#ifdef STAGGERED
      Mptr => DStagGN
!!#else
!!      print *,"STAGGERED not switched on"
!!      stop
!!#endif
      do i=1,Nferms/2
        call setGRVs(Nv,tmp)
        call IMdagMReal(tmp,phi(:,i),sigma,.false.,Mptr)
      end do

      return
      end subroutine setPhi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SF_GN(sigma,phi,hf)
      use staggeredmodule
      use axbmoduleReal
      implicit none
      real(prc) sigma(Nv),phi(Nv,Nferms/2)
      real(prc) hf
      real(prc) DR(Nv)
      procedure(),pointer :: Mptr => NULL()
      integer i,j

!!#ifdef STAGGERED
      Mptr => DStagGN
!!#else
!!      print *,"STAGGERED not switched on"
!!      stop
!!#endif
      hf=0
      do i=1,Nferms/2
c        call IMdagMReal(phi(:,i),DR,sigma,.false.,Mptr)
        do j=1,Nferms/2
c          hf=hf+sum(phi(:,j)*DR)/2
          hf=hf+sum(phi(:,j)*phi(:,i))/2
        end do
      end do

      return
      end subroutine SF_GN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine forceGN(sigma,phi,dSdsi)
      implicit none
      real(prc),intent(in) :: sigma(Nv),phi(Nv,Nferms/2)
      real(prc),intent(out) :: dSdsi(Nv)

      call dSGdsi_GN(sigma,dSdsi)
      if (DYNAMIC.gt.0) then
        call dSFdsi_GN(sigma,phi,dSdsi) 
      endif

      return
      end subroutine forceGN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dSGdsi_GN(sigma,dSdsi)
      implicit none
      real(prc) sigma(Nv),dSdsi(Nv)
      integer i

      do i=1,Nv
        dSdsi(i)=Nferms*gbeta*sigma(i)/two
      end do

      return
      end subroutine dSGdsi_GN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dSFdsi_GN(sigma,phi,dSdsi)
      use staggeredmodule
      use axbmoduleReal
      implicit none
      real(prc),intent(in) :: sigma(Nv),phi(Nv,Nferms/2)
      real(prc),intent(out) :: dSdsi(Nv)
      real(prc),dimension(Nv,Nferms/2) :: RHS
      real(prc),dimension(Nv) :: P2A,P2B,P1A,P1B,COMB
      procedure(),pointer :: Mptr => NULL()
      integer i,j,k,Nf2

c      Mptr => DStagGN
c      Nf2=Nferms/2
c      do i=1,Nf2
c        call IMdagMReal(phi(:,i),RHS(:,i),sigma,.false.,Mptr)
c      end do

      do k=1,Nv
        do i=1,Nf2
          call DStagGNderiv(phi(:,i),P2A,sigma,.false.,k)
          call DStagGN(P2A,P2B,sigma,.true.)
          call DStagGN(phi(:,i),P1A,sigma,.false.)
          call DStagGNderiv(P1A,P1B,sigma,.true.,k)
          COMB=P1B+P2B
          do j=1,Nf2
            dSdsi(k)=dSdsi(k)+sum(phi(:,j)*COMB)/two
          end do
        end do
      end do

      return
      end subroutine dSFdsi_GN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine acceptSigma(dH,sigmat,sigma)
      use gaugefield
      implicit none
      real(prc) dH
      real(prc) sigmat(Nv),sigma(Nv)
      real(prc) y
      real x
c**********************************************************************
c  Monte Carlo step: accept new fields with probability=
c              min(1,exp(H0-H1))
c**********************************************************************
      y=exp(dH)
      x=urv()
      print *,"Acceptance step:",x,y
      if (x.lt.y) then
        sigma=sigmat
        Naccept=Naccept+1
      end if
      return
      end subroutine acceptSigma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif
      end module scalargaugemodule
