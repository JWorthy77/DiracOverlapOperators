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
      complex(prc) phi3(Nv,4)
      real(prc) etime,proby,ytest,avsteps,dt,h0,h1
      integer mu,ts,tsmax,l

      etime=0.5
      tsmax=100
      dt=0.05
      avsteps=etime/dt
      proby=one/avsteps

      call coef(ut,thetat)
      call setGRVs(3*Nv,pp)  ! randomise starting momentum
!      call setCGRVs(4*Nv*Ls,ps) ! randomise pseudo-fermion field
      call setCGRVs(4*Nv,phi3) ! randomise pseudo-fermion field
      do l=1,Ls
        ps(:,:,l)=phi3
      end do

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
      subroutine fermionforce(ut,ps,dSdA) ! for Seff=1/2phi ...  phi
      use gammas
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      real(prc) :: dTmpdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls) :: eta,nu,chi
      complex(prc),dimension(Nv,4,Ls) :: TMP
      integer :: KTYPE
      procedure(),pointer :: Dptr=>NULL()
      integer baseMTYPE 
     
      if (VB_DWF) print *,"fermion force A"

      baseMTYPE=MTYPE
      call IMdagMDomWall(ps,chi,ut,baremass)
      if (VB_DWF) print *,"fermion force A - step 1 done"
      call MDomWall(chi,TMP,ut,.false.,baremass)
      if (VB_DWF) print *,"fermion force A - step 2 done"

      MTYPE=1
      call IDDW(TMP,eta,ut,.true.,one)
      if (VB_DWF) print *,"fermion force A - step 3 done"
      MTYPE=baseMTYPE

      call MDomWall(chi,nu,ut,.false.,baremass)
      call DomainWallDerivs(dSdA,eta,nu,.false.,KTYPE)

      eta=chi
      call DomainWallDerivs(dTmpdA,eta,nu,.false.,KTYPE)
      dSdA=dSdA+dTmpdA

      if (VB_DWF) print *,"fermion force A done"

      return
      end subroutine fermionforce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceB(ut,ps,dSdA) ! for Seff=1/2phi^dag (Dw^dag.Dw)^-1 phi
      use gammas
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      real(prc) :: dTmpdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls) :: eta,nu,chi
      complex(prc),dimension(Nv,4,Ls) :: TMP
      integer :: KTYPE
      procedure(),pointer :: Dptr=>NULL()
      integer baseMTYPE 
     
      if (VB_DWF) print *,"fermion force A"

      baseMTYPE=MTYPE
      call IMdagMDomWall(ps,chi,ut,baremass)
      if (VB_DWF) print *,"fermion force A - step 1 done"
      call MDomWall(chi,TMP,ut,.false.,baremass)
      if (VB_DWF) print *,"fermion force A - step 2 done"

      MTYPE=1
      call IDDW(TMP,eta,ut,.true.,one)
      if (VB_DWF) print *,"fermion force A - step 3 done"
      MTYPE=baseMTYPE

      call MDomWall(chi,nu,ut,.false.,baremass)
      call DomainWallDerivs(dSdA,eta,nu,.false.,KTYPE)

      eta=chi
      call DomainWallDerivs(dTmpdA,eta,nu,.false.,KTYPE)
      dSdA=dSdA+dTmpdA

      if (VB_DWF) print *,"fermion force A done"

      return
      end subroutine fermionforceB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine fermionforceC(ut,ps,dSdA) ! for Seff=1/2phi^dag (Dw^dag.Dw)^-1 phi
      use gammas                           ! chi^dag dM^dagdA M chi
      use WilsonDirac
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: ut(Nv,3)
      complex(prc),intent(in) :: ps(Nv,4,Ls)
      real(prc),intent(out) :: dSdA(Nv,3)
      real(prc) :: dTmpdA(Nv,3)
      complex(prc),dimension(Nv,4,Ls) :: eta,nu,chi
      complex(prc),dimension(Nv,4,Ls) :: TMP
      integer,parameter :: KTYPE=1 ! Wilson
      procedure(),pointer :: Dptr=>NULL()
      
      call IMDomWall(ps,chi,ut,.false.,baremass)
      call DDW(chi,TMP,ut,.false.,baremass)
      Dptr => DDW_Wilson
!      Dptr => DDW_OWilson
      call IMdagM_DWkernel(TMP,nu,ut,baremass,Dptr)
      eta=chi
      call DomainWallDerivs(dSdA,eta,nu,.false.,KTYPE)
      call MDomWall(chi,eta,ut,.false.,baremass)
      call DomainWallDerivs(dTmpdA,eta,nu,.false.,KTYPE)
      dSdA=dSdA+dTmpdA

      return
      end subroutine fermionforceC
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
      Dptr => MDomWall
      if (VB_DWF) print *,'energy 2 DomWall Fermions'
      call IMdagM_DWkernel(ps,tmp,ut,baremass,Dptr)
      ham2Ferms=half*dot_product(ps,tmp)

      return
      end function ham2Ferms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module hmc2domwallferms
