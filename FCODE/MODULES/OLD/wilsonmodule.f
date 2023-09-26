      module wilsonmodule
      use pacc
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DWilson(R,DR,u,DAGGER,mass,add)
      use gammas
      use indices
      use options
      implicit none
c     calculates DR = DWilson*R
      complex(prc),intent(in) :: R(Nv,Ndc)
      complex(prc),intent(out) :: DR(Nv,Ndc)
      complex(prc),intent(in) :: u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      real(prc) mult
      complex(prc) cmass
      integer i,idirac
      integer mu,igork

      mult=one
      if (DAGGER) then
        mult=-one
      end if

      cmass=mass+add
      if (DAGGER) then
        cmass=conjg(mass+add)
      end if
!     mass term 
      do idirac=1,Ndc
        do i=1,Nv
          DR(i,idirac) = cmass*R(i,idirac)
        enddo
      enddo

!     Dirac term (anti-Hermitian)
      do mu=1,NDT
        do idirac=1,Ndc
          igork=gamin(mu,idirac)
          do i=1,Nv
            DR(i,idirac)=DR(i,idirac)
     &        +mult*gamval(mu,idirac)*
     &         (u(i,mu)*R(iu(i,mu),igork)
     &         -conjg(u(id(i,mu),mu))*R(id(i,mu),igork))/two
          enddo
        enddo
      enddo

!     Wilson term (Hermitian)
      do mu=1,NDT
        do idirac=1,Ndc
          do i=1,Nv
            DR(i,idirac)=DR(i,idirac)
     &        -( u(i,mu)*R(iu(i,mu),idirac) - two*R(i,idirac)
     &         +conjg(u(id(i,mu),mu))*R(id(i,mu),idirac))/two
          enddo
        enddo
      enddo

      return
      end subroutine DWilson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DpDdag(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = (DWilson+DWilson')/2*R
      complex(prc),intent(in) :: R(Nv,Ndc)
      complex(prc),intent(out) :: DR(Nv,Ndc)
      complex(prc),intent(in) :: u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) :: TMP(Nv,4)

      call DWilson(R,DR,u,DAGGER,mass,add)
      call DWilson(R,TMP,u,.not.DAGGER,mass,add)
      DR=(DR+TMP)/2

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine G5DW(R,DR,u,DAGGER,mass,add)
      use gammas
      implicit none
c     calculates DR = G5.DWilson*R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) :: TMP(Nv,4)

      if(.not.DAGGER) then
        call DWilson(R,DR,u,DAGGER,mass,add)
        call mGmu(DR,5)
      else
        TMP=R
        call mGmu(TMP,5)
        call DWilson(TMP,DR,u,DAGGER,mass,add)
      end if
      
      return
      end subroutine G5DW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DW2(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = DWilson*R
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,4)

      call DWilson(R,TMP,u,DAGGER,mass,add)
      call DWilson(TMP,DR,u,DAGGER,mass,add)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DdagPD(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = (DWilson'+DWilson)/2*R
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,4)

      call DWilson(R,DR,u,.false.,mass,add)
      call DWilson(R,TMP,u,.true.,mass,add)
      DR=(DR+TMP)/2

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DW2pC(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = DWilson*R
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

      call DWilson(R,TMP,u,DAGGER,mass,czero)
      call DWilson(TMP,DR,u,DAGGER,mass,czero)
      DR=DR+add*R

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DdagD(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = DWilson*R
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

      call DWilson(R,TMP,u,DAGGER,mass,add)
      call DWilson(TMP,DR,u,.not.DAGGER,mass,add)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DdagDpC(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = DWilson*R
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

c      call DWilson(R,TMP,u,.false.,mass,czero)
c      call DWilson(TMP,DR,u,.true.,mass,czero)
      call DWilson(R,TMP,u,DAGGER,mass,czero)
      call DWilson(TMP,DR,u,.NOT.DAGGER,mass,czero)
      DR=DR+add*R

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDW(RR,DR,u,DAGGER,mass,add)
!     solve Dw(u) DR = RR or DwD(u) DR = RR
      use axbmodule1
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,4)
      procedure(),pointer :: Mptr

      Mptr=>DWilson
      if (.not. DAGGER) then
        call Mptr(RR,TMP,u,.true.,mass,add)
        call IMdagM(TMP,DR,u,.false.,mass,add,Mptr)
      elseif (DAGGER) then
        call IMdagM(RR,TMP,u,.true.,mass,add,Mptr)
        call Mptr(TMP,DR,u,.false.,mass,add)
      end if

      return
      end subroutine IDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDdagD(R,DR,u,DAGGER,mass,add)
!     solve DdagD DR = R or DDdag(u) DR = R
      use axbmodule1
      implicit none
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr

      Mptr=>DWilson
      call IMdagM(R,DR,u,DAGGER,mass,add,Mptr)

      return
      end subroutine IDdagD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IG5DW(RR,DR,u,DAGGER,mass,add)
!     solve Dw(u) DR = RR or DwD(u) DR = RR
      use axbmodule1
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,4)
      procedure(),pointer :: Mptr

      Mptr=>G5DW
      if (.not. DAGGER) then
        call Mptr(RR,TMP,u,.true.,mass,add)
        call IMdagM(TMP,DR,u,.false.,mass,add,Mptr)
      elseif (DAGGER) then
        call IMdagM(RR,TMP,u,.true.,mass,add,Mptr)
        call Mptr(TMP,DR,u,.false.,mass,add)
      end if

      return
      end subroutine IG5DW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDWOperators
      use numbers
      use gammas
      use indices
      use rvmodule
      use gaugefield
      use axbmodule1
      implicit none
      complex(prc) psi(Nv,4),Dpsi(Nv,4),DdDpsi(Nv,4)
      complex(prc) IDpsi(Nv,4),test(Nv,4)
      complex(prc) TMP1(Nv,4),TMP2(Nv,4),TMP3(Nv,4)
      complex(prc) urvs(Nv)
      real(prc) alpha(Nv)
      integer idirac,iv,idx
      complex(prc) add
      real(prc) mass
      procedure(),pointer :: Mptr => NULL()
      
      call setGammas
      call setIndices
      theta=zero
      call coef(u,theta)
      mass=zero/100
      add=one

      print *,'test Wilson Dirac Operators'

      do idirac=1,4
        do iv=1,Nv
          psi(iv,idirac) = iv*cmplx(one,zero) + idirac*cmplx(zero,one) 
        end do
      end do

!     check Mptr=>DWilson matches
      call DWilson(psi,Dpsi,u,.false.,mass,add)
      Mptr => DWilson
      call Mptr(psi,DdDpsi,u,.false.,mass,add)
      print *,maxval(abs(Dpsi-DdDpsi))

!     check Mptr => DdagD matches
      Mptr => DdagD
      call DWilson(psi,TMP1,u,.false.,mass,add)
      call DWilson(TMP1,TMP2,u,.true.,mass,add)
      call Mptr(psi,TMP1,u,.false.,mass,add)
      print *,maxval(abs(TMP1-TMP2))

      Mptr => DWilson
      call Mptr(psi,TMP1,u,.false.,mass,add)
      call IMnonsym(TMP1,TMP2,u,.false.,mass,add,Mptr)
      print *,maxval(realpart(TMP2-psi)),maxval(imagpart(TMP2-psi))
      Mptr => DWilson
      call Mptr(psi,TMP1,u,.true.,mass,add)
      call IMnonsym(TMP1,TMP2,u,.true.,mass,add,Mptr)
      print *,maxval(realpart(TMP2-psi)),maxval(imagpart(TMP2-psi))
     
      Mptr => DdagD
      call Mptr(psi,TMP1,u,.false.,mass,add)
      call IM(TMP1,TMP2,u,.false.,mass,add,Mptr)
      print *,maxval(realpart(TMP2-psi)),maxval(imagpart(TMP2-psi))

      Mptr => DdagDpC
      call Mptr(psi,TMP1,u,.false.,mass,add)
      call IM(TMP1,TMP2,u,.false.,mass,add,Mptr)
      print *,maxval(realpart(TMP2-psi)),maxval(imagpart(TMP2-psi))

      call G5DW(psi,TMP1,u,.false.,mass,add)
      call IG5DW(TMP1,TMP2,u,.false.,mass,add)
      print *,maxval(abs(realpart(TMP2-psi))),
     &        maxval(abs(imagpart(TMP2-psi)))

      Mptr => DWilson
      call setRVs(Nv*4,psi)
      call Mptr(psi,TMP1,u,.false.,mass,add)
      call IM_BICGSTAB(TMP1,TMP2,u,.false.,mass,add,Mptr)
      print *,'BICGSTAB error:',maxval(abs(TMP2-psi))

      print *,'Is DdagD=DDag?'
      call setRVs(Nv*4,psi)
      Mptr => DdagD
      call Mptr(psi,TMP1,u,.false.,mass,add)
      call Mptr(psi,TMP2,u,.true.,mass,add)
      print *,maxval(abs(TMP2-TMP1))

      return
      end subroutine testDWOperators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDWGaugeSymmetry
      use numbers
      use gammas
      use indices
      use rvmodule
      use gaugefield
      use options
      implicit none
      complex(prc) psi(Nv,4),Dpsi(Nv,4)
      complex(prc) psit(Nv,4),Dpsit(Nv,4)
      complex(prc) urvs(Nv)
      real(prc) alpha(Nv)
      complex(prc) pbDp,pbDpt
      integer idx
      procedure(),pointer :: Mptr => NULL()

      Mptr=>DWilson
!     set gauge transform field alpha
      call setRVs(Nv,urvs)
      alpha=real(urvs)
!     set phi
      call setRVs(Nv*4,psi)
!     make phibar D phi
      call Mptr(psi,Dpsi,u,.false.,baremass,czero)
      pbDp=czero
      do idx=1,4
        pbDp=pbDp+dot_product(psi(:,idx),Dpsi(:,idx))
      end do
!     transform u
      call gaugeTransformU(u,alpha)
!     transform phi->phit
      do idx=1,4
        psit(:,idx)=psi(:,idx)*exp(cmplx(zero,alpha))
      end do
!     make phitbar D[ut] phit
      call Mptr(psit,Dpsit,u,.false.,baremass,czero)
      pbDpt=czero
      do idx=1,4
        pbDpt=pbDpt+dot_product(psit(:,idx),Dpsit(:,idx))
      end do
      print *,pbDp,pbDpt,abs(pbDp-pbDpt)/abs(pbDp)

!     multiply Dpsit by exp(-i.alpha)
      do idx=1,4
        Dpsit(:,idx)=Dpsit(:,idx)*exp(-cmplx(zero,alpha))
      end do
      print *,maxval(abs(Dpsi-Dpsit))

      return
      end subroutine testDWGaugeSymmetry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Jacobi(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = DWilson^-1*R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) :: TMP(Nv,4)
      integer i,NJ

      NJ=200
      DR=R
      do i=1,NJ
        TMP=DR
        call JacR(TMP,DR,u,DAGGER,mass,add)
        TMP=R-DR
        call JacD(TMP,DR,u,DAGGER,mass,add)
      end do

      return
      end subroutine Jacobi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine JacR(R,DR,u,DAGGER,mass,add)
      use gammas
      use indices
      implicit none
c     calculates DR = DWilson*R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      real(prc) mult
      complex(prc) cmass
      integer i,idirac
      integer mu,igork

      mult=one
      if (DAGGER) then
        mult=-one
      end if

      cmass=mass+add
      if (DAGGER) then
        cmass=conjg(mass+add)
      end if
!     mass term 
c      do idirac=1,4
c        do i=1,Nv
c          DR(i,idirac) = cmass*R(i,idirac)
c        enddo
c      enddo

!     Dirac term (anti-Hermitian)
      do mu=1,3
        do idirac=1,4
          igork=gamin(mu,idirac)
          do i=1,Nv
            DR(i,idirac)=DR(i,idirac)
     &        +mult*gamval(mu,idirac)*
     &         (u(i,mu)*R(iu(i,mu),igork)
     &         -conjg(u(id(i,mu),mu))*R(id(i,mu),igork))/two
          enddo
        enddo
      enddo

!     Wilson term (Hermitian)
      do mu=1,3
        do idirac=1,4
          do i=1,Nv
            DR(i,idirac)=DR(i,idirac) -( 
     &        u(i,mu)*R(iu(i,mu),idirac) 
c     &        - two*R(i,idirac)
     &        +conjg(u(id(i,mu),mu))*R(id(i,mu),idirac)
     &          )/two
          enddo
        enddo
      enddo

      return
      end subroutine JacR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine JacD(R,DR,u,DAGGER,mass,add)
      use gammas
      use indices
      implicit none
c     calculates DR = DWilson*R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      real(prc) mult
      complex(prc) cmass
      integer i,idirac
      integer mu,igork

      mult=one
      if (DAGGER) then
        mult=-one
      end if

      cmass=mass+add
      if (DAGGER) then
        cmass=conjg(mass+add)
      end if
!     mass term 
      do idirac=1,4
        do i=1,Nv
          DR(i,idirac) = one/(cmass+two)*R(i,idirac)
        enddo
      enddo

!     Dirac term (anti-Hermitian)
c      do mu=1,3
c        do idirac=1,4
c          igork=gamin(mu,idirac)
c          do i=1,Nv
c            DR(i,idirac)=DR(i,idirac)
c     &        +mult*gamval(mu,idirac)*
c     &         (u(i,mu)*R(iu(i,mu),igork)
c     &         -conjg(u(id(i,mu),mu))*R(id(i,mu),igork))/two
c          enddo
c        enddo
c      enddo

!     Wilson term (Hermitian)
c      do mu=1,3
c        do idirac=1,4
c          do i=1,Nv
c            DR(i,idirac)=DR(i,idirac) -( 
c     &        u(i,mu)*R(iu(i,mu),idirac) 
c     &        - two*R(i,idirac)
c     &        +conjg(u(id(i,mu),mu))*R(id(i,mu),idirac)
c     &          )/two
c          enddo
c        enddo
c      enddo

      return
      end subroutine JacD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module wilsonmodule
