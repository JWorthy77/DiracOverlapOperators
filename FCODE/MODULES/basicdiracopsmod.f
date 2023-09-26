      module basicdiracopsmod
      use pacc
      use arraysizes
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DNaive(R,DR,u,DAGGER,mass,add)
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

      DR=czero

      mult=one
      if (DAGGER) then
        mult=-one
      end if

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

      return
      end subroutine DNaive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine WilsonTerm(R,DR,u,DAGGER,mass,add)
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

      DR=czero

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
      end subroutine WilsonTerm
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
      subroutine DShamir(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = DShamir*R = (2+Dw)^-1*Dw*R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) :: TMP(Nv,4)

      if (.not.DAGGER) then
        call DWilson(R,TMP,u,DAGGER,mass,czero)
        call IDW(TMP,DR,u,DAGGER,mass,ctwo)
        DR=DR+add*R
      else
        call DWilson(R,TMP,u,DAGGER,mass,czero)
        call IDW(TMP,DR,u,DAGGER,mass,ctwo)
        DR=DR+add*R
      endif


      return
      end subroutine DShamir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DShamirPoly(R,DR,u,DAGGER,mass,add)
      use polycoeffsmod
      implicit none
c     calculates DR = DShamir*R = (2+Dw)^-1*Dw*R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical,intent(in) :: DAGGER
      real(prc),intent(in) :: mass
      complex(prc),intent(in) :: add
      complex(prc) :: TMP(Nv,4),TMPR(Nv,4)
      integer n,j

      associate(poly => pfcs%coeffs)
      n=size(poly)

      call DWilson(R,TMPR,u,DAGGER,mass,add)
      DR=poly(1)*R+poly(2)*TMPR
      do j=3,n
        TMP=TMPR
        call DWilson(TMP,TMPR,u,DAGGER,mass,add)
        DR=DR+poly(j)*TMPR
      enddo
      TMP=DR
      call DWilson(TMP,DR,u,DAGGER,mass,add)

      end associate
      return
      end subroutine DShamirPoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DHermWilson(R,DR,u,DAGGER,mass,add)
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
      subroutine Dbasic(R,DR,u,DAGGER,mass,add)
      use options
      implicit none
c     chooses between Shamir and Wilson operator
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add

      if (dwkernel.eq.1) then
        call DShamir(R,DR,u,DAGGER,mass,add)
      elseif (dwkernel.eq.2) then
        call DWilson(R,DR,u,DAGGER,mass,add)
      end if
      return
      end subroutine Dbasic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDbasic(R,DR,u,DAGGER,mass,add)
      use options
      implicit none
c     chooses between Shamir and Wilson operator
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add

      if (dwkernel.eq.1) then
        call IDS(R,DR,u,DAGGER,mass,add)
      elseif (dwkernel.eq.2) then
        call IDW(R,DR,u,DAGGER,mass,add)
      end if
      return
      end subroutine IDbasic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Hkernel(R,DR,u,DAGGER,mass,add)
      use gammas
      implicit none
c     calculates Hw or Hs
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) :: TMP(Nv,4)

      if(.not.DAGGER) then
        call Dbasic(R,DR,u,DAGGER,mass,add)
        call mGmu(DR,4)
      else
        TMP=R
        call mGmu(TMP,4)
        call Dbasic(TMP,DR,u,DAGGER,mass,add)
      endif
      return
      end subroutine Hkernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IHkernel(R,DR,u,DAGGER,mass,add)
      use gammas
      implicit none
c     calculates Hw or Hs
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) :: TMP(Nv,4)

      if(.not.DAGGER) then
        TMP=R
        call mGmu(TMP,4)
        call IDbasic(TMP,DR,u,DAGGER,mass,add)
      else
        call IDbasic(R,DR,u,DAGGER,mass,add)
        call mGmu(DR,4)
      endif
      return
      end subroutine IHkernel
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
      subroutine DdagD(R,DR,u,DAGGER,mass,add)
      implicit none
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
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

      call DWilson(R,TMP,u,DAGGER,mass,czero)
      call DWilson(TMP,DR,u,.NOT.DAGGER,mass,czero)
      DR=DR+add*R

      return
      end subroutine DdagDpC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine scaledDdagDpC(R,DR,u,DAGGER,mass,add,alpha) 
      implicit none
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      real(prc) alpha
      complex(prc) TMP(Nv,Ndc)

      call DWilson(R,TMP,u,DAGGER,mass,czero)
      call DWilson(TMP,DR,u,.NOT.DAGGER,mass,czero)
      DR=alpha*alpha*DR+add*R

      return
      end subroutine scaledDdagDpC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SdagSpC(R,DR,u,DAGGER,mass,add) 
      implicit none
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

c      call DShamir(R,TMP,u,DAGGER,mass,czero)
c      call DShamir(TMP,DR,u,.NOT.DAGGER,mass,czero)
      call DShamir(R,TMP,u,DAGGER,mass,czero)
      call DShamir(TMP,DR,u,.NOT.DAGGER,mass,czero)
      DR=DR+add*R

      return
      end subroutine SdagSpC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine SdagSpCpoly(R,DR,u,DAGGER,mass,add) 
      implicit none
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

      call DShamirPoly(R,TMP,u,DAGGER,mass,czero)
      call DShamirPoly(TMP,DR,u,.NOT.DAGGER,mass,czero)
      DR=DR+add*R

      return
      end subroutine SdagSpCpoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine H2pC(R,DR,u,DAGGER,mass,add) 
      implicit none
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add 
      complex(prc) TMP(Nv,Ndc)
!     DAGGER options don't matter since hermitian
      call DHermWilson(R,TMP,u,DAGGER,mass,czero)
      call DHermWilson(TMP,DR,u,.NOT.DAGGER,mass,czero)
      DR=DR+add*R

      return
      end subroutine H2pC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDW(RR,DR,u,DAGGER,mass,add)
!     solve Dw.DR = RR or Dw(dag).DR = RR
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
        call Mptr(RR,TMP,u,.false.,mass,add)
        call IMdagM(TMP,DR,u,.true.,mass,add,Mptr)
c        call Mptr(TMP,DR,u,.false.,mass,add)
      end if

      return
      end subroutine IDW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDS(RR,DR,u,DAGGER,mass,add)
!     solve DS.DR = RR or DS(dag).DR = RR
!     but add must be czero
      use axbmodule1
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,4)
      procedure(),pointer :: Mptr

      call DWilson(RR,TMP,u,DAGGER,mass,ctwo)
      call IDW(TMP,DR,u,DAGGER,mass,czero)

      return
      end subroutine IDS
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
      subroutine ISdagS(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = DWilson*R
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

      call DWilson(R,DR,u,.not.DAGGER,mass,ctwo)
      call IDdagD(DR,TMP,u,DAGGER,mass,czero)
      call DWilson(TMP,DR,u,DAGGER,mass,ctwo)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine ISdagSpC(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = (DSdag(mass)*DS(mass)+add)^(-1)R
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

      call DWilson(R,DR,u,.not.DAGGER,mass,ctwo)
      call IDWB(DR,TMP,u,DAGGER,mass,add)
      call DWilson(TMP,DR,u,DAGGER,mass,ctwo)

      return
      end subroutine ISdagSpC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDWB(R,DR,u,DAGGER,mass,add)
      use axbmodule1
      implicit none
c     calculates DR = (DSdag(mass)*DS(mass)+add)^(-1)R
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)
      procedure(),pointer :: Mptr

      Mptr=>DWB
      call IM2(R,DR,u,DAGGER,mass,add,Mptr)

      return
      end subroutine IDWB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DWB(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = [(1+b)DdagD+Ddag + b.(2D +2Ddag) + 4b].R
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc),TMP2(Nv,Ndc)

      call DWilson(R,TMP,u,DAGGER,mass,czero)
      call DWilson(R,TMP2,u,.not.DAGGER,mass,czero)
      DR=2*add*(TMP+TMP2)+4*add*R
      call DWilson(TMP,TMP2,u,.not.DAGGER,mass,czero)
      DR=DR+(1+add)*TMP2

      return
      end subroutine DWB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module basicdiracopsmod
