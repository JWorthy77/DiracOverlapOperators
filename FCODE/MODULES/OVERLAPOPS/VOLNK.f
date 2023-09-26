      module VOLNKmodule
!     VOL routines for Rajamanis spectra calcs
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use axbmodule1
      use basicdiracopsmod
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLNK(RR,DR,u,DAGGER,dwmass,SRF,Vptr,Mptr,Dptr)
!     Ddag.D = (1+V)(1+Vdag) ! Alt in my plots
      use options
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Vptr
      procedure(),pointer :: Mptr
      procedure(),pointer :: Dptr
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,nn,nd

      if (VERBOSE.eq.2) then ; print *,'VOLNK' ; end if

      if (associated(Mptr)) then
        print *,"VOLNK Ddag.D Wilson"
        call Vptr(RR,TMP1,u,.false.,dwmass,SRF,Mptr,Dptr)
c        call DOverlap(RR,TMP1,u,.false.,dwmass,SRF)
        TMP2=(RR+TMP1)/two
c        print *,"RR:",RR(1,1)
c        print *,"TMP2:",TMP2(1,1)
        call Vptr(TMP2,TMP1,u,.true.,dwmass,SRF,Mptr,Dptr)
c        call Doverlap(TMP2,TMP1,u,.true.,dwmass,SRF)
        DR=(TMP1+TMP2)/two
      else
        print *,"VOLNK Ddag.D Shamir"
        call Vptr(RR,TMP1,u,.false.,dwmass,SRF)
        TMP2=(RR+TMP1)/two
c        print *,"RR:",RR(1,1)
c        print *,"TMP2:",TMP2(1,1)
        call Vptr(TMP2,TMP1,u,.true.,dwmass,SRF)
        DR=(TMP1+TMP2)/two
      end if

      return
      end subroutine VOLNK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLNK2(RR,DR,u,DAGGER,dwmass,SRF,Vptr,Mptr,Dptr)
!     Rajamani and Nikhil's (2+V+Vdag)/4
      use options
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Vptr
      procedure(),pointer :: Mptr
      procedure(),pointer :: Dptr
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,nn,nd

      if (VERBOSE.eq.2) then ; print *,'VOLNK' ; end if

      if (associated(Mptr)) then
        print *,"VOLNK 2+V+Vd Wilson"
        call Vptr(RR,DR,u,.false.,dwmass,SRF,Mptr,Dptr)
        call Vptr(RR,TMP1,u,.true.,dwmass,SRF,Mptr,Dptr)
      else
        print *,"VOLNK 2+V+Vd Shamir"
        call Vptr(RR,DR,u,.false.,dwmass,SRF)
        call Vptr(RR,TMP1,u,.true.,dwmass,SRF)
      end if
      DR=(DR+TMP1+2*RR)/4

      return
      end subroutine VOLNK2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IVNK(RR,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      implicit none
      integer,parameter :: kferm = Ndc*Nv
      integer,parameter :: niterc=10*Nv
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Vptr
      procedure(),pointer :: Mptr
      procedure(),pointer :: Dptr

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call VOLNK(RR,x2,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call VOLNK(p,x2,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        alphan=sum(conjg(r)*r)
        alphad=sum(conjg(p)*x2)
        alpha=alphan/alphad

        DR=DR+alpha*p
        r=r-alpha*x2
        betan=sum(conjg(r)*r)
        if (betan.lt.resid) goto 8
        beta=betan/alphan
        p=r+beta*p
      end do

8     continue
c      print *,itercg,niterc,betan
      return
      end subroutine IVNK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module VOLNKmodule

