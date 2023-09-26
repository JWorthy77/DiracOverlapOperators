      module VOLarchive
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use axbmodule1
      use basicdiracopsmod
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLGfac(RR,S,u,DAGGER,dwmass,SRF,Mptr,Dptr)
!     approximate Voverlap with general kernel using factored rational functions
      use options
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Mptr
      procedure(),pointer :: Dptr
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,nn,nd

      if (VERBOSE.eq.2) then ; print *,'VOLfac' ; end if

      associate(num => SRF%frf%num%zeros,
     &          denom => SRF%frf%denom%zeros)

      nn=size(num)
      nd=size(denom)

      TMP1=SRF%mult*RR
      do j=1,nn
!       numerator
        add = -num(j)
        call Mptr(TMP1,TMP2,u,DAGGER,dwmass,add)
!       denominator
        add = -denom(j)
        call IM(TMP2,TMP1,u,DAGGER,dwmass,add,Mptr)
      end do
      do j=nn+1,nd
!       denominator
        add = -denom(j)
        call IM(TMP1,TMP2,u,DAGGER,dwmass,add,Mptr)
        TMP1=TMP2
      end do
      call Dptr(TMP1,S,u,DAGGER,dwmass,czero)

      end associate
      return
      end subroutine VOLGfac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLfac(RR,S,u,DAGGER,dwmass,SRF)
!     approximate Voverlap using factored rational functions
      use options
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Mptr => NULL()
      procedure(),pointer :: Dptr => NULL()

      Mptr => DdagDpC
      Dptr => DWilson
      call VOLGfac(RR,S,u,DAGGER,dwmass,SRF,Mptr,Dptr)
      return
      end subroutine VOLfac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLrem(RR,S,u,DAGGER,dwmass,SRF)
      use options
!     approximate Voverlap using remez partial fraction
!     
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      complex(prc) add
      integer j,nn
      procedure(),pointer :: Mptr => NULL()

      associate(front => SRF%spf%front,
     &          num => SRF%spf%num,
     &          denom => SRF%spf%denom)

      nn=size(num)
      S=front*RR
      Mptr => DdagDpC
      do j=1,nn
         add = denom(j)
         call IM(RR,TMP1,u,DAGGER,dwmass,add,Mptr)
         S=S+num(j)*TMP1
      end do

      call DWilson(S,TMP1,u,DAGGER,dwmass,czero)
      S=TMP1
      end associate
      end subroutine VOLrem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLpoly(RR,S,u,DAGGER,dwmass,SRF)
      use options
!     approximate Voverlap using polynomial function
!     
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      procedure(),pointer :: Mptr => NULL()
      real(prc) mult

      associate(front => SRF%pfrf%front%coeffs)
      nn=size(front)
      
      call DdagD(RR,TMP1,u,DAGGER,dwmass,czero)
      S=front(2)*TMP1+front(1)*RR
      do j=3,nn
        call DdagD(TMP1,TMP2,u,DAGGER,dwmass,czero)
        TMP1=TMP2
        S=S+front(j)*TMP1
      end do
      call DWilson(S,TMP1,u,DAGGER,dwmass,czero)
      S=TMP1

      end associate
      return
      end subroutine VOLpoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Tpoly(RR,S,u,DAGGER,dwmass,SRF)
      use options
!     approximate Voverlap using polynomial function
!     
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      procedure(),pointer :: Mptr => NULL()
      real(prc) mult
      integer Dtype

      associate(front => SRF%pfrf%front%coeffs)
      nn=size(front)
      
      DTYPE=2

      if(DTYPE.eq.1)then
        call DWilson(RR,TMP2,u,.not.DAGGER,dwmass,czero)
        call DWilson(TMP2,TMP1,u,DAGGER,dwmass,czero)
      elseif(DTYPE.eq.2)then
        call DWilson(RR,TMP1,u,DAGGER,dwmass,czero)
      elseif(DTYPE.eq.3)then
        call DHermWilson(RR,TMP1,u,DAGGER,dwmass,czero)
      endif
      S=front(2)*TMP1+front(1)*RR
      do j=3,nn
        if(DTYPE.eq.1)then
          call DWilson(TMP1,TMP2,u,.not.DAGGER,dwmass,czero)
          call DWilson(TMP2,TMP1,u,DAGGER,dwmass,czero)
        elseif(DTYPE.eq.2)then
          call DWilson(TMP1,TMP2,u,DAGGER,dwmass,czero)
          TMP1=TMP2
        elseif(DTYPE.eq.3)then
          call DHermWilson(TMP1,TMP2,u,DAGGER,dwmass,czero)
          TMP1=TMP2
        endif
        S=S+front(j)*TMP1
      end do

      end associate
      return
      end subroutine Tpoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLepoly(RR,S,u,DAGGER,dwmass,SRF)
      use options
!     approximate Voverlap using extended polynomial function
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,pnstart,pnend
      procedure(),pointer :: Mptr => NULL()

      associate(ep => SRF%epoly, epfs => SRF%epoly%cfs)
      pnend=ep%ou-ep%ol+1
      pnstart=-ep%ol+1

c      print *,pnstart,pnend
c      print *,epfs(pnstart+1),epfs(pnstart)

      call DdagD(RR,TMP1,u,DAGGER,dwmass,czero)
      S=epfs(pnstart+1)*TMP1+epfs(pnstart)*RR

c      goto 1111

      do j=pnstart+2,pnend
        call DdagD(TMP1,TMP2,u,DAGGER,dwmass,czero)
        TMP1=TMP2
        S=S+epfs(j)*TMP1
      end do

c      goto 1111

      TMP1=RR
      do j=-ep%ol,1,-1
c        print *,j
        call IDdagD(TMP1,TMP2,u,DAGGER,dwmass,czero)
        TMP1=TMP2
        S=S+epfs(j)*TMP1
      end do


      call DWilson(S,TMP1,u,DAGGER,dwmass,czero)
      S=TMP1
1111  continue
      end associate
      return
      end subroutine VOLepoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLSpoly(RR,S,u,DAGGER,dwmass,SRF)
      use options
      use ratfuncs
!     approximate Voverlap using polynomial function
!     
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      procedure(),pointer :: Mptr => NULL()
      real(prc) mult

      associate(front => SRF%pfrf%front%coeffs)
      nn=size(front)
      
      
      call SdagSpCpoly(RR,TMP1,u,DAGGER,dwmass,czero)
      S=front(2)*TMP1+front(1)*RR
      do j=3,nn
        call SdagSpCpoly(TMP1,TMP2,u,DAGGER,dwmass,czero)
        TMP1=TMP2
        S=S+front(j)*TMP1
      end do
      call DShamirPoly(S,TMP1,u,DAGGER,dwmass,czero)
      S=TMP1

      end associate
      return
      end subroutine VOLSpoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLSepoly(RR,S,u,DAGGER,dwmass,SRF)
      use options
      use ratfuncs
!     approximate Voverlap using polynomial function
!     
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,pnstart,pnend

      associate(ep => SRF%epoly, epfs => SRF%epoly%cfs)
      pnend=ep%ou-ep%ol+1
      pnstart=-ep%ol+1
      print *,pnstart,pnend
      ! regular poly index = -ep%ou+1:nn      
      call SdagSpCpoly(RR,TMP1,u,DAGGER,dwmass,czero)
      S=epfs(pnstart+1)*TMP1+epfs(pnstart)*RR
      do j=pnstart+2,pnend
        call SdagSpCpoly(TMP1,TMP2,u,DAGGER,dwmass,czero)
        TMP1=TMP2
        S=S+epfs(j)*TMP1
      end do

      TMP1=RR
      do j=-ep%ol,1,-1
        print *,j
        call ISdagS(TMP1,TMP2,u,DAGGER,dwmass,czero)
        TMP1=TMP2
        S=S+epfs(j)*TMP1
      end do

      call DShamirPoly(S,TMP1,u,DAGGER,dwmass,czero)
      S=TMP1

      end associate
      return
      end subroutine VOLSepoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module VOLarchive
