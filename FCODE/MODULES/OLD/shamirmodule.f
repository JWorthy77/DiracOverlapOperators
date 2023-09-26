      module shamirmodule
      use pacc
      use arraysizes
      use numbers
      use wilsonmodule
      implicit none
      contains
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

      call DWilson(R,TMP,u,DAGGER,mass,add)
      call IDW(TMP,DR,u,DAGGER,mass,ctwo)
      DR=DR+add*R

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
      subroutine IHdagH(R,DR,u,DAGGER,mass,add)
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
      subroutine HdagHpC(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = DWilson*R
      complex(prc) R(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,NDT)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) TMP(Nv,Ndc)

      call DShamir(R,TMP,u,DAGGER,mass,czero)
      call DShamir(TMP,DR,u,.NOT.DAGGER,mass,czero)
      DR=DR+add*R

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine HdagHpCpoly(R,DR,u,DAGGER,mass,add)
      implicit none
c     calculates DR = DWilson*R
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
      end subroutine HdagHpCpoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLSpf(RR,S,u,DAGGER,dwmass,SRF)
      use options
      use ratfuncs
      use axbmodule1
!     approximate Voverlap using partial fraction rational functions
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
      complex(prc) cmult

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      nf=nn-nd+1
      S=zero
      if (nf.eq.1) then ! it should only ever be 0 or 1
!       b0.RR
        S=front(1)*RR
      end if

      Mptr => HdagHpC
      do j=1,nd
         add = -denom(j)
         call IM(RR,TMP1,u,DAGGER,dwmass,add,Mptr)
         TMP2 = pf(j)*TMP1
         S=S+TMP2
      end do

      call DShamir(S,TMP1,u,DAGGER,dwmass,czero)
      S=SRF%mult*TMP1

      end associate
      return
      end subroutine VOLSpf
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
      
      
      call HdagHpCpoly(RR,TMP1,u,DAGGER,dwmass,czero)
      S=front(2)*TMP1+front(1)*RR
      do j=3,nn
        call HdagHpCpoly(TMP1,TMP2,u,DAGGER,dwmass,czero)
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
      call HdagHpCpoly(RR,TMP1,u,DAGGER,dwmass,czero)
      S=epfs(pnstart+1)*TMP1+epfs(pnstart)*RR
      do j=pnstart+2,pnend
        call HdagHpCpoly(TMP1,TMP2,u,DAGGER,dwmass,czero)
        TMP1=TMP2
        S=S+epfs(j)*TMP1
      end do

      TMP1=RR
      do j=-ep%ol,1,-1
        print *,j
        call IHdagH(TMP1,TMP2,u,DAGGER,dwmass,czero)
        TMP1=TMP2
        S=S+epfs(j)*TMP1
      end do

      call DShamirPoly(S,TMP1,u,DAGGER,dwmass,czero)
      S=TMP1

      end associate
      return
      end subroutine VOLSepoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLS(R,DR,u,DAGGER,mass,SRF)
      use options
      use ratfuncs
      use axbmodule1
      use gammas
!     approximate Voverlap for Shamir kernel using partial fraction rational functions
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Sptr => NULL()
      complex(prc) :: TMP(Nv,4),VR(Nv,4),DL(Nv,4),IMG3R(Nv,4)

      if (VERBOSE.gt.2) print *,'DOLS'
      Sptr=>VOLSpf
      if (MTYPE.eq.1) then

        call Sptr(R,TMP,u,DAGGER,-MDW,SRF)
        DR=(half+mass/two)*R+(half-mass/two)*TMP

      elseif (MTYPE.eq.3) then

        if (.not.DAGGER) then
          call Sptr(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          DR=zi*mass*(R-VR)/two
          call mGmu(DR,4)
          DR=DL+DR
        elseif (DAGGER) then
          call Sptr(R,VR,u,DAGGER,-MDW,SRF)
          DL=(R+VR)/two
          TMP=R
          call mGmu(TMP,4)
          call Sptr(TMP,VR,u,DAGGER,-MDW,SRF)
          DR=-zi*mass*(TMP-VR)/two
          DR=DL+DR
        endif

      endif
      return
      end subroutine DOLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLS(R,DR,u,DAGGER,mass,SRF)
      use options
      use axbmodule2
      implicit none
!     solves DOL.DR = R 
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()
      type(ioptions)  iopts

      if (VERBOSE.eq.2) then ; print *,"IDOLS" ; endif
      iopts%SRF=SRF
      iopts%mass=mass
      Mptr => DOLS
      call IMnonsymOpts(R,DR,u,DAGGER,Mptr,iopts)

      return
      end subroutine IDOLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testDSOperators
      use rvmodule
      use gaugefield
      use options
      use axbmodule1
      use ratfuncs
      implicit none
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc) TMP1(Nv,4),TMP2(Nv,4),TMP3(Nv,4)
      integer idirac,iv
      complex(prc) add
      real(prc) mass
      procedure(),pointer :: Mptr => NULL()
      type(sgnratfunc) :: SRF
      
      print *,'test Shamir Dirac Operators'
      call setRVs(Nv*4,R)

      Mptr=>DShamir
      call Mptr(R,TMP1,u,.false.,baremass,cone)
      call IMnonsym(TMP1,DR,u,.false.,baremass,cone,Mptr)
      print *,"DShamir",maxval(abs(R-DR))
      Mptr=>DWilson
      call Mptr(R,TMP1,u,.false.,baremass,cone)
      call IMnonsym(TMP1,DR,u,.false.,baremass,cone,Mptr)
      print *,"DWilson",maxval(abs(R-DR))

      return
      end subroutine testDSOperators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testOLShamirGaugeSymmetry
      use numbers
      use gammas
      use indices
      use rvmodule
      use gaugefield
      use options
      use ratfuncs
      implicit none
      complex(prc) psi(Nv,4),Dpsi(Nv,4)
      complex(prc) psit(Nv,4),Dpsit(Nv,4)
      complex(prc) urvs(Nv)
      real(prc) alpha(Nv)
      complex(prc) pbDp,pbDpt
      integer idx,Nht
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()

      Mptr=>DOLS
      Nht=20
      call setHTcoeffs(Nht,SRF)
!     set gauge transform field alpha
      call setRVs(Nv,urvs)
      alpha=real(urvs)
!     set phi
      call setRVs(Nv*4,psi)
!     make phibar D phi
      call Mptr(psi,Dpsi,u,.false.,SRF)
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
      call Mptr(psit,Dpsit,u,.false.,SRF)
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
      end subroutine testOLShamirGaugeSymmetry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module shamirmodule
