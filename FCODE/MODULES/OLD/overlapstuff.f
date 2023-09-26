      module overlapmodule
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use axbmodule1
      use wilsonmodule
      implicit none

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLKpf(Nht,RR,S,u,DAGGER,dwmass)
!     approximate V=Dw(Dw^dag.Dw)^{-1/2}, where Dw=Dw(-M) using 
!     use Kennedy's hyperbolic tanh partial fraction formula
      implicit none
      integer Nht
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      complex(prc) rootsM1(Nht),rootsP1(Nht)
      integer j
      procedure(),pointer :: Mptr => NULL()

      if (mod(Nht,2).eq.1) then
        print *,"Kennedy SGN only for even Nht"
        stop
      end if

      print *,'VOLKpf'

      Mptr => DdagDpC
      S=zero
      do j=0,Nht/2-1
        coef=tan((j+half)*pi/Nht)**2
        call IM(RR,TMP1,u,.true.,dwmass,coef,Mptr)
        S=S+(one+coef)*TMP1
      end do
      call DWilson(S,TMP1,u,.false.,dwmass,czero)
      S=two/Nht*TMP1

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLKfac(Nht,RR,S,u,DAGGER,dwmass)
!     approximate V=Dw(Dw^dag.Dw)^{-1/2}, where Dw=Dw(-M) using 
!     use Kennedy's hyperbolic tanh partial factored formula
      implicit none
      integer Nht
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) dwmass
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      complex(prc) rootsM1(Nht),rootsP1(Nht)
      integer j
      procedure(),pointer :: Mptr => NULL()

      if (mod(Nht,2).eq.1) then
        print *,"Kennedy SGN only for even Nht"
        stop
      end if

      print *,'VOLKfac'

      Mptr => DdagDpC
      S=RR
      do j=1,Nht/2-1
        coef=tan(j*pi/Nht)**2
        call Mptr(S,TMP1,u,.true.,dwmass,coef,Mptr)
        coef=tan((j+half)*pi/Nht)**2
        call IM(TMP1,S,u,.true.,dwmass,coef,Mptr)
      end do
      coef=tan(half*pi/Nht)**2
      call IM(S,TMP1,u,.true.,dwmass,coef,Mptr)
      call DWilson(TMP1,S,u,.false.,dwmass,czero)
      S=Nht*S

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Action_DOL(R,S,u,mass,SRF)
      use options
      use gammas
      implicit none
!     calculates S = R^T.DOverlap.R 
      complex(prc),intent(in) :: R(4*Nv)
      complex(prc),intent(out) :: S
      complex(prc),intent(in) :: u(Nv,3)
      real(prc) mass
      type(sgnratfunc) :: SRF
      complex(prc) :: DR(4*Nv)
      complex(prc) S2
      
      MTYPE=1
      call DOverlap(R,DR,u,.false.,mass,SRF)
c      call DWilson(R,DR,u,.false.,mass,czero)
      S=dot_product(R,DR)
      print *,"MTYPE",MTYPE
      print *,"S:",S
      MTYPE=3
      call DOverlap(R,DR,u,.false.,mass,SRF)
c      call DWilson(R2,DR,u,.false.,mass,czero)
      S2=dot_product(R,DR)
      print *,"MTYPE",MTYPE
      print *,"S2:",S2
      print *,"S2-S:",S2-S

      return
      end subroutine Action_DOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testVOL(ZOLO,Nht)
      use gammas
      use options
      use ratfuncs
      use gaugefield
      use arpackmodule
      use sgnmodule
      use rvmodule
      implicit none
      logical ZOLO
      integer Nht
      type(sgnratfunc) :: RF
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc) DR1(Nv,4),DR2(Nv,4)
      procedure(),pointer :: Sptr => NULL()
      complex(prc) mtest(Nv,4,200)
      real ev(1)
      real(prc) lmin,lmax
      integer j,jmin,jmax

      call setRVs(Nv*4,R)
      if (.NOT.ZOLO) then
        call setHTcoeffs(Nht,RF)
      else
        call calceigs('SM',1,ev,1,Nv)
        print *,'ev min:',ev
        lmin=ev(1)
        call calceigs('LM',1,ev,1,Nv)
        print *,'ev max:',ev
        lmax=ev(1)
        call setZoloCoeffs(Nht,RF,lmin,lmax)
      endif
      call VOLpf(R,DR1,u,.false.,-MDW,RF)
      call VOLfac(R,DR2,u,.false.,-MDW,RF)
      print *,'VOLpf-VOLfac:',maxval(abs(DR1-DR2))
      call VOLpf(R,DR1,u,.true.,-MDW,RF)
      call VOLfac(R,DR2,u,.true.,-MDW,RF)
      print *,'VOLpf^dagger-VOLfac^dagger:',maxval(abs(DR1-DR2))
      call VOLfac(R,DR1,u,.false.,-MDW,RF)
      call VOLfac(R,DR2,u,.true.,-MDW,RF)
      print *,'VOLfac-VOLfac^dagger:',maxval(abs(DR1-DR2))

      return
      end subroutine testVOL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testOverlap(Nht)
      use options
      use gaugefield
      use gammas
      use rvmodule
      use axbmodule2
      use sgnmodule
      implicit none
      integer Nht
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc) TMP(Nv,4)
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()
      type(ioptions)  iopts
      complex(prc) S

      MTYPE=1     
      call setHTcoeffs(Nht,SRF)
      Mptr => VOLpf
      iopts%SRF=SRF
      iopts%mass=-MDW
      R=czero
      R(1,1)=cone
      call Mptr(R,TMP,u,.false.,-MDW,SRF)
      Mptr => VOLpf
      iopts%SRF=SRF
      iopts%mass=-MDW
      call IMnonsymOpts(TMP,DR,u,.false.,Mptr,iopts)
      print *,'VOLfac.IVOLfac',maxval(abs(R-DR))


      call SGNfactor(R,TMP,u,.false.,-MDW,SRF)
      call mGmu(TMP,5)
      TMP=(half+baremass/two)*R+(half-baremass/two)*TMP
      call DOverlap(R,DR,u,.false.,baremass,SRF)
      print *,'DOL-(1/2+1/2G5.SGN)',maxval(abs(TMP-DR))

      call DOverlap(R,TMP,u,.false.,baremass,SRF)
      Mptr => DOverlap
      iopts%SRF=SRF
      iopts%mass=baremass
      call IMnonsymOpts(TMP,DR,u,.false.,Mptr,iopts)
      print *,'DOL.IDOL',maxval(abs(R-DR))

      call setRVs(Nv*4,R)     
      call DOverlap(R,TMP,u,.false.,baremass,SRF)
      call IDOverlap(TMP,DR,u,.false.,baremass,SRF)
      print *,'DOL.IDOL',maxval(abs(R-DR))


      call setRVs(Nv*4,R)     
      call SGNfactor(R,TMP,u,.false.,-MDW,SRF)
      call mGmu(TMP,5)
      TMP=(half+baremass/two)*R+(half-baremass/two)*TMP
      call DOverlap(R,DR,u,.false.,baremass,SRF)
      print *,'DOL-(1/2+1/2G5.SGN)',maxval(abs(TMP-DR))

      call DOverlap(R,TMP,u,.true.,baremass,SRF)
      call IDOverlap(TMP,DR,u,.true.,baremass,SRF)
      print *,'DOL^dagger.IDOL^dagger',maxval(abs(R-DR))

      return
      end subroutine testOverlap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testOLAction(Nht,R)
      use options
      use gaugefield
      use rvmodule
      implicit none
      integer Nht
      complex(prc) R(Nv,4)
      type(sgnratfunc) :: SRF
      complex(prc) S

      call setHTcoeffs(Nht,SRF)
c      call setRVs(Nv*4,R)     
      call Action_DOL(R,S,u,baremass,SRF)

      return
      end subroutine testOLAction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testOLGaugeSymmetry
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
      integer idx,Nht
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()

      Mptr=>DOverlap
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
      end subroutine testOLGaugeSymmetry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcOLChiralError(R,DR,u,SRF)
!     2.g5.err = g5.D+D.g5-2.D.g5.D
      use numbers
      use gammas
      implicit none
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      type(sgnratfunc) :: SRF
      complex(prc) TMP(Nv,4),TMP1(Nv,4),TMP2(Nv,4)

!     g5.D
      call DOverlap(R,TMP1,u,.false.,zero,SRF)
      call mGmu(TMP1,4)
!     D.g5
      TMP=R
      call mGmu(TMP,4)
      call DOverlap(TMP,TMP2,u,.false.,zero,SRF)
!     2.D.g5.D
      call DOverlap(TMP1,TMP,u,.false.,zero,SRF)
      TMP=2*TMP
!     sum
      DR=TMP1+TMP2-TMP

      return
      end subroutine calcOLChiralError
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module overlapmodule
