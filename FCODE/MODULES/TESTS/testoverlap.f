      module testoverlapopsmod
      use pacc
      use arraysizes
      use numbers
      use options
      use rvmodule
      use ratfuncs
      use axbmodule1
      use gammas
      use overlapmoduledev
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testOverlapOperators
      implicit none

!      call testWilsonHTOverlap(20)
!      call testShamirHTOverlap(20)
      call testScaledHTOverlap()

      return
      end subroutine testOverlapOperators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testWilsonHTOverlap(Nht)
      use gaugefield
c      use axbmodule5
      implicit none
      integer Nht
      complex(prc) R(Nv,4),DR1(Nv,4),DR2(Nv,4)
      complex(prc) TMP(Nv,4)
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()
      procedure(),pointer :: Dptr => NULL()

      real(prc) lmin,lmax
      integer j,jmin,jmax

      call setRVs(Nv*4,R)
      call setHTcoeffs(Nht,SRF)

      dwkernel=2
      MTYPE=3
      call DOverlap(R,DR1,u,.false.,one/10,SRF)
      call DOLop2(R,DR2,u,.false.,one/10,SRF)
      print *,"DOverlap-DOL:",maxval(abs(DR2-DR1))
      call IDOLop2(DR2,DR1,u,.false.,one/10,SRF)
      print *,"DOL.IDOL:",maxval(abs(DR1-R))
      return

      print *,""
      print *,'test Wilson HT Overlap Operators'
      print *,""

      Mptr=>DdagDpC
      Dptr=>DWilson

      call VOLpf(R,DR1,u,.false.,-MDW,SRF)
      call VOLGpf(R,DR2,u,.false.,-MDW,SRF,Mptr,Dptr)
      print *,"VOLGpf-VOLpf:",maxval(abs(DR2-DR1))

      call VOLpf(R,DR1,u,.true.,-MDW,SRF)
      call VOLGpf(R,DR2,u,.true.,-MDW,SRF,Mptr,Dptr)
      print *,"VOLGpfdag-VOLpfdag:",maxval(abs(DR2-DR1))

!     G3.D.G3 = Ddag
      TMP=R
      call mGmu(TMP,4)
      call VOLpf(TMP,DR1,u,.false.,-MDW,SRF)
      call mGmu(DR1,4)
      call VOLpf(R,DR2,u,.true.,-MDW,SRF)
      print *,"G3.VOLpf.G3-VOLpf(dag):",maxval(abs(DR2-DR1))

!     G5.D.G5 = Ddag
      TMP=R
      call mGmu(TMP,5)
      call VOLpf(TMP,DR1,u,.false.,-MDW,SRF)
      call mGmu(DR1,5)
      call VOLpf(R,DR2,u,.true.,-MDW,SRF)
      print *,"G5.VOLpf.G5-VOLpf(dag):",maxval(abs(DR2-DR1))

!     M1 - G5.D.G5 = Ddag
      MTYPE=1
      TMP=R
      call mGmu(TMP,5)
      call DOverlap(TMP,DR1,u,.false.,one/10,SRF)
      call mGmu(DR1,5)
      call DOverlap(R,DR2,u,.true.,one/10,SRF)
      print *,"MTYPE=1: g5.DOL(.1).g5-DOL(.1,dag):",maxval(abs(DR2-DR1))

!     M1 - G3.D.G3 = Ddag
      MTYPE=1
      TMP=R
      call mGmu(TMP,4)
      call DOverlap(TMP,DR1,u,.false.,one/10,SRF)
      call mGmu(DR1,4)
      call DOverlap(R,DR2,u,.true.,one/10,SRF)
      print *,"MTYPE=1: g3.DOL(.1).g3-DOL(.1,dag):",maxval(abs(DR2-DR1))

!     G3.D.G3 = Ddag  
      MTYPE=3
      TMP=R
      call mGmu(TMP,4)
      call DOverlap(TMP,DR1,u,.false.,one/10,SRF)
      call mGmu(DR1,4)
      call DOverlap(R,DR2,u,.true.,one/10,SRF)
      print *,"MTYPE=3: g3.DOL(.1).g3-DOL(.1,dag):",maxval(abs(DR2-DR1))
      print *,"MTYPE=3: this is not meant to be g3-hermitian - mass part 
     & is anti-hermitian"

!     (Ddag)dag = D  
      MTYPE=3
      call DOverlap(R,TMP,u,.false.,one/10,SRF)
      call mGmu(DR1,4)
      call DOverlap(R,DR2,u,.true.,one/10,SRF)
      print *,"MTYPE=3: g3.DOL(.1).g3-DOL(.1,dag):",maxval(abs(DR2-DR1))

      return
      end subroutine testWilsonHTOverlap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testShamirHTOverlap(Nht)
      use gaugefield
      implicit none
      integer Nht
      complex(prc) R(Nv,4),DR1(Nv,4),DR2(Nv,4)
      complex(prc) DR3(Nv,4),DR4(Nv,4)
      complex(prc) TMP(Nv,4)
      complex(prc) TMP1(Nv,4),TMP2(Nv,4),TMP3(Nv,4),TMP4(Nv,4)
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()
      procedure(),pointer :: Dptr => NULL()

      real(prc) lmin,lmax,dwmass
      integer j,jmin,jmax

c      logical DAGGER
      complex(prc) add

c      DAGGER=.false.
      dwmass=-one
      add=cone/two

      call setRVs(Nv*4,R)
      call setHTcoeffs(Nht,SRF)

      print *,""
      print *,'test Shamir HT Overlap Operators'
      print *,""

c      call DWB(R,TMP1,u,.false.,dwmass,cone/two)
c      call IDWB(TMP1,DR1,u,.false.,dwmass,cone/two)
c      print *,'DWB.IDWB:',maxval(abs(DR1-R))

c      call DWB(R,TMP2,u,.true.,dwmass,cone/two)
c      call IDWB(TMP2,DR1,u,.true.,dwmass,cone/two)
c      print *,'DWBdag.IDWBdag:',maxval(abs(DR1-R))

c      print *,'DWB(dag)-DWB:',maxval(abs(TMP1-TMP2))

c      call DWilson(R,TMP1,u,.false.,dwmass,cone/two)
c      call DWilson(TMP1,TMP2,u,.true.,dwmass,cone/two)
c      call DWilson(R,TMP1,u,.true.,dwmass,cone/two)
c      call DWilson(TMP1,TMP3,u,.false.,dwmass,cone/two)
c      print *,'DdagD-DDdag:',maxval(abs(TMP2-TMP3))

!     check if DN=-DNdag
c      call DNaive(R,TMP1,u,.false.,dwmass,add)
c      call DNaive(R,TMP2,u,.true.,dwmass,add)
c      print *,"DN+DNdag:",maxval(abs(TMP1+TMP2))

!     check if W=Wdag
c      call WilsonTerm(R,TMP1,u,.false.,dwmass,add)
c      call WilsonTerm(R,TMP2,u,.true.,dwmass,add)
c      print *,"W-Wdag:",maxval(abs(TMP1-TMP2))

!     check if D=Ddag
c      call DWilson(R,TMP1,u,.false.,dwmass,add)
c      call DWilson(R,TMP2,u,.true.,dwmass,add)
c      print *,"D-Ddag:",maxval(abs(TMP1-TMP2))

!     check if DN.W=0
c      call DNaive(R,TMP1,u,.false.,dwmass,add)
c      call WilsonTerm(TMP1,TMP2,u,.false.,dwmass,add)
c      print *,"DN.W:",maxval(abs(TMP2))

!     check if W.DN=0
c      call WilsonTerm(R,TMP1,u,.false.,dwmass,add)
c      call DNaive(TMP1,TMP3,u,.false.,dwmass,add)
c      print *,"W.DN:",maxval(abs(TMP3))

c      print *,"W.DN-DN.W:",maxval(abs(TMP3-TMP2))

c      print *,"u:",maxval(abs(u)),minval(abs(u))
c      print *,"u:",u(1,:),u(2,:)

      Mptr=>SdagSpC
      call SdagSpC(R,TMP1,u,.false.,dwmass,cone/two)
      call ISdagSpC(TMP1,DR1,u,.false.,dwmass,cone/two)
      print *,'SdagSpC.ISdagSpC-I:',maxval(abs(DR1-R))

      call SdagSpC(R,TMP1,u,.true.,dwmass,cone/two)
      call ISdagSpC(TMP1,DR1,u,.true.,dwmass,cone/two)
      print *,"ISdagSpC.SdagSpC-I:",maxval(abs(DR1-R))

      Mptr => SdagSpC
      call IM(R,TMP2,u,.false.,dwmass,cone/two,Mptr)
      call IM(R,TMP4,u,.true.,dwmass,cone/two,Mptr)
      call ISdagSpC(R,TMP1,u,.false.,dwmass,cone/two)
      call ISdagSpC(R,TMP3,u,.true.,dwmass,cone/two)
      print *,"IMdag-ISdag:",maxval(abs(TMP1-TMP2))
      print *,"IM-IS:",maxval(abs(TMP3-TMP4))
      print *,"IMdag-IM:",maxval(abs(TMP2-TMP4))
      print *,"ISdag-IS:",maxval(abs(TMP1-TMP3))

      Mptr=>SdagSpC
      Dptr=>DShamir

      call VOLSpf(R,DR1,u,.false.,-MDW,SRF)
      call VOLGpf(R,DR2,u,.false.,-MDW,SRF,Mptr,Dptr)
      print *,"VOLGpf-VOLSpf:",maxval(abs(DR2-DR1))
      print *,"VOLGpf:",maxval(abs(DR2))
      print *,"VOLSpf:",maxval(abs(DR1))

      call VOLSpf(R,DR3,u,.true.,-MDW,SRF)
      call VOLGpf(R,DR4,u,.true.,-MDW,SRF,Mptr,Dptr)
      print *,"VOLGpf(dag)-VOLSpf(dag):",maxval(abs(DR3-DR4))

      print *,"VOLSpf-VOLSpf(dag):",maxval(abs(DR3-DR1))
      print *,"VOLGpf-VOLGpf(dag):",maxval(abs(DR4-DR2))

!     G3.V.G3 = Vdag
      TMP=R
      call mGmu(TMP,4)
      call VOLSpf(TMP,DR1,u,.false.,-MDW,SRF)
      call mGmu(DR1,4)
      call VOLSpf(R,DR2,u,.true.,-MDW,SRF)
      print *,"g3.VOLSpf.g3-VOLSpf(dag):",maxval(abs(DR2-DR1))

!     G5.V.G5 = Vdag
      TMP=R
      call mGmu(TMP,5)
      call VOLSpf(TMP,DR1,u,.false.,-MDW,SRF)
      call mGmu(DR1,5)
      call VOLSpf(R,DR2,u,.true.,-MDW,SRF)
      print *,"g5.VOLSpf.g5-VOLSpf(dag):",maxval(abs(DR2-DR1))

!     G5.D.G5 = Ddag
      TMP=R
      call mGmu(TMP,5)
      call DOLS(TMP,DR1,u,.false.,-MDW,SRF)
      call mGmu(DR1,5)
      call DOLS(R,DR2,u,.true.,-MDW,SRF)
      print *,"g5.DOLS.g5-DOLS(dag):",maxval(abs(DR2-DR1))

      return
      end subroutine testShamirHTOverlap
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testScaledHTOverlap()
      use gaugefield
c      use axbmodule5
      use IOmodule
      implicit none
      integer Nht
      complex(prc) R(Nv,4),DR0(Nv,4),DR1(Nv,4),DR2(Nv,4),DR3(Nv,4)
      complex(prc) DR4(Nv,4)
      complex(prc) TMP(Nv,4)
      type(sgnratfunc) :: SRF0,SRF1,SRF2,SRF3,SRF4
      type(sgnratfunc) :: ZRF0,ZRF1,ZRF2,ZRF3,ZRF4
      real(prc) alpha,lmin,lmax
      procedure(),pointer :: Vptr => NULL()
      procedure(),pointer :: Mptr => NULL()
      procedure(),pointer :: Dptr => NULL()

      call setRVs(Nv*4,R)
c      open(unit=143,file='R.dat',status='unknown',form='formatted')
c      read(143,*) R
c      close(143)

      call setHTcoeffs(10,SRF0)
      call setHTcoeffs(20,SRF1)
      call setHTcoeffs(40,SRF2)
      call setHTcoeffs(60,SRF3)
      call setHTcoeffs(80,SRF4)
      lmin=5d-4
      lmax=10d0
      call setZoloCoeffs(10,ZRF0,lmin,lmax)
      call setZoloCoeffs(30,ZRF1,5*lmin,lmax)
      call setZoloCoeffs(30,ZRF2,2*lmin,lmax)
      call setZoloCoeffs(40,ZRF3,lmin,lmax)
      call setZoloCoeffs(50,ZRF4,lmin,lmax)

      dwkernel=2
      MTYPE=3
      GAUGETYPE=2

c      call readThetaFromFile(500,theta)
      call readThetaFromFile(3000,theta)
      call coef(u,theta)

c      call DOverlap(R,DR0,u,.false.,0.005*one,ZRF3)
c      call DOLop2(R,DR1,u,.false.,0.005*one,ZRF3)
c      print *,"VOLM: ",maxval(abs(DR1-DR0))
      call DOLop2(R,DR1,u,.false.,0.005*one,ZRF3)
      call IDOLop4(DR1,DR2,u,.false.,0.005*one,ZRF3,ZRF3,ZRF3)
      print *,"D.ID: ",maxval(abs(DR2-R))
      print *,"lmin:",lmin,"lmax:",lmax
      return

      call DOverlap(R,DR1,u,.false.,0.005*one,ZRF1)
      call DOverlap(R,DR2,u,.false.,0.005*one,ZRF2)
      call DOverlap(R,DR3,u,.false.,0.005*one,ZRF3)
      call DOverlap(R,DR4,u,.false.,0.005*one,ZRF4)
      print *,"Z DOL-20-10:",maxval(abs(DR1-DR0))
      print *,"Z DOL-30-20:",maxval(abs(DR2-DR1))
      print *,"Z DOL-40-30:",maxval(abs(DR3-DR2))
      print *,"Z DOL-50-40:",maxval(abs(DR4-DR3))
      print *,"Z DOL-50-10:",maxval(abs(DR4-DR0))
      print *,"Z DOL-50-20:",maxval(abs(DR4-DR1))
      print *,"Z DOL-50-30:",maxval(abs(DR4-DR2))
      print *,"Z DOL-50-40:",maxval(abs(DR4-DR3))

      return

c      goto 2341

      call DOverlap(R,DR0,u,.false.,0.005*one,SRF0)
      call DOverlap(R,DR1,u,.false.,0.005*one,SRF1)
      call DOverlap(R,DR2,u,.false.,0.005*one,SRF2)
      call DOverlap(R,DR3,u,.false.,0.005*one,SRF3)
      call DOverlap(R,DR4,u,.false.,0.005*one,SRF4)
      print *,"DOL-20-10:",maxval(abs(DR1-DR0))
      print *,"DOL-40-20:",maxval(abs(DR2-DR1))
      print *,"DOL-60-40:",maxval(abs(DR3-DR2))
      print *,"DOL-80-60:",maxval(abs(DR4-DR3))
      print *,"DOL-80-10:",maxval(abs(DR4-DR0))
      print *,"DOL-80-20:",maxval(abs(DR4-DR1))
      print *,"DOL-80-40:",maxval(abs(DR4-DR2))
      print *,"DOL-80-60:",maxval(abs(DR4-DR3))

      alpha=one/2
      call DOLWscaled(R,DR1,u,.false.,0.005*one,SRF1,alpha)
      print *,"DOverlap-DOLWscaled(1/2):",maxval(abs(DR4-DR1))
      alpha=one/4
      call DOLWscaled(R,DR1,u,.false.,0.005*one,SRF1,alpha)
      print *,"DOverlap-DOLWscaled(1/4):",maxval(abs(DR4-DR1))
      alpha=2*one
      call DOLWscaled(R,DR1,u,.false.,0.005*one,SRF1,alpha)
      print *,"DOverlap-DOLWscaled(2):",maxval(abs(DR4-DR1))

      return

      alpha=2*one
      Vptr => VOLpf
      Mptr => DdagDpC
      Dptr => DWilson
c      call IDOLdagDOL(R,DR1,u,.false.,0.005*one,SRF1,Vptr,Mptr,Dptr)
      call IDOLWsdagDOLWs(R,DR1,u,.false.,0.005*one,SRF1,alpha)
      call DOLWscaled(DR1,DR2,u,.false.,0.005*one,SRF1,alpha)
      call DOLWscaled(DR2,DR1,u,.true.,0.005*one,SRF1,alpha)
      print *,"test:",maxval(abs(DR1-R))

      alpha=2*one
      call DOLWscaled(R,DR1,u,.false.,0.005*one,SRF1,alpha)
      call IDOLWscaled(DR1,DR2,u,.false.,0.005*one,SRF1,alpha)
      print *,"test:",maxval(abs(DR2-R))

2341  continue
      call readThetaFromFile(500,theta)
      call coef(u,theta)

      alpha=one
      open(unit=10,file="R.dat",status='unknown')
      write(10,*) R
      close(10)
      call IDOLWscaled(R,DR0,u,.false.,0.01*one,SRF0,alpha)
      open(unit=10,file="HT10.dat",status='unknown')
      write(10,*) DR0
      close(10)
      call IDOLWscaled(R,DR1,u,.false.,0.01*one,SRF1,alpha)
      open(unit=10,file="HT20.dat",status='unknown')
      write(10,*) DR1
      close(10)
      print *,"12x12 beta=0.3 HT20-HT10:",maxval(abs(DR1-DR0))
      call IDOLWscaled(R,DR2,u,.false.,0.01*one,SRF2,alpha)
      open(unit=10,file="HT40.dat",status='unknown')
      write(10,*) DR2
      close(10)
      print *,"12x12 beta=0.3 HT40-HT20:",maxval(abs(DR2-DR1))
      call IDOLWscaled(R,DR3,u,.false.,0.01*one,SRF3,alpha)
      open(unit=10,file="HT60.dat",status='unknown')
      write(10,*) DR3
      close(10)
      print *,"12x12 beta=0.3 HT60-HT40:",maxval(abs(DR3-DR2))
      print *,"12x12 beta=0.3 HT40-HT10:",maxval(abs(DR2-DR0))
      print *,"12x12 beta=0.3 HT40-HT20:",maxval(abs(DR2-DR1))
c      print *,"12x12 beta=0.3 HT:",maxval(abs(DR2-DR2))

c      alpha=one/2
c      call IDOLWscaled(R,DR2,u,.false.,0.005*one,SRF1,alpha)
c      print *,"16x16 beta=0.3 test:",maxval(abs(DR3-R))

c      alpha=one/2
c      call DOLWscaled(R,DR1,u,.false.,0.005*one,SRF1,alpha)
c      call IDOLWscaled(DR1,DR4,u,.false.,0.005*one,SRF1,alpha)
c      print *,"16x16 beta=0.3 test:",maxval(abs(DR4-R))

c      alpha=one
c      call DOLWscaled(R,DR1,u,.false.,0.005*one,SRF1,alpha)
c      call IDOLWscaled(DR1,DR5,u,.false.,0.005*one,SRF1,alpha)
c      print *,"16x16 beta=0.3 test:",maxval(abs(DR5-R))

c      call DOverlap(R,DR1,u,.false.,0.005*one,SRF1)
c      call IDOverlap(DR1,DR2,u,.false.,0.005*one,SRF1)
c      print *,"DOLW:",maxval(abs(DR2-R))

c      call DOL(R,DR1,u,.false.,0.005*one,SRF1)
c      call IDOL(DR1,DR2,u,.false.,0.005*one,SRF1)
c      print *,"DOLW:",maxval(abs(DR2-R))

      return
      end subroutine testScaledHTOverlap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testoverlapopsmod
