      module testbasicdiracopsmod
      use pacc
      use arraysizes
      use numbers
      use basicdiracopsmod
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testBasicDiracOperators
      implicit none

!      call testDWOperators
c      call testDWGaugeSymmetry
!      call testShamirDiracOperators
      call testScaledDiracOperators

      return
      end subroutine testBasicDiracOperators
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
      
      mass=-one
      add=one/2

      print *,""
      print *,'test Wilson Dirac Operators'
      print *,""
      call setRVs(Nv*4,psi)

      print *,"u:",maxval(abs(u)),minval(abs(u))
      print *,"u:",u(1,:),u(2,:)

!     check if DN=-DNdag
      call DNaive(psi,TMP1,u,.false.,mass,add)
      call DNaive(psi,TMP2,u,.true.,mass,add)
      print *,"DN+DNdag:",maxval(abs(TMP1+TMP2))

!     check if W=Wdag
      call WilsonTerm(psi,TMP1,u,.false.,mass,add)
      call WilsonTerm(psi,TMP2,u,.true.,mass,add)
      print *,"W-Wdag:",maxval(abs(TMP1-TMP2))

!     check if D=Ddag
      call DWilson(psi,TMP1,u,.false.,mass,add)
      call DWilson(psi,TMP2,u,.true.,mass,add)
      print *,"D-Ddag:",maxval(abs(TMP1-TMP2))

!     check if DN.W=0
      call DNaive(psi,TMP1,u,.false.,mass,add)
      call WilsonTerm(TMP1,TMP2,u,.false.,mass,add)
      print *,"DN.W:",maxval(abs(TMP2))

!     check if W.DN=0
      call WilsonTerm(psi,TMP1,u,.false.,mass,add)
      call DNaive(TMP1,TMP3,u,.false.,mass,add)
      print *,"W.DN:",maxval(abs(TMP3))

      print *,"W.DN-DN.W:",maxval(abs(TMP3-TMP2))

c      call DNaive(psi,TMP1,u,.true.,mass,add)
c      call WilsonTerm(TMP1,TMP3,u,.false.,mass,add)
c      print *,"DN.W:",maxval(abs(TMP2))

!     check Mptr=>DWilson matches
      call DWilson(psi,Dpsi,u,.false.,mass,add)
      Mptr => DWilson
      call Mptr(psi,DdDpsi,u,.false.,mass,add)
      print *,"D-M(D)",maxval(abs(Dpsi-DdDpsi))

!     check if DdagD=DDdag
      call DWilson(psi,TMP1,u,.false.,mass,add)
      call DWilson(TMP1,TMP2,u,.true.,mass,add)
      call DWilson(psi,TMP1,u,.true.,mass,add)
      call DWilson(TMP1,TMP3,u,.false.,mass,add)
      print *,"Ddag.D-D.Ddag:",maxval(abs(TMP2-TMP3))

!     check Mptr => DdagD matches
      Mptr => DdagD
      call DWilson(psi,TMP1,u,.false.,mass,add)
      call DWilson(TMP1,TMP2,u,.true.,mass,add)
      call Mptr(psi,TMP1,u,.false.,mass,add)
      print *,"DdagD-Ddag.D:",maxval(abs(TMP1-TMP2))

      call DWilson(psi,TMP1,u,.false.,mass,add)
      call IDW(TMP1,TMP2,u,.false.,mass,add)
      print *,"DW.IDW:",maxval(abs(TMP2-psi))

      call DWilson(psi,TMP1,u,.true.,mass,add)
      call IDW(TMP1,TMP2,u,.true.,mass,add)
      print *,"DWdag.IDWdag:",maxval(abs(TMP2-psi))

      Mptr => DWilson
      call Mptr(psi,TMP1,u,.false.,mass,add)
      call IMnonsym(TMP1,TMP2,u,.false.,mass,add,Mptr)
      print *,"D.ID-I using IMnonsym:",maxval(abs(TMP2-psi))

      Mptr => DWilson
      call Mptr(psi,TMP1,u,.true.,mass,add)
      call IMnonsym(TMP1,TMP2,u,.true.,mass,add,Mptr)
      print *,"Ddag.IDdag-I using IMnonsym:",maxval(abs(TMP2-psi))
     
      Mptr => DdagD
      call Mptr(psi,TMP1,u,.false.,mass,add)
      call IM(TMP1,TMP2,u,.false.,mass,add,Mptr)
      print *,"DdagD.IM(DdagD)-I:",maxval(abs(TMP2-psi))

      Mptr => DdagD
      call Mptr(psi,TMP1,u,.true.,mass,add)
      call IM(TMP1,TMP2,u,.true.,mass,add,Mptr)
      print *,"DdagD(dag).IM(dag,DdagD)-I:",maxval(abs(TMP2-psi))

      Mptr => DdagDpC
      call Mptr(psi,TMP1,u,.false.,mass,add)
      call Mptr(psi,TMP2,u,.true.,mass,add)
      print *,"DdagDpC-DdagDpC(dag):",maxval(abs(TMP2-TMP1))

      Mptr => DdagDpC
      call Mptr(psi,TMP1,u,.false.,mass,add)
      call IM(TMP1,TMP2,u,.false.,mass,add,Mptr)
      print *,"DdagDpC.IM(DdagDpC)-I:",maxval(abs(TMP2-psi))

      Mptr => DdagDpC
      call Mptr(psi,TMP1,u,.true.,mass,add)
      call IM(TMP1,TMP2,u,.true.,mass,add,Mptr)
      print *,"DdagDpC(dag).IM(dag,DdagDpC)-I:",maxval(abs(TMP2-psi))

      call DWilson(psi,TMP1,u,.false.,mass,czero)
      call mGmu(TMP1,5)
      TMP2=psi
      call mGmu(TMP2,5)
      call DWilson(TMP2,TMP3,u,.true.,mass,czero)
      print *,"g5.D.g5-Ddag:",maxval(abs(TMP3-TMP1))
      call DWilson(psi,TMP1,u,.false.,mass,czero)
      call mGmu(TMP1,4)
      TMP2=psi
      call mGmu(TMP2,4)
      call DWilson(TMP2,TMP3,u,.true.,mass,czero)
      print *,"g3.D.g3-Ddag:",maxval(abs(TMP3-TMP1))

!      call G5DW(psi,TMP1,u,.false.,mass,add)
!      call IG5DW(TMP1,TMP2,u,.false.,mass,add)
!      print *,maxval(abs(realpart(TMP2-psi))),
!     &        maxval(abs(imagpart(TMP2-psi)))

!      Mptr => DWilson
!      call setRVs(Nv*4,psi)
!      call Mptr(psi,TMP1,u,.false.,mass,add)
!      call IM_BICGSTAB(TMP1,TMP2,u,.false.,mass,add,Mptr)
!      print *,'BICGSTAB error:',maxval(abs(TMP2-psi))

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
      subroutine testShamirDiracOperators
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
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr => NULL()
      
      mass=-one
      add=cone/two

      print *,""
      print *,'test Shamir Dirac Operators'
      print *,""

      call setRVs(Nv*4,psi)

      TMP1=psi
      call mGmu(TMP1,5)
      call DShamir(TMP1,TMP2,u,.false.,mass,add)
      call mGmu(TMP2,5)
      call DShamir(psi,TMP3,u,.true.,mass,add)
      print *,"g5.S.g5-Sdag:",maxval(abs(TMP3-TMP2))

      TMP1=psi
      call mGmu(TMP1,4)
      call DShamir(TMP1,TMP2,u,.false.,mass,add)
      call mGmu(TMP2,4)
      call DShamir(psi,TMP3,u,.true.,mass,add)
      print *,"g3.S.g3-Sdag:",maxval(abs(TMP3-TMP2))

      call DShamir(psi,TMP1,u,.false.,mass,czero)
      call IDS(TMP1,TMP2,u,.false.,mass,czero) ! only works for add=czero
      print *,"S.IS-I:",maxval(abs(TMP2-psi))

      call DShamir(psi,TMP1,u,.true.,mass,czero) 
      call IDS(TMP1,TMP2,u,.true.,mass,czero) ! only works for add=czero
      print *,"Sdag.ISdag-I:",maxval(abs(TMP2-psi))


      call DShamir(psi,TMP1,u,.false.,mass,cone/2)
      call DShamir(TMP1,TMP2,u,.true.,mass,cone/2)
      call DShamir(psi,TMP1,u,.true.,mass,cone/2)
      call DShamir(TMP1,TMP3,u,.false.,mass,cone/2)
      print *,"Sdag.S - S.Sdag:",maxval(abs(TMP3-TMP2))

      call SdagSpC(psi,TMP1,u,.false.,mass,cone/2)
      call SdagSpC(psi,TMP2,u,.true.,mass,cone/2)
      print *,"SdagS - (SdagS)dag:",maxval(abs(TMP2-TMP1))

!     check DsDs^-1=I
      call DShamir(psi,TMP1,u,.false.,mass,czero)
      call DShamir(TMP1,TMP2,u,.true.,mass,czero)
      call ISdagS(TMP2,Dpsi,u,.false.,mass,czero)
      print *,"ISdagS.Sdag.S-I:",maxval(abs(Dpsi-psi))

      Mptr=>DShamir
      call ISdagS(psi,TMP1,u,.false.,mass,czero)
      call IMdagM2(psi,TMP2,u,.false.,mass,czero,Mptr)
      print *,"ISdagS-IMdagM2:",maxval(abs(TMP2-TMP1))

      call DWB(psi,TMP1,u,.false.,mass,cone/two)
      call DWB(psi,TMP2,u,.true.,mass,cone/two)
      print *,'DWB-DWBdag:',maxval(abs(TMP1-TMP2))

      call DWB(psi,TMP1,u,.false.,mass,cone/two)
      call IDWB(TMP1,Dpsi,u,.false.,mass,cone/two)
      print *,'DWB.IDWB:',maxval(abs(Dpsi-psi))

      call DWB(psi,TMP1,u,.true.,mass,cone/two)
      call IDWB(TMP1,Dpsi,u,.true.,mass,cone/two)
      print *,'DWBdag.IDWBdag:',maxval(abs(Dpsi-psi))

      Mptr=>SdagSpC
      call SdagSpC(psi,TMP1,u,.false.,mass,cone/two)
      call ISdagSpC(TMP1,Dpsi,u,.false.,mass,cone/two)
      print *,'SdagSpC.ISdagSpC:',maxval(abs(Dpsi-psi))

      Mptr=>SdagSpC
      call SdagSpC(psi,TMP1,u,.false.,mass,cone/two)
      call IM(TMP1,Dpsi,u,.false.,mass,cone/two,Mptr)
      print *,'SdagSpC.IM(SdagSpC):',maxval(abs(Dpsi-psi))

      Mptr=>SdagSpC
      call SdagSpC(psi,TMP1,u,.true.,mass,cone/two)
      call ISdagSpC(TMP1,Dpsi,u,.true.,mass,cone/two)
      print *,'SdagSpC(dag).ISdagSpC(dag):',maxval(abs(Dpsi-psi))

      Mptr=>SdagSpC
      call SdagSpC(psi,TMP1,u,.true.,mass,cone/two)
      call IM(TMP1,Dpsi,u,.true.,mass,cone/two,Mptr)
      print *,'SdagSpCdag.IM(SdagSpC(dag)):',maxval(abs(Dpsi-psi))

      call IM(psi,TMP1,u,.false.,mass,cone/two,Mptr)
      call IM2(psi,TMP2,u,.false.,mass,cone/two,Mptr)
      print *,'IM(SdagSpC)-IM2(SdagSpC):',maxval(abs(TMP2-TMP1))

      call IM(psi,TMP1,u,.true.,mass,cone/two,Mptr)
      call IM2(psi,TMP2,u,.true.,mass,cone/two,Mptr)
      print *,'IM(SdagSpC(dag))-IM2(SdagSpC(dag)):',
     &                                       maxval(abs(TMP2-TMP1))

      return
      end subroutine testShamirDiracOperators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testScaledDiracOperators
      use numbers
      use gammas
      use indices
      use rvmodule
      use gaugefield
      use axbmodule1
      use IOmodule
      implicit none
      complex(prc) R(Nv,4),DR(Nv,4)
      complex(prc) TMP1(Nv,4),TMP2(Nv,4),TMP3(Nv,4)
      complex(prc) urvs(Nv)
      real(prc) alpha
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr => NULL()
      
      mass=-one

      call setRVs(Nv*4,R)
      call readThetaFromFile(3000,theta) 
      call coef(u,theta)

      Mptr => DdagDpC

      add=one/1000
      call DdagDpC(R,TMP1,u,.false.,mass,add) 
      call IM(TMP1,DR,u,.false.,mass,add,Mptr) 
      print *,maxval(abs(DR-R))

      Mptr => scaledDdagDpC

      add=one/1000
      alpha=one
      print *,"alpha:",alpha,"dj:",real(add)
      call scaledDdagDpC(R,TMP1,u,.false.,mass,add,alpha) 
      call IMscaled(TMP1,DR,u,.false.,mass,add,Mptr,alpha) 
      print *,maxval(abs(DR-R))

      add=one/100000
      alpha=one
      print *,"alpha:",alpha,"dj:",real(add)
      call scaledDdagDpC(R,TMP1,u,.false.,mass,add,alpha) 
      call IMscaled(TMP1,DR,u,.false.,mass,add,Mptr,alpha) 
      print *,maxval(abs(DR-R))

      add=one/10000000
      alpha=one
      print *,"alpha:",alpha,"dj:",real(add)
      call scaledDdagDpC(R,TMP1,u,.false.,mass,add,alpha) 
      call IMscaled(TMP1,DR,u,.false.,mass,add,Mptr,alpha) 
      print *,maxval(abs(DR-R))

      add=one/100000
      alpha=100*one
      print *,"alpha:",alpha,"dj:",real(add)
      call scaledDdagDpC(R,TMP1,u,.false.,mass,add,alpha) 
      call IMscaled(TMP1,DR,u,.false.,mass,add,Mptr,alpha) 
      print *,maxval(abs(DR-R))

      add=one/100000
      alpha=one/100
      print *,"alpha:",alpha,"dj:",real(add)
      call scaledDdagDpC(R,TMP1,u,.false.,mass,add,alpha) 
      call IMscaled(TMP1,DR,u,.false.,mass,add,Mptr,alpha) 
      print *,maxval(abs(DR-R))

      return
      end subroutine testScaledDiracOperators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testbasicdiracopsmod
