      module testVOLmod
      use pacc
      use arraysizes
      use numbers
      use options
      use rvmodule
      use ratfuncs
      use axbmodule1
      use gammas
      use VOLmodule
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testVOLoperators
      implicit none

!      call testVOLW(20)
!      call testVOLW(21)
      call testVOLS(21)

      return
      end subroutine testVOLoperators
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testVOLW(Nht)
      use gaugefield
c      use axbmodule5
      implicit none
      integer Nht
      complex(prc) R(Nv,4),DR1(Nv,4),DR2(Nv,4),DR3(Nv,4)
      complex(prc) DR4(Nv,4),DR5(Nv,4),DR6(Nv,4)
      complex(prc) TMP(Nv,4)
      type(sgnratfunc) :: SRF,SRF1,SRF2,SRF3,SRF4,SRF5,SRF6
      procedure(),pointer :: Mptr => NULL()
      procedure(),pointer :: Dptr => NULL()

      real(prc) lmin,lmax
      integer j,jmin,jmax
      real(prc) alpha

      call setRVs(Nv*4,R)
      call setHTcoeffs(Nht,SRF)

      MTYPE=3

      Mptr=>DdagDpC
      Dptr=>DWilson

      call VOLpf(R,DR1,u,.false.,-MDW,SRF)
      call VOLGpf(R,DR2,u,.false.,-MDW,SRF,Mptr,Dptr)
      call VOLMpf(R,DR3,u,.false.,-MDW,SRF,Mptr,Dptr)
      call VOLMWpf(R,DR4,u,.false.,-MDW,SRF)
      print *,"VOLGpf-VOLpf:",maxval(abs(DR2-DR1))
      print *,"VOLMpf-VOLpf:",maxval(abs(DR3-DR1))
      print *,"VOLMWpf-VOLpf:",maxval(abs(DR4-DR1))


      call VOLpf(R,DR1,u,.true.,-MDW,SRF)
      call VOLGpf(R,DR2,u,.true.,-MDW,SRF,Mptr,Dptr)
      call VOLMpf(R,DR3,u,.true.,-MDW,SRF,Mptr,Dptr)
      call VOLMWpf(R,DR4,u,.true.,-MDW,SRF)
      print *,"VOLGpfdag-VOLpfdag:",maxval(abs(DR2-DR1))
      print *,"VOLMpfdag-VOMpfdag:",maxval(abs(DR3-DR1))
      print *,"VOLMWpfdag-VOLMWpfdag:",maxval(abs(DR4-DR1))

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

!     test scaled VOL
      alpha=one
      call VOLpf(R,DR1,u,.false.,-MDW,SRF)
      call VOLscaled(R,DR2,u,.false.,-MDW,SRF,alpha)
      print *,"VOLscaled-VOLpf:",maxval(abs(DR2-DR1))

      call setHTcoeffs(2,SRF1)
      call setHTcoeffs(6,SRF2)
      call setHTcoeffs(10,SRF3)
      call setHTcoeffs(14,SRF4)
      call setHTcoeffs(20,SRF5)
      call setHTcoeffs(60,SRF6)

      call VOLpf(R,DR1,u,.false.,-MDW,SRF1)
      call VOLpf(R,DR2,u,.false.,-MDW,SRF2)
      call VOLpf(R,DR3,u,.false.,-MDW,SRF3)
      call VOLpf(R,DR4,u,.false.,-MDW,SRF4)
      call VOLpf(R,DR5,u,.false.,-MDW,SRF5)
      call VOLpf(R,DR6,u,.false.,-MDW,SRF6)
      print *,"VOLpf(60)-VOLpf(2):",maxval(abs(DR6-DR1))
      print *,"VOLpf(60)-VOLpf(6):",maxval(abs(DR6-DR2))
      print *,"VOLpf(60)-VOLpf(10):",maxval(abs(DR6-DR3))
      print *,"VOLpf(60)-VOLpf(14):",maxval(abs(DR6-DR4))
      print *,"VOLpf(60)-VOLpf(20):",maxval(abs(DR6-DR5))

      alpha=0.5*one
      call VOLscaled(R,DR1,u,.false.,-MDW,SRF5,alpha)
      print *,"VOLscaled-VOLpf:",maxval(abs(DR6-DR1))

      return
      end subroutine testVOLW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testVOLS(Nht)
      use gaugefield
      implicit none
      integer Nht
      complex(prc) R(Nv,4),DR1(Nv,4),DR2(Nv,4),DR3(Nv,4)
      complex(prc) DR4(Nv,4),DR5(Nv,4),DR6(Nv,4)
      complex(prc) TMP(Nv,4)
      type(sgnratfunc) :: SRF,SRF1,SRF2,SRF3,SRF4,SRF5,SRF6
      procedure(),pointer :: Mptr => NULL()
      procedure(),pointer :: Dptr => NULL()

      real(prc) lmin,lmax
      integer j,jmin,jmax
      real(prc) alpha

      call setRVs(Nv*4,R)
      call setHTcoeffs(Nht,SRF)

      MTYPE=3

      Mptr=>SdagSpC
      Dptr=>DShamir

      call VOLSpf(R,DR1,u,.false.,-MDW,SRF)
      call VOLGpf(R,DR2,u,.false.,-MDW,SRF,Mptr,Dptr)
      call VOLMpf(R,DR3,u,.false.,-MDW,SRF,Mptr,Dptr)
      call VOLMSpf(R,DR4,u,.false.,-MDW,SRF)
      print *,"VOLGpf-VOLSpf:",maxval(abs(DR2-DR1))
      print *,"VOLMpf-VOLSpf:",maxval(abs(DR3-DR1))
      print *,"VOLMSpf-VOLSpf:",maxval(abs(DR4-DR1))


      call VOLSpf(R,DR1,u,.true.,-MDW,SRF)
      call VOLGpf(R,DR2,u,.true.,-MDW,SRF,Mptr,Dptr)
      call VOLMpf(R,DR3,u,.true.,-MDW,SRF,Mptr,Dptr)
      call VOLMSpf(R,DR4,u,.true.,-MDW,SRF)
      print *,"VOLGpfdag-VOLSpfdag:",maxval(abs(DR2-DR1))
      print *,"VOLMpfdag-VOLSpfdag:",maxval(abs(DR3-DR1))
      print *,"VOLMSpfdag-VOLSpfdag:",maxval(abs(DR4-DR1))

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

!     test scaled VOL
      alpha=one
      call VOLpf(R,DR1,u,.false.,-MDW,SRF)
      call VOLscaled(R,DR2,u,.false.,-MDW,SRF,alpha)
      print *,"VOLscaled-VOLpf:",maxval(abs(DR2-DR1))

      call setHTcoeffs(2,SRF1)
      call setHTcoeffs(6,SRF2)
      call setHTcoeffs(10,SRF3)
      call setHTcoeffs(14,SRF4)
      call setHTcoeffs(20,SRF5)
      call setHTcoeffs(60,SRF6)

      alpha=4*one
      call VOLpf(R,DR1,u,.false.,-MDW,SRF1)
      call VOLpf(R,DR2,u,.false.,-MDW,SRF2)
      call VOLpf(R,DR3,u,.false.,-MDW,SRF3)
      call VOLpf(R,DR4,u,.false.,-MDW,SRF4)
      call VOLpf(R,DR5,u,.false.,-MDW,SRF5)
      call VOLpf(R,DR6,u,.false.,-MDW,SRF6)
      print *,"VOLpf(60)-VOLpf(2):",maxval(abs(DR6-DR1))
      print *,"VOLpf(60)-VOLpf(6):",maxval(abs(DR6-DR2))
      print *,"VOLpf(60)-VOLpf(10):",maxval(abs(DR6-DR3))
      print *,"VOLpf(60)-VOLpf(14):",maxval(abs(DR6-DR4))
      print *,"VOLpf(60)-VOLpf(20):",maxval(abs(DR6-DR5))

      alpha=0.5*one
      call VOLscaled(R,DR1,u,.false.,-MDW,SRF5,alpha)
      print *,"VOLscaled-VOLpf:",maxval(abs(DR6-DR1))

      return
      end subroutine testVOLS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testVOLmod
