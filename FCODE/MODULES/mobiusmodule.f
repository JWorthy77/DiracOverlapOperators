      module mobiusmodule
      use pacc
      use arraysizes
      use numbers
      use wilsonmodule
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DWScaled(R,DR,u,DAGGER,mass,mult)
      implicit none
!     calculates DR = (2+mult.Dw)*R 
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) mult
      complex(prc) :: TMP(Nv,4)

      call DWilson(R,TMP,u,DAGGER,mass,czero)
      DR=two*R+mult*TMP

      return
      end subroutine DWScaled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DMobius(R,DR,u,DAGGER,mass,add)
      use options
      use axbmodule1
      implicit none
c     calculates DR = DShamir*R = (2+Dw)^-1*Dw*R
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      complex(prc) :: TMP(Nv,4)
      procedure(),pointer :: Mptr
      real(prc) c
      complex(prc) d

      c=amob+bmob
      d=cmplx(amob-bmob)

      call DWilson(R,TMP,u,DAGGER,mass,czero)
      TMP=c*TMP
      Mptr=>DWScaled
      call IMnonsym(TMP,DR,u,DAGGER,mass,d,Mptr)

      return
      end subroutine DMobius
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DSFreaky(R,DR,u,DAGGER,mass,mult)
      implicit none
c     calculates DR = [Dw'Dw+(2+Dw').mult.(2+Dw)].R 
c     DSF'=DSF
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) mult
      complex(prc) :: TMP(Nv,4)

      call DWilson(R,TMP,u,.false.,mass,ctwo)
      TMP=mult*TMP
      call DWilson(TMP,DR,u,.true.,mass,ctwo)

      call DdagD(R,TMP,u,.false.,mass,czero)
      DR=DR+TMP

      return
      end subroutine DSFreaky
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DMFreaky(R,DR,u,DAGGER,mass,mult)
      use options
      implicit none
c     calculates DR = [Dw'Dw+(2+Dw').mult.(2+Dw)].R 
c     DSF'=DSF
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) mult
      complex(prc) :: TMP(Nv,4),TMP2(Nv,4)
      real(prc) c
      complex(prc) d

      c=amob+bmob
      d=cmplx(amob-bmob)

      call DWScaled(R,TMP,u,.false.,mass,d)
      TMP=mult*TMP
      call DWScaled(TMP,DR,u,.true.,mass,d)

      call DdagD(R,TMP,u,.false.,mass,czero)
      DR=DR+c*c*TMP

      return
      end subroutine DMFreaky
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine G5DS(R,DR,u,DAGGER,mass,add)
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
        call DShamir(R,DR,u,DAGGER,mass,add)
        call mGmu(DR,5)
      else
        TMP=R
        call mGmu(TMP,5)
        call DShamir(TMP,DR,u,DAGGER,mass,add)
      end if
      
      return
      end subroutine G5DS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDSFreaky(RR,DR,u,DAGGER,mass,add)
      use axbmodule1
      implicit none
c     calculates DR = IDShamir*R = Dw^-1.(2+Dw)*R
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr

      Mptr=>DSFreaky
      check pointers - IMdagM cannot have a pointer which uses IMdagM (for example)
      

      call IMnonsym(RR,DR,u,DAGGER,mass,add,Mptr)

      return
      end subroutine IDSFreaky
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine VOLMpf(RR,S,u,DAGGER,SRF)
      use options
      use ratfuncs
      use axbmodule1
!     approximate Voverlap for Mobius kernel using partial fraction rational functions
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) add,coef
      complex(prc) TMP1(Nv,4),TMP2(Nv,4)
      integer j,p,nn,nd,nf
      procedure(),pointer :: Mptr => NULL()
      real(prc) mult
      complex(prc) cmult,d

      print *,'VOLMpf'

      d=cmplx(amob-bmob,zero,prc)

      associate(front => SRF%pfrf%front%coeffs,
     &          num => SRF%frf%num%zeros,
     &          denom => SRF%pfrf%denom%zeros,
     &          pf => SRF%pfrf%pf)

      nn=size(num)
      nd=size(denom)
      nf=nn-nd+1
      S=zero
      if (nf.ge.1) then
        do j=0,nf-1
          mult=front(j+1)
!         T = m.(DS'.Ds)^j.RR
          TMP1=mult*RR
          do p=1,j
            print *,"implement DMdagMS"
            stop
c            call DMdagDM(TMP1,TMP2,u,DAGGER,-MDW,czero)
            TMP1=TMP2
          end do
          S=S+TMP1
        end do
      end if

      do j=1,nd
         cmult = -denom(j)
         call DWScaled(RR,TMP1,u,.not.DAGGER,-MDW,d)
         Mptr=>DMFreaky
         call IMnonsym(TMP1,TMP2,u,DAGGER,-MDW,cmult,Mptr)
         call DWScaled(TMP2,TMP1,u,DAGGER,-MDW,d)
         S = S + pf(j)*TMP1
      end do

      call DMobius(S,TMP1,u,DAGGER,-MDW,czero)
      S=SRF%mult*TMP1

      end associate

      return
      end subroutine VOLMpf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLM(RR,S,u,DAGGER,SRF)
      use options
      use ratfuncs
      use axbmodule1
!     approximate Voverlap for Shamir kernel using partial fraction rational functions
      implicit none
      complex(prc) RR(Nv,4),S(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) TMP(Nv,4)

      print *,'DOLM'
      call VOLMpf(RR,TMP,u,DAGGER,SRF)
      S=(half+baremass/two)*RR+(half-baremass/two)*TMP

      return
      end subroutine DOLM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
      Mptr=>DMobius
      call Mptr(R,TMP1,u,.false.,baremass,cone)
      call IMnonsym(TMP1,DR,u,.false.,baremass,cone,Mptr)
      print *,"DMobius",maxval(abs(R-DR))
      Mptr=>DWilson
      call Mptr(R,TMP1,u,.false.,baremass,cone)
      call IMnonsym(TMP1,DR,u,.false.,baremass,cone,Mptr)
      print *,"DWilson",maxval(abs(R-DR))
      call DWScaled(R,TMP1,u,.false.,baremass,cone)
      call DWilson(R,DR,u,.false.,baremass,ctwo)
      print *,"DWScaled",maxval(abs(TMP1-DR))
      call DMobius(R,TMP1,u,.false.,baremass,czero)
      call DShamir(R,DR,u,.false.,baremass,czero)
      print *,"DMS",maxval(abs(TMP1-DR))

      Mptr=>DSFreaky
      call Mptr(R,TMP1,u,.false.,baremass,cone)
      call IMnonsym(TMP1,DR,u,.false.,baremass,cone,Mptr)
      print *,"DSFreaky",maxval(abs(R-DR))
      Mptr=>DMFreaky
      call Mptr(R,TMP1,u,.false.,baremass,cone)
      call IMnonsym(TMP1,DR,u,.false.,baremass,cone,Mptr)
      print *,"DMFreaky",maxval(abs(R-DR))

      call setHTcoeffs(8,SRF)
      call VOLSpf(R,DR,u,.false.,-MDW,SRF)
      call VOLMpf(R,TMP1,u,.false.,SRF)
      print *,"VMS",maxval(abs(DR-TMP1))

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
      end module mobiusmodule
