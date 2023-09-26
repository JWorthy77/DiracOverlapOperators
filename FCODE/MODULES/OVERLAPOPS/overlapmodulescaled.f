      module overlapmodulescaled
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use axbmodule1
      use basicdiracopsmod
      use VOLmodule
      use VOLNKmodule
      use VOLarchive
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine DOLWscaled(R,DR,u,DAGGER,mass,SRF,alpha)
      use options
      use gammas
      implicit none
      complex(prc),intent(in) :: R(Nv,4)
      complex(prc),intent(out) :: DR(Nv,4)
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      real(prc) alpha
      procedure(),pointer :: Sptr => NULL()
      complex(prc) :: TMP(Nv,4),VR(Nv,4),DL(Nv,4)

      Sptr => VOLscaled

      if (MTYPE.eq.1) then

        call Sptr(R,DR,u,DAGGER,-MDW,SRF,alpha)
        DR=(one+mass)/two*R + (one-mass)/two*DR

      elseif (MTYPE.eq.3) then
!       this is not the form corresponding (exactly) to the domain wall formulation
        if (.not.DAGGER) then
          call Sptr(R,VR,u,DAGGER,-MDW,SRF,alpha)
          DL=(R+VR)/two
          DR=zi*mass*(R-VR)/two
          call mGmu(DR,4)
          DR=DL+DR
        elseif (DAGGER) then
          call Sptr(R,VR,u,DAGGER,-MDW,SRF,alpha)
          DL=(R+VR)/two
          TMP=R
          call mGmu(TMP,4)
          call Sptr(TMP,VR,u,DAGGER,-MDW,SRF,alpha)
          DR=-zi*mass*(TMP-VR)/two
          DR=DL+DR
        endif

      endif

      return
      end subroutine DOLWscaled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLWscaled(RR,DR,u,DAGGER,mass,SRF,alpha) 
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      real(prc) alpha
      complex(prc) TMP(Nv,4)

      print *,"IDOLWscaled"

      if (.not. DAGGER) then
        call DOLWscaled(RR,TMP,u,.not.DAGGER,mass,SRF,alpha)
        call IDOLWsdagDOLWs(TMP,DR,u,DAGGER,mass,SRF,alpha)
      elseif (DAGGER) then
        call IDOLWsdagDOLWs(RR,TMP,u,DAGGER,mass,SRF,alpha)
        call DOLWscaled(TMP,DR,u,.not.DAGGER,mass,SRF,alpha)
      end if

      return
      end subroutine IDOLWscaled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLWsdagDOLWs(RR,DR,u,DAGGER,mass,SRF,malpha)
      implicit none
      integer,parameter :: kferm = 4*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      real(prc) malpha

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      print *,"IDOLWsdagDOLWs"

      itercg=0
c     initialise
      DR=RR
      call DOLWscaled(RR,x1,u,DAGGER,mass,SRF,malpha)
      call DOLWscaled(x1,x2,u,.not.DAGGER,mass,SRF,malpha)
      r=RR-x2
c      print *,r
c      stop
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call DOLWscaled(p,x1,u,DAGGER,mass,SRF,malpha)
        call DOLWscaled(x1,x2,u,.not.DAGGER,mass,SRF,malpha)
        alphan=sum(conjg(r)*r)
        alphad=sum(conjg(p)*x2)
        alpha=alphan/alphad

        DR=DR+alpha*p
        r=r-alpha*x2
        betan=sum(conjg(r)*r)
        print *,betan
        if (betan.lt.resid) goto 8
        beta=betan/alphan
        p=r+beta*p
      end do

8     continue
      print *,"IDOLWscaled:",itercg,niterc,betan
      return
      end subroutine IDOLWsdagDOLWs    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module overlapmodulescaled
