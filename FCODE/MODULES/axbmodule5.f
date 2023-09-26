      module axbmodule5
!     for DOL using Vptr,Mptr, and Dptr
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use options
      use basicdiracopsmod
      use overlapmoduledev
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLop2(RR,DR,u,DAGGER,mass,SRF) 
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc),intent(in) :: SRF
      complex(prc) TMP(Nv,4)
      procedure(),pointer :: Vptr ! VOLGfac, VOLGpf, VOLGrem, VOLMpf
      procedure(),pointer :: Mptr ! DdagDpC, SdagSpC, H2pC
      procedure(),pointer :: Dptr ! DWilson, DShamir, DHermWilson

      Vptr=>VOLMpf
      if (dwkernel.eq.1) then
        Mptr=>SdagSpC
        Dptr=>DShamir
      elseif (dwkernel.eq.2) then
        Mptr=>DdagDpC
        Dptr=>DWilson
      endif

      if (.not. DAGGER) then
        call DOL(RR,TMP,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        call IDdagDOL(TMP,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      elseif (DAGGER) then
        call IDdagDOL(RR,TMP,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        call DOL(TMP,DR,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      end if

      return
      end subroutine IDOLop2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDdagDOL(RR,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      implicit none
      integer,parameter :: kferm = 4*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc),intent(in) :: SRF
      procedure(),pointer :: Vptr ! VOLGfac, VOLGpf, VOLGrem, VOLMpf
      procedure(),pointer :: Mptr ! DdagDpC, SdagSpC, H2pC
      procedure(),pointer :: Dptr ! DWilson, DShamir, DHermWilson

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call DOL(RR,x1,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      call DOL(x1,x2,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call DOL(p,x1,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        call DOL(x1,x2,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
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
      print *,itercg,niterc,betan
      return
      end subroutine IDdagDOL     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module axbmodule5
