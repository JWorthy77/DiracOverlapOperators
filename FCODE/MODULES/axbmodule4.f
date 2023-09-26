      module axbmodule4
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      use overlapmoduledev
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLOpts(RR,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr) 
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Vptr,Mptr,Dptr
      complex(prc) TMP(Nv,4)

      if (.not. DAGGER) then
        call DOL(RR,TMP,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        call IDOLdagDOLOpts(TMP,DR,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      elseif (DAGGER) then
        call IDOLdagDOLOpts(RR,TMP,u,DAGGER,mass,SRF,Vptr,Mptr,Dptr)
        call DOL(TMP,DR,u,.not.DAGGER,mass,SRF,Vptr,Mptr,Dptr)
      end if

      return
      end subroutine IMnonsymOpts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IDOLdagDOL(RR,DR,u,DAGGER,mass,Vptr,Mptr,Dptr)
      implicit none
      integer,parameter :: kferm = 4*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Vptr,Mptr,Dptr
      complex(prc) add

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
      end subroutine IDOLdagDOL     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module axbmodule4
