      module axbmoduleReal
      use pacc
      use arraysizes
      use numbers
!     solves real operators
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMnonsymReal(RR,DR,sigma,DAGGER,Mptr)
!     solve Mptr.DR = RR or Mptr^dag.DR = RR 
!     RR,DR,sigma are scalar
      implicit none
      real(prc) RR(Nv,Ndc),DR(Nv,Ndc)
      real(prc) sigma(Nv)
      logical DAGGER
      real(prc) mass
      procedure(),pointer :: Mptr
      real(prc) TMP(Nv,Ndc)

      if (.not. DAGGER) then
        call Mptr(RR,TMP,sigma,.not.DAGGER)
        call IMdagMReal(TMP,DR,sigma,DAGGER,Mptr)
      elseif (DAGGER) then
        call IMdagMReal(RR,TMP,sigma,DAGGER,Mptr)
        call Mptr(TMP,DR,sigma,.not.DAGGER)
      end if

      return
      end subroutine IMnonsymReal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMdagMReal(RR,DR,u,DAGGER,Mptr)
      implicit none
      integer,parameter :: kferm = Ndc*Nv
      integer,parameter :: niterc=Nv*10
      real(prc) RR(kferm), DR(kferm)
      real(prc) u(Nv)
      logical DAGGER
      procedure(),pointer :: Mptr

      real(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call Mptr(RR,x1,u,.false.)
      call Mptr(x1,x2,u,.true.)
      r=RR-x2
      p=r
      betan=sum(r*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Mptr(p,x1,u,.false.)
        call Mptr(x1,x2,u,.true.)
        alphan=sum(r*r)
        alphad=sum(p*x2)
        alpha=alphan/alphad

        DR=DR+alpha*p
        r=r-alpha*x2
        betan=sum(r*r)
        if (betan.lt.resid) goto 8
        beta=betan/alphan
        p=r+beta*p
      end do

8     continue
c      print *,itercg,niterc,betan
      return
      end subroutine IMdagMReal    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module axbmoduleReal
