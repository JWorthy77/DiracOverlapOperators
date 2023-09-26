      module axbmodule1
      use pacc
      use arraysizes
      use numbers
!     solves for 2+1d operators
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMdagM_SJH(RR,DR,u,mass,add,Mptr)
!     solve [Dw(u)^dagger.Dw(u)] DR = RR
!     requires Dw(u)^dagger.Dw(u) to be symmetric, positive definite
      implicit none
      integer,parameter :: kferm = 4*Nv
      integer,parameter :: niterc=10*Nv
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      real(prc) mass
      complex(prc) add
      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)

      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

      DR=cone
      p=DR
      r=RR
      alpha=one

      call Mptr(p,x1,u,.false.,mass,add)
      call Mptr(x1,x2,u,.true.,mass,add)
      r=r-alpha*x2
      betan=sum(conjg(r)*r)
      betad=betan
      alphan=betan

      beta=zero
      p=r+beta*p
      if(betan.lt.resid) goto 8

      do nx=1,niterc
        itercg=itercg+1
        call Mptr(p,x1,u,.false.,mass,add)
        alphad=sum(conjg(x1)*x1)
        alpha=alphan/alphad
        DR=DR+alpha*p
        call Mptr(x1,x2,u,.true.,mass,add)
        r=r-alpha*x2
        betan=sum(conjg(r)*r) 
        beta=betan/betad
        betad=betan
        alphan=betan
        p=r+beta*p
        if(betan.lt.resid) goto 8
      end do
8     continue

      print *,itercg,niterc,betan
      return
      end subroutine IMdagM_SJH     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMnonsym(RR,DR,u,DAGGER,mass,add,Mptr)
!     solve Dw(u) DR = RR or DwD(u) DR = RR
      implicit none
      complex(prc) RR(Nv,Ndc),DR(Nv,Ndc)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr
      complex(prc) TMP(Nv,4)

c      if (.not. DAGGER) then
        call Mptr(RR,TMP,u,.not.DAGGER,mass,add)
        call IMdagM(TMP,DR,u,DAGGER,mass,add,Mptr)
c      elseif (DAGGER) then
c        call IMdagM(RR,TMP,u,DAGGER,mass,add,Mptr)
c        call Mptr(TMP,DR,u,.not.DAGGER,mass,add)
c      end if

      return
      end subroutine IMnonsym
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMdagM(RR,DR,u,DAGGER,mass,add,Mptr)
!     solve MdagM.DR = RR
!     requires Mdagger.M to be symmetric, positive definite
!     note DAGGER is for M, not for MdagM
      implicit none
      integer,parameter :: kferm = Ndc*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call Mptr(RR,x1,u,DAGGER,mass,add)
      call Mptr(x1,x2,u,.not.DAGGER,mass,add)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Mptr(p,x1,u,DAGGER,mass,add)
        call Mptr(x1,x2,u,.not.DAGGER,mass,add)
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
!      print *,itercg,niterc,betan
      return
      end subroutine IMdagM     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMdagM2(RR,DR,u,DAGGER,mass,add,Mptr)
!     solve MdagM.DR = RR
!     identical to IMdagM but allows for Mptr pointing to something using IMdagM
      implicit none
      integer,parameter :: kferm = Ndc*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call Mptr(RR,x1,u,DAGGER,mass,add)
      call Mptr(x1,x2,u,.not.DAGGER,mass,add)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Mptr(p,x1,u,DAGGER,mass,add)
        call Mptr(x1,x2,u,.not.DAGGER,mass,add)
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
!      print *,itercg,niterc,betan
      return
      end subroutine IMdagM2     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IM(RR,DR,u,DAGGER,mass,add,Mptr)
!     solve Mptr.DR = RR
!     requires M to be symmetric, positive definite
      implicit none
      integer,parameter :: kferm = Ndc*Nv
      integer,parameter :: niterc=10*Nv
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call Mptr(RR,x2,u,DAGGER,mass,add)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Mptr(p,x2,u,DAGGER,mass,add)
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
!      print *,"IM:",itercg,niterc,betan
      return
      end subroutine IM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMscaled(RR,DR,u,DAGGER,mass,add,Mptr,malpha)
!     solve Mptr.DR = RR
!     requires M to be symmetric, positive definite
      implicit none
      integer,parameter :: kferm = Ndc*Nv
      integer,parameter :: niterc=10*Nv
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr
      real(prc) :: malpha

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call Mptr(RR,x2,u,DAGGER,mass,add,malpha)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Mptr(p,x2,u,DAGGER,mass,add,malpha)
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
      print *,"IMscaled:",itercg,niterc,betan
      return
      end subroutine IMscaled
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IM2(RR,DR,u,DAGGER,mass,add,Mptr)
!     solve Mptr.DR = RR
!     same as IM, in case Mptr points to something using IM
      implicit none
      integer,parameter :: kferm = Ndc*Nv
      integer,parameter :: niterc=10*Nv
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call Mptr(RR,x2,u,DAGGER,mass,add)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Mptr(p,x2,u,DAGGER,mass,add)
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
c      print *,itercg,niterc,betan
      return
      end subroutine IM2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IM_BICGSTAB2(RR,DR,u,DAGGER,mass,add,Mptr)
!     solve Mptr.DR = RR
      use numbers
      implicit none
      integer,parameter :: kferm = 4*Nv
      integer,parameter :: niterc=10*Nv
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm),rhat0(kferm)
      complex(prc) rhat(kferm),v(kferm),best(kferm),errv(kferm)
      complex(prc) s(kferm),t(kferm),h(kferm)
      integer itercg
      integer nx,i
      real(prc) rhom1,rho,alpha,beta,omega
      real(prc) err

      itercg=0
      print *,'BICGSTAB2'

c     initialise
      DR=RR
      call Mptr(RR,x2,u,DAGGER,mass,add)
      r=RR-x2
      rhat0=r
      rhom1=one
      alpha=one
      omega=one
      v=czero
      p=czero

      do i=1,niterc
        itercg=itercg+1
        rho=dot_product(rhat0,r)
        beta=rho/rhom1*(alpha/omega)
        p=r+beta*(p-omega*v)
        call Mptr(p,v,u,DAGGER,mass,add)
        alpha=rho/dot_product(rhat0,v)
        h=DR+alpha*p
        call Mptr(h,best,u,DAGGER,mass,add)
        errv=RR-best
        err=dot_product(errv,errv)/kferm
        print *,err
        if (err .lt. resid) then
          DR=h
          return
        end if
        s=r-alpha*v
        call Mptr(s,t,u,DAGGER,mass,add)
        omega=dot_product(t,s)/dot_product(t,t)
        DR=h+omega*s
        call Mptr(DR,best,u,DAGGER,mass,add)
        errv=RR-best
        err=dot_product(errv,errv)/kferm
        print *,err
        if (err .lt. resid) then
          return
        end if
        r=s-omega*t
        rhom1=rho
      end do

      return
      end subroutine IM_BICGSTAB2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IM_BICGSTAB(RR,DR,u,DAGGER,mass,add,Mptr)
!     solve Mptr.DR = RR
      use numbers
      implicit none
      integer,parameter :: kferm = 4*Nv
      integer,parameter :: niterc=10*Nv
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm),rhat0(kferm)
      complex(prc) rhat(kferm),v(kferm),best(kferm),errv(kferm)
      complex(prc) s(kferm),t(kferm),h(kferm)
      integer itercg
      integer nx,i
      real(prc) rhom1,rho,alpha,beta,omega
      real(prc) err

      itercg=0
      print *,'BICGSTAB'

c     initialise, b=RR, let x0=DR=b
      DR=RR
      call Mptr(DR,x1,u,DAGGER,mass,add)
      r=RR-x1
      rhat0=r
      p=r
      do i=1,niterc
        itercg=itercg+1
        rhom1=dot_product(rhat0,r)
        call Mptr(p,v,u,DAGGER,mass,add)
        alpha=rhom1/dot_product(rhat0,v)
        s=r-alpha*v
        call Mptr(s,t,u,DAGGER,mass,add)
        omega=dot_product(t,s)/dot_product(t,t)
        DR=DR+alpha*p+omega*s
        r=s-omega*t
        err=dot_product(r,r)/kferm
        print *,err
        if (err .lt. resid) return
        rho=dot_product(rhat0,r)
        beta=rho/rhom1*(alpha/omega)
        p=r+beta*(p-omega*v)
      end do
      
      return
      end subroutine IM_BICGSTAB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module axbmodule1
