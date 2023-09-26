      module axbmodule3
      use pacc
      use arraysizes
      use numbers
      use options
      use ratfuncs
      use axbmodule2
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IM5nonsymOpts(RR,DR,u,DAGGER,Mptr,iopts) ! Mptr like DPF
      implicit none
      complex(prc) RR(Nv,4,Ls),DR(Nv,4,Ls)
      complex(prc) u(Nv,3)
      logical DAGGER
      procedure(),pointer :: Mptr
      type(ioptions) iopts
      complex(prc) TMP(Nv,4,Ls)

      if (.not. DAGGER) then
        call Mptr(RR,TMP,u,.not.DAGGER,iopts%SRF)
        call IM5dagM5Opts(TMP,DR,u,DAGGER,Mptr,iopts)
      elseif (DAGGER) then
        call IM5dagM5Opts(RR,TMP,u,DAGGER,Mptr,iopts)
        call Mptr(TMP,DR,u,.not.DAGGER,iopts%SRF)
      end if

      return
      end subroutine IM5nonsymOpts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IM5dagM5Opts(RR,DR,u,DAGGER,Mptr,iopts)
      implicit none
      integer,parameter :: kferm = 4*Nv*Ls
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      procedure(),pointer :: Mptr
      type(ioptions) iopts

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

c     initialise
      DR=RR
      call Mptr(RR,x1,u,DAGGER,iopts%SRF)
      call Mptr(x1,x2,u,.not.DAGGER,iopts%SRF)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Mptr(p,x1,u,DAGGER,iopts%SRF)
        call Mptr(x1,x2,u,.not.DAGGER,iopts%SRF)
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
      end subroutine IM5dagM5Opts     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IM5dagM5OptsVS(kferm,RR,DR,u,DAGGER,Mptr,iopts)
      implicit none
      integer kferm
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      procedure(),pointer :: Mptr
      type(ioptions) iopts

      integer niterc
      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      niterc=kferm*10
      itercg=0

c     initialise
      DR=RR
      call Mptr(RR,x1,u,DAGGER,iopts%SRF)
      call Mptr(x1,x2,u,.not.DAGGER,iopts%SRF)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Mptr(p,x1,u,DAGGER,iopts%SRF)
        call Mptr(x1,x2,u,.not.DAGGER,iopts%SRF)
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
      print *,itercg,niterc,betan
      return
      end subroutine IM5dagM5OptsVS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IM5_BICGSTAB2(RR,DR,u,DAGGER,Mptr,SRF)
!     solve Mptr.DR = RR
      use numbers
      implicit none
      integer,parameter :: kferm = 4*Nv*Ls
      integer,parameter :: niterc=10*Nv*Ls
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr
      type(sgnratfunc) :: SRF

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm),rhat0(kferm)
      complex(prc) rhat(kferm),v(kferm),best(kferm),errv(kferm)
      complex(prc) s(kferm),t(kferm),h(kferm)
      integer itercg
      integer nx,i
      real(prc) rhom1,rho,alpha,beta,omega
      real(prc) err

      itercg=0
      print *,'BICGSTAB'

c     initialise
      DR=RR
      call Mptr(RR,x2,u,DAGGER,SRF)
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
        print *,'rho',rho
        beta=rho/rhom1*(alpha/omega)
        p=r+beta*(p-omega*v)
        call Mptr(p,v,u,DAGGER,SRF)
        alpha=rho/dot_product(rhat0,v)
        h=DR+alpha*p
        call Mptr(h,best,u,DAGGER,SRF)
        errv=RR-best
        err=dot_product(errv,errv)/kferm
        print *,err
        if (err .lt. resid) then
          DR=h
          return
        end if
        s=r-alpha*v
        call Mptr(s,t,u,DAGGER,SRF)
        omega=dot_product(t,s)/dot_product(t,t)
        DR=h+omega*s
        call Mptr(DR,best,u,DAGGER,SRF)
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
      end subroutine IM5_BICGSTAB2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IM5_BICGSTAB(RR,DR,u,DAGGER,Mptr,SRF)
!     solve Mptr.DR = RR
      use numbers
      implicit none
      integer,parameter :: kferm = 4*Nv*Ls
      integer,parameter :: niterc=10*Nv*Ls
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      procedure(),pointer :: Mptr
      type(sgnratfunc) :: SRF

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
      call Mptr(DR,x1,u,DAGGER,SRF)
      r=RR-x1
      rhat0=r
      p=r
      do i=1,niterc
        itercg=itercg+1
        rhom1=dot_product(rhat0,r)
        call Mptr(p,v,u,DAGGER,SRF)
        alpha=rhom1/dot_product(rhat0,v)
        s=r-alpha*v
        err=dot_product(s,s)/kferm
        print *,err
        if (err .lt. resid) then
          DR=DR+alpha*p
          return
        end if
        call Mptr(s,t,u,DAGGER,SRF)
        omega=dot_product(t,s)/dot_product(t,t)
        DR=DR+alpha*p+omega*s
        r=s-omega*t
        err=dot_product(r,r)/kferm
        print *,err
        if (err .lt. resid) return
        rho=dot_product(rhat0,r)
        beta=rho/rhom1*(alpha/omega)
        p=r+beta*(p-omega*v)
        if (abs(rho) .lt. 1e-6) then
          rhat0=r
          p=r
        end if
      end do
      
      return
      end subroutine IM5_BICGSTAB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module axbmodule3
