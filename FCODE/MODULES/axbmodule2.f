      module axbmodule2
      use pacc
      use arraysizes
      use numbers
      use ratfuncs
      implicit none
      type ioptions
        type(sgnratfunc) :: SRF
        logical :: DAGGER
        real(prc) :: mass
        complex(prc) :: add
      end type
      logical,parameter :: VB_AXB2=.false.
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMnonsymOpts(RR,DR,u,DAGGER,Mptr,iopts) ! Mptr like SGNfactor
      implicit none
      complex(prc) RR(Nv,4),DR(Nv,4)
      complex(prc) u(Nv,3)
      logical DAGGER
      procedure(),pointer :: Mptr
      type(ioptions) iopts
      complex(prc) TMP(Nv,4)

      if (VB_AXB2) then ; print *,"IMnonsymOpts" ; end if ;
      if (.not. DAGGER) then
        call Mptr(RR,TMP,u,.not.DAGGER,iopts%mass,iopts%SRF)
        call IMdagMOpts(TMP,DR,u,DAGGER,Mptr,iopts)
      elseif (DAGGER) then
        call IMdagMOpts(RR,TMP,u,DAGGER,Mptr,iopts)
        call Mptr(TMP,DR,u,.not.DAGGER,iopts%mass,iopts%SRF)
      end if

      return
      end subroutine IMnonsymOpts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMdagMOpts(RR,DR,u,DAGGER,Mptr,iopts)
      use countmod
      implicit none
      integer,parameter :: kferm = 4*Nv
      integer,parameter :: niterc=Nv*10
      complex(prc) RR(kferm), DR(kferm)
      complex(prc) u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex(prc) add
      procedure(),pointer :: Mptr
      type(ioptions) iopts

      complex(prc) x1(kferm),x2(kferm),p(kferm),r(kferm)
      integer itercg
      integer nx,i
      real(prc) beta,betan,betad,alpha,alphad,alphan

      itercg=0

      if (VB_AXB2) then ; print *,"IMdagMOpts" ; end if ;
c     initialise
      DR=RR
      call Mptr(RR,x1,u,DAGGER,iopts%mass,iopts%SRF)
      call Mptr(x1,x2,u,.not.DAGGER,iopts%mass,iopts%SRF)
      r=RR-x2
      p=r
      betan=sum(conjg(r)*r)
      if (betan.lt.resid) goto 8
      
      do i=1,niterc
        itercg=itercg+1
        call Mptr(p,x1,u,DAGGER,iopts%mass,iopts%SRF)
        call Mptr(x1,x2,u,.not.DAGGER,iopts%mass,iopts%SRF)
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
      oc_idx=oc_idx+1
      outer_count=outer_count+i
!      print *,itercg,niterc,betan
      return
      end subroutine IMdagMOpts     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine IMOpts(RR,DR,u,DAGGER,mass,add,Mptr)
!     solve Mptr.DR = RR
!     requires M to be symmetric, positive definite
      use pacc
      use arraysizes
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
      print *,itercg,niterc,betan
      return
      end subroutine IMOpts
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module axbmodule2
