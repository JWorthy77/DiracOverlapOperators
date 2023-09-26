      module diracdeterminants
      use pacc
      use options
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcDWDeterminant(u,DAGGER,mass,det)
      use domainwallmod
      implicit none
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex det
      complex(prc) :: R(4*Nv*Ls)
      complex(prc) :: COL(4*Nv*Ls)
      integer i,Nv4
      integer ipiv(4*Nv*Ls),err,pow10
      complex(prc) :: DenseOL(4*Nv*Ls,4*Nv*Ls)
      complex DdD(4*Nv*Ls,4*Nv*Ls)
c      integer ierr

c      call MPI_INIT(ierr)
c      call MPI_FINALIZE()
c      return

      DenseOL=czero
      Nv4=4*Nv*Ls
      do i=1,Nv4
        R=czero
        R(i)=cone
        call DDW(R,COL,u,DAGGER,mass)
        DenseOL(1:Nv4,i)=COL(1:Nv4)
      end do
      DdD=matmul(DenseOL,conjg(transpose(DenseOL)))
      DdD=matmul(DenseOL,DenseOL)
      print *,"calc DW determinant, mass:",mass
!      call  cgetrf(Nv4,Nv4,DdD,Nv4,ipiv,err)
!      call  pcgetrf(Nv4,Nv4,DdD,1,1,Nv4,ipiv,err)
      print *,"err",err

      return
      end subroutine calcDWDeterminant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module diracdeterminants
