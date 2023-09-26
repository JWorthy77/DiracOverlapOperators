      module diracdeterminants
      use pacc
      use options
      use numbers
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcOLDeterminant(u,DAGGER,mass,SRF)
      use ratfuncs
      use overlapmoduledev
      use shamirmodule
      implicit none
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      type(sgnratfunc) :: SRF
      procedure(),pointer :: Mptr => NULL()
      complex(prc) :: R(4*Nv)
      complex(prc) :: COL(4*Nv)
      integer i,Nv4
      integer ipiv(4*Nv),err,pow10
      complex(prc) :: DenseOL(4*Nv,4*Nv)
      complex DdD(4*Nv,4*Nv),det

      if (DWkernel.eq.1) then
        Mptr => DOLS
      elseif (DWkernel.eq.2) then
        Mptr => DOverlap
      endif

      DenseOL=czero
      Nv4=4*Nv
      do i=1,Nv4
        R=czero
        R(i)=cone
        call MPtr(R,COL,u,DAGGER,baremass,SRF)
        DenseOL(1:Nv4,i)=COL(1:Nv4)
      end do
      DdD=matmul(DenseOL,conjg(transpose(DenseOL)))
      print *,"calc determinant"
      call  cgetrf(Nv4,Nv4,DdD,Nv4,ipiv,err)
      print *,"err",err

      det=1
      pow10=0
      do i=1,Nv4
        det=det*DdD(i,i)
        print *,DdD(i,i)
        if (abs(det) > 10q0) then
          det=det/10
          pow10=pow10+1
        elseif (abs(det) < 1.0q-1) then
          det=10q0*det
          pow10=pow10-1
        endif
      end do
      print *,"det",det,"pow",pow10

      return
      end subroutine calcOLDeterminant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcKDWDeterminant(u,DAGGER,mass,det)
      use domainwallmod
      implicit none
      include 'mpif.h'
      complex(prc),intent(in) :: u(Nv,3)
      logical DAGGER
      real(prc) mass
      complex det
      complex(prc) :: R(4*Nv)
      complex(prc) :: COL(4*Nv)
      integer i,Nv4
      integer ipiv(4*Nv),err,pow10
      complex(prc) :: DenseOL(4*Nv,4*Nv)
      complex DdD(4*Nv,4*Nv)
      integer ierr,rank,nprocs
      integer ictxt,npr,npc,ipr,ipc
      integer desca(9),mxllda

c      DenseOL=czero
      Nv4=4*Nv
      do i=1,Nv4
        R=czero
        R(i)=cone
        call KDDW4(R,COL,u,DAGGER,mass)
        DenseOL(1:Nv4,i)=COL(1:Nv4)
      end do
      DdD=matmul(DenseOL,conjg(transpose(DenseOL)))
      print *,"calc K determinant"
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
      print *,rank

      npr=1
      npc=1

      call SL_INIT(ictxt,npr,npc)
      print *,ictxt
      call BLACS_GRIDINFO(ictxt,npr,npc,ipr,ipc)
      print *,ipr,ipc,npr,npc
      mxllda=Nv4
      call DESCINIT(desca,Nv4,Nv4,npr,npc,0,0,ictxt,mxllda,err)
      call  pcgetrf(Nv4,Nv4,DdD,1,1,desca,ipiv,err)
      call MPI_FINALIZE()

c      call  cgetrf(Nv4,Nv4,DdD,Nv4,ipiv,err)
      print *,"err",err

      det=1
      pow10=0
      do i=1,Nv4
        det=det*DdD(i,i)
        print *,DdD(i,i)
        if (abs(det) > 10q0) then
          det=det/10
          pow10=pow10+1
        elseif (abs(det) < 1.0q-1) then
          det=10q0*det
          pow10=pow10-1
        endif
      end do
      print *,"det K",det,"pow",pow10

      return
      end subroutine calcKDWDeterminant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcDWDeterminant(u,DAGGER,mass,det)
      use domainwallmod
      implicit none
!      include 'mpif.h'
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
      integer ierr,rank,nprocs

c      call MPI_INIT(ierr)

      DenseOL=czero
      Nv4=4*Nv*Ls
      do i=1,Nv4
        R=czero
        R(i)=cone
        call DDW(R,COL,u,DAGGER,mass)
        DenseOL(1:Nv4,i)=COL(1:Nv4)
      end do
      DdD=matmul(DenseOL,conjg(transpose(DenseOL)))
      print *,"calc DW determinant, mass:",mass
      call  cgetrf(Nv4,Nv4,DdD,Nv4,ipiv,err)
!      call  pcgetrf(Nv4,Nv4,DdD,1,1,Nv4,ipiv,err)
      print *,"err",err

c      call MPI_INIT(ierr)
c      call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
c      print *,rank

c      call SL_INIT(ictxt,npr,npc)
c      call BLACS_GRIDINFO(ictxt,npr,npc,ipr,ipc)


c      call MPI_FINALIZE()
c      return

      det=1
      pow10=0
      do i=1,Nv4
        det=det*DdD(i,i)
        print *,DdD(i,i)
        if (abs(det) > 10q0) then
          det=det/10
          pow10=pow10+1
        elseif (abs(det) < 1.0q-1) then
          det=10q0*det
          pow10=pow10-1
        endif
      end do
      print *,"det DW",det,"pow",pow10

      return
      end subroutine calcDWDeterminant
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module diracdeterminants
