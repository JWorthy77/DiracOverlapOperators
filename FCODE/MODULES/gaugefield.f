      module gaugefield
      use pacc
      use arraysizes
      use numbers
      implicit none
      real(prc) theta(Nv,3) ! turn these into a derived type
      complex(prc) u(Nv,3)

      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine coef(u,theta)
c      use pacc
c      use arraysizes
      use indices
      use options
c      use numbers
      implicit none
      complex(prc) u(Nv,3)
      real(prc) theta(Nv,3)
      integer i,j,t,mu,ioffset,ind
      real(prc) a,b

      if (GAUGETYPE.eq.1) then ! QED
        do mu=1,3
          do i=1,Nv
            u(i,mu)=exp(cmplx(zero,theta(i,mu)))
          enddo
        enddo
      else if (GAUGETYPE.eq.2) then ! Thirring
        do mu=1,3
          do i=1,Nv
            u(i,mu)=cone+zi*theta(i,mu)
          enddo
        enddo
      else if (GAUGETYPE.eq.3) then ! GN
        print *,"hmm, in coef (module gaugefield)"
        stop
      elseif (GAUGETYPE.eq.4) then ! QED
        do mu=1,3
          do i=1,Nv
            a=1d0
            b=one-a
            u(i,mu)=b*exp(cmplx(zero,theta(i,mu))) 
     &                                   +a*(cone+zi*theta(i,mu))
          enddo
        enddo
      end if
!     anti-p.b.c. in timelike direction
      ioffset=(Nt-1)*Ns*Ns
      do i=1,Ns*Ns
        ind=ioffset+i
        u(ind,3)=-u(ind,3)
      enddo
      return

!     anti-p.b.c. in spacelike directions
c      j=Ns
c      do i=1,Ns
c        do t=1,Nt
c          call ia(i,j,t,ind)
c          u(ind,2)=-u(ind,2)
c        end do
c      end do
!     anti-p.b.c. in spacelike directions
c      i=Ns
c      do j=1,Ns
c        do t=1,Nt
c          call ia(i,j,t,ind)
c          u(ind,1)=-u(ind,1)
c        end do
c      end do
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine writeGaugeField(fname)
      implicit none
      character*1024 fname

      open(unit=10,file=fname,status='unknown')
      write(10,*) Nv
      write(10,*) theta
      write(10,*) u
      close(10)
      return
      end subroutine writeGaugeField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine writeGaugeSlice(fname,th3)
      implicit none
      character*1024 fname
      real(prc) th3(Ns,Ns,Nt,3)
      integer i,t

      open(unit=10,file=fname,status='unknown')
c      write(10,*) Ns,Nt
      do i=1,Ns
        do t=1,Nt
          write(10,*) i,t,th3(Ns/2,i,t,3)
        end do
      end do
      close(10)
      return
      end subroutine writeGaugeSlice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine writeGaugeLines(fname,th3)
      implicit none
      character*1024 fname
      real(prc) th3(Ns,Ns,Nt,3)
      integer i,t

      open(unit=10,file=fname,status='unknown')
c      write(10,*) Ns,Nt
      do i=1,Ns
        write(10,*) i,1,1,3,th3(i,1,1,3),i,Ns/2,1,3,th3(i,Ns/2,1,3),
     &              i,Ns/2,1,3,th3(i,1,Ns/2,3)
      end do
      close(10)
      return
      end subroutine writeGaugeLines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine gaugeTransformU(u,alpha)
!     local gauge transform of U, by alpha
      use pacc
      use arraysizes
      use indices
      use numbers
      implicit none
      complex(prc) u(Nv,3)
      real(prc) alpha(Nv)
      integer i,mu

      do mu=1,3
        do i=1,Nv
          u(i,mu)=exp(cmplx(zero,alpha(i)))*u(i,mu)*
     &                 exp(-cmplx(zero,alpha(iu(i,mu))))
        enddo
      enddo
      return
      end subroutine gaugeTransformU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module gaugefield
