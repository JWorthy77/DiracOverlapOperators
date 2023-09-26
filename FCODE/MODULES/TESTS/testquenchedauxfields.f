      module testquenchedaux
      use pacc
      use arraysizes
      use numbers
      use options
      use rvmodule
      use utilsmod
      use gaugemodule
      use gaugefield,only: fieldtheta=>theta
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine testquenchedauxfields
      implicit none
      integer,parameter :: Ncd=100
      integer,parameter :: Nbin=Ncd-1
      integer N,i,j,idb,mu
      real(prc) dpi,theta,beta1,beta2,beta3,cp,x,dx,xL,xR,cpL,cpR
      complex(prc) U
      real(prc),dimension(Ncd) :: cdist,tdist
      real(prc) th
      integer,dimension(Nbin) :: bcount

!     check cumulative distribution
      
      beta1=1*one
      call makecumcosdist(beta1,Ncd,cdist)
      call make2piedgedist(Ncd,tdist)
      bcount=0
      do j=1,Nv*3
        th=crv(Ncd,cdist)
        idb=intbin(th,Ncd,tdist)
!        print *,th,idb
        bcount(idb)=bcount(idb)+1
      end do
      print *,"bcount:",bcount
      gbeta=beta1
      call makeQuenchedCosineThirringField()
      bcount=0
      do j=1,Nv
        do mu=1,3
          idb=intbin(fieldtheta(j,mu),Ncd,tdist)
          bcount(idb)=bcount(idb)+1
        end do
      end do
      print *,"bcount:",bcount
          

      stop

      N=100
      dpi=2*pi/N
      beta1=one*2
      beta2=one
      beta3=one*0.001
!     look at pdf
      open(unit=41,file='qthetapdfB.dat',status='unknown')
      open(unit=42,file='qUpdfB.dat',status='unknown')
      do i=0,N
        theta=-pi+i*dpi
        cp=cospdf(beta3,theta)
        U=exp(zi*theta)
        print *,theta,cp
        write(41,*) theta,cp
        write(42,*) theta,real(U),dimag(U)
      end do
      close(41)
      close(42)

      N=100
      dx=pi/N
      open(unit=43,file='qthetabadpdfB.dat',status='unknown')
      open(unit=44,file='qUbadpdfB.dat',status='unknown')
      do i=0,N
        cpL=zero
        cpR=zero
        do j=5,1,-1
          xL=-j*pi+i*dx
          cpL=cpL+gausspdf(beta1,xL)
          xR=j*pi-i*dx
          cpR=cpR+gausspdf(beta1,xR)
        end do

        U=exp(zi*x)
        print *,xL,cpL
        write(43,*) xL,cpL
        write(43,*) xR,cpR
        write(44,*) x,real(U),dimag(U)
      end do
      close(43)
      close(44)

      return
      end subroutine testquenchedauxfields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module testquenchedaux
