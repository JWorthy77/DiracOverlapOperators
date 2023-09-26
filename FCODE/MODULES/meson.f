      module mesonmod
      use pacc
      use arraysizes
      use numbers
      use arpackmodule
      implicit none
      contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getThirringMesons()
      use numbers
      use options
      use IOmodule
      use gaugefield
      use statsmod
      implicit none
      integer,parameter :: Ngf=50
      real(prc) Mmin,Mmax,dM
      real(prc) lmin(Ngf),lmax(Ngf),lminav,lmaxav,lminsd,lmaxsd
      integer i,idx     

c      DWkernel=1 ! Shamir
!      open(unit=11,file='Mesons.dat',status='unknown',
!     &                                             form='formatted')
      idx=0
      do i=5,5*Ngf,5
        call readThirringConFile(i,theta)
        call coef(u,theta)
        idx=idx+1
        call calcDEigs(3,lmin(idx),lmax(idx))
!        write(11,*) lmin(idx),lmax(idx),lmax(idx)/lmin(idx)
        print *, lmin(idx),lmax(idx),lmax(idx)/lmin(idx)
      end do
      call calcVarReal(Ngf,lmin,lminav,lminsd)
      call calcVarReal(Ngf,lmax,lmaxav,lmaxsd)
      write(11,*) lminav,lminsd,lmaxav,lmaxsd
      print *,"lmin: ",lminav,lminsd
      print *,"lmax: ",lmaxav,lmaxsd
      close(11)

      return
      end subroutine getThirringMesons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine calcTimeSlice(ts)
!     calculate time slice at time ts
      use gaugefield
      use wilsonmodule
      use shamirmodule
      use rvmodule
      use options
      use axbmodule1
      use domainwallmod
      implicit none
      integer ts
      integer Nmax
      real(prc) lmax,lmin
      complex(prc),dimension(Nv,Ndc) :: R,DR,TMP
      complex(prc),dimension(Nv,Ndc,Ls) :: RL,DRL,TMPL
      real ev(1)
      integer i
      procedure(),pointer :: Mptr => NULL()
      logical DOMWALL


      
      call IDWilson(R,DR,u,.false.,-MDW,czero)
      
      do ix=1,Ns
        do iy=1,Ns
          idx=(ts-1)*Ns*Ns+(iy-1)*Ns+ix
          sval=sval+DR(idx)
        end do
      end do

      return
      end subroutine calcTimeSlice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end module mesonmod
