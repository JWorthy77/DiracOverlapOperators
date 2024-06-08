!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! must be included in order of dependence
#include "HMC_MODULES/modulesRHMC.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      program hmc
      use options
      use gaugemodule
      implicit none
      integer j,Naux

      print *,"Start program Dirac HMC"
      call init()

      Naux=50
      do j=1,Naux
        print *,"make field ",j," of ",Naux
        call makeGaugeField(.false.)
!        call measureCondensateInstance()
      end do
!      call  averageCondensate();

      call finalise()
      end program hmc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine init()
      use timer
      use countmod
      use rvmodule
      use options
      use gammas
      use paulimodule
      use indices
      use gaugefield
      implicit none

      call initRVs(.false.,.false.,0) ! uses a time based seed initialiser

      oc_idx=0
      outer_count=0
      ic_idx=0
      inner_count=0

      call init_timer()

      GAUGETYPE=2 ! 1=compact 2=non-compact
      MDW=one ! domain wall height
      MTYPE=3 ! mass term type
      DWkernel=2 ! 1 for Shamir, 2 for Wilson
      OLTYPE=1 ! 1 is direct calculation, 2 is for DW (K)-type
      baremass=0.01d0
      gbeta=1.0

      call setGammas
      call setPauliMatrices
#ifdef TWODIMENSIONS
      call setIndices2d
#else
      call setIndices
#endif

      theta=0
      call coef(u,theta)

      return
      end subroutine init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine finalise()
      implicit none

       return
       end subroutine finalise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

