!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module allmodules
      use pacc ! set precision
      use arraysizes ! array sizes
      use numbers ! parameterise one, cone, zeta, ...
#ifdef PARALLEL
      use basicparallelmod ! nprocs/rank
#endif
      use options
      use countmod
      use paulimodule ! pauli arrays
      use gammas ! gamma arrays
      use indices ! indices for lattice
      use pfmodule 
      use statsmod
      use rvmodule
      use ellipticmodule
      use zolomodule
      use polyapproxmod
      use extpolyapproxmod
      use remezapproxmod
      use ratfuncs
      use utilsmod
      use gaugefield
      use gaugefieldSU3
      use axbmodule1
      use axbmodule2
      use axbmodule3
      use axbmoduleReal
      use basicdiracopsmod
      use staggeredmodule
      use domainwallmod
      use arpackmodule
      use overlapmoduledev
      use overlapmoduleprod
      use overlapmodulescaled
      use weylmod
      use olweylmod
!      use axbmodule5
      use diracmcmodule
      use gaugemodule
      use scalargaugemodule
      use IOmodule
      use condpointtoolsmod
      use condnoisytoolsmod
      use condensatetoolsmod
      use condensatemodule
      use smearedcondensate
      use mesonmodule
      use mesonmassmodule
      use localitymodule
      use gaugesmoothness
      use multigridmod
      use kernelspectrarange
      use overlapspectrarange
      use spectra
      use gwmodule
      use testratfuncsmod
      use testbasicdiracopsmod
      use testVOLmod
      use testoverlapopsmod
      use testdomwallmod
      use testgamhermmod
      use testDequivmod
      use testDconvmod
      use testcondmod
      use validatecondmodule
      use testquenchedaux
!      use diracdeterminants ! calculates determinants of DdaggerD for different D
      use parallelmod ! sequential routines
      end module allmodules
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
