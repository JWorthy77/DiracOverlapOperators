!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! must be included in order of dependence
#include "basicmodules.f"
#include "UTILS/statsmodule.f"
#include "UTILS/rvmodule.f"
#include "UTILS/pfmodule.f"
#include "UTILS/ellipticmodule.f"
#include "UTILS/zolomodule.f"
#include "UTILS/polyapproxmod.f"
#include "UTILS/extpolyapproxmod.f"
#include "UTILS/remezapproxmod.f"
#include "UTILS/ratfuncs.f"
#include "UTILS/utils.f"
#include "DIRAC_OPS/WilsonDirac.f"
#include "DIRAC_OPS/domainwall.f"
#include "DIRAC_OPS/overlap.f"
#include "AUXFIELD/gaugefield.f"
#include "AUXFIELD/hmc2wilsonferms.f"
!#include "AUXFIELD/hmc2domwallferms.f"
!#include "gaugefieldSU3.f"
#include "AUXFIELD/gaugemodule.f"
!#include "diracmontecarlo.f"
#include "UTILS/IOmodule.f"
!#include "CONDENSATE/condpointtools.f"
!#include "CONDENSATE/condnoisytools.f"
!#include "CONDENSATE/condensatetools.f"
!#include "CONDENSATE/condensate.f"
!#include "localitymodule.f"
!#include "gaugesmoothness.f"
!#include "kernelspectra.f"
!#include "overlapspectra.f"
!#include "spectra.f"
!#include "gwmodule.f"
!#include "TESTS/testratfuncs.f"
!#include "TESTS/testbasicdiracmod.f"
!#include "TESTS/testVOL.f"
!#include "TESTS/testoverlap.f"
!#include "TESTS/testdomwall.f"
!#include "TESTS/testgamherm.f"
!#include "TESTS/testopconv.f"
!#include "TESTS/testopequiv.f"
!#include "TESTS/testcondensates.f"
!#include "TESTS/validatecondensates.f"
!#include "TESTS/testquenchedauxfields.f"
!#include "determinants.f"
!#include "parallelmodule.f"
#include "allhmcmodules.f"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
