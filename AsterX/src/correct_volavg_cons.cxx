#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include "aster_utils.hxx"

namespace AsterX {
using namespace Loop;
using namespace AsterUtils;

extern "C" void AsterX_CorrectVolAvgCons(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CorrectVolAvgCons;
  DECLARE_CCTK_PARAMETERS;

  // Loop over interior
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Reduce order in presence of shock
        const CCTK_INT corr_ordL = LOflag(p.I) ? correction_order : 2;
        dens(p.I) = correct_conservs(dens, p, corr_ordL);
        momx(p.I) = correct_conservs(momx, p, corr_ordL);
        momy(p.I) = correct_conservs(momy, p, corr_ordL);
        momz(p.I) = correct_conservs(momz, p, corr_ordL);
        tau(p.I) = correct_conservs(tau, p, corr_ordL);
        DYe(p.I) = correct_conservs(DYe, p, corr_ordL);
        DEnt(p.I) = correct_conservs(DEnt, p, corr_ordL);
      });
}

} // namespace AsterX
