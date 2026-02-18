#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "loop_device.hxx"
#include "aster_utils.hxx"

namespace AsterX {
using namespace Loop;
using namespace AsterUtils;

extern "C" void AsterX_Add_Dissipation_AvecPsi(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Add_Dissipation_AvecPsi;
  DECLARE_CCTK_PARAMETERS;

  // KO-5 requires >=3 ghost zones
  const int need = 3;
  if (cctk_nghostzones[0] < need || cctk_nghostzones[1] < need ||
      cctk_nghostzones[2] < need) {
    CCTK_VERROR("KO-5 requires >=3 ghost zones; have (%d,%d,%d).",
                cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]);
  }

  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_x_rhs(p.I) += eps_diss * KO_dissipation(p, Avec_x);
      });

  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_y_rhs(p.I) += eps_diss * KO_dissipation(p, Avec_y);
      });

  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Avec_z_rhs(p.I) += eps_diss * KO_dissipation(p, Avec_z);
      });

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Psi_rhs(p.I) += eps_diss * KO_dissipation(p, Psi);
      });
}

} // namespace AsterX
