#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "setup_eos.hxx"
#include "seeds_utils.hxx"

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;

extern "C" void Atmosphere_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Atmosphere_Initialize;
  DECLARE_CCTK_PARAMETERS;

  // For all the tests, the initial data EOS is ideal gas
  // Get local eos object
  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  const CCTK_REAL dummy_ye = 0.5;
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        rho(p.I) = rho_atmosphere;
        velx(p.I) = 0.;
        vely(p.I) = 0.;
        velz(p.I) = 0.;
        press(p.I) = press_atmosphere;
        eps(p.I) = eos_3p_ig->eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                          dummy_ye);
      });

  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.; });

  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.; });

  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.; });
}

} // namespace AsterSeeds
