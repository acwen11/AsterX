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

extern "C" void Tests3D_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Tests3D_Initialize;
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

  if (CCTK_EQUALS(test_case, "spherical shock")) {

    const CCTK_REAL rho_in    = 1e-2;
    const CCTK_REAL rho_out   = 1e-4;
    const CCTK_REAL press_in  = 1.0;
    const CCTK_REAL press_out = 3e-5;

    const CCTK_REAL rin  = shock_radius * (1.0 - 0.2); // 0.8 * shock_radius
    const CCTK_REAL rout = shock_radius;

    const CCTK_REAL rin2  = pow2(rin);
    const CCTK_REAL rout2 = pow2(rout);

    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

      	  const CCTK_REAL r2 = pow2(p.x) + pow2(p.y) + pow2(p.z);
          if (r2 <= rin2) {
	    // inner region
            rho(p.I) = rho_in;
            press(p.I) = press_in;
          } else if (r2 >= rout2){
	    // outer region
            rho(p.I) = rho_out;
            press(p.I) = press_out;
          } else {
	   // transition region
            const CCTK_REAL r = sqrt(r2);
            const CCTK_REAL denom = (rout - rin);

            const CCTK_REAL w_in  = (denom > 0.0) ? (rout - r) / denom : 1.0;
            const CCTK_REAL w_out = (denom > 0.0) ? (r - rin) / denom : 0.0;

            const CCTK_REAL log_rho   = w_in * log(rho_in)   + w_out * log(rho_out);
            const CCTK_REAL log_press = w_in * log(press_in) + w_out * log(press_out);

            rho(p.I) = exp(log_rho);
            press(p.I) = exp(log_press);
          }

	  velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          eps(p.I) = eos_3p_ig->eps_from_valid_rho_press_ye(
              rho(p.I), press(p.I), dummy_ye);
        });

    grid.loop_all<1, 0, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_x(p.I) = 0.0;
                                                 });

    grid.loop_all<0, 1, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_y(p.I) = amplitude*p.x;
                                                 });

    grid.loop_all<0, 0, 1>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_z(p.I) = 0.0;
                                                 });

  } else {
    CCTK_ERROR("Test case not defined");
  }
}

} // namespace AsterSeeds
