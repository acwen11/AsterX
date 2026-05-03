#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <seeds_utils.hxx>

#include <setup_eos.hxx>

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace EOSX;

extern "C" void TabEOSTests_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_TabEOSTests_Initialize;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(test_case, "Isotropic Gas")) {

    grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                  [=] CCTK_DEVICE(const PointDesc &p)
                                      CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                        rho(p.I) = 1e-4;
                                        velx(p.I) = 0.0;
                                        vely(p.I) = 0.0;
                                        velz(p.I) = 0.0;
                                      });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

  } else if (CCTK_EQUALS(test_case, "TabEOS shock")) {

    const CCTK_REAL rho_in    = 1e-4;
    const CCTK_REAL rho_out   = 1e-6;

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
          } else if (r2 >= rout2){
      // outer region
            rho(p.I) = rho_out;
          } else {
     // transition region
            const CCTK_REAL r = sqrt(r2);
            const CCTK_REAL denom = (rout - rin);

            const CCTK_REAL w_in  = (denom > 0.0) ? (rout - r) / denom : 1.0;
            const CCTK_REAL w_out = (denom > 0.0) ? (r - rin) / denom : 0.0;

            const CCTK_REAL log_rho   = w_in * log(rho_in)   + w_out * log(rho_out);

            rho(p.I) = exp(log_rho);
          }

          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
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
