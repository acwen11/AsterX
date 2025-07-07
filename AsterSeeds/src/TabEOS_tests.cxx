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

  } else {
    CCTK_ERROR("Test case not defined");
  }
}

} // namespace AsterSeeds
