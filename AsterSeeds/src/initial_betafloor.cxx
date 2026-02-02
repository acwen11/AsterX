#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "seeds_utils.hxx"
#include "setup_eos.hxx"

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace AsterUtils;
using namespace EOSX;

extern "C" void AsterSeeds_SetInitialBetaFloor(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_SetInitialBetaFloor;
  DECLARE_CCTK_PARAMETERS;

  auto eos_3p_tab3d = global_eos_3p_tab3d;
  if (not CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"Tabulated3d\" in your parameter file.",
                evolution_eos);
  }

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        const CCTK_REAL pressL = press(p.I);
        const CCTK_REAL rhoL = rho(p.I);
        const CCTK_REAL tempL = temperature(p.I);
        const CCTK_REAL YeL = Ye(p.I);
        const CCTK_REAL epsL = eps(p.I);
        const CCTK_REAL entL = entropy(p.I);

        const CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        // Grading rho
        // CCTK_REAL rho_atm = (radial_distance > r_atmo)
        //     ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
        //     : rho_abs_min;
        // const CCTK_REAL rho_atmo_cut = rho_atm * (1 + atmo_tol);

        // Compute b^2
        /* Get covariant metric */
        const smat<CCTK_REAL, 3> glo(
            [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

        vec<CCTK_REAL, 3> B_up{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
        vec<CCTK_REAL, 3> B_low = calc_contraction(glo, B_up);

        vec<CCTK_REAL, 3> v_up{velx(p.I), vely(p.I), velz(p.I)};
        vec<CCTK_REAL, 3> v_low = calc_contraction(glo, v_up);

        const CCTK_REAL wlor = calc_wlorentz(v_up, v_low);
        const CCTK_REAL alp_b0 = wlor * calc_contraction(B_up, v_low);

        const CCTK_REAL B2 = calc_contraction(B_up, B_low);
        const CCTK_REAL bsq = ( B2 + alp_b0 * alp_b0 ) / ( wlor*wlor );

        // Increase P if necessary
        CCTK_REAL press_lim = initial_beta_min * bsq * 0.5;

        if ((pressL >= press_lim) || (rhoL > initial_beta_rhocut)) {
          press(p.I) = pressL;
          
          rho(p.I) = rhoL;
          eps(p.I) = epsL;
          entropy(p.I) = entL;
        }
        else {
          // Recalculate primitives
          press(p.I) = press_lim;
          rho(p.I) = eos_3p_tab3d->rho_from_valid_press_temp_ye(press_lim, tempL, YeL);
          eps(p.I) = eos_3p_tab3d->eps_from_valid_rho_temp_ye(rho(p.I), tempL, YeL);
          entropy(p.I) = eos_3p_tab3d->entropy_from_valid_rho_temp_ye(rho(p.I), tempL, YeL);
        }
        
      });

  if (unmagnetized_test)
  {
    // Reset everything to 0 again
    grid.loop_int<1, 0, 0>(grid.nghostzones,
                           [=] CCTK_HOST(const Loop::PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 Avec_x(p.I) = 0.0;
                               });

    grid.loop_int<0, 1, 0>(grid.nghostzones,
                           [=] CCTK_HOST(const Loop::PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 Avec_y(p.I) = 0.0;
                               });

    grid.loop_int<0, 0, 1>(grid.nghostzones,
                           [=] CCTK_HOST(const Loop::PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 Avec_z(p.I) = 0.0;
                               });
  }

}

} // namespace AsterSeeds
