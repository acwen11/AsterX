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

extern "C" void Tests2D_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Tests2D_Initialize;
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

  // See Cipolletta et al (2020) and Del Zanna, Bucciantini, Londrillo (2003)
  if (CCTK_EQUALS(test_case, "magnetic rotor")) {
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL r_cyl = sqrt(pow2(p.x) + pow2(p.y));

          if (r_cyl < 0.1) {
            const CCTK_REAL omega = 9.95;
            rho(p.I) = 10.;
            velx(p.I) = -omega * p.y;
            vely(p.I) = omega * p.x;
            velz(p.I) = 0.;
          }

          else {
            rho(p.I) = 1.;
            velx(p.I) = 0.;
            vely(p.I) = 0.;
            velz(p.I) = 0.;
          }

          press(p.I) = 1.;
          eps(p.I) = eos_3p_ig->eps_from_valid_rho_press_ye(
              rho(p.I), press(p.I), dummy_ye);
        });

    grid.loop_all<1, 0, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_x(p.I) = 0.;
                                                 });

    grid.loop_all<0, 1, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_y(p.I) = 0.;
                                                 });

    grid.loop_all<0, 0, 1>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_z(p.I) = p.y;
                                                 });
  }

  // See Cipolletta et al (2020), and Beckwith, Stone (2011)
  else if (CCTK_EQUALS(test_case, "magnetic loop advection")) {
    CCTK_REAL axial_vel;

    if (CCTK_EQUALS(mag_loop_adv_type, "2D")) {
      if (CCTK_EQUALS(mag_loop_adv_axial_vel, "zero"))
        axial_vel = 0.;
      else if (CCTK_EQUALS(mag_loop_adv_axial_vel, "non-zero"))
        axial_vel = 1. / 24.;
      else {
        CCTK_VERROR("Invalid loop advection case");
      }
    }

    // TODO: implement this (see Cipolletta et al (2020) and code)
    else if (CCTK_EQUALS(mag_loop_adv_type, "3D")) {
      CCTK_VERROR("Sorry, 3D loop advection hasn't been implemented yet");
    }

    else {
      CCTK_VERROR("Invalid loop advection type");
    }

    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.;
          press(p.I) = 3.;
          eps(p.I) = eos_3p_ig->eps_from_valid_rho_press_ye(
              rho(p.I), press(p.I), dummy_ye);
          velx(p.I) = 1. / 12.0;
          vely(p.I) = 1. / 24.;
          velz(p.I) = axial_vel;
        });

    grid.loop_all<1, 0, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_x(p.I) = 0.;
                                                 });

    grid.loop_all<0, 1, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_y(p.I) = 0.;
                                                 });

    grid.loop_all<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL r_cyl = sqrt(pow2(p.x) + pow2(p.y));
          const CCTK_REAL R_loop = 0.3;
          const CCTK_REAL A_loop = 0.001;

          if (r_cyl < R_loop)
            Avec_z(p.I) = A_loop * (R_loop - r_cyl);
          else
            Avec_z(p.I) = 0.;
        });
  }

  // See Cipolletta et al (2020)
  else if (CCTK_EQUALS(test_case, "cylindrical blast")) {
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL r_cyl = sqrt(pow2(p.x) + pow2(p.y));
          const CCTK_REAL r_in = 0.8;
          const CCTK_REAL r_out = 1.0;
          const CCTK_REAL rho_in = 1e-2;
          const CCTK_REAL rho_out = 1e-4;
          const CCTK_REAL p_in = 1.0;
          const CCTK_REAL p_out = 3.0e-5;

          if (r_cyl < r_in) {
            rho(p.I) = rho_in;
            press(p.I) = p_in;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
          }

          else if (r_cyl > r_out) {
            rho(p.I) = rho_out;
            press(p.I) = p_out;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
          } else {
            rho(p.I) = exp(((r_out - r_cyl) * log(rho_in) +
                            (r_cyl - r_in) * log(rho_out)) /
                           (r_out - r_in));
            press(p.I) = exp(
                ((r_out - r_cyl) * log(p_in) + (r_cyl - r_in) * log(p_out)) /
                (r_out - r_in));
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
          }
          eps(p.I) = eos_3p_ig->eps_from_valid_rho_press_ye(
              rho(p.I), press(p.I), dummy_ye);
        });

    grid.loop_all<1, 0, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_x(p.I) = 0.0;
                                                 });

    grid.loop_all<0, 1, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_y(p.I) = 0.0;
                                                 });

    grid.loop_all<0, 0, 1>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_z(p.I) = 0.1 * p.y;
                                                 });
  }

  else if (CCTK_EQUALS(test_case, "KHI")) {
    printf("Setting KHI data: \n");
    printf(" rhoUp = %16.8e \n", rhoUp);
    printf(" rhoLow = %16.8e \n", rhoLow);

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // set velocities in opposite directions in x along the slip surface
          // set slightly different densities along the slip surface
          // location of slip surfaces hardcoded at (abs(p.y) = 0.25)
          using std::abs;
          if (abs(p.y) >= 0.25) {
            rho(p.I) = rhoUp;
            velx(p.I) = vxUp;
          } else {
            rho(p.I) = rhoLow;
            velx(p.I) = vxLow;
          }
        //  printf("  rho, p.y = %16.8e, %16.8e \n", rho(p.I), p.y); 
          // excite the instability by peturbing v^y
          using std::exp, std::pow, std::sin;
          vely(p.I) = w0 * sin(4 * M_PI * p.x) *
                      (exp(-pow(p.y - 0.25, 2) / (2 * pow(sigma, 2))) +
                       exp(-pow(p.y + 0.25, 2) / (2 * pow(sigma, 2))));

          velz(p.I) = 0.0;

          // set constant initial pressure throughout the domain
          press(p.I) = p_val;

          // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          eps(p.I) = eos_3p_ig->eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
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
  }

  else if (CCTK_EQUALS(test_case, "magKHI")) {

    const CCTK_REAL alpha = 0.01;
    const CCTK_REAL beta = 0.1;
    const CCTK_REAL v0 = 0.25;
    const CCTK_REAL epsilon = 0.01*v0;
    const CCTK_REAL pr = 20.0;

    const CCTK_REAL sigma_pol = 0.01;
    const CCTK_REAL sigma_tor = 1;
    const CCTK_REAL Bx = sqrt(2.0 * pr * sigma_pol);
    const CCTK_REAL Bz = sqrt(2.0 * pr * sigma_tor);

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          const CCTK_REAL kx = 2*M_PI*p.x;
          const CCTK_REAL arg = (p.y*p.y)/(beta*beta);

          rho(p.I) = 1.0;
          press(p.I) = pr;
          velx(p.I) = -v0 * tanh(p.y/alpha);
          // excite the instability by peturbing v^y
          vely(p.I) = 0.5 * epsilon * (sin(kx) - sin(-kx) )*exp(-arg);
          velz(p.I) = 0.0;

          eps(p.I) = eos_3p_ig->eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = -Bz * p.y; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = Bz * p.x - Bx * p.y; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = Bx * p.y;; });
  }

  else {
    CCTK_ERROR("Test case not defined");
  }
}

} // namespace AsterSeeds
