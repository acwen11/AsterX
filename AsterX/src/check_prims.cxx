#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include "eos.hxx"
#include "eos_idealgas.hxx"
#include "aster_utils.hxx"

namespace AsterX {
using namespace AsterUtils;
using namespace Loop;
using namespace EOSX;
using namespace std;

extern "C" void AsterX_CheckPrims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CheckPrims;
  DECLARE_CCTK_PARAMETERS;

  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
       rgye(ye_min, ye_max);
  
  const eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);

  // Loop over the entire grid (0 to n-1 cells in each direction)
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Interpolate metric terms from vertices to center
        const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                   calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                   calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};

        vec<CCTK_REAL, 3> v_up{velx(p.I), vely(p.I), velz(p.I)};

        CCTK_REAL rhoL     = rho(p.I);
        CCTK_REAL epsL     = eps(p.I);
        CCTK_REAL pressL   = press(p.I);
        CCTK_REAL entropyL = entropy(p.I);

        // Lower velocity
        vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);

        // Check validity of primitives ----------
        // Mostly copy of code from c2p::prims_floors_and_ceilings

        const CCTK_REAL w_lim = sqrt(1.0 + vw_lim * vw_lim);
        const CCTK_REAL v_lim = vw_lim / w_lim;
      
        // ----------
        // Floor and ceiling for rho and velocity
        // Keeps pressure the same and changes eps
        // ----------
      
        // check if computed velocities are within the specified limit
        CCTK_REAL vsq_Sol = calc_contraction(v_low, v_up);
        CCTK_REAL sol_v = sqrt(vsq_Sol);
        if (sol_v > v_lim) {

          v_up *= v_lim / sol_v;
          v_low = v_low * v_lim / sol_v;
     
        }
      
        if (rhoL > eos_th.rgrho.max) {
      
          // remove mass
          rhoL = eos_th.rgrho.max;
          epsL = eos_th.eps_from_valid_rho_press_ye(eos_th.rgrho.max, pressL, 0.5);
          entropyL =
              eos_th.kappa_from_valid_rho_eps_ye(eos_th.rgrho.max, epsL, 0.5);

        }

        if (rhoL < rho_abs_min * (1 + atmo_tol)) {
      
          // add mass
          rhoL = rho_abs_min;
          epsL = eos_th.eps_from_valid_rho_press_ye(rho_abs_min, pressL, 0.5);
          entropyL =
              eos_th.kappa_from_valid_rho_eps_ye(rho_abs_min, epsL, 0.5);

        }
      
        // ----------
        // Floor and ceiling for eps
        // Keeps rho the same and changes press
        // ----------
      
        // check the validity of the computed eps
        auto rgeps = eos_th.range_eps_from_valid_rho_ye(rhoL, 0.5);
        if (epsL > rgeps.max) {

          epsL   = rgeps.max;
          pressL = eos_th.press_from_valid_rho_eps_ye(rhoL, epsL, 0.5);
          entropyL = eos_th.kappa_from_valid_rho_eps_ye(rhoL, epsL, 0.5);  

        } else if (epsL < rgeps.min) {

          epsL = rgeps.min;
          pressL = eos_th.press_from_valid_rho_eps_ye(rhoL, epsL, 0.5);
          entropyL = eos_th.kappa_from_valid_rho_eps_ye(rhoL, epsL, 0.5); 

        }
 
        // ---------- End of validity check

        rho(p.I)     = rhoL;
        velx(p.I)    = v_up(0);
        vely(p.I)    = v_up(1);
        velz(p.I)    = v_up(2);
        eps(p.I)     = epsL;
        press(p.I)   = pressL;
        entropy(p.I) = entropyL;

        saved_rho(p.I) = rhoL;
        saved_velx(p.I) = v_up(0);
        saved_vely(p.I) = v_up(1);
        saved_velz(p.I) = v_up(2);
        saved_eps(p.I) = epsL;

        CCTK_REAL wlor = calc_wlorentz(v_low, v_up);

	zvec_x(p.I) = wlor * v_up(0);
	zvec_y(p.I) = wlor * v_up(1);
	zvec_z(p.I) = wlor * v_up(2);

	svec_x(p.I) = (rhoL+rhoL*epsL+pressL)*wlor*wlor*v_up(0);
	svec_y(p.I) = (rhoL+rhoL*epsL+pressL)*wlor*wlor*v_up(1);
	svec_z(p.I) = (rhoL+rhoL*epsL+pressL)*wlor*wlor*v_up(2);

      });

}

} // namespace AsterX
