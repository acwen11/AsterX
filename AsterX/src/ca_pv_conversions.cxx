#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "aster_utils.hxx"
#include "prim2con.hxx"
#include "setup_eos.hxx"

#include "../../../CarpetX/CarpetX/src/schedule.hxx"

namespace AsterX {
using namespace Loop;
using namespace EOSX;
using namespace Arith;
using namespace AsterUtils;

extern "C" void AsterX_SetCellAverage(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetCellAverage;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1) / CCTK_REAL(24);

  grid.loop_allmn_device<1, 1, 1>(
      grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL densL = dens(p.I);
        const CCTK_REAL momxL = momx(p.I);
        const CCTK_REAL momyL = momy(p.I);
        const CCTK_REAL momzL = momz(p.I);
        const CCTK_REAL tauL = tau(p.I);
        const CCTK_REAL DYeL = DYe(p.I);
        const CCTK_REAL DEntL = DEnt(p.I);

        const bool c2pflag_stencil = (con2prim_flag(p.I) >= 6) ||
                                     (con2prim_flag(p.I + p.DI[0]) >= 6) ||
                                     (con2prim_flag(p.I - p.DI[0]) >= 6) ||
                                     (con2prim_flag(p.I + p.DI[1]) >= 6) ||
                                     (con2prim_flag(p.I - p.DI[1]) >= 6) ||
                                     (con2prim_flag(p.I + p.DI[2]) >= 6) ||
                                     (con2prim_flag(p.I - p.DI[2]) >= 6);
        if ((con2prim_flag(p.I) >= 6) &&
            ((LOflag(p.I) > 0.0) && shock_pv_fallback)) {
          dens(p.I) = dens_pv(p.I);
          momx(p.I) = momx_pv(p.I);
          momy(p.I) = momy_pv(p.I);
          momz(p.I) = momz_pv(p.I);
          tau(p.I) = tau_pv(p.I);
          DYe(p.I) = DYe_pv(p.I);
          DEnt(p.I) = DEnt_pv(p.I);
        } else if (c2pflag_stencil || (cctk_iteration == 0)) {
          // Eq. (17) from https://arxiv.org/pdf/2310.11831
          bool thetac = ((LOflag(p.I) == 0.0) || !shock_pv_fallback);
          dens(p.I) =
              dens_pv(p.I) + thetac * one_over_24 * laplace_3d(dens_pv, p);
          momx(p.I) =
              momx_pv(p.I) + thetac * one_over_24 * laplace_3d(momx_pv, p);
          momy(p.I) =
              momy_pv(p.I) + thetac * one_over_24 * laplace_3d(momy_pv, p);
          momz(p.I) =
              momz_pv(p.I) + thetac * one_over_24 * laplace_3d(momz_pv, p);
          tau(p.I) = tau_pv(p.I) + thetac * one_over_24 * laplace_3d(tau_pv, p);
          DYe(p.I) = DYe_pv(p.I) + thetac * one_over_24 * laplace_3d(DYe_pv, p);
          DEnt(p.I) =
              DEnt_pv(p.I) + thetac * one_over_24 * laplace_3d(DEnt_pv, p);
        } else {
          dens(p.I) = densL;
          momx(p.I) = momxL;
          momy(p.I) = momyL;
          momz(p.I) = momzL;
          tau(p.I) = tauL;
          DYe(p.I) = DYeL;
          DEnt(p.I) = DEntL;
        }
      });

  grid.loop_outer_n_device<1, 1, 1>(grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            // Use 2nd order accurate conversion
            // at boundary
            dens(p.I) = dens_pv(p.I);
            momx(p.I) = momx_pv(p.I);
            momy(p.I) = momy_pv(p.I);
            momz(p.I) = momz_pv(p.I);
            tau(p.I) = tau_pv(p.I);
            DYe(p.I) = DYe_pv(p.I);
            DEnt(p.I) = DEnt_pv(p.I);
      });
}

extern "C" void AsterX_InitPointValues(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_InitPointValues;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            dens_pv(p.I) = 0.0;
            momx_pv(p.I) = 0.0;
            momy_pv(p.I) = 0.0;
            momz_pv(p.I) = 0.0;
            tau_pv(p.I) = 0.0;
            DYe_pv(p.I) = 0.0;
            DEnt_pv(p.I) = 0.0;
      });
}

/*
extern "C" void AsterX_CopySavedPrims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CopySavedPrims;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            saved_rho(p.I) = rho(p.I);
            saved_velx(p.I) = velx(p.I);
            saved_vely(p.I) = vely(p.I);
            saved_velz(p.I) = velz(p.I);
            // Use eps to store T
            saved_eps(p.I) = temperature(p.I);
            saved_Ye(p.I) = Ye(p.I);
      });
}

extern "C" void AsterX_AdjustCA(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_AdjustCA;
  DECLARE_CCTK_PARAMETERS;

  if (cctk_iteration > 0) {
    constexpr CCTK_REAL one_over_24 = CCTK_REAL(1) / CCTK_REAL(24);

    grid.loop_all_device<1, 1, 1>(grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (LOflag(p.I) > LOflag_p(p.I)) {
            dens(p.I) = pvdens_frac * dens_pv(p.I)  + (1.0 - pvdens_frac) * dens(p.I) ;
            momx(p.I) = momx_pv(p.I);
            momy(p.I) = momy_pv(p.I);
            momz(p.I) = momz_pv(p.I);
            tau(p.I) = tau_pv(p.I);
            DYe(p.I) =  pvdens_frac * DYe_pv(p.I)  + (1.0 - pvdens_frac) * DYe(p.I) ;
            DEnt(p.I) =  pvdens_frac * DEnt_pv(p.I)  + (1.0 - pvdens_frac) * DEnt(p.I) ;
          }
        });
  }
}
*/

extern "C" void AsterX_AdjustPV(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_AdjustPV;
  DECLARE_CCTK_PARAMETERS;

  if (cctk_iteration > 0) {
    constexpr CCTK_REAL one_over_24 = CCTK_REAL(1) / CCTK_REAL(24);

    grid.loop_int_device<1, 1, 1>(grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // bool reavg = false;
          // for (int dir = 0; dir < 3; dir++) {
          //   for (int step = -1; step <= 1; step++) {
          //     const auto idx = p.I + step * p.DI[dir];
          //     reavg = reavg || (LOflag(idx) < LOflag_p(idx));
          //   }
          // }
          // reavg = reavg && (LOflag(p.I) == 0.0);
          if (LOflag(p.I) > LOflag_p(p.I)) {
            dens_pv(p.I) = dens(p.I);
            momx_pv(p.I) = momx(p.I);
            momy_pv(p.I) = momy(p.I);
            momz_pv(p.I) = momz(p.I);
            tau_pv(p.I) = tau(p.I);
            DYe_pv(p.I) = DYe(p.I);
            DEnt_pv(p.I) = DEnt(p.I);
          }
          else if (LOflag(p.I) < LOflag_p(p.I)) {
            dens_pv(p.I) = dens(p.I) - one_over_24 * laplace_3d(dens, p);
            momx_pv(p.I) = momx(p.I) - one_over_24 * laplace_3d(momx, p);
            momy_pv(p.I) = momy(p.I) - one_over_24 * laplace_3d(momy, p);
            momz_pv(p.I) = momz(p.I) - one_over_24 * laplace_3d(momz, p);
            tau_pv(p.I) = tau(p.I) - one_over_24 * laplace_3d(tau, p);
            DYe_pv(p.I) = DYe(p.I) - one_over_24 * laplace_3d(DYe, p);
            DEnt_pv(p.I) = DEnt(p.I) - one_over_24 * laplace_3d(DEnt, p);
          }
            
        });
  }
}

/* KEEP KE = rho*h*(Wv)^2 CONSTANT */
extern "C" void AsterX_AdjustPrimsPostStep(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_AdjustPrimsPostStep;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1) / CCTK_REAL(24);
  // Get local eos objects
  // auto eos_1p = global_eos_1p_poly;
  auto eos_3p = global_eos_3p_tab3d;

  grid.loop_allmn_device<1, 1, 1>(
      grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        bool reavg = false;
        for (int dir = 0; dir < 3; dir++) {
          for (int step = -1; step <= 1; step++) {
            const auto idx = p.I + step * p.DI[dir];
            reavg = reavg || (LOflag(idx) < LOflag_p(idx));
          }
        }
        reavg = reavg && (LOflag(p.I) == 0.0);
        if (reavg || (LOflag(p.I) > LOflag_p(p.I))) {
          // Recalculate point value cons from prims in neighbors of newly unflagged cells
          const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                     calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                     calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};
          const CCTK_REAL sqrt_detg = sqrt(calc_det(g));

          prim pv;
          pv.Bvec(0) = Bvecx(p.I);
          pv.Bvec(1) = Bvecy(p.I);
          pv.Bvec(2) = Bvecz(p.I);

          /*
          // Rescale v^i to preserve v^2, W
          vec<CCTK_REAL, 3> v_up{velx(p.I), vely(p.I), velz(p.I)};
          const vec<CCTK_REAL, 3> s_v_up{saved_velx(p.I), saved_vely(p.I), saved_velz(p.I)};
          const CCTK_REAL vnorm = sqrt(calc_norm(v_up, g));
          const CCTK_REAL s_vnorm = sqrt(calc_norm(s_v_up, g));
          const CCTK_REAL fac = s_vnorm / (vnorm + 1e-15);
          if (vnorm > 0.0) {
            for (int ii = 0; ii < 3; ii++) {
              v_up(ii) *= fac;
            }
          }
          else {
            for (int ii = 0; ii < 3; ii++) {
              v_up(ii) = s_v_up(ii);
            }
          }
          const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);
          const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);
          const CCTK_REAL vsq = calc_contraction(v_low, v_up);
          pv.vel(0) = v_up(0);
          pv.vel(1) = v_up(1);
          pv.vel(2) = v_up(2);
          */
          pv.vel(0) = velx(p.I);
          pv.vel(1) = vely(p.I);
          pv.vel(2) = velz(p.I);
          vec<CCTK_REAL, 3> v_up{velx(p.I), vely(p.I), velz(p.I)};
          const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);
          const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

          // Set prims from T00 or kinetic energy
          // const CCTK_REAL s_rho = rho(p.I);
          // const CCTK_REAL s_eps = eps(p.I);
          // const CCTK_REAL s_press = press(p.I);
          // const CCTK_REAL s_rhoh = s_rho + s_rho * s_eps + s_press;

          // const CCTK_REAL s_temp = saved_eps(p.I);
          // const CCTK_REAL s_Ye = Ye(p.I);

          // Set prims from tau_pv

          /* Computing B_j */
          const vec<CCTK_REAL, 3> &B_up = pv.Bvec;
          const vec<CCTK_REAL, 3> B_low = calc_contraction(g, B_up);

          /* Computing b^t : this is b^0 * alp */
          const CCTK_REAL bst = w_lorentz * calc_contraction(B_up, v_low);

          /* Computing b^mu b_mu */
          const CCTK_REAL bs2 =
              (calc_contraction(B_up, B_low) + bst * bst) / (w_lorentz * w_lorentz);

          const CCTK_REAL tau_target = (tau_pv(p.I) / sqrt_detg) - bs2 * w_lorentz * w_lorentz + 0.5 * bs2 + bst * bst;

          // Option 1: Fix T, Ye; adjust rho
          const CCTK_REAL tempL = temperature(p.I);
          const CCTK_REAL YeL = Ye(p.I);
          const CCTK_REAL rhoL = eos_3p->rho_from_tautarget_temp_ye(tau_target, w_lorentz, tempL, YeL);
          // Option 2: Fix Dens, Ye; adjust temp
          // const CCTK_REAL rhoL = dens(p.I) / (w_lorentz * sqrt_detg);
          // const CCTK_REAL YeL = s_Ye;
          // const CCTK_REAL tempL = eos_3p->temp_from_rhoh_rho_ye(s_rhoh, rhoL, YeL);
          // Option 3: Fix rho, Ye; adjust temp
          // const CCTK_REAL rhoL = 0.5 * (s_rho + dens(p.I) / (w_lorentz * sqrt_detg));
          // const CCTK_REAL YeL = s_Ye;
          // const CCTK_REAL tempL = eos_3p->temp_from_rhoh_rho_ye(s_rhoh, rhoL, YeL);
          // CCTK_VINFO("saved temp = %g; new temp = %g; adjusted temp = %g", s_temp, temperature(p.I), tempL);
          // CCTK_VINFO("saved rho = %g; new rho = %g; adjusted rho = %g", s_rho, rho(p.I), rhoL);

          pv.rho = rhoL;
          pv.Ye = YeL;
          pv.eps = eos_3p->eps_from_rho_temp_ye(rhoL, tempL, YeL);
          pv.press = eos_3p->press_from_rho_temp_ye(rhoL, tempL, YeL);
          pv.entropy = eos_3p->entropy_from_rho_temp_ye(rhoL, tempL, YeL);

          // Update prims
          rho(p.I) = pv.rho;
          velx(p.I) = pv.vel(0);
          vely(p.I) = pv.vel(1);
          velz(p.I) = pv.vel(2);
          temperature(p.I) = tempL;
          // Ye(p.I) = YeL;
          press(p.I) = pv.press;
          eps(p.I) = pv.eps;
          entropy(p.I) = pv.entropy;

          // Recalculate cons from prims
          cons cv;
          prim2con(g, pv, cv);

          dens_pv(p.I) = cv.dens;
          momx_pv(p.I) = cv.mom(0);
          momy_pv(p.I) = cv.mom(1);
          momz_pv(p.I) = cv.mom(2);
          tau_pv(p.I) = cv.tau;
          DYe_pv(p.I) = cv.DYe;
          DEnt_pv(p.I) = cv.DEnt;
        } 
      });
}

/* KEEP ALL C2P OUTPUT, RESET TEMPERATURE
extern "C" void AsterX_AdjustPrimsPostStep(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_AdjustPrimsPostStep;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1) / CCTK_REAL(24);
  // Get local eos objects
  // auto eos_1p = global_eos_1p_poly;
  auto eos_3p = global_eos_3p_tab3d;
  grid.loop_allmn_device<1, 1, 1>(
      grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        bool reavg = false;
        for (int dir = 0; dir < 3; dir++) {
          for (int step = -1; step <= 1; step++) {
            const auto idx = p.I + step * p.DI[dir];
            reavg = reavg || (LOflag(idx) < LOflag_p(idx));
          }
        }
        reavg = reavg && (LOflag(p.I) == 0.0);
        if (reavg || (LOflag(p.I) > LOflag_p(p.I))) {
          const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                     calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                     calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};
          // // Recalculate point value cons from prims in newly unflagged cells
          // // while keeping dens fixed
          // const CCTK_REAL sqrt_detg = sqrt(calc_det(g));

          prim pv;
          pv.vel(0) = velx(p.I);
          pv.vel(1) = vely(p.I);
          pv.vel(2) = velz(p.I);
          pv.Bvec(0) = Bvecx(p.I);
          pv.Bvec(1) = Bvecy(p.I);
          pv.Bvec(2) = Bvecz(p.I);
          pv.rho = rho(p.I);
          pv.Ye = Ye(p.I);
          // Recalculate prims from saved temperature
          const CCTK_REAL tempL = exp(saved_T_frac * log(saved_eps(p.I)) + (1.0 - saved_T_frac) * temperature(p.I));
          // const CCTK_REAL tempL = saved_eps(p.I);
          pv.eps = eos_3p->eps_from_rho_temp_ye(pv.rho, tempL, pv.Ye);
          pv.press = eos_3p->press_from_rho_temp_ye(pv.rho, tempL, pv.Ye);
          pv.entropy = eos_3p->entropy_from_rho_temp_ye(pv.rho, tempL, pv.Ye);

          temperature(p.I) = tempL;
          eps(p.I) = pv.eps;
          press(p.I) = pv.press;
          entropy(p.I) = pv.entropy;

          // Recalculate PV cons from prims
          cons cv;
          prim2con(g, pv, cv);
          dens_pv(p.I) = cv.dens;
          momx_pv(p.I) = cv.mom(0);
          momy_pv(p.I) = cv.mom(1);
          momz_pv(p.I) = cv.mom(2);
          tau_pv(p.I) = cv.tau;
          DYe_pv(p.I) = cv.DYe;
          DEnt_pv(p.I) = cv.DEnt;
        } 
      });
}
*/

extern "C" void AsterX_ReAvgConsPostStep(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ReAvgConsPostStep;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1) / CCTK_REAL(24);
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // bool c2pflag_stencil = false;
        bool reavg = false;
        for (int dk = -2; dk <= 2; dk++) {
          for (int dj = -2; dj <= 2; dj++) {
            for (int di = -2; di <= 2; di++) {
              auto idx = p.I + di * p.DI[0] + dj * p.DI[1] + dk * p.DI[2];
              // if (abs(dk) + abs(dj) + abs(di) <= 1)
              //   c2pflag_stencil = c2pflag_stencil || (con2prim_flag(idx) >= 6);
              if (abs(dk) + abs(dj) + abs(di) <= 2) 
                reavg = reavg || (LOflag(idx) < LOflag_p(idx));
              else
                continue;
            }
          }
        }
        reavg = reavg && (LOflag(p.I) == 0.0);
        // c2pflag_stencil = c2pflag_stencil && (LOflag(p.I) == 0.0);
        if (reavg) {
          // Re-average cells that are no longer flagged
          dens(p.I) = dens_pv(p.I) + one_over_24 * laplace_3d(dens_pv, p);
          momx(p.I) = momx_pv(p.I) + one_over_24 * laplace_3d(momx_pv, p);
          momy(p.I) = momy_pv(p.I) + one_over_24 * laplace_3d(momy_pv, p);
          momz(p.I) = momz_pv(p.I) + one_over_24 * laplace_3d(momz_pv, p);
          tau(p.I) = tau_pv(p.I) + one_over_24 * laplace_3d(tau_pv, p);
          DYe(p.I) = DYe_pv(p.I) + one_over_24 * laplace_3d(DYe_pv, p);
          DEnt(p.I) = DEnt_pv(p.I) + one_over_24 * laplace_3d(DEnt_pv, p);
        } 
        else if (LOflag(p.I) > LOflag_p(p.I)) {
          dens(p.I) = dens_pv(p.I);
          momx(p.I) = momx_pv(p.I);
          momy(p.I) = momy_pv(p.I);
          momz(p.I) = momz_pv(p.I);
          tau(p.I) = tau_pv(p.I);
          DYe(p.I) = DYe_pv(p.I);
          DEnt(p.I) = DEnt_pv(p.I);
        }

      });
}

/* BEGIN ITERATIVE POINT VALUED CONSERVATIVES CALCULATION */
extern "C" void AsterX_SetConsIter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetConsIter;
  DECLARE_CCTK_PARAMETERS;

  *cons_pv_iter = n_conspv_iters;
}

extern "C" void AsterX_SetPVConsnMinus1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetPVConsnMinus1;
  DECLARE_CCTK_PARAMETERS;

  if (*cons_pv_iter == n_conspv_iters) {
    grid.loop_all_device<1, 1, 1>(grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
              // Use RHS GFs as temporary
              // helper
              densrhs(p.I) = dens(p.I);
              momxrhs(p.I) = momx(p.I);
              momyrhs(p.I) = momy(p.I);
              momzrhs(p.I) = momz(p.I);
              taurhs(p.I) = tau(p.I);
              DYe_rhs(p.I) = DYe(p.I);
              DEntrhs(p.I) = DEnt(p.I);
        });
  } else {
    grid.loop_all_device<1, 1, 1>(grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
              densrhs(p.I) = dens_pv(p.I);
              momxrhs(p.I) = momx_pv(p.I);
              momyrhs(p.I) = momy_pv(p.I);
              momzrhs(p.I) = momz_pv(p.I);
              taurhs(p.I) = tau_pv(p.I);
              DYe_rhs(p.I) = DYe_pv(p.I);
              DEntrhs(p.I) = DEnt_pv(p.I);
        });
  }
}

extern "C" void AsterX_SetPointValues(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetPointValues;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1) / CCTK_REAL(24);

  grid.loop_allmn_device<1, 1, 1>(
      grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        bool thetac = ((LOflag(p.I) == 0.0) || !shock_pv_fallback);
        // Use RHS GFs as temporary helper
        dens_pv(p.I) = dens(p.I) - thetac * one_over_24 * laplace_3d(densrhs, p);
        momx_pv(p.I) = momx(p.I) - thetac * one_over_24 * laplace_3d(momxrhs, p);
        momy_pv(p.I) = momy(p.I) - thetac * one_over_24 * laplace_3d(momyrhs, p);
        momz_pv(p.I) = momz(p.I) - thetac * one_over_24 * laplace_3d(momzrhs, p);
        tau_pv(p.I) = tau(p.I) - thetac * one_over_24 * laplace_3d(taurhs, p);
        DYe_pv(p.I) = DYe(p.I) - thetac * one_over_24 * laplace_3d(DYe_rhs, p);
        DEnt_pv(p.I) =
            DEnt(p.I) - thetac * one_over_24 * laplace_3d(DEntrhs, p);
      });

  grid.loop_outer_n_device<1, 1, 1>(grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            // Use 2nd order accurate conversion at boundary
            dens_pv(p.I) = dens(p.I);
            momx_pv(p.I) = momx(p.I);
            momy_pv(p.I) = momy(p.I);
            momz_pv(p.I) = momz(p.I);
            tau_pv(p.I) = tau(p.I);
            DYe_pv(p.I) = DYe(p.I);
            DEnt_pv(p.I) = DEnt(p.I);
      });
}

extern "C" void AsterX_DecConsIter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_DecConsIter;
  DECLARE_CCTK_PARAMETERS;

  *cons_pv_iter -= 1;

  int minghosts =
      min(cctk_nghostzones[0], min(cctk_nghostzones[1], cctk_nghostzones[2]));
  if (*cons_pv_iter == 0 ||
      ((n_conspv_iters - *cons_pv_iter) % minghosts == 0)) {
    static const std::vector<int> groups = {
        CCTK_GroupIndex("AsterX::cons_vector_pv")};

    if (use_subcycling)
      SyncGroupsByDirISubcycling(cctkGH, groups.size(), groups.data(), nullptr);
    else
      SyncGroupsByDirI(cctkGH, groups.size(), groups.data(), nullptr);
  }
}
/* END ITERATIVE POINT VALUED CONSERVATIVE VECTORY CALCULATION */

} // namespace AsterX
