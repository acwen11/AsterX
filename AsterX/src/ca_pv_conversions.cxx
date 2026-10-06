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

extern "C" void AsterX_CA2PVPostStep(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CA2PVPostStep;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1) / CCTK_REAL(24);
  // Get local eos objects
  // auto eos_1p = global_eos_1p_poly;
  auto eos_3p = global_eos_3p_tab3d;
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if ((LOflag(p.I) > 0.0) && (LOflag_p(p.I) == 0.0)) {
          // Recalculate point value cons from prims in newly flagged cells
          // while keeping dens fixed
          const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                     calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                     calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};
          const CCTK_REAL sqrt_detg = sqrt(calc_det(g));

          prim pv;
          pv.vel(0) = velx(p.I);
          pv.vel(1) = vely(p.I);
          pv.vel(2) = velz(p.I);
          pv.Bvec(0) = Bvecx(p.I);
          pv.Bvec(1) = Bvecy(p.I);
          pv.Bvec(2) = Bvecz(p.I);
          const vec<CCTK_REAL, 3> &v_up = pv.vel;
          const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);
          const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

          // Recalculate prims from dens
          const CCTK_REAL rhoL = dens(p.I) / (sqrt_detg * w_lorentz);
          const CCTK_REAL tempL = temperature(p.I);
          const CCTK_REAL YeL = Ye(p.I);
          pv.rho = rhoL;
          pv.eps = eos_3p->eps_from_rho_temp_ye(rhoL, tempL, YeL);
          pv.press = eos_3p->press_from_rho_temp_ye(rhoL, tempL, YeL);
          pv.entropy = eos_3p->entropy_from_rho_temp_ye(rhoL, tempL, YeL);
          pv.Ye = Ye(p.I);
          rho(p.I) = pv.rho;
          press(p.I) = pv.press;
          eps(p.I) = pv.eps;
          entropy(p.I) = pv.entropy;

          // Recalculate cons from prims
          cons cv;
          prim2con(g, pv, cv);
          dens(p.I) = cv.dens;
          momx(p.I) = cv.mom(0);
          momy(p.I) = cv.mom(1);
          momz(p.I) = cv.mom(2);
          tau(p.I) = cv.tau;
          DYe(p.I) = cv.DYe;
          DEnt(p.I) = cv.DEnt;
        } 
      });
}

extern "C" void AsterX_PV2CAPostStep(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_PV2CAPostStep;
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
            reavg = reavg || (LOflag(idx) != LOflag_p(idx));
          }
        }
        reavg = reavg && (LOflag(p.I) == 0.0);
        if (reavg) {
          const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                     calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                     calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};
          // Recalculate point value cons from prims in newly unflagged cells
          // while keeping dens fixed
          const CCTK_REAL sqrt_detg = sqrt(calc_det(g));

          prim pv;
          pv.vel(0) = velx(p.I);
          pv.vel(1) = vely(p.I);
          pv.vel(2) = velz(p.I);
          pv.Bvec(0) = Bvecx(p.I);
          pv.Bvec(1) = Bvecy(p.I);
          pv.Bvec(2) = Bvecz(p.I);
          const vec<CCTK_REAL, 3> &v_up = pv.vel;
          const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);
          const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

          // Recalculate prims from dens_pv
          const CCTK_REAL rhoL = dens_pv(p.I) / (sqrt_detg * w_lorentz);
          const CCTK_REAL tempL = temperature(p.I);
          const CCTK_REAL YeL = Ye(p.I);
          pv.rho = rhoL;
          pv.eps = eos_3p->eps_from_rho_temp_ye(rhoL, tempL, YeL);
          pv.press = eos_3p->press_from_rho_temp_ye(rhoL, tempL, YeL);
          pv.entropy = eos_3p->entropy_from_rho_temp_ye(rhoL, tempL, YeL);
          pv.Ye = Ye(p.I);
          rho(p.I) = pv.rho;
          press(p.I) = pv.press;
          eps(p.I) = pv.eps;
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

extern "C" void AsterX_ReAvgConsPostStep(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ReAvgConsPostStep;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1) / CCTK_REAL(24);
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        bool reavg = false;
        for (int dk = -2; dk <= 2; dk++) {
          for (int dj = -2; dj <= 2; dj++) {
            for (int di = -2; di <= 2; di++) {
              if (abs(dk) + abs(dj) + abs(di) > 2) 
                continue;
              else {
                const auto idx = p.I + di * p.DI[0] + dj * p.DI[1] + dk * p.DI[2];
                reavg = reavg || (LOflag(idx) != LOflag_p(idx));
              }
            }
          }
        }
        reavg = reavg && (LOflag(p.I) == 0.0);
        if (reavg) {
          // Re-average cells that are no longer flagged
          // dens(p.I) = dens(p.I) + one_over_24 * laplace_3d(dens, p);
          momx(p.I) = momx_pv(p.I) + one_over_24 * laplace_3d(momx_pv, p);
          momy(p.I) = momy_pv(p.I) + one_over_24 * laplace_3d(momy_pv, p);
          momz(p.I) = momz_pv(p.I) + one_over_24 * laplace_3d(momz_pv, p);
          tau(p.I) = tau_pv(p.I) + one_over_24 * laplace_3d(tau_pv, p);
          DYe(p.I) = DYe_pv(p.I) + one_over_24 * laplace_3d(DYe_pv, p);
          DEnt(p.I) = DEnt_pv(p.I) + one_over_24 * laplace_3d(DEnt_pv, p);
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
