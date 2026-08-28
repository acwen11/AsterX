#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "aster_utils.hxx"

#include "../../../CarpetX/CarpetX/src/schedule.hxx"

namespace AsterX {
using namespace Loop;
using namespace AsterUtils;

extern "C" void AsterX_SetCellAverage(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetCellAverage;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1)/CCTK_REAL(24);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, 
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        const CCTK_REAL densL =  dens(p.I);
        const CCTK_REAL momxL =  momx(p.I);
        const CCTK_REAL momyL =  momy(p.I);
        const CCTK_REAL momzL =  momz(p.I);
        const CCTK_REAL tauL =  tau(p.I);
        const CCTK_REAL DYeL =  DYe(p.I);
        const CCTK_REAL DEntL =  DEnt(p.I);

        const bool c2pflag_stencil = (con2prim_flag(p.I) >= 6) 
          || (con2prim_flag(p.I + p.DI[0]) >= 6) || (con2prim_flag(p.I - p.DI[0]) >= 6)
          || (con2prim_flag(p.I + p.DI[1]) >= 6) || (con2prim_flag(p.I - p.DI[1]) >= 6)
          || (con2prim_flag(p.I + p.DI[2]) >= 6) || (con2prim_flag(p.I - p.DI[2]) >= 6);
        if ((con2prim_flag(p.I) >= 6) && ((LOflag(p.I) > 0.0) && shock_pv_fallback)) {
          dens(p.I) =  dens_pv(p.I);
          momx(p.I) =  momx_pv(p.I);
          momy(p.I) =  momy_pv(p.I);
          momz(p.I) =  momz_pv(p.I);
          tau(p.I)  =  tau_pv(p.I);
          DYe(p.I)  =  DYe_pv(p.I);
          DEnt(p.I) =  DEnt_pv(p.I);
        } else if (c2pflag_stencil || (cctk_iteration == 0)) {
          // Eq. (17) from https://arxiv.org/pdf/2310.11831 
          bool thetac = ((LOflag(p.I) == 0.0) || !shock_pv_fallback);
          dens(p.I) =  dens_pv(p.I) + thetac * one_over_24*laplace_3d(dens_pv,p);
          momx(p.I) =  momx_pv(p.I) + thetac * one_over_24*laplace_3d(momx_pv,p);
          momy(p.I) =  momy_pv(p.I) + thetac * one_over_24*laplace_3d(momy_pv,p);
          momz(p.I) =  momz_pv(p.I) + thetac * one_over_24*laplace_3d(momz_pv,p);
          tau(p.I)  =  tau_pv(p.I) + thetac * one_over_24*laplace_3d(tau_pv,p);
          DYe(p.I)  =  DYe_pv(p.I) + thetac * one_over_24*laplace_3d(DYe_pv,p);
          DEnt(p.I) =  DEnt_pv(p.I) + thetac * one_over_24*laplace_3d(DEnt_pv,p);
        } else {
          dens(p.I) =  densL;
          momx(p.I) =  momxL;
          momy(p.I) =  momyL;
          momz(p.I) =  momzL;
          tau(p.I)  =  tauL;
          DYe(p.I)  =  DYeL;
          DEnt(p.I) =  DEntL;
        }
      });

  /*
  grid.loop_outer_n_device<1, 1, 1>(
      grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        // Use 2nd order accurate conversion at boundary
        dens(p.I) =  dens_pv(p.I);
        momx(p.I) =  momx_pv(p.I);
        momy(p.I) =  momy_pv(p.I);
        momz(p.I) =  momz_pv(p.I);
        tau(p.I)  =  tau_pv(p.I);
        DYe(p.I)  =  DYe_pv(p.I);
        DEnt(p.I) =  DEnt_pv(p.I);

      });
  */
}

extern "C" void AsterX_InitPointValues(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_InitPointValues;
  DECLARE_CCTK_PARAMETERS;

  // if (cctk_iteration != 0) { // Point values are available at iteration 0
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
       
        // Set to zero to make them valid for CA C2P,
        // where they will not be used.
        dens_pv(p.I) = 0.0;
        momx_pv(p.I) = 0.0;
        momy_pv(p.I) = 0.0;
        momz_pv(p.I) = 0.0;
        tau_pv(p.I)  = 0.0;
        DYe_pv(p.I)  = 0.0;
        DEnt_pv(p.I) = 0.0;

      });
  // }
}

extern "C" void AsterX_AdjustConsPostStep(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_AdjustConsPostStep;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1)/CCTK_REAL(24);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
       
        if (LOflag(p.I) > 0.0) {
          // Use point values in flagged cells
          dens(p.I) = dens_pv(p.I);
          momx(p.I) = momx_pv(p.I);
          momy(p.I) = momy_pv(p.I);
          momz(p.I) = momz_pv(p.I);
          tau(p.I) = tau_pv(p.I);
          DYe(p.I) = DYe_pv(p.I);
          DEnt(p.I) = DEnt_pv(p.I);
        }
        else if ((LOflag(p.I) == 0.0) && (LOflag_p(p.I) > 0.0)) {
          // Re-average cells that are no longer flagged
          dens(p.I) = dens_pv(p.I) + one_over_24*laplace_3d(dens_pv,p);
          momx(p.I) = momx_pv(p.I) + one_over_24*laplace_3d(momx_pv,p);
          momy(p.I) = momy_pv(p.I) + one_over_24*laplace_3d(momy_pv,p);
          momz(p.I) = momz_pv(p.I) + one_over_24*laplace_3d(momz_pv,p);
          tau(p.I) = tau_pv(p.I) + one_over_24*laplace_3d(tau_pv,p);
          DYe(p.I) = DYe_pv(p.I) + one_over_24*laplace_3d(DYe_pv,p);
          DEnt(p.I) = DEnt_pv(p.I) + one_over_24*laplace_3d(DEnt_pv,p);
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
    grid.loop_all_device<1, 1, 1>(
       grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // Use HydroBaseX GFs as temporary helper
          rho(p.I) =  dens(p.I);
          velx(p.I) =  momx(p.I);
          vely(p.I) =  momy(p.I);
          velz(p.I) =  momz(p.I);
          eps(p.I)  =  tau(p.I);
          Ye(p.I)  =  DYe(p.I);
          entropy(p.I) =  DEnt(p.I);
        });
  } else {
    grid.loop_all_device<1, 1, 1>(
       grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) =  dens_pv(p.I);
          velx(p.I) =  momx_pv(p.I);
          vely(p.I) =  momy_pv(p.I);
          velz(p.I) =  momz_pv(p.I);
          eps(p.I)  =  tau_pv(p.I);
          Ye(p.I)  =  DYe_pv(p.I);
          entropy(p.I) =  DEnt_pv(p.I);
        });
  }
}

extern "C" void AsterX_SetPointValues(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetPointValues;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1)/CCTK_REAL(24);

  //if (cctk_iteration != 0) { // use ID values at iteration 0
  grid.loop_allmn_device<1, 1, 1>(
      grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        bool thetac = ((LOflag(p.I) == 0.0) || !shock_pv_fallback);
        // Use HydroBaseX GFs as temporary helper
        dens_pv(p.I) =  dens(p.I) - thetac * one_over_24*laplace_3d(rho,p);
        momx_pv(p.I) =  momx(p.I) - thetac * one_over_24*laplace_3d(velx,p);
        momy_pv(p.I) =  momy(p.I) - thetac * one_over_24*laplace_3d(vely,p);
        momz_pv(p.I) =  momz(p.I) - thetac * one_over_24*laplace_3d(velz,p);
        tau_pv(p.I)  =  tau(p.I) - thetac * one_over_24*laplace_3d(eps,p);
        DYe_pv(p.I)  =  DYe(p.I) - thetac * one_over_24*laplace_3d(Ye,p);
        DEnt_pv(p.I) =  DEnt(p.I) - thetac * one_over_24*laplace_3d(entropy,p);

      });

  grid.loop_outer_n_device<1, 1, 1>(
      grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

  // Use 2nd order accurate conversion at boundary
        dens_pv(p.I) =  dens(p.I);
        momx_pv(p.I) =  momx(p.I);
        momy_pv(p.I) =  momy(p.I);
        momz_pv(p.I) =  momz(p.I);
        tau_pv(p.I)  =  tau(p.I);
        DYe_pv(p.I)  =  DYe(p.I);
        DEnt_pv(p.I) =  DEnt(p.I);

      });
  //}
}

extern "C" void AsterX_DecConsIter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_DecConsIter;
  DECLARE_CCTK_PARAMETERS;

  *cons_pv_iter -= 1;

  int minghosts = min(cctk_nghostzones[0], min(cctk_nghostzones[1], cctk_nghostzones[2]));
  if (*cons_pv_iter==0 || ((n_conspv_iters - *cons_pv_iter) % minghosts == 0))  {
    static const std::vector<int> groups = {CCTK_GroupIndex("AsterX::cons_vector_pv")};

    if (use_subcycling)
      SyncGroupsByDirISubcycling(cctkGH, groups.size(), groups.data(), nullptr);
    else
      SyncGroupsByDirI(cctkGH, groups.size(), groups.data(), nullptr);
  }
}
/* END ITERATIVE POINT VALUED CONSERVATIVE VECTORY CALCULATION */

} // namespace
