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

  grid.loop_allmn_device<1, 1, 1>(
      grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const bool c2pflag_stencil = (con2prim_flag(p.I) >= 6) 
          || (con2prim_flag(p.I + p.DI[0]) >= 6) || (con2prim_flag(p.I - p.DI[0]) >= 6)
          || (con2prim_flag(p.I + p.DI[1]) >= 6) || (con2prim_flag(p.I - p.DI[1]) >= 6)
          || (con2prim_flag(p.I + p.DI[2]) >= 6) || (con2prim_flag(p.I - p.DI[2]) >= 6);
        if ((con2prim_flag(p.I) >= 6) && (LOflag(p.I) && shock_pv_fallback)) {
          dens(p.I) =  dens_pv(p.I);
          momx(p.I) =  momx_pv(p.I);
          momy(p.I) =  momy_pv(p.I);
          momz(p.I) =  momz_pv(p.I);
          tau(p.I)  =  tau_pv(p.I);
          DYe(p.I)  =  DYe_pv(p.I);
          DEnt(p.I) =  DEnt_pv(p.I);
        } else if (c2pflag_stencil || (cctk_iteration == 0)) {
          // Eq. (17) from https://arxiv.org/pdf/2310.11831 
          bool thetac = (!LOflag(p.I) || !shock_pv_fallback);
          dens(p.I) =  dens_pv(p.I) + thetac * one_over_24*laplace_3d(dens_pv,p);
          momx(p.I) =  momx_pv(p.I) + thetac * one_over_24*laplace_3d(momx_pv,p);
          momy(p.I) =  momy_pv(p.I) + thetac * one_over_24*laplace_3d(momy_pv,p);
          momz(p.I) =  momz_pv(p.I) + thetac * one_over_24*laplace_3d(momz_pv,p);
          tau(p.I)  =  tau_pv(p.I) + thetac * one_over_24*laplace_3d(tau_pv,p);
          DYe(p.I)  =  DYe_pv(p.I) + thetac * one_over_24*laplace_3d(DYe_pv,p);
          DEnt(p.I) =  DEnt_pv(p.I) + thetac * one_over_24*laplace_3d(DEnt_pv,p);
        }
        // else do nothing, leave cons as is
      });

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
}

extern "C" void AsterX_InitPointValues(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_InitPointValues;
  DECLARE_CCTK_PARAMETERS;

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
}

extern "C" void AsterX_AdjustConsPostStep(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_AdjustConsPostStep;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1)/CCTK_REAL(24);

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
       
        if (LOflag(p.I)) {
          // Use point values in flagged cells
          dens(p.I) = dens_pv(p.I);
          momx(p.I) = momx_pv(p.I);
          momy(p.I) = momy_pv(p.I);
          momz(p.I) = momz_pv(p.I);
          tau(p.I) = tau_pv(p.I);
          DYe(p.I) = DYe_pv(p.I);
          DEnt(p.I) = DEnt_pv(p.I);
        }
        else if (!LOflag(p.I) && LOflag_p(p.I)) {
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
          dens_aux(p.I) =  dens(p.I);
          momx_aux(p.I) =  momx(p.I);
          momy_aux(p.I) =  momy(p.I);
          momz_aux(p.I) =  momz(p.I);
          tau_aux(p.I)  =  tau(p.I);
          DYe_aux(p.I)  =  DYe(p.I);
          DEnt_aux(p.I) =  DEnt(p.I);
        });
  } else {
    grid.loop_all_device<1, 1, 1>(
       grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          dens_aux(p.I) =  dens_pv(p.I);
          momx_aux(p.I) =  momx_pv(p.I);
          momy_aux(p.I) =  momy_pv(p.I);
          momz_aux(p.I) =  momz_pv(p.I);
          tau_aux(p.I)  =  tau_pv(p.I);
          DYe_aux(p.I)  =  DYe_pv(p.I);
          DEnt_aux(p.I) =  DEnt_pv(p.I);
        });
  }
}

extern "C" void AsterX_SetPointValues(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetPointValues;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1)/CCTK_REAL(24);

  if (cctk_iteration != 0) { // use ID values at iteration 0
    grid.loop_allmn_device<1, 1, 1>(
        grid.nghostzones, 1,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          bool thetac = (!LOflag(p.I) || !shock_pv_fallback);
          dens_pv(p.I) =  dens(p.I) - thetac * one_over_24*laplace_3d(dens_aux,p);
          momx_pv(p.I) =  momx(p.I) - thetac * one_over_24*laplace_3d(momx_aux,p);
          momy_pv(p.I) =  momy(p.I) - thetac * one_over_24*laplace_3d(momy_aux,p);
          momz_pv(p.I) =  momz(p.I) - thetac * one_over_24*laplace_3d(momz_aux,p);
          tau_pv(p.I)  =  tau(p.I) - thetac * one_over_24*laplace_3d(tau_aux,p);
          DYe_pv(p.I)  =  DYe(p.I) - thetac * one_over_24*laplace_3d(DYe_aux,p);
          DEnt_pv(p.I) =  DEnt(p.I) - thetac * one_over_24*laplace_3d(DEnt_aux,p);

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
  }
}

extern "C" void AsterX_DecConsIter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_DecConsIter;
  DECLARE_CCTK_PARAMETERS;

  *cons_pv_iter -= 1;

  int minghosts = min(cctk_nghostzones[0], min(cctk_nghostzones[1], cctk_nghostzones[2]));
  if (*cons_pv_iter && ((n_conspv_iters - *cons_pv_iter) % minghosts == 0))  {
    static const std::vector<int> groups = {CCTK_GroupIndex("AsterX::cons_vector_pv")};

    SyncGroupsByDirISubcycling(cctkGH, groups.size(), groups.data(), nullptr);
  }
}
/* END ITERATIVE POINT VALUED CONSERVATIVE VECTORY CALCULATION */

} // namespace
