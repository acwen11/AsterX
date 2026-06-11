#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "aster_utils.hxx"

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

extern "C" void AsterX_SetPointValues(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetPointValues;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1)/CCTK_REAL(24);

  if (cctk_iteration != 0) { // use ID values at iteration 0
         
    grid.loop_allmn_device<1, 1, 1>(
        grid.nghostzones, 1,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

	  // Eq. (16) from https://arxiv.org/pdf/2310.11831 
          bool thetac = (!LOflag(p.I) || !shock_pv_fallback);
          dens_pv(p.I) =  dens(p.I) - thetac * one_over_24*laplace_3d(dens,p);
          momx_pv(p.I) =  momx(p.I) - thetac * one_over_24*laplace_3d(momx,p);
          momy_pv(p.I) =  momy(p.I) - thetac * one_over_24*laplace_3d(momy,p);
          momz_pv(p.I) =  momz(p.I) - thetac * one_over_24*laplace_3d(momz,p);
          tau_pv(p.I)  =  tau(p.I) - thetac * one_over_24*laplace_3d(tau,p);
          DYe_pv(p.I)  =  DYe(p.I) - thetac * one_over_24*laplace_3d(DYe,p);
          DEnt_pv(p.I) =  DEnt(p.I) - thetac * one_over_24*laplace_3d(DEnt,p);

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

  } else {
    // Do nothing
    return;
  }
}

extern "C" void AsterX_PVConsFallback(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_PVConsFallback;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
       
        if (LOflag(p.I)) {
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

} // namespace
