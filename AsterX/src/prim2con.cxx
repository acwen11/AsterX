#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include "prim2con.hxx"

namespace AsterX {
using namespace Loop;
using namespace std;

extern "C" void AsterX_Prim2Con_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Prim2Con_Initial;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the entire grid (0 to n-1 cells in each direction)
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Interpolate metric terms from vertices to center
        const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                   calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                   calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};

        prim pv;
        pv.rho = rho(p.I);
        pv.vel(0) = velx(p.I);
        pv.vel(1) = vely(p.I);
        pv.vel(2) = velz(p.I);
        pv.eps = eps(p.I);
        pv.press = press(p.I);
        pv.entropy = entropy(p.I);
        pv.Ye = Ye(p.I);
        pv.Bvec(0) = Bvecx(p.I);
        pv.Bvec(1) = Bvecy(p.I);
        pv.Bvec(2) = Bvecz(p.I);

        cons cv;
        prim2con(g, pv, cv);

        dens(p.I) = cv.dens;
        momx(p.I) = cv.mom(0);
        momy(p.I) = cv.mom(1);
        momz(p.I) = cv.mom(2);
        tau(p.I) = cv.tau;
        DYe(p.I) = cv.DYe;
        dBx(p.I) = cv.dBvec(0);
        dBy(p.I) = cv.dBvec(1);
        dBz(p.I) = cv.dBvec(2);
        DEnt(p.I) = cv.DEnt;

        if (use_ho_fv) {
          dens_pv(p.I) = cv.dens;
          momx_pv(p.I) = cv.mom(0);
          momy_pv(p.I) = cv.mom(1);
          momz_pv(p.I) = cv.mom(2);
          tau_pv(p.I) = cv.tau;
          DYe_pv(p.I) = cv.DYe;
          DEnt_pv(p.I) = cv.DEnt;
	      } else {
          dens_pv(p.I) = 0.0;
          momx_pv(p.I) = 0.0;
          momy_pv(p.I) = 0.0;
          momz_pv(p.I) = 0.0;
          tau_pv(p.I) = 0.0;
          DYe_pv(p.I) = 0.0;
          DEnt_pv(p.I) = 0.0;
        }

      });
}

extern "C" void AsterX_LineAvgAvec_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_LineAvgAvec_Initial;
  DECLARE_CCTK_PARAMETERS;

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1)/CCTK_REAL(24);
  // Use E GFs as temporary storage
  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { 
          Ex(p.I) = Avec_x(p.I) + one_over_24 * laplace_1d<0>(Avec_x, p);
      });
  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { 
          Ey(p.I) = Avec_y(p.I) + one_over_24 * laplace_1d<1>(Avec_y, p);
      });
  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { 
          Ez(p.I) = Avec_z(p.I) + one_over_24 * laplace_1d<2>(Avec_z, p);
      });
}

extern "C" void AsterX_CopyLineAvgAvec_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CopyLineAvgAvec_Initial;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { 
          Avec_x(p.I) = Ex(p.I);
      });
  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { 
          Avec_y(p.I) = Ey(p.I);
      });
  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { 
          Avec_z(p.I) = Ez(p.I);
      });
}

extern "C" void AsterX_PsiZero_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_PsiZero_Initial;
  DECLARE_CCTK_PARAMETERS;

  /* Initilaize Psi to 0.0 */
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Psi(p.I) = 0.0; });
}

extern "C" void AsterX_FluxAuxZero_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_FluxAuxZero_Initial;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { Fx_stag(p.I) = 0.0; });

  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { Fy_stag(p.I) = 0.0; });

  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { Fz_stag(p.I) = 0.0; });
}

} // namespace AsterX
