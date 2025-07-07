#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <stdlib.h>
#include <time.h>

#include "seeds_utils.hxx"

namespace AsterSeeds {
using namespace std;
using namespace amrex;
using namespace Loop;
using namespace AsterUtils;

extern "C" void AsterSeeds_InitializeCenteredAvec_TOV(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_InitializeCenteredAvec_TOV;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(Afield_config, "internal dipole")) {

    /* computing cell centered vector potential components */
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL x_local = p.x - dipole_x[0];
          CCTK_REAL y_local = p.y - dipole_y[0];
          CCTK_REAL Pcut = press_max * press_cut;
          CCTK_REAL Pdiff = std::max(press(p.I) - Pcut, 0.0);
          CCTK_REAL Aphi_local = Ab * pow(Pdiff, Avec_kappa);
          Avec_x_cent(p.I) = -y_local * Aphi_local;
          Avec_y_cent(p.I) = x_local * Aphi_local;
          Avec_z_cent(p.I) = 0.0;
        });

  } else if (CCTK_EQUALS(Afield_config, "external dipole")) {

    /* computing cell centered vector potential components */
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL x_local = p.x - dipole_x[0];
          CCTK_REAL y_local = p.y - dipole_y[0];
          CCTK_REAL z_local = p.z - dipole_z[0];
          CCTK_REAL cylrad2 = x_local * x_local + y_local * y_local;
          CCTK_REAL rsph =
              sqrt(x_local * x_local + y_local * y_local + z_local * z_local);
          CCTK_REAL rsph3 = pow(rsph, 3.0);
          CCTK_REAL r03 = pow(r0, 3.0);
          CCTK_REAL Aphi_local = A0 * (r03 / (r03 + rsph3));
          Avec_x_cent(p.I) = -y_local * Aphi_local;
          Avec_y_cent(p.I) = x_local * Aphi_local;
          Avec_z_cent(p.I) = 0.0;
        });

  } else if (CCTK_EQUALS(Afield_config, "random")) {

    // Init RNG
    srand (75447);
    CCTK_REAL* as = (CCTK_REAL*)The_Managed_Arena()->alloc(sizeof(CCTK_REAL) * 3 * nk * nk * nk);
    CCTK_REAL* bs = (CCTK_REAL*)The_Managed_Arena()->alloc(sizeof(CCTK_REAL) * 3 * nk * nk * nk);
    CCTK_REAL* cs = (CCTK_REAL*)The_Managed_Arena()->alloc(sizeof(CCTK_REAL) * 3 * nk * nk * nk);
    CCTK_REAL* ds = (CCTK_REAL*)The_Managed_Arena()->alloc(sizeof(CCTK_REAL) * 3 * nk * nk * nk);

    const int npts = nk * nk * nk;
    for (int ll=0; ll<nk; ll++) {
      for (int mm=0; mm<nk; mm++) {
        for (int nn=0; nn<nk; nn++) {
          const int lmn = nn + nk * mm + nk * nk * ll;
          CCTK_REAL ax_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 
          CCTK_REAL bx_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 
          CCTK_REAL cx_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 
          CCTK_REAL dx_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 

          CCTK_REAL ay_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 
          CCTK_REAL by_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 
          CCTK_REAL cy_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 
          CCTK_REAL dy_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 

          CCTK_REAL az_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 
          CCTK_REAL bz_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 
          CCTK_REAL cz_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 
          CCTK_REAL dz_lmn = ((CCTK_REAL)rand()) / ((CCTK_REAL)RAND_MAX); 

          as[         lmn] = ax_lmn;
          as[  npts + lmn] = ay_lmn;
          as[2*npts + lmn] = az_lmn;
          bs[         lmn] = bx_lmn;
          bs[  npts + lmn] = by_lmn;
          bs[2*npts + lmn] = bz_lmn;
          cs[         lmn] = cx_lmn;
          cs[  npts + lmn] = cy_lmn;
          cs[2*npts + lmn] = cz_lmn;
          ds[         lmn] = dx_lmn;
          ds[  npts + lmn] = dy_lmn;
          ds[2*npts + lmn] = dz_lmn;
        }
      }
    }

    const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

    CCTK_REAL PI = 4.0 * atan(1.0);

    /* computing cell centered vector potential components */
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL x_local = p.x - dipole_x[0];
          CCTK_REAL y_local = p.y - dipole_y[0];
          CCTK_REAL z_local = p.z - dipole_z[0];
          CCTK_REAL rsph =
              sqrt(x_local * x_local + y_local * y_local + z_local * z_local);
          CCTK_REAL rsph3 = pow(rsph, 3.0);
          CCTK_REAL r03 = pow(r0, 3.0);

          const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
            return calc_avg_v2c(gf_g(i, j), p);
          });
          const CCTK_REAL detg = calc_det(g_avg);
          const CCTK_REAL sqrt_detg = sqrt(detg);

          const CCTK_REAL scale = A0 * sqrt_detg * (r03 / (r03 + rsph3));

          // Generate random Avec
          CCTK_REAL Ax_modes = 0;
          CCTK_REAL Ay_modes = 0;
          CCTK_REAL Az_modes = 0;
          for (int ll=0; ll<nk; ll++) {
            for (int mm=0; mm<nk; mm++) {
              for (int nn=0; nn<nk; nn++) {
                const int lmn = nn + nk * mm + nk * nk * ll;
                CCTK_REAL ax_lmn = as[           lmn];
                CCTK_REAL ay_lmn = as[    npts + lmn];
                CCTK_REAL az_lmn = as[2 * npts + lmn];
                CCTK_REAL bx_lmn = bs[           lmn];
                CCTK_REAL by_lmn = bs[    npts + lmn];
                CCTK_REAL bz_lmn = bs[2 * npts + lmn];
                CCTK_REAL cx_lmn = cs[           lmn];
                CCTK_REAL cy_lmn = cs[    npts + lmn];
                CCTK_REAL cz_lmn = cs[2 * npts + lmn];
                CCTK_REAL dx_lmn = ds[           lmn];
                CCTK_REAL dy_lmn = ds[    npts + lmn];
                CCTK_REAL dz_lmn = ds[2 * npts + lmn];

                CCTK_REAL kl = ll == 0 ? 0.0 : 2 * PI / (lambda_min + (ll - 1) * delta_lambda);
                CCTK_REAL km = mm == 0 ? 0.0 : 2 * PI / (lambda_min + (mm - 1) * delta_lambda);
                CCTK_REAL kn = nn == 0 ? 0.0 : 2 * PI / (lambda_min + (nn - 1) * delta_lambda);

                CCTK_REAL X_k = p.x * kl + p.y * km + p.z * kn;
                Ax_modes += ax_lmn * cos(X_k + 2 * PI * bx_lmn) + cx_lmn * sin(X_k + 2 * PI * dx_lmn);
                Ay_modes += ay_lmn * cos(X_k + 2 * PI * by_lmn) + cy_lmn * sin(X_k + 2 * PI * dy_lmn);
                Az_modes += az_lmn * cos(X_k + 2 * PI * bz_lmn) + cz_lmn * sin(X_k + 2 * PI * dz_lmn);
              }
            }
          }
          Avec_x_cent(p.I) = scale * Ax_modes;
          Avec_y_cent(p.I) = scale * Ay_modes;
          Avec_z_cent(p.I) = scale * Az_modes;
        });

    The_Managed_Arena()->free(as);
    The_Managed_Arena()->free(bs);
    The_Managed_Arena()->free(cs);
    The_Managed_Arena()->free(ds);

  } else {
    CCTK_ERROR("Vector potential configuration not defined");
  }
}

extern "C" void AsterSeeds_InitializeStagAvec_TOV(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterSeeds_InitializeStagAvec_TOV;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<1, 0, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               Avec_x(p.I) = calc_avg_c2e(Avec_x_cent, p, 0);
                             });

  grid.loop_int<0, 1, 0>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               Avec_y(p.I) = calc_avg_c2e(Avec_y_cent, p, 1);
                             });

  grid.loop_int<0, 0, 1>(grid.nghostzones,
                         [=] CCTK_HOST(const Loop::PointDesc &p)
                             CCTK_ATTRIBUTE_ALWAYS_INLINE {
                               Avec_z(p.I) = calc_avg_c2e(Avec_z_cent, p, 2);
                             });
}

} // namespace AsterSeeds
