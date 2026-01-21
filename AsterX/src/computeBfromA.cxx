#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

#include "aster_utils.hxx"

struct metric {
  CCTK_REAL gxx, gxy, gxz, gyy, gyz, gzz;
};

namespace AsterX {
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;

template <int dir> void ComputeStaggeredB(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBstagFromA;
  DECLARE_CCTK_PARAMETERS;

  static_assert(dir >= 0 && dir < 3, "");

  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};
  const int nloop = (mag_correction_order - 2) / 2;
  grid.loop_allmn_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones, nloop,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (dir == 0) {
          /* dBx is curl(A) at (i-1/2,j,k) */
          dBx_stag(p.I) =
              calc_fd_forward_midpoint<1>(Avec_z, p, mag_correction_order) -
              calc_fd_forward_midpoint<2>(Avec_y, p, mag_correction_order);
        } else if (dir == 1) {
          /* dBy is curl(A) at (i,j-1/2,k) */
          dBy_stag(p.I) =
              calc_fd_forward_midpoint<2>(Avec_x, p, mag_correction_order) -
              calc_fd_forward_midpoint<0>(Avec_z, p, mag_correction_order);
        } else if (dir == 2) {
          /* dBz is curl(A) at (i,j,z-1/2) */
          dBz_stag(p.I) =
              calc_fd_forward_midpoint<0>(Avec_y, p, mag_correction_order) -
              calc_fd_forward_midpoint<1>(Avec_x, p, mag_correction_order);
        }

        // TODO: need to implement copy conditions?
      });

  // Calculate Bstag in boundaries/ghosts with lower order
  if (nloop != 0) {
    grid.loop_outer_n_device<face_centred[0], face_centred[1], face_centred[2]>(
        grid.nghostzones, nloop,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (dir == 0) {
            /* dBx is curl(A) at (i-1/2,j,k) */
            dBx_stag(p.I) = calc_fd_forward_midpoint<1>(Avec_z, p, 2) -
                            calc_fd_forward_midpoint<2>(Avec_y, p, 2);
          } else if (dir == 1) {
            /* dBy is curl(A) at (i,j-1/2,k) */
            dBy_stag(p.I) = calc_fd_forward_midpoint<2>(Avec_x, p, 2) -
                            calc_fd_forward_midpoint<0>(Avec_z, p, 2);
          } else if (dir == 2) {
            /* dBz is curl(A) at (i,j,z-1/2) */
            dBz_stag(p.I) = calc_fd_forward_midpoint<0>(Avec_y, p, 2) -
                            calc_fd_forward_midpoint<1>(Avec_x, p, 2);
          }
        });
  }
}

extern "C" void AsterX_ComputedBstagFromA(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBstagFromA;
  DECLARE_CCTK_PARAMETERS;

  ComputeStaggeredB<0>(cctkGH);
  ComputeStaggeredB<1>(cctkGH);
  ComputeStaggeredB<2>(cctkGH);
}

extern "C" void AsterX_ComputedBFromdBstag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBFromdBstag;
  DECLARE_CCTK_PARAMETERS;

  const int nloop = (mag_correction_order - 2) / 2;
  grid.loop_allmn_device<1, 1, 1>(
      grid.nghostzones, nloop,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* Interpolation of staggered B components to cell center
         */
        dBx(p.I) = calc_avg_f2c(dBx_stag, p, 0, mag_correction_order);
        dBy(p.I) = calc_avg_f2c(dBy_stag, p, 1, mag_correction_order);
        dBz(p.I) = calc_avg_f2c(dBz_stag, p, 2, mag_correction_order);
      });

  // Interpolate dB in boundaries/ghosts at lower order
  if (nloop != 0) {
    // Lower bound
    grid.loop_outer_n_device<1, 1, 1>(
        grid.nghostzones, nloop,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          /* Interpolation of staggered B components to cell center
           */
          dBx(p.I) = calc_avg_f2c(dBx_stag, p, 0, 2);
          dBy(p.I) = calc_avg_f2c(dBy_stag, p, 1, 2);
          dBz(p.I) = calc_avg_f2c(dBz_stag, p, 2, 2);
        });
  }
}

extern "C" void AsterX_ComputeBFromdB(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputeBFromdB;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* Interpolate metric terms from vertices to center */
        metric g;
        g.gxx = calc_avg_v2c(gxx, p);
        g.gxy = calc_avg_v2c(gxy, p);
        g.gxz = calc_avg_v2c(gxz, p);
        g.gyy = calc_avg_v2c(gyy, p);
        g.gyz = calc_avg_v2c(gyz, p);
        g.gzz = calc_avg_v2c(gzz, p);

        /* Determinant of spatial metric */
        const smat<CCTK_REAL, 3> gmat{g.gxx, g.gxy, g.gxz, g.gyy, g.gyz, g.gzz};
        const CCTK_REAL sqrt_detg = sqrt(calc_det(gmat));

        /* Second order interpolation of staggered B components to cell center
         */
        Bvecx(p.I) = dBx(p.I) / sqrt_detg;
        Bvecy(p.I) = dBy(p.I) / sqrt_detg;
        Bvecz(p.I) = dBz(p.I) / sqrt_detg;
      });
}

} // namespace AsterX
