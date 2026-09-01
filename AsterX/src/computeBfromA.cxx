#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

#include "aster_utils.hxx"

#include "../../../CarpetX/CarpetX/src/schedule.hxx"

struct metric {
  CCTK_REAL gxx, gxy, gxz, gyy, gyz, gzz;
};

namespace AsterX {
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;

template <int dir> void ComputeStaggeredFaceAvgB(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputeAvgdBstagFromAvgA;
  DECLARE_CCTK_PARAMETERS;

  static_assert(dir >= 0 && dir < 3, "");

  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};
  grid.loop_all_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (dir == 0) {
          /* <dB>x is curl(<A>) at (i-1/2,j,k) */
          dBx_stag_fa(p.I) =
              calc_fd_forward_midpoint<1>(Avec_z, p, 2) -
              calc_fd_forward_midpoint<2>(Avec_y, p, 2);
          // Also init PV, will be overwritten
          dBx_stag(p.I) = 0.0;
        } else if (dir == 1) {
          /* <dB>y is curl(<A>) at (i,j-1/2,k) */
          dBy_stag_fa(p.I) =
              calc_fd_forward_midpoint<2>(Avec_x, p, 2) -
              calc_fd_forward_midpoint<0>(Avec_z, p, 2);
          dBy_stag(p.I) = 0.0;
        } else if (dir == 2) {
          /* <dB>z is curl(<A>) at (i,j,z-1/2) */
          dBz_stag_fa(p.I) =
              calc_fd_forward_midpoint<0>(Avec_y, p, 2) -
              calc_fd_forward_midpoint<1>(Avec_x, p, 2);
          dBz_stag(p.I) = 0.0;
        }
      });
}

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

template <int dir> void SetdBstagnMinus1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetdBstagnMinus1;
  DECLARE_CCTK_PARAMETERS;

  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};

  if (*dBstag_pv_iter == n_dBstagpv_iters) {
    grid.loop_all_device<face_centred[0], face_centred[1], face_centred[2]>(
       grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          // Temporarily use flux GF as helper
          if (dir == 0) {
            fxdens(p.I) = dBx_stag_fa(p.I);
          } else if (dir == 1) {
            fydens(p.I) = dBy_stag_fa(p.I);
          } else if (dir == 2) {
            fzdens(p.I) = dBz_stag_fa(p.I);
          }
        });
  } else {
    grid.loop_all_device<face_centred[0], face_centred[1], face_centred[2]>(
       grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (dir == 0) {
            fxdens(p.I) = dBx_stag(p.I);
          } else if (dir == 1) {
            fydens(p.I) = dBy_stag(p.I);
          } else if (dir == 2) {
            fzdens(p.I) = dBz_stag(p.I);
          }
        });
  }
}

template <int dir> void ComputeStaggeredPointValB(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBstagIter;
  DECLARE_CCTK_PARAMETERS;

  static_assert(dir >= 0 && dir < 3, "");

  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};

  // Here, we have already computed the face averaged <dBi_stag> from <Avec>.
  // Now, we use Equation 21 of https://arxiv.org/pdf/2310.11831 to compute the point value
  // dBi_stag.
  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1)/CCTK_REAL(24);
  grid.loop_allmn_device<face_centred[0], face_centred[1], face_centred[2]>(
     grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Temporarily use flux GF as helper
        if (dir == 0) {
          dBx_stag(p.I) = dBx_stag_fa(p.I) - one_over_24 * laplace_perp<0>(fxdens, p);
        } else if (dir == 1) {
          dBy_stag(p.I) = dBy_stag_fa(p.I) - one_over_24 * laplace_perp<1>(fydens, p);
        } else if (dir == 2) {
          dBz_stag(p.I) = dBz_stag_fa(p.I) - one_over_24 * laplace_perp<2>(fzdens, p);
        }
      });

  // Calculate Bstag in boundaries/ghosts with lower order
  grid.loop_outer_n_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones, 1,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (dir == 0) {
          dBx_stag(p.I) = dBx_stag_fa(p.I);
        } else if (dir == 1) {
          dBy_stag(p.I) = dBy_stag_fa(p.I);
        } else if (dir == 2) {
          dBz_stag(p.I) = dBz_stag_fa(p.I);
        }
      });
}

template <int dir> void ComputeStaggeredPointValB_Int(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBstagIter;
  DECLARE_CCTK_PARAMETERS;

  static_assert(dir >= 0 && dir < 3, "");

  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};

  // Here, we have already computed the face averaged <dBi_stag> from <Avec>.
  // Now, we use Equation 21 of https://arxiv.org/pdf/2310.11831 to compute the point value
  // dBi_stag.
  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1)/CCTK_REAL(24);
  grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
     grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Temporarily use flux GF as helper
        if (dir == 0) {
          dBx_stag(p.I) = dBx_stag_fa(p.I) - one_over_24 * laplace_perp<0>(fxdens, p);
        } else if (dir == 1) {
          dBy_stag(p.I) = dBy_stag_fa(p.I) - one_over_24 * laplace_perp<1>(fydens, p);
        } else if (dir == 2) {
          dBz_stag(p.I) = dBz_stag_fa(p.I) - one_over_24 * laplace_perp<2>(fzdens, p);
        }
      });
}

extern "C" void AsterX_ComputedBstagFromA(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBstagFromA;
  DECLARE_CCTK_PARAMETERS;

  ComputeStaggeredB<0>(cctkGH);
  ComputeStaggeredB<1>(cctkGH);
  ComputeStaggeredB<2>(cctkGH);
}

/* Point Value dBstag Calculation */
extern "C" void AsterX_ComputeAvgdBstagFromAvgA(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputeAvgdBstagFromAvgA;
  DECLARE_CCTK_PARAMETERS;

  ComputeStaggeredFaceAvgB<0>(cctkGH);
  ComputeStaggeredFaceAvgB<1>(cctkGH);
  ComputeStaggeredFaceAvgB<2>(cctkGH);
}

extern "C" void AsterX_SetdBstagIter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetdBstagIter;
  DECLARE_CCTK_PARAMETERS;

  *dBstag_pv_iter = n_dBstagpv_iters;
}

extern "C" void AsterX_SetdBstagnMinus1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetdBstagnMinus1;
  DECLARE_CCTK_PARAMETERS;

  SetdBstagnMinus1<0>(cctkGH);
  SetdBstagnMinus1<1>(cctkGH);
  SetdBstagnMinus1<2>(cctkGH);
}

extern "C" void AsterX_ComputedBstagIter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBstagIter;
  DECLARE_CCTK_PARAMETERS;

  if (n_dBstagpv_iters - *dBstag_pv_iter >= dBstag_iter_loopswitch) {
    ComputeStaggeredPointValB_Int<0>(cctkGH);
    ComputeStaggeredPointValB_Int<1>(cctkGH);
    ComputeStaggeredPointValB_Int<2>(cctkGH);
  }
  else {
    ComputeStaggeredPointValB<0>(cctkGH);
    ComputeStaggeredPointValB<1>(cctkGH);
    ComputeStaggeredPointValB<2>(cctkGH);
  }
}

extern "C" void AsterX_DecdBstagIter(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_DecdBstagIter;
  DECLARE_CCTK_PARAMETERS;

  *dBstag_pv_iter -= 1;

  int minghosts = min(cctk_nghostzones[0], min(cctk_nghostzones[1], cctk_nghostzones[2]));
  if (*dBstag_pv_iter==0 || ((n_dBstagpv_iters - *dBstag_pv_iter) % minghosts == 0))  {
    static const std::vector<int> groups = {CCTK_GroupIndex("AsterX::dBx_stag"),
                                        CCTK_GroupIndex("AsterX::dBy_stag"),
                                        CCTK_GroupIndex("AsterX::dBz_stag")};

    SyncGroupsByDirIGhostOnly(cctkGH, groups.size(), groups.data(), nullptr);
  }
}
/* End Point Value dBstag Calculation */

// Fall back to 2nd order interpolation from faces to cell centers if requested
extern "C" void AsterX_ComputedBFromdBstag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_ComputedBFromdBstag;
  DECLARE_CCTK_PARAMETERS;

  const int interp_order = use_ho_fv ? 4 : mag_correction_order;
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* Interpolation of staggered B components to cell center
         */
        const int ordL = ((LOflag(p.I) > 0.0) && shock_Bstag_fallback) ? 2 : interp_order;
        dBx(p.I) = calc_avg_f2c(dBx_stag, p, 0, ordL);
        dBy(p.I) = calc_avg_f2c(dBy_stag, p, 1, ordL);
        dBz(p.I) = calc_avg_f2c(dBz_stag, p, 2, ordL);
      });

  // Interpolate dB in boundaries/ghosts at lower order
  // if (nloop != 0) {
  //   grid.loop_outer_n_device<1, 1, 1>(
  //       grid.nghostzones, nloop,
  //       [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
  //         /* Interpolation of staggered B components to cell center
  //          */
  //         dBx(p.I) = calc_avg_f2c(dBx_stag, p, 0, 2);
  //         dBy(p.I) = calc_avg_f2c(dBy_stag, p, 1, 2);
  //         dBz(p.I) = calc_avg_f2c(dBz_stag, p, 2, 2);
  //       });
  // }
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
