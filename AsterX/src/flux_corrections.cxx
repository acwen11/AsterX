#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "aster_utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;

template <int dir> void CalcFluxCorrections(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Flux_Corrections;
  DECLARE_CCTK_PARAMETERS;

  /* grid functions for fluxes */
  const vec<GF3D2<CCTK_REAL>, dim> fluxdenss{fxdens, fydens, fzdens};
  const vec<GF3D2<CCTK_REAL>, dim> fluxDEnts{fxDEnt, fyDEnt, fzDEnt};
  const vec<GF3D2<CCTK_REAL>, dim> fluxDYes{fxDYe, fyDYe, fzDYe};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomxs{fxmomx, fymomx, fzmomx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomys{fxmomy, fymomy, fzmomy};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomzs{fxmomz, fymomz, fzmomz};
  const vec<GF3D2<CCTK_REAL>, dim> fluxtaus{fxtau, fytau, fztau};
  // const vec<GF3D2<CCTK_REAL>, dim> fluxBxs{fxBx, fyBx, fzBx};
  // const vec<GF3D2<CCTK_REAL>, dim> fluxBys{fxBy, fyBy, fzBy};
  // const vec<GF3D2<CCTK_REAL>, dim> fluxBzs{fxBz, fyBz, fzBz};

  // const vec<GF3D2<CCTK_REAL>, dim> fluxdenss_HO{fHOxdens, fHOydens, fHOzdens};
  // const vec<GF3D2<CCTK_REAL>, dim> fluxDEnts_HO{fHOxDEnt, fHOyDEnt, fHOzDEnt};
  // const vec<GF3D2<CCTK_REAL>, dim> fluxmomxs_HO{fHOxmomx, fHOymomx, fHOzmomx};
  // const vec<GF3D2<CCTK_REAL>, dim> fluxmomys_HO{fHOxmomy, fHOymomy, fHOzmomy};
  // const vec<GF3D2<CCTK_REAL>, dim> fluxmomzs_HO{fHOxmomz, fHOymomz, fHOzmomz};
  // const vec<GF3D2<CCTK_REAL>, dim> fluxtaus_HO{fHOxtau, fHOytau, fHOztau};

  // Face-centred grid functions (in direction `dir`)
  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};

  grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Update correction grid functions:

        if ((correction_order != 2) && (correction_order != 4) &&
            (correction_order != 6)) {
          printf("Incorrect correction order for the fluxes! \n");
          assert(0);
        }

        fluxdenss(dir)(p.I) =
            higher_order_correction(fluxdenss(dir), p, dir, correction_order);
        fluxDEnts(dir)(p.I) =
            higher_order_correction(fluxDEnts(dir), p, dir, correction_order);
        fluxDYes(dir)(p.I) =
            higher_order_correction(fluxDYes(dir), p, dir, correction_order);
        fluxmomxs(dir)(p.I) =
            higher_order_correction(fluxmomxs(dir), p, dir, correction_order);
        fluxmomys(dir)(p.I) =
            higher_order_correction(fluxmomys(dir), p, dir, correction_order);
        fluxmomzs(dir)(p.I) =
            higher_order_correction(fluxmomzs(dir), p, dir, correction_order);
        fluxtaus(dir)(p.I) =
            higher_order_correction(fluxtaus(dir), p, dir, correction_order);

        if (isnan(fluxdenss(dir)(p.I)) || isnan(fluxDEnts(dir)(p.I)) ||
            isnan(fluxmomxs(dir)(p.I)) || isnan(fluxmomys(dir)(p.I)) ||
            isnan(fluxmomzs(dir)(p.I)) || isnan(fluxtaus(dir)(p.I))) {
          printf("NaN in HO Flux corr. at cctk_iteration = %i,  dir = %i,  ijk = %i, %i, %i, "
                 "x, y, z = %16.8e, %16.8e, %16.8e.\n",
                 cctk_iteration, dir, p.i, p.j, p.k, p.x, p.y, p.z);
          printf("  fluxHOdenss = %16.8e,\n", fluxdenss(dir)(p.I));
          printf("  fluxHODEnts = %16.8e,\n", fluxDEnts(dir)(p.I));
          printf("  fluxHOmoms  = %16.8e, %16.8e, %16.8e,\n",
                 fluxmomxs(dir)(p.I), fluxmomys(dir)(p.I),
                 fluxmomzs(dir)(p.I));
          printf("  fluxHOtaus  = %16.8e,\n", fluxtaus(dir)(p.I));
          assert(0);
        }
      });

  // Update flux grid functions with computed corrections
  // grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
  //     grid.nghostzones,
  //     [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
  //       fluxdenss(dir)(p.I) = fluxdenss_HO(dir)(p.I);
  //       fluxDEnts(dir)(p.I) = fluxDEnts_HO(dir)(p.I);
  //       fluxmomxs(dir)(p.I) = fluxmomxs_HO(dir)(p.I);
  //       fluxmomys(dir)(p.I) = fluxmomys_HO(dir)(p.I);
  //       fluxmomzs(dir)(p.I) = fluxmomzs_HO(dir)(p.I);
  //       fluxtaus(dir)(p.I) = fluxtaus_HO(dir)(p.I);
  //     });
}

extern "C" void AsterX_Flux_Corrections(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Flux_Corrections;
  DECLARE_CCTK_PARAMETERS;

  CalcFluxCorrections<0>(cctkGH);
  CalcFluxCorrections<1>(cctkGH);
  CalcFluxCorrections<2>(cctkGH);
}

} // namespace AsterX
