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
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomxs{fxmomx, fymomx, fzmomx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomys{fxmomy, fymomy, fzmomy};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomzs{fxmomz, fymomz, fzmomz};
  const vec<GF3D2<CCTK_REAL>, dim> fluxtaus{fxtau, fytau, fztau};
  const vec<GF3D2<CCTK_REAL>, dim> fluxDYes{fxDYe, fyDYe, fzDYe};
  const vec<GF3D2<CCTK_REAL>, dim> fluxDEnts{fxDEnt, fyDEnt, fzDEnt};

  const vec<GF3D2<CCTK_REAL>, dim> fluxdenss_HO{fxdensHO, fydensHO, fzdensHO};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomxs_HO{fxmomxHO, fymomxHO, fzmomxHO};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomys_HO{fxmomyHO, fymomyHO, fzmomyHO};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomzs_HO{fxmomzHO, fymomzHO, fzmomzHO};
  const vec<GF3D2<CCTK_REAL>, dim> fluxtaus_HO{fxtauHO, fytauHO, fztauHO};
  const vec<GF3D2<CCTK_REAL>, dim> fluxDYes_HO{fxDYeHO, fyDYeHO, fzDYeHO};
  const vec<GF3D2<CCTK_REAL>, dim> fluxDEnts_HO{fxDEntHO, fyDEntHO, fzDEntHO};
 
  // Face-centred grid functions (in direction `dir`)
  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};

  grid.loop_int_device<
      face_centred[0], face_centred[1],
      face_centred
          [2]>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    // Update correction grid functions:
    const int ordL = (LOflag(p.I) || LOflag(p.I - p.DI[dir])) ? 2 : correction_order;

    fluxdenss_HO(dir)(p.I) = higher_order_correction<dir>(fluxdenss(dir), p, ordL);
    fluxmomxs_HO(dir)(p.I) = higher_order_correction<dir>(fluxmomxs(dir), p, ordL);
    fluxmomys_HO(dir)(p.I) = higher_order_correction<dir>(fluxmomys(dir), p, ordL);
    fluxmomzs_HO(dir)(p.I) = higher_order_correction<dir>(fluxmomzs(dir), p, ordL);
    fluxtaus_HO(dir)(p.I) = higher_order_correction<dir>(fluxtaus(dir), p, ordL);
    fluxDYes_HO(dir)(p.I) = higher_order_correction<dir>(fluxDYes(dir), p, ordL);
    fluxDEnts_HO(dir)(p.I) = higher_order_correction<dir>(fluxDEnts(dir), p, ordL);

  });

  // Update flux grid functions with computed corrections
  grid.loop_int_device<
      face_centred[0], face_centred[1],
      face_centred
          [2]>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

    fluxdenss(dir)(p.I) = fluxdenss_HO(dir)(p.I);
    fluxmomxs(dir)(p.I) = fluxmomxs_HO(dir)(p.I);
    fluxmomys(dir)(p.I) = fluxmomys_HO(dir)(p.I);
    fluxmomzs(dir)(p.I) = fluxmomzs_HO(dir)(p.I);
    fluxtaus(dir)(p.I) = fluxtaus_HO(dir)(p.I);
    fluxDYes(dir)(p.I) = fluxDYes_HO(dir)(p.I);
    fluxDEnts(dir)(p.I) = fluxDEnts_HO(dir)(p.I);

  });

}

extern "C" void AsterX_Flux_Corrections(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Flux_Corrections;
  DECLARE_CCTK_PARAMETERS;

  CalcFluxCorrections<0>(cctkGH);
  CalcFluxCorrections<1>(cctkGH);
  CalcFluxCorrections<2>(cctkGH);
}

} // namespace AsterX
