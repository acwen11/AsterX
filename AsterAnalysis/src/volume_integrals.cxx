#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <limits>
#include <cassert>

#include "aster_utils.hxx"
#include "../../../CarpetX/CarpetX/src/driver.hxx"
#include "../../../CarpetX/CarpetX/src/reduction.hxx"

namespace AsterAnalysis {
using namespace Loop;
using namespace AsterUtils;

static CCTK_REAL reduce_sum_carpetx(const char *varname, int tl = 0) {
  const int varindex = CCTK_VarIndex(varname);
  assert(varindex >= 0 && "Unknown variable name for reduction");
  const int gi = CCTK_GroupIndexFromVarI(varindex);
  assert(gi >= 0);
  const int v0 = CCTK_FirstVarIndexI(gi);
  assert(v0 >= 0);
  const int vi = varindex - v0;
  const CarpetX::reduction<CCTK_REAL, 3> red = CarpetX::reduce(gi, vi, tl);
  return red.sum; // includes dx*dy*dz
}

extern "C" void AsterAnalysis_VolumeIntegrals(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterAnalysis_VolumeIntegrals;
  DECLARE_CCTK_PARAMETERS;

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // metric & proper-volume weight
        const smat<CCTK_REAL, 3> g([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });
        const CCTK_REAL sqrt_detg = std::sqrt(calc_det(g));

        // v^2
        const CCTK_REAL Wlor = w_lorentz(p.I);
        const CCTK_REAL v2 = ((Wlor - 1.0) * (Wlor + 1.0)) / (Wlor * Wlor);

        const CCTK_REAL Bnorm = B_norm(p.I);
        const CCTK_REAL B2big = Bnorm * Bnorm; // B^j B_j
        const CCTK_REAL b2 = b2small(p.I);     // b^μ b_μ

        // weights: W * √γ for “fluid-weighted” integrals; √γ otherwise
        const CCTK_REAL Wcell = sqrt_detg;
        const CCTK_REAL Wcell_fluid = Wlor * sqrt_detg;

        // ---- EM energy integrand: (B^2 - 0.5 b^2) * √γ
        EMenergy_temps(p.I) = (B2big - 0.5 * b2) * Wcell;

        // ---- Kinetic energy: 0.5*(ρ + ρ ε + p) v^2 * W * √γ
        const CCTK_REAL rhoH = rho(p.I) + rho(p.I) * eps(p.I) + press(p.I);
        kinetic_energy_temps(p.I) = 0.5 * rhoH * v2 * Wcell_fluid;

        // ---- ⟨|B|⟩
        mean_B_numerator_temps(p.I) = rho(p.I) * Bnorm * Wcell_fluid;
        mean_B_denominator_temps(p.I) = rho(p.I) * Wcell_fluid;

        // ---- ⟨|B_tor|⟩, ⟨|B_pol|⟩ using prior spherical analysis
        mean_B_Tor_numerator_temps(p.I) =
            rho(p.I) * B_Toroidal_norm(p.I) * Wcell_fluid;
        mean_B_Pol_numerator_temps(p.I) =
            rho(p.I) * B_Poloidal_norm(p.I) * Wcell_fluid;
        mean_B_Tor_denominator_temps(p.I) = rho(p.I) * Wcell_fluid;
      });

  // ---- global reductions ----
  const CCTK_REAL EM_sum = reduce_sum_carpetx("AsterAnalysis::EMenergy_temps");
  const CCTK_REAL KE_sum =
      reduce_sum_carpetx("AsterAnalysis::kinetic_energy_temps");
  const CCTK_REAL meanB_num =
      reduce_sum_carpetx("AsterAnalysis::mean_B_numerator_temps");
  const CCTK_REAL meanB_den =
      reduce_sum_carpetx("AsterAnalysis::mean_B_denominator_temps");
  const CCTK_REAL tor_num =
      reduce_sum_carpetx("AsterAnalysis::mean_B_Tor_numerator_temps");
  const CCTK_REAL pol_num =
      reduce_sum_carpetx("AsterAnalysis::mean_B_Pol_numerator_temps");
  const CCTK_REAL tor_den =
      reduce_sum_carpetx("AsterAnalysis::mean_B_Tor_denominator_temps");

  const CCTK_REAL tiny = 1e-300;

  // ---- scalars ----
  *total_magnetic_energy = EM_sum;
  *total_kinetic_energy = KE_sum;
  *mean_B = meanB_num / (meanB_den + tiny);
  *mean_B_Tor = tor_num / (tor_den + tiny);
  *mean_B_Pol = pol_num / (tor_den + tiny);

  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "EM=%.6e  KE=%.6e  <|B|>=%.6e  <|B_tor|>=%.6e  <|B_pol|>=%.6e",
               *total_magnetic_energy, *total_kinetic_energy,
               *mean_B, *mean_B_Tor, *mean_B_Pol);
  }
}
} // namespace AsterAnalysis
