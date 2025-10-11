#ifndef RECONX_WENOZP_HXX
#define RECONX_WENOZP_HXX

#include <cctk.h>

#include <loop_device.hxx>

#include <array>
#include <cmath>

namespace ReconX {

using std::array;

/**
 * WENO-Z: Performs the reconstruction of a given variable using 5th order
 * WENO-Z method based on Borges et al., "An improved weighted essentially
 * non-oscillatory scheme for hyperbolic conservation laws", 2008). Also, see
 * the Spritz code for details
 */
template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<T, 2>
wenozp(T gf_Imm, T gf_Im, T gf_I, T gf_Ip, T gf_Ipp, const CCTK_REAL dx, T weno_eps) {

  using Arith::vec;
  using std::abs;

  // Computing the smoothness indicators (Borges et al. 2008)
  const vec<T, 3> betaZ{
      (13.0 / 12.0) * pow2(gf_Imm - 2.0 * gf_Im + gf_I) +
          (1.0 / 4.0) * pow2(gf_Imm - 4.0 * gf_Im + 3.0 * gf_I),

      (13.0 / 12.0) * pow2(gf_Im - 2.0 * gf_I + gf_Ip) +
          (1.0 / 4.0) * pow2(gf_Im - gf_Ip),

      (13.0 / 12.0) * pow2(gf_I - 2.0 * gf_Ip + gf_Ipp) +
          (1.0 / 4.0) * pow2(3.0 * gf_I - 4.0 * gf_Ip + gf_Ipp)};

  // Defining tau5 based on eq. 25 of (Borges et al. 2008)
  const T tau5{abs(betaZ(0) - betaZ(2))};

  // Unnormalized weights according to Eq. 10 of Acker+2016 "An improved WENO-Z
  // scheme", Journal of Computational Physics, 313, 726–753. 

  vec<T, 3> aux_alphaZ;
  aux_alphaZ(0) = (tau5 + weno_eps) / (betaZ(0) + weno_eps);
  aux_alphaZ(1) = (tau5 + weno_eps) / (betaZ(1) + weno_eps);
  aux_alphaZ(2) = (tau5 + weno_eps) / (betaZ(2) + weno_eps);

  // Optimal weights according to Del Zanna et al. 2007
  const vec<T, 3> wt{5.0 / 16.0, 10.0 / 16.0, 1.0 / 16.0};

  // Original weights as suggested in (Borges et al. 2008)
  // const vec<T, 3> wt{3.0 / 10.0, 3.0 / 5.0, 1.0 / 10.0};

  vec<vec<T, 2>, 3> alphaZ;

  // for minus side
  alphaZ(0)(0) = wt(0) * (1.0 + pow2(aux_alphaZ(0)) +
                          pow(dx, 2.0 / 3.0) * (1.0 / aux_alphaZ(0)));
  alphaZ(1)(0) = wt(1) * (1.0 + pow2(aux_alphaZ(1)) +
                          pow(dx, 2.0 / 3.0) * (1.0 / aux_alphaZ(1)));
  alphaZ(2)(0) = wt(2) * (1.0 + pow2(aux_alphaZ(2)) +
                          pow(dx, 2.0 / 3.0) * (1.0 / aux_alphaZ(2)));

  // for plus side
  alphaZ(0)(1) = wt(2) * (1.0 + pow2(aux_alphaZ(0)) +
                          pow(dx, 2.0 / 3.0) * (1.0 / aux_alphaZ(0)));
  alphaZ(1)(1) = wt(1) * (1.0 + pow2(aux_alphaZ(1)) +
                          pow(dx, 2.0 / 3.0) * (1.0 / aux_alphaZ(1)));
  alphaZ(2)(1) = wt(0) * (1.0 + pow2(aux_alphaZ(2)) +
                          pow(dx, 2.0 / 3.0) * (1.0 / aux_alphaZ(2)));

  // Normalized weights for the reconstruction (Del Zanna et al. 2007)
  const vec<T, 2> omega_denom([&](int j) ARITH_INLINE {
    return alphaZ(0)(j) + alphaZ(1)(j) + alphaZ(2)(j);
  });

  const vec<vec<T, 2>, 3> omegaZ([&](int j) ARITH_INLINE {
    return vec<T, 2>(
        [&](int f) ARITH_INLINE { return 0.125 * alphaZ(j)(f) / omega_denom(f); });
  });

  // Reconstruct cell-centered variable to left (minus) and right (plus) cell
  // interfaces

  // Spritz Weights:
  const T var_m =
      omegaZ(2)(0) * (3.0 * gf_Ipp - 10.0 * gf_Ip + 15.0 * gf_I) +
      omegaZ(1)(0) * (-1.0 * gf_Ip + 6.0 * gf_I + 3.0 * gf_Im) +
      omegaZ(0)(0) * (3.0 * gf_I + 6.0 * gf_Im - 1.0 * gf_Imm);

  const T var_p =
      omegaZ(0)(1) * (3.0 * gf_Imm - 10.0 * gf_Im + 15.0 * gf_I) +
      omegaZ(1)(1) * (-1.0 * gf_Im + 6.0 * gf_I + 3.0 * gf_Ip) +
      omegaZ(2)(1) * (3.0 * gf_I + 6.0 * gf_Ip - 1.0 * gf_Ipp);

  /* GRHydro Weights:
  const T var_m{
      (omegaZ(2)(0) / 6.0) * (2.0 * gf_Ipp - 7.0 * gf_Ip + 11.0 * gf_I) +
      (omegaZ(1)(0) / 6.0) * (-1.0 * gf_Ip + 5.0 * gf_I + 2.0 * gf_Im) +
      (omegaZ(0)(0) / 6.0) * (2.0 * gf_I + 5.0 * gf_Im - 1.0 * gf_Imm)};

  const T var_p{
      (omegaZ(0)(1) / 6.0) * (2.0 * gf_Imm - 7.0 * gf_Im + 11.0 * gf_I) +
      (omegaZ(1)(1) / 6.0) * (-1.0 * gf_Im + 5.0 * gf_I + 2.0 * gf_Ip) +
      (omegaZ(2)(1) / 6.0) * (2.0 * gf_I + 5.0 * gf_Ip - 1.0 * gf_Ipp)};
  */

  return {var_m, var_p};
}

template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<T, 2>
wenozp_reconstruct(T gf_Immm, T gf_Imm, T gf_Im, T gf_Ip, T gf_Ipp, T gf_Ippp,
                  const CCTK_REAL dx, T weno_eps) {

  const auto rc_Im{wenozp(gf_Immm, gf_Imm, gf_Im, gf_Ip, gf_Ipp, dx, weno_eps)};
  const auto rc_Ip{wenozp(gf_Imm, gf_Im, gf_Ip, gf_Ipp, gf_Ippp, dx, weno_eps)};

  return {rc_Im[1], rc_Ip[0]};
}

} // namespace ReconX

#endif // RECONX_WENOZP_HXX
