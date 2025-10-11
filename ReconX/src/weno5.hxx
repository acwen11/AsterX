#ifndef RECONX_WENO5_HXX
#define RECONX_WENO5_HXX

#include <cctk.h>

#include <loop_device.hxx>

#include <array>
#include <cmath>

namespace ReconX {

using std::array;

/**
 * WENO5:
 */
template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<T, 2>
weno5(T gf_Imm, T gf_Im, T gf_I, T gf_Ip, T gf_Ipp, T weno_eps) {

  using Arith::vec;
  using std::abs;

	// GRHydro Weights
	const CCTK_REAL beta_shu[3][6] = { { 4.0/3.0,  -19.0/3.0, 25.0/3.0, 11.0/3.0, -31.0/3.0, 10.0/3.0 },
                                      { 4.0/3.0,  -13.0/3.0, 13.0/3.0, 5.0/3.0,  -13.0/3.0, 4.0/3.0 },
                                      { 10.0/3.0, -31.0/3.0, 25.0/3.0, 11.0/3.0, -19.0/3.0, 4.0/3.0 } };
  const CCTK_REAL weno_coeffs[3][5] = { { 3.0/8.0, -5.0/4.0, 15.0/8.0, 0,      0 },
                                         { 0,       -1.0/8.0, 3.0/4.0,  3.0/8.0, 0 },
                                         { 0,       0,        3.0/8.0,  3.0/4.0, -1.0/8.0 } };

/*
  constexpr array<array<CCTK_REAL, 6>, 3> beta_shu = {{ { 4.0/3.0,  -19.0/3.0, 25.0/3.0, 11.0/3.0, -31.0/3.0, 10.0/3.0 },
		{ 4.0/3.0,  -13.0/3.0, 13.0/3.0, 5.0/3.0,  -13.0/3.0, 4.0/3.0 },
		{ 10.0/3.0, -31.0/3.0, 25.0/3.0, 11.0/3.0, -19.0/3.0, 4.0/3.0 } }};

	constexpr array<array<CCTK_REAL, 5>, 3> weno_coeffs = {{ { 3.0/8.0, -5.0/4.0, 15.0/8.0, 0,      0 },
		{ 0,       -1.0/8.0, 3.0/4.0,  3.0/8.0, 0 },
		{ 0,       0,        3.0/8.0,  3.0/4.0, -1.0/8.0 } }};
*/

  const vec<T, 3> beta{
			beta_shu[0][0] * pow2(gf_Imm)
			+ beta_shu[0][1] * gf_Imm * gf_Im
			+ beta_shu[0][2] * pow2(gf_Im)
			+ beta_shu[0][3] * gf_Imm * gf_I
			+ beta_shu[0][4] * gf_Im * gf_I
			+ beta_shu[0][5] * pow2(gf_I),

	 		beta_shu[1][0] * pow2(gf_Im)
			+ beta_shu[1][1] * gf_Im * gf_I
			+ beta_shu[1][2] * pow2(gf_I)
			+ beta_shu[1][3] * gf_Im * gf_Ip
			+ beta_shu[1][4] * gf_I * gf_Ip
			+ beta_shu[1][5] * pow2(gf_Ip),

			beta_shu[2][0] * pow2(gf_I)
			+ beta_shu[2][1] * gf_I * gf_Ip
			+ beta_shu[2][2] * pow2(gf_Ip)
			+ beta_shu[2][3] * gf_I * gf_Ipp
			+ beta_shu[2][4] * gf_Ip * gf_Ipp
			+ beta_shu[2][5] * pow2(gf_Ipp) };

  vec<T, 3> wbarplus;
	wbarplus(0) = 1.0/16.0 / pow2(weno_eps + beta(0));
	wbarplus(1) = 5.0/8.0 / pow2(weno_eps + beta(1));
	wbarplus(2) = 5.0/16.0 / pow2(weno_eps + beta(2));

	const CCTK_REAL iwbarplussum = 1.0 / (wbarplus(0) + wbarplus(1) + wbarplus(2));

  vec<T, 3> wplus;
  wplus(0) = wbarplus(0) * iwbarplussum;
  wplus(1) = wbarplus(1) * iwbarplussum;
  wplus(2) = wbarplus(2) * iwbarplussum;

  vec<T, 3> wbarminus;
	wbarminus(0) = 5.0/16.0 / pow2(weno_eps + beta(0));
	wbarminus(1) = 5.0/8.0 / pow2(weno_eps + beta(1));
	wbarminus(2) = 1.0/16.0 / pow2(weno_eps + beta(2));

	const CCTK_REAL iwbarminussum = 1.0 / (wbarminus(0) + wbarminus(1) + wbarminus(2));

  vec<T, 3> wminus;
  wminus(0) = wbarminus(0) * iwbarminussum;
  wminus(1) = wbarminus(1) * iwbarminussum;
  wminus(2) = wbarminus(2) * iwbarminussum;

  // Reconstruct cell-centered variable to left (minus) and right (plus) cell
  // interfaces
  const T var_m =
		(wminus(0) * weno_coeffs[2][4] + wminus(1) * weno_coeffs[1][4] + wminus(2) * weno_coeffs[0][4]) * gf_Imm
		+ (wminus(0) * weno_coeffs[2][3] + wminus(1) * weno_coeffs[1][3] + wminus(2) * weno_coeffs[0][3]) * gf_Im
		+ (wminus(0) * weno_coeffs[2][2] + wminus(1) * weno_coeffs[1][2] + wminus(2) * weno_coeffs[0][2]) * gf_I
		+ (wminus(0) * weno_coeffs[2][1] + wminus(1) * weno_coeffs[1][1] + wminus(2) * weno_coeffs[0][1]) * gf_Ip
		+ (wminus(0) * weno_coeffs[2][0] + wminus(1) * weno_coeffs[1][0] + wminus(2) * weno_coeffs[0][0]) * gf_Ipp;

  const T var_p =
		(wplus(0) * weno_coeffs[0][0] + wplus(1) * weno_coeffs[1][0] + wplus(2) * weno_coeffs[2][0]) * gf_Imm
		+ (wplus(0) * weno_coeffs[0][1] + wplus(1) * weno_coeffs[1][1] + wplus(2) * weno_coeffs[2][1]) * gf_Im
		+ (wplus(0) * weno_coeffs[0][2] + wplus(1) * weno_coeffs[1][2] + wplus(2) * weno_coeffs[2][2]) * gf_I
		+ (wplus(0) * weno_coeffs[0][3] + wplus(1) * weno_coeffs[1][3] + wplus(2) * weno_coeffs[2][3]) * gf_Ip
		+ (wplus(0) * weno_coeffs[0][4] + wplus(1) * weno_coeffs[1][4] + wplus(2) * weno_coeffs[2][4]) * gf_Ipp;

  return {var_m, var_p};
}

template <typename T = CCTK_REAL>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<T, 2>
weno5_reconstruct(T gf_Immm, T gf_Imm, T gf_Im, T gf_Ip, T gf_Ipp, T gf_Ippp,
                  T weno_eps) {

  const auto rc_Im{weno5(gf_Immm, gf_Imm, gf_Im, gf_Ip, gf_Ipp, weno_eps)};
  const auto rc_Ip{weno5(gf_Imm, gf_Im, gf_Ip, gf_Ipp, gf_Ippp, weno_eps)};

  return {rc_Im[1], rc_Ip[0]};
}

} // namespace ReconX

#endif // RECONX_WENO5_HXX
