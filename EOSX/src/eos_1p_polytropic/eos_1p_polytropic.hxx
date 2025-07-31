#ifndef EOS_1P_POLYTROPIC_HXX
#define EOS_1P_POLYTROPIC_HXX

#include "../eos_1p.hxx"
#include <cmath>

namespace EOSX {

/// Polytropic EOS
/**
The polytropic EOS is written (using units with \f$ c=1 \f$) as
\f[ P = \rho_p \left(\frac{\rho}{\rho_p}\right)^\Gamma ,\qquad
\Gamma= 1+\frac{1}{n}  \f]
Here we use a <em>polytropic density scale</em> \f$ \rho_p \f$ to specify
the EOS instead of the usual form \f$ P = K \rho^\Gamma \f$
because it has simpler units than
\f$ K = \rho_p^{-1/n} \f$.

See eos_cold for notation used and eos_cold_api for a description of the member
functions.
*/

class eos_1p_polytropic : public eos_1p {
  CCTK_REAL poly_k;     ///< Polytropic density scale \f$ K \f$
  CCTK_REAL poly_gamma; ///< Polytropic exponent \f$ \Gamma \f$
  CCTK_REAL n;          ///< Polytropic index \f$ n = 1/(\Gamma-1) \f$
  CCTK_REAL np1;        ///< \f$ n+1 \f$
  CCTK_REAL invn;       ///< \f$ \frac{1}{n} \f$
  CCTK_REAL rho_p;      ///< Weird polytropic density scal \f$ rho_p \f$

public:
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  init(CCTK_REAL poly_gamma_, ///< Adiabatic index \f$ n \f$
       CCTK_REAL poly_k_,     ///< Density scale \f$ K \f$
       CCTK_REAL rho_max_     ///< Max valid density
  );

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  gm1_from_valid_rho(const CCTK_REAL rho) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  gm1_from_valid_p(const CCTK_REAL p) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  sed_from_valid_gm1(const CCTK_REAL gm1) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  ied_from_valid_gm1(const CCTK_REAL gm1) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  p_from_valid_gm1(const CCTK_REAL gm1) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  rho_from_valid_gm1(const CCTK_REAL gm1) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  hm1_from_valid_gm1(const CCTK_REAL gm1) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd2_from_valid_gm1(const CCTK_REAL gm1) const;
};

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
eos_1p_polytropic::init(CCTK_REAL poly_gamma_, CCTK_REAL poly_k_,
                        CCTK_REAL rho_max_) {
  poly_k = poly_k_;
  poly_gamma = poly_gamma_;
  assert(poly_gamma_ > 1.0);
  n = 1.0 / (poly_gamma - 1.0);
  np1 = n + 1;
  invn = 1.0 / n;
  rho_p = pow(poly_k, -n);
  // set_ranges(range(0, rho_max_));
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_1p_polytropic::gm1_from_valid_rho(const CCTK_REAL rho) const {
  return np1 * pow(rho / rho_p, invn);
}

/**
\return \f$ g-1 = h-1 = (n+1) \left(\frac{P}{\rho_p} \right)^\frac{1}{n+1} \f$
*/
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_1p_polytropic::gm1_from_valid_p(const CCTK_REAL p) const {
  return np1 * pow(p / rho_p, 1.0 / np1);
}

/**
\return Specific internal energy \f$ \epsilon = \frac{g-1}{\Gamma} \f$
*/
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_1p_polytropic::sed_from_valid_gm1(const CCTK_REAL gm1 ///< \f$ g-1 \f$
                                      ) const {
  return gm1 / poly_gamma;
}

/**
\return Internal energy density
\f$ \rho_I = n \rho_p \left(\frac{g-1}{n+1}\right)^{n+1} \f$
*/
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_1p_polytropic::ied_from_valid_gm1(const CCTK_REAL gm1 ///< \f$ g-1 \f$
                                      ) const {
  // return sed_from_gm1(gm1)*rho_from_gm1(gm1);
  return n * rho_p * pow(gm1 / np1, np1);
}

/**
\return Pressure \f$ P = \rho_p \left( \frac{g-1}{1+n} \right)^{1+n} \f$
*/
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_1p_polytropic::p_from_valid_gm1(const CCTK_REAL gm1 ///< \f$ g-1 \f$
                                    ) const {
  return rho_p * pow(gm1 / np1, np1);
}

/**
\return Rest mass density \f$ \rho = \rho_p \left( \frac{g-1}{1+n} \right)^n \f$
*/
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_1p_polytropic::rho_from_valid_gm1(const CCTK_REAL gm1 ///< \f$ g-1 \f$
                                      ) const {
  return rho_p * pow(gm1 / np1, n);
}

/**
\return specific enthalpy \f$ h-1 = g-1 \f$
*/
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_1p_polytropic::hm1_from_valid_gm1(const CCTK_REAL gm1 ///< \f$ g-1 \f$
                                      ) const {
  return gm1;
}

/**
\return Soundspeed squared \f$ c_s^2 = \frac{g-1}{ng} \f$
*/
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eos_1p_polytropic::csnd2_from_valid_gm1(const CCTK_REAL gm1) const {
  return gm1 / (n * (gm1 + 1));
}

} // namespace EOSX

#endif
