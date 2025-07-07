/*! \file eos_3p.hxx
\brief Defines an EOS

EOS is effectively an interface that describes how to create an equation of
state, where the independent quantities are rest mass density, specific energy,
and electron fraction.

*/

#ifndef EOS_3P_HXX
#define EOS_3P_HXX
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "utils/eos_utils.hxx"
#include "utils/eos_constants.hxx"

namespace EOSX {

/// Abstract class eos_3p

class eos_3p {
public:
  typedef eos_status status;
  typedef eos_range range;

  range rgrho;  ///< Valid range for density \f$ \rho \f$
  range rgye;   ///< Valid range for electron fraction \f$ Y_e \f$
  range rgtemp; ///< Valid range for temperature \f$ T \f$

protected:
  /// Set the density range. Has to be called in the constructor of any
  /// implementation.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_range_rho(const range &r) {
    rgrho = r;
  }
  /// Set the electron fraction range. Has to be called in the constructor of
  /// any implementation.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_range_ye(const range &r) {
    rgye = r;
  }
  /// Set the temperature range. Has to be called in the constructor of any
  /// implementation.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_range_temp(const range &r) {
    rgtemp = r;
  }

public:
  CCTK_DEVICE CCTK_HOST eos_3p() {}

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps,
                        const CCTK_REAL ye, status &stat) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps,
                       const CCTK_REAL ye, status &stat) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps,
                       const CCTK_REAL ye, status &stat) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_rho_eps_ye(CCTK_REAL &press, CCTK_REAL &dpdrho,
                               CCTK_REAL &dpdeps, const CCTK_REAL rho,
                               const CCTK_REAL eps, const CCTK_REAL ye,
                               status &stat) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                           const CCTK_REAL ye, status &stat) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps,
                          const CCTK_REAL ye, status &stat) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_csnd_from_rho_eps_ye(CCTK_REAL &press, CCTK_REAL &csnd,
                             const CCTK_REAL rho, const CCTK_REAL eps,
                             const CCTK_REAL ye, status &stat) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                       const CCTK_REAL ye, status &stat) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline range
  range_eps(const CCTK_REAL rho, const CCTK_REAL ye) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline const range &
  range_rho() const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline const range &
  range_ye() const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline const range &
  range_temp() const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_rho_valid(CCTK_REAL rho) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_ye_valid(CCTK_REAL ye) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_temp_valid(CCTK_REAL temp) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_eps_valid(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_rho_eps_ye_valid(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL ye) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_rho_temp_ye_valid(CCTK_REAL rho, CCTK_REAL temp, CCTK_REAL ye) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  check_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps, const CCTK_REAL ye,
                   status &stat) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  check_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                    const CCTK_REAL ye, status &stat) const;
  CCTK_DEVICE CCTK_HOST static CCTK_REAL nan() { return 0.0 / 0.0; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps,
                        const CCTK_REAL ye) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_rho_press_ye(const CCTK_REAL rho, const CCTK_REAL press,
                        const CCTK_REAL ye) const;
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                       const CCTK_REAL ye) const;
};

} // namespace EOSX

#endif
