#ifndef EOS_1P_HXX
#define EOS_1P_HXX
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "utils/eos_utils.hxx"
#include "utils/eos_constants.hxx"

namespace EOSX {

class eos_1p {
public:
  typedef eos_range range;

private:
  bool isentropic;  ///< Whether EOS is isentropic
  bool temp_avail;  ///< Whether EOS has temperature information
  bool efrac_avail; ///< Whether EOS has electron fraction information
  range rg_rho;     ///< Density range in which EOS is valid.
  range rg_gm1;     ///< Range of \f$ g-1 \f$ in which EOS is valid.
  range rg_p;       ///< Pressure range in which EOS is valid.

public:
  /// Set valid density range. Has to be called by constructor of each
  /// implementation.
  CCTK_DEVICE CCTK_HOST void set_ranges(const range &rg_rho_);
  /// Check if density is in valid range, else throw exception.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  check_rho(CCTK_REAL rho) const;
  /// Check if \f$ g-1 \f$ is in valid range, else throw exception.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  check_gm1(CCTK_REAL gm1) const;
  /// Check if pressure is in valid range, else throw exception.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  check_p(CCTK_REAL p) const;

  /// Whether EOS is isentropic
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_isentropic() const {
    return isentropic;
  }
  /// Whether EOS can compute temperature
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  has_temp() const {
    return temp_avail;
  }
  /// Whether EOS can compute electron fraction
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  has_efrac() const {
    return efrac_avail;
  }
  /// Returns range of validity for density
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline const range &
  range_rho() const {
    return rg_rho;
  }
  /// Returns range of validity for \f$ g-1 \f$
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline const range &
  range_gm1() const {
    return rg_gm1;
  }
  /// Returns range of validity for pressure
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline const range &
  range_p() const {
    return rg_p;
  }
  /// Whether given density is in the valid range of the EOS
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_rho_valid(CCTK_REAL rho) const {
    return rg_rho.contains(rho);
  }
  /// Whether given \f$ g-1 \f$ is in the valid range of the EOS
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_gm1_valid(CCTK_REAL gm1) const {
    return rg_gm1.contains(gm1);
  }
  /// Whether given pressure is in the valid range of the EOS
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  is_p_valid(CCTK_REAL p) const {
    return rg_p.contains(p);
  }
};

} // namespace EOSX

#endif
