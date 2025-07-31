/* file eos_hybrid.hxx
 * child class of eos, implements hybrid EoS given barotropic EoS from
 * eos_1p.hxx added by Johnny Tsao btsao@utexas.edu hybrid EoS following the
 * convention in thcextra/EOS_Thermal_Hybrid/src/eos_hybrid.h by Wolfgang
 * Kastaun physik@fangwolg.de
 */

#ifndef EOS_3P_HYBRID_HXX
#define EOS_3P_HYBRID_HXX

#include <stdexcept>
#include <algorithm>
#include <cmath>

#include "../eos_3p.hxx"
#include "../eos_1p_polytropic/eos_1p_polytropic.hxx"

using namespace std;

namespace EOSX {

class eos_3p_hybrid : public eos_3p {
public:
  // CCTK_REAL gamma_th, gm1, temp_over_eps;
  CCTK_REAL gamma;            // dummy
  CCTK_REAL gamma_th, gm1_th; // TODO: temp_over_eps
  range rgeps;
  eos_1p_polytropic *eos_c;

  // constructor
  CCTK_HOST
  CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_3p_hybrid(
      eos_1p_polytropic *eos_c_, CCTK_REAL gamma_th_, range &rgrho_,
      range &rgeps_, range &rgye_)
      : eos_c(eos_c_), gamma_th(gamma_th_), rgeps(rgeps_) {

    if (gamma_th < 1) {
      assert(0);
      // runtime_error("EOS_Hybrid: initialized with gamma_th < 1");
    }
    if (gamma_th > 2) { // Ensure subluminal Soundspeed and P < E
      rgeps.max = min(rgeps.max, 1 / (gamma_th * (gamma_th - 2)));
    }

    // temp_over_eps = gm1 * umass_; // TODO: temp_over_eps = (gamma - 1) *
    // u_mass, write in no gamma form

    gm1_th = gamma_th - 1.0;
    set_range_rho(rgrho_);
    set_range_ye(rgye_);
    // set_range_temp(range(temp_over_eps * rgeps.min, temp_over_eps *
    // rgeps.max));
    set_range_temp(range(0, 0)); // TODO: Implement temperature;
  }
  // ranges found with {rho, eps, ye} instead of {rho, temp ye}
  // range for rho can be replaced with {0, rho_max}

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_cold(const CCTK_REAL rho) const {
    CCTK_REAL gm1 =
        eos_c->gm1_from_valid_rho(rho); // TODO: name change -> gm1 here is hm1
    return eos_c->sed_from_valid_gm1(gm1);
  }

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  p_cold(const CCTK_REAL rho) const {
    CCTK_REAL gm1 =
        eos_c->gm1_from_valid_rho(rho); // TODO: name change -> gm1 here is hm1
    return eos_c->p_from_valid_gm1(gm1);
  }

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  hm1_cold(const CCTK_REAL rho) const {
    CCTK_REAL gm1 =
        eos_c->gm1_from_valid_rho(rho); // TODO: name change -> gm1 here is hm1
    return eos_c->hm1_from_valid_gm1(gm1);
  }

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  cs2_cold(const CCTK_REAL rho) const {
    CCTK_REAL gm1 =
        eos_c->gm1_from_valid_rho(rho); // TODO: name change -> gm1 here is hm1
    return eos_c->csnd2_from_valid_gm1(gm1);
  }

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                              const CCTK_REAL ye) const {
    CCTK_REAL p_c = p_cold(rho);
    CCTK_REAL eps_c = eps_cold(rho);
    CCTK_REAL p_th = gm1_th * rho * (eps - eps_c);
    return p_c + p_th;
  }

  // this is not in thc, TODO: check if correct
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_press_ye(const CCTK_REAL rho, const CCTK_REAL press,
                              const CCTK_REAL ye) const {
    CCTK_REAL p_c = p_cold(rho);
    CCTK_REAL eps_c = eps_cold(rho);
    CCTK_REAL p_th = press - p_c;
    CCTK_REAL eps_th = p_th / gm1_th / rho;
    return eps_c + eps_th;
  }

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                             const CCTK_REAL ye) const {
    CCTK_REAL cs2_c = cs2_cold(rho);
    CCTK_REAL eps_c = eps_cold(rho);
    CCTK_REAL h_c = 1.0 + hm1_cold(rho);
    CCTK_REAL eps_th = eps - eps_c;
    CCTK_REAL h_th = gamma_th * eps_th;
    CCTK_REAL w = h_th / (h_c + h_th);
    CCTK_REAL cs2 = (1.0 - w) * cs2_c + w * gm1_th;
    return sqrt(cs2);
  }

  // temperature is not yet implemented in thc
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                             const CCTK_REAL ye) const {
    // return temp_over_eps * eps;
    printf("eos_3p_hybrid: temperature not implemented");
    return 0.0;
  }

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_valid_rho_eps_ye(CCTK_REAL &press, CCTK_REAL &dpdrho,
                                     CCTK_REAL &dpdeps, const CCTK_REAL rho,
                                     const CCTK_REAL eps,
                                     const CCTK_REAL ye) const {
    CCTK_REAL p_c = p_cold(rho);
    CCTK_REAL cs2_c = cs2_cold(rho);
    CCTK_REAL eps_c = eps_cold(rho);
    CCTK_REAL h_c = 1.0 + hm1_cold(rho);
    CCTK_REAL eps_th = eps - eps_c;
    CCTK_REAL p_th = gm1_th * rho * eps_th;
    press = p_c + p_th;
    dpdrho = h_c * cs2_c + gm1_th * (eps_th - p_c / rho);
    dpdeps = gm1_th * rho;
  }

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  entropy_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                                 const CCTK_REAL ye) const {
    // return log(temp * pow(rho, -gm1_th) / temp_over_eps);
    printf("EOS: entropy from temperature not implemented for eos_3p_hybrid.");
  }

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_valid_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps,
                                const CCTK_REAL ye) const {
    CCTK_REAL eps_th = eps - eps_cold(rho);
    return log(eps_th * pow(rho, -gm1_th));
  }

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                             const CCTK_REAL ye) const {
    // return temp / temp_over_eps;
    printf("EOS: eps from temperature not implemented for eos_3p_hybrid.");
    return 0.0;
  }

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                               const CCTK_REAL ye) const {
    // return temp / temp_over_eps;
    printf("EOS: press from temperature not implemented for eos_3p_hybrid.");
    return 0.0;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_kappa_ye(const CCTK_REAL rho,
                                const CCTK_REAL kappa, // p/rho^gamma_th
                                const CCTK_REAL ye) const {
    return kappa * pow(rho, gamma_th);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_kappa_ye(const CCTK_REAL rho,
                              const CCTK_REAL kappa, // p/rho^gamma
                              const CCTK_REAL ye) const {
    return kappa * pow(rho, gm1_th) / gm1_th;
  };

  // Note that kappa implements a generic thermodynamic quantity
  // meant to describe the "evolved" entropy by an evolution/application
  // thorn.
  // The notion of the "evolved" entropy (kappa) might differ from the
  // definition of the actual entropy (entropy_from..., see above) for different
  // EOS, e.g. for the thermal part of hybrid EOS we have kappa = p *
  // (rho)^(-gamma_th), where gamma is the adiabatic index.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  kappa_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                              const CCTK_REAL ye) const {
    CCTK_REAL eps_th = eps - eps_cold(rho);
    return gm1_th * eps_th * pow(rho, gm1_th);
  };

  // edited
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_3p::range
  range_eps_from_valid_rho_ye(const CCTK_REAL rho, const CCTK_REAL ye) const {
    return rgeps;
  }
};
} // namespace EOSX

#endif
