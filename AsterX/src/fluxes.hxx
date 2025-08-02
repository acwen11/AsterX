#ifndef ASTERX_FLUXES_HXX
#define ASTERX_FLUXES_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Arith;

// Lax-Friedrichs solver
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
laxf(vec<vec<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
     vec<CCTK_REAL, 2> flux) {
  const CCTK_REAL charmax =
      max({CCTK_REAL(0), fabs(lam(0)(0)), fabs(lam(0)(1)), fabs(lam(0)(2)),
           fabs(lam(0)(3)), fabs(lam(1)(0)), fabs(lam(1)(1)), fabs(lam(1)(2)),
           fabs(lam(1)(3))});

  return 0.5 * ((flux(0) + flux(1)) - charmax * (var(1) - var(0)));
}

// HLLE solver
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
hlle(vec<vec<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
     vec<CCTK_REAL, 2> flux) {
  const CCTK_REAL charmax =
      max({CCTK_REAL(0), lam(0)(0), lam(0)(1), lam(0)(2), lam(0)(3), lam(1)(0),
           lam(1)(1), lam(1)(2), lam(1)(3)});

  // Note that charmin is just the minimum, not with the minus sign
  const CCTK_REAL charmin =
      min({CCTK_REAL(0), lam(0)(0), lam(0)(1), lam(0)(2), lam(0)(3), lam(1)(0),
           lam(1)(1), lam(1)(2), lam(1)(3)});

  const CCTK_REAL charpm = charmax - charmin;

  return (charmax * flux(0) - charmin * flux(1) +
          charmax * charmin * (var(1) - var(0))) /
         charpm;
}

struct recon_prims {
  vec<CCTK_REAL, 2> rho_rc;
  vec<CCTK_REAL, 2> entropy_rc;
  vec<CCTK_REAL, 2> Ye_rc;
  vec<CCTK_REAL, 2> eps_rc;
  vec<CCTK_REAL, 2> press_rc;
  vec<CCTK_REAL, 2> temp_rc;
  vec<vec<CCTK_REAL, 2>, 3> Bs_rc;
  vec<vec<CCTK_REAL, 2>, 3> vels_rc;
  vec<vec<CCTK_REAL, 2>, 3> vlows_rc;
  vec<CCTK_REAL, 2> w_lorentz_rc;

  // Default constructor. Leaves all members uninitialized.
  recon_prims() = default;
  
  /// Construct from single variables.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline recon_prims(
    vec<CCTK_REAL, 2> rho_rc_, vec<CCTK_REAL, 2> entropy_rc_, vec<CCTK_REAL, 2> Ye_rc_, vec<CCTK_REAL, 2> eps_rc_,
    vec<CCTK_REAL, 2> press_rc_, vec<CCTK_REAL, 2> temp_rc_, vec<vec<CCTK_REAL, 2>, 3> Bs_rc_, vec<vec<CCTK_REAL, 2>, 3> vels_rc_,
    vec<vec<CCTK_REAL, 2>, 3> vlows_rc_, vec<CCTK_REAL, 2> w_lorentz_rc_)
    : rho_rc(rho_rc_), entropy_rc(entropy_rc_), Ye_rc(Ye_rc_), eps_rc(eps_rc_), press_rc(press_rc_), temp_rc(temp_rc_),
    Bs_rc(Bs_rc_), vels_rc(vels_rc_), vlows_rc(vlows_rc_), w_lorentz_rc(w_lorentz_rc_){};
};
  
struct recon_cons {
  vec<CCTK_REAL, 2> dens_rc;
  vec<CCTK_REAL, 2> DEnt_rc;
  vec<vec<CCTK_REAL, 2>, 3> moms_rc;
  vec<CCTK_REAL, 2> tau_rc;
  vec<CCTK_REAL, 2> DYe_rc;
  vec<vec<CCTK_REAL, 2>, 3> Btildes_rc;
  vec<vec<CCTK_REAL, 2>, 3> vtildes_rc;
  vec<vec<CCTK_REAL, 4>, 2> lambda;
  vec<CCTK_REAL, 2> flux_dens;
  vec<CCTK_REAL, 2> flux_DEnt;
  vec<vec<CCTK_REAL, 2>, 3> flux_moms;
  vec<CCTK_REAL, 2> flux_tau;
  vec<CCTK_REAL, 2> flux_DYe;
  vec<vec<CCTK_REAL, 2>, 3> flux_Btildes;

  // Default constructor. Leaves all members uninitialized.
  recon_cons() = default;
  
  /// Construct from single variables.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline recon_cons(
    vec<CCTK_REAL, 2> dens_rc_, vec<CCTK_REAL, 2> DEnt_rc_, vec<vec<CCTK_REAL, 2>, 3> moms_rc_, vec<CCTK_REAL, 2> tau_rc_, vec<CCTK_REAL, 2> DYe_rc_,
    vec<vec<CCTK_REAL, 2>, 3> Btildes_rc_, vec<vec<CCTK_REAL, 2>, 3> vtildes_rc_, vec<vec<CCTK_REAL, 4>, 2> lambda_,
    vec<CCTK_REAL, 2> flux_dens_, vec<CCTK_REAL, 2> flux_DEnt_, vec<vec<CCTK_REAL, 2>, 3> flux_moms_, vec<CCTK_REAL, 2> flux_tau_, vec<CCTK_REAL, 2> flux_DYe_,
    vec<vec<CCTK_REAL, 2>, 3> flux_Btildes_)
    : dens_rc(dens_rc_), DEnt_rc(DEnt_rc_), moms_rc(moms_rc_), tau_rc(tau_rc_), DYe_rc(DYe_rc_), Btildes_rc(Btildes_rc_), vtildes_rc(vtildes_rc_),
    lambda(lambda_),
    flux_dens(flux_dens_), flux_DEnt(flux_DEnt_), flux_moms(flux_moms_),
    flux_tau(flux_tau_), flux_DYe(flux_DYe_), flux_Btildes(flux_Btildes_){};

  /// Copy assignment
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline recon_cons &
    operator=(const recon_cons &other) {
      if (this == &other)
        return *this; // Handle self-assignment
    // Copy data members from 'other' to 'this'
    dens_rc = other.dens_rc;
    DEnt_rc = other.DEnt_rc;
    moms_rc = other.moms_rc;
    tau_rc = other.tau_rc;
    DYe_rc = other.DYe_rc;
    Btildes_rc = other.Btildes_rc;
    vtildes_rc = other.vtildes_rc;
    lambda = other.lambda;
    flux_dens = other.flux_dens;
    flux_DEnt = other.flux_DEnt;
    flux_moms = other.flux_moms;
    flux_tau = other.flux_tau;
    flux_DYe = other.flux_DYe;
    flux_Btildes = other.flux_Btildes;
    return *this;
  }
  
};


} // namespace AsterX

#endif // ASTERX_FLUXES_HXX
