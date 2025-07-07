/*! \file cons.hxx
\brief Class definitions representing conserved mhd variables.
*/

#ifndef CONS_HXX
#define CONS_HXX

#include "prims.hxx"

namespace Con2PrimFactory {

using namespace AsterUtils;

// Structure to represent conserved variables for mhd
struct cons_vars {
  CCTK_REAL dens;
  vec<CCTK_REAL, 3> mom;
  CCTK_REAL tau;
  CCTK_REAL DYe;
  CCTK_REAL DEnt;
  vec<CCTK_REAL, 3> dBvec;

  // Default constructor, no initialization.
  cons_vars() = default;

  // Construct from single variables.
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline cons_vars(
      CCTK_REAL dens_, vec<CCTK_REAL, 3> mom_, CCTK_REAL tau_, CCTK_REAL DYe_,
      CCTK_REAL DEnt_, vec<CCTK_REAL, 3> dBvec_)
      : dens{dens_}, mom{mom_}, tau{tau_}, DYe{DYe_}, DEnt{DEnt_},
        dBvec{dBvec_} {}

  /// Copy assignment
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline cons_vars &
  operator=(const cons_vars &other) {
    if (this == &other)
      return *this; // Handle self-assignment
    // Copy data members from 'other' to 'this'
    dens = other.dens;
    mom = other.mom;
    tau = other.tau;
    DYe = other.DYe;
    DEnt = other.DEnt;
    dBvec = other.dBvec;
    return *this;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  from_prim(const prim_vars &pv, const smat<CCTK_REAL, 3> &g) {

    // determinant of spatial metric
    const CCTK_REAL sqrt_detg = sqrt(calc_det(g));

    /* Computing v_j */
    const vec<CCTK_REAL, 3> &v_up = pv.vel;
    const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);

    const CCTK_REAL w_lorentz = pv.w_lor;
    ;
    // const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

    /* Computing B_j */
    const vec<CCTK_REAL, 3> &B_up = pv.Bvec;
    const vec<CCTK_REAL, 3> B_low = calc_contraction(g, B_up);

    /* Computing b^t : this is b^0 * alp */
    const CCTK_REAL bst = w_lorentz * calc_contraction(B_up, v_low);

    /* Computing b_j */
    const vec<CCTK_REAL, 3> b_low = B_low / w_lorentz + bst * v_low;

    /* Computing b^mu b_mu */
    const CCTK_REAL bs2 =
        (calc_contraction(B_up, B_low) + bst * bst) / (w_lorentz * w_lorentz);

    // computing conserved from primitives
    dens = sqrt_detg * pv.rho * w_lorentz;

    mom = sqrt_detg * (w_lorentz * w_lorentz *
                           (pv.rho * (1.0 + pv.eps) + pv.press + bs2) * v_low -
                       bst * b_low);

    tau = sqrt_detg * (w_lorentz * w_lorentz *
                           (pv.rho * (1.0 + pv.eps) + pv.press + bs2) -
                       (pv.press + 0.5 * bs2) - bst * bst) -
          dens;

    dBvec = sqrt_detg * pv.Bvec;
    DYe = dens * pv.Ye;
    DEnt = dens * pv.entropy;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  scatter(CCTK_REAL &dens_, CCTK_REAL &momx_, CCTK_REAL &momy_,
          CCTK_REAL &momz_, CCTK_REAL &tau_, CCTK_REAL &DYe_, CCTK_REAL &DEnt_,
          CCTK_REAL &dBvecx_, CCTK_REAL &dBvecy_, CCTK_REAL &dBvecz_) const {

    dens_ = dens;
    momx_ = mom(0);
    momy_ = mom(1);
    momz_ = mom(2);
    tau_ = tau;
    DYe_ = DYe;
    DEnt_ = DEnt;
    dBvecx_ = dBvec(0);
    dBvecy_ = dBvec(1);
    dBvecz_ = dBvec(2);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void set_to_nan() {
    dens = mom(0) = mom(1) = mom(2) = tau = DYe = DEnt = dBvec(0) = dBvec(1) =
        dBvec(2) = std::numeric_limits<CCTK_REAL>::quiet_NaN();
  }
};

} // namespace Con2PrimFactory
#endif
