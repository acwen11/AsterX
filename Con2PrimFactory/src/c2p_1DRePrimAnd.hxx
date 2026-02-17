#ifndef C2P_1DREPRIMAND_HXX
#define C2P_1DREPRIMAND_HXX

#include <algorithm>
#include <cmath>

#include "c2p.hxx"
#include "c2p_report.hxx"
#include "prims.hxx"
#include "cons.hxx"
#include "roots.hxx"

#include "c2p_1DRePrimAnd_intervals.hxx"
#include "c2p_1DRePrimAnd_rootfinder.hxx"
#include "c2p_1DRePrimAnd_internals.hxx"

namespace Con2PrimFactory {

class c2p_1DRePrimAnd : public c2p {
public:
  CCTK_REAL GammaIdealFluid;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DRePrimAnd(
      const EOSType *eos_3p, const atmosphere &atm_in, CCTK_INT maxIter,
      CCTK_REAL tol, CCTK_REAL alp_thresh_in, CCTK_REAL consError,
      CCTK_REAL vwlim, CCTK_REAL B_lim, CCTK_REAL rho_BH_in,
      CCTK_REAL eps_BH_in, CCTK_REAL vwlim_BH_in, bool ye_len, bool use_z,
      bool use_temperature, bool use_pressure_atmo) {

    atmo = atm_in;
    maxIterations = maxIter;
    tolerance = tol;
    alp_thresh = alp_thresh_in;
    cons_error = consError;
    vw_lim = vwlim;
    w_lim = sqrt(1.0 + vw_lim * vw_lim);
    v_lim = vw_lim / w_lim;
    Bsq_lim = B_lim * B_lim;
    rho_BH = rho_BH_in;
    eps_BH = eps_BH_in;
    vwlim_BH = vwlim_BH_in;
    ye_lenient = ye_len;
    use_zprim = use_z;
    use_temp = use_temperature;
    use_press_atmo = use_pressure_atmo;

    GammaIdealFluid = eos_3p->gamma;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_to_nan(prim_vars &pv, cons_vars &cv) const {
    pv.set_to_nan();
    cv.set_to_nan();
  }

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(const EOSType *eos_3p, prim_vars &pv, cons_vars &cv,
        const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const {
    rep.iters = 0;
    rep.adjust_cons = false;
    rep.set_atmo = false;
    rep.status = c2p_report::SUCCESS;

    const CCTK_REAL spatial_detg = calc_det(glo);
    const CCTK_REAL sqrt_detg = sqrt(spatial_detg);
    const bool minor1{glo(X, X) > 0.0};
    const bool minor2{glo(X, X) * glo(Y, Y) - glo(X, Y) * glo(X, Y) > 0.0};
    const bool minor3{spatial_detg > 0.0};
    if (!(minor1 && minor2 && minor3)) {
      rep.set_invalid_detg(sqrt_detg);
      set_to_nan(pv, cv);
      return;
    }

    const smat<CCTK_REAL, 3> gup = calc_inv(glo, spatial_detg);

    const cons_vars cv_const = cv;

    cv.dens /= sqrt_detg;
    cv.tau /= sqrt_detg;
    cv.mom /= sqrt_detg;
    cv.dBvec /= sqrt_detg;
    cv.DYe /= sqrt_detg;
    cv.DEnt /= sqrt_detg;

    const CCTK_REAL Ssq =
        calc_contraction(calc_contraction(gup, cv.mom), cv.mom);
    const CCTK_REAL Bsq =
        calc_contraction(calc_contraction(glo, cv.dBvec), cv.dBvec);
    const CCTK_REAL BiSi = cv.dBvec(X) * cv.mom(X) + cv.dBvec(Y) * cv.mom(Y) +
                           cv.dBvec(Z) * cv.mom(Z);

    if ((!isfinite(cv.dens)) || (!isfinite(Ssq)) || (!isfinite(Bsq)) ||
        (!isfinite(BiSi)) || (!isfinite(cv.DYe)) || (!isfinite(cv.DEnt))) {
      rep.set_nans_in_cons(cv.dens, Ssq, Bsq, BiSi, cv.DYe);
      set_to_nan(pv, cv);
      return;
    }

    if (Bsq > Bsq_lim) {
      rep.set_B_limit(Bsq);
      set_to_nan(pv, cv);
      return;
    }

    // Ye needs to be a non-const lvalue for EOS calls
    CCTK_REAL Ye_raw = cv.DYe / cv.dens;
    CCTK_REAL Ye = fmin(fmax(eos_3p->rgye.min, Ye_raw), eos_3p->rgye.max);

    typename RePrimAnd::froot<EOSType>::cache cache{};
    RePrimAnd::froot<EOSType> f(
        eos_3p, Ye, cv.dens, cv.tau / cv.dens, Ssq / (cv.dens * cv.dens),
        (BiSi * BiSi) / (cv.dens * cv.dens * cv.dens), Bsq / cv.dens, cache);

    interval<CCTK_REAL> mu_br = f.initial_bracket();

    interval<CCTK_REAL> rgrho_iv{eos_3p->rgrho.min, eos_3p->rgrho.max};
    RePrimAnd::rarecase<EOSType> rc(mu_br, rgrho_iv, f);
    mu_br = rc.bracket;

    const CCTK_REAL log2 = std::log(2.0);
    const CCTK_INT minbits = int(std::abs(std::log(tolerance)) / log2);
    const CCTK_REAL tolerance_0 = std::ldexp(double(1.0), -minbits);
    const CCTK_INT maxiters = maxIterations;

    auto fn = [&](CCTK_REAL mu) { return f(mu); };
    if (fn(mu_br.min()) * fn(mu_br.max()) > 0) {
      const CCTK_REAL qtot = cv.tau / cv.dens;
      const CCTK_REAL sPal = Bsq / cv.dens;
      CCTK_REAL new_hi = 3.0 + 3.0 * qtot - 1.5 * sPal;
      mu_br = interval<CCTK_REAL>{mu_br.min(), new_hi};
    }

    auto result =
        Algo::brent(fn, mu_br.min(), mu_br.max(), minbits, maxiters, rep.iters);

    const CCTK_REAL a_root = result.first;
    const CCTK_REAL b_root = result.second;
    const CCTK_REAL fa = fn(a_root);
    const CCTK_REAL fb = fn(b_root);
    const CCTK_REAL mu =
        (fb == (CCTK_REAL)0 || std::abs(fb) < std::abs(fa)) ? b_root
        : (std::abs(fa) < std::abs(fb))                     ? a_root
                                        : CCTK_REAL(0.5) * (a_root + b_root);

    // ------------------------------------------------------------------
    // IMPORTANT: ensure cache corresponds to the final chosen mu
    // ------------------------------------------------------------------
    (void)f(mu);

    // ------------------------------------------------------------------
    // Use the EOS-consistent RePrimAnd cached primitives
    // ------------------------------------------------------------------
    pv.rho = cache.rho;
    pv.Ye = cache.ye;
    pv.eps = cache.eps;
    pv.press = cache.press;
    pv.w_lor = cache.w;

    pv.temperature = eos_3p->temp_from_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.entropy = eos_3p->kappa_from_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    // ------------------------------------------------------------------
    // Velocity: RePrimAnd magnetic-aware reconstruction
    //   v^i = mu * x * ( r^i + (r.b) * mu * b^i )
    // where r^i = S^i / D, b^i = B^i / sqrt(D), (r.b) = (B^i S_i)/(D^(3/2)).
    // ------------------------------------------------------------------

    const vec<CCTK_REAL, 3> mom_up = calc_contraction(gup, cv.mom);

    const CCTK_REAL D = cv.dens;
    const CCTK_REAL sqD = std::sqrt(D);

    const vec<CCTK_REAL, 3> r_u = mom_up * (1.0 / (D + 1e-300));
    const vec<CCTK_REAL, 3> b_u = cv.dBvec * (1.0 / (sqD + 1e-300));
    const CCTK_REAL rb = BiSi / (D * (sqD + 1e-300));

    pv.vel = (mu * cache.x) * (r_u + (rb * mu) * b_u);

    vec<CCTK_REAL, 3> v_low = calc_contraction(glo, pv.vel);
    CCTK_REAL vsq = calc_contraction(v_low, pv.vel);
    if (vsq >= v_lim * v_lim) {
      const CCTK_REAL scale = v_lim / (std::sqrt(vsq) + 1e-300);
      pv.vel *= scale;
      vsq = v_lim * v_lim;
    }
    pv.w_lor = CCTK_REAL(1) / std::sqrt(std::max(1e-32, 1.0 - vsq));

    pv.Bvec = cv.dBvec;
    const vec<CCTK_REAL, 3> Elow = calc_cross_product(pv.Bvec, pv.vel);
    pv.E = calc_contraction(gup, Elow);

    if (std::abs(result.first - result.second) >
        tolerance_0 *
            std::min(std::abs(result.first), std::abs(result.second))) {

      cons_vars cv_check;
      cv_check.from_prim(pv, glo);

      cv_check.dens /= sqrt_detg;
      cv_check.tau /= sqrt_detg;
      cv_check.mom /= sqrt_detg;
      cv_check.dBvec /= sqrt_detg;
      cv_check.DYe /= sqrt_detg;
      cv_check.DEnt /= sqrt_detg;

      const CCTK_REAL small = 1e-50;
      const CCTK_REAL max_error = sqrt(
          max({pow((cv_check.dens - cv.dens) / (cv.dens + small), 2.0),
               pow((cv_check.mom(0) - cv.mom(0)) / (cv.mom(0) + small), 2.0),
               pow((cv_check.mom(1) - cv.mom(1)) / (cv.mom(1) + small), 2.0),
               pow((cv_check.mom(2) - cv.mom(2)) / (cv.mom(2) + small), 2.0),
               pow((cv_check.tau - cv.tau) / (cv.tau + small), 2.0)}));

      if (max_error > cons_error) {
        rep.set_root_conv();
        cv = cv_const;
        return;
      }
    }

    if (pv.rho < atmo.rho_cut) {
      rep.set_atmo_set();
      atmo.set(pv, cv, glo);
      return;
    }

    c2p::prims_floors_and_ceilings(eos_3p, pv, cv, glo, rep);

    if (rep.adjust_cons) {
      cv.from_prim(pv, glo);
      cv.dBvec = cv_const.dBvec;
    } else {
      cv = cv_const;
      cv.DEnt = cv.dens * pv.entropy;
    }
  }

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  operator()(const EOSType *eos_3p, prim_vars &pv, cons_vars &cv,
             const smat<CCTK_REAL, 3> & /*gup_unused*/,
             const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const {
    solve(eos_3p, pv, cv, glo, rep);
  }
};

} // namespace Con2PrimFactory

#endif
