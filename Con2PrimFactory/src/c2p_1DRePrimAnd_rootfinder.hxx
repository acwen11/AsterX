#ifndef C2P_1DREPRIMAND_ROOTFINDER_HXX
#define C2P_1DREPRIMAND_ROOTFINDER_HXX

#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

#include "roots.hxx" // Algo::brent
#include "c2p_1DRePrimAnd_intervals.hxx"
#include "c2p_utils.hxx"

namespace Con2PrimFactory {

// Derivative-assisted adapter (we still just use Brent on the residual)
template <class F, class T = typename F::value_t>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
findroot_using_deriv(F &f, interval<T> bracket, ROOTSTAT &errs, int digits,
                     unsigned int max_calls = 256) {
  const int minbits = digits;
  int iters = 0;

  auto fn = [&](T x) -> T {
    auto val = f(x);
    if constexpr (std::is_same<decltype(val), std::pair<T, T> >::value)
      return val.first;
    else
      return val;
  };

  const T fa = fn(bracket.min()), fb = fn(bracket.max());
  if (fa * fb > T(0)) {
    errs = ROOTSTAT::NOT_BRACKETED;
    return std::numeric_limits<T>::quiet_NaN();
  }

  auto ab =
      Algo::brent(fn, bracket.min(), bracket.max(), minbits, max_calls, iters);
  errs =
      (iters >= (int)max_calls) ? ROOTSTAT::NOT_CONVERGED : ROOTSTAT::SUCCESS;

  const T a_root = ab.first, b_root = ab.second;
  const T fa2 = fn(a_root), fb2 = fn(b_root);

  if (fb2 == T(0) || std::abs(fb2) < std::abs(fa2))
    return b_root;
  if (std::abs(fa2) < std::abs(fb2))
    return a_root;
  return (a_root + b_root) * T(0.5);
}

// Derivative-free Brent adapter
template <class F, class T = typename F::value_t>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval<T>
findroot_no_deriv(F &f, interval<T> bracket, T /*acc*/, unsigned int max_calls,
                  ROOTSTAT &errs) {
  const int minbits = std::numeric_limits<T>::digits - 4;
  int iters = 0;

  auto fn = [&](T x) -> T { return f(x); };

  const T fa = fn(bracket.min()), fb = fn(bracket.max());
  if (fa * fb > T(0)) {
    errs = ROOTSTAT::NOT_BRACKETED;
    return {bracket.min(), bracket.max()};
  }

  auto ab =
      Algo::brent(fn, bracket.min(), bracket.max(), minbits, max_calls, iters);
  errs =
      (iters >= (int)max_calls) ? ROOTSTAT::NOT_CONVERGED : ROOTSTAT::SUCCESS;

  return {ab.first, ab.second};
}

template <class F, class T = typename F::value_t>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
findroot_using_deriv(F &f, ROOTSTAT &errs, int digits,
                     unsigned int max_calls = 256) {
  return findroot_using_deriv(f, f.initial_bracket(), errs, digits, max_calls);
}

template <class F, class T = typename F::value_t>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline interval<T>
findroot_no_deriv(F &f, T tol, unsigned int max_calls, ROOTSTAT &errs) {
  return findroot_no_deriv(f, f.initial_bracket(), tol, max_calls, errs);
}

} // namespace Con2PrimFactory
#endif
