#ifndef REDUCTION_HXX
#define REDUCTION_HXX

#include <cctk.h>
#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#include <mpi.h>

#include <cmath>

namespace CarpetX {
using namespace std;

template <typename T> T fmax1(T x, T y) {
  if (CCTK_isnan(x))
    return x;
  if (CCTK_isnan(y))
    return y;
  return fmax(x, y);
}
template <typename T> T fmin1(T x, T y) {
  if (CCTK_isnan(x))
    return x;
  if (CCTK_isnan(y))
    return y;
  return fmin(x, y);
}

template <typename T> struct reduction {
  // TODO: contains_inf, contains_nan?
  // TODO: minloc, maxloc
  T min, max, sum, sum2;
  T vol, maxabs, sumabs, sum2abs;
  reduction();
  reduction(const T &V, const T &x);
  reduction(const reduction &x, const reduction &y);
  reduction operator+(const reduction &y) const { return reduction(*this, y); }
  reduction &operator+=(const reduction &y) { return *this = *this + y; }

  T avg() const { return sum / vol; }
  T sdv() const { return sqrt(fmax1(T(0.0), vol * sum2 - pow(sum, 2))); }
  T norm0() const { return vol; }
  T norm1() const { return sumabs / vol; }
  T norm2() const { return sqrt(sum2abs / vol); }
  T norm_inf() const { return maxabs; }
};

template <typename T>
reduction<T>::reduction()
    : min(1.0 / 0.0), max(-1.0 / 0.0), sum(0.0), sum2(0.0), vol(0.0),
      maxabs(0.0), sumabs(0.0), sum2abs(0.0) {}

template <typename T>
reduction<T>::reduction(const T &V, const T &x)
    : min(x), max(x), sum(V * x), sum2(pow(V * x, 2)), vol(V), maxabs(fabs(x)),
      sumabs(fabs(V * x)), sum2abs(pow(fabs(V * x), 2)) {}

template <typename T>
reduction<T>::reduction(const reduction &x, const reduction &y)
    : min(fmin1(x.min, y.min)), max(fmax1(x.max, y.max)), sum(x.sum + y.sum),
      sum2(x.sum2 + y.sum2), vol(x.vol + y.vol),
      maxabs(fmax1(x.maxabs, y.maxabs)), sumabs(x.sumabs + y.sumabs),
      sum2abs(x.sum2abs + y.sum2abs) {}

typedef reduction<CCTK_REAL> reduction_CCTK_REAL;
#pragma omp declare reduction(reduction:reduction_CCTK_REAL : omp_out += omp_in)

MPI_Op reduction_mpi_op();

reduction<CCTK_REAL> reduce(int gi, int vi, int tl);

} // namespace CarpetX

#endif // #ifndef REDUCTION_HXX