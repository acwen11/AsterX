#include "prolongate_3d_rf2.hxx"

#include <array>
#include <cassert>
#include <cmath>
#include <vector>

namespace AMReX {
using namespace amrex;
using namespace std;

namespace {
constexpr int VC = 0; // centering
constexpr int CC = 1;
constexpr bool NOCONS = false; // conservative
constexpr bool CONS = true;
} // namespace

// 1D interpolation coefficients

template <int CENTERING, bool CONSERVATIVE, int ORDER, typename T>
struct coeffs1d;
//  static constexpr const T coeffs[];

template <typename T> struct coeffs1d<VC, NOCONS, /*order*/ 1, T> {
  static constexpr array<T, 2> coeffs = {
      +1 / T(2),
      +1 / T(2),
  };
};
template <typename T> struct coeffs1d<VC, NOCONS, /*order*/ 3, T> {
  static constexpr array<T, 4> coeffs = {
      -1 / T(16),
      +9 / T(16),
      +9 / T(16),
      -1 / T(16),
  };
};
template <typename T> struct coeffs1d<VC, NOCONS, /*order*/ 5, T> {
  static constexpr array<T, 6> coeffs = {
      +3 / T(256),  -25 / T(256), +75 / T(128),
      +75 / T(128), -25 / T(256), +3 / T(256),
  };
};
template <typename T> struct coeffs1d<VC, NOCONS, /*order*/ 7, T> {
  static constexpr array<T, 8> coeffs = {
      -5 / T(2048),    +49 / T(2048),  -245 / T(2048), +1225 / T(2048),
      +1225 / T(2048), -245 / T(2048), +49 / T(2048),  -5 / T(2048),
  };
};

template <typename T> struct coeffs1d<CC, NOCONS, /*order*/ 0, T> {
  static constexpr array<T, 1> coeffs = {
      +1 / T(1),
  };
};
template <typename T> struct coeffs1d<CC, NOCONS, /*order*/ 1, T> {
  static constexpr array<T, 2> coeffs = {
      +1 / T(4),
      +3 / T(4),
  };
};
template <typename T> struct coeffs1d<CC, NOCONS, /*order*/ 2, T> {
  static constexpr array<T, 3> coeffs = {
      +5 / T(32),
      +15 / T(16),
      -3 / T(32),
  };
};
template <typename T> struct coeffs1d<CC, NOCONS, /*order*/ 3, T> {
  static constexpr array<T, 4> coeffs = {
      -5 / T(128),
      +35 / T(128),
      +105 / T(128),
      -7 / T(128),
  };
};
template <typename T> struct coeffs1d<CC, NOCONS, /*order*/ 4, T> {
  static constexpr array<T, 5> coeffs = {
      -45 / T(2048), +105 / T(512), +945 / T(1024), -63 / T(512), +35 / T(2048),
  };
};

template <typename T> struct coeffs1d<VC, CONS, /*order*/ 2, T> {
  static constexpr array<T, 3> coeffs0 = {
      -1 / T(32),
      +17 / T(16),
      -1 / T(32),
  };
  static constexpr array<T, 3> coeffs1 = {
      +13 / T(16),
      -5 / T(32),
      +11 / T(32),
  };
};
template <typename T> struct coeffs1d<VC, CONS, /*order*/ 4, T> {
  static constexpr array<T, 5> coeffs0 = {
      +7 / T(2048), -23 / T(512), +1109 / T(1024), -23 / T(512), +7 / T(2048),
  };
  static constexpr array<T, 5> coeffs1 = {
      +63 / T(2048), -103 / T(512), +781 / T(1024),
      +233 / T(512), -97 / T(2048),
  };
};

template <typename T> struct coeffs1d<CC, CONS, /*order*/ 0, T> {
  static constexpr array<T, 1> coeffs = {
      +1 / T(1),
  };
};
template <typename T> struct coeffs1d<CC, CONS, /*order*/ 2, T> {
  static constexpr array<T, 3> coeffs = {
      +1 / T(8),
      +1 / T(1),
      -1 / T(8),
  };
};
template <typename T> struct coeffs1d<CC, CONS, /*order*/ 4, T> {
  static constexpr array<T, 5> coeffs = {
      -3 / T(128), +11 / T(64), +1 / T(1), -11 / T(64), +3 / T(128),
  };
};

// 1D interpolation operators

template <int CENTERING, bool CONSERVATIVE, int ORDER> struct interp1d;
// template <typename T>
// inline T operator()(const T *restrict const crseptr, const ptrdiff_t di,
//                     const int off) const;

// off=0: on coarse point
// off=1: between coarse points
template <int ORDER> struct interp1d<VC, NOCONS, ORDER> {
  static_assert(ORDER % 2 == 1);
  template <typename T>
  inline T operator()(const T *restrict const crseptr, const ptrdiff_t di,
                      const int off) const {
#ifdef CCTK_DEBUG
    assert(off == 0 || off == 1);
#endif
    if (off == 0)
      return crseptr[0];
    constexpr array<T, ORDER + 1> cs = coeffs1d<VC, NOCONS, ORDER, T>::coeffs;
    constexpr int i0 = (ORDER + 1) / 2;
    T y = 0;
    // for (int i = 0; i < ORDER + 1; ++i)
    //   y += cs[i] * crseptr[(i - i0) * di];
    // Make use of symmetry in coefficients
    for (int i = 0; i < i0; ++i)
      y += cs[i] * (crseptr[(i - i0) * di] + crseptr[(i0 - 1 - i) * di]);
#ifdef CCTK_DEBUG
    assert(!isnan(y));
#endif
    return y;
  }
};

// off=0: left sub-cell
// off=1: right sub-cell
template <int ORDER> struct interp1d<CC, NOCONS, ORDER> {
  template <typename T>
  inline T operator()(const T *restrict const crseptr, const ptrdiff_t di,
                      const int off) const {
#ifdef CCTK_DEBUG
    assert(off == 0 || off == 1);
#endif
    constexpr array<T, ORDER + 1> cs = coeffs1d<CC, NOCONS, ORDER, T>::coeffs;
    constexpr int i0 = (ORDER + 1) / 2;
    T y = 0;
    if (off == 0)
      for (int i = 0; i < ORDER + 1; ++i)
        y += cs[i] * crseptr[(i - i0) * di];
    else
      for (int i = 0; i < ORDER + 1; ++i)
        y += cs[i] * crseptr[-(i - i0) * di];
#ifdef CCTK_DEBUG
    assert(!isnan(y));
#endif
    return y;
  }
};

// off=0: on coarse point
// off=1: between coarse points
template <int ORDER> struct interp1d<VC, CONS, ORDER> {
  static_assert(ORDER % 2 == 0);
  template <typename T>
  inline T operator()(const T *restrict const crseptr, const ptrdiff_t di,
                      const int off) const {
#ifdef CCTK_DEBUG
    assert(off == 0 || off == 1);
#endif
    constexpr int i0 = (ORDER + 1) / 2;
    T y = 0;
    if (off == 0) {
      constexpr array<T, ORDER + 1> cs = coeffs1d<VC, CONS, ORDER, T>::coeffs0;
      for (int i = 0; i < ORDER + 1; ++i)
        y += cs[i] * crseptr[(i - i0) * di];
    } else {
      constexpr array<T, ORDER + 1> cs = coeffs1d<VC, CONS, ORDER, T>::coeffs1;
      for (int i = 0; i < ORDER + 1; ++i)
        y += cs[i] * crseptr[(i - i0) * di];
    }
    return y;
  }
};

// off=0: left sub-cell
// off=1: right sub-cell
template <int ORDER> struct interp1d<CC, CONS, ORDER> {
  static_assert(ORDER % 2 == 0, "");
  template <typename T>
  inline T operator()(const T *restrict const crseptr, const ptrdiff_t di,
                      const int off) const {
#ifdef CCTK_DEBUG
    assert(off == 0 || off == 1);
#endif
    constexpr array<T, ORDER + 1> cs = coeffs1d<CC, CONS, ORDER, T>::coeffs;
    constexpr int i0 = (ORDER + 1) / 2;
    T y = 0;
    if (off == 0)
      for (int i = 0; i < ORDER + 1; ++i)
        y += cs[i] * crseptr[(i - i0) * di];
    else
      for (int i = 0; i < ORDER + 1; ++i)
        y += cs[i] * crseptr[-(i - i0) * di];
    return y;
  }
};

// Test 1d interpolators

template <int CENTERING, bool CONSERVATIVE, int ORDER, typename T>
struct test_interp1d;

template <int CENTERING, int ORDER, typename T>
struct test_interp1d<CENTERING, NOCONS, ORDER, T> {
  void operator()() const {
    for (int order = 0; order < ORDER; ++order) {
      auto f = [&](T x) { return order == 0 ? T(1) : pow(x, order); };
      constexpr int n = (ORDER + 1) / 2 * 2 + 1;
      array<T, n + 2> ys;
      ys[0] = ys[n + 1] = 0 / T(0);
      const int i0 = n / 2;
      for (int i = 0; i < n; ++i) {
        T x = (i - i0) + CENTERING / T(2);
        T y = f(x);
        ys[i + 1] = y;
      }
      for (int off = 0; off < 2; ++off) {
        T x = CENTERING / T(4) + off / T(2);
        T y = f(x);
        T y1 = interp1d<CENTERING, NOCONS, ORDER>()(&ys[i0 + 1], 1, off);
        assert(!isnan(y1));
        assert(y1 == y);
      }
    }
  }
};

////////////////////////////////////////////////////////////////////////////////

template <int CENTERING, bool CONSERVATIVE, int ORDER, int D, typename T>
void interp3d(const T *restrict const crseptr, const Box &restrict crsebox,
              T *restrict const fineptr, const Box &restrict finebox,
              const Box &restrict targetbox) {
  test_interp1d<CENTERING, CONSERVATIVE, ORDER, T>()();

  assert(crseptr);
  assert(crsebox.ok());
  assert(fineptr);
  assert(finebox.ok());
  assert(targetbox.ok());

  const IntVect first_crseind(finebox.loVect());
  IntVect next_crseind = first_crseind;
  ++next_crseind.getVect()[D];
  const ptrdiff_t di =
      crsebox.index(next_crseind) - crsebox.index(first_crseind);
  assert(di > 0);

  for (int k = targetbox.loVect()[2]; k <= targetbox.hiVect()[2]; ++k) {
    for (int j = targetbox.loVect()[1]; j <= targetbox.hiVect()[1]; ++j) {
#pragma omp simd
      for (int i = targetbox.loVect()[0]; i <= targetbox.hiVect()[0]; ++i) {
        const IntVect fineind(i, j, k);
        IntVect crseind = fineind;
        crseind.getVect()[D] = coarsen(crseind.getVect()[D], 2);
        const int off = fineind.getVect()[D] - crseind.getVect()[D] * 2;
        fineptr[finebox.index(fineind)] =
            interp1d<CENTERING, CONSERVATIVE, ORDER>()(
                &crseptr[crsebox.index(crseind)], di, off);
#ifdef CCTK_DEBUG
        assert(!isnan(fineptr[finebox.index(fineind)]));
#endif
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

template <int CENTI, int CENTJ, int CENTK, bool CONSI, bool CONSJ, bool CONSK,
          int ORDERI, int ORDERJ, int ORDERK>
prolongate_3d_rf2<CENTI, CENTJ, CENTK, CONSI, CONSJ, CONSK, ORDERI, ORDERJ,
                  ORDERK>::~prolongate_3d_rf2() {}

template <int CENTI, int CENTJ, int CENTK, bool CONSI, bool CONSJ, bool CONSK,
          int ORDERI, int ORDERJ, int ORDERK>
Box prolongate_3d_rf2<CENTI, CENTJ, CENTK, CONSI, CONSJ, CONSK, ORDERI, ORDERJ,
                      ORDERK>::CoarseBox(const Box &fine, int ratio) {
  return CoarseBox(fine, IntVect(ratio));
}

template <int CENTI, int CENTJ, int CENTK, bool CONSI, bool CONSJ, bool CONSK,
          int ORDERI, int ORDERJ, int ORDERK>
Box prolongate_3d_rf2<CENTI, CENTJ, CENTK, CONSI, CONSJ, CONSK, ORDERI, ORDERJ,
                      ORDERK>::CoarseBox(const Box &fine,
                                         const IntVect &ratio) {
  for (int d = 0; d < dim; ++d)
    assert(fine.type(d) ==
           (indextype()[d] == 0 ? IndexType::NODE : IndexType::CELL));
  for (int d = 0; d < dim; ++d)
    assert(ratio.getVect()[d] == 2);
  Box crse = amrex::coarsen(fine, 2);
  for (int d = 0; d < dim; ++d)
    crse = crse.grow(d, (order()[d] + 1) / 2);
  return crse;
}

template <int CENTI, int CENTJ, int CENTK, bool CONSI, bool CONSJ, bool CONSK,
          int ORDERI, int ORDERJ, int ORDERK>
void prolongate_3d_rf2<CENTI, CENTJ, CENTK, CONSI, CONSJ, CONSK, ORDERI, ORDERJ,
                       ORDERK>::interp(const FArrayBox &crse, int crse_comp,
                                       FArrayBox &fine, int fine_comp,
                                       int ncomp, const Box &fine_region,
                                       const IntVect &ratio,
                                       const Geometry &crse_geom,
                                       const Geometry &fine_geom,
                                       Vector<BCRec> const &bcr,
                                       int actual_comp, int actual_state,
                                       RunOn gpu_or_cpu) {
  for (int d = 0; d < dim; ++d)
    assert(ratio.getVect()[d] == 2);
  // ??? assert(gpu_or_cpu == RunOn::Cpu);

  const BCRec bcrec(BCType::int_dir, BCType::int_dir, BCType::int_dir,
                    BCType::int_dir, BCType::int_dir, BCType::int_dir);
  for (const auto &bc : bcr)
    assert(bc == bcrec);

  assert(actual_comp == 0);  // ???
  assert(actual_state == 0); // ???

  // Target box is intersection of fine_region and domain of fine
  const Box target_region = fine_region & fine.box();
  assert(target_region == fine_region);

  // We prolongate first in the x, then y, then the z direction. Each
  // direction changes the target from corse-plus-ghosts to fine.
  const Box source_region = CoarseBox(target_region, 2);
  array<Box, dim> targets;
  for (int d = 0; d < dim; ++d) {
    targets[d] = d == 0 ? source_region : targets[d - 1];
    targets[d].setRange(d, target_region.loVect()[d], target_region.length(d));
  }
  assert(targets[dim - 1] == target_region);
  array<vector<CCTK_REAL>, dim - 1> tmps;
  for (int d = 0; d < dim - 1; ++d)
    tmps[d].resize(targets[d].numPts());

  for (int comp = 0; comp < ncomp; ++comp) {
    const CCTK_REAL *restrict crseptr = crse.dataPtr(crse_comp + comp);
    CCTK_REAL *restrict fineptr = fine.dataPtr(fine_comp + comp);

#ifdef CCTK_DEBUG
    for (int k = source_region.loVect()[2]; k <= source_region.hiVect()[2];
         ++k) {
      for (int j = source_region.loVect()[1]; j <= source_region.hiVect()[1];
           ++j) {
#pragma omp simd
        for (int i = source_region.loVect()[0]; i <= source_region.hiVect()[0];
             ++i) {
          const IntVect ind(i, j, k);
          assert(crse.box().contains(ind));
          assert(!isnan(crseptr[crse.box().index(ind)]));
        }
      }
    }
#endif
#ifdef CCTK_DEBUG
    for (auto &x : tmps[0])
      x = 0.0 / 0.0;
#endif
    interp3d<CENTI, CONSI, ORDERI, /*D*/ 0>(crseptr, crse.box(), tmps[0].data(),
                                            targets[0], targets[0]);
#ifdef CCTK_DEBUG
    for (auto x : tmps[0])
      assert(!isnan(x));
#endif

#ifdef CCTK_DEBUG
    for (auto &x : tmps[1])
      x = 0.0 / 0.0;
#endif
    interp3d<CENTJ, CONSJ, ORDERJ, /*D*/ 1>(
        tmps[0].data(), targets[0], tmps[1].data(), targets[1], targets[1]);
#ifdef CCTK_DEBUG
    // for (auto x : tmps[1])
    //   assert(!isnan(x));
    for (size_t i = 0; i < tmps[1].size(); ++i)
      assert(!isnan(tmps[1][i]));
#endif

#ifdef CCTK_DEBUG
    for (int k = target_region.loVect()[2]; k <= target_region.hiVect()[2];
         ++k) {
      for (int j = target_region.loVect()[1]; j <= target_region.hiVect()[1];
           ++j) {
#pragma omp simd
        for (int i = target_region.loVect()[0]; i <= target_region.hiVect()[0];
             ++i) {
          const IntVect ind(i, j, k);
          assert(fine.box().contains(ind));
          fineptr[fine.box().index(ind)] = 0.0 / 0.0;
        }
      }
    }
#endif
    interp3d<CENTK, CONSK, ORDERK, /*D*/ 2>(tmps[1].data(), targets[1], fineptr,
                                            fine.box(), target_region);
#ifdef CCTK_DEBUG
    for (int k = target_region.loVect()[2]; k <= target_region.hiVect()[2];
         ++k) {
      for (int j = target_region.loVect()[1]; j <= target_region.hiVect()[1];
           ++j) {
#pragma omp simd
        for (int i = target_region.loVect()[0]; i <= target_region.hiVect()[0];
             ++i) {
          const IntVect ind(i, j, k);
          assert(fine.box().contains(ind));
          assert(!isnan(fineptr[fine.box().index(ind)]));
        }
      }
    }
#endif
  }
}

prolongate_3d_rf2<0, 0, 0, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c000_o1;
prolongate_3d_rf2<0, 0, 1, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c001_o1;
prolongate_3d_rf2<0, 1, 0, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c010_o1;
prolongate_3d_rf2<0, 1, 1, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c011_o1;
prolongate_3d_rf2<1, 0, 0, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c100_o1;
prolongate_3d_rf2<1, 0, 1, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c101_o1;
prolongate_3d_rf2<1, 1, 0, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c110_o1;
prolongate_3d_rf2<1, 1, 1, false, false, false, 1, 1, 1>
    prolongate_3d_rf2_c111_o1;

} // namespace AMReX
