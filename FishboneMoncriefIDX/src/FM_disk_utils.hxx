#ifndef FM_UTILS_HXX
#define FM_UTILS_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "aster_utils.hxx"

namespace FMdisk {

using namespace std;

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(T x) {
  return x * x;
}

enum class atmosphere_t {isentropic_graded,
	                 free_graded,
		         constant };

} // namespace FMdisk

#endif // FM_UTILS_HXX
