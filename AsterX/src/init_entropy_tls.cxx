#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

#include "aster_utils.hxx"

namespace AsterX {
using namespace Loop;

extern "C" void AsterX_InitEntropyTLs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_InitEntropyTLs;
  DECLARE_CCTK_PARAMETERS;
  
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
		// Simple copy
		ent_m1(p.I) = entropy(p.I);
		ent_m2(p.I) = entropy(p.I);
      });
}

} // namespace AsterX
