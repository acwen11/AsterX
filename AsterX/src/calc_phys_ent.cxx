#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

#include "setup_eos.hxx"

namespace AsterX {
using namespace Loop;
using namespace EOSX;

enum class eos_3param { IdealGas, Hybrid, Tabulated };

template <typename EOSType>
void AsterX_CalcPhysicalEntropy_typeEoS(CCTK_ARGUMENTS, EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CalcPhysicalEntropy;
  DECLARE_CCTK_PARAMETERS;

  // Calculate "physical" entropy from "evolved" entropy
  // and cycle time levels?
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // phys_ent_p_p(p.I) = phys_ent_p(p.I);
        // phys_ent_p(p.I) = phys_ent(p.I);
        // ent_m2(p.I) = ent_m1(p.I);
        // ent_m1(p.I) = phys_ent(p.I);
        phys_ent(p.I) = eos_3p->physical_from_evolved_ent(entropy(p.I));
      });
}

extern "C" void AsterX_CalcPhysicalEntropy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_CalcPhysicalEntropy;
  DECLARE_CCTK_PARAMETERS;

  // defining EOS objects
  eos_3param eos_3p_type;
  
  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eos_3p_type = eos_3param::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eos_3p_type = eos_3param::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    eos_3p_type = eos_3param::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eos_3p_type) {
    case eos_3param::IdealGas: {
      auto eos_3p_ig = global_eos_3p_ig;

      AsterX_CalcPhysicalEntropy_typeEoS(CCTK_PASS_CTOC, eos_3p_ig);
      break;
    }
    case eos_3param::Hybrid: {
      auto eos_3p_hyb = global_eos_3p_hyb;

      AsterX_CalcPhysicalEntropy_typeEoS(CCTK_PASS_CTOC, eos_3p_hyb);
      break;
    }
    case eos_3param::Tabulated: {
      auto eos_3p_tab3d = global_eos_3p_tab3d;

      AsterX_CalcPhysicalEntropy_typeEoS(CCTK_PASS_CTOC, eos_3p_tab3d);
      break;
    }
    default:
      assert(0);
  }
}
  
} // namespace AsterX
