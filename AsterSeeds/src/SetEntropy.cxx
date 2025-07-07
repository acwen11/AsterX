#include <loop.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdio>
#include <cstdbool>
#include <cmath>

#include <setup_eos.hxx>

using namespace EOSX;
enum class eos_3param { IdealGas, Hybrid, Tabulated };

template <typename EOSType>
void SetEntropy_typeEoS(CCTK_ARGUMENTS, EOSType *eos_3p) {

  DECLARE_CCTK_ARGUMENTSX_SetEntropy;
  DECLARE_CCTK_PARAMETERS;

  if (set_entropy_postinitial) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          entropy(p.I) =
              eos_3p->kappa_from_valid_rho_eps_ye(rho(p.I), eps(p.I), Ye(p.I));
        });
  }
}

extern "C" void SetEntropy(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTSX_SetEntropy;
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
    SetEntropy_typeEoS(CCTK_PASS_CTOC, eos_3p_ig);
    break;
  }
  case eos_3param::Hybrid: {
    auto eos_3p_hyb = global_eos_3p_hyb;
    SetEntropy_typeEoS(CCTK_PASS_CTOC, eos_3p_hyb);
    break;
  }
  case eos_3param::Tabulated: {
    auto eos_3p_tab3d = global_eos_3p_tab3d;
    SetEntropy_typeEoS(CCTK_PASS_CTOC, eos_3p_tab3d);
    break;
  }
  default:
    assert(0);
  }
}
