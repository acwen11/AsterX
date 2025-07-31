#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>

#include <setup_eos.hxx>

namespace EOSX {

using namespace amrex;

enum class eos_1param { Polytropic, PWPolytropic };
enum class eos_3param { IdealGas, Hybrid, Tabulated };

// initial data EOS
eos_1p_polytropic *global_eos_1p_poly = nullptr;

// evolution EOS
eos_3p_idealgas *global_eos_3p_ig = nullptr;
eos_3p_hybrid *global_eos_3p_hyb = nullptr;
eos_3p_tabulated3d *global_eos_3p_tab3d = nullptr;

extern "C" void EOSX_Setup_EOSID(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  eos_1param eos_1p_type;
  if (CCTK_EQUALS(initial_data_eos, "Polytropic")) {
    eos_1p_type = eos_1param::Polytropic;
  } else if (CCTK_EQUALS(initial_data_eos, "PWPolytropic")) {
    eos_1p_type = eos_1param::PWPolytropic;
  } else {
    CCTK_ERROR("Unknown value for parameter \"initial_data_eos\"");
  }

  switch (eos_1p_type) {
  case eos_1param::Polytropic: {
    CCTK_INFO("Setting initial data EOS to Polytropic");
    global_eos_1p_poly = (eos_1p_polytropic *)The_Managed_Arena()->alloc(
        sizeof *global_eos_1p_poly);
    assert(global_eos_1p_poly);
    new (global_eos_1p_poly) eos_1p_polytropic;
    global_eos_1p_poly->init(poly_gamma, poly_k, rho_max);
    break;
  }
  case eos_1param::PWPolytropic: {
    CCTK_ERROR("Piecewise Polytrope EOS is not supported yet!");
    break;
  }
  default:
    assert(0);
  }
}

extern "C" void EOSX_Setup_EOS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  eos_3param eos_3p_type;
  eos_3p::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

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
    CCTK_INFO("Setting evolution EOS to Ideal Gas");
    global_eos_3p_ig =
        (eos_3p_idealgas *)The_Managed_Arena()->alloc(sizeof *global_eos_3p_ig);
    assert(global_eos_3p_ig);
    new (global_eos_3p_ig) eos_3p_idealgas;
    global_eos_3p_ig->init(gl_gamma, particle_mass, rgeps, rgrho, rgye);
    break;
  }
  case eos_3param::Hybrid: {
    CCTK_INFO("Setting evolution EOS to Hybrid");
    global_eos_3p_hyb =
        (eos_3p_hybrid *)The_Managed_Arena()->alloc(sizeof *global_eos_3p_hyb);
    assert(global_eos_3p_hyb);
    new (global_eos_3p_hyb)
        eos_3p_hybrid(global_eos_1p_poly, gamma_th, rgeps, rgrho, rgye);
    break;
  }
  case eos_3param::Tabulated: {
    CCTK_INFO("Setting evolution EOS to Tabulated3D");
    const string eos_filename = EOSTable_filename;
    global_eos_3p_tab3d = (eos_3p_tabulated3d *)The_Managed_Arena()->alloc(
        sizeof *global_eos_3p_tab3d);
    assert(global_eos_3p_tab3d);
    new (global_eos_3p_tab3d) eos_3p_tabulated3d;
    global_eos_3p_tab3d->init(eos_filename, rgeps, rgrho, rgye);
    break;
  }
  default:
    assert(0);
  }
}

} // namespace EOSX
