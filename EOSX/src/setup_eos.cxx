#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>

#include <setup_eos.hxx>


namespace EOSX {

using namespace amrex;

enum class eos_id {Polytropic, PWPolytropic};
enum class eos_evol { IdealGas, Hybrid, Tabulated };

// initial data EOS
AMREX_GPU_MANAGED eos_polytrope *eos_poly = nullptr;

// evolution EOS
AMREX_GPU_MANAGED eos_idealgas    *eos_ig    = nullptr;
AMREX_GPU_MANAGED eos_hybrid      *eos_hyb    = nullptr;
AMREX_GPU_MANAGED eos_tabulated3d *eos_tab3d = nullptr;


extern "C" void EOSX_Setup_EOSID(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  eos_id eos_id_type;
  if (CCTK_EQUALS(initial_data_eos, "Polytropic")) {
    eos_id_type = eos_id::Polytropic;
  } else if (CCTK_EQUALS(initial_data_eos, "PWPolytropic")) {
    eos_id_type = eos_id::PWPolytropic;
  } else {
    CCTK_ERROR("Unknown value for parameter \"initial_data_eos\"");
  }

  switch (eos_id_type) {
    case eos_id::Polytropic: {
      CCTK_INFO("Setting initial data EOS to Polytropic");
      eos_poly = (eos_polytrope*)The_Managed_Arena()->alloc(sizeof *eos_poly);
      new (eos_poly) eos_polytrope;
      assert(eos_poly);
      eos_poly->init(poly_gamma, poly_k, rho_max);
      break;
    }
    case eos_id::PWPolytropic: {
      CCTK_ERROR("Piecewise Polytrope EOS is not supported yet!");
      break;
    }
    default:
      assert(0);
  }
}

extern "C" void EOSX_Setup_EOS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  eos_evol eos_evol_type;
  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eos_evol_type = eos_evol::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eos_evol_type = eos_evol::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    eos_evol_type = eos_evol::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eos_evol_type) {
  case eos_evol::IdealGas: {
    CCTK_INFO("Setting evolution EOS to Ideal Gas");
    eos_ig = (eos_idealgas*)The_Managed_Arena()->alloc(sizeof *eos_ig);
    new (eos_ig) eos_idealgas;
    assert(eos_ig);
    eos_ig->init(gl_gamma, particle_mass, rgeps, rgrho, rgye);
    break;
  }
  case eos_evol::Hybrid: {
    CCTK_INFO("Setting evolution EOS to Hybrid");
    eos_hyb = (eos_hybrid*)The_Managed_Arena()->alloc(sizeof *eos_hyb);
    new (eos_hyb) eos_hybrid(eos_poly, gamma_th, rgeps, rgrho, rgye);
    assert(eos_hyb);
    break;
  }
  case eos_evol::Tabulated: {
    CCTK_INFO("Setting evolution EOS to Tabulated3D");
    const string eos_filename = EOSTable_filename;
    eos_tab3d = (eos_tabulated3d*)The_Managed_Arena()->alloc(sizeof *eos_tab3d);
    new (eos_tab3d) eos_tabulated3d;
    assert(eos_tab3d);
    eos_tab3d->init(rgeps, rgrho, rgye);
    eos_tab3d->read_eos_table(eos_filename);
    break;
  }
  default:
    assert(0);
  }

}

} // namespace EOSX

