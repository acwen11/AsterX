#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "setup_eos.hxx"
#include "aster_utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;
using namespace EOSX;

enum class eos_3param { IdealGas, Hybrid, Tabulated };

CCTK_DEVICE CCTK_HOST inline CCTK_REAL 
LOFlagVar(const GF3D2<const CCTK_REAL> &gf, const CCTK_REAL ref, const PointDesc &p) {
    CCTK_REAL etac_sq = 0.0;
    const CCTK_REAL epslo = 1e-22;

    // Calculate derivatives, indicator func in each direction
    for (int dir=0; dir<3; dir++) {
      const CCTK_REAL qimm = gf(p.I - 2*p.DI[dir]);
      const CCTK_REAL qim = gf(p.I - p.DI[dir]);
      const CCTK_REAL qi = gf(p.I);
      const CCTK_REAL qip = gf(p.I + p.DI[dir]);
      const CCTK_REAL qipp = gf(p.I + 2*p.DI[dir]);

      const CCTK_REAL d0Q = abs(ref);
      const CCTK_REAL d1Q = abs(0.5 * (qip - qim));
      const CCTK_REAL d2Q = abs(qip - 2.0 * qi + qim);
      const CCTK_REAL d3Q = abs(0.5 * (qipp - 2.0 * qip + 2.0 * qim - qimm));
      const CCTK_REAL d4Q = abs(qipp - 4.0 * qip + 6.0 * qi - 4.0 * qim + qimm);

      const CCTK_REAL eta_o = d3Q / (d0Q + d1Q + d3Q + epslo);
      const CCTK_REAL eta_e = d4Q / (d0Q + d2Q + d4Q + epslo);
      const CCTK_REAL eta_ci = max(eta_o, eta_e);

      etac_sq += eta_ci * eta_ci;
    }

    return sqrt(etac_sq);  
}

// Calculate low-order flag for a particular gridfunction
template<typename EOSType>
void CalcLOFlag(CCTK_ARGUMENTS, EOSType *eos_3p, const bool havetemp) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetLOFlag;
  DECLARE_CCTK_PARAMETERS;

  const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_vels{velx, vely, velz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_Bvecs{Bvecx, Bvecy, Bvecz};

  // Loop over the grid
  grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
  
    if (rho(p.I) < LO_rhothresh) {
      LOflag(p.I) = 1.0;
      return;
    }

    // Store largest etac as diagnostic. Note that this is only truly the largest
    // if all checks pass.
    CCTK_REAL etac_tot = 0.0;

    // Check density
    CCTK_REAL etacL = LOFlagVar(rho, rho(p.I), p);
    if (etacL > etac_tot)
      etac_tot = etacL;
    if (etacL > eta_thresh){
      etac(p.I) = etac_tot;
      LOflag(p.I) = 1.0;
      return;
    }

    // Check pressure
    etacL = LOFlagVar(press, press(p.I), p);
    if (etacL > etac_tot)
      etac_tot = etacL;
    if (etacL > eta_thresh){
      etac(p.I) = etac_tot;
      LOflag(p.I) = 1.0;
      return;
    }

    // Calculate c_sound
    CCTK_REAL csL = 0.0;
    if (havetemp) {
      csL = eos_3p->csnd_from_valid_rho_temp_ye(rho(p.I), temperature(p.I), Ye(p.I));
    } else {
      CCTK_REAL epsL = eps(p.I);
      csL = eos_3p->csnd_from_valid_rho_eps_ye(rho(p.I), epsL, Ye(p.I));
    }
    const CCTK_REAL cs = csL;

    // Check velocity 
    for (int dir=0; dir<3; dir++) {
      etacL = LOFlagVar(gf_vels(dir), cs, p);
      if (etacL > etac_tot)
        etac_tot = etacL;
      if (etacL > eta_thresh){
        etac(p.I) = etac_tot;
        LOflag(p.I) = 1.0;
        return;
      }
    }
      
    /* Disable this for now
    // Calculate |B|
    const vec<CCTK_REAL, 3> BvecsL{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
    const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
      return calc_avg_v2c(gf_g(i, j), p);
    });
    const CCTK_REAL normBL = calc_norm(BvecsL, g_avg);
    
    // Check B^i's
    for (int dir=0; dir<3; dir++) {
      etacL = LOFlagVar(gf_Bvecs(dir), normBL, p);
      if (etacL > etac_tot)
        etac_tot = etacL;
      if (etacL > eta_thresh){
        etac(p.I) = etac_tot;
        LOflag(p.I) = 0.0;
        return;
      }
    }
    */

    // All checks pass
    etac(p.I) = etac_tot;
    LOflag(p.I) = 0.0;

  });
}

extern "C" void AsterX_SetLOFlag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetLOFlag;
  DECLARE_CCTK_PARAMETERS;

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
    // Get local eos object
    auto eos_3p_ig = global_eos_3p_ig;
    CalcLOFlag(cctkGH, eos_3p_ig, false);
    break;
  }
  case eos_3param::Hybrid: {
    auto eos_3p_hyb = global_eos_3p_hyb;
    CalcLOFlag(cctkGH, eos_3p_hyb, false);
    break;
  }
  case eos_3param::Tabulated: {
    auto eos_3p_tab3d = global_eos_3p_tab3d;
    CalcLOFlag(cctkGH, eos_3p_tab3d, true);
    break;
  }
  default:
    assert(0);
  }
}

extern "C" void AsterX_InitLOFlag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_InitLOFlag;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the grid
  grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
		LOflag(p.I) = 0.0;
	});
}

} // namespace AsterX
