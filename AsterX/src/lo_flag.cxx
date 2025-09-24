#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

/*
#include <mat.hxx>
#include <simd.hxx>
#include <sum.hxx>
#include <vec.hxx>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
*/

#include "setup_eos.hxx"
#include "aster_utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;
using namespace EOSX;

enum class eos_3param { IdealGas, Hybrid, Tabulated };

/*
// Calculate sound speed
template <typename EOSType>
void cs_LOAuxGF(CCTK_ARGUMENTS, EOSType *eos_3p, bool havetemp) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetLOFlag;
  //DECLARE_CCTK_PARAMETERS;

  // Loop over the entire grid (0 to n-1 cells in each direction)
  grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    if (havetemp) {
      lo_auxgf(p.I) = eos_3p->csnd_from_valid_rho_temp_ye(rho(p.I), temperature(p.I), Ye(p.I));
    }
    else {
      lo_auxgf(p.I) = eos_3p->csnd_from_valid_rho_eps_ye(rho(p.I), eps(p.I), Ye(p.I));
    }
  });
}
*/

// Calculate low-order flag for a particular gridfunction
template<typename T>
void CalcLOFlag(CCTK_ARGUMENTS, const GF3D2<const CCTK_REAL> &gf, const GF3D2<T> &refgf) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetLOFlag;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the entire grid (0 to n-1 cells in each direction)
  grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
  
    if (LOflag(p.I) == 0)
      return;

    CCTK_REAL etac_sq = 0;
    const CCTK_REAL epslo = 1e-22;

    for (int dir=0; dir<3; dir++) {
      const CCTK_REAL qimm = gf(p.I - 2*p.DI[dir]);
      const CCTK_REAL qim = gf(p.I - p.DI[dir]);
      const CCTK_REAL qi = gf(p.I);
      const CCTK_REAL qip = gf(p.I + p.DI[dir]);
      const CCTK_REAL qipp = gf(p.I + 2*p.DI[dir]);

      const CCTK_REAL d0Q = abs(refgf(p.I));
      const CCTK_REAL d1Q = abs(0.5 * (qip - qim));
      const CCTK_REAL d2Q = abs(qip - 2.0 * qi + qim);
      const CCTK_REAL d3Q = abs(0.5 * (qipp - 2.0 * qip + 2.0 * qim - qimm));
      const CCTK_REAL d4Q = abs(qipp - 4.0 * qip + 6.0 * qi - 4.0 * qim + qimm);

      const CCTK_REAL eta_o = d3Q / (d0Q + d1Q + d3Q + epslo);
      const CCTK_REAL eta_e = d4Q / (d0Q + d2Q + d4Q + epslo);
      const CCTK_REAL eta_ci = max(eta_o, eta_e);

      etac_sq += eta_ci * eta_ci;
    }

    const CCTK_REAL etacL = sqrt(etac_sq);  
    if (etacL > etac(p.I))
      etac(p.I) = etacL;
    
    LOflag(p.I) *= etacL < eta_thresh;
  });
}

extern "C" void AsterX_SetLOFlag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetLOFlag;
  DECLARE_CCTK_PARAMETERS;

  // Initialize flags
  grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    etac(p.I) = 0.0;
    LOflag(p.I) = 1;
  });

  // Detect shocks in rho, press
  CalcLOFlag(cctkGH, rho, rho);
  CalcLOFlag(cctkGH, press, press);

  // Set reference gridfunction to c_s
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
    grid.loop_int_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      CCTK_REAL epsL = eps(p.I);
      lo_auxgf(p.I) = eos_3p_ig->csnd_from_valid_rho_eps_ye(rho(p.I), epsL, Ye(p.I));
    });
    break;
  }
  case eos_3param::Hybrid: {
    auto eos_3p_hyb = global_eos_3p_hyb;
    grid.loop_int_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      CCTK_REAL epsL = eps(p.I);
      lo_auxgf(p.I) = eos_3p_hyb->csnd_from_valid_rho_eps_ye(rho(p.I), epsL, Ye(p.I));
    });
    break;
  }
  case eos_3param::Tabulated: {
    auto eos_3p_tab3d = global_eos_3p_tab3d;
    grid.loop_int_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      lo_auxgf(p.I) = eos_3p_tab3d->csnd_from_valid_rho_temp_ye(rho(p.I), temperature(p.I), Ye(p.I));
    });
    break;
  }
  default:
    assert(0);
  }

  // Detect shocks in velocity
  CalcLOFlag(cctkGH, velx, lo_auxgf);
  CalcLOFlag(cctkGH, vely, lo_auxgf);
  CalcLOFlag(cctkGH, velz, lo_auxgf);

  // Set reference gridfunction to |B|
  const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

    const vec<CCTK_REAL, 3> Bvecs{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
    const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
      return calc_avg_v2c(gf_g(i, j), p);
    });

    lo_auxgf(p.I) = calc_norm(Bvecs, g_avg);
  });

  // Detect shocks in Bvec 
  CalcLOFlag(cctkGH, Bvecx, lo_auxgf);
  CalcLOFlag(cctkGH, Bvecy, lo_auxgf);
  CalcLOFlag(cctkGH, Bvecz, lo_auxgf);
}

} // namespace AsterX
