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

// Calculate low-order flag for a particular gridfunction
CCTK_DEVICE CCTK_HOST inline CCTK_REAL
LOFlagVar(const GF3D2<const CCTK_REAL> &gf, const CCTK_REAL ref,
          const PointDesc &p) {
  CCTK_REAL etac_sq = 0.0;
  const CCTK_REAL epslo = 1e-22;

  // Calculate derivatives, indicator func in each direction
  for (int dir = 0; dir < 3; dir++) {
    const CCTK_REAL qimm = gf(p.I - 2 * p.DI[dir]);
    const CCTK_REAL qim = gf(p.I - p.DI[dir]);
    const CCTK_REAL qi = gf(p.I);
    const CCTK_REAL qip = gf(p.I + p.DI[dir]);
    const CCTK_REAL qipp = gf(p.I + 2 * p.DI[dir]);

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

template <typename EOSType> void CalcLOFlag(CCTK_ARGUMENTS, EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetLOeta;
  DECLARE_CCTK_PARAMETERS;

  const vec<GF3D2<const CCTK_REAL>, dim> gf_vels{velx, vely, velz};
  // const vec<GF3D2<const CCTK_REAL>, dim> gf_Bvecs{Bvecx, Bvecy, Bvecz};
  // const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  // Loop over the grid
  grid.loop_allmn_device<1, 1, 1>(
      grid.nghostzones, 2,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        // Grading rho
        const CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        CCTK_REAL rho_atm = (radial_distance > r_atmo)
                ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
                : rho_abs_min;
        rho_atm = std::max(eos_3p->rgrho.min, rho_atm);
        if (rho(p.I) <= rho_atm * (1.0 + atmo_tol)) {
          etac(p.I) = eta_thresh * 1.001;
          return;
        }

        // Store largest etac as diagnostic. Note that this is only truly the
        // largest if all checks pass.
        CCTK_REAL etac_tot = 0.0;

        // Check density
        CCTK_REAL etacL = LOFlagVar(rho, rho(p.I), p);
        if (etacL > etac_tot) {
          etac_tot = etacL;
          if (etac_tot > eta_thresh) {
            etac(p.I) = etac_tot;
            return;
          }
        }

        // Check pressure
        etacL = LOFlagVar(press, press(p.I), p);
        if (etacL > etac_tot) {
          etac_tot = etacL;
          if (etac_tot > eta_thresh) {
            etac(p.I) = etac_tot;
            return;
          }
        }

        // Calculate c_sound
        const CCTK_REAL cs =
            eos_3p->csnd_from_rho_temp_ye(rho(p.I), temperature(p.I), Ye(p.I));

        // Check velocity
        for (int dir = 0; dir < 3; dir++) {
          etacL = LOFlagVar(gf_vels(dir), cs, p);
          if (etacL > etac_tot) {
            etac_tot = etacL;
            if (etac_tot > eta_thresh) {
              etac(p.I) = etac_tot;
              return;
            }
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
            if (etac_tot > eta_thresh) {
              etac(p.I) = etac_tot;
              return;
            }
        }
        */

        // All checks pass
        etac(p.I) = etac_tot;
      });

    grid.loop_outer_n_device<1, 1, 1>(
        grid.nghostzones, 2,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        etac(p.I) = 0.0; // set eta_c = 0 here to mark as valid
      });
}

extern "C" void AsterX_SetLOeta(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetLOeta;
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
    CalcLOFlag(cctkGH, eos_3p_ig);
    break;
  }
  case eos_3param::Hybrid: {
    if (global_eos_3p_hyb_pwpoly) {
      // Get local eos object
      auto eos_3p_hyb = global_eos_3p_hyb_pwpoly;
      CalcLOFlag(cctkGH, eos_3p_hyb);

    } else if (global_eos_3p_hyb_poly) {
      // Get local eos object
      auto eos_3p_hyb = global_eos_3p_hyb_poly;
      CalcLOFlag(cctkGH, eos_3p_hyb);

    } else {
      CCTK_ERROR(
          "Hybrid EOS selected but no hybrid EOS object was initialized");
    }

    break;
  }
  case eos_3param::Tabulated: {
    auto eos_3p_tab3d = global_eos_3p_tab3d;
    CalcLOFlag(cctkGH, eos_3p_tab3d);
    break;
  }
  default:
    assert(0);
  }
}

extern "C" void AsterX_SetLOFlag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SetLOFlag;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the grid
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { 
        if (rho(p.I) < loworder_rho_thresh) {
          LOflag(p.I) = 1.0;
          return;
        }
        // Check nearest neighbors
        bool flag = false;
        for (int kk=-1; kk<2; kk++) { 
          for (int jj=-1; jj<2; jj++) { 
            for (int ii=-1; ii<2; ii++) { 
              const auto II = p.I + ii * p.DI[0] + jj * p.DI[1] + kk * p.DI[2];
              flag = flag || (etac(II) > eta_thresh);
            }
          }
        }
        LOflag(p.I) = flag;
  });
}

extern "C" void AsterX_InitLOFlag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_InitLOFlag;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { LOflag(p.I) = 0.0; });
}

extern "C" void AsterX_LOFlagCopyTLs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_LOFlagCopyTLs;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { LOflag_p(p.I) = LOflag(p.I); });
}

extern "C" void AsterX_LOFlagCopyCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_LOFlagCopyCheck;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { flagcheck(p.I) = LOflag_p(p.I); });
}

// extern "C" void AsterX_AdjustLOFlagPostRegrid(CCTK_ARGUMENTS) {
//   DECLARE_CCTK_ARGUMENTSX_AsterX_AdjustLOFlagPostRegrid;
//   DECLARE_CCTK_PARAMETERS;
// 
//   grid.loop_all_device<1, 1, 1>(
//         grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
//             CCTK_ATTRIBUTE_ALWAYS_INLINE {
//           LOflag(p.I) = (LOflag(p.I) > 0.5) ? 1.0 : 0.0;
//           LOflag_p(p.I) = LOflag(p.I);
//         });
// }

} // namespace AsterX
