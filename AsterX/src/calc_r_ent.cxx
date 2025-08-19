#include <loop_device.hxx>

#include <derivs.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>

#include "aster_utils.hxx"
#include "setup_eos.hxx"

namespace AsterX {
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;
using namespace EOSX;

enum class eos_3param { IdealGas, Hybrid, Tabulated };

template <typename EOSType>
void AsterX_CalcEntropyResidual_typeEoS(CCTK_ARGUMENTS, EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CalcEntropyResidual;
  DECLARE_CCTK_PARAMETERS;

  // Compute d_i (phys entropy)
  // Derivs boilerplate
  Arith::vect<int, dim> imin, imax;
  const std::array<int, dim> nghostzones = {
      cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]};
  GridDescBase(cctkGH).box_int<1, 1, 1>(nghostzones, imin, imax);
  const GF3D5layout layout5(imin, imax);
  const Arith::vect<CCTK_REAL, dim> dx(std::array<CCTK_REAL, dim>{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  });

  constexpr int nvars = 4;
  int ivar = 0;
  GF3D5vector<CCTK_REAL> vars(layout5, nvars);
  const auto make_gf = [&]() { return GF3D5<CCTK_REAL>(vars(ivar++)); };
  const auto make_vec = [&](const auto &f) {
    return Arith::vec<std::result_of_t<decltype(f)()>, dim>(
        [&](int) { return f(); });
  };
  const auto make_vec_gf = [&]() { return make_vec(make_gf); };

  const GF3D5<CCTK_REAL> t5_s(make_gf());
  const Arith::vec<GF3D5<CCTK_REAL>, dim> t5_ds(make_vec_gf());

  Derivs::calc_derivs<1, 1, 1>(t5_s, t5_ds, layout5, grid, phys_ent, dx,
                               efl_deriv_order);

  // Prep spacetime GFs
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        // Calculate spatial contraction
        const GF3D5index index5(layout5, p.I);
        const vec<CCTK_REAL, 3> di_s = t5_ds(index5);

        
        const CCTK_REAL alp_avg = calc_avg_v2c(alp, p);
        const vec<CCTK_REAL, 3> betas_avg(
            [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });
        const vec<CCTK_REAL, 3> vels{velx(p.I), vely(p.I), velx(p.I)};

        const CCTK_REAL v_dis =
            calc_contraction(alp_avg * vels - betas_avg, di_s);

        // Calculate d_t s
        const CCTK_REAL i2dt = 1 / (2 * cctk_delta_time);
        const CCTK_REAL dts =
            i2dt * (3 * phys_ent(p.I) - 4 * ent_m1(p.I) + ent_m2(p.I));
        
        // Calculate R
        r_ent(p.I) = std::abs(dts + v_dis);
      });
}

extern "C" void AsterX_CalcEntropyResidual(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_CalcEntropyResidual;
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

      AsterX_CalcEntropyResidual_typeEoS(CCTK_PASS_CTOC, eos_3p_ig);
      break;
    }
    case eos_3param::Hybrid: {
      auto eos_3p_hyb = global_eos_3p_hyb;

      AsterX_CalcEntropyResidual_typeEoS(CCTK_PASS_CTOC, eos_3p_hyb);
      break;
    }
    case eos_3param::Tabulated: {
      auto eos_3p_tab3d = global_eos_3p_tab3d;

      AsterX_CalcEntropyResidual_typeEoS(CCTK_PASS_CTOC, eos_3p_tab3d);
      break;
    }
    default:
      assert(0);
  }
}
  
} // namespace AsterX
