#include <array>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include "aster_utils.hxx"

namespace AsterX {
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;

enum class vector_potential_gauge_t { algebraic, generalized_lorenz };

template <int i, vector_potential_gauge_t gauge>
void CalcRHSofAvec_impl(CCTK_ARGUMENTS, const int order) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_RHS;

  const vec<GF3D2<const CCTK_REAL>, dim> gf_E{Ex, Ey, Ez};
  const vec<GF3D2<CCTK_REAL>, dim> gf_Avec_rhs{Avec_x_rhs, Avec_y_rhs,
                                               Avec_z_rhs};

  if constexpr (gauge == vector_potential_gauge_t::algebraic) {

    grid.loop_int_device<i == 0, i == 1, i == 2>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          gf_Avec_rhs(i)(p.I) = -gf_E(i)(p.I);
        });

  } else if constexpr (gauge == vector_potential_gauge_t::generalized_lorenz) {

    grid.loop_int_device<i == 0, i == 1, i == 2>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          gf_Avec_rhs(i)(p.I) = -gf_E(i)(p.I) - calc_fd2_v2e<i>(G, p, order);
        });
  }
}

template <vector_potential_gauge_t gauge>
void CalcRHSofPsi_impl(CCTK_ARGUMENTS, const CCTK_REAL damp_fac) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_RHS;

  const vec<GF3D2<const CCTK_REAL>, dim> gf_Fstag{Fx_stag, Fy_stag, Fz_stag};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_Fbeta{Fbetax, Fbetay, Fbetaz};

  if constexpr (gauge == vector_potential_gauge_t::algebraic) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Psi_rhs(p.I) = 0.0; });
  } else if constexpr (gauge == vector_potential_gauge_t::generalized_lorenz) {
    grid.loop_int_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL dF = 0.0;
          for (int i = 0; i < dim; i++) {
            dF += calc_fd2_e2v(gf_Fstag(i), p, i) -
                  (gf_beta(i)(p.I) < 0
                       ? calc_fd2_v2v_oneside<-1>(gf_Fbeta(i), p, i)
                       : calc_fd2_v2v_oneside<+1>(gf_Fbeta(i), p, i));
          }
          Psi_rhs(p.I) = -dF - damp_fac * alp(p.I) * Psi(p.I);
        });
  }
}

template <int i>
void CalcRHSofAvec(CCTK_ARGUMENTS, const vector_potential_gauge_t gauge, const int order) {
  switch (gauge) {
  case vector_potential_gauge_t::algebraic: {
    CalcRHSofAvec_impl<i, vector_potential_gauge_t::algebraic>(CCTK_PASS_CTOC);
    break;
  }
  case vector_potential_gauge_t::generalized_lorenz: {
    CalcRHSofAvec_impl<i, vector_potential_gauge_t::generalized_lorenz>(
        CCTK_PASS_CTOC, order);
    break;
  }
  default:
    assert(0);
  }
}

void CalcRHSofPsi(CCTK_ARGUMENTS, const vector_potential_gauge_t gauge,
                  const CCTK_REAL damp_fac) {
  switch (gauge) {
  case vector_potential_gauge_t::algebraic: {
    CalcRHSofPsi_impl<vector_potential_gauge_t::algebraic>(CCTK_PASS_CTOC,
                                                           damp_fac);
    break;
  }
  case vector_potential_gauge_t::generalized_lorenz: {
    CalcRHSofPsi_impl<vector_potential_gauge_t::generalized_lorenz>(
        CCTK_PASS_CTOC, damp_fac);
  } break;
  default:
    assert(0);
  }
}

extern "C" void AsterX_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_RHS;
  DECLARE_CCTK_PARAMETERS;

  vector_potential_gauge_t gauge;
  if (CCTK_EQUALS(vector_potential_gauge, "algebraic"))
    gauge = vector_potential_gauge_t::algebraic;
  else if (CCTK_EQUALS(vector_potential_gauge, "generalized Lorenz"))
    gauge = vector_potential_gauge_t::generalized_lorenz;
  else
    CCTK_ERROR("Unknown value for parameter \"vector_potential_gauge\"");

  const vec<CCTK_REAL, dim> idx{1 / CCTK_DELTA_SPACE(0),
                                1 / CCTK_DELTA_SPACE(1),
                                1 / CCTK_DELTA_SPACE(2)};

  const vec<GF3D2<const CCTK_REAL>, dim> gf_fdens{fxdens, fydens, fzdens};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_fDEnt{fxDEnt, fyDEnt, fzDEnt};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_fmomx{fxmomx, fymomx, fzmomx};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_fmomy{fxmomy, fymomy, fzmomy};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_fmomz{fxmomz, fymomz, fzmomz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_ftau{fxtau, fytau, fztau};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_fDYe{fxDYe, fyDYe, fzDYe};

  const auto calcupdate_hydro =
      [=] CCTK_DEVICE(const vec<GF3D2<const CCTK_REAL>, dim> &gf_fluxes,
                      const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        vec<CCTK_REAL, 3> dfluxes([&](int i) ARITH_INLINE {
          return compute_flux_derivatives(gf_fluxes(i), p, i, correction_order);
        });
        return -calc_contraction(idx, dfluxes);
      };

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        densrhs(p.I) += calcupdate_hydro(gf_fdens, p);
        DEntrhs(p.I) += calcupdate_hydro(gf_fDEnt, p);
        momxrhs(p.I) += calcupdate_hydro(gf_fmomx, p);
        momyrhs(p.I) += calcupdate_hydro(gf_fmomy, p);
        momzrhs(p.I) += calcupdate_hydro(gf_fmomz, p);
        taurhs(p.I) += calcupdate_hydro(gf_ftau, p);
        DYe_rhs(p.I) += calcupdate_hydro(gf_fDYe, p);

        // Diagnostic only, save min(theta)
        theta_tot(p.I) = 1.0;
        for (int ii=0; ii<3; ii++)
          theta_tot(p.I) = min({theta_tot(p.I),
              theta_x(p.I), theta_x(p.I + p.DI[ii]), theta_y(p.I), theta_y(p.I + p.DI[ii]),
              theta_z(p.I), theta_z(p.I + p.DI[ii])});

#ifdef CCTK_DEBUG
        if (isnan(densrhs(p.I))) {
          printf("calcupdate = %f, ", calcupdate_hydro(gf_fdens, p));
          printf("densrhs = %f, gf_fdens = %f, %f, %f, %f, %f, %f \n",
                 densrhs(p.I), gf_fdens(0)(p.I), gf_fdens(1)(p.I),
                 gf_fdens(2)(p.I), gf_fdens(0)(p.I + p.DI[0]),
                 gf_fdens(1)(p.I + p.DI[1]), gf_fdens(2)(p.I + p.DI[2]));
        }
        assert(!isnan(densrhs(p.I)));
#endif
      });

  CalcRHSofAvec<0>(CCTK_PASS_CTOC, gauge);
  CalcRHSofAvec<1>(CCTK_PASS_CTOC, gauge);
  CalcRHSofAvec<2>(CCTK_PASS_CTOC, gauge);

  CalcRHSofPsi(CCTK_PASS_CTOC, gauge, lorenz_damp_fac);
}

} // namespace AsterX
