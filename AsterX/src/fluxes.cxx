#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <mat.hxx>
#include <simd.hxx>
#include <sum.hxx>
#include <vec.hxx>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include "aster_utils.hxx"
#include "eigenvalues.hxx"
#include "fluxes.hxx"
#include "reconstruct.hxx"
#include "setup_eos.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;
using namespace EOSX;
using namespace ReconX;
using namespace AsterUtils;

enum class flux_t { LxF, HLLE };
enum class eos_3param { IdealGas, Hybrid, Tabulated };
enum class rec_var_t { v_vec, z_vec, s_vec };

template <int dir, typename EOSType>
void CalcConsFluxesFromPrims(auto eos_3p, const bool use_temp_flag,
 CCTK_REAL alp_avg, vec<CCTK_REAL, 3> betas_avg, smat<CCTK_REAL, 3> g_avg, CCTK_REAL sqrtg,
 recon_prims& rcprims, recon_cons& rccons) {

    /* vtilde^i = alpha * v^i - beta^i */
    const vec<vec<CCTK_REAL, 2>, 3> vtildes_rc([&](int i) ARITH_INLINE {
      return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
        return alp_avg * rcprims.vels_rc(i)(f) - betas_avg(i);
      });
    });

    /* alpha * b0 = W * B^i * v_i */
    const vec<CCTK_REAL, 2> alp_b0_rc([&](int f) ARITH_INLINE {
      return rcprims.w_lorentz_rc(f) * calc_contraction(rcprims.Bs_rc, rcprims.vlows_rc)(f);
    });
    /* covariant magnetic field measured by the Eulerian observer */
    const vec<vec<CCTK_REAL, 2>, 3> Blows_rc = calc_contraction(g_avg, rcprims.Bs_rc);
    /* B^2 = B^i * B_i */
    const vec<CCTK_REAL, 2> B2_rc = calc_contraction(rcprims.Bs_rc, Blows_rc);
    /* covariant magnetic field measured by the comoving observer:
     *  b_i = B_i/W + alpha*b^0*v_i */
    const vec<vec<CCTK_REAL, 2>, 3> blows_rc([&](int i) ARITH_INLINE {
      return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
        return Blows_rc(i)(f) / rcprims.w_lorentz_rc(f) + alp_b0_rc(f) * rcprims.vlows_rc(i)(f);
      });
    });
    /* b^2 = b^{\mu} * b_{\mu} */
    const vec<CCTK_REAL, 2> bsq_rc([&](int f) ARITH_INLINE {
      return (B2_rc(f) + pow2(alp_b0_rc(f))) / pow2(rcprims.w_lorentz_rc(f));
    });

    /* componets correspond to the dir we are considering */
    const CCTK_REAL beta_avg = betas_avg(dir);
    const vec<CCTK_REAL, 2> vel_rc{rcprims.vels_rc(dir)};
    const vec<CCTK_REAL, 2> B_rc{rcprims.Bs_rc(dir)};
    const vec<CCTK_REAL, 2> vtilde_rc{vtildes_rc(dir)};

    vec<CCTK_REAL, 2> cs2_rc;
    if (use_temp_flag) {
      const vec<CCTK_REAL, 2> cs_dummy([&](int f) ARITH_INLINE {
        return eos_3p->csnd_from_valid_rho_temp_ye(rcprims.rho_rc(f), rcprims.eps_rc(f),
                                                  rcprims.Ye_rc(f));
      });
      for (int i=0; i<=1; i++)
        cs2_rc(i) = cs_dummy(i) * cs_dummy(i);
    }
    else {
      const vec<CCTK_REAL, 2> cs2_rc([&](int f) ARITH_INLINE {
        return eos_3p->csnd_from_valid_rho_eps_ye(rcprims.rho_rc(f), rcprims.eps_rc(f),
                                                  rcprims.Ye_rc(f)) *
               eos_3p->csnd_from_valid_rho_eps_ye(rcprims.rho_rc(f), rcprims.eps_rc(f), rcprims.Ye_rc(f));
      });
    }

    const vec<CCTK_REAL, 2> h_rc([&](int f) ARITH_INLINE {
      return 1 + rcprims.eps_rc(f) + rcprims.press_rc(f) / rcprims.rho_rc(f);
    });

    /* Computing conservatives from primitives: */

    /* dens = sqrt(g) * D = sqrt(g) * (rho * W) */
    const vec<CCTK_REAL, 2> dens_rc([&](int f) ARITH_INLINE {
      return sqrtg * rcprims.rho_rc(f) * rcprims.w_lorentz_rc(f);
    });

    /* DEnt = sqrt(g) * D * s  = sqrt(g) * (rho * W) * s */
    /*    s = entropy */
    const vec<CCTK_REAL, 2> DEnt_rc([&](int f) ARITH_INLINE {
      return sqrtg * rcprims.rho_rc(f) * rcprims.w_lorentz_rc(f) * rcprims.entropy_rc(f);
    });

    /* auxiliary: dens * h * W = sqrt(g) * rho * h * W^2 */
    const vec<CCTK_REAL, 2> dens_h_W_rc([&](int f) ARITH_INLINE {
      return dens_rc(f) * h_rc(f) * rcprims.w_lorentz_rc(f);
    });
    /* auxiliary: sqrt(g) * (rho*h + b^2)*W^2 */
    const vec<CCTK_REAL, 2> dens_h_W_plus_sqrtg_W2b2_rc =
        dens_h_W_rc + sqrtg * (pow2(alp_b0_rc) + B2_rc);
    /* auxiliary: (pgas + pmag) */
    const vec<CCTK_REAL, 2> press_plus_pmag_rc = rcprims.press_rc + 0.5 * bsq_rc;

    /* mom_i = sqrt(g)*S_i = sqrt(g)((rho*h+b^2)*W^2*v_i - alpha*b^0*b_i) */
    const vec<vec<CCTK_REAL, 2>, 3> moms_rc([&](int i) ARITH_INLINE {
      return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
        return dens_h_W_plus_sqrtg_W2b2_rc(f) * rcprims.vlows_rc(i)(f) -
               sqrtg * alp_b0_rc(f) * blows_rc(i)(f);
      });
    });

    /* tau = sqrt(g)*t =
     *  sqrt(g)((rho*h + b^2)*W^2 - (pgas+pmag) - (alpha*b^0)^2 - D) */
    const vec<CCTK_REAL, 2> tau_rc =
        dens_h_W_rc - dens_rc + sqrtg * (B2_rc - press_plus_pmag_rc);

    /* Btildes^i = sqrt(g) * B^i */
    const vec<vec<CCTK_REAL, 2>, 3> Btildes_rc(
        [&](int i) ARITH_INLINE { return sqrtg * rcprims.Bs_rc(i); });

    /* Computing fluxes of conserved variables: */

    /* auxiliary: unit in 'dir' */
    const vec<CCTK_REAL, 3> unit_dir{vec<int, 3>::unit(dir)};
    /* auxiliary: alpha * sqrt(g) */
    const CCTK_REAL alp_sqrtg = alp_avg * sqrtg;
    /* auxiliary: B^i / W */
    const vec<CCTK_REAL, 2> B_over_w_lorentz_rc(
        [&](int f) ARITH_INLINE { return B_rc(f) / rcprims.w_lorentz_rc(f); });

    /* flux(dens) = sqrt(g) * D * vtilde^i = sqrt(g) * rho * W * vtilde^i */
    const vec<CCTK_REAL, 2> flux_dens(
        [&](int f) ARITH_INLINE { return dens_rc(f) * vtilde_rc(f); });

    /* flux(DEnt) = sqrt(g) * D * s * vtilde^i = sqrt(g) * rho * W * s *
     * vtilde^i */
    const vec<CCTK_REAL, 2> flux_DEnt(
        [&](int f) ARITH_INLINE { return DEnt_rc(f) * vtilde_rc(f); });

    /* flux(mom_j)^i = sqrt(g)*(
     *  S_j*vtilde^i + alpha*((pgas+pmag)*delta^i_j - b_jB^i/W) ) */
    const vec<vec<CCTK_REAL, 2>, 3> flux_moms([&](int j) ARITH_INLINE {
      return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
        return moms_rc(j)(f) * vtilde_rc(f) +
               alp_sqrtg * (press_plus_pmag_rc(f) * unit_dir(j) -
                            blows_rc(j)(f) * B_over_w_lorentz_rc(f));
      });
    });

    /* flux(tau) = sqrt(g)*(
     *  t*vtilde^i + alpha*((pgas+pmag)*v^i-alpha*b0*B^i/W) ) */
    const vec<CCTK_REAL, 2> flux_tau([&](int f) ARITH_INLINE {
      return tau_rc(f) * vtilde_rc(f) +
             alp_sqrtg * (press_plus_pmag_rc(f) * vel_rc(f) -
                          alp_b0_rc(f) * B_over_w_lorentz_rc(f));
    });

    /* flux(DYe) = sqrt(g) * (D * Ye * vtilde^i) */
    const vec<CCTK_REAL, 2> DYe_rc(
        [&](int f) ARITH_INLINE { return dens_rc(f) * rcprims.Ye_rc(f); });
    const vec<CCTK_REAL, 2> flux_DYe(
        [&](int f) ARITH_INLINE { return DYe_rc(f) * vtilde_rc(f); });

    /* electric field E_i = \tilde\epsilon_{ijk} Btilde_j * vtilde_k */
    const vec<vec<CCTK_REAL, 2>, 3> Es_rc =
        calc_cross_product(Btildes_rc, vtildes_rc);
    /* flux(Btildes) = {{0, -Ez, Ey}, {Ez, 0, -Ex}, {-Ey, Ex, 0}} */
    const vec<vec<CCTK_REAL, 2>, 3> flux_Btildes =
        calc_cross_product(unit_dir, Es_rc);

    /* Calculate eigenvalues: */

    /* variable for either g^xx, g^yy or g^zz depending on the direction */
    const CCTK_REAL u_avg = calc_inv(g_avg, sqrtg*sqrtg)(dir, dir);
    /* eigenvalues */
    vec<vec<CCTK_REAL, 4>, 2> lambda =
        eigenvalues(alp_avg, beta_avg, u_avg, vel_rc, rcprims.rho_rc, cs2_rc,
                    rcprims.w_lorentz_rc, h_rc, bsq_rc);

    // Save cons and fluxes
    recon_cons rccons_tmp{dens_rc, DEnt_rc, moms_rc, tau_rc, DYe_rc, Btildes_rc, vtildes_rc,
     lambda, flux_dens, flux_DEnt, flux_moms, flux_tau, flux_DYe, flux_Btildes};
    rccons = rccons_tmp; 
}


// Calculate the fluxes in direction `dir`. This function is more
// complex because it has to handle any direction, but as reward,
// there is only one function, not three.
template <int dir, typename EOSType>
void CalcFlux(CCTK_ARGUMENTS, EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  /* grid functions for fluxes */
  const vec<GF3D2<CCTK_REAL>, dim> fluxdenss{fxdens, fydens, fzdens};
  const vec<GF3D2<CCTK_REAL>, dim> fluxDEnts{fxDEnt, fyDEnt, fzDEnt};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomxs{fxmomx, fymomx, fzmomx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomys{fxmomy, fymomy, fzmomy};
  const vec<GF3D2<CCTK_REAL>, dim> fluxmomzs{fxmomz, fymomz, fzmomz};
  const vec<GF3D2<CCTK_REAL>, dim> fluxtaus{fxtau, fytau, fztau};
  const vec<GF3D2<CCTK_REAL>, dim> fluxDYes{fxDYe, fyDYe, fzDYe};
  const vec<GF3D2<CCTK_REAL>, dim> fluxBxs{fxBx, fyBx, fzBx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxBys{fxBy, fyBy, fzBy};
  const vec<GF3D2<CCTK_REAL>, dim> fluxBzs{fxBz, fyBz, fzBz};
  /* grid functions for PP flux limiter */
  const vec<GF3D2<CCTK_REAL>, dim> thetagf{theta_x, theta_y, theta_z};
  const vec<GF3D2<CCTK_REAL>, dim> fluxLOdenss{fLOxdens, fLOydens, fLOzdens};
  const vec<GF3D2<CCTK_REAL>, dim> fluxLODEnts{fLOxDEnt, fLOyDEnt, fLOzDEnt};
  const vec<GF3D2<CCTK_REAL>, dim> fluxLODYes{fLOxDYe, fLOyDYe, fLOzDYe};
  const vec<GF3D2<CCTK_REAL>, dim> fluxLOmomxs{fLOxmomx, fLOymomx, fLOzmomx};
  const vec<GF3D2<CCTK_REAL>, dim> fluxLOmomys{fLOxmomy, fLOymomy, fLOzmomy};
  const vec<GF3D2<CCTK_REAL>, dim> fluxLOmomzs{fLOxmomz, fLOymomz, fLOzmomz};
  const vec<GF3D2<CCTK_REAL>, dim> fluxLOtaus{fLOxtau, fLOytau, fLOztau};
  /* grid functions */
  const vec<GF3D2<const CCTK_REAL>, dim> gf_vels{velx, vely, velz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_zvec{zvec_x, zvec_y, zvec_z};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_svec{svec_x, svec_y, svec_z};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_Bvecs{Bvecx, Bvecy, Bvecz};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_dBstags{dBx_stag, dBy_stag,
                                                    dBz_stag};
  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  /* grid functions for Upwind CT */
  const vec<GF3D2<CCTK_REAL>, dim> vtildes_one{vtilde_y_xface, vtilde_z_yface,
                                               vtilde_x_zface};
  const vec<GF3D2<CCTK_REAL>, dim> vtildes_two{vtilde_z_xface, vtilde_x_yface,
                                               vtilde_y_zface};
  const vec<GF3D2<CCTK_REAL>, dim> amax{amax_xface, amax_yface, amax_zface};
  const vec<GF3D2<CCTK_REAL>, dim> amin{amin_xface, amin_yface, amin_zface};

  static_assert(dir >= 0 && dir < 3, "");

  rec_var_t rec_var;
  if (CCTK_EQUALS(recon_type, "v_vec")) {
    rec_var = rec_var_t::v_vec;
  } else if (CCTK_EQUALS(recon_type, "z_vec")) {
    rec_var = rec_var_t::z_vec;
  } else if (CCTK_EQUALS(recon_type, "s_vec")) {
    rec_var = rec_var_t::s_vec;
  } else {
    CCTK_ERROR("Unknown value for parameter \"recon_type\"");
  }

  // Primary reconstruction method
  reconstruction_t reconstruction;
  if (CCTK_EQUALS(reconstruction_method, "Godunov"))
    reconstruction = reconstruction_t::Godunov;
  else if (CCTK_EQUALS(reconstruction_method, "minmod"))
    reconstruction = reconstruction_t::minmod;
  else if (CCTK_EQUALS(reconstruction_method, "monocentral"))
    reconstruction = reconstruction_t::monocentral;
  else if (CCTK_EQUALS(reconstruction_method, "ppm"))
    reconstruction = reconstruction_t::ppm;
  else if (CCTK_EQUALS(reconstruction_method, "eppm"))
    reconstruction = reconstruction_t::eppm;
  else if (CCTK_EQUALS(reconstruction_method, "wenoz"))
    reconstruction = reconstruction_t::wenoz;
  else if (CCTK_EQUALS(reconstruction_method, "mp5"))
    reconstruction = reconstruction_t::mp5;
  else
    CCTK_ERROR("Unknown value for parameter \"reconstruction_method\"");

  bool rc_temp = reconstruct_with_temperature;

  // Lower-order fallback for negative values
  reconstruction_t reconstruction_LO;
  if (CCTK_EQUALS(loworder_method, "Godunov"))
    reconstruction_LO = reconstruction_t::Godunov;
  else if (CCTK_EQUALS(loworder_method, "minmod"))
    reconstruction_LO = reconstruction_t::minmod;
  else if (CCTK_EQUALS(loworder_method, "monocentral"))
    reconstruction_LO = reconstruction_t::monocentral;
  else
    CCTK_ERROR("Unknown value for parameter \"loworder_method\"");

  flux_t fluxtype;
  if (CCTK_EQUALS(flux_type, "LxF")) {
    fluxtype = flux_t::LxF;
  } else if (CCTK_EQUALS(flux_type, "HLLE")) {
    fluxtype = flux_t::HLLE;
  } else {
    CCTK_ERROR("Unknown value for parameter \"flux_type\"");
  }

  switch (reconstruction) {
  case reconstruction_t::Godunov:
    assert(cctk_nghostzones[dir] >= 1);
    break;
  case reconstruction_t::minmod:
    assert(cctk_nghostzones[dir] >= 2);
    break;
  case reconstruction_t::monocentral:
    assert(cctk_nghostzones[dir] >= 2);
    break;
  case reconstruction_t::ppm:
    assert(cctk_nghostzones[dir] >= 3);
    break;
  case reconstruction_t::eppm:
    assert(cctk_nghostzones[dir] >= 3);
    break;
  case reconstruction_t::wenoz:
    assert(cctk_nghostzones[dir] >= 3);
  case reconstruction_t::mp5:
    assert(cctk_nghostzones[dir] >= 3);
    break;
  }

  // reconstruction parameters struct
  reconstruct_params_t reconstruct_params;

  // ppm parameters
  reconstruct_params.ppm_shock_detection = ppm_shock_detection;
  reconstruct_params.ppm_zone_flattening = ppm_zone_flattening;
  reconstruct_params.poly_k = poly_k;
  reconstruct_params.poly_gamma = poly_gamma;
  reconstruct_params.ppm_eta1 = ppm_eta1;
  reconstruct_params.ppm_eta2 = ppm_eta2;
  reconstruct_params.ppm_eps = ppm_eps;
  reconstruct_params.ppm_eps_shock = ppm_eps_shock;
  reconstruct_params.ppm_small = ppm_small;
  reconstruct_params.ppm_omega1 = ppm_omega1;
  reconstruct_params.ppm_omega2 = ppm_omega2;
  reconstruct_params.enhanced_ppm_C2 = enhanced_ppm_C2;
  // wenoz parameters
  reconstruct_params.weno_eps = weno_eps;
  // mp5 parameters
  reconstruct_params.mp5_alpha = mp5_alpha;

  const auto reconstruct_pt =
      [=] CCTK_DEVICE(const GF3D2<const CCTK_REAL> &var, const PointDesc &p,
                      const bool &gf_is_rho,
                      const bool &gf_is_press) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        return reconstruct(var, p, reconstruction, dir, gf_is_rho, gf_is_press,
                           press, gf_vels(dir), reconstruct_params);
      };
  const auto reconstruct_loworder =
      [=] CCTK_DEVICE(const GF3D2<const CCTK_REAL> &var, const PointDesc &p,
                      const bool &gf_is_rho,
                      const bool &gf_is_press) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        return reconstruct(var, p, reconstruction_LO, dir, gf_is_rho, gf_is_press,
                           press, gf_vels(dir), reconstruct_params);
      };
  const auto calcflux =
      [=] CCTK_DEVICE(vec<vec<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
                      vec<CCTK_REAL, 2> flux) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL flx;
        switch (fluxtype) {

        case flux_t::LxF: {
          flx = laxf(lam, var, flux);
          break;
        }

        case flux_t::HLLE: {
          flx = hlle(lam, var, flux);
          break;
        }

        default:
          assert(0);
        }

        return flx;
      };

  // Face-centred grid functions (in direction `dir`)
  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};

  // Precompute mapping table once for clarity and efficiency
  constexpr array<array<int, 3>, 3> dir_arr_table = {{
      {0, 1, 2}, // dir == 0: x, y, z
      {1, 2, 0}, // dir == 1: y, z, x
      {2, 0, 1}  // dir == 2: z, x, y
  }};
  constexpr auto dir_arr = dir_arr_table[dir];

  grid.loop_int_device<
      face_centred[0], face_centred[1],
      face_centred
          [2]>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    /* Reconstruct primitives from the cells on left (indice 0) and right
     * (indice 1) side of this face rc = reconstructed variables or
     * computed from reconstructed variables */

    /* Interpolate metric components from vertices to faces */
    const CCTK_REAL alp_avg = calc_avg_v2f(alp, p, dir);
    const vec<CCTK_REAL, 3> betas_avg(
        [&](int i) ARITH_INLINE { return calc_avg_v2f(gf_beta(i), p, dir); });
    const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
      return calc_avg_v2f(gf_g(i, j), p, dir);
    });

    /* determinant of spatial metric */
    const CCTK_REAL detg_avg = calc_det(g_avg);
    const CCTK_REAL sqrtg = sqrt(detg_avg);

    // ---------------------------- BEGIN LO FOR FLUX LIM -----------------------------------
    // Reconstruct density
    vec<CCTK_REAL, 2> rhoLO_rc{reconstruct_loworder(rho, p, true, true)};

    // Reconstruct entropy
    vec<CCTK_REAL, 2> entropyLO_rc{reconstruct_loworder(entropy, p, false, false)};

    // Reconstruct Ye
    vec<CCTK_REAL, 2> YeLO_rc{reconstruct_loworder(Ye, p, false, false)};

    // Initialize variables for eps, pressure, and temperature
    vec<CCTK_REAL, 2> epsLO_rc;
    vec<CCTK_REAL, 2> pressLO_rc;
    vec<CCTK_REAL, 2> tempLO_rc;

    if (reconstruct_with_temperature) {

      // Reconstruct temperature
      vec<CCTK_REAL, 2> tempLO_rc{reconstruct_loworder(temperature, p, false, false)};

      // Compute eps_rc and press_rc using lambdas
      for (int f = 0; f < 2; ++f) {
        epsLO_rc(f) =
            eos_3p->eps_from_valid_rho_temp_ye(rhoLO_rc(f), tempLO_rc(f), YeLO_rc(f));
        pressLO_rc(f) = eos_3p->press_from_valid_rho_temp_ye(
            rhoLO_rc(f), tempLO_rc(f), YeLO_rc(f));
      }

    } else {
      printf("Flux limiters only supported for temperature reconstruction for now.");
      assert(0);
    }

    const vec<CCTK_REAL, 2> rhohLO_rc([&](int f) ARITH_INLINE {
      return rhoLO_rc(f) + rhoLO_rc(f) * epsLO_rc(f) + pressLO_rc(f);
    });

    // Introduce reconstructed Bs
    // Use staggered dB for i == dir
    vec<vec<CCTK_REAL, 2>, 3> BsLO_rc;

    // Assign the value for the primary direction
    const CCTK_REAL val = gf_dBstags(dir)(p.I) / sqrtg;
    BsLO_rc(dir)(0) = val;
    BsLO_rc(dir)(1) = val;

    // Lambda to assign the reconstructed values
    auto assign_reconstructedLO = [&](int d) {
      auto tmp = reconstruct_loworder(gf_Bvecs(d), p, false, false);
      BsLO_rc(d)(0) = tmp[0];
      BsLO_rc(d)(1) = tmp[1];
    };

    // Assign reconstructed values for the two perpendicular directions
    assign_reconstructedLO(dir_arr[1]);
    assign_reconstructedLO(dir_arr[2]);
    // End of setting Bs
    
    vec<vec<CCTK_REAL, 2>, 3> velsLO_rc;
    vec<vec<CCTK_REAL, 2>, 3> vlowsLO_rc;
    vec<CCTK_REAL, 2> w_lorentzLO_rc;
    switch (rec_var) {
    case rec_var_t::v_vec: {

      array<CCTK_REAL, 2>
          vels_rc_dummy; // note: can't copy array<,2> to vec<,2>, only construct

      for (int i = 0; i <= 2; ++i) { // loop over components
        vels_rc_dummy = reconstruct_loworder(gf_vels(i), p, false, false);
        velsLO_rc(i)(0) = vels_rc_dummy[0];
        velsLO_rc(i)(1) = vels_rc_dummy[1];
      }

      /* co-velocity measured by Eulerian observer: v_j */
      vlowsLO_rc = calc_contraction(g_avg, velsLO_rc);

      /* Lorentz factor: W = 1 / sqrt(1 - v^2) */
      w_lorentzLO_rc(0) = 1 / sqrt(1 - calc_contraction(vlowsLO_rc, velsLO_rc)(0));
      w_lorentzLO_rc(1) = 1 / sqrt(1 - calc_contraction(vlowsLO_rc, velsLO_rc)(1));
      break;
    };
    case rec_var_t::z_vec: {

      vec<vec<CCTK_REAL, 2>, 3> zvec_rc([&](int i) ARITH_INLINE {
        return vec<CCTK_REAL, 2>{reconstruct_loworder(gf_zvec(i), p, false, false)};
      });

      const vec<vec<CCTK_REAL, 2>, 3> zveclow_rc =
          calc_contraction(g_avg, zvec_rc);

      w_lorentzLO_rc(0) = sqrt(1 + calc_contraction(zveclow_rc, zvec_rc)(0));
      w_lorentzLO_rc(1) = sqrt(1 + calc_contraction(zveclow_rc, zvec_rc)(1));

      for (int i = 0; i <= 2; ++i) {   // loop over components
        for (int j = 0; j <= 1; ++j) { // loop over left and right state
          velsLO_rc(i)(j) = zvec_rc(i)(j) / w_lorentzLO_rc(j);
          vlowsLO_rc(i)(j) = zveclow_rc(i)(j) / w_lorentzLO_rc(j);
        }
      }
      break;
    };
    case rec_var_t::s_vec: {

      vec<vec<CCTK_REAL, 2>, 3> svec_rc([&](int i) ARITH_INLINE {
        return vec<CCTK_REAL, 2>{reconstruct_loworder(gf_svec(i), p, false, false)};
      });

      const vec<vec<CCTK_REAL, 2>, 3> sveclow_rc =
          calc_contraction(g_avg, svec_rc);

      w_lorentzLO_rc(0) =
          sqrt(0.5 + sqrt(0.25 + calc_contraction(sveclow_rc, svec_rc)(0) /
                                     rhohLO_rc(0) / rhohLO_rc(0)));
      w_lorentzLO_rc(1) =
          sqrt(0.5 + sqrt(0.25 + calc_contraction(sveclow_rc, svec_rc)(1) /
                                     rhohLO_rc(1) / rhohLO_rc(1)));

      for (int i = 0; i <= 2; ++i) {   // loop over components
        for (int j = 0; j <= 1; ++j) { // loop over left and right state
          velsLO_rc(i)(j) =
              svec_rc(i)(j) / w_lorentzLO_rc(j) / w_lorentzLO_rc(j) / rhohLO_rc(j);
          vlowsLO_rc(i)(j) =
              sveclow_rc(i)(j) / w_lorentzLO_rc(j) / w_lorentzLO_rc(j) / rhohLO_rc(j);
        }
      }
      break;
    };
    }

    recon_prims rcprimsLO{rhoLO_rc, entropyLO_rc, YeLO_rc, epsLO_rc, pressLO_rc, tempLO_rc, BsLO_rc,
      velsLO_rc, vlowsLO_rc, w_lorentzLO_rc};
    recon_cons rcconsLO;

    CalcConsFluxesFromPrims(eos_3p, rc_temp,
     alp_avg, betas_avg, g_avg, sqrtg,
     rcprimsLO, rcconsLO);
    // ---------------------------- END LO FOR FLUX LIM -----------------------------------
    
    
    
    // Booleans to decide whether to use low order
    // reconstruction
    bool useLO_0 = false;
    bool useLO_1 = false;

    // Reconstruct density
    vec<CCTK_REAL, 2> rho_rc{reconstruct_pt(rho, p, true, true)};

    // Reconstruct entropy
    vec<CCTK_REAL, 2> entropy_rc{reconstruct_pt(entropy, p, false, false)};

    // Reconstruct Ye
    vec<CCTK_REAL, 2> Ye_rc{reconstruct_pt(Ye, p, false, false)};

    // Initialize variables for eps, pressure, and temperature
    vec<CCTK_REAL, 2> eps_rc;
    vec<CCTK_REAL, 2> press_rc;
    vec<CCTK_REAL, 2> temp_rc;

    if (reconstruct_with_temperature) {

      // Reconstruct temperature
      array<CCTK_REAL, 2> temp_rc_dummy;
      temp_rc_dummy = reconstruct_pt(temperature, p, false, false);

      // Use lower-order if reconstructed rho, entropy, Ye or T is <= 0
      if ((rho_rc(0) <= 0.0) || (entropy_rc(0) <= 0.0) || (Ye_rc(0) <= 0.0) || (temp_rc_dummy[0] <= 0.0)) {
        useLO_0 = true;
      }
      if ((rho_rc(1) <= 0.0) || (entropy_rc(1) <= 0.0) || (Ye_rc(1) <= 0.0) || (temp_rc_dummy[1] <= 0.0)) {
        useLO_1 = true;
      }

      // Lower-order
      if (useLO_0) {
        vec<CCTK_REAL, 2> rhoLO_rc{reconstruct_loworder(rho, p, true, true)};
        vec<CCTK_REAL, 2> entropyLO_rc{reconstruct_loworder(entropy, p, false, false)};
        vec<CCTK_REAL, 2> YeLO_rc{reconstruct_loworder(Ye, p, false, false)};
        vec<CCTK_REAL, 2> tempLO_rc{reconstruct_loworder(temperature, p, false, false)};

        rho_rc(0) = rhoLO_rc(0);
        entropy_rc(0) = entropyLO_rc(0);
        Ye_rc(0) = YeLO_rc(0);
        temp_rc_dummy[0] = tempLO_rc(0);
      }

      if (useLO_1) {
        vec<CCTK_REAL, 2> rhoLO_rc{reconstruct_loworder(rho, p, true, true)};
        vec<CCTK_REAL, 2> entropyLO_rc{reconstruct_loworder(entropy, p, false, false)};
        vec<CCTK_REAL, 2> YeLO_rc{reconstruct_loworder(Ye, p, false, false)};
        vec<CCTK_REAL, 2> tempLO_rc{reconstruct_loworder(temperature, p, false, false)};

        rho_rc(1) = rhoLO_rc(1);
        entropy_rc(1) = entropyLO_rc(1);
        Ye_rc(1) = YeLO_rc(1);
        temp_rc_dummy[1] = tempLO_rc(1);
      }
      // End lower-order

      // Compute eps_rc and press_rc using lambdas
      for (int f = 0; f < 2; ++f) {
        temp_rc(f) = temp_rc_dummy[f];
        eps_rc(f) =
            eos_3p->eps_from_valid_rho_temp_ye(rho_rc(f), temp_rc(f), Ye_rc(f));
        press_rc(f) = eos_3p->press_from_valid_rho_temp_ye(
            rho_rc(f), temp_rc(f), Ye_rc(f));
      }

    } else {

      // Reconstruct pressure
      array<CCTK_REAL, 2> press_rc_dummy;
      press_rc_dummy = reconstruct_pt(press, p, false, true);

      // Use lower-order if reconstructed rho, entropy, Ye or pressure is <= 0
      if ((rho_rc(0) <= 0.0) || (entropy_rc(0) <= 0.0) || (Ye_rc(0) <= 0.0) || (press_rc_dummy[0] <= 0.0)) {
        useLO_0 = true;
      }
      if ((rho_rc(1) <= 0.0) || (entropy_rc(1) <= 0.0) || (Ye_rc(1) <= 0.0) || (press_rc_dummy[1] <= 0.0)) {
        useLO_1 = true;
      }

      // Lower-order
      if (useLO_0) {
        vec<CCTK_REAL, 2> rhoLO_rc{reconstruct_loworder(rho, p, true, true)};
        vec<CCTK_REAL, 2> entropyLO_rc{reconstruct_loworder(entropy, p, false, false)};
        vec<CCTK_REAL, 2> YeLO_rc{reconstruct_loworder(Ye, p, false, false)};
        vec<CCTK_REAL, 2> pressLO_rc{reconstruct_loworder(press, p, false, true)};

        rho_rc(0) = rhoLO_rc(0);
        entropy_rc(0) = entropyLO_rc(0);
        Ye_rc(0) = YeLO_rc(0);
        press_rc_dummy[0] = pressLO_rc(0);
      }

      if (useLO_1) {
        vec<CCTK_REAL, 2> rhoLO_rc{reconstruct_loworder(rho, p, true, true)};
        vec<CCTK_REAL, 2> entropyLO_rc{reconstruct_loworder(entropy, p, false, false)};
        vec<CCTK_REAL, 2> YeLO_rc{reconstruct_loworder(Ye, p, false, false)};
        vec<CCTK_REAL, 2> pressLO_rc{reconstruct_loworder(press, p, false, true)};

        rho_rc(1) = rhoLO_rc(1);
        entropy_rc(1) = entropyLO_rc(1);
        Ye_rc(1) = YeLO_rc(1);
        press_rc_dummy[1] = pressLO_rc(1);
      }
      // End lower-order

      // Compute eps_rc and temp_rc using lambdas
      for (int f = 0; f < 2; ++f) {
        press_rc(f) = press_rc_dummy[f];
        eps_rc(f) = eos_3p->eps_from_valid_rho_press_ye(rho_rc(f), press_rc(f),
                                                        Ye_rc(f));
        temp_rc(f) =
            eos_3p->temp_from_valid_rho_eps_ye(rho_rc(f), eps_rc(f), Ye_rc(f));
      }

    }

    const vec<CCTK_REAL, 2> rhoh_rc([&](int f) ARITH_INLINE {
      return rho_rc(f) + rho_rc(f) * eps_rc(f) + press_rc(f);
    });

    // Introduce reconstructed Bs
    // Use staggered dB for i == dir
    vec<vec<CCTK_REAL, 2>, 3> Bs_rc;

    // Assign the value for the primary direction
    // const CCTK_REAL val = gf_dBstags(dir)(p.I) / sqrtg;
    Bs_rc(dir)(0) = val;
    Bs_rc(dir)(1) = val;

    // Lambda to assign the reconstructed values
    auto assign_reconstructed = [&](int d) {
      auto tmp = reconstruct_pt(gf_Bvecs(d), p, false, false);
      Bs_rc(d)(0) = tmp[0];
      Bs_rc(d)(1) = tmp[1];

      // Lower-order
      if (useLO_0) {
        tmp = reconstruct_loworder(gf_Bvecs(d), p, false, false);
        Bs_rc(d)(0) = tmp[0];
      }
      if (useLO_1) {
        tmp = reconstruct_loworder(gf_Bvecs(d), p, false, false);
        Bs_rc(d)(1) = tmp[1];
      }
      // End lower-order
    };

    // Assign reconstructed values for the two perpendicular directions
    assign_reconstructed(dir_arr[1]);
    assign_reconstructed(dir_arr[2]);
    // End of setting Bs

    vec<vec<CCTK_REAL, 2>, 3> vels_rc;
    vec<vec<CCTK_REAL, 2>, 3> vlows_rc;
    vec<CCTK_REAL, 2> w_lorentz_rc;
    switch (rec_var) {
    case rec_var_t::v_vec: {

      array<CCTK_REAL, 2>
          vels_rc_dummy; // note: can't copy array<,2> to vec<,2>, only construct

      for (int i = 0; i <= 2; ++i) { // loop over components
        vels_rc_dummy = reconstruct_pt(gf_vels(i), p, false, false);
        vels_rc(i)(0) = vels_rc_dummy[0];
        vels_rc(i)(1) = vels_rc_dummy[1];

	// Lower-order
        if (useLO_0) {
          vels_rc_dummy = reconstruct_loworder(gf_vels(i), p, false, false);
          vels_rc(i)(0) = vels_rc_dummy[0];
        }
        if (useLO_1) {
          vels_rc_dummy = reconstruct_loworder(gf_vels(i), p, false, false);
          vels_rc(i)(1) = vels_rc_dummy[1];
        }
	// End lower-order
      }

      /* co-velocity measured by Eulerian observer: v_j */
      vlows_rc = calc_contraction(g_avg, vels_rc);

      /* Lorentz factor: W = 1 / sqrt(1 - v^2) */
      w_lorentz_rc(0) = 1 / sqrt(1 - calc_contraction(vlows_rc, vels_rc)(0));
      w_lorentz_rc(1) = 1 / sqrt(1 - calc_contraction(vlows_rc, vels_rc)(1));
      break;
    };
    case rec_var_t::z_vec: {

      vec<vec<CCTK_REAL, 2>, 3> zvec_rc([&](int i) ARITH_INLINE {
        return vec<CCTK_REAL, 2>{reconstruct_pt(gf_zvec(i), p, false, false)};
      });

      // Lower-order
      if (useLO_0) {
        vec<vec<CCTK_REAL, 2>, 3> zvecLO_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>{reconstruct_loworder(gf_zvec(i), p, false, false)};
        });
        
	for (int i = 0; i <= 2; ++i) { // loop over components
          zvec_rc(i)(0) = zvecLO_rc(i)(0);
	}
      }
      if (useLO_1) {
        vec<vec<CCTK_REAL, 2>, 3> zvecLO_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>{reconstruct_loworder(gf_zvec(i), p, false, false)};
        });
        
	for (int i = 0; i <= 2; ++i) { // loop over components
          zvec_rc(i)(1) = zvecLO_rc(i)(1);
	}
      }
      // End lower-order

      const vec<vec<CCTK_REAL, 2>, 3> zveclow_rc =
          calc_contraction(g_avg, zvec_rc);

      w_lorentz_rc(0) = sqrt(1 + calc_contraction(zveclow_rc, zvec_rc)(0));
      w_lorentz_rc(1) = sqrt(1 + calc_contraction(zveclow_rc, zvec_rc)(1));

      for (int i = 0; i <= 2; ++i) {   // loop over components
        for (int j = 0; j <= 1; ++j) { // loop over left and right state
          vels_rc(i)(j) = zvec_rc(i)(j) / w_lorentz_rc(j);
          vlows_rc(i)(j) = zveclow_rc(i)(j) / w_lorentz_rc(j);
        }
      }
      break;
    };
    case rec_var_t::s_vec: {

      vec<vec<CCTK_REAL, 2>, 3> svec_rc([&](int i) ARITH_INLINE {
        return vec<CCTK_REAL, 2>{reconstruct_pt(gf_svec(i), p, false, false)};
      });

      // Lower-order
      if (useLO_0) {
        vec<vec<CCTK_REAL, 2>, 3> svecLO_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>{reconstruct_loworder(gf_svec(i), p, false, false)};
        });
        
	for (int i = 0; i <= 2; ++i) { // loop over components
          svec_rc(i)(0) = svecLO_rc(i)(0);
	}
      }
      if (useLO_1) {
        vec<vec<CCTK_REAL, 2>, 3> svecLO_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>{reconstruct_loworder(gf_svec(i), p, false, false)};
        });
        
	for (int i = 0; i <= 2; ++i) { // loop over components
          svec_rc(i)(1) = svecLO_rc(i)(1);
	}
      }
      // End lower-order

      const vec<vec<CCTK_REAL, 2>, 3> sveclow_rc =
          calc_contraction(g_avg, svec_rc);

      w_lorentz_rc(0) =
          sqrt(0.5 + sqrt(0.25 + calc_contraction(sveclow_rc, svec_rc)(0) /
                                     rhoh_rc(0) / rhoh_rc(0)));
      w_lorentz_rc(1) =
          sqrt(0.5 + sqrt(0.25 + calc_contraction(sveclow_rc, svec_rc)(1) /
                                     rhoh_rc(1) / rhoh_rc(1)));

      // printf("  wlor = %16.8e, %16.8e\n", w_lorentz_rc(0), w_lorentz_rc(1));

      for (int i = 0; i <= 2; ++i) {   // loop over components
        for (int j = 0; j <= 1; ++j) { // loop over left and right state
          vels_rc(i)(j) =
              svec_rc(i)(j) / w_lorentz_rc(j) / w_lorentz_rc(j) / rhoh_rc(j);
          vlows_rc(i)(j) =
              sveclow_rc(i)(j) / w_lorentz_rc(j) / w_lorentz_rc(j) / rhoh_rc(j);
        }
      }
      break;
    };
    }

    // /* vtilde^i = alpha * v^i - beta^i */
    // const vec<vec<CCTK_REAL, 2>, 3> vtildes_rc([&](int i) ARITH_INLINE {
    //   return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
    //     return alp_avg * vels_rc(i)(f) - betas_avg(i);
    //   });
    // });

    // /* alpha * b0 = W * B^i * v_i */
    // const vec<CCTK_REAL, 2> alp_b0_rc([&](int f) ARITH_INLINE {
    //   return w_lorentz_rc(f) * calc_contraction(Bs_rc, vlows_rc)(f);
    // });
    // /* covariant magnetic field measured by the Eulerian observer */
    // const vec<vec<CCTK_REAL, 2>, 3> Blows_rc = calc_contraction(g_avg, Bs_rc);
    // /* B^2 = B^i * B_i */
    // const vec<CCTK_REAL, 2> B2_rc = calc_contraction(Bs_rc, Blows_rc);
    // /* covariant magnetic field measured by the comoving observer:
    //  *  b_i = B_i/W + alpha*b^0*v_i */
    // const vec<vec<CCTK_REAL, 2>, 3> blows_rc([&](int i) ARITH_INLINE {
    //   return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
    //     return Blows_rc(i)(f) / w_lorentz_rc(f) + alp_b0_rc(f) * vlows_rc(i)(f);
    //   });
    // });
    // /* b^2 = b^{\mu} * b_{\mu} */
    // const vec<CCTK_REAL, 2> bsq_rc([&](int f) ARITH_INLINE {
    //   return (B2_rc(f) + pow2(alp_b0_rc(f))) / pow2(w_lorentz_rc(f));
    // });

    // /* componets correspond to the dir we are considering */
    // const CCTK_REAL beta_avg = betas_avg(dir);
    // const vec<CCTK_REAL, 2> vel_rc{vels_rc(dir)};
    // const vec<CCTK_REAL, 2> B_rc{Bs_rc(dir)};
    // const vec<CCTK_REAL, 2> vtilde_rc{vtildes_rc(dir)};

    // // TODO: Compute pressure based on user-specified EOS.
    // // Currently, computing press for classical ideal gas from reconstructed
    // // vars

    // const vec<CCTK_REAL, 2> cs2_rc([&](int f) ARITH_INLINE {
    //   return eos_3p->csnd_from_valid_rho_eps_ye(rho_rc(f), eps_rc(f),
    //                                             Ye_rc(f)) *
    //          eos_3p->csnd_from_valid_rho_eps_ye(rho_rc(f), eps_rc(f), Ye_rc(f));
    // });

    // const vec<CCTK_REAL, 2> h_rc([&](int f) ARITH_INLINE {
    //   return 1 + eps_rc(f) + press_rc(f) / rho_rc(f);
    // });

    // /* Computing conservatives from primitives: */

    // /* dens = sqrt(g) * D = sqrt(g) * (rho * W) */
    // const vec<CCTK_REAL, 2> dens_rc([&](int f) ARITH_INLINE {
    //   return sqrtg * rho_rc(f) * w_lorentz_rc(f);
    // });

    // /* DEnt = sqrt(g) * D * s  = sqrt(g) * (rho * W) * s */
    // /*    s = entropy */
    // const vec<CCTK_REAL, 2> DEnt_rc([&](int f) ARITH_INLINE {
    //   return sqrtg * rho_rc(f) * w_lorentz_rc(f) * entropy_rc(f);
    // });

    // /* auxiliary: dens * h * W = sqrt(g) * rho * h * W^2 */
    // const vec<CCTK_REAL, 2> dens_h_W_rc([&](int f) ARITH_INLINE {
    //   return dens_rc(f) * h_rc(f) * w_lorentz_rc(f);
    // });
    // /* auxiliary: sqrt(g) * (rho*h + b^2)*W^2 */
    // const vec<CCTK_REAL, 2> dens_h_W_plus_sqrtg_W2b2_rc =
    //     dens_h_W_rc + sqrtg * (pow2(alp_b0_rc) + B2_rc);
    // /* auxiliary: (pgas + pmag) */
    // const vec<CCTK_REAL, 2> press_plus_pmag_rc = press_rc + 0.5 * bsq_rc;

    // /* mom_i = sqrt(g)*S_i = sqrt(g)((rho*h+b^2)*W^2*v_i - alpha*b^0*b_i) */
    // const vec<vec<CCTK_REAL, 2>, 3> moms_rc([&](int i) ARITH_INLINE {
    //   return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
    //     return dens_h_W_plus_sqrtg_W2b2_rc(f) * vlows_rc(i)(f) -
    //            sqrtg * alp_b0_rc(f) * blows_rc(i)(f);
    //   });
    // });

    // /* tau = sqrt(g)*t =
    //  *  sqrt(g)((rho*h + b^2)*W^2 - (pgas+pmag) - (alpha*b^0)^2 - D) */
    // const vec<CCTK_REAL, 2> tau_rc =
    //     dens_h_W_rc - dens_rc + sqrtg * (B2_rc - press_plus_pmag_rc);

    // /* Btildes^i = sqrt(g) * B^i */
    // const vec<vec<CCTK_REAL, 2>, 3> Btildes_rc(
    //     [&](int i) ARITH_INLINE { return sqrtg * Bs_rc(i); });

    // /* Computing fluxes of conserved variables: */

    // /* auxiliary: unit in 'dir' */
    // const vec<CCTK_REAL, 3> unit_dir{vec<int, 3>::unit(dir)};
    // /* auxiliary: alpha * sqrt(g) */
    // const CCTK_REAL alp_sqrtg = alp_avg * sqrtg;
    // /* auxiliary: B^i / W */
    // const vec<CCTK_REAL, 2> B_over_w_lorentz_rc(
    //     [&](int f) ARITH_INLINE { return B_rc(f) / w_lorentz_rc(f); });

    // /* flux(dens) = sqrt(g) * D * vtilde^i = sqrt(g) * rho * W * vtilde^i */
    // const vec<CCTK_REAL, 2> flux_dens(
    //     [&](int f) ARITH_INLINE { return dens_rc(f) * vtilde_rc(f); });

    // /* flux(DEnt) = sqrt(g) * D * s * vtilde^i = sqrt(g) * rho * W * s *
    //  * vtilde^i */
    // const vec<CCTK_REAL, 2> flux_DEnt(
    //     [&](int f) ARITH_INLINE { return DEnt_rc(f) * vtilde_rc(f); });

    // /* flux(mom_j)^i = sqrt(g)*(
    //  *  S_j*vtilde^i + alpha*((pgas+pmag)*delta^i_j - b_jB^i/W) ) */
    // const vec<vec<CCTK_REAL, 2>, 3> flux_moms([&](int j) ARITH_INLINE {
    //   return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
    //     return moms_rc(j)(f) * vtilde_rc(f) +
    //            alp_sqrtg * (press_plus_pmag_rc(f) * unit_dir(j) -
    //                         blows_rc(j)(f) * B_over_w_lorentz_rc(f));
    //   });
    // });

    // /* flux(tau) = sqrt(g)*(
    //  *  t*vtilde^i + alpha*((pgas+pmag)*v^i-alpha*b0*B^i/W) ) */
    // const vec<CCTK_REAL, 2> flux_tau([&](int f) ARITH_INLINE {
    //   return tau_rc(f) * vtilde_rc(f) +
    //          alp_sqrtg * (press_plus_pmag_rc(f) * vel_rc(f) -
    //                       alp_b0_rc(f) * B_over_w_lorentz_rc(f));
    // });

    // /* flux(DYe) = sqrt(g) * (D * Ye * vtilde^i) */
    // const vec<CCTK_REAL, 2> DYe_rc(
    //     [&](int f) ARITH_INLINE { return dens_rc(f) * Ye_rc(f); });
    // const vec<CCTK_REAL, 2> flux_DYe(
    //     [&](int f) ARITH_INLINE { return DYe_rc(f) * vtilde_rc(f); });

    // /* electric field E_i = \tilde\epsilon_{ijk} Btilde_j * vtilde_k */
    // const vec<vec<CCTK_REAL, 2>, 3> Es_rc =
    //     calc_cross_product(Btildes_rc, vtildes_rc);
    // /* flux(Btildes) = {{0, -Ez, Ey}, {Ez, 0, -Ex}, {-Ey, Ex, 0}} */
    // const vec<vec<CCTK_REAL, 2>, 3> flux_Btildes =
    //     calc_cross_product(unit_dir, Es_rc);

    // /* Calculate eigenvalues: */

    // /* variable for either g^xx, g^yy or g^zz depending on the direction */
    // const CCTK_REAL u_avg = calc_inv(g_avg, detg_avg)(dir, dir);
    // /* eigenvalues */
    // vec<vec<CCTK_REAL, 4>, 2> lambda =
    //     eigenvalues(alp_avg, beta_avg, u_avg, vel_rc, rho_rc, cs2_rc,
    //                 w_lorentz_rc, h_rc, bsq_rc);

    recon_prims rcprimsHO{rho_rc, entropy_rc, Ye_rc, eps_rc, press_rc, temp_rc, Bs_rc,
      vels_rc, vlows_rc, w_lorentz_rc};
    recon_cons rcconsHO;

    CalcConsFluxesFromPrims(eos_3p, rc_temp,
     alp_avg, betas_avg, g_avg, sqrtg,
     rcprimsHO, rcconsHO);

    /* Calculate numerical fluxes */
    fluxdenss(dir)(p.I) = calcflux(rcconsHO.lambda, rcconsHO.dens_rc,       rcconsHO.flux_dens);
    fluxDEnts(dir)(p.I) = calcflux(rcconsHO.lambda, rcconsHO.DEnt_rc,       rcconsHO.flux_DEnt);
    fluxmomxs(dir)(p.I) = calcflux(rcconsHO.lambda, rcconsHO.moms_rc(0),    rcconsHO.flux_moms(0));
    fluxmomys(dir)(p.I) = calcflux(rcconsHO.lambda, rcconsHO.moms_rc(1),    rcconsHO.flux_moms(1));
    fluxmomzs(dir)(p.I) = calcflux(rcconsHO.lambda, rcconsHO.moms_rc(2),    rcconsHO.flux_moms(2));
    fluxtaus(dir)(p.I) =  calcflux(rcconsHO.lambda, rcconsHO.tau_rc,        rcconsHO.flux_tau);
    fluxDYes(dir)(p.I) =  calcflux(rcconsHO.lambda, rcconsHO.DYe_rc,        rcconsHO.flux_DYe);
    fluxBxs(dir)(p.I) =
             (dir != 0) * calcflux(rcconsHO.lambda, rcconsHO.Btildes_rc(0), rcconsHO.flux_Btildes(0));
    fluxBys(dir)(p.I) =
             (dir != 1) * calcflux(rcconsHO.lambda, rcconsHO.Btildes_rc(1), rcconsHO.flux_Btildes(1));
    fluxBzs(dir)(p.I) =
             (dir != 2) * calcflux(rcconsHO.lambda, rcconsHO.Btildes_rc(2), rcconsHO.flux_Btildes(2));

    /* Limit Fluxes */
    if (use_pplim || use_efl) {
      const auto Ip = p.I;
      const auto Im = p.I - p.DI[dir];

      CCTK_REAL theta = 1.0;
      CCTK_REAL theta_m = 1.0;
      CCTK_REAL theta_p = 1.0;

      // Calc LO flux, currently only LLF
      fluxLOdenss(dir)(Ip) = laxf(rcconsLO.lambda, rcconsLO.dens_rc,    rcconsLO.flux_dens);
      fluxLODEnts(dir)(Ip) = laxf(rcconsLO.lambda, rcconsLO.DEnt_rc,    rcconsLO.flux_DEnt);
      fluxLODYes(dir)(Ip) =  laxf(rcconsLO.lambda, rcconsLO.DYe_rc,     rcconsLO.flux_DYe);
      fluxLOmomxs(dir)(Ip) = laxf(rcconsLO.lambda, rcconsLO.moms_rc(0), rcconsLO.flux_moms(0));
      fluxLOmomys(dir)(Ip) = laxf(rcconsLO.lambda, rcconsLO.moms_rc(1), rcconsLO.flux_moms(1));
      fluxLOmomzs(dir)(Ip) = laxf(rcconsLO.lambda, rcconsLO.moms_rc(2), rcconsLO.flux_moms(2));
      fluxLOtaus(dir)(Ip) =  laxf(rcconsLO.lambda, rcconsLO.tau_rc,     rcconsLO.flux_tau);

      if (use_pplim) {
        // Get atmosphere
        CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        CCTK_REAL rho_atm = (radial_distance > r_atmo)
                      ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
                      : rho_abs_min;
        rho_atm = std::max(eos_3p->rgrho.min, rho_atm);

        // Get 2 * \alpha * CFL with \alpha = 3
        const CCTK_REAL a2cfl = 6 * cctk_delta_time / p.DX[dir];
        printf("DX = %e, a2cfl = %e", p.DX[dir], a2cfl);

        // Calc dens floor
        const CCTK_REAL densmin_p = volform(Ip) * w_lorentz(Ip) * rho_atm;
        const CCTK_REAL densmin_m = volform(Im) * w_lorentz(Im) * rho_atm;

        // Calc theta
        const CCTK_REAL newdens_p = dens(Ip) + a2cfl * fluxdenss(dir)(Ip);
        const CCTK_REAL newdens_m = dens(Im) - a2cfl * fluxdenss(dir)(Ip);

        const CCTK_REAL newdensLO_p = dens(Ip) + a2cfl * fluxLOdenss(dir)(Ip);
        const CCTK_REAL newdensLO_m = dens(Im) - a2cfl * fluxLOdenss(dir)(Ip);

        if (newdens_p < densmin_p)
          theta_p = std::min(
              theta, std::max(0.0, (densmin_p - newdensLO_p) /
                                       (a2cfl * (fluxdenss(dir)(Ip) -
                                                 fluxLOdenss(dir)(Ip)))));

        if (newdens_m < densmin_m)
          theta_m = std::min(
              theta, std::max(0.0, -(densmin_m - newdensLO_m) /
                                       (a2cfl * (fluxdenss(dir)(Ip) -
                                                 fluxLOdenss(dir)(Ip)))));

        theta = std::min(theta_m, theta_p);

        if (isnan(theta) || isnan(theta_m) || isnan(theta_p) || isnan(a2cfl) ||
            isnan(fluxLOdenss(dir)(Ip)) || isnan(fluxLODEnts(dir)(Ip)) ||
            isnan(fluxLOmomxs(dir)(Ip)) || isnan(fluxLOmomys(dir)(Ip)) ||
            isnan(fluxLOmomzs(dir)(Ip)) || isnan(fluxLOtaus(dir)(Ip)) ||
            isnan(fluxLOmomys(dir)(Ip)) || isnan(fluxLOmomzs(dir)(Ip)) ||
            isnan(densmin_p) || isnan(densmin_m)) {
          printf("cctk_iteration = %i,  dir = %i,  ijk = %i, %i, %i, "
                 "x, y, z = %16.8e, %16.8e, %16.8e.\n",
                 cctk_iteration, dir, p.i, p.j, p.k, p.x, p.y, p.z);
          printf("  fluxLOdenss = %16.8e,\n", fluxLOdenss(dir)(p.I));
          printf("  fluxLOmoms  = %16.8e, %16.8e, %16.8e,\n",
                 fluxLOmomxs(dir)(p.I), fluxLOmomys(dir)(p.I),
                 fluxLOmomzs(dir)(p.I));
          printf("  fluxLOtaus  = %16.8e,\n", fluxLOtaus(dir)(p.I));
          printf("  theta = %16.8e,\n", theta);
          printf("  a2cfl = %16.8e,\n", a2cfl);
          printf("  densmins = %16.8e, %16.8e\n", densmin_m, densmin_p);
          assert(0);
        }
      }

      if (use_efl) {
        // If use_ppl, theta = theta_ppl. Else, theta = 1.
        // This is the limiter from Doulis et al. 2022. efl_ce is typically 1.
        const CCTK_REAL nu_m = std::min(1.0, efl_ce * r_ent(Im));
        const CCTK_REAL nu_p = std::min(1.0, efl_ce * r_ent(Ip));
        theta = std::min(theta, 1 - 0.5 * (nu_m + nu_p));
      }

      // Update flux GF
      fluxdenss(dir)(Ip) =
          (1 - theta) * fluxLOdenss(dir)(Ip) + theta * fluxdenss(dir)(Ip);
      fluxDEnts(dir)(Ip) =
          (1 - theta) * fluxLODEnts(dir)(Ip) + theta * fluxDEnts(dir)(Ip);
      fluxDYes(dir)(Ip) =
          (1 - theta) * fluxLODYes(dir)(Ip) + theta * fluxDYes(dir)(Ip);
      fluxmomxs(dir)(Ip) =
          (1 - theta) * fluxLOmomxs(dir)(Ip) + theta * fluxmomxs(dir)(Ip);
      fluxmomys(dir)(Ip) =
          (1 - theta) * fluxLOmomys(dir)(Ip) + theta * fluxmomys(dir)(Ip);
      fluxmomzs(dir)(Ip) =
          (1 - theta) * fluxLOmomzs(dir)(Ip) + theta * fluxmomzs(dir)(Ip);
      fluxtaus(dir)(Ip) =
          (1 - theta) * fluxLOtaus(dir)(Ip) + theta * fluxtaus(dir)(Ip);
      thetagf(dir)(Ip) = theta;

      if (isnan(theta) || isnan(fluxLOdenss(dir)(Ip)) ||
          isnan(fluxLODEnts(dir)(Ip)) || isnan(fluxLOmomxs(dir)(Ip)) ||
          isnan(fluxLOmomys(dir)(Ip)) || isnan(fluxLOmomzs(dir)(Ip)) ||
          isnan(fluxLOtaus(dir)(Ip)) || isnan(fluxLOmomys(dir)(Ip)) ||
          isnan(fluxLOmomzs(dir)(Ip))) {
        printf("cctk_iteration = %i,  dir = %i,  ijk = %i, %i, %i, "
               "x, y, z = %16.8e, %16.8e, %16.8e.\n",
               cctk_iteration, dir, p.i, p.j, p.k, p.x, p.y, p.z);
        printf("  fluxLOdenss = %16.8e,\n", fluxLOdenss(dir)(p.I));
        printf("  fluxLOmoms  = %16.8e, %16.8e, %16.8e,\n",
               fluxLOmomxs(dir)(p.I), fluxLOmomys(dir)(p.I),
               fluxLOmomzs(dir)(p.I));
        printf("  fluxLOtaus  = %16.8e,\n", fluxLOtaus(dir)(p.I));
        printf("  theta = %16.8e,\n", theta);
        assert(0);
      }
    } else {
      fluxLOdenss(dir)(p.I) = 0.0;
      fluxLODEnts(dir)(p.I) = 0.0;
      fluxLOmomxs(dir)(p.I) = 0.0;
      fluxLOmomys(dir)(p.I) = 0.0;
      fluxLOmomzs(dir)(p.I) = 0.0;
      fluxLOtaus(dir)(p.I) = 0.0;
      thetagf(dir)(p.I) = 1.0;
    }

    /* TODO: def NaN check in recon structs
    if (isnan(dens_rc(0)) || isnan(dens_rc(1)) || isnan(moms_rc(0)(0)) ||
        isnan(moms_rc(0)(1)) || isnan(moms_rc(1)(0)) || isnan(moms_rc(1)(1)) ||
        isnan(moms_rc(2)(0)) || isnan(moms_rc(2)(1)) || isnan(tau_rc(0)) ||
        isnan(tau_rc(1)) || isnan(Btildes_rc(0)(0)) ||
        isnan(Btildes_rc(0)(1)) || isnan(Btildes_rc(1)(0)) ||
        isnan(Btildes_rc(1)(1)) || isnan(Btildes_rc(2)(0)) ||
        isnan(Btildes_rc(2)(1)) || isnan(DYe_rc(0)) || isnan(DYe_rc(1)) ||
        isnan(flux_dens(0)) || isnan(flux_dens(1)) || isnan(flux_moms(0)(0)) ||
        isnan(flux_moms(0)(1)) || isnan(flux_moms(1)(0)) ||
        isnan(flux_moms(1)(1)) || isnan(flux_moms(2)(0)) ||
        isnan(flux_moms(2)(1)) || isnan(flux_DYe(0)) || isnan(flux_DYe(1)) ||
        isnan(flux_tau(0)) || isnan(flux_tau(1)) || isnan(flux_Btildes(0)(0)) ||
        isnan(flux_Btildes(0)(1)) || isnan(flux_Btildes(1)(0)) ||
        isnan(flux_Btildes(1)(1)) || isnan(flux_Btildes(2)(0)) ||
        isnan(flux_Btildes(2)(1)) || isnan(fluxdenss(dir)(p.I)) ||
        isnan(fluxmomxs(dir)(p.I)) || isnan(fluxmomys(dir)(p.I)) ||
        isnan(fluxmomzs(dir)(p.I)) || isnan(fluxtaus(dir)(p.I)) ||
        isnan(fluxBxs(dir)(p.I)) || isnan(fluxBys(dir)(p.I)) ||
        isnan(fluxBzs(dir)(p.I)) || rho_rc(0) < 0.0 || rho_rc(1) < 0.0 ||
        press_rc(0) < 0.0 || press_rc(1) < 0.0) {
      */
    if (isnan(rcconsHO.dens_rc(0)) || isnan(rcconsHO.dens_rc(1)) || isnan(rcconsHO.moms_rc(0)(0))) {
      printf("cctk_iteration = %i,  dir = %i,  ijk = %i, %i, %i, "
             "x, y, z = %16.8e, %16.8e, %16.8e.\n",
             cctk_iteration, dir, p.i, p.j, p.k, p.x, p.y, p.z);
      // printf("  fluxdenss = %16.8e,\n", fluxdenss(dir)(p.I));
      // printf("  fluxmoms  = %16.8e, %16.8e, %16.8e,\n", fluxmomxs(dir)(p.I),
      //        fluxmomys(dir)(p.I), fluxmomzs(dir)(p.I));
      // printf("  fluxtaus  = %16.8e,\n", fluxtaus(dir)(p.I));
      // printf("  fluxBs    = %16.8e, %16.8e, %16.8e\n", fluxBxs(dir)(p.I),
      //        fluxBys(dir)(p.I), fluxBzs(dir)(p.I));
      // printf("  flux_denss = %16.8e, %16.8e,\n", flux_dens(0), flux_dens(1));
      // printf("  flux_moms  = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e,\n",
      //        flux_moms(0)(0), flux_moms(0)(1), flux_moms(1)(0), flux_moms(1)(1),
      //        flux_moms(2)(0), flux_moms(2)(1));
      // printf("  flux_taus  = %16.8e, %16.8e,\n", flux_tau(0), flux_tau(1));
      // printf("  flux_DYes  = %16.8e, %16.8e,\n", flux_DYe(0), flux_DYe(1));
      // printf("  flux_Bts   = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e,\n",
      //        flux_Btildes(0)(0), flux_Btildes(0)(1), flux_Btildes(1)(0),
      //        flux_Btildes(1)(1), flux_Btildes(2)(0), flux_Btildes(2)(1));
      // printf("  dens_rc = %16.8e, %16.8e,\n", dens_rc(0), dens_rc(1));
      // printf("  moms_rc = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e,\n",
      //        moms_rc(0)(0), moms_rc(0)(1), moms_rc(1)(0), moms_rc(1)(1),
      //        moms_rc(2)(0), moms_rc(2)(1));
      // printf("  tau_rc  = %16.8e, %16.8e,\n", tau_rc(0), tau_rc(1));
      // printf("  DYe_rc  = %16.8e, %16.8e,\n", DYe_rc(0), DYe_rc(1));
      // printf("  Bs_rc  = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e,\n",
      //        Bs_rc(0)(0), Bs_rc(0)(1), Bs_rc(1)(0), Bs_rc(1)(1), Bs_rc(2)(0),
      //        Bs_rc(2)(1));
      // printf("  Bts_rc  = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e,\n",
      //        Btildes_rc(0)(0), Btildes_rc(0)(1), Btildes_rc(1)(0),
      //        Btildes_rc(1)(1), Btildes_rc(2)(0), Btildes_rc(2)(1));
      // printf("  lam = %16.8e, %16.8e, %16.8e, %16.8e,\n"
      //        "        %16.8e, %16.8e, %16.8e, %16.8e.\n",
      //        lambda(0)(0), lambda(0)(1), lambda(0)(2), lambda(0)(3),
      //        lambda(1)(0), lambda(1)(1), lambda(1)(2), lambda(1)(3));
      // printf("  alp_avg = %16.8e, beta_avg = %16.8e, u_avg = %16.8e \n",
      //        alp_avg, beta_avg, u_avg);
      // printf("  vel_rc  = %16.8e, %16.8e \n", vel_rc(0), vel_rc(1));
      // printf("  rho_rc  = %16.8e, %16.8e \n", rho_rc(0), rho_rc(1));
      // printf("  cs2_rc  = %16.8e, %16.8e \n", cs2_rc(0), cs2_rc(1));
      // printf("  wlor_rc = %16.8e, %16.8e \n", w_lorentz_rc(0), w_lorentz_rc(1));
      // printf("  h_rc    = %16.8e, %16.8e \n", h_rc(0), h_rc(1));
      // printf("  bsq_rc  = %16.8e, %16.8e \n", bsq_rc(0), bsq_rc(1));
      // printf("  press_rc = %16.8e, %16.8e \n", press_rc(0), press_rc(1));
      // printf("  eps_rc   = %16.8e, %16.8e \n", eps_rc(0), eps_rc(1));
      printf("  rho = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e;\n",
             rho(p.I - p.DI[dir] * 3), rho(p.I - p.DI[dir] * 2),
             rho(p.I - p.DI[dir]), rho(p.I), rho(p.I + p.DI[dir]),
             rho(p.I + p.DI[dir] * 2));
      printf("  press = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e;\n",
             press(p.I - p.DI[dir] * 3), press(p.I - p.DI[dir] * 2),
             press(p.I - p.DI[dir]), press(p.I), press(p.I + p.DI[dir]),
             press(p.I + p.DI[dir] * 2));
      printf("  eps   = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e;\n",
             eps(p.I - p.DI[dir] * 3), eps(p.I - p.DI[dir] * 2),
             eps(p.I - p.DI[dir]), eps(p.I), eps(p.I + p.DI[dir]),
             eps(p.I + p.DI[dir] * 2));
      printf("  alp_avg, beta_avg = %16.8e, %16.8e, %16.8e, %16.8e,\n", alp_avg,
             betas_avg(0), betas_avg(1), betas_avg(2));
      printf("  g_avg = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e.\n",
             g_avg(0, 0), g_avg(0, 1), g_avg(0, 2), g_avg(1, 1), g_avg(1, 2),
             g_avg(2, 2));
      printf("  sqrtg = %16.8e,\n", sqrtg);
      // printf("  vlows_rc  = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e.\n",
      //        vlows_rc(0)(0), vlows_rc(0)(1), vlows_rc(1)(0), vlows_rc(1)(1),
      //        vlows_rc(2)(0), vlows_rc(2)(1));
      // printf("  vups_rc   = %16.8e, %16.8e, %16.8e, %16.8e, %16.8e, %16.8e.\n",
      //        vels_rc(0)(0), vels_rc(0)(1), vels_rc(1)(0), vels_rc(1)(1),
      //        vels_rc(2)(0), vels_rc(2)(1));
      // printf("  vtilde_rc = %16.8e, %16.8e.\n", vtilde_rc(0), vtilde_rc(1));
      assert(0);
    }

    /* Begin code for upwindCT */
    // if dir==0: dir1=1, dir2=2 | dir==1: dir1=2, dir2=0 | dir==2; dir1=0,
    // dir2=1

    const int dir1 = (dir == 0) ? 1 : ((dir == 1) ? 2 : 0);
    const int dir2 = (dir == 0) ? 2 : ((dir == 1) ? 0 : 1);

    vec<vec<CCTK_REAL, 4>, 2> lambda = rcconsHO.lambda;
    vec<vec<CCTK_REAL, 2>, 3> vtildes_rc = rcconsHO.vtildes_rc;

    amax(dir)(p.I) = max({CCTK_REAL(0), lambda(0)(0), lambda(0)(1),
                          lambda(0)(2), lambda(0)(3), lambda(1)(0),
                          lambda(1)(1), lambda(1)(2), lambda(1)(3)});

    amin(dir)(p.I) = -1 * (min({CCTK_REAL(0), lambda(0)(0), lambda(0)(1),
                                lambda(0)(2), lambda(0)(3), lambda(1)(0),
                                lambda(1)(1), lambda(1)(2), lambda(1)(3)}));

    vtildes_one(dir)(p.I) = (amax(dir)(p.I) * vtildes_rc(dir1)(0) +
                             amin(dir)(p.I) * vtildes_rc(dir1)(1)) /
                            (amax(dir)(p.I) + amin(dir)(p.I));
    vtildes_two(dir)(p.I) = (amax(dir)(p.I) * vtildes_rc(dir2)(0) +
                             amin(dir)(p.I) * vtildes_rc(dir2)(1)) /
                            (amax(dir)(p.I) + amin(dir)(p.I));

    /* End code for upwindCT */
  });
}

void CalcAuxForAvecPsi(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  const vec<GF3D2<const CCTK_REAL>, dim> gf_Avecs{Avec_x, Avec_y, Avec_z};
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* interpolate A to vertices */
        const vec<CCTK_REAL, 3> A_vert([&](int i) ARITH_INLINE {
          return calc_avg_e2v(gf_Avecs(i), p, i);
        });
        const smat<CCTK_REAL, 3> g{gxx(p.I), gxy(p.I), gxz(p.I),
                                   gyy(p.I), gyz(p.I), gzz(p.I)};
        const vec<CCTK_REAL, 3> betas{betax(p.I), betay(p.I), betaz(p.I)};
        const CCTK_REAL detg = calc_det(g);
        const CCTK_REAL sqrtg = sqrt(detg);
        const smat<CCTK_REAL, 3> ug = calc_inv(g, detg);
        const vec<CCTK_REAL, 3> Aup = calc_contraction(ug, A_vert);

        Fx(p.I) = alp(p.I) * sqrtg * Aup(0);
        Fy(p.I) = alp(p.I) * sqrtg * Aup(1);
        Fz(p.I) = alp(p.I) * sqrtg * Aup(2);
        Fbetax(p.I) = betas(0) * Psi(p.I);
        Fbetay(p.I) = betas(1) * Psi(p.I);
        Fbetaz(p.I) = betas(2) * Psi(p.I);
        G(p.I) = alp(p.I) * Psi(p.I) / sqrtg - calc_contraction(betas, A_vert);
      });
}

extern "C" void AsterX_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Fluxes;
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

    CalcFlux<0>(cctkGH, eos_3p_ig);
    CalcFlux<1>(cctkGH, eos_3p_ig);
    CalcFlux<2>(cctkGH, eos_3p_ig);
    break;
  }
  case eos_3param::Hybrid: {
    // Get local eos object
    auto eos_3p_hyb = global_eos_3p_hyb;

    CalcFlux<0>(cctkGH, eos_3p_hyb);
    CalcFlux<1>(cctkGH, eos_3p_hyb);
    CalcFlux<2>(cctkGH, eos_3p_hyb);
    break;
  }
  case eos_3param::Tabulated: {
    // Get local eos object
    auto eos_3p_tab3d = global_eos_3p_tab3d;

    CalcFlux<0>(cctkGH, eos_3p_tab3d);
    CalcFlux<1>(cctkGH, eos_3p_tab3d);
    CalcFlux<2>(cctkGH, eos_3p_tab3d);
    break;
  }
  default:
    assert(0);
  }

  /* Set auxiliary variables for the rhs of A and Psi  */
  CalcAuxForAvecPsi(cctkGH);
}

} // namespace AsterX
