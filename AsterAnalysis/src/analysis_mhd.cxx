#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <limits>

#include "aster_utils.hxx"

namespace AsterAnalysis {
using namespace std;
using namespace Loop;
using namespace AsterUtils;

extern "C" void AsterAnalysis_MHD(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterAnalysis_MHD;
  DECLARE_CCTK_PARAMETERS;

  const bool use_v_vec = CCTK_EQUALS(recon_type, "v_vec");

  /* handy aliases for grid functions */
  const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  const vec<GF3D2<const CCTK_REAL>, 3> gf_Avec{Avec_x, Avec_y, Avec_z};
  const vec<GF3D2<const CCTK_REAL>, 3> gf_Bvec{Bvecx, Bvecy, Bvecz};
  const vec<GF3D2<const CCTK_REAL>, 3> gf_vel{velx, vely, velz};

  /* inverse grid spacings (assumed uniform) */
  const CCTK_REAL idx = 1.0 / CCTK_DELTA_SPACE(0);
  const CCTK_REAL idy = 1.0 / CCTK_DELTA_SPACE(1);
  const CCTK_REAL idz = 1.0 / CCTK_DELTA_SPACE(2);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* ---- metric averages at cell center ---- */
        const CCTK_REAL alp_avg = calc_avg_v2c(alp, p);
        const vec<CCTK_REAL, 3> beta_avg(
            [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });
        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });
        const CCTK_REAL detg = calc_det(g_avg);
        const CCTK_REAL sqrt_detg = sqrt(detg);

        const smat<CCTK_REAL, 3> ug_avg = calc_inv(g_avg, detg);

        /* compute w_lorentz */
        vec<CCTK_REAL, 3> v_up;
        vec<CCTK_REAL, 3> v_low;
        if (use_v_vec) {

          v_up(0) = velx(p.I);
          v_up(1) = vely(p.I);
          v_up(2) = velz(p.I);
          v_low = calc_contraction(g_avg, v_up);

        } else {

          const vec<CCTK_REAL, 3> z_up{zvec_x(p.I), zvec_y(p.I), zvec_z(p.I)};
          const vec<CCTK_REAL, 3> z_low = calc_contraction(g_avg, z_up);
          v_up = z_up / w_lorentz(p.I);
          v_low = z_low / w_lorentz(p.I);
        }

        const CCTK_REAL v2 = calc_contraction(v_low, v_up);

        /* v^i - beta^i / alpha */
        const vec<CCTK_REAL, 3> velshift = v_up - beta_avg / alp_avg;

        /* ---- magnetic field ---- */
        const vec<CCTK_REAL, 3> B_up{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
        const vec<CCTK_REAL, 3> B_low = calc_contraction(g_avg, B_up);
        const CCTK_REAL B2big = calc_contraction(B_up, B_low);

        /* cell-centered A^i and A_i */
        const vec<CCTK_REAL, 3> A_low{calc_avg_e2c<0>(gf_Avec(0), p),
                                      calc_avg_e2c<1>(gf_Avec(1), p),
                                      calc_avg_e2c<2>(gf_Avec(2), p)};

        const vec<CCTK_REAL, 3> A_up = calc_contraction(ug_avg, A_low);

        A_norm(p.I) = sqrt(calc_contraction(A_low, A_up));

        /* b^mu: b^0 and b^i */
        const CCTK_REAL bs0 =
            w_lorentz(p.I) * calc_contraction(B_up, v_low) / alp_avg;
        const vec<CCTK_REAL, 3> bs =
            (B_up / w_lorentz(p.I) + alp_avg * bs0 * velshift);

        /* inverse Beta and magnetization */
        mhd_press_ratio(p.I) = b2small(p.I) / (2.0 * press(p.I));
        mhd_energy_ratio(p.I) = b2small(p.I) / (2.0 * rho(p.I));

        /* EM energy integrand: (B^2 - 0.5 b^2) * √γ */
        magnetic_energy(p.I) = (B2big - 0.5 * b2small(p.I)) * sqrt_detg;

        /* Kinetic energy: 0.5*(ρ + ρ ε + p) v^2 * W * √γ */
        const CCTK_REAL rhoH = rho(p.I) + rho(p.I) * eps(p.I) + press(p.I);
        kinetic_energy(p.I) = 0.5 * rhoH * v2 * w_lorentz(p.I) * sqrt_detg;

        /* ⟨|B|⟩ */
        mean_B_numerator(p.I) =
            rho(p.I) * B_norm(p.I) * w_lorentz(p.I) * sqrt_detg;
        mean_B_denominator(p.I) = rho(p.I) * w_lorentz(p.I) * sqrt_detg;

        // NOTE: Computing volume integral of mean_B must be done via
        // post-processing To compute VI of mean_B, use sum column values in the
        // norms files of mean_B_num and mean_B_den e.g.: mean_B = mean_B_num /
        // mean_B_den

        /* ---- shift & four‑velocity ---- */
        const vec<CCTK_REAL, 3> beta_avg_low =
            calc_contraction(g_avg, beta_avg);
        const CCTK_REAL beta2 = calc_contraction(beta_avg, beta_avg_low);

        const vec<CCTK_REAL, 3> u = velshift * w_lorentz(p.I); /* u^i */
        const CCTK_REAL utlow =
            (-alp_avg + calc_contraction(v_up, beta_avg_low)) *
            w_lorentz(p.I); /* u_0 */

        /* b_i & b_0 */
        const vec<CCTK_REAL, 3> bs_low =
            bs0 * beta_avg_low + calc_contraction(g_avg, bs);
        const CCTK_REAL bstlow =
            bs0 * (-pow2(alp_avg) + beta2) + calc_contraction(bs, beta_avg_low);

        /* ---- Poynting flux S^i ---- */
        const vec<CCTK_REAL, 3> S =
            -alp_avg * (b2small(p.I) * u * utlow - bs * bstlow);
        poynting_flux_x(p.I) = S(0);
        poynting_flux_y(p.I) = S(1);
        poynting_flux_z(p.I) = S(2);
        const vec<CCTK_REAL, 3> S_low = calc_contraction(g_avg, S);
        poynting_flux_norm(p.I) = calc_contraction(S, S_low);

        /* ---- alternative Poynting (E x B) ---- */
        const vec<CCTK_REAL, 3> vtilde_low = alp_avg * v_low - beta_avg_low;
        const vec<CCTK_REAL, 3> E{
            -(vtilde_low(1) * B_low(2) - vtilde_low(2) * B_low(1)),
            -(vtilde_low(2) * B_low(0) - vtilde_low(0) * B_low(2)),
            -(vtilde_low(0) * B_low(1) - vtilde_low(1) * B_low(0))};
        const vec<CCTK_REAL, 3> S_alt{E(1) * B_up(2) - E(2) * B_up(1),
                                      E(2) * B_up(0) - E(0) * B_up(2),
                                      E(0) * B_up(1) - E(1) * B_up(0)};
        poynting_vector_x(p.I) = S_alt(0);
        poynting_vector_y(p.I) = S_alt(1);
        poynting_vector_z(p.I) = S_alt(2);

        /* normal to sphere */
        const CCTK_REAL rmag =
            max(sqrt(p.x * p.x + p.y * p.y + p.z * p.z), 1e-20);
        const vec<CCTK_REAL, 3> n{p.x / rmag, p.y / rmag, p.z / rmag};
        poynting_scalar(p.I) = calc_contraction(S_alt, n);

        /* ---- Lambda_MRI ---- */
        const CCTK_REAL omega_x = alp_avg * velx(p.I) - beta_avg(0);
        const CCTK_REAL omega_y = alp_avg * vely(p.I) - beta_avg(1);
        CCTK_REAL local_omega = 0.0;
        if (pow2(p.x) + pow2(p.y) > 0.0) {
          local_omega =
              (p.x * omega_y - p.y * omega_x) / (pow2(p.x) + pow2(p.y));
        }

        const CCTK_REAL Wtot =
            rho(p.I) * (1.0 + eps(p.I)) + press(p.I) + b2small(p.I);
        auto lambda_MRI = [&](CCTK_REAL bi, CCTK_REAL bi_low,
                              CCTK_REAL &lambda) {
          const CCTK_REAL va = sqrt(bi * bi_low / Wtot);
          lambda = (va > numeric_limits<CCTK_REAL>::epsilon() &&
                    fabs(local_omega) > 1.0e-15)
                       ? 2.0 * M_PI * fabs(va / local_omega)
                       : 0.0;
        };
        lambda_MRI(bs(0), bs_low(0), lambda_MRI_x(p.I));
        lambda_MRI(bs(1), bs_low(1), lambda_MRI_y(p.I));
        lambda_MRI(bs(2), bs_low(2), lambda_MRI_z(p.I));

        /* ---- divB (CT, densitized faces) ---- */
        const CCTK_REAL divB_val =
            idx * (dBx_stag(p.I + p.DI[0]) - dBx_stag(p.I)) +
            idy * (dBy_stag(p.I + p.DI[1]) - dBy_stag(p.I)) +
            idz * (dBz_stag(p.I + p.DI[2]) - dBz_stag(p.I));
        divB(p.I) = divB_val / sqrt_detg;

        /* ---- divA (flat divergence of edge-centered A) ---- */
        const CCTK_REAL divA_val = idx * (Avec_x(p.I) - Avec_x(p.I - p.DI[0])) +
                                   idy * (Avec_y(p.I) - Avec_y(p.I - p.DI[1])) +
                                   idz * (Avec_z(p.I) - Avec_z(p.I - p.DI[2]));
        divA(p.I) = divA_val;
      }); /* end grid loop */
}

extern "C" void AsterAnalysis_HOFV(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterAnalysis_HOFV;
  DECLARE_CCTK_PARAMETERS;

  /* handy aliases for grid functions */
  // const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{betax, betay, betaz};
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  // const vec<GF3D2<const CCTK_REAL>, 3> gf_Avec{Avec_x, Avec_y, Avec_z};
  // const vec<GF3D2<const CCTK_REAL>, 3> gf_Bvec{Bvecx, Bvecy, Bvecz};
  // const vec<GF3D2<const CCTK_REAL>, 3> gf_vel{velx, vely, velz};

  /* inverse grid spacings (assumed uniform) */
  const CCTK_REAL idx = 1.0 / CCTK_DELTA_SPACE(0);
  const CCTK_REAL idy = 1.0 / CCTK_DELTA_SPACE(1);
  const CCTK_REAL idz = 1.0 / CCTK_DELTA_SPACE(2);

  constexpr CCTK_REAL one_over_24 = CCTK_REAL(1)/CCTK_REAL(24);

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* ---- metric averages at cell center ---- */
        // const CCTK_REAL alp_avg = calc_avg_v2c(alp, p);
        // const vec<CCTK_REAL, 3> beta_avg(
        //     [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });
        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });
        const CCTK_REAL detg = calc_det(g_avg);
        const CCTK_REAL sqrt_detg = sqrt(detg);

        /* ---- divB (CT, densitized faces) ---- */
        const CCTK_REAL divB_val =
            idx * (dBx_stag_fa(p.I + p.DI[0]) - dBx_stag_fa(p.I)) +
            idy * (dBy_stag_fa(p.I + p.DI[1]) - dBy_stag_fa(p.I)) +
            idz * (dBz_stag_fa(p.I + p.DI[2]) - dBz_stag_fa(p.I));
        divB_FA(p.I) = divB_val / sqrt_detg;

        divB_fromPV(p.I) = calc_FV_flux<0>(dBx_stag, p) 
          + calc_FV_flux<1>(dBy_stag, p) 
          + calc_FV_flux<2>(dBz_stag, p);

        bool thetac = (!LOflag(p.I) || !shock_pv_fallback);
        const CCTK_REAL tau_reavg = tau_pv(p.I) + thetac * one_over_24*laplace_3d(tau_pv,p);
        tau_avg_err(p.I) = abs(tau(p.I) - tau_reavg) / tau(p.I);

        flag_weight(p.I) = (!LOflag(p.I)) * dens(p.I);

      }); /* end grid loop */
}
} /* namespace AsterAnalysis */
