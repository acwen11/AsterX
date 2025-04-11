#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "eos.hxx"
#include "eos_idealgas.hxx"

#include "c2p.hxx"
#include "c2p_1DPalenzuela.hxx"
#include "c2p_2DNoble.hxx"
#include "c2p_1DEntropy.hxx"

#include "c2p_utils.hxx"

namespace Con2PrimFactory {

using namespace Arith;
using namespace EOSX;
using namespace AsterUtils;

extern "C" void Con2PrimFactory_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  // Initializing the IG EOS
  eos::range rgeps(0, 100), rgrho(1e-13, 20), rgye(0.5, 0.5);
  const eos_idealgas eos_th(2.0, 938.985, rgeps, rgrho, rgye);

  // Set atmo values
  const CCTK_REAL rho_atmo = 1e-10;
  const CCTK_REAL eps_atmo = 1e-8;
  const CCTK_REAL Ye_atmo = 0.5;
  const CCTK_REAL press_atmo = eos_th.press_from_valid_rho_eps_ye(rho_atmo,eps_atmo,Ye_atmo);
  const CCTK_REAL entropy_atmo = eos_th.kappa_from_valid_rho_eps_ye(rho_atmo,eps_atmo,Ye_atmo);

  // Setting up atmosphere
  const CCTK_REAL rho_atmo_cut = rho_atmo * (1 + 1.0e-3);
  atmosphere atmo(rho_atmo, eps_atmo, Ye_atmo, press_atmo, entropy_atmo, rho_atmo_cut);

  // Metric
  const smat<CCTK_REAL, 3> g{1.0, 0.0, 0.0,
                             1.0, 0.0, 1.0}; // xx, xy, xz, yy, yz, zz

  // Set BH limiters
  const CCTK_REAL alp_thresh = -1.;
  const CCTK_REAL rho_BH = 1e20;
  const CCTK_REAL eps_BH = 1e20;
  const CCTK_REAL vwlim_BH = 1e20;

  // Con2Prim objects
  // (eos_th, atmo, max_iter, c2p_tol
  //  alp_thresh, cons_error_limit,
  //  vw_lim, B_lim, rho_BH, eps_BH, vwlim_BH,
  //  Ye_lenient, use_z)
  c2p_2DNoble c2p_Noble(eos_th, atmo, 100, 1e-8, alp_thresh, -1., 1, 1, rho_BH, eps_BH, vwlim_BH, true, false);
  c2p_1DPalenzuela c2p_Pal(eos_th, atmo, 100, 1e-8, alp_thresh, -1., 1, 1, rho_BH, eps_BH, vwlim_BH, true, false);
  c2p_1DEntropy c2p_Ent(eos_th, atmo, 100, 1e-8, alp_thresh, -1., 1, 1, rho_BH, eps_BH, vwlim_BH, true, false);

  // Construct error report object:
  c2p_report rep_Noble;
  c2p_report rep_Pal;
  c2p_report rep_Ent;

  // Set primitive seeds
  const CCTK_REAL rho_in = 0.125;
  const CCTK_REAL eps_in = 0.8;
  const CCTK_REAL Ye_in = 0.5;
  const CCTK_REAL press_in = eos_th.press_from_valid_rho_eps_ye(rho_in,eps_in,Ye_in);
  const CCTK_REAL entropy_in = eos_th.kappa_from_valid_rho_eps_ye(rho_in,eps_in,Ye_in);
  const vec<CCTK_REAL, 3> vup_in = {0.0,0.0,0.0};
  const vec<CCTK_REAL, 3> Bup_in = {0.5,-1.0,0.0};
  const vec<CCTK_REAL, 3> vdown_in = calc_contraction(g,vup_in);
  const CCTK_REAL wlor_in = calc_wlorentz(vdown_in,vup_in); 

  prim_vars pv;
  // rho(p.I), eps(p.I), dummy_Ye, press(p.I), entropy, v_up, wlor, Bup
  prim_vars pv_seeds{rho_in, eps_in, Ye_in, press_in, entropy_in, vup_in, wlor_in, Bup_in};

  // cons_vars cv{dens(p.I), {momx(p.I), momy(p.I), momz(p.I)}, tau(p.I),
  //  dummy_dYe, DEnt, {dBx(p.I), dBy(p.I), dBz(p.I)}};

  cons_vars cv_Noble;
  cons_vars cv_Pal;
  cons_vars cv_Ent;
  cons_vars cv_all;

  // C2P solve may modify cons_vars input
  // Use the following to test C2Ps independently
  cv_Noble.from_prim(pv_seeds, g);
  cv_Pal.from_prim(pv_seeds, g);
  cv_Ent.from_prim(pv_seeds, g);

  // Use the following to test C2Ps in sequence
  // (As done in AsterX, the evolution thorn)
  cv_all.from_prim(pv_seeds, g);

  // Testing C2P Noble
  CCTK_VINFO("Testing C2P Noble...");
  c2p_Noble.solve(eos_th, pv, pv_seeds, cv_all, g, rep_Noble);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.entropy, pv.entropy, 
         pv_seeds.vel(0), pv.vel(0), pv_seeds.vel(1),
         pv.vel(1), pv_seeds.vel(2), pv.vel(2), pv_seeds.Bvec(0), pv.Bvec(0),
         pv_seeds.Bvec(1), pv.Bvec(1), pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "dYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n", 
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2), cv_all.dYe, cv_all.dBvec(0),
         cv_all.dBvec(1), cv_all.dBvec(2), cv_all.DEnt);
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */

  rep_Noble.debug_message();

  // Testing C2P Palenzuela
  CCTK_VINFO("Testing C2P Palenzuela...");
  //c2p_Pal.solve(eos_th, pv, pv_seeds, cv, g, rep_Pal);
  c2p_Pal.solve(eos_th, pv, cv_all, g, rep_Pal);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.entropy, pv.entropy, 
         pv_seeds.vel(0), pv.vel(0), pv_seeds.vel(1),
         pv.vel(1), pv_seeds.vel(2), pv.vel(2), pv_seeds.Bvec(0), pv.Bvec(0),
         pv_seeds.Bvec(1), pv.Bvec(1), pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "dYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n", 
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2), cv_all.dYe, cv_all.dBvec(0),
         cv_all.dBvec(1), cv_all.dBvec(2), cv_all.DEnt);
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */
  rep_Pal.debug_message();

  // Testing C2P Entropy
  CCTK_VINFO("Testing C2P Entropy...");
  c2p_Ent.solve(eos_th, pv, cv_all, g, rep_Ent);

  printf("pv_seeds, pv: \n"
         "rho: %f, %f \n"
         "eps: %f, %f \n"
         "Ye: %f, %f \n"
         "press: %f, %f \n"
         "entropy: %f, %f \n"
         "velx: %f, %f \n"
         "vely: %f, %f \n"
         "velz: %f, %f \n"
         "Bx: %f, %f \n"
         "By: %f, %f \n"
         "Bz: %f, %f \n",
         pv_seeds.rho, pv.rho, pv_seeds.eps, pv.eps, pv_seeds.Ye, pv.Ye,
         pv_seeds.press, pv.press, pv_seeds.entropy, pv.entropy, 
         pv_seeds.vel(0), pv.vel(0), pv_seeds.vel(1),
         pv.vel(1), pv_seeds.vel(2), pv.vel(2), pv_seeds.Bvec(0), pv.Bvec(0),
         pv_seeds.Bvec(1), pv.Bvec(1), pv_seeds.Bvec(2), pv.Bvec(2));
  printf("cv: \n"
         "dens: %f \n"
         "tau: %f \n"
         "momx: %f \n"
         "momy: %f \n"
         "momz: %f \n"
         "dYe: %f \n"
         "dBx: %f \n"
         "dBy: %f \n"
         "dBz: %f \n"
         "DEnt: %f \n", 
         cv_all.dens, cv_all.tau, cv_all.mom(0), cv_all.mom(1), cv_all.mom(2), cv_all.dYe, cv_all.dBvec(0),
         cv_all.dBvec(1), cv_all.dBvec(2), cv_all.DEnt);
  /*
    assert(pv.rho == pv_seeds.rho);
    assert(pv.eps == pv_seeds.eps);
    assert(pv.press == pv_seeds.press);
    assert(pv.vel == pv_seeds.vel);
    assert(pv.Bvec == pv_seeds.Bvec);
  */
  rep_Ent.debug_message();
}

} // namespace Con2PrimFactory
