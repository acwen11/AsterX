#ifndef RECONX_RECONSTRUCT_HXX
#define RECONX_RECONSTRUCT_HXX

#include <loop_device.hxx>

#include <cctk.h>
// #include <cctk_Arguments.h>
// #include <cctk_Parameters.h>

#include "monocentral.hxx"
#include "minmod.hxx"
#include "ppm.hxx"
#include "eppm.hxx"
#include "weno5.hxx"
#include "wenoz.hxx"
#include "wenozp.hxx"
#include "mp5.hxx"

#include <type_traits>
#include <array>

namespace ReconX {

using std::array;
using namespace Arith;

// enum class for different reconstruction routines

enum class reconstruction_t {
  Godunov,
  minmod,
  monocentral,
  ppm,
  eppm,
  weno5,
  wenoz,
  wenozp,
  mp5
};

template <typename Container = std::array<CCTK_REAL, 2>>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE Container
reconstruct(const GF3D2<const CCTK_REAL> &gf_var, const PointDesc &p,
            reconstruction_t reconstruction, int dir,
            bool gf_is_rho, bool gf_is_press,
            const GF3D2<const CCTK_REAL> &gf_press,
            const GF3D2<const CCTK_REAL> &gf_vel_dir,
            const reconstruct_params_t &reconstruct_params) {

  // Neighbouring "plus" and "minus" cell indices
  const auto Immm = p.I - 3 * p.DI[dir];
  const auto Imm = p.I - 2 * p.DI[dir];
  const auto Im = p.I - p.DI[dir];
  const auto Ip = p.I;
  const auto Ipp = p.I + p.DI[dir];
  const auto Ippp = p.I + 2 * p.DI[dir];

  // canonical compute buffer
  std::array<CCTK_REAL,2> tmp{};

  switch (reconstruction) {

  case reconstruction_t::Godunov: {
    tmp = {gf_var(Im), gf_var(Ip)};
    break;
  }

  case reconstruction_t::minmod: {
    tmp = minmod_reconstruct(gf_var(Imm), gf_var(Im), gf_var(Ip), gf_var(Ipp));
    break;
  }

  case reconstruction_t::monocentral: {
    tmp = monocentral_reconstruct(gf_var(Imm), gf_var(Im), gf_var(Ip),
                                   gf_var(Ipp));
    break;
  }

  case reconstruction_t::ppm: {
    tmp = ppm_reconstruct(
        gf_var(Immm), gf_var(Imm), gf_var(Im), gf_var(Ip), gf_var(Ipp),
        gf_var(Ippp), gf_press(Immm), gf_press(Imm), gf_press(Im), gf_press(Ip),
        gf_press(Ipp), gf_press(Ippp), gf_vel_dir(Imm), gf_vel_dir(Im),
        gf_vel_dir(Ip), gf_vel_dir(Ipp), gf_is_rho, reconstruct_params);
    break;
  }

  case reconstruction_t::weno5: {
    tmp = weno5_reconstruct(gf_var(Immm), gf_var(Imm), gf_var(Im), gf_var(Ip),
                             gf_var(Ipp), gf_var(Ippp),
                             reconstruct_params.weno_eps);
    break;
  }

  case reconstruction_t::wenoz: {
    tmp = wenoz_reconstruct(gf_var(Immm), gf_var(Imm), gf_var(Im), gf_var(Ip),
                             gf_var(Ipp), gf_var(Ippp),
                             reconstruct_params.weno_eps);
    break;
  }

  case reconstruction_t::wenozp: {
		// const CCTK_REAL dx = p.DX[dir];
		const CCTK_REAL dx = 0.0;
    tmp = wenozp_reconstruct(gf_var(Immm), gf_var(Imm), gf_var(Im), gf_var(Ip),
                             gf_var(Ipp), gf_var(Ippp), dx,
                             reconstruct_params.weno_eps, reconstruct_params.weno_mp);
    break;
  }

  case reconstruction_t::mp5: {
    tmp = mp5_reconstruct(gf_var(Immm), gf_var(Imm), gf_var(Im), gf_var(Ip),
                           gf_var(Ipp), gf_var(Ippp),
                           reconstruct_params.mp5_alpha);
    break;
  }

  case reconstruction_t::eppm: {
    const array<const vect<int, dim>, 5> cells_Im = {Immm, Imm, Im, Ip, Ipp};
    const array<const vect<int, dim>, 5> cells_Ip = {Imm, Im, Ip, Ipp, Ippp};

    const array<CCTK_REAL, 2> rc_Im =
        eppm(gf_var, cells_Im, gf_is_press, gf_press, gf_vel_dir,
             reconstruct_params);
    const array<CCTK_REAL, 2> rc_Ip =
        eppm(gf_var, cells_Ip, gf_is_press, gf_press, gf_vel_dir,
             reconstruct_params);

    tmp = {rc_Im[1], rc_Ip[0]};
    break;
  }

  default:
    assert(0);
  }

  // In the following we assume that Container
  // can be constructed from an std::array!

  if constexpr (std::is_same<Container, std::array<CCTK_REAL,2>>::value) {
    return tmp; // zero-cost
  } else {
    return Container{tmp};
  }

}

} // namespace ReconX

#endif // RECONX_RECONSTRUCT_HXX
