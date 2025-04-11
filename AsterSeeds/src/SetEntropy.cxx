#include <loop.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdio>
#include <cstdbool>
#include <cmath> 

#include "eos.hxx"
#include "eos_idealgas.hxx"

extern "C" void SetEntropy(CCTK_ARGUMENTS)
{

  using namespace EOSX;

  DECLARE_CCTK_ARGUMENTSX_SetEntropy;
  DECLARE_CCTK_PARAMETERS;

  if (set_entropy_postinitial) {

    eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
         rgye(ye_min, ye_max);
  
    const eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);
  
    grid.loop_all_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
  
            entropy(p.I) = eos_th.kappa_from_valid_rho_eps_ye(rho(p.I),eps(p.I),1.0);
  
    });
  }
}


