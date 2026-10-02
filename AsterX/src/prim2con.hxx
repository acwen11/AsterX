#ifndef ASTERX_PRIM2CON_HXX
#define ASTERX_PRIM2CON_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "aster_utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;
using namespace AsterUtils;

struct prim {
  CCTK_REAL rho;
  vec<CCTK_REAL, 3> vel;
  CCTK_REAL eps, press, entropy;
  CCTK_REAL Ye;
  vec<CCTK_REAL, 3> Bvec;
};

struct cons {
  CCTK_REAL dens;
  vec<CCTK_REAL, 3> mom;
  CCTK_REAL tau, DEnt;
  CCTK_REAL DYe;
  vec<CCTK_REAL, 3> dBvec;
};

CCTK_DEVICE CCTK_HOST void prim2con(const smat<CCTK_REAL, 3> &g, const prim &pv,
                                    cons &cv);

} // namespace AsterX

#endif // ASTERX_PRIM2CON_HXX
