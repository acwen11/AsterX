#ifndef EOSX_SETUP_EOS_HXX
#define EOSX_SETUP_EOS_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>

#include "eos_1p.hxx"
#include "eos_1p_polytropic/eos_1p_polytropic.hxx"
#include "eos_1p_polytropic/eos_1p_piecewise_polytropic.hxx"

#include "eos_3p.hxx"
#include "eos_3p_idealgas/eos_3p_idealgas.hxx"
#include "eos_3p_hybrid/eos_3p_hybrid.hxx"
#include "eos_3p_tabulated3d/eos_3p_tabulated3d.hxx"

namespace EOSX {

// initial data EOS
extern eos_1p_polytropic *global_eos_1p_poly;
extern eos_1p_piecewise_polytropic *global_eos_1p_pwpoly;

// evolution EOS
extern eos_3p_idealgas *global_eos_3p_ig;

using eos_3p_hybrid_poly   = eos_3p_hybrid<eos_1p_polytropic>;
using eos_3p_hybrid_pwpoly = eos_3p_hybrid<eos_1p_piecewise_polytropic>;

extern eos_3p_hybrid_poly   *global_eos_3p_hyb_poly;
extern eos_3p_hybrid_pwpoly *global_eos_3p_hyb_pwpoly;

extern eos_3p_tabulated3d *global_eos_3p_tab3d;

} // namespace EOSX

#endif
