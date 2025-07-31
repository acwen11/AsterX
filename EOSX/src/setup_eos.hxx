#ifndef EOSX_SETUP_EOS_HXX
#define EOSX_SETUP_EOS_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>

#include "eos_1p.hxx"
#include "eos_1p_polytropic/eos_1p_polytropic.hxx"

#include "eos_3p.hxx"
#include "eos_3p_idealgas/eos_3p_idealgas.hxx"
#include "eos_3p_hybrid/eos_3p_hybrid.hxx"
#include "eos_3p_tabulated3d/eos_3p_tabulated3d.hxx"

namespace EOSX {

// initial data EOS
extern eos_1p_polytropic *global_eos_1p_poly;

// evolution EOS
extern eos_3p_idealgas *global_eos_3p_ig;
extern eos_3p_hybrid *global_eos_3p_hyb;
extern eos_3p_tabulated3d *global_eos_3p_tab3d;

} // namespace EOSX

#endif
