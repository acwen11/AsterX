#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "../../../CarpetX/CarpetX/src/fillpatch.hxx"
#include "../../../CarpetX/CarpetX/src/schedule.hxx"
#include "../../../CarpetX/CarpetX/src/task_manager.hxx"

namespace AsterX {
using namespace CarpetX;

////////////////////////////////////////////////////////////////////////////////

extern "C" void AsterX_Sync(CCTK_ARGUMENTS) {
  // do nothing
}

void ApplyOuterBC(CCTK_ARGUMENTS, const std::vector<int> &groups) {
  task_manager tasks1;
  task_manager tasks2;

  for (const int gi : groups) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);

      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;

      if (leveldata.level == 0) {
        // Copy from adjacent boxes on same level and apply boundary conditions
        for (int tl = 0; tl < sync_tl; ++tl) {
          tasks1.submit_serially([&tasks2, &leveldata, &groupdata, tl]() {
            FillPatch_Sync(tasks2, groupdata, *groupdata.mfab.at(tl),
                           ghext->patchdata.at(leveldata.patch)
                               .amrcore->Geom(leveldata.level));
          });
        } // for tl
      }
    });
  } // for gi

  tasks1.run_tasks_serially();
  synchronize();
  tasks2.run_tasks_serially();
  synchronize();

  assert(ghext->num_patches() == 1);
}

extern "C" void AsterX_ApplyOuterBCOnPrim(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  static const std::vector<int> groups = {
      CCTK_GroupIndex("HydroBaseX::rho"),
      CCTK_GroupIndex("HydroBaseX::vel"),
      CCTK_GroupIndex("HydroBaseX::eps"),
      CCTK_GroupIndex("HydroBaseX::press"),
      CCTK_GroupIndex("HydroBaseX::Bvec"),
      CCTK_GroupIndex("HydroBaseX::temperature"),
      CCTK_GroupIndex("HydroBaseX::entropy"),
      CCTK_GroupIndex("HydroBaseX::Ye"),
      CCTK_GroupIndex("AsterX::zvec"),
      CCTK_GroupIndex("AsterX::svec"),
      CCTK_GroupIndex("AsterX::dBx_stag"),
      CCTK_GroupIndex("AsterX::dBy_stag"),
      CCTK_GroupIndex("AsterX::dBz_stag")};

  ApplyOuterBC(CCTK_PASS_CTOC, groups);
}

extern "C" void AsterX_RestrictAndApplyOuterBCOnFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  static const std::vector<int> restrict_groups = {
      CCTK_GroupIndex("AsterX::flux_x"), CCTK_GroupIndex("AsterX::flux_y"),
      CCTK_GroupIndex("AsterX::flux_z")};

  static const std::vector<int> bc_groups = {
      CCTK_GroupIndex("AsterX::vbar_xface"),
      CCTK_GroupIndex("AsterX::vbar_yface"),
      CCTK_GroupIndex("AsterX::vbar_zface"),
      CCTK_GroupIndex("AsterX::a_xface"),
      CCTK_GroupIndex("AsterX::a_yface"),
      CCTK_GroupIndex("AsterX::a_zface")};

  active_levels->loop_fine_to_coarse([&](const auto &leveldata) {
    if (leveldata.level < ghext->num_levels() - 1)
      RestrictNoPoison(cctkGH, leveldata.level, restrict_groups);
  });

  ApplyOuterBC(CCTK_PASS_CTOC, bc_groups);
}

extern "C" void AsterX_RestrictAuxTermsForAvecPsiRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;

  static const std::vector<int> restrict_groups = {
      CCTK_GroupIndex("AsterX::G"), CCTK_GroupIndex("AsterX::Ex"),
      CCTK_GroupIndex("AsterX::Ey"), CCTK_GroupIndex("AsterX::Ez")};

  static const std::vector<int> bc_groups = {CCTK_GroupIndex("AsterX::G")};

  active_levels->loop_fine_to_coarse([&](const auto &leveldata) {
    if (leveldata.level < ghext->num_levels() - 1)
      RestrictNoPoison(cctkGH, leveldata.level, restrict_groups);
  });

  ApplyOuterBC(CCTK_PASS_CTOC, bc_groups);
}

} // namespace AsterX
