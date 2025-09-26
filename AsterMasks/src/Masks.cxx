#include <loop.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C" void AsterMasks_One(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterMasks_One;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
  
          aster_mask_vc(p.I) = 1.0;

  });

  grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          aster_mask_cc(p.I) = 1.0;

  });
}

extern "C" void AsterMasks_OnePostStep(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTSX_AsterMasks_OnePostStep;
  DECLARE_CCTK_PARAMETERS;

  if (*actually_set == 1.0) {

    grid.loop_all_device<0, 0, 0>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    
            aster_mask_vc(p.I) = 1.0;

    });

    grid.loop_all_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

            aster_mask_cc(p.I) = 1.0;

    });

  }
}

extern "C" void AsterMasks_ActuallySetZero(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTSX_AsterMasks_ActuallySetZero;
  DECLARE_CCTK_PARAMETERS;

  *actually_set = 0.0;

}
