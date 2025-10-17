#include <loop.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C" void AsterMasks_SetToOne(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterMasks_SetToOne;
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

extern "C" void AsterMasks_SetToOnePostStep(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTSX_AsterMasks_SetToOnePostStep;
  DECLARE_CCTK_PARAMETERS;

  if (*MaskBool == 1.0) {

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

extern "C" void AsterMasks_SetMaskBoolZero(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTSX_AsterMasks_SetMaskBoolZero;
  DECLARE_CCTK_PARAMETERS;

  *MaskBool = 0.0;

}
