/* Adopted from Margherita EOS framework of FIL */
#ifndef EOS_CONSTANTS
#define EOS_CONSTANTS

namespace eos_constants {

constexpr double LENGTHGF = 6.77269222552442e-06;
constexpr double TIMEGF = 2.03040204956746e05;
constexpr double RHOGF = 1.61887093132742e-18;
constexpr double PRESSGF = 1.80123683248503e-39;
constexpr double EPSGF = 1.11265005605362e-21;
constexpr double INVRHOGF = 6.17714470405638e17;
constexpr double INVEPSGF = 8.98755178736818e20;
constexpr double INVPRESSGF = 5.55174079257738e38;
// for Temperature in polytropic EOS
constexpr double c2_cgs = 8.9875517873681764e+20;
constexpr double mnuc_MeV = 931.494061;      // Atomic mass unit
constexpr double mnuc_cgs = 1.660539040e-24; // Atomic mass unit cgs
constexpr double mnuc_Msun =
    mnuc_cgs * PRESSGF * (LENGTHGF * LENGTHGF * LENGTHGF) * c2_cgs;
constexpr double MeV_to_erg = 1.60217733e-6;
constexpr double cm3_to_fm3 = 1.0e39;
constexpr double avogadro = 6.0221367e23;
constexpr double m_neutron_MeV = 939.565379;

constexpr double hbarc_MeV_fm = 197.326978812;
} // namespace eos_constants
#endif
