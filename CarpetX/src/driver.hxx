#ifndef DRIVER_HXX
#define DRIVER_HXX

#include "loop.hxx"
#include "valid.hxx"

#include <rational.hxx>

#include <cctk.h>

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <vector>

namespace CarpetX {
using namespace std;
using namespace Arith;

using Loop::dim;

using rat64 = rational<int64_t>;

static_assert(AMREX_SPACEDIM == dim,
              "AMReX's AMREX_SPACEDIM must be the same as Cactus's cctk_dim");

static_assert(is_same<amrex::Real, CCTK_REAL>::value,
              "AMReX's Real type must be the same as Cactus's CCTK_REAL");

class CactusAmrCore final : public amrex::AmrCore {
  int patch;

public:
  vector<bool> level_modified;

  CactusAmrCore();
  CactusAmrCore(int patch, const amrex::RealBox *rb, int max_level_in,
                const amrex::Vector<int> &n_cell_in, int coord = -1,
                amrex::Vector<amrex::IntVect> ref_ratios =
                    amrex::Vector<amrex::IntVect>(),
                const int *is_per = nullptr);
  CactusAmrCore(int patch, const amrex::RealBox &rb, int max_level_in,
                const amrex::Vector<int> &n_cell_in, int coord,
                amrex::Vector<amrex::IntVect> const &ref_ratios,
                amrex::Array<int, AMREX_SPACEDIM> const &is_per);
  CactusAmrCore(const amrex::AmrCore &rhs) = delete;
  CactusAmrCore &operator=(const amrex::AmrCore &rhs) = delete;

  virtual ~CactusAmrCore() override;

  virtual void ErrorEst(int level, amrex::TagBoxArray &tags, amrex::Real time,
                        int ngrow) override;
  void SetupLevel(int level, const amrex::BoxArray &ba,
                  const amrex::DistributionMapping &dm,
                  const function<string()> &why);
  virtual void
  MakeNewLevelFromScratch(int level, amrex::Real time,
                          const amrex::BoxArray &ba,
                          const amrex::DistributionMapping &dm) override;
  virtual void
  MakeNewLevelFromCoarse(int level, amrex::Real time, const amrex::BoxArray &ba,
                         const amrex::DistributionMapping &dm) override;
  virtual void RemakeLevel(int level, amrex::Real time,
                           const amrex::BoxArray &ba,
                           const amrex::DistributionMapping &dm) override;
  virtual void ClearLevel(int level) override;
};

// Cactus grid hierarchy extension
struct GHExt {

  struct CommonGroupData {
    int groupindex;
    int firstvarindex;
    int numvars;

    bool do_checkpoint; // whether to checkpoint
    bool do_restrict;   // whether to restrict

    vector<vector<why_valid_t> > valid; // [time level][var index]

    // TODO: add poison_invalid and check_valid functions

    friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                     const CommonGroupData &commongroupdata);
  };

  struct GlobalData {
    // all data that exists on all levels

    // For subcycling in time, there really should be one copy of each
    // integrated grid scalar per level. We don't do that yet; instead,
    // we assume that grid scalars only hold "analysis" data.

    struct ArrayGroupData : public CommonGroupData {
      vector<vector<CCTK_REAL> > data; // [time level][var index]
      int array_size;
      int dimension;
      int activetimelevels;
      int lsh[dim];
      int ash[dim];
      int gsh[dim];
      int lbnd[dim];
      int ubnd[dim];
      int bbox[2 * dim];
      int nghostzones[dim];

      ArrayGroupData() {
        dimension = -1;
        activetimelevels = -1;
        for (int d = 0; d < dim; d++) {
          lsh[d] = -1;
          ash[d] = -1;
          gsh[d] = -1;
          lbnd[d] = -1;
          ubnd[d] = -1;
          bbox[2 * d] = bbox[2 * d + 1] = -1;
          nghostzones[d] = -1;
        }
      }

      friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                       const ArrayGroupData &arraygroupdata);
    };
    // TODO: right now this is sized for the total number of groups
    vector<unique_ptr<ArrayGroupData> > arraygroupdata; // [group index]

    friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                     const GlobalData &globaldata);
  };
  GlobalData globaldata;

  struct PatchData {
    int patch;

    // AMReX grid structure
    // TODO: convert this from unique_ptr to optional
    unique_ptr<CactusAmrCore> amrcore;

    struct LevelData {
      int patch, level;
      const PatchData &patchdata() const;
      PatchData &patchdata();

      // This level uses subcycling with respect to the next coarser
      // level. (Ignored for the coarsest level.)
      bool is_subcycling_level;

      // Iteration and time at which this cycle level is valid
      rat64 iteration, delta_iteration;

      // Fabamrex::ArrayBase object holding a cell-centred BoxArray for
      // iterating over grid functions. This stores the grid structure
      // and its distribution over all processes, but holds no data.
      unique_ptr<amrex::FabArrayBase> fab;

      struct GroupData : public CommonGroupData {
        GroupData(int patch, int level, int gi);

        int patch, level;

        const LevelData &leveldata() const;
        LevelData &leveldata();

        array<int, dim> indextype;
        array<int, dim> nghostzones;
        vector<array<int, dim> > parities;

        // each amrex::MultiFab has numvars components
        vector<unique_ptr<amrex::MultiFab> > mfab; // [time level]

        // flux register between this and the next coarser level
        unique_ptr<amrex::FluxRegister> freg;
        // associated flux group indices
        array<int, dim> fluxes; // [dir]

        friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                         const GroupData &groupdata);
      };
      // TODO: right now this is sized for the total number of groups
      vector<unique_ptr<GroupData> > groupdata; // [group index]

      friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                       const LevelData &leveldata);
    };
    vector<LevelData> leveldata; // [reflevel]

    friend YAML::Emitter &operator<<(YAML::Emitter &yaml,
                                     const PatchData &patchdata);
  };
  vector<PatchData> patchdata; // [patch]

  int num_patches() const { return patchdata.size(); }
  int num_levels() const {
    int nlevels = 0;
    using std::max;
    for (const auto &pd : patchdata)
      nlevels = max(nlevels, int(pd.leveldata.size()));
    return nlevels;
  }

  friend YAML::Emitter &operator<<(YAML::Emitter &yaml, const GHExt &ghext);
  friend ostream &operator<<(ostream &os, const GHExt &ghext);
};

extern unique_ptr<GHExt> ghext;

inline const GHExt::PatchData &GHExt::PatchData::LevelData::patchdata() const {
  return ghext->patchdata.at(patch);
}
inline GHExt::PatchData &GHExt::PatchData::LevelData::patchdata() {
  return ghext->patchdata.at(patch);
}

inline const GHExt::PatchData::LevelData &
GHExt::PatchData::LevelData::GroupData::leveldata() const {
  return ghext->patchdata.at(patch).leveldata.at(level);
}
inline GHExt::PatchData::LevelData &
GHExt::PatchData::LevelData::GroupData::leveldata() {
  return ghext->patchdata.at(patch).leveldata.at(level);
}

amrex::Interpolater *get_interpolator(const array<int, dim> indextype);

typedef void apply_physbcs_t(const amrex::Box &, const amrex::FArrayBox &, int,
                             int, const amrex::Geometry &, CCTK_REAL,
                             const amrex::Vector<amrex::BCRec> &, int, int);
typedef amrex::PhysBCFunct<apply_physbcs_t *> CarpetXPhysBCFunct;
tuple<CarpetXPhysBCFunct, amrex::Vector<amrex::BCRec> >
get_boundaries(const GHExt::PatchData::LevelData::GroupData &groupdata);

} // namespace CarpetX

#endif // #ifndef DRIVER_HXX
