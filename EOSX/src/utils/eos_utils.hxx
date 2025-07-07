#ifndef EOS_UTILS_HXX
#define EOS_UTILS_HXX

#include <cmath>
#include <array>
#include <string>
#include <unordered_map>
#include <algorithm>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <hdf5.h>
#include <mpi.h>

namespace EOSX {

using namespace std;
using namespace Loop;
using namespace Arith;

// Macro checking for errors coming from routines returning error codes
#define CHECK_ERROR(routine)                                                   \
  do {                                                                         \
    const auto err = routine;                                                  \
    if (err < 0) {                                                             \
      CCTK_VERROR("Routine '%s' returned error code %d", #routine, err);       \
    }                                                                          \
  } while (0)

/// Class representing a range
struct eos_range {
  CCTK_REAL min; ///< Minimum
  CCTK_REAL max; ///< Maximum
  /// Default constructor: empty range.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_range()
      : min(0), max(0) {}
  /// Construct from minimum and maximum
  CCTK_DEVICE
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_range(CCTK_REAL min_,
                                                          CCTK_REAL max_)
      : min(min_), max(max_) {}
  /// Check if value is contained in [min,max]. False for NAN or INF.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
  contains(CCTK_REAL x) const {
    return (x >= min) && (x <= max);
  }
};

// TODO: add enums for error messages
/// Class representing error conditions in EOS calls.
struct eos_status {
  bool failed; ///< Set to true if parameters are out of range/NAN/INF
  //  std::string  err_msg; ///< Error description in case of failure, else
  //  undefined.
  /// Default constructor: Set to no failure.
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline eos_status()
      : failed(false) {}
  /// Set fail flag and error message.
  //  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  //  fail(std::string msg) {
  //    failed = true;
  //    err_msg = msg;
  //  }
};

// A map linking HDF5 actual I/O modes to appropriate descriptive strings
#ifdef H5_HAVE_PARALLEL
const std::unordered_map<H5D_mpio_actual_io_mode_t, std::string>
    H5D_mpio_actual_io_mode_map{
        {H5D_MPIO_NO_COLLECTIVE, "H5D_MPIO_NO_COLLECTIVE"},
        {H5D_MPIO_CHUNK_INDEPENDENT, "H5D_MPIO_CHUNK_INDEPENDENT"},
        {H5D_MPIO_CHUNK_COLLECTIVE, "H5D_MPIO_CHUNK_COLLECTIVE"},
        {H5D_MPIO_CHUNK_MIXED, "H5D_MPIO_CHUNK_MIXED"},
        {H5D_MPIO_CONTIGUOUS_COLLECTIVE, "H5D_MPIO_CONTIGUOUS_COLLECTIVE"}};

/* A map linking HDF5 actual chunk optimization modes to appropriate descriptive
 * strings */
const std::unordered_map<H5D_mpio_actual_chunk_opt_mode_t, std::string>
    H5D_mpio_actual_chunk_opt_mode_map{
        {H5D_MPIO_NO_CHUNK_OPTIMIZATION, "H5D_MPIO_NO_CHUNK_OPTIMIZATION"},
        {H5D_MPIO_LINK_CHUNK, "H5D_MPIO_LINK_CHUNK"},
        {H5D_MPIO_MULTI_CHUNK, "H5D_MPIO_MULTI_CHUNK"}};

/* A map linking HDF5 non-collective I/O mode causes to appropriate descriptive
 * strings */
const std::unordered_map<uint32_t, std::string>
    H5Pget_mpio_no_collective_cause_map{
        {H5D_MPIO_COLLECTIVE, "H5D_MPIO_COLLECTIVE"},
        {H5D_MPIO_SET_INDEPENDENT, "H5D_MPIO_SET_INDEPENDENT"},
        {H5D_MPIO_DATATYPE_CONVERSION, "H5D_MPIO_DATATYPE_CONVERSION"},
        {H5D_MPIO_DATA_TRANSFORMS, "H5D_MPIO_DATA_TRANSFORMS"},
        {H5D_MPIO_MPI_OPT_TYPES_ENV_VAR_DISABLED,
         "H5D_MPIO_MPI_OPT_TYPES_ENV_VAR_DISABLED"},
        {H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES,
         "H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES"},
        {H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET,
         "H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET"},
        {H5D_MPIO_PARALLEL_FILTERED_WRITES_DISABLED,
         "H5D_MPIO_PARALLEL_FILTERED_WRITES_DISABLED"},
        {H5D_MPIO_ERROR_WHILE_CHECKING_COLLECTIVE_POSSIBLE,
         "H5D_MPIO_ERROR_WHILE_CHECKING_COLLECTIVE_POSSIBLE"},
        {H5D_MPIO_NO_COLLECTIVE_MAX_CAUSE, "H5D_MPIO_NO_COLLECTIVE_MAX_CAUSE"}};
#endif // H5_HAVE_PARALLEL

CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
get_hdf5_int_dset(const hid_t &file_id, const string &dset_name,
                  const int npoints, int *var) {

  const auto dset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
  assert(dset_id >= 0);

  const auto dtype_id = H5Dget_type(dset_id);
  assert(dtype_id >= 0);

  const auto dtypeclass = H5Tget_class(dtype_id);
  assert(dtypeclass == H5T_INTEGER);
  CHECK_ERROR(H5Tclose(dtype_id));

  auto dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  assert(dxpl_id >= 0);

#ifdef H5_HAVE_PARALLEL
  CHECK_ERROR(H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE));
#else
  dxpl_id = H5P_DEFAULT;
#endif

  CHECK_ERROR(H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, dxpl_id, var));
  CHECK_ERROR(H5Dclose(dset_id));

#ifdef H5_HAVE_PARALLEL
  H5D_mpio_actual_io_mode_t actual_io_mode;
  CHECK_ERROR(H5Pget_mpio_actual_io_mode(dxpl_id, &actual_io_mode));

  H5D_mpio_actual_chunk_opt_mode_t actual_chunk_opt_mode;
  CHECK_ERROR(
      H5Pget_mpio_actual_chunk_opt_mode(dxpl_id, &actual_chunk_opt_mode));

  uint32_t no_collective_cause_local, no_collective_cause_global;
  CHECK_ERROR(H5Pget_mpio_no_collective_cause(
      dxpl_id, &no_collective_cause_local, &no_collective_cause_global));

  if (actual_io_mode == H5D_MPIO_NO_COLLECTIVE or
      // actual_chunk_opt_mode      == H5D_MPIO_NO_CHUNK_OPTIMIZATION or  //
      // In general, input files are not chunked
      no_collective_cause_local != H5D_MPIO_COLLECTIVE or
      no_collective_cause_global != H5D_MPIO_COLLECTIVE) {
    CCTK_VWARN(
        1,
        "Actual I/O mode, chunk optimization and local and global "
        "non-collective I/O causes when reading data from dataset '%s': "
        "'%s', '%s', '%s', '%s'",
        H5D_mpio_actual_io_mode_map.at(actual_io_mode).c_str(),
        H5D_mpio_actual_chunk_opt_mode_map.at(actual_chunk_opt_mode).c_str(),
        H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_local)
            .c_str(),
        H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_global)
            .c_str(),
        dset_name.c_str());
  }
#endif // H5_HAVE_PARALLEL

  CHECK_ERROR(H5Pclose(dxpl_id));
}

// Routine reading an HDF5 real number dataset
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
get_hdf5_real_dset(const hid_t &file_id, const string &dset_name,
                   const int npoints, CCTK_REAL *var) {

  const auto dset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
  assert(dset_id >= 0);

  auto dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  assert(dxpl_id >= 0);

#ifdef H5_HAVE_PARALLEL
  CHECK_ERROR(H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE));
#else
  dxpl_id = H5P_DEFAULT;
#endif

  CHECK_ERROR(
      H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl_id, var));
  CHECK_ERROR(H5Dclose(dset_id));

#ifdef H5_HAVE_PARALLEL
  H5D_mpio_actual_io_mode_t actual_io_mode;
  CHECK_ERROR(H5Pget_mpio_actual_io_mode(dxpl_id, &actual_io_mode));

  H5D_mpio_actual_chunk_opt_mode_t actual_chunk_opt_mode;
  CHECK_ERROR(
      H5Pget_mpio_actual_chunk_opt_mode(dxpl_id, &actual_chunk_opt_mode));

  uint32_t no_collective_cause_local, no_collective_cause_global;
  CHECK_ERROR(H5Pget_mpio_no_collective_cause(
      dxpl_id, &no_collective_cause_local, &no_collective_cause_global));

  if (actual_io_mode == H5D_MPIO_NO_COLLECTIVE or
      // actual_chunk_opt_mode      == H5D_MPIO_NO_CHUNK_OPTIMIZATION or  //
      // In general, input files are not chunked
      no_collective_cause_local != H5D_MPIO_COLLECTIVE or
      no_collective_cause_global != H5D_MPIO_COLLECTIVE) {
    CCTK_VWARN(
        1,
        "Actual I/O mode, chunk optimization and local and global "
        "non-collective I/O causes when reading data from dataset '%s': "
        "'%s', '%s', '%s', '%s'",
        H5D_mpio_actual_io_mode_map.at(actual_io_mode).c_str(),
        H5D_mpio_actual_chunk_opt_mode_map.at(actual_chunk_opt_mode).c_str(),
        H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_local)
            .c_str(),
        H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_global)
            .c_str(),
        dset_name.c_str());
  }
#endif // H5_HAVE_PARALLEL

  CHECK_ERROR(H5Pclose(dxpl_id));
}

} // namespace EOSX

#endif // #ifndef UTILS_HXX
