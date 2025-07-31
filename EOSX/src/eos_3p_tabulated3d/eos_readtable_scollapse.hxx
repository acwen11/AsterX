#ifndef EOS_READTABLE_SCOLLAPSE_HXX
#define EOS_READTABLE_SCOLLAPSE_HXX

#define NTABLES 19

#include <cctk.h>
#include <string>
#include <mpi.h>
#include <hdf5.h>
#include "../eos_3p.hxx"
#include "../utils/eos_linear_interp_ND.hxx"

namespace EOSX {
using namespace std;
using namespace eos_constants;

CCTK_REAL *energy_shift;
linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> *interptable;

CCTK_HOST void eos_readtable_scollapse(const string &filename) {
  CCTK_VINFO("Reading Stellar Collapse EOS table '%s'", filename.c_str());
  CCTK_INT ntemp, nrho, nye;

  int rank_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_id);

  hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  assert(fapl_id >= 0);

#ifdef H5_HAVE_PARALLEL
  CHECK_ERROR(H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL));
  CHECK_ERROR(H5Pset_all_coll_metadata_ops(fapl_id, true));
#endif

  hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl_id);
  assert(file_id >= 0);

  if (rank_id == 0) {
    get_hdf5_int_dset(file_id, "pointstemp", 1, &ntemp);
    get_hdf5_int_dset(file_id, "pointsrho", 1, &nrho);
    get_hdf5_int_dset(file_id, "pointsye", 1, &nye);
  }

  MPI_Bcast(&ntemp, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nrho, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nye, 1, MPI_INT, 0, MPI_COMM_WORLD);

  CCTK_INT npoints = ntemp * nrho * nye;
  CCTK_VINFO("EOS dimensions: ntemp=%d, nrho=%d, nye=%d", ntemp, nrho, nye);

  energy_shift =
      (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(sizeof(CCTK_REAL));
  CCTK_REAL *logrho =
      (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(nrho * sizeof(CCTK_REAL));
  CCTK_REAL *logtemp =
      (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(ntemp * sizeof(CCTK_REAL));
  CCTK_REAL *yes =
      (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(nye * sizeof(CCTK_REAL));
  CCTK_REAL *alltables = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
      npoints * NTABLES * sizeof(CCTK_REAL));
  CCTK_REAL *epstable = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
      npoints * sizeof(CCTK_REAL));

  if (rank_id == 0) {
    const char *datasets[NTABLES] = {
        "logpress", "logenergy", "entropy", "munu", "cs2",  "dedt",
        "dpdrhoe",  "dpderho",   "muhat",   "mu_e", "mu_p", "mu_n",
        "Xa",       "Xh",        "Xn",      "Xp",   "Abar", "Zbar"};

    for (int i = 0; i < NTABLES; i++)
      get_hdf5_real_dset(file_id, datasets[i], npoints,
                         &alltables[i * npoints]);

    get_hdf5_real_dset(file_id, "logrho", nrho, logrho);
    get_hdf5_real_dset(file_id, "logtemp", ntemp, logtemp);
    get_hdf5_real_dset(file_id, "ye", nye, yes);
    get_hdf5_real_dset(file_id, "energy_shift", 1, energy_shift);
  }

  MPI_Bcast(alltables, npoints * NTABLES, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(logrho, nrho, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(logtemp, ntemp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(yes, nye, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(energy_shift, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  H5Fclose(file_id);
  H5Pclose(fapl_id);

  *energy_shift *= EPSGF;
  const CCTK_REAL ln10 = log(10.0);
  const CCTK_REAL inv_time2 = 1 / (TIMEGF * TIMEGF);

  for (int i = 0; i < nrho; i++)
    logrho[i] = logrho[i] * ln10 + log(RHOGF);

  for (int i = 0; i < ntemp; i++)
    logtemp[i] *= ln10;

  for (int i = 0; i < npoints; i++) {
    alltables[i] = alltables[i] * ln10 + log(PRESSGF); // PRESS
    alltables[npoints + i] = alltables[npoints + i] * ln10 + log(EPSGF); // EPS
    epstable[i] = exp(alltables[i]);                                     // EPS
    alltables[4 * npoints + i] *= LENGTHGF * LENGTHGF * inv_time2;       // CS2
    alltables[5 * npoints + i] *= EPSGF;                                 // DEDT
    alltables[6 * npoints + i] *= PRESSGF / RHOGF; // DPDRHOE
    alltables[7 * npoints + i] *= PRESSGF / EPSGF; // DPDERHO
  }

  interptable = new (amrex::The_Managed_Arena()->alloc(sizeof(*interptable)))
      linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES>(
          alltables, {(size_t)nrho, (size_t)ntemp, (size_t)nye}, logrho,
          logtemp, yes);
}
} // namespace EOSX
#endif
