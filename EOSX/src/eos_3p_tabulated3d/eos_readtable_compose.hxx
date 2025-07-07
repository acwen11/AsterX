#ifndef EOS_READTABLE_COMPOSE_HXX
#define EOS_READTABLE_COMPOSE_HXX

#define NTABLES 19

#include <cctk.h>
#include <string>
#include "../eos_3p.hxx"
#include "../utils/eos_linear_interp_ND.hxx"

namespace EOSX {
using namespace std;
using namespace eos_constants;

static linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> interptable;
static CCTK_INT ntemp, nrho, nye;
static CCTK_REAL energy_shift;

// Routine reading the EOS table and filling the corresponding object
CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
eos_readtable_compose(const string &filename) {
  CCTK_VINFO("Reading EOS table '%s'", filename.c_str());

  auto fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  assert(fapl_id >= 0);
  hid_t file_id = 0;
  int rank_id;
  CHECK_ERROR(MPI_Comm_rank(MPI_COMM_WORLD, &rank_id));

#ifdef H5_HAVE_PARALLEL
  CHECK_ERROR(H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL));
  CHECK_ERROR(H5Pset_all_coll_metadata_ops(fapl_id, true));
  file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl_id);
#else
  fapl_id = H5P_DEFAULT;
  if (rank_id == 0) {
    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl_id);
  }
// TODO: bcast table info to all ranks
#endif

  assert(file_id >= 0);

  // Get number of points
  get_hdf5_int_dset(file_id, "pointstemp", 1, &ntemp);
  get_hdf5_int_dset(file_id, "pointsrho", 1, &nrho);
  get_hdf5_int_dset(file_id, "pointsye", 1, &nye);

  const int npoints = ntemp * nrho * nye;

  CCTK_VINFO("EOS table dimensions: ntemp = %d, nrho = %d, nye = %d", ntemp,
             nrho, nye);

  // Allocate memory for tables

  CCTK_REAL *logrho, *logtemp, *yes;
  CCTK_REAL *epstable;
  CCTK_REAL *alltables;

  CCTK_REAL *alltables_temp;
  if (!(alltables_temp = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
            npoints * NTABLES * sizeof(CCTK_REAL)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }
  if (!(logrho = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
            nrho * sizeof(CCTK_REAL)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }
  if (!(logtemp = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
            ntemp * sizeof(CCTK_REAL)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }
  if (!(yes = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
            nye * sizeof(CCTK_REAL)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims[2] = {NTABLES, (hsize_t)npoints};
  hsize_t var3[2] = {1, (hsize_t)npoints};
  hid_t mem3 = H5Screate_simple(2, table_dims, NULL);

  // Read additional tables and variables
  get_hdf5_real_dset(file_id, "nb", nrho, logrho);
  get_hdf5_real_dset(file_id, "t", ntemp, logtemp);
  get_hdf5_real_dset(file_id, "yq", nye, yes);

  // Thermo Table
  // Number of variables in the thermo table

  hid_t thermo_id;
  HDF5_ERROR(thermo_id = H5Gopen(file, "/Thermo_qty"));
  int nthermo;
  get_hdf5_int_dset(thermo_id, "pointsqty", 1, &nthermo);

  // Read thermo index array
  int *thermo_index = new int[nthermo];
  get_hdf5_real_dset(thermo_id, "index_thermo", thermo_index, H5T_NATIVE_INT,
                     H5S_ALL);

  // Allocate memory and read table
  double *thermo_table = new double[nthermo * nrho * ntemp * nye];
  READ_EOS_HDF5_COMPOSE(thermo_id, "thermo", thermo_table, H5T_NATIVE_DOUBLE,
                        H5S_ALL);

  /////////////////

  // hydro (and munu)
  get_hdf5_real_dset(file_id, "logpress", npoints,
                     &alltables_temp[0 * npoints]);
  get_hdf5_real_dset(file_id, "logenergy", npoints,
                     &alltables_temp[1 * npoints]);
  get_hdf5_real_dset(file_id, "entropy", npoints, &alltables_temp[2 * npoints]);
  get_hdf5_real_dset(file_id, "munu", npoints, &alltables_temp[3 * npoints]);
  get_hdf5_real_dset(file_id, "cs2", npoints, &alltables_temp[4 * npoints]);
  get_hdf5_real_dset(file_id, "dedt", npoints, &alltables_temp[5 * npoints]);
  get_hdf5_real_dset(file_id, "dpdrhoe", npoints, &alltables_temp[6 * npoints]);
  get_hdf5_real_dset(file_id, "dpderho", npoints, &alltables_temp[7 * npoints]);

  // chemical potentials
  get_hdf5_real_dset(file_id, "muhat", npoints, &alltables_temp[8 * npoints]);
  get_hdf5_real_dset(file_id, "mu_e", npoints, &alltables_temp[9 * npoints]);
  get_hdf5_real_dset(file_id, "mu_p", npoints, &alltables_temp[10 * npoints]);
  get_hdf5_real_dset(file_id, "mu_n", npoints, &alltables_temp[11 * npoints]);

  // compositions
  get_hdf5_real_dset(file_id, "Xa", npoints, &alltables_temp[12 * npoints]);
  get_hdf5_real_dset(file_id, "Xh", npoints, &alltables_temp[13 * npoints]);
  get_hdf5_real_dset(file_id, "Xn", npoints, &alltables_temp[14 * npoints]);
  get_hdf5_real_dset(file_id, "Xp", npoints, &alltables_temp[15 * npoints]);

  // average nucleus
  get_hdf5_real_dset(file_id, "Abar", npoints, &alltables_temp[16 * npoints]);
  get_hdf5_real_dset(file_id, "Zbar", npoints, &alltables_temp[17 * npoints]);

  // Gamma
  get_hdf5_real_dset(file_id, "gamma", npoints, &alltables_temp[18 * npoints]);

  // Read additional tables and variables
  get_hdf5_real_dset(file_id, "logrho", nrho, logrho);
  get_hdf5_real_dset(file_id, "logtemp", ntemp, logtemp);
  get_hdf5_real_dset(file_id, "ye", nye, yes);
  get_hdf5_real_dset(file_id, "energy_shift", 1, &energy_shift);

  CHECK_ERROR(H5Pclose(fapl_id));
  CHECK_ERROR(H5Sclose(mem3));

#ifdef H5_HAVE_PARALLEL
  CHECK_ERROR(H5Fclose(file_id));
#else
  if (rank_id == 0) {
    CHECK_ERROR(H5Fclose(file_id));
  }
#endif

  // Fill actual table
  if (!(alltables = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
            npoints * NTABLES * sizeof(CCTK_REAL)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }
  for (int iv = 0; iv < NTABLES; iv++)
    for (int k = 0; k < nye; k++)
      for (int j = 0; j < ntemp; j++)
        for (int i = 0; i < nrho; i++) {
          int indold = i + nrho * (j + ntemp * (k + nye * iv));
          int indnew = iv + NTABLES * (i + nrho * (j + ntemp * k));

          // Maybe swap temp axis?
          // int indnew = iv + NTABLES*(j + ntemp*(i + nrho*k));
          alltables[indnew] = alltables_temp[indold];
        }

  // free memory of temporary array
  amrex::The_Managed_Arena()->free(alltables_temp);

  // allocate epstable; a linear-scale eps table
  // that allows us to extrapolate to negative eps
  if (!(epstable = (CCTK_REAL *)amrex::The_Managed_Arena()->alloc(
            npoints * sizeof(CCTK_REAL)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Cannot allocate memory for EOS table");
  }

  // convert units, convert logs to natural log
  // The latter is great, because exp() is way faster than pow()
  // pressure
  energy_shift = energy_shift * EPSGF;
  for (int i = 0; i < nrho; i++) {
    // rewrite:
    // logrho[i] = log(pow(10.0,logrho[i]) * RHOGF);
    // by using log(a^b*c) = b*log(a)+log(c)
    logrho[i] = logrho[i] * log(10.) + log(RHOGF);
  }

  for (int i = 0; i < ntemp; i++) {
    // logtemp[i] = log(pow(10.0,logtemp[i]));
    logtemp[i] = logtemp[i] * log(10.0);
  }

  // convert units
  for (int i = 0; i < npoints; i++) {

    { // pressure
      int idx = 0 + NTABLES * i;
      alltables[idx] = alltables[idx] * log(10.0) + log(PRESSGF);
    }

    { // eps
      int idx = 1 + NTABLES * i;
      alltables[idx] = alltables[idx] * log(10.0) + log(EPSGF);
      epstable[i] = exp(alltables[idx]);
    }

    { // cs2
      int idx = 4 + NTABLES * i;
      alltables[idx] *= LENGTHGF * LENGTHGF / TIMEGF / TIMEGF;
    }

    { // dedT
      int idx = 5 + NTABLES * i;
      alltables[idx] *= EPSGF;
    }

    { // dpdrhoe
      int idx = 6 + NTABLES * i;
      alltables[idx] *= PRESSGF / RHOGF;
    }

    { // dpderho
      int idx = 7 + NTABLES * i;
      alltables[idx] *= PRESSGF / EPSGF;
    }
  }

  auto num_points =
      std::array<size_t, 3>{size_t(nrho), size_t(ntemp), size_t(nye)};

  auto logrho_ptr = std::unique_ptr<CCTK_REAL[]>(new CCTK_REAL[nrho]);
  auto logtemp_ptr = std::unique_ptr<CCTK_REAL[]>(new CCTK_REAL[ntemp]);
  auto ye_ptr = std::unique_ptr<CCTK_REAL[]>(new CCTK_REAL[nye]);
  auto alltables_ptr =
      std::unique_ptr<CCTK_REAL[]>(new CCTK_REAL[npoints * NTABLES]);

  for (int i = 0; i < nrho; ++i)
    logrho_ptr[i] = logrho[i];
  for (int i = 0; i < ntemp; ++i)
    logtemp_ptr[i] = logtemp[i];
  for (int i = 0; i < nye; ++i)
    ye_ptr[i] = yes[i];
  for (int i = 0; i < npoints * NTABLES; ++i)
    alltables_ptr[i] = alltables[i];

  amrex::The_Managed_Arena()->free(logrho);
  amrex::The_Managed_Arena()->free(logtemp);
  amrex::The_Managed_Arena()->free(yes);
  amrex::The_Managed_Arena()->free(alltables);
  amrex::The_Managed_Arena()->free(epstable);

  interptable = linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES>(
      std::move(alltables_ptr), std::move(num_points), std::move(logrho_ptr),
      std::move(logtemp_ptr), std::move(ye_ptr));

  // set up steps, mins, maxes here?
  return;
};
} // namespace EOSX
#endif
