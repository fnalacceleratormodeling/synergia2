#include "deposit.h"

#include "space_charge_3d_fd.h"
#include "space_charge_3d_fd_alias.h"
#include "space_charge_3d_fd_impl.h"
#include "space_charge_3d_fd_utils.h"

#include "space_charge_3d_kernels.h"

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/simple_timer.h"

namespace {
    double
    get_smallest_non_tiny(double val, double other1, double other2, double tiny)
    {
        double retval;
        if (val > tiny) {
            retval = val;
        } else {
            if ((other1 > tiny) && (other2 > tiny)) {
                retval = std::min(other1, other2);
            } else {
                retval = std::max(other1, other2);
            }
        }
        return retval;
    }
} // namespace

// Constructor
Space_charge_3d_fd::Space_charge_3d_fd(Space_charge_3d_fd_options const& ops)
    : Collective_operator("sc_3d_fd", 1.0)
    , options(ops)
    , bunch_sim_id()
    , domain(ops.shape, {1.0, 1.0, 1.0})
    , use_fixed_domain(false)
    , allocated(false)
{}

// Destructor
Space_charge_3d_fd::~Space_charge_3d_fd()
{
    if (allocated) {
        destroy_sc3d_fd();
        allocated = false;
    }
}

void
Space_charge_3d_fd::set_fixed_domain(std::array<double, 3> offset,
                                     std::array<double, 3> size)
{

    /* todo check for physics here! */
    domain = Rectangular_grid_domain(options.shape, size, offset, false);

    use_fixed_domain = true;

    gctx.Lx = static_cast<PetscReal>(domain.get_physical_size()[0]);
    gctx.Ly = static_cast<PetscReal>(domain.get_physical_size()[1]);
    gctx.Lz = static_cast<PetscReal>(domain.get_physical_size()[2]);
}

void
Space_charge_3d_fd::get_local_charge_density(Bunch const& bunch)
{
    scoped_simple_timer timer("sc3d_fd_local_rho");

    deposit_charge_rectangular_3d_kokkos_scatter_view(
        lctx.seqrho_view, domain, domain.get_grid_shape(), bunch);
}

void
Space_charge_3d_fd::apply_impl(Bunch_simulator& sim,
                               double time_step,
                               Logger& logger)
{

    PetscErrorCode ierr;

    logger << "    Space charge 3d finite difference\n";
    scoped_simple_timer timer("sc3d_fd_total");

    // count number of bunches
    int num_bunches_in_bunch_sim = 0;
    int num_bunches_in_train_0 = 0;
    for (size_t t = 0; t < 2; ++t) {
        for (size_t b = 0; b < sim[t].get_bunch_array_size(); ++b) {
            num_bunches_in_bunch_sim += 1;
            if (t == 0) num_bunches_in_train_0 += 1;
        }
    }
    if (num_bunches_in_bunch_sim != num_bunches_in_train_0) {
        throw std::runtime_error(
            "sc3d-fd only works on a single bunch train, work is ongoing \
		    to make it work on multiple bunch trains!");
    }

    // construct the workspace for a new bunch simulator
    if (bunch_sim_id != sim.id()) {
        bunch_sim_id = sim.id();
        // Assumption: When an MPI rank has more than 1 bunch (within the same
        // train), all the bunches on the rank have the same communicator and
        // the same distribution. Reason: when constructing a bunch train, one
        // can have one the two scenarios: [1] number of bunches > number of MPI
        // ranks, each rank always has only 1 bunch [2] number of bunches <
        // number of MPI ranks, we require number of bunches to be divisible by
        // number of  MPI ranks and each bunch has a MPI communicator of size 1.
        // Update this bit when enabling bunch sim with two bunch trains!
        allocate_sc3d_fd(sim[0][0]);
        allocated = true;

        /* Functionality that is present in update_domain that must be called
          where a static domain is used ! */
        if (use_fixed_domain) {

            auto left_x = static_cast<PetscReal>(domain.get_left()[0]);
            auto left_y = static_cast<PetscReal>(domain.get_left()[1]);
            auto left_z = static_cast<PetscReal>(domain.get_left()[2]);

            /* Public API that cannot return PetscErrorCode, abort on failure */
            PetscCallAbort(gctx.bunch_comm,
                           DMDASetUniformCoordinates(sctx.da,
                                                     left_x,
                                                     left_x + gctx.Lx,
                                                     left_y,
                                                     left_y + gctx.Ly,
                                                     left_z,
                                                     left_z + gctx.Lz));
            PetscCallAbort(gctx.bunch_comm, compute_mat(sctx, gctx));
        }
    }

    // apply to bunches
    for (size_t t = 0; t < 2; ++t) {
        for (size_t b = 0; b < sim[t].get_bunch_array_size(); ++b) {
            // Using PetscCallAbort instead of PetscCall as this fucntion cannot
            // return a PetscErrorCode. Since we can't do much if the following
            // fails, we just abort instead of try/catch/recover.
            PetscCallAbort(gctx.bunch_comm,
                           apply_bunch(sim[t][b], time_step, logger));
        }
    }
}

PetscErrorCode
Space_charge_3d_fd::apply_bunch(Bunch& bunch, double time_step, Logger& logger)
{
    PetscFunctionBeginUser;

    // update domain only when not using fixed
    if (!use_fixed_domain) update_domain(bunch);

    // charge density
    get_local_charge_density(bunch); // [C/m^3]

    // DEBUGGING!
    if (gctx.dumps) {
        PetscViewer hdf5_viewer;
        PetscCall(PetscPrintf(gctx.bunch_comm,
                              "Dumping rho vector on all ranks!\n "));
        std::string filename = "rho_on_rank_";
        filename.append(std::to_string(gctx.global_rank));
        filename.append(".h5");
        PetscCall(PetscViewerHDF5Open(
            PETSC_COMM_SELF, filename.c_str(), FILE_MODE_WRITE, &hdf5_viewer));
        PetscCall(VecView(lctx.seqrho, hdf5_viewer));
        PetscCall(PetscViewerDestroy(&hdf5_viewer));
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       global to subcomm scatters
       - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    {
        scoped_simple_timer timer("sc3d_fd_rho_local_to_subcomms");

        /* Begin global (alias of local) to subcomm scatters! */
        for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
            PetscCall(VecScatterBegin(gctx.scat_glocal_to_subcomms[i],
                                      gctx.rho_global_local,
                                      gctx.rho_global_subcomm[i],
                                      ADD_VALUES,
                                      SCATTER_FORWARD));
        }

        /* Hopefully there is some unrelated work that can occur here! */

        /* End global (alias of local) to subcomm scatters! */
        for (PetscInt i = 0; i < gctx.nsubcomms; i++) {
            PetscCall(VecScatterEnd(gctx.scat_glocal_to_subcomms[i],
                                    gctx.rho_global_local,
                                    gctx.rho_global_subcomm[i],
                                    ADD_VALUES,
                                    SCATTER_FORWARD));
        }
    }

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         concurrent operations on each solver subcommunicator
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    // DEBUGGING!
    if (gctx.dumps) {
        PetscViewer hdf5_viewer;
        PetscCall(PetscPrintf(gctx.bunch_comm,
                              "Dumping rho vector on all subcomms!\n"));
        std::string filename = "rho_on_subcomm";
        filename.append(std::to_string(sctx.solversubcommid));
        filename.append(".h5");
        PetscCall(PetscViewerHDF5Open(sctx.solversubcomm,
                                      filename.c_str(),
                                      FILE_MODE_WRITE,
                                      &hdf5_viewer));
        PetscCall(VecView(sctx.rho_subcomm, hdf5_viewer));
        PetscCall(PetscViewerDestroy(&hdf5_viewer));
    }

    // DEBUGGING!
    if (gctx.dumps) {
        PetscViewer ascii_viewer;
        PetscCall(
            PetscPrintf(gctx.bunch_comm, "Dumping matrix on all subcomms!\n"));
        std::string filename = "mat_on_subcomm";
        filename.append(std::to_string(sctx.solversubcommid));
        filename.append(".ascii");
        PetscCall(PetscObjectSetName((PetscObject)(sctx.A), "A_on_sctx"));
        PetscCall(PetscViewerASCIIOpen(
            sctx.solversubcomm, filename.c_str(), &ascii_viewer));
        PetscCall(MatView(sctx.A, ascii_viewer));
        PetscCall(PetscViewerDestroy(&ascii_viewer));
    }

    {
        scoped_simple_timer timer("sc3d_fd_linear_solver");

        /* Solve for phi on each subcomm! */
        PetscCall(solve(sctx, gctx));
    }

    // DEBUGGING!
    if (gctx.dumps) {
        PetscViewer hdf5_viewer;
        PetscCall(PetscPrintf(gctx.bunch_comm,
                              "Dumping phi vector on all subcomms!\n"));
        std::string filename = "phi_on_subcomm_";
        filename.append(std::to_string(sctx.solversubcommid));
        filename.append(".h5");
        PetscCall(PetscViewerHDF5Open(sctx.solversubcomm,
                                      filename.c_str(),
                                      FILE_MODE_WRITE,
                                      &hdf5_viewer));
        PetscCall(VecView(sctx.phi_subcomm, hdf5_viewer));
        PetscCall(PetscViewerDestroy(&hdf5_viewer));
    }

    {
        scoped_simple_timer timer("sc3d_fd_phi_subcomms_to_local");

        /* Begin subcomm (alias of local) to local scatters! */
        PetscCall(VecScatterBegin(sctx.scat_subcomm_to_local,
                                  sctx.phi_subcomm,
                                  sctx.phi_subcomm_local,
                                  INSERT_VALUES,
                                  SCATTER_REVERSE));

        /* Hopefully there is some unrelated work that can occur here! */

        /* End subcomm (alias of local) to local scatters! */
        PetscCall(VecScatterEnd(sctx.scat_subcomm_to_local,
                                sctx.phi_subcomm,
                                sctx.phi_subcomm_local,
                                INSERT_VALUES,
                                SCATTER_REVERSE));
    }

    // DEBUGGING!
    if (gctx.dumps) {
        PetscViewer hdf5_viewer;
        PetscCall(
            PetscPrintf(gctx.bunch_comm, "Dumping phi vector on all ranks!\n"));
        std::string filename = "phi_on_rank_";
        filename.append(std::to_string(gctx.global_rank));
        filename.append(".h5");
        PetscCall(PetscViewerHDF5Open(
            PETSC_COMM_SELF, filename.c_str(), FILE_MODE_WRITE, &hdf5_viewer));
        PetscCall(VecView(lctx.seqphi, hdf5_viewer));
        PetscCall(PetscViewerDestroy(&hdf5_viewer));
    }

    get_force();

    // DEBUGGING!
    if (gctx.dumps) {
        /*
          Disable this temporarily until a commxx of type MPI_COMM_SELF
          can be created

    PetscCall(PetscPrintf(gctx.bunch_comm,
                        "Dumping enx/eny/enz vector on all ranks!\n"));

    std::string filename;

    filename = "enx_on_rank_";
    filename.append(std::to_string(gctx.global_rank));
    filename.append(".h5");
    Hdf5_file file_x(filename, Hdf5_file::Flag::read_write, );
    file_x.write("enx", lctx.enx.data(), lctx.enx.size(), true);

    filename = "eny_on_rank_";
    filename.append(std::to_string(gctx.global_rank));
    filename.append(".h5");
    Hdf5_file file_y(filename, Hdf5_file::Flag::read_write, Commxx());
    file_y.write("eny", lctx.eny.data(), lctx.eny.size(), true);

    filename = "enz_on_rank_";
    filename.append(std::to_string(gctx.global_rank));
    filename.append(".h5");
    Hdf5_file file_z(filename, Hdf5_file::Flag::read_write, Commxx());
    file_z.write("enz", lctx.enz.data(), lctx.enz.size(), true);
    */
    }

    apply_kick(bunch, time_step);

    PetscFunctionReturn(0);
}

// get force
void
Space_charge_3d_fd::get_force()
{
    scoped_simple_timer timer("sc3d_fd_force");

    sc3d_kernels::zyx::alg_force_extractor alg(lctx.seqphi_view,
                                               lctx.enx,
                                               lctx.eny,
                                               lctx.enz,
                                               domain.get_grid_shape(),
                                               domain.get_cell_size());
    Kokkos::parallel_for(gctx.nsize, alg);
    Kokkos::fence();
}

// apply kick
void
Space_charge_3d_fd::apply_kick(Bunch& bunch, double time_step)
{
    scoped_simple_timer timer("sc3d_fd_kick");

    auto ref = bunch.get_reference_particle();

    double q = bunch.get_particle_charge() * pconstants::e;
    double m = bunch.get_mass();

    double gamma = ref.get_gamma();
    double beta = ref.get_beta();
    double pref = ref.get_momentum();

    auto g = domain.get_grid_shape();
    auto h = domain.get_cell_size();
    auto l = domain.get_left();

    double fn_norm = (1.0 / (pconstants::epsilon0));

    //                   (1.0 / (4.0 * mconstants::pi * pconstants::epsilon0));

    double unit_conversion = pconstants::c / (1e9 * pconstants::e);
    double factor = options.kick_scale * unit_conversion * q * time_step *
                    fn_norm / (pref * gamma * gamma * beta);

    auto parts = bunch.get_local_particles();
    auto masks = bunch.get_local_particle_masks();

    sc3d_kernels::zyx::alg_kicker kicker(
        parts, masks, lctx.enx, lctx.eny, lctx.enz, g, h, l, factor, pref, m);

    Kokkos::parallel_for(bunch.size(), kicker);
    Kokkos::fence();
}

// update_domain
PetscErrorCode
Space_charge_3d_fd::update_domain(Bunch const& bunch)
{

    PetscFunctionBeginUser;

    scoped_simple_timer timer("sc3d_fd_domain");

    auto spatial_mean_stddev =
        Core_diagnostics::calculate_spatial_mean_stddev(bunch);
    auto mean_x = spatial_mean_stddev(0);
    auto mean_y = spatial_mean_stddev(1);
    auto mean_z = spatial_mean_stddev(2);
    auto stddev_x = spatial_mean_stddev(3);
    auto stddev_y = spatial_mean_stddev(4);
    auto stddev_z = spatial_mean_stddev(5);

    const double tiny = 1.0e-10;

    if ((stddev_x < tiny) && (stddev_y < tiny) && (stddev_z < tiny)) {
        throw std::runtime_error(
            "Space_charge_3d_open_hockney_eigen::update_domain: "
            "all three spatial dimensions have neglible extent");
    }

    std::array<double, 3> offset{mean_x, mean_y, mean_z};

    std::array<double, 3> size{
        options.n_sigma *
            get_smallest_non_tiny(stddev_x, stddev_y, stddev_z, tiny),
        options.n_sigma *
            get_smallest_non_tiny(stddev_y, stddev_x, stddev_z, tiny),
        options.n_sigma *
            get_smallest_non_tiny(stddev_z, stddev_x, stddev_y, tiny)};

    domain = Rectangular_grid_domain(options.shape, size, offset, false);

    gctx.Lx = static_cast<PetscReal>(domain.get_physical_size()[0]);
    gctx.Ly = static_cast<PetscReal>(domain.get_physical_size()[1]);
    gctx.Lz = static_cast<PetscReal>(domain.get_physical_size()[2]);

    auto left_x = static_cast<PetscReal>(domain.get_left()[0]);
    auto left_y = static_cast<PetscReal>(domain.get_left()[1]);
    auto left_z = static_cast<PetscReal>(domain.get_left()[2]);

    PetscCall(DMDASetUniformCoordinates(sctx.da,
                                        left_x,
                                        left_x + gctx.Lx,
                                        left_y,
                                        left_y + gctx.Ly,
                                        left_z,
                                        left_z + gctx.Lz));

    PetscCall(compute_mat(sctx, gctx));

    PetscFunctionReturn(0);
}

PetscErrorCode
Space_charge_3d_fd::allocate_sc3d_fd(const Bunch& bunch)
{
    PetscFunctionBeginUser;

    scoped_simple_timer timer("sc3d_fd_allocations");

    /* size of seqphi/seqrho vectors/views is size of domain! */
    gctx.nsize = options.shape[0] * options.shape[1] * options.shape[2];
    gctx.nsize_x = options.shape[0];
    gctx.nsize_y = options.shape[1];
    gctx.nsize_z = options.shape[2];

    MPI_Comm bunch_comm = MPI_Comm(bunch.get_comm());

    /* store MPI communicator of bunch in gctx */
    PetscCall(PetscCommDuplicate(bunch_comm, &gctx.bunch_comm, NULL));
    PetscCallMPI(MPI_Comm_rank(gctx.bunch_comm, &gctx.global_rank));
    PetscCallMPI(MPI_Comm_size(gctx.bunch_comm, &gctx.global_size));

    /* Initialize task subcomms, display task-subcomm details */
    PetscCall(init_solver_subcomms(sctx, gctx));

    /* Local rho and phi vectors on each MPI rank */
    PetscCall(init_local_vecs(lctx, gctx));

    /* rho and phi vectors on each subcomm */
    PetscCall(init_subcomm_vecs(sctx, gctx));

    /* create global aliases of local vectors */
    PetscCall(init_global_local_aliases(lctx, gctx));

    /* create global aliases of subcomm vectors */
    PetscCall(init_global_subcomm_aliases(sctx, gctx));

    /* create subcomm aliases of local vectors */
    PetscCall(init_subcomm_local_aliases(lctx, sctx, gctx));

    /* create DM and Matrix on subcomms */
    PetscCall(init_subcomm_mat(sctx, gctx));

    /* Initialize global (alias of local) to subcomm scatters */
    PetscCall(init_global_subcomm_scatters(sctx, gctx));

    /* Initialize subcomm (alias of local) to local scatters */
    PetscCall(init_subcomm_local_scatters(lctx, sctx, gctx));

    PetscFunctionReturn(0);
}

PetscErrorCode
Space_charge_3d_fd::destroy_sc3d_fd()
{
    PetscFunctionBeginUser;

    PetscCall(finalize(lctx, sctx, gctx));

    PetscFunctionReturn(0);
}
