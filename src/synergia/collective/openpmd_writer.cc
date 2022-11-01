#include "openpmd_writer.h"

Space_charge_openPMD_writer::Space_charge_openPMD_writer(
    std::string const& file)
    : series(file, openPMD::Access::CREATE)
    , iteration(0)
    , write_interval(0) // every iteration
    , write(false)
{}

bool
Space_charge_openPMD_writer::start_iteration()
{
    ++iteration;

    if (write_interval < 0) {
        write = false;
    } else if (write_interval == 0) {
        write = true;
    } else if ((iteration % write_interval) == 1) {
        write = true;
    } else {
        write = false;
    }

    return write;
}

void
Space_charge_openPMD_writer::write_particles(Bunch const& bunch)
{
    if (!write) return;

    size_t npart = bunch.get_local_num();
    auto h_parts = bunch.get_host_particles();
    Reference_particle ref = bunch.get_reference_particle();
    double p_ref = ref.get_momentum() * (1.0 + ref.get_state()[Bunch::dpop]);

    auto parts_x = Kokkos::subview(h_parts, Kokkos::ALL, Bunch::x);
    auto parts_y = Kokkos::subview(h_parts, Kokkos::ALL, Bunch::y);
    auto parts_z = Kokkos::subview(h_parts, Kokkos::ALL, Bunch::cdt);

    auto moments_x = Kokkos::subview(h_parts, Kokkos::ALL, Bunch::xp);
    auto moments_y = Kokkos::subview(h_parts, Kokkos::ALL, Bunch::yp);
    auto moments_dpop = Kokkos::subview(h_parts, Kokkos::ALL, Bunch::dpop);

    Kokkos::View<double*, Kokkos::DefaultHostExecutionSpace> moments_z(
        "moments_z_openpmd", npart);

    Kokkos::parallel_for(
        "OpenPMD-writer-dpop-pz", npart, KOKKOS_LAMBDA(const int& i) {
            double dpop = moments_dpop(i);
            double p = p_ref * (dpop + 1.0);
            moments_z(i) =
                sqrt(p * p - parts_x(i) * parts_x(i) - parts_y(i) * parts_y(i));
        });

    openPMD::ParticleSpecies& protons =
        series.iterations[iteration].particles["bunch"];
    openPMD::Datatype datatype =
        openPMD::determineDatatype(openPMD::shareRaw(h_parts.data()));
    openPMD::Dataset dataset(datatype, {npart});

    protons["position"]["x"].resetDataset(dataset);
    protons["position"]["y"].resetDataset(dataset);
    protons["position"]["z"].resetDataset(dataset);

    protons["moments"]["x"].resetDataset(dataset);
    protons["moments"]["y"].resetDataset(dataset);
    protons["moments"]["z"].resetDataset(dataset);

    series.flush();

    protons["position"]["x"].storeChunk(
        openPMD::shareRaw(parts_x.data()), {0}, {npart});
    protons["position"]["y"].storeChunk(
        openPMD::shareRaw(parts_y.data()), {0}, {npart});
    protons["position"]["z"].storeChunk(
        openPMD::shareRaw(parts_z.data()), {0}, {npart});

    protons["moments"]["x"].storeChunk(
        openPMD::shareRaw(moments_x.data()), {0}, {npart});
    protons["moments"]["y"].storeChunk(
        openPMD::shareRaw(moments_y.data()), {0}, {npart});
    protons["moments"]["z"].storeChunk(
        openPMD::shareRaw(moments_z.data()), {0}, {npart});

    series.flush();
}
