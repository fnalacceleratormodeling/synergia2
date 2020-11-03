
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/core_diagnostics.h"

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <sstream>


#if 0
void
Bunch::assign_spectator_ids(int local_offset)
{

#if 0
    int global_offset, request_num;
    if (comm_sptr->get_rank() == 0) {
        request_num = total_s_num;
    } else {
        request_num = 0;
    }
    global_offset = particle_id_offset.get(request_num, *comm_sptr);
#endif

    int global_offset = 0;
    for (int i = 0; i < local_s_num; ++i) {
        (*local_s_particles)[i][id] = i + local_offset + global_offset;
    }
}
#endif


#if 0
std::string
Bunch::get_local_particles_serialization_path() const
{
    std::stringstream sstream;
    sstream << "local_particles_";
    sstream << bucket_index;
    sstream << ".h5";
    return get_serialization_path(sstream.str());
}
#endif


template<>
Bunch::bunch_t(
        Reference_particle const& reference_particle, 
        int total_num, 
        double real_num, 
        Commxx bunch_comm,
        int total_spectator_num,
        int bunch_index,
        int bucket_index,
        int array_index )
    : comm(std::make_shared<Commxx>(bunch_comm))
    , boundary(LB::open)
    , boundary_param(0.0)
    , ref_part(reference_particle)
    , design_ref_part(reference_particle)
    , particle_charge(reference_particle.get_charge())
    , real_num(real_num)
    , parts{ BunchParticles(PG::regular, total_num, -1, *comm),
             BunchParticles(PG::spectator, total_spectator_num, -1, *comm) }
    , bunch_index(bunch_index)
    , bucket_index(bucket_index)
    , array_index(array_index)
{
}

template<>
Bunch::bunch_t()
    : comm(new Commxx())
    , boundary(LB::open)
    , boundary_param(0.0)
    , ref_part()
    , design_ref_part()
    , particle_charge(ref_part.get_charge())
    , real_num(1.0)
    , parts{ BunchParticles(PG::regular, 0, 0, *comm), 
             BunchParticles(PG::spectator, 0, 0, *comm) }
    , bunch_index(0)
    , bucket_index(0)
    , array_index(0)
{
}

template<>
void
Bunch::inject(Bunch const& o)
{
    const double weight_tolerance = 1.0e-14;
    const double particle_tolerance = 1.0e-14;

    auto const& ref = ref_part;
    auto const& oref = o.ref_part;

    // The charge and mass of the bunch particles must match
    if (particle_charge != o.particle_charge) 
    {
        throw std::runtime_error(
                "Bunch.inject: bunch particle charges do not match.");
    }

    if (std::abs(ref.get_mass()/oref.get_mass() - 1.0) > particle_tolerance) 
    {
        throw std::runtime_error(
                "Bunch:inject: bunch particle masses do not match.");
    }

    // total num for regualr particles
    int total_num = parts[0].num_total();

    // can only check particle weight if total_num is nonzero
    if (total_num != 0) 
    {
        double wgt1 = real_num/total_num;
        double wgt2 = o.get_real_num()/o.get_total_num();

        if (std::abs(wgt1/wgt2 - 1.0) > weight_tolerance) 
        {
            throw std::runtime_error(
                "Bunch.inject: macroparticle weight of injected bunch does not match.");
        }
    }

    // prepare
    double target_momentum = ref.get_momentum();
    double injected_momentum = oref.get_momentum();
    double pdiff = injected_momentum / target_momentum;

    karray1d_dev ref_st_diff("", 6);
    karray1d_dev tgt_st("", 6);
    karray1d_dev inj_st("", 6);

    auto h_ref_st_diff = Kokkos::create_mirror_view(ref_st_diff);
    auto h_tgt_st      = Kokkos::create_mirror_view(tgt_st);
    auto h_inj_st      = Kokkos::create_mirror_view(inj_st);

    for (int i = 0; i < 6; ++i) 
    {
        h_tgt_st(i) = ref.get_state()[i];
        h_inj_st(i) = oref.get_state()[i];
        h_ref_st_diff(i) = h_inj_st(i) - h_tgt_st(i);
    }

    Kokkos::deep_copy(ref_st_diff, h_ref_st_diff);
    Kokkos::deep_copy(tgt_st, h_tgt_st);
    Kokkos::deep_copy(inj_st, h_inj_st);

    // regular
    get_bunch_particles(PG::regular).inject(
            o.get_bunch_particles(PG::regular),
            ref_st_diff, tgt_st, inj_st, pdiff );

    // spectator
    get_bunch_particles(PG::spectator).inject(
            o.get_bunch_particles(PG::spectator),
            ref_st_diff, tgt_st, inj_st, pdiff );

    // update total number, for both real and spectator particles
    int old_total = update_total_num();

    // target bunch is empty.  Set the weights from the injected bunch
    if (old_total == 0) real_num = o.get_real_num();
}

#if 0
Diagnostics_worker & 
Bunch::get_diag(std::string const & name)
{ 
    auto it = diags.find(name);
    if (it == diags.end()) throw std::runtime_error("cannot find diagnostics " + name);
    return it->second;
}
#endif

template<>
void
Bunch::print_statistics(Logger& logger) const
{
    using LV = LoggerV;

    logger(LV::DEBUG)
        << "Bunch statistics: "
        << "num_valid = " << get_local_num()
        << ", size = " << size()
        << ", capacity = " << capacity()
        << ", total_num = " << get_total_num()
        <<"\nMean and std: ";

    // print particles after propagate
    auto mean = Core_diagnostics::calculate_mean(*this);
    auto std  = Core_diagnostics::calculate_std(*this, mean);

    logger(LV::DEBUG) 
        << std::resetiosflags(std::ios::fixed)
        << std::setprecision(16)
        << std::setiosflags(std::ios::showpos | std::ios::scientific)
        << "\n"
        //<< "\nmean\tstd\n"
        ;

    for (int i=0; i<6; ++i) 
        logger(LV::DEBUG) << mean[i] << ", " << std[i] << "\n";

    logger(LV::DEBUG)
        << std::resetiosflags(std::ios::showpos | std::ios::scientific)
        << "\n";

#if 0
    for (int p=0; p<4; ++p) print_particle(p, logger);
    logger << "\n";
#endif
}



