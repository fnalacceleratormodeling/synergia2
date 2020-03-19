
#include "bunch.h"

#include <iostream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <sstream>

#if 0
// particle padding based on GSVector settings
#if defined(GSV_SSE)
  const int Bunch::particle_alignment = 4;
#elif defined(GSV_AVX)
  const int Bunch::particle_alignment = 4;
#elif defined(GSV_AVX512)
  const int Bunch::particle_alignment = 8;
#elif defined(GSV_QPX)
  const int Bunch::particle_alignment = 4;
#else
  const int Bunch::particle_alignment = 4;
#endif
#endif


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

Bunch::Bunch(
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
    , parts{ BunchParticles("particles", total_num, -1, *comm),
             BunchParticles("spectators", total_spectator_num, -1, *comm) }
    , bunch_index(bunch_index)
    , bucket_index(bucket_index)
    , array_index(array_index)
{
}

Bunch::Bunch()
    : comm(new Commxx())
    , boundary(LB::open)
    , boundary_param(0.0)
    , ref_part()
    , design_ref_part()
    , particle_charge(ref_part.get_charge())
    , real_num(1.0)
    , parts{ BunchParticles("particles", 0, 0, *comm), 
             BunchParticles("spectators", 0, 0, *comm) }
    , bunch_index(0)
    , bucket_index(0)
    , array_index(0)
{
}

void
Bunch::inject(Bunch const& bunch)
{
    throw std::runtime_error("Bunch::inject not implemented");

#if 0
    const double weight_tolerance = 1.0e-14;
    const double particle_tolerance = 1.0e-14;

    // The charge and mass of the bunch particles must match
    if (particle_charge != bunch.get_particle_charge()) 
    {
        throw std::runtime_error(
                "Bunch.inject: bunch particle charges do not match.");
    }

    if (std::abs(reference_particle.get_four_momentum().get_mass()/
                 bunch.get_reference_particle().get_four_momentum().get_mass() - 1.0) > particle_tolerance) 
    {
        throw std::runtime_error(
                "Bunch:inject: bunch particle masses do not match.");
    }

    // can only check particle weight if total_num is nonzero
    if (total_num == 0) 
    {
        // target bunch is empty.  Set the weights from the injected bunch
        real_num = bunch.get_real_num();
        total_num = bunch.get_total_num();
    } 
    else 
    {
        double wgt1 = real_num/total_num;
        double wgt2 = bunch.get_real_num()/bunch.get_total_num();

        if (std::abs(wgt1/wgt2 - 1.0) > weight_tolerance) 
        {
            throw std::runtime_error(
                "Bunch.inject: macroparticle weight of injected bunch does not match.");
        }
    }

    ConstParticles injected_particles(bunch.get_local_particles());
    ConstParticles injected_spectator_particles(bunch.get_local_spectator_particles());

    double target_momentum = reference_particle.get_momentum();
    double injected_momentum = bunch.get_reference_particle().get_momentum();

    MArray1d ref_state_diff(boost::extents[6]);
    MArray1d target_state(boost::extents[6]);
    MArray1d injected_state(boost::extents[6]);

    for (int i = 0; i < 6; ++i) 
    {
        ref_state_diff[i] = bunch.get_reference_particle().get_state()[i]
                - reference_particle.get_state()[i];
    }

    for (int i = 0; i < 6; ++i) 
    {
        target_state[i] = reference_particle.get_state()[i];
        injected_state[i] = bunch.get_reference_particle().get_state()[i];
    }

    // real particles
    int old_local_num = local_num;
    int old_local_num_lost = get_local_num_lost();

    // reallocate the particle array
    expand_local_num(
            old_local_num + bunch.get_local_num(), 
            bunch.get_local_num_lost() );

    for (int part = 0; part < bunch.get_local_num(); ++part) 
    {
        // space-like coordinates
        for (int i = 0; i < 6; i += 2) 
        {
            (*local_particles)[old_local_num + part][i]
                    = injected_particles[part][i] + ref_state_diff[i];
        }

        // npx and npy coordinates are scaled with p_ref which can be different
        // for different bunches
        for (int i = 1; i < 4; i += 2) 
        {
            (*local_particles)[old_local_num + part][i] =
                    (injected_momentum/target_momentum) *
                    (injected_particles[part][i] - injected_state[i]) + target_state[i];
        }

        // ndp coordinate is delta-p scaled with pref
        (*local_particles)[old_local_num + part][5] =
                (injected_momentum/target_momentum) *
                (1.0 + injected_particles[part][5] - injected_state[5]) + target_state[5] - 1.0;

        (*local_particles)[old_local_num + part][Bunch::id]
                = injected_particles[part][Bunch::id];
    }

    // copy the lost particles from bunch
    for (int p = 0; p < bunch.get_local_num_lost(); ++p)
    {
        for (int i = 0; i < 6; i += 2) 
        {
            (*local_particles)[local_num_padded + old_local_num_lost + p][i] =
                injected_particles[bunch.get_local_num_padded() + p][i];
        }
    }

    // spectator particles
    int old_local_s_num = local_s_num;
    int old_local_s_num_lost = get_local_spectator_num_lost();

    // reallocate the spectator particle array
    expand_local_spectator_num(
            old_local_s_num + bunch.get_local_spectator_num(), 
            bunch.get_local_spectator_num_lost() );

    for (int part = 0; part < bunch.get_local_spectator_num(); ++part) 
    {
        // space-like coordinates
        for (int i = 0; i < 6; i += 2) 
        {
            (*local_s_particles)[old_local_s_num + part][i]
                    = injected_spectator_particles[part][i] + ref_state_diff[i];
        }

        // npx and npy coordinates are scaled with p_ref which can be different
        // for different bunches
        for (int i = 1; i < 4; i += 2) 
        {
            (*local_s_particles)[old_local_s_num + part][i] =
                    (injected_momentum/target_momentum) *
                    (injected_spectator_particles[part][i] - injected_state[i]) + target_state[i];
        }

        // ndp coordinate is delta-p scaled with pref
        (*local_s_particles)[old_local_s_num + part][5] =
                (injected_momentum/target_momentum) *
                (1.0 + injected_spectator_particles[part][5] - injected_state[5]) + target_state[5] - 1.0;

        (*local_s_particles)[old_local_s_num + part][Bunch::id]
                = injected_spectator_particles[part][Bunch::id];
    }

    // copy the lost spectator particles from bunch
    for (int p = 0; p < bunch.get_local_spectator_num_lost(); ++p)
    {
        for (int i = 0; i < 6; i += 2) 
        {
            (*local_s_particles)[local_s_num_padded + old_local_s_num_lost + p][i] =
                injected_spectator_particles[bunch.get_local_spectator_num_padded() + p][i];
        }
    }

    // update total number, for both real and spectator particles
    update_total_num();
#endif
}

Diagnostics_worker & 
Bunch::get_diag(std::string const & name)
{ 
    auto it = diags.find(name);
    if (it == diags.end()) throw std::runtime_error("cannot find diagnostics " + name);
    return it->second;
}

#if 0
template<class Archive>
void
Bunch::save(Archive & ar, const unsigned int version) const
{
    ar << CEREAL_NVP(longitudinal_extent)
       << CEREAL_NVP(z_periodic)
       << CEREAL_NVP(longitudinal_aperture)

       << CEREAL_NVP(reference_particle)
       << CEREAL_NVP(design_reference_particle)
       << CEREAL_NVP(particle_charge)

       << CEREAL_NVP(total_num)
       << CEREAL_NVP(total_s_num)
       << CEREAL_NVP(real_num)

       << CEREAL_NVP(bucket_index)
       << CEREAL_NVP(bucket_index_assigned)

       << CEREAL_NVP(comm_sptr)

       //<< CEREAL_NVP(state)
       //<< CEREAL_NVP(default_converter)
       //<< CEREAL_NVP(converter_ptr)
       ;

    if (comm_sptr->has_this_rank()) 
    {
        int attempts = 0;
        bool fail = true;

        while ((attempts<5) && fail)
        {
            try 
            {
                boost::filesystem::remove(get_local_particles_serialization_path());
                Hdf5_file file(get_local_particles_serialization_path(), Hdf5_file::truncate);

                file.write(local_num, "local_num");
                file.write(local_num_aligned, "local_num_aligned");
                file.write(local_num_padded, "local_num_padded");
                file.write(local_num_slots, "local_num_slots");
                file.write(storage, local_num_slots*7, "local_storage");

                file.write(local_s_num, "local_s_num");
                file.write(local_s_num_aligned, "local_s_num_aligned");
                file.write(local_s_num_padded, "local_s_num_padded");
                file.write(local_s_num_slots, "local_s_num_slots");
                file.write(s_storage, local_s_num_slots*7, "local_s_storage");

                file.close();
                fail=false;
            }
            catch(Hdf5_exception & he) 
            {
                ++attempts;
                fail=true;
                std::cout<<"bunch.cc: H5 Exception thrown, attempts number="
                    <<attempts<<" on rank="<<Commxx().get_rank()<<std::endl;
                sleep(3);
            }
        }
    }
}


template<class Archive>
void
Bunch::load(Archive & ar, const unsigned int version)
{
    ar >> CEREAL_NVP(longitudinal_extent)
       >> CEREAL_NVP(z_periodic)
       >> CEREAL_NVP(longitudinal_aperture)

       >> CEREAL_NVP(reference_particle)
       >> CEREAL_NVP(design_reference_particle)
       >> CEREAL_NVP(particle_charge)

       >> CEREAL_NVP(total_num)
       >> CEREAL_NVP(total_s_num)
       >> CEREAL_NVP(real_num)

       >> CEREAL_NVP(bucket_index)
       >> CEREAL_NVP(bucket_index_assigned)

       >> CEREAL_NVP(comm_sptr)

       //>> CEREAL_NVP(state)
       //>> CEREAL_NVP(default_converter)
       //>> CEREAL_NVP(converter_ptr)
       ;

    if (comm_sptr->has_this_rank()) 
    {
        Hdf5_file file(get_local_particles_serialization_path(), Hdf5_file::read_only);

        local_num = file.read<int> ("local_num");
        local_num_aligned = file.read<int> ("local_num_aligned");
        local_num_padded = file.read<int> ("local_num_padded");
        local_num_slots = file.read<int> ("local_num_slots");

        local_s_num = file.read<int> ("local_s_num");
        local_s_num_aligned = file.read<int> ("local_s_num_aligned");
        local_s_num_padded = file.read<int> ("local_s_num_padded");
        local_s_num_slots = file.read<int> ("local_s_num_slots");

        storage = file.read<double *>("local_storage");
        local_particles = new MArray2d_ref(storage, boost::extents[local_num_slots][7], boost::fortran_storage_order());

        s_storage = file.read<double *>("local_s_storage");
        local_s_particles = new MArray2d_ref(s_storage, boost::extents[local_s_num_slots][7], boost::fortran_storage_order());
    } 
    else 
    {
        local_num = 0;
        local_num_aligned = 0;
        local_num_padded = 0;
        local_num_slots = 0;

        local_s_num = 0;
        local_s_num_aligned = 0;
        local_s_num_padded = 0;
        local_s_num_slots = 0;

        storage = NULL;
        local_particles = new MArray2d_ref(storage, boost::extents[0][7], boost::fortran_storage_order());

        s_storage = NULL;
        local_s_particles = new MArray2d_ref(s_storage, boost::extents[0][7], boost::fortran_storage_order());
    }
}

template
void
Bunch::save<cereal::BinaryOutputArchive >(
        cereal::BinaryOutputArchive & ar, const unsigned int version) const;

template
void
Bunch::save<cereal::XMLOutputArchive >(
        cereal::XMLOutputArchive & ar, const unsigned int version) const;

template
void
Bunch::load<cereal::BinaryInputArchive >(
        cereal::BinaryInputArchive & ar, const unsigned int version);

template
void
Bunch::load<cereal::XMLInputArchive >(
        cereal::XMLInputArchive & ar, const unsigned int version);

#endif
