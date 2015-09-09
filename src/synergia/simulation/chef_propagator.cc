#include "chef_propagator.h"
#include "synergia/lattice/chef_utils.h"

Chef_propagator::Chef_propagator(
        Chef_lattice_section_sptr chef_lattice_section_sptr) :
    chef_lattice_section_sptr(chef_lattice_section_sptr)
{
}

Chef_propagator::Chef_propagator()
{
}

// jfa: This routine is incorrect when passing through an accelerating element.
// Please fix it.
void
Chef_propagator::apply(Bunch & bunch, int verbosity, Logger & logger)
{

    Particle
            particle(
                    reference_particle_to_chef_particle(
                            bunch.get_reference_particle()));
    double length = 0.0;
    if (verbosity > 4) {
        logger << "Chef_propagator: begin_index = "
                << chef_lattice_section_sptr->get_begin_index()
                << ", end_index = "
                << chef_lattice_section_sptr->get_end_index() << std::endl;
    }

    bool has_rf_cavity = false;
    double initial_reference_energy = particle.ReferenceEnergy();
    if (verbosity > 4) {
        logger << "reference particle energy: " << initial_reference_energy << std::endl;
    }
    for (Chef_lattice_section::iterator it = chef_lattice_section_sptr->begin(); it
            != chef_lattice_section_sptr->end(); ++it) {
        if (verbosity > 5) {
			logger << "-->" << (*it)->Type() << "<--" << std::endl;
		}
        if (strcmp((*it)->Type(), "thinrfcavity") == 0) {
            has_rf_cavity = true;
        }
        (*it)->propagate(particle);
        double this_length = (*it)->OrbitLength(particle);
        length += this_length;
        if (verbosity > 4) {
            logger << "Chef_propagator: name = " << (*it)->Name() << ", type = " <<
                    (*it)->Type() << ", length = " << this_length << ", ReferenceTime: " << (*it)->getReferenceTime() << std::endl;
        }
    }
    bunch.get_reference_particle().increment_trajectory(length);
    double final_reference_energy = particle.ReferenceEnergy();

    if ((verbosity > 4) && has_rf_cavity) {
        logger << "Chef_propagator, lattice section has rf cavity, energy gain: " << final_reference_energy-initial_reference_energy << std::endl;
    }

    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    #pragma omp parallel for firstprivate(particle)
    for (int part = 0; part < local_num; ++part) {
        // propagation through an RF cavity cavity section changes the ReferenceEnergy of particle,
        // so I have to reset it each loop
        particle.SetReferenceEnergy(initial_reference_energy);
        particle.set_x(particles[part][Bunch::x]);
        particle.set_npx(particles[part][Bunch::xp]);
        particle.set_y(particles[part][Bunch::y]);
        particle.set_npy(particles[part][Bunch::yp]);
        particle.set_cdt(particles[part][Bunch::cdt]);
        particle.set_ndp(particles[part][Bunch::dpop]);

        for (Chef_lattice_section::iterator it =
                chef_lattice_section_sptr->begin(); it
                != chef_lattice_section_sptr->end(); ++it) {
            (*it)->propagate(particle);
        }

        particles[part][Bunch::x] = particle.get_x();
        particles[part][Bunch::xp] = particle.get_npx();
        particles[part][Bunch::y] = particle.get_y();
        particles[part][Bunch::yp] = particle.get_npy();
        particles[part][Bunch::cdt] = particle.get_cdt();
        particles[part][Bunch::dpop] = particle.get_ndp();
    }

    // update the reference particle if there was an energy change
    if (has_rf_cavity &&
        (final_reference_energy-initial_reference_energy != 0.0)) {
        bunch.get_reference_particle().set_total_energy(final_reference_energy);
    }
}

template<class Archive>
    void
    Chef_propagator::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(chef_lattice_section_sptr);
    }

template
void
Chef_propagator::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Chef_propagator::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Chef_propagator::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Chef_propagator::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
