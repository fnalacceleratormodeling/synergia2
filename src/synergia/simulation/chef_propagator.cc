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
    for (Chef_lattice_section::iterator it = chef_lattice_section_sptr->begin(); it
            != chef_lattice_section_sptr->end(); ++it) {
        (*it)->propagate(particle);
        double this_length = (*it)->OrbitLength(particle);
        length += this_length;
        if (verbosity > 4) {
            logger << "Chef_propagator: name = " << (*it)->Name() << ", type = " <<
                    (*it)->Type() << ", length = " << this_length << std::endl;
        }
    }
    // std::cout<<"chef operate apply with length= "<<length<<std::endl;
    bunch.get_reference_particle().increment_trajectory(length);

    std::vector<double > u = chef_unit_conversion(
            bunch.get_reference_particle());
    int local_num = bunch.get_local_num();
    MArray2d_ref particles = bunch.get_local_particles();

    Vector chef_state(6);
    for (int part = 0; part < local_num; ++part) {
        for (int synergia_index = 0; synergia_index < 6; ++synergia_index) {
            int chef_idx = get_chef_index(synergia_index);
            chef_state[chef_idx] = particles[part][synergia_index]
                    / u[synergia_index];
        }

        particle.State() = chef_state;

        for (Chef_lattice_section::iterator it =
                chef_lattice_section_sptr->begin(); it
                != chef_lattice_section_sptr->end(); ++it) {
            (*it)->propagate(particle);
        }
        chef_state = particle.State();

        for (int synergia_index = 0; synergia_index < 6; ++synergia_index) {
            int chef_idx = get_chef_index(synergia_index);
            particles[part][synergia_index] = chef_state[chef_idx]
                    * u[synergia_index];
        }

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
