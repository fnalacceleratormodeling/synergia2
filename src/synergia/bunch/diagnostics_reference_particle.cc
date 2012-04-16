#include "diagnostics_reference_particle.h"

const char Diagnostics_reference_particle::name[] = "diagnostics_reference_particle";

Diagnostics_reference_particle::Diagnostics_reference_particle(
        std::string const& filename) :
            Diagnostics_reference_particle::Diagnostics(
                    Diagnostics_reference_particle::name, filename),
            have_writers(false), writer_beta(0), writer_gamma(0),
            writer_state(0), writer_s(0)
{
}

Diagnostics_reference_particle::Diagnostics_reference_particle()
{
}

bool
Diagnostics_reference_particle::is_serial() const
{
    return true;
}

void
Diagnostics_reference_particle::update()
{
}

void
Diagnostics_reference_particle::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) {
        writer_beta = new Hdf5_serial_writer<double > (file_sptr, "beta");
        writer_gamma = new Hdf5_serial_writer<double > (file_sptr, "gamma");
        writer_state = new Hdf5_serial_writer<MArray1d_ref > (file_sptr,
                "state");
        writer_s = new Hdf5_serial_writer<double > (file_sptr, "s");
        have_writers = true;
    }
}

void
Diagnostics_reference_particle::write()
{
    if (get_write_helper().write_locally()) {
        init_writers(get_write_helper().get_hdf5_file_sptr());
        double beta = get_bunch().get_reference_particle().get_beta();
        writer_beta->append(beta);
        double gamma = get_bunch().get_reference_particle().get_gamma();
        writer_gamma->append(gamma);
        MArray1d state(get_bunch().get_reference_particle().get_state());
        writer_state->append(state);
        double s = get_bunch().get_reference_particle().get_s();
        writer_s->append(s);
        get_write_helper().finish_write();
    }
}

template<class Archive>
    void
    Diagnostics_reference_particle::serialize(Archive & ar,
            const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics)
                & BOOST_SERIALIZATION_NVP(have_writers)
                & BOOST_SERIALIZATION_NVP(writer_beta)
                & BOOST_SERIALIZATION_NVP(writer_gamma)
                & BOOST_SERIALIZATION_NVP(writer_state)
                & BOOST_SERIALIZATION_NVP(writer_s);
    }

template
void
Diagnostics_reference_particle::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_reference_particle::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_reference_particle::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_reference_particle::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_reference_particle::~Diagnostics_reference_particle()
{
    if (have_writers) {
        delete writer_beta;
        delete writer_gamma;
        delete writer_state;
        delete writer_s;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_reference_particle)

