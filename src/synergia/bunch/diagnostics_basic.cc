#include "diagnostics_basic.h"
#include "core_diagnostics.h"
#include <cmath>

const char Diagnostics_basic::name[] = "diagnostics_basic";

Diagnostics_basic::Diagnostics_basic(std::string const& filename, std::string const& local_dir) : 
    Diagnostics_basic::Diagnostics(Diagnostics_basic::name, filename, local_dir), 
    have_writers(false), 
    writer_s_n(0), 
    writer_repetition(0), 
    writer_s(0), 
    writer_num_particles(0), 
    writer_real_num_particles(0), 
    writer_pz(0), 
    mean(boost::extents[6]), 
    writer_mean(0), 
    std(boost::extents[6]), 
    writer_std(0), 
    min(boost::extents[3]), 
    writer_min(0), 
    max(boost::extents[3]), 
    writer_max(0)
{
}

Diagnostics_basic::Diagnostics_basic() : 
    have_writers(false)
{
}

bool
Diagnostics_basic::is_serial() const
{
    return true;
}

void
Diagnostics_basic::update()
{
    if (get_bunch().get_comm().has_this_rank())
    {
	    get_bunch().convert_to_state(get_bunch().fixed_z_lab);

	    s_n = get_bunch().get_reference_particle().get_s_n();
	    repetition = get_bunch().get_reference_particle().get_repetition();
	    s = get_bunch().get_reference_particle().get_s();
        pz = get_bunch().get_reference_particle().get_momentum();

	    num_particles = get_bunch().get_total_num();
	    real_num_particles = get_bunch().get_real_num();

	    mean = Core_diagnostics::calculate_mean(get_bunch());
	    std = Core_diagnostics::calculate_std(get_bunch(), mean);
	    min = Core_diagnostics::calculate_min(get_bunch());
	    max = Core_diagnostics::calculate_max(get_bunch());
    }
}

double
Diagnostics_basic::get_s_n() const
{
    return s_n;
}

int
Diagnostics_basic::get_repetition() const
{
    return repetition;
}

double
Diagnostics_basic::get_s() const
{
    return s;
}

int
Diagnostics_basic::get_num_particles() const
{
    return num_particles;
}

double
Diagnostics_basic::get_real_num_particles() const
{
    return real_num_particles;
}

double
Diagnostics_basic::get_pz() const
{
    return pz;
}

Const_MArray1d_ref
Diagnostics_basic::get_mean() const
{
    return mean;
}

Const_MArray1d_ref
Diagnostics_basic::get_std() const
{
    return std;
}

const MArray1d
Diagnostics_basic::get_min() const
{
    return min;
}

const MArray1d
Diagnostics_basic::get_max() const
{
    return max;
}

void
Diagnostics_basic::init_writers(Hdf5_file_sptr file_sptr)
{
    if (!have_writers) 
    {
        Four_momentum fourp( get_bunch().get_reference_particle().get_four_momentum() );

        int chg = get_bunch().get_reference_particle().get_charge();
        file_sptr->write(chg, "charge");

        double pmass = fourp.get_mass();
        file_sptr->write(pmass, "mass");

        writer_s_n = new Hdf5_serial_writer<double> (file_sptr, "s_n");
        writer_repetition = new Hdf5_serial_writer<int> (file_sptr, "repetition");
        writer_s = new Hdf5_serial_writer<double> (file_sptr, "s");
        writer_num_particles = new Hdf5_serial_writer<int> (file_sptr, "num_particles");
        writer_real_num_particles = new Hdf5_serial_writer<double> (file_sptr, "real_num_particles");
        writer_pz = new Hdf5_serial_writer<double> (file_sptr, "pz");
        writer_mean = new Hdf5_serial_writer<MArray1d_ref> (file_sptr, "mean");
        writer_std = new Hdf5_serial_writer<MArray1d_ref> (file_sptr, "std");
        writer_min = new Hdf5_serial_writer<MArray1d_ref> (file_sptr, "min");
        writer_max = new Hdf5_serial_writer<MArray1d_ref> (file_sptr, "max");

        have_writers = true;
    }
}

void
Diagnostics_basic::write()
{
    if (get_bunch().get_comm().has_this_rank())
    {
        if (get_write_helper().write_locally()) 
        {
            init_writers(get_write_helper().get_hdf5_file_sptr());

            writer_s_n->append(s_n);
            writer_repetition->append(repetition);
            writer_s->append(s);
            writer_num_particles->append(num_particles);
            writer_real_num_particles->append(real_num_particles);
            writer_pz->append(pz);
            writer_mean->append(mean);
            writer_std->append(std);
            writer_min->append(min);
            writer_max->append(max);

            get_write_helper().finish_write();
        }
    }
}

template<class Archive>
void
Diagnostics_basic::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Diagnostics);
    ar & BOOST_SERIALIZATION_NVP(have_writers);
    ar & BOOST_SERIALIZATION_NVP(s_n);
    ar & BOOST_SERIALIZATION_NVP(writer_s_n);
    ar & BOOST_SERIALIZATION_NVP(repetition);
    ar & BOOST_SERIALIZATION_NVP(writer_repetition);
    ar & BOOST_SERIALIZATION_NVP(s);
    ar & BOOST_SERIALIZATION_NVP(writer_s);
    ar & BOOST_SERIALIZATION_NVP(num_particles);
    ar & BOOST_SERIALIZATION_NVP(writer_num_particles);
    ar & BOOST_SERIALIZATION_NVP(real_num_particles);
    ar & BOOST_SERIALIZATION_NVP(writer_real_num_particles);
    ar & BOOST_SERIALIZATION_NVP(writer_pz);
    ar & BOOST_SERIALIZATION_NVP(mean);
    ar & BOOST_SERIALIZATION_NVP(writer_mean);
    ar & BOOST_SERIALIZATION_NVP(std);
    ar & BOOST_SERIALIZATION_NVP(writer_std);
    ar & BOOST_SERIALIZATION_NVP(min);
    ar & BOOST_SERIALIZATION_NVP(writer_min);
    ar & BOOST_SERIALIZATION_NVP(max);
    ar & BOOST_SERIALIZATION_NVP(writer_max);
}

template
void
Diagnostics_basic::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics_basic::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics_basic::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics_basic::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Diagnostics_basic::~Diagnostics_basic()
{
    if (have_writers) {
        delete writer_max;
        delete writer_min;
        delete writer_std;
        delete writer_mean;
        delete writer_pz;
        delete writer_real_num_particles;
        delete writer_num_particles;
        delete writer_s;
        delete writer_repetition;
        delete writer_s_n;
    }
}

BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics_basic)


