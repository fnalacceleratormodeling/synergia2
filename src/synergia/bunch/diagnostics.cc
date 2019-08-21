#include "diagnostics.h"
#include "synergia/bunch/bunch.h"
#include <stdexcept>

Diagnostics ::Diagnostics(
        std::string const& name, bool serial,
        std::string const& filename,
        std::string const& local_dir) 
    : name(name)
    , filename(filename)
    , local_dir(local_dir)
    , serial(serial)
    , bunch(nullptr)
    , writers()
{
}

void
Diagnostics::delete_write_helper(std::string const & name)
{
    writers.erase(name);
}

bool
Diagnostics::have_write_helper(std::string const & name) const
{
    return (writers.find(name) != writers.end());
}

Diagnostics_write_helper &
Diagnostics::get_write_helper(std::string const & name)
{
    if (have_write_helper(name))
    {
        return writers.find(name)->second;
    }
    else
    {
        return writers.emplace( 
                name,
                Diagnostics_write_helper(
                    get_filename(),
                    is_serial(), 
                    get_bunch().get_comm(), 
                    local_dir,
                    (name == DEFAULT_WRITER_NAME) ? "" : name )
                ).first->second;
    }
}

#if 0
template<class Archive>
    void
    Diagnostics::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(name);
        ar & BOOST_SERIALIZATION_NVP(filename);
        ar & BOOST_SERIALIZATION_NVP(local_dir);
        ar & BOOST_SERIALIZATION_NVP(bunch_sptr);
        ar & BOOST_SERIALIZATION_NVP(have_bunch_);
        ar & BOOST_SERIALIZATION_NVP(write_helper_ptr);
        ar & BOOST_SERIALIZATION_NVP(have_write_helper_);
        ar & BOOST_SERIALIZATION_NVP(extra_writers);
    }

template
void
Diagnostics::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Diagnostics::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Diagnostics::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Diagnostics::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
#endif


