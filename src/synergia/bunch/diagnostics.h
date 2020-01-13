#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include <string>
#include <map>

#include "synergia/foundation/diagnostics_write_helper.h"
#include "synergia/utils/hdf5_serial_writer.h"

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
//#include <cereal/archives/binary.hpp>
//#include <cereal/archives/xml.hpp>

class Bunch;

class Diagnostics
{
private:

    std::string type_;
    bool serial_;
    bool first_write = true;

    virtual void do_update(Bunch const&) = 0;
    virtual void do_collect(Commxx, int) = 0;
    virtual void do_write(Hdf5_file&, bool) = 0;

public:

    Diagnostics(std::string const& type = "Diagnostics", bool serial = true)
        : type_(type), serial_(serial)
    { }

    virtual ~Diagnostics() = default;

    std::string type() const { return type_; }

    void update(Bunch const& bunch)
    { do_update(bunch); }

    void collect(Commxx comm, int root)
    { do_collect(comm, root); }

    void write(Hdf5_file & file)
    { do_write(file, first_write); first_write = false; }

    bool serial() const { return serial_; }

    template<class Archive>
    void serialize(Archive & ar)
    {
        ar(CEREAL_NVP(type_));
        ar(CEREAL_NVP(serial_));
        ar(CEREAL_NVP(first_write));
    }
};

class Diagnostics_dummy : public Diagnostics
{

public:

    Diagnostics_dummy() : Diagnostics("diag_dummy", true)
    { }

private:

    void do_update(Bunch const& bunch) override { }
    void do_collect(Commxx comm, int writer_rank) override { }
    void do_write(Hdf5_file & file, bool first_write) override { }

    friend class cereal::access;
    
    template<class Archive>
    void serialize(Archive & ar)
    { ar(cereal::base_class<Diagnostics>(this)); }
};

CEREAL_REGISTER_TYPE(Diagnostics_dummy)

#if 0
/// Diagnostics is an abstract base class for bunch diagnostics classes
class Diagnostics
{

private:

    constexpr static const char * DEFAULT_WRITER_NAME = "__default_writer__";

    std::string name;
    std::string filename;
    std::string local_dir;

    bool serial;

    std::map<std::string, Diagnostics_write_helper> writers;

private:

    virtual void do_update(Bunch const&) = 0;
    virtual void do_write (Bunch const&) = 0;
    virtual std::unique_ptr<Diagnostics> do_pilfer() = 0;

public:

    Diagnostics( std::string const& name, bool serial,
                 std::string const& filename, 
                 std::string const& local_dir="" );

    virtual ~Diagnostics() = default;
    Diagnostics(Diagnostics &&) noexcept = default;

    std::string const& get_filename()  const { return filename; }
    std::string const& get_local_dir() const { return local_dir; }

    // extra write helpers
    void delete_write_helper(std::string const & name = DEFAULT_WRITER_NAME);
    bool have_write_helper(std::string const& name = DEFAULT_WRITER_NAME) const;

    Diagnostics_write_helper &
    get_write_helper(Bunch const& bunch, std::string const& name = DEFAULT_WRITER_NAME);

    /// Multiple serial diagnostics can be written to a single file.
    bool is_serial() const
    { return serial; }

    /// Update the diagnostics
    void update(Bunch const& bunch)
    { do_update(bunch); }

    /// Write the diagnostics to the file
    void write(Bunch const& bunch)
    { do_write(bunch); }

    /// Update the diagnostics and write them to the file
    void update_and_write(Bunch const& bunch)
    { do_update(bunch); do_write(bunch); }

    /// Create a unique_ptr by moving the current diagnostics.
    /// the current diagnostics after pilfer() will be invalid
    std::unique_ptr<Diagnostics> pilfer()
    { return do_pilfer(); }

    template<class Archive>
    void serialize(Archive & ar)
    {
        ar(CEREAL_NVP(name));
        ar(CEREAL_NVP(filename));
        ar(CEREAL_NVP(local_dir));
        ar(CEREAL_NVP(serial));
        //ar(CEREAL_NVP(writers));
    }
};
#endif

#endif /* DIAGNOSTICS_H_ */
