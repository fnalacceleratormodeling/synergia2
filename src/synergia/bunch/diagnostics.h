#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include <string>

#include <cereal/archives/json.hpp>
#include <cereal/types/polymorphic.hpp>

#include "synergia/utils/hdf5_file.h"

template <class PART>
class bunch_t;

using Bunch = bunch_t<double>;

class Diagnostics {
  private:
    std::string type_;
    std::string filename_;
    bool single_file_;
    bool first_write = true;

    virtual void do_update(Bunch const&) = 0;
    virtual void do_reduce(Commxx const&, int) = 0;
    virtual void do_write(Hdf5_file&) = 0;
    virtual void do_first_write(Hdf5_file&) = 0;

  public:
    Diagnostics(std::string const& type = "Diagnostics",
                std::string const& filename = "diag.h5",
                bool single_file = true)
        : type_(type), filename_(filename), single_file_(single_file)
    {}

    virtual ~Diagnostics() = default;

    std::string
    type() const
    {
        return type_;
    }
    std::string
    filename() const
    {
        return filename_;
    }

    void
    update(Bunch const& bunch)
    {
        do_update(bunch);
    }

    void
    reduce(Commxx const& comm, int root)
    {
        do_reduce(comm, root);
    }

    void
    write(Hdf5_file& file)
    {
        if (first_write) {
            do_first_write(file);
            first_write = false;
        }
        do_write(file);
    }

    bool
    serial() const
    {
        return single_file();
    }
    bool
    single_file() const
    {
        return single_file_;
    }

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(CEREAL_NVP(type_));
        ar(CEREAL_NVP(filename_));
        ar(CEREAL_NVP(single_file_));
        ar(CEREAL_NVP(first_write));
    }
};

class Diagnostics_dummy : public Diagnostics {

  public:
    Diagnostics_dummy() : Diagnostics("diag_dummy", "diag_dummy.h5", true) {}

  private:
    void
    do_update(Bunch const& bunch) override
    {}
    void
    do_reduce(Commxx const& comm, int writer_rank) override
    {}
    void
    do_first_write(Hdf5_file& file) override
    {}
    void
    do_write(Hdf5_file& file) override
    {}

    friend class cereal::access;

    template <class Archive>
    void
    serialize(Archive& ar)
    {
        ar(cereal::base_class<Diagnostics>(this));
    }
};

CEREAL_REGISTER_TYPE(Diagnostics_dummy)

#endif /* DIAGNOSTICS_H_ */
