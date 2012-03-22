#include "diagnostics.h"
#include "synergia/utils/hdf5_writer.h"
#include "synergia/utils/hdf5_chunked_array2d_writer.h"
#include <cmath>
#include "synergia/utils/eigen2/Eigen/Core"
#include "synergia/utils/eigen2/Eigen/LU"
#include <stdexcept>
#include "synergia/utils/simple_timer.h"

// import most common Eigen types
USING_PART_OF_NAMESPACE_EIGEN

Diagnostics::Diagnostics(std::string const& name, std::string const& filename) :
    name(name), filename(filename), have_bunch_(false), write_helper_ptr(0),
            have_write_helper_(0)
{
}

std::string const&
Diagnostics::get_filename() const
{
    return filename;
}

void
Diagnostics::set_bunch_sptr(Bunch_sptr bunch_sptr)
{
    this->bunch_sptr = bunch_sptr;
    have_bunch_ = true;
}

bool
Diagnostics::have_bunch() const
{
    return have_bunch_;
}

Bunch &
Diagnostics::get_bunch()
{
    if (!have_bunch_) {
        throw std::runtime_error(name + ": bunch not set");
    }
    return *bunch_sptr;
}

void
Diagnostics::delete_write_helper_ptr()
{
    if (have_write_helper_) {
        delete write_helper_ptr;
        have_write_helper_ = false;
    }
}

Diagnostics_write_helper *
Diagnostics::new_write_helper_ptr()
{
    delete_write_helper_ptr();
    return new Diagnostics_write_helper(get_filename(),
            is_serial(), get_bunch().get_comm());
}

bool
Diagnostics::have_write_helper() const
{
    return have_write_helper_;
}

Diagnostics_write_helper &
Diagnostics::get_write_helper()
{
    if (!have_write_helper_) {
        write_helper_ptr = new_write_helper_ptr();
        have_write_helper_ = true;
    }
    return *write_helper_ptr;
}

Diagnostics::Diagnostics()
{
}

Diagnostics::~Diagnostics()
{
    if (have_write_helper_) {
        delete write_helper_ptr;
    }
}
BOOST_CLASS_EXPORT_IMPLEMENT(Diagnostics)
