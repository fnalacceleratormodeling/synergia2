#include "diagnostics_writer.h"

Diagnostics_writer::Diagnostics_writer(std::string const& filename,
        Diagnostics_sptr const& diagnostics_sptr) :
    diagnostics_sptr(diagnostics_sptr), dummy(false)
{
    file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    diagnostics_sptr->init_writers(file);
}

Diagnostics_writer::Diagnostics_writer() :
    dummy(true)
{
}

bool
Diagnostics_writer::is_dummy() const
{
    return dummy;
}

Diagnostics_sptr &
Diagnostics_writer::get_diagnostics_sptr()
{
    return diagnostics_sptr;
}

void
Diagnostics_writer::write()
{
    if (!dummy) {
        diagnostics_sptr->write_hdf5();
    }
}

void
Diagnostics_writer::update_and_write(Bunch const& bunch)
{
    if (!dummy) {
        diagnostics_sptr->update(bunch);
        diagnostics_sptr->write_hdf5();
    }
}

Diagnostics_writer::~Diagnostics_writer()
{
    if (!dummy) {
        H5Fclose(file);
    }
}

Diagnostics_writer
no_diagnostics()
{
    return Diagnostics_writer();
}
