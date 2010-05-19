#ifndef DIAGNOSTICS_WRITER_H_
#define DIAGNOSTICS_WRITER_H_
#include "components/bunch/diagnostics.h"
#include <string>

class Diagnostics_writer
{
    Diagnostics_sptr diagnostics_sptr;
    hid_t file;
public:
    Diagnostics_writer(std::string const& filename, Diagnostics_sptr const& diagnostics_sptr);
    Diagnostics_sptr &
    get_diagnostics_sptr();
    void
    write();
    void
    update_and_write(Bunch const& bunch);
    ~Diagnostics_writer();
};

#endif /* DIAGNOSTICS_WRITER_H_ */
