#ifndef DIAGNOSTICS_WRITER_H_
#define DIAGNOSTICS_WRITER_H_
#include "components/bunch/diagnostics.h"
#include <string>

class Diagnostics_writer
{
private:
    Diagnostics_sptr diagnostics_sptr;
    hid_t file;
    bool dummy;
    int count;
    std::string filename_base, filename_suffix;
    void
    open_file_and_init();
public:
    Diagnostics_writer(std::string const& filename,
            Diagnostics_sptr diagnostics_sptr);
    Diagnostics_writer();
    bool
    is_dummy() const;
    Diagnostics_sptr
    get_diagnostics_sptr();
    int
    get_count() const;
    void
    set_count(int count);
    void
    write();
    void
    update_and_write(Bunch const& bunch);
    ~Diagnostics_writer();
};

typedef boost::shared_ptr<Diagnostics_writer > Diagnostics_writer_sptr;

Diagnostics_writer
no_diagnostics();

#endif /* DIAGNOSTICS_WRITER_H_ */
