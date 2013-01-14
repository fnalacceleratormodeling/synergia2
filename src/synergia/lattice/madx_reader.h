#ifndef MADX_READER_CC_
#define MADX_READER_CC_
#include "synergia/lattice/madx.h"
#include "synergia/lattice/lattice.h"

class MadX_reader
{
private:
    boost::shared_ptr<synergia::MadX > madx_sptr;
public:
    MadX_reader();
    void
    parse(std::string const& string);
    void
    parse_file(std::string const& filename);
    std::vector<std::string >
    get_line_names();
    Lattice
    get_lattice(std::string const& line_name);
    Lattice
    get_lattice(std::string const& line_name, std::string const& filename);
    ~MadX_reader();
};
#endif /* MADX_READER_CC_ */
