#ifndef MADX_READER_CC_
#define MADX_READER_CC_
#include "synergia/lattice/madx.h"
#include "synergia/lattice/lattice.h"
#include "synergia/lattice/element_adaptor_map.h"

class MadX_reader
{
private:
    boost::shared_ptr<synergia::MadX > madx_sptr;
    Element_adaptor_map_sptr element_adaptor_map_sptr;
    void
    extract_reference_particle(Lattice & lattice);
public:
    MadX_reader();
    MadX_reader(Element_adaptor_map_sptr element_adaptor_map_sptr);
    void
    parse(std::string const& string);
    void
    parse_file(std::string const& filename);
    std::vector<std::string >
    get_line_names();
    Lattice_sptr
    get_lattice_sptr(std::string const& line_name);
    Lattice_sptr
    get_lattice_sptr(std::string const& line_name, std::string const& filename);
    ~MadX_reader();
};
#endif /* MADX_READER_CC_ */
