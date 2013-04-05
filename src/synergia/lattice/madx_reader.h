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
    get_line_names() const;
    std::vector<std::string >
    get_sequence_names() const;
    std::vector<std::string >
    get_all_names() const;
    double
    get_double_variable(std::string const& name) const;
    std::string
    get_string_variable(std::string const& name) const;
    Lattice_sptr
    get_lattice_sptr(std::string const& line_name);
    Lattice_sptr
    get_lattice_sptr(std::string const& line_name, std::string const& filename);
    ~MadX_reader();
};
#endif /* MADX_READER_CC_ */
