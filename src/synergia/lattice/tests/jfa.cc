#include <iostream>
#include "synergia/lattice/madx_reader.h"

int
main()
{
    MadX_reader madx_reader;
    madx_reader.parse_file("ps_lattice.madx");
    Lattice_sptr lattice_sptr(madx_reader.get_lattice_sptr("ps"));
    lattice_sptr->print();
}
