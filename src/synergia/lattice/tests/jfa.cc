#include <iostream>
#include "synergia/lattice/madx_reader.h"

int
main()
{
    MadX_reader madx_reader;
    madx_reader.parse_file("ps_lattice.madx");
    Lattice lattice(madx_reader.get_lattice("ps"));
    lattice.print();
}
