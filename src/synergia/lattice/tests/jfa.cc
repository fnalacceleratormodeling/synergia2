#include <iostream>
#include "synergia/lattice/madx_reader.h"

int
main()
{
    MadX_reader madx_reader;
    madx_reader.parse_file("lattices/fodo2work.madx");
    Lattice lattice(madx_reader.get_lattice("fodo"));
    lattice.print();
}