#include <iostream>

#include "lattice.h"
#include "four_momentum.h"

int
main(int argc, char **argv)
{
	if (argc != 3) {
		std::cout << "usage: demo <lattice file> <line name>\n";
		std::exit(1);
	}

	std::cout << "starting demo...\n";

	Lattice lattice(argv[1],argv[2]);
	std::list<Lattice_element> elements(lattice.get_elements());
	for(std::list<Lattice_element>::iterator it = elements.begin();
		it != elements.end();
		++it) {
		std::cout << it->get_name() << ": " << it->get_type() << std::endl;
	}

	Four_momentum mom4(0.938,0.0);
	mom4.set_momentum(100.0);
	std::cout << "mom4.energy = " << mom4.get_total_energy() << std::endl;
	std::cout << "success!\n";
	return 0;
}
