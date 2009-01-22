#include "lattice.h"

#include <iostream>
#include <algorithm>
#include <cctype>
#include <stdlib.h>

#include <parsers/xsif/XSIFFactory.h>
#include <bmlfactory/MAD8Factory.h>

#include <basic_toolkit/PhysicsConstants.h>
#include <beamline/Particle.h>

extern beamline* DriftsToSlots(beamline const& beamline_drifts);

std::string
upcase(const std::string str)
{
	std::string retval(str);
	// explicit cast needed to resolve ambiguity in overloaded versions
	// of std::toupper
	std::transform(retval.begin(), retval.end(), retval.begin(),
			static_cast<int(*)(int)>(std::toupper));
	return retval;
}

Lattice_element::Lattice_element(ElmPtr chef_elmptr)
{
	this->chef_elmptr = chef_elmptr;
}

std::string
Lattice_element::get_name()
{
	return std::string(chef_elmptr->Name());
}

Lattice::Lattice(std::string const &file_name,
		std::string const &line_name)
{
	XSIFFactory factory(file_name);
	chef_bmlptr = BmlPtr(DriftsToSlots(*factory.create_beamline(
			upcase(line_name).c_str())));
	int nelm = chef_bmlptr->countHowManyDeeply();
		std::cout << "nelm = " << nelm << std::endl;

		for (beamline::deep_iterator it = chef_bmlptr->deep_begin(); it != chef_bmlptr->deep_end(); ++it ) {
		   elements.push_back(*it);
		 }
}

std::list<Lattice_element> &
Lattice::get_elements()
{
	return elements;
}
