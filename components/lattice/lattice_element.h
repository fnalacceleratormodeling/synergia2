#ifndef LATTICE_ELEMENT_H_
#define LATTICE_ELEMENT_H_

#include <string>
#include <beamline/bmlnElmnt.h>

class Lattice_element
{
private:
	ElmPtr chef_elmptr;

public:
	Lattice_element(ElmPtr chef_elmptr);
	std::string get_name();
	std::string get_type();
//	double length();
//	double quadrupole_strength();
};

#endif /* LATTICE_ELEMENT_H_ */
