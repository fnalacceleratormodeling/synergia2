#ifndef LATTICE_H_
#define LATTICE_H_

#include <string>
#include <list>

#include <beamline/beamline.h>

class Lattice_element
{
private:
	ElmPtr chef_elmptr;

public:
	Lattice_element(ElmPtr chef_elmptr);
	std::string get_name();
//	std::string type();
//	double length();
//	double quadrupole_strength();
};

class Lattice
{
private:
//	std::string file_name;
	BmlPtr chef_bmlptr;
	std::list<Lattice_element> elements;

public:
	Lattice(std::string const &file_name,
			std::string const &line_name);
	std::list<Lattice_element> &get_elements();
};



#endif /* LATTICE_H_ */
