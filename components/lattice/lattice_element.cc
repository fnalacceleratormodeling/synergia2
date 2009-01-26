#include "lattice_element.h"

Lattice_element::Lattice_element(ElmPtr chef_elmptr)
{
	this->chef_elmptr = chef_elmptr;
}

std::string
Lattice_element::get_name()
{
	return std::string(chef_elmptr->Name());
}

std::string
Lattice_element::get_type()
{
	return std::string(chef_elmptr->Type());
}
