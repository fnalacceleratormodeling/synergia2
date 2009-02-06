#include "macro_bunch_store.h"
#include <numpy/arrayobject.h>
#include "array_nd/array_nd_python.h"
#include <iostream>
using namespace boost::python;

Macro_bunch_store::Macro_bunch_store(numeric::array& numeric_local_particles,
                                     int local_num, int total_num,
                                     double mass, int charge, double total_current,
                                     numeric::array& numeric_units,
                                     numeric::array& numeric_ref_particle,
                                     bool is_fixedz)
{
    this->numeric_local_particles = &numeric_local_particles;
    local_particles = Array_nd_from_PyObject<double>(numeric_local_particles.ptr());
    this->numeric_units = &numeric_units;
    units = Array_nd_from_PyObject<double>(numeric_units.ptr());
    this->numeric_ref_particle = &numeric_ref_particle;
    ref_particle = Array_nd_from_PyObject<double>(numeric_ref_particle.ptr());
    this->local_num = local_num;
    this->total_num = total_num;
    this->mass = mass;
    this->charge = charge;
    this->total_current = total_current;
    this->is_fixedz = is_fixedz;
}

numeric::array Macro_bunch_store::get_local_particles()
{
    return *(this->numeric_local_particles);
}

numeric::array Macro_bunch_store::get_units()
{
    return *(this->numeric_units);
}

numeric::array Macro_bunch_store::get_ref_particle()
{
    return *(this->numeric_ref_particle);
}

void
Macro_bunch_store::check(numeric::array &array)
{
    std::cout << "stored value: " << numeric_ref_particle << std::endl;
    std::cout << "passed value: " << &array << std::endl;
}

void Macro_bunch_store::convert_to_fixedt()
{
    if (is_fixedz) {
        double gamma = -ref_particle(5);
        for (int i = 0; i < local_num; ++i) {
            local_particles(0, i) /= units(0);
            local_particles(2, i) /= units(2);
            double xp = local_particles(1, i);
            double yp = local_particles(3, i);
            double rcp_gammai = 1.0 / (gamma - local_particles(5, i));
            double betai = sqrt(1.0 - rcp_gammai * rcp_gammai * (1 + xp * xp + yp * yp));
            local_particles(4, i) *= -gamma * betai / units(0); // units(0)
            // is not
            // an error!
        }
        is_fixedz = false;
    }
}

void Macro_bunch_store::convert_to_fixedz()
{
    if (! is_fixedz) {
        double gamma = -ref_particle(5);
        for (int i = 0; i < local_num; ++i) {
            local_particles(0, i) *= units(0);
            local_particles(2, i) *= units(2);
            double xp = local_particles(1, i);
            double yp = local_particles(3, i);
            double rcp_gammai = 1.0 / (gamma - local_particles(5, i));
            double betai = sqrt(1.0 - rcp_gammai * rcp_gammai * (1 + xp * xp + yp * yp));
            local_particles(4, i) /= -gamma * betai / units(0); // units(0)
            // is not
            // an error!
        }
        is_fixedz = true;
    }
}

double Macro_bunch_store::get_coord(int coord_index, int particle_index)
{
    return local_particles(coord_index, particle_index);
}

void Macro_bunch_store::set_local_particles(numeric::array &local_particles)
{
    this->numeric_local_particles = &local_particles;
    this->local_particles = Array_nd_from_PyObject<double>(local_particles.ptr());
}

Macro_bunch_store::~Macro_bunch_store()
{}
