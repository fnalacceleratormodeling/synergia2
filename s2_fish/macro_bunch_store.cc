#include "macro_bunch_store.h"
#include <numpy/arrayobject.h>
#include "array_nd/array_nd_python.h"
#include <iostream>
using namespace boost::python;

Macro_bunch_store::Macro_bunch_store(numeric::array& numeric_local_particles,
                                     int local_num, int total_num,
                                     double mass, int charge, double total_current, double bunch_np,
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
    this->bunch_np = bunch_np;
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
            double gammai=gamma - local_particles(5, i);
            double rcp_gammai = 1.0 /gammai;
            double sum1=1.0 + xp * xp + yp * yp;
            double betazi;
            if ((1.0 - rcp_gammai * rcp_gammai * sum1)<0.0){  
                std::cout<<"beta i ^2=" <<1.0 - rcp_gammai * rcp_gammai * (1.0 + xp * xp + yp * yp)<<std::endl;
                std::cout<<"rcp_gammai="<<rcp_gammai<<" delta gamma_i="<<local_particles(5, i)<<" gamma="<< gamma<<std::endl;
                std::cout<<"xp="<<xp<<"   yp="<<yp<<std::endl;
                throw std::runtime_error (" error in convert_to_fixedt,"
                        "probably the approximation beta ~betaz broke down, unstable beam ");
                abort();        
            } 
            else {
                betazi = sqrt(1.0 - rcp_gammai * rcp_gammai * sum1);
            }
            local_particles(4, i) *= -gamma * betazi / units(0); // units(0)
            // is not
            // an error!
            if (gammai*gammai-sum1<0.0){
                std::cout<<"gammai*gammai-sum1="<<gammai*gammai-sum1<<std::endl;
                std::cout<<"gammai="<<gammai<<std::endl;
                std::cout<<"xp="<<xp<<"   yp="<<yp<<std::endl;
                throw std::runtime_error (" error in convert_to_fixedt,"
                        "gamma*gamma-1-px^2-py^2 cannot be negative");
                abort();        
            }
            local_particles(5, i)=sqrt(gammai*gammai-sum1);

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
            double zp = local_particles(5, i);
            double gammai=sqrt(1.0 + xp * xp + yp * yp+zp * zp);            
            double betazi=zp/gammai;
            local_particles(4, i) /= -gamma * betazi / units(0); // units(0)
            // is not
            // an error!
            local_particles(5, i) =gamma-gammai;
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
