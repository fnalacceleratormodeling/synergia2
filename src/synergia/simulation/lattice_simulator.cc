#include "lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"

void
Lattice_simulator::construct_extractor_map()
{
    Operation_extractor_sptr chef_mixed_operation_extractor(
            new Chef_mixed_operation_extractor(chef_lattice_sptr, map_order));

    extractor_map_sptr->set_extractor(default_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map_sptr->set_extractor(chef_mixed_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map_sptr->set_extractor(chef_propagate_operation_extractor_name,
            Operation_extractor_sptr(new Chef_propagate_operation_extractor(
                    chef_lattice_sptr, map_order)));
    extractor_map_sptr->set_extractor(chef_map_operation_extractor_name,
            Operation_extractor_sptr(new Chef_map_operation_extractor(
                    chef_lattice_sptr, map_order)));
}

Lattice_simulator::Lattice_simulator(Lattice_sptr lattice_sptr,
        int map_order) :
    lattice_sptr(lattice_sptr), chef_lattice_sptr(new Chef_lattice(
            lattice_sptr)), extractor_map_sptr(new Operation_extractor_map),
            map_order(map_order)
{
    construct_extractor_map();
    set_bucket_length(); 
}

void
Lattice_simulator::construct_sliced_chef_beamline(
        Lattice_element_slices const& slices)
{
    chef_lattice_sptr->construct_sliced_beamline(slices);
}

int
Lattice_simulator::get_map_order() const
{
    return map_order;
}

Operation_extractor_map_sptr
Lattice_simulator::get_operation_extractor_map_sptr()
{
    return extractor_map_sptr;
}

Lattice_sptr
Lattice_simulator::get_lattice_sptr()
{
    return lattice_sptr;
}

Chef_lattice_sptr
Lattice_simulator::get_chef_lattice_sptr()
{
    return chef_lattice_sptr;
}

void
Lattice_simulator::set_bucket_length()
{
    double freq(0.), freq2(0.);
    int isw=0;
    double eps=1e-6;
    for (Lattice_elements::const_iterator it= this->lattice_sptr->get_elements().begin();
           it != this->lattice_sptr->get_elements().end(); ++it){
        
        if ((*it)->has_double_attribute("freq")) {
                 freq=(*it)->get_double_attribute("freq");
                 if ((isw==1) && (fabs(freq-freq2) > eps)) {
                    throw std::runtime_error(
                         "set_bucket_length: rf elements with different frequencies found!!");
                  }
                 freq2=freq;
                 isw=1;
        }
        if (isw==1) {
        double beta=this->get_lattice_sptr()->get_reference_particle().get_beta();
                 this->bucket_length=pconstants::c*beta/freq;
                 }
        else {
                this->bucket_length=0.0;
        }
    } 
}

double 
Lattice_simulator::get_bucket_length()
{
    return this->bucket_length;
}

int
Lattice_simulator::get_number_buckets()
{
double eps=1e-5;
int number_buckets;
double bl=get_bucket_length();
double ol=this->get_lattice_sptr()->get_length();
bl>eps ? number_buckets=int(ol/bl):number_buckets=1;
return number_buckets;
}



Lattice_simulator::~Lattice_simulator()
{
}

