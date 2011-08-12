#include "lattice_simulator.h"
#include <stdexcept>

Lattice_functions::Lattice_functions() :
    alpha_x(0.0), alpha_y(0.0), beta_x(0.0), beta_y(0.0), psi_x(0.0),
            psi_y(0.0), D_x(0.0), D_y(0.0), Dprime_x(0.0), Dprime_y(0.0)
{
}

Lattice_functions::Lattice_functions(LattFuncSage::lattFunc const& latt_func) :
    alpha_x(latt_func.alpha.hor), alpha_y(latt_func.alpha.ver),
            beta_x(latt_func.beta.hor), beta_y(latt_func.beta.ver),
            psi_x(latt_func.psi.hor), psi_y(latt_func.psi.ver),
            D_x(latt_func.dispersion.hor), D_y(latt_func.dispersion.ver),
            Dprime_x(latt_func.dPrime.hor), Dprime_y(latt_func.dPrime.ver)
{
}

void
Lattice_simulator::construct_extractor_map()
{
    Operation_extractor_sptr chef_mixed_operation_extractor(
            new Chef_mixed_operation_extractor(chef_lattice_sptr, map_order));

    extractor_map_sptr->set_extractor(default_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map_sptr->set_extractor(chef_mixed_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map_sptr->set_extractor(
            chef_propagate_operation_extractor_name,
            Operation_extractor_sptr(
                    new Chef_propagate_operation_extractor(chef_lattice_sptr,
                            map_order)));
    extractor_map_sptr->set_extractor(
            chef_map_operation_extractor_name,
            Operation_extractor_sptr(
                    new Chef_map_operation_extractor(chef_lattice_sptr,
                            map_order)));
}

Lattice_simulator::Lattice_simulator(Lattice_sptr lattice_sptr, int map_order) :
    lattice_sptr(lattice_sptr),
            chef_lattice_sptr(new Chef_lattice(lattice_sptr)),
            extractor_map_sptr(new Operation_extractor_map),
            map_order(map_order)
{
    construct_extractor_map();
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
Lattice_simulator::calculate_lattice_functions()
{
    LattFuncSage latt_func_sage(chef_lattice_sptr->get_sliced_beamline_sptr());
    const int map_order = 1;
    JetParticle jet_particle = reference_particle_to_chef_jet_particle(
            lattice_sptr->get_reference_particle(), map_order);
    int cslf_retval =
            latt_func_sage.CourantSnyderLatticeFunctions(jet_particle);
    if (cslf_retval != LattFuncSage::OK) {
        std::string
                message(
                        "Lattice_simulator::calculate_lattice_functions failed with message: ");
        if (cslf_retval == LattFuncSage::SLOTS_DETECTED) {
            message += "slots detected";
        } else if (cslf_retval == LattFuncSage::UNSTABLE) {
            message += "unstable";
        } else if (cslf_retval == LattFuncSage::INTEGER_TUNE) {
            message += "integer tune";
        } else if (cslf_retval == LattFuncSage::PHASE_ERROR) {
            message += "phase error";
        } else if (cslf_retval == LattFuncSage::WRONG_COUNT) {
            message += "wrong count";
        } else if (cslf_retval == LattFuncSage::NOT_WRITTEN) {
            message += "not written";
        } else if (cslf_retval == LattFuncSage::TOO_MANY_VECTORS) {
            message += "too many vectors";
        } else {
            message += "(unknown message)";
        }
        throw std::runtime_error(message);
    }
    std::vector<LattFuncSage::lattFunc > latt_func(
            latt_func_sage.getTwissArray());
    for (int i = 0; i < latt_func.size(); ++i) {
        Lattice_functions latttice_functions(latt_func.at(i));
        std::cout << latt_func.at(i).arcLength << " "
                << latttice_functions.beta_x << std::endl;
    }
}

Lattice_simulator::~Lattice_simulator()
{
}

