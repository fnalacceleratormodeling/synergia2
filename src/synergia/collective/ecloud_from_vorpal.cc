#include "ecloud_from_vorpal.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/simple_timer.h"
#
Ecloud_from_vorpal::Ecloud_from_vorpal():
file_name_archive("none"),
field_name("blank")
{
;
}
Ecloud_from_vorpal::Ecloud_from_vorpal(Commxx_sptr comm_sptr_in, 
                                       const std::string &f_name_archive, const std::string aDeviceName):
file_name_archive(f_name_archive),
comm_sptr(comm_sptr_in),
e_field(comm_sptr_in, f_name_archive.c_str()),
field_name("")
{
  field_name += e_field.getVORPALJobName();
  subjectedDevices.push_back(aDeviceName);
}

Ecloud_from_vorpal::~Ecloud_from_vorpal() { ; }

void Ecloud_from_vorpal::apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger) {
//
// Find if we have an electron cloud in this device. 
//
// Actually, we are already in a device where the e-cloud is present. 
// 
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    double beta = bunch.get_reference_particle().get_beta();
    double lStep = step.get_length();
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);
    double p_ref = bunch.get_reference_particle().get_momentum();
    
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        const double x = bunch.get_local_particles()[part][Bunch::x];
        const double y = bunch.get_local_particles()[part][Bunch::y];
        const double dz = bunch.get_local_particles()[part][Bunch::z];
	double px = p_ref*bunch.get_local_particles()[part][Bunch::xp]; // in GeV/c
	const double delta_px = e_field.GetFieldEX(x, y, dz) * q * lStep / (beta*pconstants::c); // in SI units. 
	px += unit_conversion*delta_px; // adding in GeV/c 
	double py = p_ref*bunch.get_local_particles()[part][Bunch::yp]; // in GeV/c
	const double delta_py = e_field.GetFieldEX(x, y, dz) * q * lStep / (beta*pconstants::c); // in SI units. 
	py += unit_conversion*delta_py; // adding in GeV/c 
	bunch.get_local_particles()[part][Bunch::xp] = px/p_ref;
	bunch.get_local_particles()[part][Bunch::yp] = py/p_ref;
   }  	    
	    
	    
}
