#include "ecloud_from_vorpal.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/commxx.h"

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
// Find if we have an electron cloud in this device. Done above!! 

   bunch.convert_to_state(Bunch::fixed_z_lab); 
//
// Actually, we are already in a device where the e-cloud is present. 
// 
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    double beta = bunch.get_reference_particle().get_beta();
    double lStep = step.get_length();
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e); // convert a momentum from SI to GeV/c
    double fact = q * lStep / (beta*pconstants::c) ; // dp = d(mv) = F * dt = F *lStep/(beta c) = q*E*lStep/(beta c) = fact * E
     // So, fact = q*lStep/(beta*c).  In SI units. 
    double p_ref = bunch.get_reference_particle().get_momentum();
    Commxx myComm; // The communicator 
    int my_rank= myComm.get_rank();
    int partDump = 3*bunch.get_local_num()/4;
    if ((my_rank == 0) && (verbosity > 3) ) std::cerr << " Ecloud_from_vorpal::apply p_ref = " << p_ref << " lStep " << lStep 
                                 << " kinetics at particle "  << partDump << " verbosity " << verbosity << std::endl;
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        const double x = bunch.get_local_particles()[part][Bunch::x];
        const double y = bunch.get_local_particles()[part][Bunch::y];
        const double dz = bunch.get_local_particles()[part][Bunch::z];
	double px = p_ref*bunch.get_local_particles()[part][Bunch::xp]; // in GeV/c
	double py = p_ref*bunch.get_local_particles()[part][Bunch::yp]; // in GeV/c
	if ((part == partDump) && (my_rank == 0) && (verbosity > 3)) 
	  std::cerr << " At x= " << x << " y " << y << " dz " << dz << " before kick px " << px << " py " << py << std::endl;
	const double delta_px = e_field.GetFieldEX(x, y, dz) * fact; // in SI units. 
	px += unit_conversion*delta_px; // adding in GeV/c 
	const double delta_py = e_field.GetFieldEX(x, y, dz) * fact; // in SI units. 
	py += unit_conversion*delta_py; // adding in GeV/c 
	bunch.get_local_particles()[part][Bunch::xp] = px/p_ref;
	bunch.get_local_particles()[part][Bunch::yp] = py/p_ref;
	if ((part == partDump) && (my_rank == 0) && (verbosity > 3)) 
	  std::cerr << " Ex= " << e_field.GetFieldEX(x, y, dz) << " y " << e_field.GetFieldEY(x, y, dz) 
	            << " after kick px " << px << " py " << py << std::endl;
   }  	    
	    
	    
}
