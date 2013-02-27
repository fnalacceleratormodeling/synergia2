#include "ecloud_from_vorpal.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/commxx.h"

Ecloud_from_vorpal::Ecloud_from_vorpal():
file_name_archive("none"),
field_name("blank"),
enhanceFactor(1.0)
{
;
}
Ecloud_from_vorpal::Ecloud_from_vorpal(Commxx_sptr comm_sptr_in, 
                                       const std::string &f_name_archive, const std::string aDeviceName):
file_name_archive(f_name_archive),
comm_sptr(comm_sptr_in),
e_field(comm_sptr_in, f_name_archive.c_str()),
field_name(""),
enhanceFactor(1.0)
{
  field_name += e_field.getVORPALJobName();
  subjectedDevices.push_back(aDeviceName);
//  std::cerr << "   Ecloud_from_vorpal::Ecloud_from_vorpal, defined, for device type " << aDeviceName << std::endl;
  
} 

Ecloud_from_vorpal::~Ecloud_from_vorpal() { ; }

Ecloud_from_vorpal *
Ecloud_from_vorpal::clone()
{
    return new Ecloud_from_vorpal(*this);
}

void Ecloud_from_vorpal::apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger) {
//
// Actually, we are already in a device where the e-cloud is present.
// Not really ! The clever element stepper is not yet implemented yet. 
//

    Commxx myComm; // The communicator 
    int my_rank= myComm.get_rank();
   if (!this->checkElementType(step)) {
      if ((my_rank == 0) && (verbosity > 3) ) std::cerr << " Rank 0; Ecloud_from_vorpal, element with no e-cloud , skip " << std::endl; 
      return;
   } 
//
// Find if we have an electron cloud in this device. Done above!!  But we learn how to check where we are.. 
//   
   bunch.convert_to_state(Bunch::fixed_t_lab); 
//
//    this->getElementBoudaries(step); 
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    double beta = bunch.get_reference_particle().get_beta();
    double lStep = step.get_length();
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e); // convert a momentum from SI to GeV/c
    double fact = q * lStep / (beta*pconstants::c) ; // dp = d(mv) = F * dt = F *lStep/(beta c) = q*E*lStep/(beta c) = fact * E
     // So, fact = q*lStep/(beta*c).  In SI units. 
    double p_ref = bunch.get_reference_particle().get_momentum();
    int partDump = 3*bunch.get_local_num()/4;
    if ((my_rank == 0) && (verbosity > 3) ) std::cerr << " Rank 0; Ecloud_from_vorpal::apply p_ref = " << p_ref << " lStep " << lStep 
                    << "  Num Particle " <<   bunch.get_local_num() << " kinetics at particle "  << partDump << " verbosity " << verbosity << std::endl;
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        const double x = bunch.get_local_particles()[part][Bunch::x];
        const double y = bunch.get_local_particles()[part][Bunch::y];
        const double dz = bunch.get_local_particles()[part][Bunch::z];
	double px = p_ref*bunch.get_local_particles()[part][Bunch::xp]; // in GeV/c
	double py = p_ref*bunch.get_local_particles()[part][Bunch::yp]; // in GeV/c
	if ((part == partDump) && (my_rank == 0) && (verbosity > 3)) 
	  std::cerr << " At x= " << x << " y " << y << " dz " << dz << " before kick px " << px << " py " << py << std::endl;
	const double delta_px = e_field.GetFieldEX(x, y, dz) * fact; // in SI units. 
	px += enhanceFactor*unit_conversion*delta_px; // adding in GeV/c 
	const double delta_py = e_field.GetFieldEX(x, y, dz) * fact; // in SI units. 
	py += enhanceFactor*unit_conversion*delta_py; // adding in GeV/c 
	bunch.get_local_particles()[part][Bunch::xp] = px/p_ref;
	bunch.get_local_particles()[part][Bunch::yp] = py/p_ref;
	if ((part == partDump) && (my_rank == 0) && (verbosity > 3)) 
	  std::cerr << " Ex= " << e_field.GetFieldEX(x, y, dz) << " y " << e_field.GetFieldEY(x, y, dz) 
	            << " after kick px " << px << " py " << py << std::endl;
   }  	    
	    
	    
}
template<class Archive>
        void
        Ecloud_from_vorpal::save(Archive & ar, const unsigned int version) const
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
            ar & BOOST_SERIALIZATION_NVP(comm_sptr)
	       & BOOST_SERIALIZATION_NVP(file_name_archive)
	       & BOOST_SERIALIZATION_NVP(field_name);
	       // We do not archive the field map, already archived. And it can't change. 
        }
template<class Archive>
        void
        Ecloud_from_vorpal::load(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
            ar & BOOST_SERIALIZATION_NVP(comm_sptr)
	       & BOOST_SERIALIZATION_NVP(file_name_archive)
	       & BOOST_SERIALIZATION_NVP(field_name);
	       std::cerr << " Reloading e_field...from file  " <<file_name_archive << std::endl;
	       e_field = ECloudEFieldVORPAL2D(comm_sptr, file_name_archive.c_str()); 
	       std::cerr << " Reloading e_field.Ex at relevant pts " 
	                 <<  e_field.GetFieldEX(-0.00877598, -0.00401517, -0.489653) << std::endl;
//		std::cerr << " And quit for now ... " << std::endl;
//		MPI_Abort(MPI_COMM_WORLD, 111);	 exit(2);
        }
	
	
void Ecloud_from_vorpal::getElementBoudaries( const Step & step ) const {

    Operators operators(step.get_operators());
    Commxx myComm; // The communicator 
    int my_rank= myComm.get_rank();
    if (my_rank == 0) std::cerr << " Studying ElementBoudaries for step of length " << step.get_length() << std::endl;
    for (Operators::iterator oit = operators.begin(); oit!= operators.end(); ++oit) {
        if((*oit)->get_type() != "independent") continue;
	if (my_rank == 0) std::cerr << " .....  got an independent operator .. " << std::endl;
        Lattice_element_slices
                    slices(
                            boost::static_pointer_cast<Independent_operator >(
                                    *oit)->get_slices());
	for (Lattice_element_slices::iterator slit = slices.begin(); slit
                    != slices.end(); ++slit) {
                double l = (*slit)->get_right() - (*slit)->get_left();
		if (my_rank == 0) {
		   if ((*slit)->is_whole()) std::cerr << " Whole slice, left  ";
		   else std::cerr << " fractional slice, left  ";
		   std::cerr << (*slit)->get_left() 
		                            << " right " << (*slit)->get_right() << " length " << l << std::endl;
                   Lattice_element tmpElem = (*slit)->get_lattice_element();
                   std::cerr << " .............. Element type " 
		             << (*slit)->get_lattice_element().get_type() << " name " <<  (*slit)->get_lattice_element().get_name() << std::endl;
		   
		}	
	}		    
   }

}

bool Ecloud_from_vorpal::checkElementType( const Step & step ) const {

    Operators operators(step.get_operators());
    Commxx myComm; // The communicator 
    int my_rank= myComm.get_rank();
    for (Operators::iterator oit = operators.begin(); oit!= operators.end(); ++oit) {
        if((*oit)->get_type() != "independent") continue;
        Lattice_element_slices
                    slices(
                            boost::static_pointer_cast<Independent_operator >(
                                    *oit)->get_slices());
	for (Lattice_element_slices::iterator slit = slices.begin(); slit
                    != slices.end(); ++slit) {
                std::string type = (*slit)->get_lattice_element().get_type();
		for(std::vector<std::string>::const_iterator its=subjectedDevices.begin(); its != subjectedDevices.end(); its++) {
		  if (type.find((*its).c_str()) != std::string::npos) return true;
		}
	}	
    }
    return false;
}


template
void
Ecloud_from_vorpal::save<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version) const;
template
void
Ecloud_from_vorpal::save<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version) const;

template
void
Ecloud_from_vorpal::load<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);
template
void
Ecloud_from_vorpal::load<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);


BOOST_CLASS_EXPORT_IMPLEMENT(Ecloud_from_vorpal);
