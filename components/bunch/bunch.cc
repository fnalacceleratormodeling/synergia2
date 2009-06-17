#include "bunch.h"

//	Reference_particle reference_particle;
//	int particle_charge;
//	Array_2d<double> local_particles;
//	int local_num, total_num;
//	double total_current;
//	Bunch_state state;
//	bool particles_valid;
//	bool total_num_valid;


Bunch::Bunch(Reference_particle& reference_particle, int particle_charge,
		int total_num, double real_num):
			reference_particle(reference_particle)
{
}

//	void set_particle_charge(int particle_charge);
//	void set_local_num(int local_num);
//	void set_total_num(int total_num, bool update_current=true);
//	void set_total_current(double total_current);
//
//	Array_2d<double>& get_local_particles();
//	int get_particle_charge();
//	int get_local_num();
//	int get_total_num();
//	double get_total_current();
//	Bunch_state get_state();
//
//	int convert_to_state(Bunch_state state);
//	void update_total_num();
//
//	void inject(Bunch bunch);
