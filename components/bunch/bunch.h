#ifndef BUNCH_H_
#define BUNCH_H_

#include "array_nd/array_2d.h"
#include "components/foundation/reference_particle.h"

enum Bunch_state {fixed_z = 1, fixed_t = 2};

class Bunch {
private:
	Reference_particle reference_particle;
	int particle_charge;
	Array_2d<double> local_particles;
	int local_num, total_num;
	double real_num;
	Bunch_state state;
	bool particles_valid;
	bool total_num_valid;
public:
	Bunch(Reference_particle& reference_particle, int particle_charge,
			int total_num, double real_num);
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
};

#endif /* BUNCH_H_ */
