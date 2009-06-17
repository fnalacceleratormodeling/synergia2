#ifndef REFERENCE_PARTICLE_H_
#define REFERENCE_PARTICLE_H_

class Reference_particle {
private:
	double total_energy;
	double state[6];
	double units[6];
public:
	Reference_particle(double total_energy, double units[6]);
	Reference_particle(double total_energy, double units[6], double state[6]);

	void set_total_energy(double total_energy);
	void set_state(double state[6]);

	double get_total_energy();
	double * get_units();
	double * get_state();
};

#endif /* REFERENCE_PARTICLE_H_ */
