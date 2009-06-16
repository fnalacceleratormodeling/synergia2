#ifndef REFERENCE_PARTICLE_H_
#define REFERENCE_PARTICLE_H_

class Reference_particle {
private:
	double energy;
	double state[6];
public:
	Reference_particle(double energy);
	Reference_particle(double energy, double state[6]);

	void set_energy(double energy);
	void set_state(double state[6]);

	double get_energy();
	double * get_state();
};

#endif /* REFERENCE_PARTICLE_H_ */
