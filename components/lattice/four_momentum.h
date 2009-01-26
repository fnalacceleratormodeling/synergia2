#ifndef FOUR_MOMENTUM_H_
#define FOUR_MOMENTUM_H_

class Four_momentum {
private:
	double mass,energy,momentum,gamma,beta;
	void update_from_gamma();
public:
	Four_momentum(double mass, double total_energy=0.0);
	void set_total_energy(double total_energy);
	void set_kinetic_energy(double kinetic_energy);
	void set_momentum(double momentum);
	void set_gamma(double gamma);
	void set_beta(double beta);

	double get_mass();
	double get_total_energy();
	double get_kinetic_energy();
	double get_momentum();
	double get_gamma();
	double get_beta();
};

#endif /* FOUR_MOMENTUM_H_ */
