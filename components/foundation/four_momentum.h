#ifndef FOUR_MOMENTUM_H_
#define FOUR_MOMENTUM_H_

///
/// Provide conversion between various kinematic parameters
///

class Four_momentum {
private:
	double mass,energy,momentum,gamma,beta;
	void update_from_gamma();
public:
    /// Construct a Four_momentum in the rest frame
    /// \param mass in :math: GeV/c^2
	Four_momentum(double mass); 
	void set_total_energy(double total_energy); ///< total energy in GeV
	void set_kinetic_energy(double kinetic_energy); ///< kinetic energy in GeV
	void set_momentum(double momentum); ///< momentum in GeV/c
	void set_gamma(double gamma); ///< relativistic gamma
	void set_beta(double beta); ///< relativistic beta

	double get_mass(); ///
	double get_total_energy();
	double get_kinetic_energy();
	double get_momentum();
	double get_gamma();
	double get_beta();
};

#endif /* FOUR_MOMENTUM_H_ */
