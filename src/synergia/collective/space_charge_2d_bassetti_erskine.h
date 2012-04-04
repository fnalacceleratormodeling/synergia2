#ifndef SPACE_CHARGE_2D_BASSETTI_ERSKINE_H_
#define SPACE_CHARGE_2D_BASSETTI_ERSKINE_H_

#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"

#include "synergia/utils/serialization.h"

class Space_charge_2d_bassetti_erskine : public Collective_operator
{
private:
    double sigma[3];
    bool use_round;
    // By default = 1
    // If 1: then round beam approximation
    // used when horizontal and vertical
    // sigmas approximately equal.

public:
    Space_charge_2d_bassetti_erskine();
    // pointer to an array containing sigma_x, sigma_y and sigma_cdt[m]
    void
    set_sigma(double* = 0);
    // returns the "normalized" electric field in the rest frame of the bunch,
    // in inverse meters.  To get the field [V/m], this must be multiplied
    // by Q/(2 pi epsilon_o), where Q is the line density of charge [C/m]
    // (in rest frame).
    std::vector<double >
    normalized_efield(double x, double y);
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
    virtual
    ~Space_charge_2d_bassetti_erskine();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
	    ar & BOOST_SERIALIZATION_NVP(sigma)
		& BOOST_SERIALIZATION_NVP(use_round);
        }
};
BOOST_CLASS_EXPORT_KEY(Space_charge_2d_bassetti_erskine)

typedef boost::shared_ptr<Space_charge_2d_bassetti_erskine >
        Space_charge_2d_bassetti_erskine_sptr;

#endif /* SPACE_CHARGE_2D_BASSETTI_ERSKINE_H_ */
