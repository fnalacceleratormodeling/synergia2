#ifndef SPACE_CHARGE_2D_BASSETTI_ERSKINE_H_
#define SPACE_CHARGE_2D_BASSETTI_ERSKINE_H_

#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"

#include "synergia/utils/serialization.h"

class Space_charge_2d_bassetti_erskine : public Collective_operator
{
private:
    double sigma_x, sigma_y, sigma_cdt;
    bool is_round;
public:
    Space_charge_2d_bassetti_erskine();
    virtual Space_charge_2d_bassetti_erskine *
    clone();
    void
    set_sigma(double sigma_x, double sigma_y, double sigma_cdt);
    // returns the "normalized" electric field in the rest frame of the bunch,
    // in inverse meters.  To get the field [V/m], this must be multiplied
    // by Q/(2 pi epsilon_o), where Q is the line density of charge [C/m]
    // (in rest frame).
    std::vector<double >
    normalized_efield(double x, double y);
    void
    normalized_efield(double x, double y, double & E_x, double & E_y);
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
    virtual
    ~Space_charge_2d_bassetti_erskine();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};
BOOST_CLASS_EXPORT_KEY(Space_charge_2d_bassetti_erskine)

typedef boost::shared_ptr<Space_charge_2d_bassetti_erskine >
        Space_charge_2d_bassetti_erskine_sptr; // syndoc:include

#endif /* SPACE_CHARGE_2D_BASSETTI_ERSKINE_H_ */
