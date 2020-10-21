#ifndef SPACE_CHARGE_2D_KV_H_
#define SPACE_CHARGE_2D_KV_H_

#include "synergia/simulation/operator.h"
#include "synergia/simulation/collective_operator_options.h"

class Space_charge_2d_kv;

struct Space_charge_2d_kv_options 
    : public CO_base_options<Space_charge_2d_kv_options, Space_charge_2d_kv>
{
    enum class LongitudinalDistribution
    {
        gaussian,
        uniform,
    };

    using LD = LongitudinalDistribution;

    // switch to control whether the linear charge density is assumed to be
    // distributed gaussian or uniform over the bunch length.
    LD longitudinal_distribution = LD::uniform;

    bool strictly_linear = true;
    bool strictly_centered = false;

    template<class Archive>
    void serialize(Archive & ar)
    { 
        ar(cereal::base_class<CO_base_options>(this)); 
        ar(longitudinal_distribution);
        ar(strictly_linear);
        ar(strictly_centered);
    }
};

CEREAL_REGISTER_TYPE(Space_charge_2d_kv_options)


class Space_charge_2d_kv : public Collective_operator
{

private:

    const Space_charge_2d_kv_options opts;

    void apply_impl(
            Bunch_simulator& simulator, 
            double time_step, 
            Logger& logger) override;

    void apply_bunch( 
            Bunch& bunch, 
            double time_step, 
            Logger& logger);

public:

    Space_charge_2d_kv(Space_charge_2d_kv_options const& opts);
};

#endif /* SPACE_CHARGE_2D_KV_H_ */
