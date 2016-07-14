#ifndef FF_RFCAVITY_H
#define FF_RFCAVITY_H

#include "ff_element.h"

class FF_rfcavity : public FF_element
{
public:
    FF_rfcavity();

    template <typename T>
    inline static void rfcavity_unit(T & x, T & xp,
                                      T & y, T & yp,
                                      T & cdt, T const& dpop,
                                      double length, 
                                      double freq, 
                                      double reference_momentum,
                                      double m, 
                                      double reference_brho );

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_rfcavity();
};

typedef boost::shared_ptr<FF_rfcavity > FF_rfcavity_sptr;

template <typename T>
inline void FF_rfcavity::rfcavity_unit(T & x, T & xp,
                                         T & y, T & yp,
                                         T & cdt, T const& dpop,
                                         double length, 
                                         double freq,
                                         double reference_momentum,
                                         double m, 
                                         double reference_brho )
{
}

#endif // FF_RFCAVITY_H
