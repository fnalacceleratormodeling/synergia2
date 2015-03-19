#ifndef FF_MULTIPOLE_H
#define FF_MULTIPOLE_H

#include "ff_element.h"

class FF_multipole : public FF_element
{
public:
    FF_multipole();

    template <typename T>
    inline static void multipole_unit(T & x, T & xp,
                                      T & y, T & yp,
                                      T & cdt, T const& dpop,
                                      double length, 
                                      double strength, 
                                      double reference_momentum,
                                      double m, 
                                      double reference_brho );

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_multipole();
};

typedef boost::shared_ptr<FF_multipole > FF_multipole_sptr;

template <typename T>
inline void FF_multipole::multipole_unit(T & x, T & xp,
                                         T & y, T & yp,
                                         T & cdt, T const& dpop,
                                         double length, 
                                         double strength,
                                         double reference_momentum,
                                         double m, 
                                         double reference_brho )
{
}

#endif // FF_MULTIPOLE_H
