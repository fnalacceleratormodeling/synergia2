#ifndef PHYSICAL_CONSTANTS_H_
#define PHYSICAL_CONSTANTS_H_
#include "synergia/foundation/math_constants.h"

namespace pconstants
{
    const double mp = 0.93827203; // Mass of proton [GeV/c^2]
    const double me = 0.51099892e-3; // Mass of electron [GeV/c^2]
    const double mmu = 0.105658369; // Mass of muon [GeV/c^2]

    const double e = 1.6021892e-19; // Charge of proton [C]

    const double c = 299792458.0; // Speed of light [m/s]

    const double mu0 = 4*mconstants::pi*1.0e-7; // Permittivity of free space [H/m]
    const double epsilon0 = 1.0/(c*c*mu0); // Permeability of free space [F/m]

    // Classical radius of a particle = e^2/(4 pi epsilon0 m c^2)
    const double re = e/(4*mconstants::pi*epsilon0*me*1.0e9); // Classical
                                                   // radius of electron [m]
    const double rp = e/(4*mconstants::pi*epsilon0*mp*1.0e9); // Classical
                                                   // radius of proton [m]
    const double rmu = e/(4*mconstants::pi*epsilon0*mmu*1.0e9); // Classical
                                                   // radius of muon [m]

    const int proton_charge = 1; // Charge in units of e
    const int antiproton_charge = -1; // Charge in units of e
    const int electron_charge = -1; // Charge in units of e
    const int positron_charge = 1; // Charge in units of e
    const int muon_charge = -1; // Charge in units of e
    const int antimuon_charge = 1; // Charge in units of e
}

#endif /* PHYSICAL_CONSTANTS_H_ */
