#ifndef PHYSICAL_CONSTANTS_H_
#define PHYSICAL_CONSTANTS_H_

namespace pconstants
{
    const double mp = 0.93827203; // Mass of proton [GeV/c^2]
    const double me = 0.51099892e-3; // Mass of electron [GeV/c^2]
    const double mmu = 0.105658369; // Mass of muon [GeV/c^2]

    const double e = 1.6021892e-19; // Charge of proton [C]

    const int proton_charge = 1; // Charge in units of e
    const int antiproton_charge = -1; // Charge in units of e
    const int electron_charge = -1; // Charge in units of e
    const int positron_charge = 1; // Charge in units of e
    const int muon_charge = -1; // Charge in units of e
    const int antimuon_charge = 1; // Charge in units of e
}

#endif /* PHYSICAL_CONSTANTS_H_ */
