#ifndef PHYSICAL_CONSTANTS_H_
#define PHYSICAL_CONSTANTS_H_
#include "synergia/foundation/math_constants.h"
#include <string>

#define PDG 2010 // use 2010 pdg values for now

namespace pconstants
{
#if PDG == 2012
    const std::string pdg_year("2012");
//  J. Beringer et al. (Particle Data Group), PR D86, 010001 (2012) and 2013 partial
//     update for the 2014 edition (URL: http://pdg.lbl.gov)
    const double mp = 0.938272046; // Mass of proton [GeV/c^2]
    const double proton_mass = 0.938272046; // Mass of proton [GeV/c^2]
    const double me = 0.510998928e-3; // Mass of electron [GeV/c^2]
    const double electron_mass = 0.510998928e-3; // Mass of electron [GeV/c^2]
    const double mmu = 0.1056583715; // Mass of muon [GeV/c^2]
    const double muon_mass = 0.1056583715; // Mass of muon [GeV/c^2]
    const double e = 1.602176565e-19; // Charge of proton [C]
#else 
#if PDG == 2010
//  K. Nakamura et al. (Particle Data Group), J. Phys. G 37, 075021 (2010) 
    const std::string pdg_year("2010");
    const double mp = 0.938272013; // Mass of proton [GeV/c^2]
    const double proton_mass = 0.938272013; // Mass of proton [GeV/c^2]
    const double me = 0.510998910e-3; // Mass of electron [GeV/c^2]
    const double electron_mass = 0.510998910e-3; // Mass of electron [GeV/c^2]
    const double mmu = 0.105658367; // Mass of muon [GeV/c^2]
    const double muon_mass = 0.105658367; // Mass of muon [GeV/c^2]
    const double e = 1.602176487e-19; // Charge of proton [C]
#else
#if PDG == 2008 // PDG2008
//  C. Amsler et al. (Particle Data Group), Physics Letters B667, 1 (2008) 
    const std::string pdg_year("2008");
    const double mp = 0.938272013; // Mass of proton [GeV/c^2]
    const double proton_mass = 0.938272013; // Mass of proton [GeV/c^2]
    const double me = 0.510998910e-3; // Mass of electron [GeV/c^2]
    const double electron_mass = 0.510998910e-3; // Mass of electron [GeV/c^2]
    const double mmu = 0.105658367; // Mass of muon [GeV/c^2]
    const double muon_mass = 0.105658367; // Mass of muon [GeV/c^2]
    const double e = 1.602176487e-19; // Charge of proton [C]
#if PDG == 2006 // PDG2006
//  S. Eidelman et al. (Particle Data Group), Phys. Lett. B 592, 1 (2004)
    const std::string pdg_year("2006");
    const double mp = 0.938272029; // Mass of proton [GeV/c^2]
    const double proton_mass = 0.938272029; // Mass of proton [GeV/c^2]
    const double me = 0.510998918e-3; // Mass of electron [GeV/c^2]
    const double electron_mass = 0.510998918e-3; // Mass of electron [GeV/c^2]
    const double mmu = 0.105658369; // Mass of muon [GeV/c^2]
    const double muon_mass = 0.105658369; // Mass of muon [GeV/c^2]
    const double e = 1.60217653e-19; // Charge of proton [C]
#else
#if PDG == 2004
//  S. Eidelman et al. (Particle Data Group), Phys. Lett. B 592, 1 (2004)
    const std::string pdg_year("2004");
    const double mp = 0.938272029; // Mass of proton [GeV/c^2]
    const double proton_mass = 0.938272029; // Mass of proton [GeV/c^2]
    const double me = 0.510998918e-3; // Mass of electron [GeV/c^2]
    const double electron_mass = 0.510998918e-3; // Mass of electron [GeV/c^2]
    const double mmu = 0.105658369; // Mass of muon [GeV/c^2]
    const double muon_mass = 0.105658369; // Mass of muon [GeV/c^2]
    const double e = 1.60217653e-19; // Charge of proton [C]
#else
#if PDG == -1 // legacy
    const std::string pdg_year("legacy");
    const double mp = 0.93827203; // Mass of proton [GeV/c^2]
    const double proton_mass = 0.93827203; // Mass of proton [GeV/c^2]
    const double me = 0.51099892e-3; // Mass of electron [GeV/c^2]
    const double electron_mass = 0.51099892e-3; // Mass of electron [GeV/c^2]
    const double mmu = 0.105658369; // Mass of muon [GeV/c^2]
    const double muon_mass = 0.105658369; // Mass of muon [GeV/c^2]
    const double e = 1.6021892e-19; // Charge of proton [C]
#else
#error "No selection for physicsl_constants PDG set"
#endif // legacy
#endif // PDG 2004
#endif // PDG 2006
#endif // PDG 2008
#endif // PDG 2010
#endif // PDG 2012
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
