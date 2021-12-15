var data = {lines:[
{"lineNum":"    1","line":"#ifndef PHYSICAL_CONSTANTS_H_"},
{"lineNum":"    2","line":"#define PHYSICAL_CONSTANTS_H_"},
{"lineNum":"    3","line":"//#include <basic_toolkit/PhysicsConstants.h>"},
{"lineNum":"    4","line":"#include \"synergia/foundation/math_constants.h\""},
{"lineNum":"    5","line":"#include <string>"},
{"lineNum":"    6","line":""},
{"lineNum":"    7","line":"#ifndef PDG_VERSION"},
{"lineNum":"    8","line":"#define PDG_VERSION 2012"},
{"lineNum":"    9","line":"#endif"},
{"lineNum":"   10","line":""},
{"lineNum":"   11","line":"namespace pconstants"},
{"lineNum":"   12","line":"{"},
{"lineNum":"   13","line":"#if PDG_VERSION == 2012"},
{"lineNum":"   14","line":"    const std::string pdg_year(\"2012\");","class":"lineCov","hits":"1","order":"706","possible_hits":"1",},
{"lineNum":"   15","line":"//  J. Beringer et al. (Particle Data Group), PR D86, 010001 (2012) and 2013 partial"},
{"lineNum":"   16","line":"//     update for the 2014 edition (URL: http://pdg.lbl.gov)"},
{"lineNum":"   17","line":"    const double mp = 0.938272046; // Mass of proton [GeV/c^2]"},
{"lineNum":"   18","line":"    const double proton_mass = 0.938272046; // Mass of proton [GeV/c^2]"},
{"lineNum":"   19","line":"    const double me = 0.510998928e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   20","line":"    const double electron_mass = 0.510998928e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   21","line":"    const double mmu = 0.1056583715; // Mass of muon [GeV/c^2]"},
{"lineNum":"   22","line":"    const double muon_mass = 0.1056583715; // Mass of muon [GeV/c^2]"},
{"lineNum":"   23","line":"    const double e = 1.602176565e-19; // Charge of proton [C]"},
{"lineNum":"   24","line":"#else"},
{"lineNum":"   25","line":"#if PDG_VERSION == 2010"},
{"lineNum":"   26","line":"//  K. Nakamura et al. (Particle Data Group), J. Phys. G 37, 075021 (2010)"},
{"lineNum":"   27","line":"    const std::string pdg_year(\"2010\");"},
{"lineNum":"   28","line":"    const double mp = 0.938272013; // Mass of proton [GeV/c^2]"},
{"lineNum":"   29","line":"    const double proton_mass = 0.938272013; // Mass of proton [GeV/c^2]"},
{"lineNum":"   30","line":"    const double me = 0.510998910e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   31","line":"    const double electron_mass = 0.510998910e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   32","line":"    const double mmu = 0.105658367; // Mass of muon [GeV/c^2]"},
{"lineNum":"   33","line":"    const double muon_mass = 0.105658367; // Mass of muon [GeV/c^2]"},
{"lineNum":"   34","line":"    const double e = 1.602176487e-19; // Charge of proton [C]"},
{"lineNum":"   35","line":"#else"},
{"lineNum":"   36","line":"#if PDG_VERSION == 2008 // PDG_VERSION2008"},
{"lineNum":"   37","line":"//  C. Amsler et al. (Particle Data Group), Physics Letters B667, 1 (2008)"},
{"lineNum":"   38","line":"    const std::string pdg_year(\"2008\");"},
{"lineNum":"   39","line":"    const double mp = 0.938272013; // Mass of proton [GeV/c^2]"},
{"lineNum":"   40","line":"    const double proton_mass = 0.938272013; // Mass of proton [GeV/c^2]"},
{"lineNum":"   41","line":"    const double me = 0.510998910e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   42","line":"    const double electron_mass = 0.510998910e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   43","line":"    const double mmu = 0.105658367; // Mass of muon [GeV/c^2]"},
{"lineNum":"   44","line":"    const double muon_mass = 0.105658367; // Mass of muon [GeV/c^2]"},
{"lineNum":"   45","line":"    const double e = 1.602176487e-19; // Charge of proton [C]"},
{"lineNum":"   46","line":"#if PDG_VERSION == 2006 // PDG_VERSION2006"},
{"lineNum":"   47","line":"//  S. Eidelman et al. (Particle Data Group), Phys. Lett. B 592, 1 (2004)"},
{"lineNum":"   48","line":"    const std::string pdg_year(\"2006\");"},
{"lineNum":"   49","line":"    const double mp = 0.938272029; // Mass of proton [GeV/c^2]"},
{"lineNum":"   50","line":"    const double proton_mass = 0.938272029; // Mass of proton [GeV/c^2]"},
{"lineNum":"   51","line":"    const double me = 0.510998918e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   52","line":"    const double electron_mass = 0.510998918e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   53","line":"    const double mmu = 0.105658369; // Mass of muon [GeV/c^2]"},
{"lineNum":"   54","line":"    const double muon_mass = 0.105658369; // Mass of muon [GeV/c^2]"},
{"lineNum":"   55","line":"    const double e = 1.60217653e-19; // Charge of proton [C]"},
{"lineNum":"   56","line":"#else"},
{"lineNum":"   57","line":"#if PDG_VERSION == 2004"},
{"lineNum":"   58","line":"//  S. Eidelman et al. (Particle Data Group), Phys. Lett. B 592, 1 (2004)"},
{"lineNum":"   59","line":"    const std::string pdg_year(\"2004\");"},
{"lineNum":"   60","line":"    const double mp = 0.938272029; // Mass of proton [GeV/c^2]"},
{"lineNum":"   61","line":"    const double proton_mass = 0.938272029; // Mass of proton [GeV/c^2]"},
{"lineNum":"   62","line":"    const double me = 0.510998918e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   63","line":"    const double electron_mass = 0.510998918e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   64","line":"    const double mmu = 0.105658369; // Mass of muon [GeV/c^2]"},
{"lineNum":"   65","line":"    const double muon_mass = 0.105658369; // Mass of muon [GeV/c^2]"},
{"lineNum":"   66","line":"    const double e = 1.60217653e-19; // Charge of proton [C]"},
{"lineNum":"   67","line":"#else"},
{"lineNum":"   68","line":"#if PDG_VERSION == -1 // legacy"},
{"lineNum":"   69","line":"    const std::string pdg_year(\"legacy\");"},
{"lineNum":"   70","line":"    const double mp = 0.93827203; // Mass of proton [GeV/c^2]"},
{"lineNum":"   71","line":"    const double proton_mass = 0.93827203; // Mass of proton [GeV/c^2]"},
{"lineNum":"   72","line":"    const double me = 0.51099892e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   73","line":"    const double electron_mass = 0.51099892e-3; // Mass of electron [GeV/c^2]"},
{"lineNum":"   74","line":"    const double mmu = 0.105658369; // Mass of muon [GeV/c^2]"},
{"lineNum":"   75","line":"    const double muon_mass = 0.105658369; // Mass of muon [GeV/c^2]"},
{"lineNum":"   76","line":"    const double e = 1.6021892e-19; // Charge of proton [C]"},
{"lineNum":"   77","line":"#else"},
{"lineNum":"   78","line":"#error \"No selection for physicsl_constants PDG_VERSION set\""},
{"lineNum":"   79","line":"#endif // legacy"},
{"lineNum":"   80","line":"#endif // PDG_VERSION 2004"},
{"lineNum":"   81","line":"#endif // PDG_VERSION 2006"},
{"lineNum":"   82","line":"#endif // PDG_VERSION 2008"},
{"lineNum":"   83","line":"#endif // PDG_VERSION 2010"},
{"lineNum":"   84","line":"#endif // PDG_VERSION 2012"},
{"lineNum":"   85","line":"    const double c = 299792458.0; // Speed of light [m/s]"},
{"lineNum":"   86","line":""},
{"lineNum":"   87","line":"    const double mu0 = 4*mconstants::pi*1.0e-7; // Permittivity of free space [H/m]"},
{"lineNum":"   88","line":"    const double epsilon0 = 1.0/(c*c*mu0); // Permeability of free space [F/m]"},
{"lineNum":"   89","line":""},
{"lineNum":"   90","line":"    // Classical radius of a particle = e^2/(4 pi epsilon0 m c^2)"},
{"lineNum":"   91","line":"    const double re = e/(4*mconstants::pi*epsilon0*me*1.0e9); // Classical"},
{"lineNum":"   92","line":"                                                   // radius of electron [m]"},
{"lineNum":"   93","line":"    const double rp = e/(4*mconstants::pi*epsilon0*mp*1.0e9); // Classical"},
{"lineNum":"   94","line":"                                                   // radius of proton [m]"},
{"lineNum":"   95","line":"    const double rmu = e/(4*mconstants::pi*epsilon0*mmu*1.0e9); // Classical"},
{"lineNum":"   96","line":"                                                   // radius of muon [m]"},
{"lineNum":"   97","line":""},
{"lineNum":"   98","line":"    const int proton_charge = 1; // Charge in units of e"},
{"lineNum":"   99","line":"    const int antiproton_charge = -1; // Charge in units of e"},
{"lineNum":"  100","line":"    const int electron_charge = -1; // Charge in units of e"},
{"lineNum":"  101","line":"    const int positron_charge = 1; // Charge in units of e"},
{"lineNum":"  102","line":"    const int muon_charge = -1; // Charge in units of e"},
{"lineNum":"  103","line":"    const int antimuon_charge = 1; // Charge in units of e"},
{"lineNum":"  104","line":"}"},
{"lineNum":"  105","line":""},
{"lineNum":"  106","line":"#endif /* PHYSICAL_CONSTANTS_H_ */"},
]};
var percent_low = 25;var percent_high = 75;
var header = { "command" : "fodo_cxx", "date" : "2021-12-15 09:03:42", "instrumented" : 1, "covered" : 1,};
var merged_data = [];
