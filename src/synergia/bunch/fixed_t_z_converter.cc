#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/foundation/reference_particle.h"
#include "fixed_t_z_converter.h"
#include "bunch.h"

#include <iostream>
#include <cmath>

void
Fixed_t_z_zeroth::fixed_t_to_fixed_z(Bunch &bunch)
{
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        // ct in accelerator frame
        particles[part][Bunch::cdt] = -particles[part][Bunch::z]
                / (gamma * beta);

        // p'_{x,y,z} in beam frame
        double pxp = particles[part][Bunch::xp] * p_ref;
        double pyp = particles[part][Bunch::yp] * p_ref;
        double pzp = particles[part][Bunch::zp] * p_ref;
        double p_perp2 = pxp * pxp + pyp * pyp;
        // E'/c in beam frame
        double Epoc = std::sqrt(p_perp2 + pzp * pzp + m * m);
        double pz = gamma * (pzp + beta * Epoc);
        // dpop = (p - p_ref)/p_ref
        double p = std::sqrt(p_perp2 + pz * pz);
        particles[part][Bunch::dpop] = (p - p_ref) / p_ref;
    }
}

void
Fixed_t_z_zeroth::fixed_z_to_fixed_t(Bunch &bunch)
{
    double gamma = bunch.get_reference_particle().get_gamma();
    double beta = bunch.get_reference_particle().get_beta();
    double m = bunch.get_mass();
    double p_ref = bunch.get_reference_particle().get_momentum();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        // z in beam rest frame
        particles[part][Bunch::z] = -1.0*gamma * beta * particles[part][Bunch::cdt];

        // total momentum in accelerator frame
        double p = p_ref + particles[part][Bunch::dpop] * p_ref;
        // E/c in accelerator frame
        double Eoc = std::sqrt(p * p + m * m);
        // p_{x,y,z} in accelerator frame
        double px = particles[part][Bunch::xp] * p_ref;
        double py = particles[part][Bunch::yp] * p_ref;
        double pz2 = p * p - px * px - py * py;
        if (pz2 < 0.0) {
            throw std::runtime_error(
                    "Fixed_t_z_zeroth::fixed_z_to_fixed_t: particle has negative pz^2");
        }
        double pz = std::sqrt(pz2);
        // zp = pz/p_{ref}^{total}
        particles[part][Bunch::zp] = gamma * (pz - beta * Eoc) / p_ref;

        // n.b. in the zeroth approximation, the transformation from
        //      t' = gamma cdt to t' = 0
        //      is a no-op.
    }
}

void
Fixed_t_z_ballistic::fixed_t_to_fixed_z(Bunch &bunch)
{
    std::cout << "stub: ballistic fixed_t_to_fixed_z\n";
}

void
Fixed_t_z_ballistic::fixed_z_to_fixed_t(Bunch &bunch)
{
    std::cout << "stub: ballistic fixed_z_to_fixed_t\n";
}
