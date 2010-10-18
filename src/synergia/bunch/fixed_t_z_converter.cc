#include "utils/multi_array_typedefs.h"
#include "components/foundation/reference_particle.h"
#include "fixed_t_z_converter.h"
#include "bunch.h"

void
Fixed_t_z_zeroth::fixed_t_to_fixed_z(Bunch &bunch)
{
    // This routine just copied/translated from older branch.
    // Check/fix me please.
    double gamma_ref = bunch.get_reference_particle().get_gamma();
    double mass = bunch.get_mass();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double f_x = particles[part][Bunch::xp] / mass;
        double f_y = particles[part][Bunch::yp] / mass;
        double gamma = gamma_ref - particles[part][Bunch::tp] / mass;
        double beta = sqrt(1.0 - (1 + f_x * f_x + f_y * f_y) / (gamma * gamma));
        particles[part][Bunch::z] /= -gamma * beta;
    }
}

void
Fixed_t_z_zeroth::fixed_z_to_fixed_t(Bunch &bunch)
{
    // This routine just copied/translated from older branch.
    // Check/fix me please.
    double gamma_ref = bunch.get_reference_particle().get_gamma();
    double mass = bunch.get_mass();
    MArray2d_ref particles = bunch.get_local_particles();
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double f_x = particles[part][Bunch::xp] / mass;
        double f_y = particles[part][Bunch::yp] / mass;
        double gamma = gamma_ref - particles[part][Bunch::tp] / mass;
        double beta = sqrt(1.0 - (1 + f_x * f_x + f_y * f_y) / (gamma * gamma));
        particles[part][Bunch::z] *= -gamma * beta;
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
