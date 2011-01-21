#include "aperture.h"

void
apply_circular_aperture(Bunch & bunch, Lattice_element_slices & slices)
{
    const double large = 1.0e30;
    double min_radius(large);
    bool found_radius(false);
    for (Lattice_element_slices::const_iterator it = slices.begin(); it
            != slices.end(); ++it) {
        if ((*it)->get_lattice_element().has_double_attribute("aperture_radius")) {
            double radius((*it)->get_lattice_element().get_double_attribute(
                    "aperture_radius"));
            found_radius = true;
            if (radius < min_radius) {
                min_radius = radius;
            }
        }
    }
    if (min_radius > 0.0) {
        if (!found_radius) {
            min_radius = default_aperture_radius;
        }
        double radius2 = min_radius * min_radius;
        MArray2d_ref particles(bunch.get_local_particles());
        int kept = 0;
        int discarded = 0;
        int local_num = bunch.get_local_num();
        for (int part = 0; part < local_num; ++part) {
            bool try_discard = true;
            while (try_discard) {
                double r2 = particles[part][Bunch::x]
                        * particles[part][Bunch::x] + particles[part][Bunch::y]
                        * particles[part][Bunch::y];
                if (r2 > radius2) {
                    ++discarded;
                    --local_num;
                    if (part == local_num) {
                        // No more particles left
                        try_discard = false;
                    } else {
                        // Move the last particle into this newly empty position
                        int last = local_num;
                        particles[part][0] = particles[last][0];
                        particles[part][1] = particles[last][1];
                        particles[part][2] = particles[last][2];
                        particles[part][3] = particles[last][3];
                        particles[part][4] = particles[last][4];
                        particles[part][5] = particles[last][5];
                        particles[part][6] = particles[last][6];
                    }
                } else {
                    ++kept;
                    try_discard = false;
                }
            }
        }
//        std::cout << "kept = " << kept << ", discarded = " << discarded
//                << std::endl;
        bunch.set_local_num(local_num);
        bunch.update_total_num();
    }
}

