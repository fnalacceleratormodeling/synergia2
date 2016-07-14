#include "aperture.h"

#include <omp.h>

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
        int npart = bunch.get_local_num();

        int nt;

        #pragma omp parallel
        { nt = omp_get_num_threads(); }

        char * discard = new char[npart];
        int * discard_count = new int[nt];

        #pragma omp parallel shared(nt, npart, particles, radius2, discard, discard_count)
        {
            int it = omp_get_thread_num();

            int l = npart / nt;
            int s = it * l;
            int e = (it==nt-1) ? npart : (s+l);

            discard_count[it] = 0;

            for (int part = s; part < e; ++part) 
            {
                double r2 = particles[part][Bunch::x] * particles[part][Bunch::x] 
                          + particles[part][Bunch::y] * particles[part][Bunch::y];

                if (r2 > radius2 )
                {
                    discard[s + discard_count[it]] = part;
                    ++ discard_count[it];
                }
            }
        }

        int total_discarded = discard_count[0];

        for (int t = 1; t < nt; ++t)
        {
            std::memcpy( discard + discard_count[t-1], discard + t*npart/nt, discard_count[t] );
            total_discarded += discard_count[t];
        }

        int last = npart-1;

        for (int n = 0; n < total_discarded; ++n, --last)
        {
            while ( last == discard[total_discarded-1] )
            {
                --last;
                --total_discarded;

                if (n == total_discarded) goto done;
            }

            int idx = discard[n];

            particles[idx][0] = particles[last][0];
            particles[idx][1] = particles[last][1];
            particles[idx][2] = particles[last][2];
            particles[idx][3] = particles[last][3];
            particles[idx][4] = particles[last][4];
            particles[idx][5] = particles[last][5];
            particles[idx][6] = particles[last][6];
        }

done:
        bunch.set_local_num(npart - total_discarded);
        bunch.update_total_num();
    }
}

