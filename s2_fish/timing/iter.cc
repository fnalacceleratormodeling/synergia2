#include <iostream>
#include "array_1d.h"
#include "array_2d.h"
#include "array_3d.h"
#include "array_nd.h"
#include "vector_helper.h"
#include "mytimer.h"
#include <mpi.h>

int main(int argc, char **argv)
{
    int xmax = 1000000;
    int ymax = 7;
    int zmax = 8192;
    //~ int xmax = 4;
    //~ int ymax = 4;
    //~ int zmax = 4;
    MPI_Init(&argc, &argv);

    reset_timer();
    Array_1d<double> a1(xmax*ymax);
    timer("a1 construct");

    a1.set_all(0.0);
    timer("a1 set_all");
    a1.set_all(0.0);
    timer("a1 set_all2");
    for (Array_1d<double>::Iterator it = a1.begin(); it != a1.end(); ++it) {
        *it = 0.0;
    }
    timer("a1 Iterator");
    for (Array_1d<double>::Iterator it = a1.begin(); it != a1.end(); ++it) {
        *it = 0.0;
    }
    timer("a1 Iterator1");

    for (int i = 0 ; i < xmax*ymax; ++i) {
        a1(i) = 0.0;
    }
    timer("a1 dumb set all");
    for (int i = 0 ; i < xmax*ymax; ++i) {
        a1(i) = 0.0;
    }
    timer("a1 dumb set all2");

    //~ reset_timer();
    //~ Array_2d<double> a2(xmax,ymax);
    //~ timer("a2 construct");

    //~ a2.set_all(0.0);
    //~ timer("a2 set_all");
    //~ a2.set_all(0.0);
    //~ timer("a2 set_all2");

    //~ for(Array_2d<double>::Iterator it = a2.begin(); it != a2.end(); ++it) {
    //~ *it = 0.0;
    //~ }
    //~ timer("a2 Iterator");
    //~ for(Array_2d<double>::Iterator it = a2.begin(); it != a2.end(); ++it) {
    //~ *it = 0.0;
    //~ }
    //~ timer("a2 Iterator2");

    //~ for(int i = 0 ; i < xmax; ++i) {
    //~ for(int j = 0; j<ymax; ++j) {
    //~ a2(i,j) = 0.0;
    //~ }
    //~ }
    //~ timer("a2 dumb set all");
    //~ for(int i = 0 ; i < xmax; ++i) {
    //~ for(int j = 0; j<ymax; ++j) {
    //~ a2(i,j) = 0.0;
    //~ }
    //~ }
    //~ timer("a2 dumb set all2");

    //~ reset_timer();
    //~ Array_2d<double> at2(ymax,xmax);
    //~ timer("at2 construct");

    //~ at2.set_all(0.0);
    //~ timer("at2 set_all");
    //~ at2.set_all(0.0);
    //~ timer("at2 set_all2");

    //~ for(Array_2d<double>::Iterator it = at2.begin(); it != at2.end(); ++it) {
    //~ *it = 0.0;
    //~ }
    //~ timer("at2 Iterator");
    //~ for(Array_2d<double>::Iterator it = at2.begin(); it != at2.end(); ++it) {
    //~ *it = 0.0;
    //~ }
    //~ timer("at2 Iterator2");

    //~ for(int i = 0 ; i < ymax; ++i) {
    //~ for(int j = 0; j<xmax; ++j) {
    //~ at2(i,j) = 0.0;
    //~ }
    //~ }
    //~ timer("at2 dumb set all");
    //~ for(int i = 0 ; i < ymax; ++i) {
    //~ for(int j = 0; j<xmax; ++j) {
    //~ at2(i,j) = 0.0;
    //~ }
    //~ }
    //~ timer("at2 dumb set all2");

    //~ xmax = 2;
    //~ ymax = 256;

    //~ reset_timer();
    //~ Array_3d<double> a3(xmax,ymax,zmax);
    //~ timer("a3 construct");

    //~ a3.set_all(0.0);
    //~ timer("a3 set_all");
    //~ a3.set_all(0.0);
    //~ timer("a3 set_all2");

    //~ for(Array_3d<double>::Iterator it = a3.begin(); it != a3.end(); ++it) {
    //~ *it = 0.0;
    //~ }
    //~ timer("a3 Iterator");
    //~ for(Array_3d<double>::Iterator it = a3.begin(); it != a3.end(); ++it) {
    //~ *it = 0.0;
    //~ }
    //~ timer("a3 Iterator2");

    //~ for(int i = 0 ; i < xmax; ++i) {
    //~ for(int j = 0; j<ymax; ++j) {
    //~ for(int k = 0; k<zmax; ++k) {
    //~ a3(i,j,k) = 0.0;
    //~ }
    //~ }
    //~ }
    //~ timer("a3 dumb set all ");
    //~ for(int i = 0 ; i < xmax; ++i) {
    //~ for(int j = 0; j<ymax; ++j) {
    //~ for(int k = 0; k<zmax; ++k) {
    //~ a3(i,j,k) = 0.0;
    //~ }
    //~ }
    //~ }
    //~ timer("a3 dumb set all2");

    reset_timer();
    xmax = 2;
    ymax = 256;
    Array_nd<double> an(vector3(xmax, ymax, zmax));
    timer("an construct");

    an.set_all(0.0);
    timer("an set_all");
    an.set_all(0.0);
    timer("an set_all2");

    for (Array_nd<double>::Iterator it = an.begin(); it != an.end(); ++it) {
        *it = 0.0;
    }
    timer("an Iterator");
    for (Array_nd<double>::Iterator it = an.begin(); it != an.end(); ++it) {
        *it = 0.0;
    }
    timer("an Iterator2");

    for (int i = 0 ; i < xmax; ++i) {
        for (int j = 0; j < ymax; ++j) {
            for (int k = 0; k < zmax; ++k) {
                an(vector3(i, j, k)) = 0.0;
            }
        }
    }
    timer("an dumb set all ");
    for (int i = 0 ; i < xmax; ++i) {
        for (int j = 0; j < ymax; ++j) {
            for (int k = 0; k < zmax; ++k) {
                an(vector3(i, j, k)) = 0.0;
            }
        }
    }
    timer("an dumb set all2");

    std::vector<int> indices(3);
    for (indices[0] = 0 ; indices[0] < xmax; ++indices[0]) {
        for (indices[1] = 0; indices[1] < ymax; ++indices[1]) {
            for (indices[2] = 0; indices[2] < zmax; ++indices[2]) {
                an(indices) = 0.0;
            }
        }
    }
    timer("an less dumb set all ");
    for (indices[0] = 0 ; indices[0] < xmax; ++indices[0]) {
        for (indices[1] = 0; indices[1] < ymax; ++indices[1]) {
            for (indices[2] = 0; indices[2] < zmax; ++indices[2]) {
                an(indices) = 0.0;
            }
        }
    }
    timer("an less dumb set all2");

    reset_timer();
    Array_nd<double> an_yslice(an.slice(
                                   vector3(Range(0),
                                           Range(),
                                           Range(0))));
    std::cout << "an_yslice rank = " << an_yslice.get_rank() << std::endl;
    timer("an_yslice construct");

    an_yslice.set_all(0.0);
    timer("an_yslice set_all");
    an_yslice.set_all(0.0);
    timer("an_yslice set_all2");

    Array_nd<double>::Iterator  end = an_yslice.end();
    for (Array_nd<double>::Iterator it = an_yslice.begin(); it != end; ++it) {
        *it = 0.0;
    }
    timer("an_yslice Iterator");
    for (Array_nd<double>::Iterator it = an_yslice.begin(); it != an_yslice.end(); ++it) {
        *it = 0.0;
    }
    timer("an_yslice Iterator2");

    for (int j = 0; j < ymax; ++j) {
        an_yslice(vector1(j)) = 0.0;
    }
    timer("an_yslice dumb set all ");

    for (int j = 0; j < ymax; ++j) {
        an_yslice(vector1(j)) = 0.0;
    }
    timer("an_yslice dumb set all2 ");

    reset_timer();
    Array_nd<double> an_halfslice(an.slice(
                                   vector3(Range(Range::begin,Range::end,2),
                                           Range(Range::begin,Range::end,2),
                                           Range(Range::begin,Range::end,2))));
    timer("an_halfslice construct");

    an_halfslice.set_all(0.0);
    timer("an_halfslice set_all");
    an_halfslice.set_all(0.0);
    timer("an_halfslice set_all2");

    end = an_halfslice.end();
    for (Array_nd<double>::Iterator it = an_halfslice.begin();
        it != end;
        ++it) {
        *it = 0.0;
    }
    timer("an_halfslice Iterator");
    for (Array_nd<double>::Iterator it = an_halfslice.begin(); it != an_halfslice.end(); ++it) {
        *it = 0.0;
    }
    timer("an_halfslice Iterator2");

    std::vector<int> shape = an_halfslice.get_shape();
    for (int i = 0 ; i < shape[0]; ++i) {
        for (int j = 0; j < shape[1]; ++j) {
            for (int k = 0; k < shape[2]; ++k) {
                an_halfslice(vector3(i, j, k)) = 0.0;
            }
        }
    }
    timer("an_halfslice dumb set all ");
    for (int i = 0 ; i < shape[0]; ++i) {
        for (int j = 0; j < shape[1]; ++j) {
            for (int k = 0; k < shape[2]; ++k) {
                an_halfslice(vector3(i, j, k)) = 0.0;
            }
        }
    }
    timer("an_halfslice dumb set all2");

    return 0;
}
