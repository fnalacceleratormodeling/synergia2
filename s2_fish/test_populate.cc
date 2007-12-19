#include "populate.h"

extern void
array_2d_to_octave_file(const Array_2d<double> &array, const std::string filename); 

int main()
{
    int num_particles=1000;
    Array_2d<double> p(7,num_particles);
    Array_2d<double> covs(6,6);
    covs.set_all(0.0);
    Array_1d<double> means(6);
    means.set_all(0.0);

    for (int i=0; i<6; ++i) {
        means.at(i) = 0.01*i;
        for (int j=0; j<6; ++j) {
            if (i==j) {
                covs.at(i,j) = i+1.0;
            }
            if ((i+1 == j) && (i%2 == 0)) {
                covs.at(i,j) = 0.1*(j+1.0);
                covs.at(j,i) = 0.1*(j+1.0);
            }
        }
    }
    covs.print("requested covs");
    means.print("requested means");
    
    populate_6d_gaussian(p,means,covs,0);
    array_2d_to_octave_file(p,"r.dat");
    std::cout << "particles dumped to r.dat\n";
    return 0;
}
