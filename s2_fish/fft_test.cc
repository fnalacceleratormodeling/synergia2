#include <iostream>
#include <fftw3.h>
#include "array_nd.h"
#include "array_3d.h"
#include <complex>

int 
main()
{
    int mr = 2;
    int mphi = 4;
    int mz = 4;
    Array_3d<double> b(mr,mphi,mz);
    Array_3d<double> bnew(mr,mphi,mz);
    for(int i=0; i<mr; ++i){
        for(int j=0; j<mphi; ++j) {
            for(int k=0; k<mphi; ++k) {
                b(i,j,k) = b(i,j,k) = (i+1)*100.0+(j+1)*10.0+(k+1)*.011;
            }
        }
    }
    // the shape of the FFT'd array (shape_lm) is halved in the third 
    // dimension because of the peculiar (but efficient) way FFTW does
    // real-to-complex transforms.
    std::vector<int> shape = b.get_shape();
    std::vector<int> shape_lm = vector3(shape[0],shape[1],shape[2]/2+1);
    Array_3d<std::complex<double> > blm(shape_lm);
        
    fftw_plan plan = fftw_plan_many_dft_r2c(2,
        &shape[1], shape[0],
        b.get_data_ptr(),
        NULL, 1, shape[1]*shape[2],
        reinterpret_cast<double (*)[2]>(blm.get_data_ptr()),
        NULL, 1, shape_lm[1]*shape_lm[2],
        FFTW_ESTIMATE);
    fftw_execute(plan);

    plan = fftw_plan_many_dft_c2r(2,
                    &shape[1], shape_lm[0],
                    reinterpret_cast<double (*)[2]>(blm.get_data_ptr()),
                    NULL, 1, shape_lm[1]*shape_lm[2],
                    bnew.get_data_ptr(),
                    NULL, 1, shape[1]*shape[2],
                    FFTW_ESTIMATE);
    
    fftw_execute(plan);
    bnew.scale(1.0/(shape[1]*shape[2]));

    b.print("b");
    bnew.print("bnew");

    for(int i=0; i<mr; ++i){
        for(int j=0; j<mphi; ++j) {
            for(int k=0; k<mphi; ++k) {
        
    return 0;
}
