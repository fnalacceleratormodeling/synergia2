#ifndef HAVE_TRANSVERSE_FFTW_H
#define HAVE_TRANSVERSE_FFTW_H

#include "macro_bunch_store.h"
#include "array_nd/array_2d.h"
#include "array_nd/array_3d.h"
#include "field_domain.h"

// This solver is based on ORBIT program transverse solver
// check the corresponding paper for details

struct Extremum {
    double max;
    double min;
    Extremum() : max(0), min(0) { }
};

class TransverseSolver 
{
    //Private data &  methods
    
    int nXBins, nYBins, nZBins; // Grid size
    bool includeLocalDensity;
    std::vector<int> indicies; // Indicies inside of the distribution area ( [ mean - size / 2, mean + size / 2] )
    Array_2d<std::complex<double> > greensFunc; // Green function
    Array_1d<int> xBin, yBin, zBin; // grid charge bins
    Array_1d<double> xFractBin, yFractBin; // fractioned bin (the shift)
    Array_2d<double> rho; // Density
    Array_1d<double> xGrid, yGrid; // 2D grid
    Array_2d<double> fscx, fscy; // Forces at grid cells
    Array_1d<double> stds, means;
    Array_1d<double> z_density;
    double xGridMin, xGridMax, yGridMin, yGridMax; // Grid extremas
    double dx, dy, dz; // Grid dx, dy
    double eps; // smoothing parameter
    void form_grid(Macro_bunch_store& mbs); // form grid for charge distribution via max and min (not recommended)
    void form_grid(Macro_bunch_store& mbs, const std::vector<double>& means, 
		    const std::vector<double>& sizes); // form grid via means and stds
    void deposit_charge(Macro_bunch_store& mbs); // distribue charge on celss (PIC tecnique)
    void calc_local_z_density(Macro_bunch_store& mbs); // density along z direction
     
public:
   // Constructor
    TransverseSolver(int nx, int ny, double e, bool includeZDensity = true) 
	    : nXBins(nx), nYBins(ny), nZBins(nx), eps(e), includeLocalDensity(includeZDensity)  {  }
    // Method includes fftw and convolution
    void kick_transverse_charge(Macro_bunch_store &mbs, double tau, const std::vector<double>& means, const std::vector<double>& sizes);
    // Method implements direct charge calculation on grid cells, used for testing
    void kick_transverse_charge_direct(Macro_bunch_store &mbs);
    // Find min and max value for selected coordinate of the bunch
    static void findBunchExtrema(const Macro_bunch_store& mbs, int index, Extremum& e);
    const Array_nd<double>& getXGrid() const { return xGrid; }
    const Array_nd<double>& getYGrid() const { return yGrid; }
    const Array_nd<double>& getRho() const { return rho; }
    const Array_nd<double>& getFscx() const { return fscx; }
    const Array_nd<double>& getFscy() const { return fscy; }
    
};


#endif // HAVE_TRANSVERSE_FFTW_H

