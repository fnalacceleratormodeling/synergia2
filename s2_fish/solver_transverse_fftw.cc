#include <iostream>
#include "solver_transverse_fftw.h"
#include <fftw3.h>
#include "math_constants.h"
#include "basic_toolkit/PhysicsConstants.h"

static std::ofstream logFile("log_transverse_solver.txt");

template<class T> void dbg_print_array_2d(const std::string& name, const Array_2d<T>& array) {
#ifdef DBG
    logFile << name << ":" << std::endl;
    std::vector<int> shapes = array.get_shape();
    for (int i = 0; i < shapes[0]; ++i) {
	for (int j = 0; j < shapes[1]; ++j) {
	    logFile << array(i,j) << "   ";
	}
	logFile << std::endl;
    }
    logFile << std::endl;
#endif
}

void TransverseSolver::findBunchExtrema(const Macro_bunch_store& mbs, int index, Extremum& e)
{
    e.min = e.max = mbs.local_particles(index,0);
    int local_num = mbs.local_num;   
    for (int i = 0; i < local_num; ++i) {
    	double value = mbs.local_particles(index,i);
	if (e.min > value) {
	    e.min = value;
	}
	if (e.max < value) {
	    e.max = value;
	}
               	
    }
}

inline std::complex<double> calcGreensFunc(double rTransX, double rTransY, double eps)
{
    double  rTot2 = rTransX*rTransX + rTransY*rTransY + eps;
    // 1/2*PI factor is from Gauss law
    return (1/(2.*pi))*std::complex<double>( rTransX / rTot2, rTransY / rTot2); 
}

void TransverseSolver::form_grid(Macro_bunch_store& mbs, const std::vector<double>& means, const std::vector<double>& sizes ) {

    double xMax, xMin, yMax, yMin;
    xMax = means[0] + sizes[0] / 2.;
    xMin = means[0] - sizes[0] / 2.;
    yMax = means[1] + sizes[1] / 2.;
    yMin = means[1] - sizes[1] / 2.;
            
    double dxExtra = 0.5*(xMax - xMin);
    double dyExtra = 0.5*(yMax - yMin);

    xGridMin = xMin - dxExtra;
    xGridMax = xMax + dxExtra;
    yGridMin = yMin - dyExtra;
    yGridMax = yMax + dyExtra;

    dx = (xGridMax - xGridMin) / double(nXBins);
    dy = (yGridMax - yGridMin) / double(nYBins);

    xGrid.reshape(nXBins + 1);
    yGrid.reshape(nYBins + 1);
    
    for (int iX = 0; iX <=  nXBins; ++iX) {
	xGrid(iX) =  xGridMin + double(iX) * dx;
    }
    
    for (int iY = 0; iY<= nYBins; ++iY){
	yGrid(iY) = yGridMin + double(iY) * dy;
    }
}


void TransverseSolver::form_grid(Macro_bunch_store& mbs) {

    double xMax, xMin, yMax, yMin;
    Extremum eX, eY;
    findBunchExtrema(mbs, 0, eX);
    findBunchExtrema(mbs, 2, eY);
    xMax = eX.max;
    xMin = eX.min;
    yMax = eY.max;
    yMin = eY.min;
            
    double dxExtra = 0.5*(xMax - xMin);
    double dyExtra = 0.5*(yMax - yMin);

    xGridMin = xMin - dxExtra;
    xGridMax = xMax + dxExtra;
    yGridMin = yMin - dyExtra;
    yGridMax = yMax + dyExtra;

    dx = (xGridMax - xGridMin) / double(nXBins);
    dy = (yGridMax - yGridMin) / double(nYBins);

    xGrid.reshape(nXBins + 1);
    yGrid.reshape(nYBins + 1);
    
    for (int iX = 0; iX <=  nXBins; ++iX) {
	xGrid(iX) =  xGridMin + double(iX) * dx;
    }
    
    for (int iY = 0; iY<= nYBins; ++iY){
	yGrid(iY) = yGridMin + double(iY) * dy;
    }
}


void TransverseSolver::deposit_charge(Macro_bunch_store& mbs) {
    
    xBin.reshape(mbs.local_num);
    yBin.reshape(mbs.local_num);
    xFractBin.reshape(mbs.local_num);
    yFractBin.reshape(mbs.local_num);
    rho.reshape(nXBins+2,nYBins+2);
    rho.set_all(0.);
    greensFunc.reshape(nXBins, nYBins);
    
    for(int i=0; i < mbs.local_num; ++i) {
	int iX = int( (mbs.local_particles(0,i) - xGridMin)/ dx);
	int iY = int( (mbs.local_particles(2,i) - yGridMin)/ dy);
	//Ignore particles that are out of grid size
	if (iX < 0  || iX > nXBins || iY < 0 || iY > nYBins)	
	     continue;
	indicies.push_back(i);
	xBin(i) = iX; 
	yBin(i) = iY;
	xFractBin(i) = (mbs.local_particles(0,i) - xGrid(iX)) / dx;
	yFractBin(i) = (mbs.local_particles(2,i) - yGrid(iY)) / dy;
    
	// Bilinear binning:
												        
	rho(iX,iY) += (1.- xFractBin(i)) * (1. - yFractBin(i));
	rho(iX,iY+1) += (1.- xFractBin(i)) * yFractBin(i);
	rho(iX+1,iY) += xFractBin(i) * (1.-yFractBin(i));
	rho(iX+1,iY+1) += xFractBin(i) * yFractBin(i);
    }
         
}

void TransverseSolver::calc_local_z_density(Macro_bunch_store& mbs) {
    
    // Init grid along Z direction to calc local line density
    Extremum eZ;
    findBunchExtrema(mbs, 4, eZ);
    double zMax = eZ.max;
    double zMin = eZ.min;
            
    double dzExtra = 0.0; //0.001*(zMax - zMin);
    double zGridMin = zMin - dzExtra;
    double zGridMax = zMax + dzExtra;

    dz = (zGridMax - zGridMin) / double(nZBins);
    
    size_t num = mbs.local_num;
    zBin.reshape(num);
    z_density.reshape(nZBins + 1);
    z_density.set_all(0.0);
    
    //logFile << "zMin: " << zGridMin << std::endl;
    //logFile << "zMax: " << zGridMax << std::endl;
    //logFile << "dz: " << dz  << std::endl;
    
    // Calc density (linear binning)
    double sum = 0.0;
    for (int i = 0; i < num; ++i) {
	int iZ = int( (mbs.local_particles(4,i) - zGridMin) / dz );
	zBin(i) = iZ;
	assert(iZ < nZBins + 1);
	z_density(iZ) += 1.;
	sum += 1.;
    }	
    assert(int(sum) == num); 
    //double avrDensity = sum / (double) (nZBins + 1);
    z_density.scale(1. / ((double)num * dz ) );

}

void TransverseSolver::kick_transverse_charge(Macro_bunch_store &mbs, double tau, const std::vector<double>& means, const std::vector<double>& sizes) {
    int index, iX, iY, iXS, iYS, nPoints;
    double x1, x2, y1, y2;
    double epsSq, rTransX, rTransY, rTot2, rr;
    int nMacros;
   
    // Form grid (2 times larger than real particle extent, check Hockney convolution method for details)
    logFile << "******************** FFTW TRANSVERSE METHOD **************************************" << std::endl;
    
    logFile << "\nForming the grid...\n";
    form_grid(mbs, means, sizes);
    
    // Bin the particles (set rho)
    
    logFile << "Binning the particles...\n";
    deposit_charge(mbs);
        
    // Calculate the Greens funtion grid 
    
    dbg_print_array_2d("rho", rho);
    
    logFile << "Calculating Greens function...\n"; 
    
    epsSq = dx*dy*eps;
        
    for (iY = 0; iY < nYBins/2; iY++) { // assign middle to top-1 rows
	rTransY = dy * iY;
        for (iX = 0; iX < nXBins/2; iX++) {
	    rTransX = iX * dx;
            greensFunc(iX, iY) = calcGreensFunc(rTransX, rTransY, epsSq);
	}

        greensFunc(nXBins/2, iY) = std::complex<double>(0,0); // end point

        for (iX = nXBins/2+1; iX < nXBins; iX++) {
	    rTransX = (iX -  nXBins) * dx;
	    greensFunc(iX, iY) = calcGreensFunc(rTransX, rTransY, epsSq);
	}
    }
    
    for(iX = 0; iX < nXBins; iX++)   // Null the top row:
    {
	greensFunc(iX, nYBins/2) = std::complex<double>(0,0);
    }

    for (iY = nYBins/2+1; iY < nYBins; iY++)  // Bottom rows:
    {
	rTransY = dy * (iY - nYBins);
        for (iX = 0; iX < nXBins/2; iX++)
        {
	    rTransX = (iX) * dx;
       	    greensFunc(iX, iY) = calcGreensFunc(rTransX, rTransY, epsSq);
	}
	
	greensFunc(nXBins/2, iY) = std::complex<double>(0,0);  // end point
        
	for (iX = nXBins/2 + 1; iX < nXBins; iX++) 
	{
	    rTransX = (iX  - nXBins) * dx;
       	    greensFunc(iX, iY) = calcGreensFunc(rTransX, rTransY, epsSq);
	}
	    
    }

    dbg_print_array_2d("greens func", greensFunc);

    //  Calculate the FFT of the Greens Function:

    Array_2d<std::complex<double> > out(nXBins, nYBins), out2(nXBins, nYBins);

    fftw_plan forward_plan_greens =  fftw_plan_dft_2d(nYBins, nXBins, 
					    (fftw_complex*) greensFunc.get_data_ptr(), (fftw_complex*) out.get_data_ptr(), 
					    FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(forward_plan_greens);
    dbg_print_array_2d("fft of greens", out);
    
    // Calculate the FFT of the binned charge distribution:
    
    Array_2d<std::complex<double> > rho_c(nXBins,nYBins); // Complex array rho
    
    for (int i = 0; i < nXBins; ++i) {
	for (int j = 0; j < nYBins; ++j) { 
	    rho_c(i,j) = std::complex<double>(rho(i,j), 0.0);
	}
    }
	    
    fftw_plan forward_plan_charge =  fftw_plan_dft_2d(nYBins, nXBins, 
					    (fftw_complex*) rho_c.get_data_ptr(), (fftw_complex*) out2.get_data_ptr(), 
					    FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(forward_plan_charge);
    dbg_print_array_2d("fft of rho", out2);

    // Do Convolution:
    
    for(int i = 0; i < nXBins; ++i) {
	for (int j = 0; j < nYBins; ++j)
        {
	    out(i,j) = out(i,j)*out2(i,j);
	 }
    }


    dbg_print_array_2d("convolution", out);
    
    // Do Inverse FFT to get the Force:
    
    fftw_plan backward_plan =  fftw_plan_dft_2d(nYBins, nXBins, 
					    (fftw_complex*) out.get_data_ptr(), (fftw_complex*) out2.get_data_ptr(), 
					    FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(backward_plan);
    dbg_print_array_2d("result", out2);
    
    fscx.reshape(nXBins, nYBins);
    fscy.reshape(nXBins, nYBins);
    
    for (int i = 0; i < nXBins; ++i) {
	for (int j = 0; j < nYBins; ++j) {
	    double real = std::abs( out2(i,j).real() ) > 0.00001 ? out2(i,j).real() : 0;
	    double imag = std::abs( out2(i,j).imag() ) > 0.00001 ? out2(i,j).imag() : 0;
	    double factor = 1./ double(nXBins*nYBins*mbs.local_num);
	    fscx(i,j) = factor*real;
	    fscy(i,j) = factor*imag;
	}
    }
    
    // Caclulate perveance (see e_field calculation for detailed comments)
    
    double gamma = -1.*mbs.ref_particle(5);
    double beta = sqrt(gamma*gamma- 1.) / gamma;
    double eps0 = PH_MKS_eps0;
    double c = PH_MKS_c;
    double qe=PH_MKS_e;
    
    double length = 2.0*pi*beta/mbs.units(0); // Bunch length in lab frame
   // double charge  = (PH_CNV_brho_to_p/eps0)*length*mbs.total_current / (beta*c);
    double charge  = (PH_CNV_brho_to_p/eps0)*mbs.bunch_np*qe*mbs.charge;
    double factor = charge* ( 1./ (beta*c) )*(1. / gamma)*mbs.units(1);		
    
    /*
    logFile << "Charge is " << charge << std::endl;
    logFile << "Number of particles is " << mbs.local_num << std::endl;
    logFile << "Length is " << length << std::endl; 
    logFile << "Gamma is " << gamma << std::endl;
    logFile << "Beta is " << beta << std::endl;
    logFile << "Beta*c is " << beta*c << std::endl;
    logFile << "MBS.units(0) is " << mbs.units(0) << std::endl;
    logFile << "MBS.units(1) is " << mbs.units(1) << std::endl;
    */
    
    // Interpolate forces on particles
    
    //FOR DEBUG PRINT
    static int step = 0;
    int num = 20;
    std::ios_base::openmode mode;  
    step == 0 ? mode = std::ios::out : mode = std::ios::app;
    //std::ofstream outf("out_orbit.txt", mode);

    //outf << "\n$$$STEP " << step << " $$$\n" << std::endl;
    ++step;
    
    if (includeLocalDensity) {
	logFile << "Calculate line density along z direction..." << std::endl;   
	calc_local_z_density(mbs);
	//for (int iZ = 0; iZ < nZBins + 1; ++iZ) 
	//    logFile << iZ << " : " << z_density(iZ) << std::endl;
    } else {
	logFile << "Assume uniform charge distribution..." << std::endl;
    }

    // Local line density (if assumed uniform charge distribution) 
    double lambda = 1. / (gamma * length);

    for (std::vector<int>::iterator it = indicies.begin(); it != indicies.end(); ++it)
    {
	int j = *it;
	// Find horizontal force:
        double f1 = fscx(xBin(j), yBin(j));
	double f2 = fscx(xBin(j)+1, yBin(j));
	double f3 = fscx(xBin(j)+1, yBin(j)+1);
	double f4 = fscx(xBin(j), yBin(j)+1);

	double fx = (1. - xFractBin(j)) * (1. - yFractBin(j))*f1 +
         xFractBin(j) * (1. - yFractBin(j))*f2 +
         (1. - xFractBin(j)) * yFractBin(j)*f4 +       
         xFractBin(j) * yFractBin(j) * f3;

	// Find Vertical Force:
        f1 = fscy(xBin(j), yBin(j));
	f2 = fscy(xBin(j)+1, yBin(j));
	f3 = fscy(xBin(j)+1, yBin(j)+1);
	f4 = fscy(xBin(j), yBin(j)+1);
	
	double fy = (1. - xFractBin(j)) * (1. - yFractBin(j))*f1 +
         yFractBin(j) * (1. - yFractBin(j))*f2 +
         (1. - xFractBin(j)) * yFractBin(j)*f4 +       
	 xFractBin(j) * yFractBin(j) * f3;
	
	//DEBUG PRINT
	/*if ( j < num) {
	    double x = mbs.local_particles(0,j);
	    double y = mbs.local_particles(2,j);
	    outf << std::endl << j << ": (" << x << "," << y << ") : ";
	    outf << " " << fx*tau*factor;
	    outf << " " << fy*tau*factor;

	}*/
	
	if (includeLocalDensity) {
	    // Include local line density effect
	    lambda = z_density(zBin(j)); 
	}
		    		
	mbs.local_particles(1,j) += lambda*tau*factor*fx;
	mbs.local_particles(3,j) += lambda*tau*factor*fy;

    }
    
    fftw_destroy_plan(forward_plan_greens);
    fftw_destroy_plan(forward_plan_charge);
    fftw_destroy_plan(backward_plan);
    
    logFile << "Calculation finished" << std::endl;

}


void TransverseSolver::kick_transverse_charge_direct(Macro_bunch_store &mbs)
{
    
    int i, index, iX, iY, j, nPoints;
    double x1, x2, y1, y2;
    double epsSq, rTransX, rTransY, rTot2, rr;
    int nMacros;
   
    // Form grid (2 times larger than real particle extent, check Hockney convolution method for details)
    
    logFile << "******************** DIRECT METHOD **************************************" << std::endl;
    
    logFile << "\nForming the grid...\n";
    
    form_grid(mbs);
    
    // Calculate perveance 
    /*
    double lambda = Real(mp._globalNMacros) * Injection::nReals_Macro 
              * 2. * Consts::pi * Ring::harmonicNumber / 
	       (Ring::lRing * (mp._phiMax - mp._phiMin) );
				   
    double perveance =  Sqr(mp._syncPart._charge)  *  _lambda * rClassical / 
	               (2. * Sqr(mp._syncPart._betaSync) * Cube(mp._syncPart._gammaSync)* mp._syncPart._mass); 
					     
    */
    
    // Bin the particles (set rho)
    
    logFile << "Binning the particles...\n";
    deposit_charge(mbs);
    dbg_print_array_2d("rho", rho);
     
    // Calculate forces
    
    fscx.reshape(nXBins+1, nYBins+1);
    fscy.reshape(nXBins+1, nYBins+1);
    
    fscx.set_all(0.);
    fscy.set_all(0.);
    
    for (iX = 0; iX <= nXBins; ++iX)
    {
	double a = xGrid(iX);
	for (iY = 0; iY <= nYBins; ++iY)
	{
	    double b = yGrid(iY);
	    logFile << "Calculating forces for point: (" << a << "," << b << ")" << std::endl;
	    for(int iXS = 0; iXS <= nXBins; iXS++)
	    {
		rTransX = a - xGrid(iXS);
		for(int iYS = 0; iYS <= nYBins; iYS++)
		{
		    rTransY = b - yGrid(iYS);
		    logFile << "\t...contribution of point (" << xGrid(iXS) << "," << yGrid(iYS) << ") is (";
		    logFile << rTransX << "," << rTransY << ")" << std::endl;
		    rTot2 = rTransX*rTransX + rTransY*rTransY + epsSq;
	            rr = rho(iXS, iYS) / rTot2;
		    
		    fscx(iX,iY) += rr * rTransX;
		    fscy(iX,iY) += rr * rTransY;
		}
	    }
	    logFile << "The force is : (" << fscx(iX,iY) << "," << fscy(iX,iY) << ")" << std::endl;
	
	}
    }
    
    //std::ofstream outfile("direct_results");
    std::ofstream& outfile = logFile;
    outfile << "RESULTS: " << std::endl;
    for (int i = 0; i <= nYBins; ++i) {
	for (int j = 0; j <= nXBins; ++j) {
	    outfile << "(" << fscx(i,j)  << "," << fscy(i,j) << ") ";
	}
	outfile << std::endl;
    }
	
    
/*
    // Include the perveance: (factor of 1.e3 to get F in MKS.
    //                           since our distances are in mm)

        Real factor = 4.0e3 * _perveance/
             (Real(mp._globalNMacros) * Real(nXBins) * Real(nYBins));

        for (iX = 1; iX <=nXBins; iX++)
          {
          for (iY = 1; iY <= nYBins; iY++)
            {
               index = iY-1 + nYBins * (iX-1);
               _fscx(iX, iY) = c_re(_out1[index]) * factor;
               _fscy(iX, iY) = c_im(_out1[index]) * factor;
	    }
          }

    */

    logFile << "Calculation finished" << std::endl;
}
 


   
