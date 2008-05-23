#include <iostream>
#include <fstream>
#include <cmath>
#include "biGaussFunc.h"

// BiGaussian, used for testing..
 
biGaussFunc::biGaussFunc(double *sigsIn, double dens):
density(dens),
sigsSqx(sigsIn[0]*sigsIn[0]),
sigsSqy(sigsIn[1]*sigsIn[1])
{
  
}

gsl_integration_workspace *biGaussFunc::workspace =  
          gsl_integration_workspace_alloc(biGaussFunc::maxLim);

double biGaussFunc::val(double x, double y) const {
  double parVal[4];
  parVal[0] = x; parVal[1] = y; parVal[2] = sigsSqx;  parVal[3] = sigsSqy; 
  double currentVal;
  double a = 0.;
  double b = 20.0*std::max(sigsSqx, sigsSqy);
  double absErr;
  gsl_function myGFunc;
  myGFunc.function = this->fbiGaussIntegralFunc;
  myGFunc.params = (void*) parVal;
  double absErrIn = 0.1*std::sqrt(sigsSqx*sigsSqy);
  int nn = gsl_integration_qag (&myGFunc, a, b, absErrIn, 1.0e-4, maxLim, 5, 
                       workspace, &currentVal, &absErr); 
//  std::cerr << " nn " << nn << " currentVal " 
//            << currentVal << " absErr " << absErr << std::endl;		       
  currentVal *= density;
  return currentVal;
}

double biGaussFunc::fbiGaussIntegralFunc(double t, void * p){
  // From H. Weideman book, ISBN 3-540-56550-7, p. 386, 
  double *pars = (double *) p; 
  double arg1 = pars[0]*pars[0]/(2.0*(pars[2] + t));
  double arg2 = pars[1]*pars[1]/(2.0*(pars[3] + t));
  double aa = 1.0 - std::exp(-1.0*(arg1+arg2));
  aa /= std::sqrt((pars[2] + t)*(pars[3] + t));
  return aa;
}
