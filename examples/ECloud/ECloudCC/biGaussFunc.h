#ifndef HAVE_ECLOUD_BIGAUSSFUNC_H
#define HAVE_ECLOUD_BIGAUSSFUNC_H
#include <cmath>
#include <gsl/gsl_integration.h>
class biGaussFunc {

  public:
    biGaussFunc(double *sigsIn, double dens);
    double density;
    double sigsSqx; double sigsSqy;
    
    double val(double x, double y) const;
    
  private:
    static  gsl_integration_workspace * workspace;
    static  double fbiGaussIntegralFunc(double t, void *params);
    static const size_t maxLim = 10000;
    double myParams[4];
};

#endif
