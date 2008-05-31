#include <boost/python.hpp>
#include "util1.h"
#include "RKIntegrators.h"

using namespace boost::python;

  BOOST_PYTHON_MODULE(ECloudPy)
  {
    // set the docstring of the current module scope
    scope().attr("__doc__") = "The ECloud C++ to Python interface.";
    def("initGslRan", &setGSLRandom);
    def("getLandauEnergyDist", &getLandauEnergyDist);
    def("getTechXEnergyDist", &getTechXEnergyDist);
    // The RK Integrator class
    class_<RKIntegrator>("RKIntegrator", init<bool>())
      .def("propagateF", &RKIntegrator::propagateFPy)
      .def("propagateV", &RKIntegrator::propagateVPy)
      .def("propThroughBunch", &RKIntegrator::propThroughBunchPy)
      .def("propBetweenBunches", &RKIntegrator::propBetweenBunchesPy)
      .def("setUnits", &RKIntegrator::setUnits)
      .def("closeTrajectoryFile", &RKIntegrator::closeTrajectoryFile)
      .def("reOpenTrajectoryFile", &RKIntegrator::reOpenTrajectoryFilePy)
      .def("setToRelativstic", &RKIntegrator::setToRelativstic)
      .def("isRelativistic",  &RKIntegrator::isRelativistic)
      .def("setPrecisionStep", &RKIntegrator::setPrecisionStep)
      .def("setTrajectoryFileName", &RKIntegrator::setTrajectoryFileNamePy)
      .def("getTrajectoryFileName", &RKIntegrator::getTrajectoryFileName)
      .def("setGamProtonBunch", &RKIntegrator::setGamProtonBunch)
      .def("getGamProtonBunch", &RKIntegrator::getGamProtonBunch)
      .def("setBFieldStaticCmp", &RKIntegrator::setBFieldStaticCmp)
      .def("getBFieldStaticCmp", &RKIntegrator::getBFieldStaticCmp)
      .def("setStepRatio", &RKIntegrator::setStepRatio)
      .def("getStepRatio", &RKIntegrator::getStepRatio)
      .def("setDynamicRelativistic", &RKIntegrator::setDynamicRelativistic)
      .def("getDynamicRelativistic", &RKIntegrator::getDynamicRelativistic)
      .def("setRelEFieldChange", &RKIntegrator::setRelEFieldChange)
      .def("getRelEFieldChange", &RKIntegrator::getRelEFieldChange)
      .def("setMaximumXBeamPipe", &RKIntegrator::setMaximumXBeamPipe)
      .def("getMaximumXBeamPipe", &RKIntegrator::getMaximumXBeamPipe)
      .def("setMaximumYBeamPipe", &RKIntegrator::setMaximumYBeamPipe)
      .def("getMaximumYBeamPipe", &RKIntegrator::getMaximumYBeamPipe)
      .def("reachedBeamPipe", &RKIntegrator::reachedBeamPipe)
      .def("resetClock", &RKIntegrator::resetClock)
      .def("getTime", &RKIntegrator::getTime)
      .def("setDebugOn", &RKIntegrator::setDebugOn)
      .def("setDebugOff", &RKIntegrator::setDebugOff)
      .def("gotPropagationError", &RKIntegrator::gotPropagationError)
    ;   
  }


