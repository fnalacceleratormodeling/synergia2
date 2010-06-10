#include <iostream>
#include <vector>
#include <cmath>
#include "VorpalElectrons.h"
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

const double  aVorpalElectron::MASS=510998.918;
const double  aVorpalElectron::speedOfLight=2.99792458e8;
const std::string VorpalElectrons::DATASET_NAME = std::string("electrons");
aVorpalElectron::aVorpalElectron():
x(0.),
betax (0.),
y(0.),
betay (0.),
z(0.),
betaz (0.),
gamVx(0.),
gamVy(0.),
gamVz(0.),
energy(MASS),
kineticEnergy(0.) 
 {
}

aVorpalElectron::aVorpalElectron(const double *ee):
x(ee[0]),
betax (0.),
y(ee[1]),
betay (0.),
z(ee[2]),
betaz (0.),
gamVx(ee[3]),
gamVy(ee[4]),
gamVz(ee[5])
 {
  double gamVN = std::sqrt(ee[3]*ee[3] + ee[4]*ee[4] + ee[5]*ee[5]); // in m/sec
  double gamBeta = gamVN/speedOfLight;
  double pMom = MASS*gamBeta;  // in eV
  energy = std::sqrt(MASS*MASS + pMom*pMom);
  kineticEnergy = energy-MASS;
  double gam = energy/MASS;
  betax = ee[3]/(gam*speedOfLight);
  betay = ee[4]/(gam*speedOfLight);
  betaz = ee[5]/(gam*speedOfLight);
}

VorpalElectrons::VorpalElectrons(const char *fName):
gotIt(false),
data(0),
time(0.),
stepNum(0)
 {

  //   debugIsOn = true;
   try
   {
      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
      Exception::dontPrint();
//      if (debugIsOn) std::cout << " VorpalElectrons constructor, on file " 
//                               << std::string(fName) << std::endl;
      /*
       * Open the specified file and the specified dataset in the file.
       */
      H5File file( fName, H5F_ACC_RDONLY );
      
      DataSet dataset = file.openDataSet( DATASET_NAME.c_str() );
      hsize_t dims_out[2];
      DataSpace dataspace = dataset.getSpace();
      int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
      size_t nTot = dims_out[0]*dims_out[1];
      std::vector<double> allData(nTot);
      dataset.read(&allData[0], PredType::NATIVE_DOUBLE);
      for (int k=0; k < nTot; k+=6) data.push_back(aVorpalElectron(&allData[k]));
      
      Group timeGrp = file.openGroup("time");
      int locIdTime =  timeGrp.getLocId(); 
      hid_t hid_time = H5Aopen_name(locIdTime, "time");
      H5Aread(hid_time, H5Aget_type(hid_time),  &time);
      hid_t hid_step = H5Aopen_name(locIdTime, "step");
      H5Aread(hid_step, H5Aget_type(hid_step),  &stepNum);
      
   }
   // catch failure caused by the H5File operations
   catch( FileIException error )
   {
//      error.printError();
      return ;
   }

   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
//      error.printError();
      return ;
   }

   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
//      error.printError();
      
      return ;
   }

   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
//      error.printError();
      return ;
   }
   gotIt=true;
}
