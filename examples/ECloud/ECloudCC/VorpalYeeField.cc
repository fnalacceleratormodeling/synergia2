#include <iostream>
#include <vector>
#include <cmath>
#include "VorpalYeeField.h"
#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

const std::string VorpalYeeField::DATASET_ENAME = std::string("YeeStaticElecFld");
const std::string VorpalYeeField::DATASET_BNAME = std::string("YeeStaticMagnFld");

VorpalYeeField::VorpalYeeField(const char* fNameField, bool isE, bool yC):
EFieldDefined(false),
BFieldDefined(false),
yeeCorr(yC), 
debugIsOn(false),
Ex(0),  Ey(0), Ez(0),
Bx(0),  By(0), Bz(0)
 {

  if (isE) constructEField(fNameField);
  else constructBField(fNameField);

}


VorpalYeeField::VorpalYeeField(const char* fNameEField, const char* fNameBField, bool yC):
EFieldDefined(false),
BFieldDefined(false),
yeeCorr(yC), 
debugIsOn(false),
Ex(0),  Ey(0), Ez(0),
Bx(0),  By(0), Bz(0)
{
  constructEField(fNameEField);
  constructBField(fNameBField);
}
VorpalYeeField::~VorpalYeeField() { 
  delete Ex; delete Ey; delete Ez;
  delete Bx; delete By; delete Bz;
}

void VorpalYeeField::addCstEField(const std::vector<double> eField) {
  int indexNd[3];
  std::vector<int> BNumPts= Ex->get_points().get_shape();  
  for (int i=0; i!=BNumPts[0]; i++) {
    indexNd[0]=i;
    for (int j=0; j!=BNumPts[1]; j++) {
      indexNd[1]=j;
      for (int k=0; k!=BNumPts[2]; k++) {
        indexNd[2]=k;
        double eG = Ex->get_points().get(indexNd);
	Ex->get_points().set(indexNd, eG+eField[0] );
        eG = Ey->get_points().get(indexNd);
	Ey->get_points().set(indexNd, eG+eField[1] );
        eG = Ez->get_points().get(indexNd);
	Ez->get_points().set(indexNd, eG+eField[2] );
      }
    }
  }
}

void VorpalYeeField::constructEField(const char* fName) {
  
//   debugIsOn = true;
   try
   {
      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
      Exception::dontPrint();
      if (debugIsOn) std::cout << " VorpalYeeField::constructEField, on file " 
                               << std::string(fName) << std::endl;
      /*
       * Open the specified file and the specified dataset in the file.
       */
      H5File file( fName, H5F_ACC_RDONLY );
       /*
      * Start by looking at the mesh. 
      */
      Group meshGrp = file.openGroup("/mesh");
//      std::cerr << " Number of object in group mesh " << meshGrp.getNumObjs() << std::endl;
      int locIdMesh =  meshGrp.getLocId(); 
//      std::cerr << " Loc Id for this group " << locIdMesh << std::endl;
      
//      std::cerr << "  Number of attributes in this group " 
//                << H5Aget_num_attrs(locIdMesh) << std::endl;
		
      hid_t hid_numPhysCells = H5Aopen_name(locIdMesh, "totalNumPhysCells");
//      std::cerr << " hid_numPhysCells " << hid_numPhysCells << std::endl;
      DataSpace numPhysCellsSpace =  H5Aget_space(hid_numPhysCells);
       hsize_t  nnDimsnumPhysCells = H5Aget_storage_size(hid_numPhysCells);
       int physCellDims[nnDimsnumPhysCells];
       H5Aread(hid_numPhysCells, H5Aget_type(hid_numPhysCells),  physCellDims);
       if (debugIsOn) std::cout << " Number of cells, X, Y , Z  " << physCellDims[0] << " / "  
                <<  physCellDims[1] << " / " << physCellDims[2] << std::endl;
      
      hid_t hid_numLBounds = H5Aopen_name(locIdMesh, "lowerBounds");
      double physCellLB[nnDimsnumPhysCells];
      H5Aread(hid_numLBounds, H5Aget_type(hid_numLBounds),  physCellLB);
      if (debugIsOn) std::cout << " Lower of cells, X, Y , Z  " << physCellLB[0] << " / "  
                <<  physCellLB[1] << " / " << physCellLB[2] << std::endl;
      
       hid_t hid_numUBounds = H5Aopen_name(locIdMesh, "upperBounds");
      double physCellUB[nnDimsnumPhysCells];
      H5Aread(hid_numUBounds, H5Aget_type(hid_numUBounds),  physCellUB);
      if (debugIsOn) std::cout << " Upper of cells, X, Y , Z  " << physCellUB[0] << " / "  
                <<  physCellUB[1] << " / " << physCellUB[2] << std::endl;
      
      //
      // Define the field.. 
      //
      double physicalSize[3];
      for (int k=0; k != 3; k++) physicalSize[k] = physCellUB[k] - physCellLB[k];
      double physicalOffset[3];
      
      for (int k=0; k != 3; k++) physicalOffset[k] = (physCellUB[k] + physCellLB[k])/2.;
      //
      // Before we define this field, check for consistency with the magnetic field 
      //
      if (BFieldDefined) {
        std::vector<int> BNumPts= Bx->get_points().get_shape();
        std::vector<double> BPhysSize = Bx->get_physical_size();
        std::vector<double> BPhysOffset = Bx->get_physical_offset();
        for (int k=0; k!=3; k++) {
          if (BNumPts[k] != physCellDims[k]) {
             std::cerr <<  
   " VorpalYeeField::constructEField Inconsistency between pre-defined Magnetic and Electric: " << std::endl;
             std::cerr << " Direction " << k << " Magnetic array size " 
	               << BNumPts[k] << " in Vorpal file " << physCellDims[k] << std::endl;
		       
	     std::cerr << " Stop re-defining Electric Field ... " <<std::endl;
	     return;
	  }
          if (std::abs(BPhysSize[k] - physicalSize[k]) > 1.0e-10)  {
             std::cerr <<  
   " VorpalYeeField::constructEField Inconsistency between Magnetic and Electric phys. size grid : " << std::endl;
             std::cerr << " Direction " << k << " Magnetic phys. size " 
	               << BPhysSize[k] << " in Vorpal file " << physicalSize[k] << std::endl;
		       
	     std::cerr << " Stop re-defining Electric Field ... " <<std::endl;
	     return;
	  }
          if (std::abs(BPhysOffset[k] - physicalOffset[k]) > 1.0e-10)  {
             std::cerr <<  
   " VorpalYeeField::constructEField Inconsistency between Magnetic and Electric phys. size grid : " << std::endl;
             std::cerr << " Direction " << k << " Magnetic phys. offset " 
	               << BPhysOffset[k] << " in Vorpal file " << physicalOffset[k] << std::endl;
		       
	     std::cerr << " Stop re-defining Electric Field ... " <<std::endl;
	     return;
	  }
        }
      }
      
      Ex = new Real_scalar_field(physCellDims, physicalSize, physicalOffset);
      Ey = new Real_scalar_field(physCellDims, physicalSize, physicalOffset);
      Ez = new Real_scalar_field(physCellDims, physicalSize, physicalOffset);
      Ex->get_points().zero_all(); Ey->get_points().zero_all(); Ez->get_points().zero_all();
      /*
       * Get dataspace of the dataset.
       */
      DataSet dataset = file.openDataSet( DATASET_ENAME.c_str() );
      DataSpace dataspace = dataset.getSpace();
      if (!dataspace.isSimple()) {
         std::cerr << " Data space for the EField in file " << 
	  fName << " is NOT simple .. " << std::endl;
         return;
      }

      /*
       * Get the number of dimensions in the dataspace.
       */
      int rank = dataspace.getSimpleExtentNdims();
      if (rank != 4) {
          std::cerr << " Rank associated to the HDF5 Data space for the EField in file " << 
	  fName << " is NOT 4 .. " << std::endl;
	  return;
      }
      size_t nTot=3;
      for (int k=0; k !=3; k++) nTot *= (size_t) physCellDims[k];
      
      std::vector<double> allData(nTot);
      /*
       * Read all the data in one statement... may be a waste of time, as one could map things out. 
      */
      if (debugIsOn) std::cerr << " ... Ntot " << nTot << std::endl;
      dataset.read(&allData[0], PredType::NATIVE_DOUBLE);
      
      if (debugIsOn) {
        std::cout << " Values at middle ";
	for (int kk=0; kk !=2; kk++) { 
	  for (int k=0; k!=6; k++) std::cout << " " << allData[nTot/2 + k + 21*42*3 + 42*3 ];
	  std::cout << std::endl;
	}
	int nnD = 0;
        for (int k=0; k < allData.size(); k +=3) {
	  if (std::abs(allData[k]) > 100.) {
	       std::cout << " Some Ex values at index " 
	                                            << k << " val " << allData[k] << std::endl;
	       nnD++; 
	  }
	  if (nnD > 10) break;
	}
      }
      int indexXShift, indexYShift, indexZShift;
      if (yeeCorr) {
        indexXShift = 3*physCellDims[1]*physCellDims[2];
        indexYShift = 3*physCellDims[2];
	indexZShift = 3;
      }
      // Now load the Nd_array.. The slow way.... 
      int indexNd[3];  for (int k=0; k!=3; k++) indexNd[k]=0;
      for (int k=0; k < allData.size(); k +=3) {
        double Exi = allData[k];
        double Eyi = allData[k+1];
        double Ezi = allData[k+2];
	if (yeeCorr) { 
	  if (indexNd[0] < physCellDims[0]/2) Exi = 0.5 * (allData[k] + allData[k+indexXShift]);
	  else Exi = 0.5 * (allData[k] + allData[k-indexXShift]);
	  if (indexNd[1] < physCellDims[1]/2) Eyi = 0.5 * (allData[k+1] + allData[k+1+indexYShift]);
	  else Eyi = 0.5 * (allData[k+1] + allData[k+1-indexYShift]);
	  if (indexNd[2] < physCellDims[2]/2) Ezi = 0.5 * (allData[k+2] + allData[k+2+indexZShift]);
	  else Ezi = 0.5 * (allData[k+2] + allData[k+2-indexZShift]);
	}
        Ex->get_points().set(indexNd, Exi );
        Ey->get_points().set(indexNd, Eyi);
        Ez->get_points().set(indexNd, Ezi);
	indexNd[2]++;
	if (indexNd[2] == physCellDims[2]) {
	  indexNd[2]=0;
	  indexNd[1]++;
//	  if (debugIsOn) std::cout <<  "... at k index 1 bumped " << k << " Index[1] " 
//	                           << indexNd[1] << " [0] " << indexNd[0] << std::endl;
	  if (indexNd[1] == physCellDims[1]) {
	    indexNd[1]=0;
	    indexNd[0]++;
//	    if (debugIsOn) std::cout <<  "... at k index 0 bumped " << k << " Index[1] " 
//	                           << indexNd[1] << " [0] " << indexNd[0] << std::endl;
	  }
	}
      }
      
      if (indexNd[0] != physCellDims[0]) {
        std::cerr << " VorpalYeeField::constructEField Inconsistent array size, indexNd: ";
	for (int k=0; k!=3; k++) std::cerr << indexNd[k] << " ";
	std::cerr << std::endl;
	std::cerr << " Fatal, quit ! " << std::endl;
	exit(2);  
      }
      EFieldDefined = true;
      // test... 
      if (debugIsOn) {
         std::vector<double> loc(3); loc[0] = 0.01, loc[1] = 0.000034; loc[2] = .0000007;
         std::cerr << " Ex val at 0.01, 0.000034, .0000007 =  " << Ex->get_val(loc) << std::endl;
      }
   }  // end of try block

   // catch failure caused by the H5File operations
   catch( FileIException error )
   {
      error.printError();
      return ;
   }

   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
      error.printError();
      return ;
   }

   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
      error.printError();
      
      return ;
   }

   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
      error.printError();
      return ;
   }

  EFieldDefined=true;
  return;  // successfully terminated
}
void VorpalYeeField::constructBField(const char* fName) {

//   debugIsOn = true;
   try
   {
      /*
       * Turn off the auto-printing when failure occurs so that we can
       * handle the errors appropriately
       */
      Exception::dontPrint();
      if (debugIsOn) std::cout << " VorpalYeeField::constructBField, on file " 
                               << std::string(fName) << std::endl;
      /*
       * Open the specified file and the specified dataset in the file.
       */
      H5File file( fName, H5F_ACC_RDONLY );
       /*
      * Start by looking at the mesh. 
      */
      Group meshGrp = file.openGroup("/mesh");
//      std::cerr << " Number of object in group mesh " << meshGrp.getNumObjs() << std::endl;
      int locIdMesh =  meshGrp.getLocId(); 
//      std::cerr << " Loc Id for this group " << locIdMesh << std::endl;
      
//      std::cerr << "  Number of attributes in this group " 
//                << H5Aget_num_attrs(locIdMesh) << std::endl;
		
      hid_t hid_numPhysCells = H5Aopen_name(locIdMesh, "totalNumPhysCells");
//      std::cerr << " hid_numPhysCells " << hid_numPhysCells << std::endl;
      DataSpace numPhysCellsSpace =  H5Aget_space(hid_numPhysCells);
       hsize_t  nnDimsnumPhysCells = H5Aget_storage_size(hid_numPhysCells);
       int physCellDims[nnDimsnumPhysCells];
       H5Aread(hid_numPhysCells, H5Aget_type(hid_numPhysCells),  physCellDims);
       if (debugIsOn) std::cout << " Number of cells, X, Y , Z  " << physCellDims[0] << " / "  
                <<  physCellDims[1] << " / " << physCellDims[2] << std::endl;
      
      hid_t hid_numLBounds = H5Aopen_name(locIdMesh, "lowerBounds");
      double physCellLB[nnDimsnumPhysCells];
      H5Aread(hid_numLBounds, H5Aget_type(hid_numLBounds),  physCellLB);
      if (debugIsOn) std::cout << " Lower of cells, X, Y , Z  " << physCellLB[0] << " / "  
                <<  physCellLB[1] << " / " << physCellLB[2] << std::endl;
      
       hid_t hid_numUBounds = H5Aopen_name(locIdMesh, "upperBounds");
      double physCellUB[nnDimsnumPhysCells];
      H5Aread(hid_numUBounds, H5Aget_type(hid_numUBounds),  physCellUB);
      if (debugIsOn) std::cout << " Upper of cells, X, Y , Z  " << physCellUB[0] << " / "  
                <<  physCellUB[1] << " / " << physCellUB[2] << std::endl;
      //
      // Define the field.. 
      //
      double physicalSize[3];
      for (int k=0; k != 3; k++) physicalSize[k] = physCellUB[k] - physCellLB[k];
      double physicalOffset[3];
      for (int k=0; k != 3; k++) physicalOffset[k] = (physCellUB[k] + physCellLB[k])/2.;
      //
      // Before we define this field, check for consistency with the magnetic field 
      //
      if (EFieldDefined) {
        std::vector<int> ENumPts= Ex->get_points().get_shape();
        std::vector<double> EPhysSize = Ex->get_physical_size();
        std::vector<double> EPhysOffset = Ex->get_physical_offset();
        for (int k=0; k!=3; k++) {
          if (ENumPts[k] != physCellDims[k]) {
             std::cerr <<  
   " VorpalYeeField::constructBField Inconsistency between pre-defined Magnetic and Electric: " << std::endl;
             std::cerr << " Direction " << k << " Electric array size " 
	               << ENumPts[k] << " in Vorpal file " << physCellDims[k] << std::endl;
		       
	     std::cerr << " Stop re-defining Magnetic Field ... " <<std::endl;
	     return;
	  }
          if (std::abs(EPhysSize[k] - physicalSize[k]) > 1.0e-10)  {
             std::cerr <<  
   " VorpalYeeField::constructBField Inconsistency between Magnetic and Electric phys. size grid : " << std::endl;
             std::cerr << " Direction " << k << " Magnetic phys. size " 
	               << EPhysSize[k] << " in Vorpal file " << physicalSize[k] << std::endl;
		       
	     std::cerr << " Stop re-defining Magnetic Field ... " <<std::endl;
	     return;
	  }
          if (std::abs(EPhysOffset[k] - physicalOffset[k]) > 1.0e-10)  {
             std::cerr <<  
   " VorpalYeeField::constructBField Inconsistency between Magnetic and Electric phys. size grid : " << std::endl;
             std::cerr << " Direction " << k << " Magnetic phys. offset " 
	               << EPhysOffset[k] << " in Vorpal file " << physicalOffset[k] << std::endl;
		       
	     std::cerr << " Stop re-defining Magnetic Field ... " <<std::endl;
	     return;
	  }
        }
      }
      
      Bx = new Real_scalar_field(physCellDims, physicalSize, physicalOffset);
      By = new Real_scalar_field(physCellDims, physicalSize, physicalOffset);
      Bz = new Real_scalar_field(physCellDims, physicalSize, physicalOffset);
      Bx->get_points().zero_all(); By->get_points().zero_all(); Bz->get_points().zero_all();
      /*
       * Get dataspace of the dataset.
       */
      DataSet dataset = file.openDataSet( DATASET_BNAME.c_str() );
      DataSpace dataspace = dataset.getSpace();
      if (!dataspace.isSimple()) {
         std::cerr << " Data space for the BField in file " << 
	  fName << " is NOT simple .. " << std::endl;
         return;
      }

      /*
       * Get the number of dimensions in the dataspace.
       */
      int rank = dataspace.getSimpleExtentNdims();
      if (rank != 4) {
          std::cerr << " Rank associated to the HDF5 Data space for the BField in file " << 
	  fName << " is NOT 4 .. " << std::endl;
	  return;
      }
      size_t nTot=3;
      for (int k=0; k !=3; k++) nTot *= (size_t) physCellDims[k];
      
      std::vector<double> allData(nTot);
      /*
       * Read all the data in one statement... may be a waste of time, as one could map things out. 
      */
      if (debugIsOn) std::cerr << " ... Ntot " << nTot << std::endl;
      dataset.read(&allData[0], PredType::NATIVE_DOUBLE);
      
      if (debugIsOn) {
        std::cout << " Values at middle ";
	for (int kk=0; kk !=2; kk++) { 
	  for (int k=0; k!=6; k++) std::cout << " " << allData[nTot/2 + k + 21*42*3 + 42*3 ];
	  std::cout << std::endl;
	}
	int nnD = 0;
        for (int k=0; k < allData.size(); k +=3) {
	  if (std::abs(allData[k]) > 100.) {
	       std::cout << " Some Ex values at index " 
	                                            << k << " val " << allData[k] << std::endl;
	       nnD++; 
	  }
	  if (nnD > 10) break;
	}
      }
// Not needed for the BField.. I think.... 
//      int indexXShift, indexYShift, indexZShift;
//      if (yeeCorr) {
//        indexXShift = 3*physCellDims[1]*physCellDims[2];
//        indexYShift = 3*physCellDims[2];
//	indexZShift = 3;
//      }
      // Now load the Nd_array.. The slow way.... 
      int indexNd[3];  for (int k=0; k!=3; k++) indexNd[k]=0;
      for (int k=0; k < allData.size(); k +=3) {
        double Bxi = allData[k];
        double Byi = allData[k+1];
        double Bzi = allData[k+2];
//	if (yeeCorr) { 
//	  if (indexNd[0] < physCellDims[0]/2) Exi = 0.5 * (allData[k] + allData[k+indexXShift]);
//	  else Exi = 0.5 * (allData[k] + allData[k-indexXShift]);
//	  if (indexNd[1] < physCellDims[1]/2) Eyi = 0.5 * (allData[k+1] + allData[k+1+indexYShift]);
//	  else Eyi = 0.5 * (allData[k+1] + allData[k+1-indexYShift]);
//	  if (indexNd[2] < physCellDims[2]/2) Ezi = 0.5 * (allData[k+2] + allData[k+2+indexZShift]);
//	  else Ezi = 0.5 * (allData[k+2] + allData[k+2-indexZShift]);
//	}
        Bx->get_points().set(indexNd, Bxi );
        By->get_points().set(indexNd, Byi);
        Bz->get_points().set(indexNd, Bzi);
	indexNd[2]++;
	if (indexNd[2] == physCellDims[2]) {
	  indexNd[2]=0;
	  indexNd[1]++;
//	  if (debugIsOn) std::cout <<  "... at k index 1 bumped " << k << " Index[1] " 
//	                           << indexNd[1] << " [0] " << indexNd[0] << std::endl;
	  if (indexNd[1] == physCellDims[1]) {
	    indexNd[1]=0;
	    indexNd[0]++;
//	    if (debugIsOn) std::cout <<  "... at k index 0 bumped " << k << " Index[1] " 
//	                           << indexNd[1] << " [0] " << indexNd[0] << std::endl;
	  }
	}
      }
      
      if (indexNd[0] != physCellDims[0]) {
        std::cerr << " VorpalYeeField::constructBField Inconsistent array size, indexNd: ";
	for (int k=0; k!=3; k++) std::cerr << indexNd[k] << " ";
	std::cerr << std::endl;
	std::cerr << " Fatal, quit ! " << std::endl;
	exit(2);  
      }
      BFieldDefined = true;
      // test... 
      if (debugIsOn) {
         std::vector<double> loc(3); loc[0] = 0.01, loc[1] = 0.000034; loc[2] = .0000007;
         std::cerr << " Ex val at 0.01, 0.000034, .0000007 =  " << Ex->get_val(loc) << std::endl;
      }
   }  // end of try block

   // catch failure caused by the H5File operations
   catch( FileIException error )
   {
      error.printError();
      return ;
   }

   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
      error.printError();
      return ;
   }

   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
      error.printError();
      
      return ;
   }

   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
      error.printError();
      return ;
   }

  BFieldDefined=true;
  return;  // successfully terminated

}

std::vector<int> VorpalYeeField::getDims() const {
  if (EFieldDefined) {
    return Ex->get_points().get_shape();
  } else if (BFieldDefined) {
    return Bx->get_points().get_shape();
  } else { 
    std::vector<int> empty;
    return empty; 
  } 
} 

std::vector<double> VorpalYeeField::getPhysSize() const  {
  if (EFieldDefined) {
    return Ex->get_physical_size();
  } else if (BFieldDefined) {
    return Bx->get_physical_size();
  } else { 
    std::vector<double> empty;
    return empty; 
  } 
} 

std::vector<double> VorpalYeeField::getPhysOffset() const  {
  if (EFieldDefined) {
    return Ex->get_physical_offset();
  } else if (BFieldDefined) {
    return Bx->get_physical_offset();
  } else { 
    std::vector<double> empty;
    return empty; 
  } 
}

double VorpalYeeField::getEField(const std::vector<double> location, int iComp) const { // x=0, y=1, z=2

   switch (iComp) {
     case 0: 
       return Ex->get_val(location);
     case 1: 
       return Ey->get_val(location);
     case 2: 
       return Ez->get_val(location);
     default:
        std::cerr << " VorpalYeeField::getEField Incorrect component index " << iComp << std::endl;
       return 0.;
  }
  return 0.; // should never be executed.   
}

double VorpalYeeField::getBField(const std::vector<double> location, int iComp) const { // x=0, y=1, z=2

   switch (iComp) {
     case 0: 
       return Bx->get_val(location);
     case 1: 
       return By->get_val(location);
     case 2: 
       return Bz->get_val(location);
     default:
        std::cerr << " VorpalYeeField::getEField Incorrect component index " << iComp << std::endl;
       return 0.;
  }
  return 0.; // should never be executed.   
}
