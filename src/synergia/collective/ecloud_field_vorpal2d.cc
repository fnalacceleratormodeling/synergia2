#include <iostream>
#include <sstream>
#include <fstream>
#include "ecloud_field_vorpal2d.h"
#include <boost/serialization/vector.hpp>
#include <stdexcept>
#include "mpi.h"

ECloudEFieldVORPAL2D::YScanAtX::YScanAtX():
x(0.),
coefChebEX(0),
coefChebEY(0)
{ ; }

ECloudEFieldVORPAL2D::YXScanAtdZ::YXScanAtdZ():
dz(0.),
data(0)
{
;
}

ECloudEFieldVORPAL2D::ECloudEFieldVORPAL2D():
// myComm() // default constructor, takes over the whole MPI world
version("1.0"),
VORPALJobName("none"),
verticalChebychevOrder(16),
yLow(-2.2),
yUp(2.2),
data(0),
tCheb(0)
{
  tCheb = gsl_cheb_alloc ((size_t) verticalChebychevOrder);
  free(tCheb->c); tCheb->c=0; // We will always point to our std::vector
  tCheb->a = yLow;
  tCheb->b = yUp;
//  std::cerr << " tCheb created in blank constructor .. " 
//            << tCheb->order <<  " Ptr " << tCheb << " Ptr f " << tCheb->f << std::endl;
  if (tCheb->f != 0) free(tCheb->f); tCheb->f = 0;
	    
 }
  
ECloudEFieldVORPAL2D::ECloudEFieldVORPAL2D(const char *archiveName):
// myComm() // default constructor, takes over the whole MPI world
version("1.0"),
VORPALJobName("none"),
verticalChebychevOrder(16),
yLow(-2.2),
yUp(2.2),
data(0),
tCheb(0)
{
   MPI_Barrier(myComm_sptr->get()); // Might not be needed.
   int my_rank= myComm_sptr->get_rank();
   if (my_rank == 0)  loadFromFile(archiveName); 
   int numProc=1;
   MPI_Comm_size(myComm_sptr->get(), &numProc);     
   if (numProc > 1) broadcastIt();
     
}
ECloudEFieldVORPAL2D::ECloudEFieldVORPAL2D(Commxx_sptr comm_sptr, const char *archiveName):
myComm_sptr(comm_sptr),
version("1.0"),
VORPALJobName("none"),
verticalChebychevOrder(16),
yLow(-2.2),
yUp(2.2),
data(0),
tCheb(0)
{
   MPI_Barrier(myComm_sptr->get()); // Might not be needed.
   int my_rank= myComm_sptr->get_rank();
   if (my_rank == 0)  loadFromFile(archiveName); 
   int numProc=1;
   MPI_Comm_size(myComm_sptr->get(), &numProc);     
   if (numProc > 1) broadcastIt();
     
}    
ECloudEFieldVORPAL2D::ECloudEFieldVORPAL2D( const ECloudEFieldVORPAL2D &c) {
     myComm_sptr = c.getComm(); // The communicator 
     version = c.getVersion();  // I'll probably change my mind, so there it is..
     VORPALJobName = c.getVORPALJobName(); // strictly ofr experts..
     verticalChebychevOrder = c.getVerticalChebychevOrder(); // the order for the Chebyshev
     yLow = c.getYLow(); // the lower limit of the map
     yUp = c.getYUp();  // the upper limit of the map.
     data = c.getData();  // vector along z. 
     if (tCheb != 0) {
       tCheb->c=0; // Pointing to our vector 
       if (tCheb->f != 0) free(tCheb->f); tCheb->f = 0;
       free(tCheb);
       tCheb=0;
     }
     tCheb = gsl_cheb_alloc ((size_t) verticalChebychevOrder);
//  std::cerr << " verticalChebychevOrder from file " << verticalChebychevOrder << std::endl;
//  std::cerr << " Freeing coefficients after upload file, with order   .. " << tCheb->order << std::endl;
     free(tCheb->c); tCheb->c=0; // We will always point to our std::vector
     if (tCheb->f != 0) free(tCheb->f); tCheb->f = 0;
//  std::cerr << " Done  .. " << tCheb->order << std::endl;
     tCheb->a = yLow;
     tCheb->b = yUp;
}
ECloudEFieldVORPAL2D& ECloudEFieldVORPAL2D::operator=(const ECloudEFieldVORPAL2D &c) {

     myComm_sptr = c.getComm(); // The communicator 
     version = c.getVersion();  // I'll probably change my mind, so there it is..
     VORPALJobName = c.getVORPALJobName(); // strictly ofr experts..
     verticalChebychevOrder = c.getVerticalChebychevOrder(); // the order for the Chebyshev
     yLow = c.getYLow(); // the lower limit of the map
     yUp = c.getYUp();  // the upper limit of the map.
     data = c.getData();  // vector along z. 
     if (tCheb != 0) {
       tCheb->c=0; // Pointing to our vector 
       if (tCheb->f != 0) free(tCheb->f); tCheb->f = 0;
       free(tCheb);
       tCheb=0;
     }
     tCheb = gsl_cheb_alloc ((size_t) verticalChebychevOrder);
//  std::cerr << " verticalChebychevOrder from file " << verticalChebychevOrder << std::endl;
//  std::cerr << " Freeing coefficients after upload file, with order   .. " << tCheb->order << std::endl;
     free(tCheb->c); tCheb->c=0; // We will always point to our std::vector
     if (tCheb->f != 0) free(tCheb->f); tCheb->f = 0;
//  std::cerr << " Done  .. " << tCheb->order << std::endl;
     tCheb->a = yLow;
     tCheb->b = yUp;
}

ECloudEFieldVORPAL2D::~ECloudEFieldVORPAL2D() {
  if (tCheb != 0) {
     tCheb->c=0; // Pointing to our vector 
     if (tCheb->f != 0) free(tCheb->f); tCheb->f = 0;
     free(tCheb);
     tCheb=0;
  }
}

bool ECloudEFieldVORPAL2D::loadFromFile(const char *archiveName) {

  if (tCheb != 0) {
      tCheb->c = 0;
      if (tCheb->f != 0) free(tCheb->f); tCheb->f = 0;
      free(tCheb);
      tCheb=0;
  }
  std::ifstream ifs(archiveName, std::ios::binary);
// Error handling to be done, always return true. 
  boost::archive::binary_iarchive ia(ifs);
  ia >> (*this);
  tCheb = gsl_cheb_alloc ((size_t) verticalChebychevOrder);
//  std::cerr << " verticalChebychevOrder from file " << verticalChebychevOrder << std::endl;
//  std::cerr << " Freeing coefficients after upload file, with order   .. " << tCheb->order << std::endl;
  free(tCheb->c); tCheb->c=0; // We will always point to our std::vector
  if (tCheb->f != 0) free(tCheb->f); tCheb->f = 0;
//  std::cerr << " Done  .. " << tCheb->order << std::endl;
  tCheb->a = yLow;
  tCheb->b = yUp;
  return true;
}  

void ECloudEFieldVORPAL2D::archiveIt(const char *archiveName) {

  std::ofstream ofs(archiveName, std::ios::binary);
  boost::archive::binary_oarchive oa(ofs);
  oa << (*this);

}

double ECloudEFieldVORPAL2D::GetFieldEX(double x, double y, double dz) const {

  if (data.size() == 0) return 0.; // dah.. 
  // Locate the relevant slice in dz. Not the fatest algorithm, but general..
  std::vector<YXScanAtdZ>::const_iterator itZLow=data.begin();
  if (dz < itZLow->dz) return 0.; // Outside. head end.. 
  std::vector<YXScanAtdZ>::const_iterator itZUp=data.end();
  itZUp--; if (dz > itZUp->dz) return 0.;
  // Assume the are sorted... 
  for(std::vector<YXScanAtdZ>::const_iterator itZ=data.begin(); itZ!=data.end(); itZLow=itZ, itZ++) {
    if (dz > itZ->dz) continue;
    itZUp = itZ;  
    break;  
  }
  const double ddz=itZUp->dz - itZLow->dz;
  const double udz = (dz - itZLow->dz)/ddz;
  double ExLow = GetFieldEAtdz(true, itZLow, x, y);
  double ExUp = GetFieldEAtdz(true, itZUp, x, y); 
  return ExLow*(1.0-udz) + ExUp*udz; 
}

double ECloudEFieldVORPAL2D::GetFieldEY(double x, double y, double dz) const {

//  std::cerr << " ECloudEFieldVORPAL2D::GetFieldEY, Num dz Slices " << data.size() << std::endl;
  if (data.size() == 0) return 0.; // dah.. 
  // Locate the relevant slice in dz. Not the fatest algorithm, but general..
  std::vector<YXScanAtdZ>::const_iterator itZLow=data.begin();
//  std::cerr << " Min dZ " <<  itZLow->dz << " Num dz Slices " << data.size() << std::endl;
  if (dz < itZLow->dz) return 0.; // Outside. head end.. 
  std::vector<YXScanAtdZ>::const_iterator itZUp=data.end();
  itZUp--; if (dz > itZUp->dz) return 0.;
  // Assume the are sorted... 
  for(std::vector<YXScanAtdZ>::const_iterator itZ=data.begin(); itZ!=data.end(); itZLow=itZ, itZ++) {
    if (dz > itZ->dz) continue;
    itZUp = itZ;
    break;  
  }
//  std::cerr << " for dz " << dz << " between " << itZLow->dz << " and " << itZUp->dz << std::endl;
  const double ddz=itZUp->dz - itZLow->dz;
  const double udz = (dz - itZLow->dz)/ddz;
  double EyLow = GetFieldEAtdz(false, itZLow, x, y);
  double EyUp = GetFieldEAtdz(false, itZUp,  x, y); 
  return EyLow*(1.0-udz) + EyUp*udz; 
}

double ECloudEFieldVORPAL2D::GetFieldEAtdz(bool isX, std::vector<YXScanAtdZ>::const_iterator itZ, double x, double y) const {

  std::vector<YScanAtX>::const_iterator itXLow=itZ->data.begin();
  if (x < itXLow->x) return 0.; // Outside. head end.. 
  std::vector<YScanAtX>::const_iterator itXUp=itZ->data.end();
  itXUp--;
//  for (std::vector<YScanAtX>::const_iterator itX=itZ->data.begin(); itX != itZ->data.end(); itX++) 
//      std::cerr << " xb " << itX->x << std::endl;
  // Assume uniform grid. (dictated by VORPAL )
  if (x > itXUp->x) return 0.;
  const double ldx = (itXUp->x - itXLow->x)/(itZ->data.size()-1);
  size_t ix = (int) ((x - itXLow->x)/ldx);
  itXLow += ix; 
  itXUp = itXLow; itXUp++;
//  std::cerr << " at cell x " << itXLow->x << " / " << itXUp->x << " ix " << ix <<  std::endl;
  const double udx = (x - itXLow->x)/ldx;
  double ELow = GetFieldEAtdzX(isX, itXLow, y);
  double EUp = GetFieldEAtdzX(isX, itXUp, y); 
  return ELow*(1.0-udx) + EUp*udx; 
}
  
double ECloudEFieldVORPAL2D::GetFieldEAtdzX(bool isX, std::vector<YScanAtX>::const_iterator itXZ, double y) const {
   
//   std::cerr << " ECloudEFieldVORPAL2D::GetFieldEAtdzX  y =" << y << std::endl;
   if (y < this->yLow) return 0.;
   if (y > this->yUp) return 0.;
   if (isX) {
     tCheb->c = (double *) &(itXZ->coefChebEX[0]);
   } else {
     tCheb->c = (double *) &(itXZ->coefChebEY[0]);
   }
//   double *pTmp = tCheb->c;
//   for (int k=0; k!=5; k++) pTmp++;
//   if (isX) std::cerr << " tCheb->a " << tCheb->a << " b " << tCheb->b << " orders " 
//             << tCheb->order << "/" << tCheb->order_sp << ", tCheb-5 " << itXZ->coefChebEX[5] << std::endl;
//   else std::cerr << " tCheb->a " << tCheb->a << " b " << tCheb->b << " orders " 
//             << tCheb->order << "/" << tCheb->order_sp << ", tCheb-5 " << itXZ->coefChebEY[5] << std::endl;
   const double res = gsl_cheb_eval(tCheb, y);	
//   std::cerr << " result " << res << std::endl;     	     
   return res;
}  

void ECloudEFieldVORPAL2D::cleanup() { //to save memory. 

   for (std::vector<YXScanAtdZ>::iterator itz=data.begin(); itz!=data.end(); itz++) 
     itz->data.clear();
   data.clear();
   // Leave the damm gsl pointer as is, done in loadArchive
}
 
bool ECloudEFieldVORPAL2D::loadSynergia(const char *archiveName, 
                      boost::shared_ptr<Rectangular_grid> sGridEx, boost::shared_ptr<Rectangular_grid> sGridEy) {
		      
   this->cleanup();
   this->loadFromFile(archiveName);
//   std::cerr << " Returned from loading from file, name " << VORPALJobName << " num ZSlice " << data.size() << std::endl;
   return this->loadSynergia(sGridEx, sGridEy);
}
  
bool ECloudEFieldVORPAL2D::loadSynergia( 
                      boost::shared_ptr<Rectangular_grid> sGridEx, boost::shared_ptr<Rectangular_grid> sGridEy) {
   
   std::cerr << " Loading Synergia Grid " << std::endl; 
   {
     // Ex grid 
     const Rectangular_grid_domain_sptr	sGridExD = sGridEx->get_domain_sptr();
     const std::vector<int> shapeEx = sGridExD->get_grid_shape();
     const std::vector<double> sizeEx = sGridExD->get_physical_size();
     const std::vector<double> offsetEx = sGridExD->get_physical_offset();
     MArray3d_ref a(sGridEx->get_grid_points());
     for (int ix=0; ix!=shapeEx[0]; ix++) {
       const double x = offsetEx[0] - sizeEx[0]/2. + ((double) ix + 0.5)*sizeEx[0]/shapeEx[0];
       for (int iy=0; iy!=shapeEx[1]; iy++) {
         const double y = offsetEx[1] - sizeEx[1]/2. + ((double) iy + 0.5)*sizeEx[1]/shapeEx[1];
         for (int iz=0; iz!=shapeEx[2]; iz++) {
           const double z = offsetEx[2] - sizeEx[2]/2. + ((double) iz + 0.5)*sizeEx[2]/shapeEx[2];
           const double Ex = this->GetFieldEX(x, y, z);
//	   std::cerr << " At ix = " << ix << " iy " << iy << " z " << z << " iz " << iz << " Ex " << Ex << std::endl;
	   a[ix][iy][iz] = Ex;
         }
       }
     }
   }
   {
     // Ey grid 
     const Rectangular_grid_domain_sptr	sGridEyD = sGridEy->get_domain_sptr();
     const std::vector<int> shapeEy = sGridEyD->get_grid_shape();
     const std::vector<double> sizeEy = sGridEyD->get_physical_size();
     const std::vector<double> offsetEy = sGridEyD->get_physical_offset();
     MArray3d_ref a(sGridEy->get_grid_points());
     for (int ix=0; ix!=shapeEy[0]; ix++) {
       const double x = offsetEy[0] - sizeEy[0]/2.  + ((double) ix + 0.5)*sizeEy[0]/shapeEy[0];
       for (int iy=0; iy!=shapeEy[1]; iy++) {
         const double y = offsetEy[1]  - sizeEy[1]/2. + ((double) iy + 0.5)*sizeEy[1]/shapeEy[1];
         for (int iz=0; iz!=shapeEy[2]; iz++) {
           const double z = offsetEy[2] - sizeEy[2]/2.  + ((double) iz + 0.5)*sizeEy[2]/shapeEy[2];
           const double Ey = this->GetFieldEY(x, y, z);
//	   std::cerr << " At x = " << x << " y " << y << " z " << z << " iz " << iz << " Ey " << Ey << std::endl;
	   a[ix][iy][iz] = Ey;
         }
       }
     }
   }
   return true; 
}

void ECloudEFieldVORPAL2D::loadOneScan(bool isEX, double dz, double x, double yLow, double yUp, 
                      size_t chebOrder,  gsl_cheb_series *chebSerie) {
		      
  std::vector<YXScanAtdZ>::iterator itzSel = data.end();
  for (std::vector<YXScanAtdZ>::iterator itz=data.begin(); itz!=data.end(); itz++) 
      if (std::abs(dz - itz->dz) < 1.0e-12) {itzSel = itz; break;}
  if (itzSel != data.end()) {
    if ((std::abs(yLow - chebSerie->a) > 1.0e-12) || (std::abs(yUp - chebSerie->b) > 1.0e-12)) {
        std::ostringstream ostr; 
	ostr << "ECloudEFieldVORPAL2D::loadOneScan Chebychev ranges must consistant over the entire grid \n" <<
	        " Existing range " << yLow << " to " << yUp 
	      << " submited " << chebSerie->a << " to " << chebSerie->b << std::endl;
	std::string ostrstr(ostr.str());      
        throw std::runtime_error(ostrstr.c_str());
    }
  } else {
    yLow = chebSerie->a; yUp = chebSerie->b;
  }   
  YScanAtX tDat;
  tDat.x = x;
  int n = chebSerie->order + 1;
  double *cE = chebSerie->c;
  if (isEX) {
    tDat.coefChebEX.resize(n);
    for (size_t i=0; i != n; i++, cE++) tDat.coefChebEX[i] = *cE;
  } else {
    tDat.coefChebEY.resize(n);
    for (size_t i=0; i != n; i++, cE++) tDat.coefChebEY[i] = *cE;
  }   
  // Look if Y Scan at fixed X has not been already loaded. 
  if (itzSel != data.end()) {
    for(std::vector<YScanAtX>::iterator itx = itzSel->data.begin(); itx != itzSel->data.end(); itx++) {
      if (std::abs(itx->x - x) < 1.0e-24) { // dangerous, but will work for most Synergia problems.
        if (isEX) itx->coefChebEX = tDat.coefChebEX; // Deep copy of stl vector should take place 
	else itx->coefChebEY = tDat.coefChebEY;
	return;
      }
    }
    itzSel->data.push_back(tDat);
  } else {
   YXScanAtdZ tYX;
   tYX.dz = dz;
   tYX.data.push_back(tDat);
   this->data.push_back(tYX);
  }
  
  return;
		      
}		      
void ECloudEFieldVORPAL2D::broadcastIt() {
   
   int my_rank= myComm_sptr->get_rank();
//   if (my_rank == 0) std::cerr << " got to broadcast.... " << std::endl;
   
   // We do not broadacast the name of the VORPAL job, version, since it should not be used in the calculation 
   if (my_rank != 0) { version = std::string("na"); VORPALJobName = version; }
   MPI_Bcast((void *) &verticalChebychevOrder, 1, MPI_INTEGER, 0, myComm_sptr->get());
   MPI_Bcast((void *) &yLow, 1, MPI_DOUBLE, 0, myComm_sptr->get());
   MPI_Bcast((void *) &yUp, 1, MPI_DOUBLE, 0, myComm_sptr->get());
   int numYXScanAtdZ = 0;
   if (my_rank == 0) numYXScanAtdZ = data.size();
   MPI_Bcast((void *) &numYXScanAtdZ, 1, MPI_INTEGER, 0, myComm_sptr->get());
   if (my_rank != 0) data.clear();
   for (size_t iYX = 0; iYX != numYXScanAtdZ; ++iYX) {
      if (my_rank != 0) {
        YXScanAtdZ aYxScan;
	data.push_back(aYxScan);
      }
      MPI_Bcast((void *) &(data[iYX].dz), 1, MPI_DOUBLE, 0, myComm_sptr->get());
      int numYScanAtX = 0;
      if (my_rank == 0) numYScanAtX = data[iYX].data.size();
      MPI_Bcast((void *) &numYScanAtX, 1, MPI_INTEGER, 0, myComm_sptr->get());
      if (my_rank != 0) data[iYX].data.clear();
      for (size_t iY = 0; iY != numYScanAtX; ++iY) {
        if (my_rank != 0) {
          YScanAtX aYScan;
	  data[iYX].data.push_back(aYScan);
	}
        MPI_Bcast((void *) &(data[iYX].data[iY].x), 1, MPI_DOUBLE, 0, myComm_sptr->get());
	int nOrderEx = 0;
	int nOrderEy = 0;
	if (my_rank == 0) {
	    nOrderEx = data[iYX].data[iY].coefChebEX.size();
	    nOrderEy = data[iYX].data[iY].coefChebEY.size();
	}
        MPI_Bcast((void *) &nOrderEx, 1, MPI_INTEGER, 0, myComm_sptr->get());
        MPI_Bcast((void *) &nOrderEy, 1, MPI_INTEGER, 0, myComm_sptr->get());
        if (my_rank != 0) {
	  data[iYX].data[iY].coefChebEX.resize(nOrderEx);
	  data[iYX].data[iY].coefChebEY.resize(nOrderEy);
	} 
        MPI_Bcast((void *) &(data[iYX].data[iY].coefChebEX[0]), nOrderEx, MPI_DOUBLE, 0, myComm_sptr->get());
        MPI_Bcast((void *) &(data[iYX].data[iY].coefChebEY[0]), nOrderEy, MPI_DOUBLE, 0, myComm_sptr->get());
     } // on Y scans
   } // on YX scans 
   // prepare for using GSL Cheb polynoms. 
   //
   if (my_rank != 0) {
       if (tCheb != 0) {
         tCheb->c=0; // Pointing to our vector 
         if (tCheb->f != 0) free(tCheb->f); tCheb->f = 0;
         free(tCheb);
       }
       tCheb = gsl_cheb_alloc ((size_t) verticalChebychevOrder);
       free(tCheb->c); tCheb->c=0; // We will always point to our std::vector
       tCheb->a = yLow;
       tCheb->b = yUp;
//  std::cerr << " tCheb created in blank constructor .. " 
//            << tCheb->order <<  " Ptr " << tCheb << " Ptr f " << tCheb->f << std::endl;
       if (tCheb->f != 0) free(tCheb->f); tCheb->f = 0;
    }

}
		      
