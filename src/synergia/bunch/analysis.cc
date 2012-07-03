#include "analysis.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <fftw3.h>
#include <fstream>
#include <cfloat>

Analysis::Analysis(std::string const& fname, size_t maxTurn):
    maxNumberOfTurn(maxTurn),
    seqCount(0),
    tokenname(fname),
    filename(fname),
    numTurns(0),
    numParticles1rstBunch(0),
    selectedParticleId(0),
    gotHorizontalTunes(false),
    gotVerticalTunes(false),
    numXTunesFound(0),
    XTunes(boost::extents[1]),
    numYTunesFound(0),
    YTunes(boost::extents[1]),
    minimumRequiredTurnNum(50),
    coords(boost::extents[maxNumberOfTurn][1][2])
{
   // Load the coordinates, all of it, hoping we have enough memory. 
   uploadCoords();
   XTunes.resize(boost::extents[coords.shape()[1]]);
   YTunes.resize(boost::extents[coords.shape()[1]]);
   compute_betatron_tunes(true);   
   compute_betatron_tunes(false);   
}

Analysis::~Analysis() 
{ 
//  std::cerr << " Trying to close the file .. " << file.getFileName() << std::endl;
//  file.close();
//  std::cerr << " ...... closed it ....  " << std::endl;
}

void Analysis::uploadCoords() {

  numTurns=0;
  hsize_t      offset[2]={0,0};   // hyperslab offset in the file
  hsize_t      dims_out[2]={0,0};   // the whole thing.
  if (maxNumberOfTurn < minimumRequiredTurnNum) {
    std::cerr << " Not enough turn to measure decent tunes, do nothing " << std::endl;
    return;
  }
  while (numTurns < maxNumberOfTurn) { 
     std::ostringstream sstream;
     sstream << tokenname << "_" << std::setw(4) << std::setfill('0') << numTurns << ".h5";
     std::string fullname(sstream.str());
     // Check if the file exists..
     std::ifstream inFileTmp;
     inFileTmp.open(fullname.c_str());
     if (!inFileTmp.is_open()) {
         break; // We do not increment the turn number, ad stop, no further turns. 
     } 
     inFileTmp.close();
     inFileTmp.clear();
     H5::H5File file=H5::H5File(fullname.c_str(), H5F_ACC_RDONLY);
     H5::DataSet bunch_ds = file.openDataSet("particles");
     H5::FloatType datatypeDble(H5::PredType::NATIVE_DOUBLE);
     H5::DataSpace spaceBunch = bunch_ds.getSpace();
     size_t nTotalEntries = (size_t) dims_out[0];
//     std::cerr << " Number of particle at turn " << numTurns << " is " << nTotalEntries << std::endl;
     hsize_t      countVect[2]={1,7};    // size of the hyperslab in the file
     H5::DataSpace  memspace(2, countVect);
     MArray1d aPartCoord(boost::extents[7]); // coordinate values.
     if (coords.shape()[1] == 1)  {
         coords.resize(boost::extents[maxNumberOfTurn][nTotalEntries][2]);
         for(size_t i=0; i!=maxNumberOfTurn; i++) { 
	   for(size_t j=0; j!=nTotalEntries; j++) { 
	     coords[i][j][0] = ((double) FLT_MAX) ;
	     coords[i][j][1] = ((double) FLT_MAX) ;
	   }
	 }
     }
     size_t nEntry=0;
     int nGoodParticles=0;
     while (nEntry < nTotalEntries) { // Slow and tedious get every X or y position, 
       offset[0] = nEntry;
       spaceBunch.selectHyperslab( H5S_SELECT_SET, countVect, offset );
       bunch_ds.read(&aPartCoord[0], datatypeDble, memspace, spaceBunch );
       size_t partNum = (size_t) ( (int) aPartCoord[Bunch::id]);
       if (partNum >= nTotalEntries) {nEntry++; continue; } // Should not happen, but depend on algorithms to define 
                                                // Selected particles. 
       coords[numTurns][partNum][0] = aPartCoord[Bunch::x];
       coords[numTurns][partNum][1] = aPartCoord[Bunch::y];
       nEntry++; nGoodParticles++;
     } // on Entries within a turn... 
    if (numTurns == 0) {
       numParticles1rstBunch = nGoodParticles;
     }
     numTurns++;
     file.close();
  } // on Files, one file per turn. 
}


std::vector<double> Analysis::get_betatron_tunes(bool isH) {
    std::vector<double> tunes(numParticles1rstBunch, 0.);
    if (isH && (!gotHorizontalTunes)) this->compute_betatron_tunes(isH);
    if (!isH && (!gotVerticalTunes))  this->compute_betatron_tunes(isH);
    for (size_t k=0; k != numParticles1rstBunch; k++) {
      if (isH) tunes[k] = XTunes[k];
      else tunes[k] = YTunes[k];
    }
    return tunes;
}
    
std::vector<double> Analysis::get_XCoords_forTunes(size_t selectedParticle) const {
    std::vector<double> cVals(coords.shape()[0], FLT_MAX);
    if (selectedParticle >= coords.shape()[1]) return cVals;
    for(size_t i=0; i!= coords.shape()[0]; i++) cVals[i] = coords[i][selectedParticle][0];
    return cVals;
}

std::vector<double> Analysis::get_YCoords_forTunes(size_t selectedParticle) const {
    std::vector<double> cVals(coords.shape()[0], FLT_MAX);
    if (selectedParticle >= coords.shape()[1]) return cVals;
    for(size_t i=0; i!= coords.shape()[0]; i++) cVals[i] = coords[i][selectedParticle][1];
    return cVals;
}
//
// Get the transver action, turn after turn, for a selected particles
//
// need to loop again over all files. 
// 
std::vector<double> Analysis::get_transverse_action_for_particle(bool isH,  
                                   size_t partNumber, double alpha, double beta) const {
  std::vector<double> cVals;
  size_t nTurns=0;
  hsize_t      offset[2]={0,0};   // hyperslab offset in the file
  hsize_t      dims_out[2]={0,0};   // the whole thing.
  double gamma = (1. + alpha*alpha)/beta;
  while (nTurns < maxNumberOfTurn) { 
     std::ostringstream sstream;
     sstream << tokenname << "_" << std::setw(4) << std::setfill('0') << nTurns << ".h5";
     std::string fullname(sstream.str());
     // Check if the file exists..
     std::ifstream inFileTmp;
     inFileTmp.open(fullname.c_str());
     if (!inFileTmp.is_open()) break; // We do not increment the turn number. 
     inFileTmp.close();
     inFileTmp.clear();
     H5::H5File file=H5::H5File(fullname.c_str(), H5F_ACC_RDONLY);
     H5::DataSet bunch_ds = file.openDataSet("particles");
     H5::FloatType datatypeDble(H5::PredType::NATIVE_DOUBLE);
     H5::DataSpace spaceBunch = bunch_ds.getSpace();
     size_t nTotalEntries = (size_t) dims_out[0];
     hsize_t      countVect[2]={1,7};    // size of the hyperslab in the file
     H5::DataSpace  memspace(2, countVect);
     MArray1d aPartCoord(boost::extents[7]); // coordinate values.
     size_t nEntry=0;
     while (nEntry < nTotalEntries) { // Slow and tedious get every X or y position, 
       offset[0] = nEntry;
       spaceBunch.selectHyperslab( H5S_SELECT_SET, countVect, offset );
       bunch_ds.read(&aPartCoord[0], datatypeDble, memspace, spaceBunch );
       size_t partNum = (size_t) ( (int) aPartCoord[Bunch::id]);
       if (partNumber != partNum) {nEntry++; continue;}
       double action=0;
       if (isH) action = beta*aPartCoord[Bunch::xp]*aPartCoord[Bunch::xp] + 
	                           2.0*alpha*aPartCoord[Bunch::x]*aPartCoord[Bunch::xp] +
				   gamma*aPartCoord[Bunch::x]*aPartCoord[Bunch::x];
       else action = beta*aPartCoord[Bunch::yp]*aPartCoord[Bunch::yp] + 
	                           2.0*alpha*aPartCoord[Bunch::y]*aPartCoord[Bunch::yp] +
				   gamma*aPartCoord[Bunch::y]*aPartCoord[Bunch::y];
       cVals.push_back(action);
       break;
      } // on Entries within a turn...  
    nTurns++;
    file.close();
  } // on Files, one file per turn. 
  return cVals;
} 
//
// Get the transverse action, for a selected turn, all particles 
//
// Got to the right file. 
// 
std::vector<double> Analysis::get_transverse_action_for_bunch(bool isH,  
                                   size_t iTurn, double alpha, double beta) const {
  std::vector<double> cVals;
  std::ostringstream sstream;
  sstream << tokenname << "_" << std::setw(4) << std::setfill('0') << iTurn << ".h5";
  std::string fullname(sstream.str());
   // Check if the file exists..
  std::ifstream inFileTmp;
  inFileTmp.open(fullname.c_str());
  if (!inFileTmp.is_open()) return cVals; // Empty, no such turns 
  inFileTmp.close();
  inFileTmp.clear();
  H5::H5File file=H5::H5File(fullname.c_str(), H5F_ACC_RDONLY);
  H5::DataSet bunch_ds = file.openDataSet("particles");
  H5::FloatType datatypeDble(H5::PredType::NATIVE_DOUBLE);
  H5::DataSpace spaceBunch = bunch_ds.getSpace();
  hsize_t	offset[2]={0,0};   // hyperslab offset in the file
  hsize_t      dims_out[2]={0,0};   // the whole thing.
  size_t nTotalEntries = (size_t) dims_out[0];
  hsize_t	countVect[2]={1,7};    // size of the hyperslab in the file
  H5::DataSpace  memspace(2, countVect);
  size_t nEntry=0;
  double gamma = (1. + alpha*alpha)/beta;
  MArray1d aPartCoord(boost::extents[7]); // coordinate values.
  while (nEntry < nTotalEntries) { // Slow and tedious get every X or y position, 
    offset[0] = nEntry;
    spaceBunch.selectHyperslab( H5S_SELECT_SET, countVect, offset );
    bunch_ds.read(&aPartCoord[0], datatypeDble, memspace, spaceBunch );
    double action=0;
    if (isH) action = beta*aPartCoord[Bunch::xp]*aPartCoord[Bunch::xp] + 
        			 2.0*alpha*aPartCoord[Bunch::x]*aPartCoord[Bunch::xp] +
        			 gamma*aPartCoord[Bunch::x]*aPartCoord[Bunch::x];
    else action = beta*aPartCoord[Bunch::yp]*aPartCoord[Bunch::yp] + 
        			 2.0*alpha*aPartCoord[Bunch::y]*aPartCoord[Bunch::yp] +
        			 gamma*aPartCoord[Bunch::y]*aPartCoord[Bunch::y];
    cVals.push_back(action);
    nEntry++;
  } // on Entries within a turn...  
  return cVals;
} 
 
double Analysis::get_betatron_averaged_tune(bool isH) {

    if (isH && (!gotHorizontalTunes)) this->compute_betatron_tunes(isH);
    if (!isH && (!gotVerticalTunes)) this->compute_betatron_tunes(isH);
    double tune=0;
    if (isH) {
      size_t nT = numXTunesFound;
      for (size_t i=0; i != nT; i++)  tune += XTunes[i]; tune /= nT;
    } else {
      size_t nT = numYTunesFound;
      for (size_t i=0; i != nT; i++)  tune += YTunes[i]; tune /= nT;
    }
    return tune;
}

double Analysis::get_betatron_tune_spread(bool isH) {

    if (isH && (!gotHorizontalTunes)) this->compute_betatron_tunes(isH);
    if (!isH && (!gotVerticalTunes)) this->compute_betatron_tunes(isH);
    double tune=0;
    double tune2=0;
    if (isH) {
      size_t nT = numXTunesFound;
      if (nT < 2) return 0.;
      for (size_t i=0; i != nT; i++)  {
         tune += XTunes[i];  tune2 += XTunes[i]*XTunes[i];
      } 
      tune/=nT;
      return std::sqrt((tune2 - nT*tune*tune)/(nT - 1));
    } else {
      size_t nT = numYTunesFound;
      for (size_t i=0; i != nT; i++)  {
         tune += YTunes[i];  tune2 += YTunes[i]*YTunes[i];
      } 
      tune/=nT;
      return std::sqrt((tune2 - nT*tune*tune)/(nT - 1));
    }
    return 0.; // Should not happen 
}

void Analysis::compute_betatron_tunes(bool isH) {

    const int memory_fudge_factor = 2;
    fftw_complex *out=(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*memory_fudge_factor*numParticles1rstBunch);
    if (isH) numXTunesFound=0;
    else numYTunesFound=0;
    size_t numTunes=0;
    for (size_t partNum=0; partNum != numParticles1rstBunch; partNum++) {
      size_t lastTurn=0;
      for(size_t iTurn=0; iTurn!=coords.shape()[0]; iTurn++) {
        if ((isH) && (coords[iTurn][partNum][0] == FLT_MAX)) break;
        if ((!isH) && (coords[iTurn][partNum][1] == FLT_MAX)) break;
	lastTurn=iTurn;
      }
      lastTurn++;
      if (lastTurn < minimumRequiredTurnNum) continue; // Lost too early, not entered in Tune list. 
      std::vector<double> x(lastTurn, 0.);
      for(size_t iTurn=0; iTurn!=lastTurn; iTurn++) {
        if (isH) x[iTurn] = coords[iTurn][partNum][0];
        else x[iTurn] = coords[iTurn][partNum][1];
      }
      fftw_plan plan = fftw_plan_dft_r2c_1d(lastTurn+1, &x[0], out, FFTW_ESTIMATE);
      fftw_execute(plan);
      fftw_destroy_plan(plan);
      double tune=0.;
      double ampMax = -1.0;
      fftw_complex* outPtr = out;
      for (size_t i=0; i != lastTurn/2; i++, outPtr++) { // Only through 1/2 of the array 
        double  aCReal=*outPtr[0];
        double aCImag=*outPtr[1];
        std::complex<double> aC(aCReal, aCImag);
        double a = std::norm(aC);
        if (a > ampMax) {
          ampMax = a;
	  tune = ((double) i + 1.0)/lastTurn;
        }
      }
      if (isH) XTunes[numTunes]= tune; 
      else YTunes[numTunes]= tune;
      numTunes++; 
    }  
    if (isH) numXTunesFound=numTunes;
    else numYTunesFound=numTunes;      
    if (isH) gotHorizontalTunes=true;
    else gotVerticalTunes=true;
 }   

