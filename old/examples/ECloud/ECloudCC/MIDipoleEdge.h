#ifndef HAVE_ECLOUD_MIDIPOLEEDGE_H
#define HAVE_ECLOUD_MIDIPOLEEDGE_H
#include <cmath>
#include <vector>

class MIDipoleEdge {
  public: 
    MIDipoleEdge();
    double getFieldYFact(double z);
  
  private:
    const static size_t numPts=100;
    static bool initDone;
    static double zs[numPts];
    static double vals[numPts];
    static double zMax;
};

// Combine quadrupole/Dipole edge.. 

class MIDipQuadEdge {
  public:
    MIDipQuadEdge();
    double getFieldY(const double *location); // dimensioned to 3 
    double getFieldX(const double *location);
    double getFieldYFact(double z);
    static double BDipoleOverMIEnergy;
    static double BQuadStrOverMIEnergy;
    
   private:
    MIDipoleEdge miDipE;
    double MIEnergy;
    static double zOffDipEdge; // With respect to the ziptrack data..
    
   public:
     inline double setMIEnergy(double e) {MIEnergy=e;}
     // Not negotiable... 
//      inline double setZOffDipEdge(double z) {zOffDipEdge=z;}
    
};

#endif
