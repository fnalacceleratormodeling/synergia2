#!/usr/bin/env bwpython
# Electron flock Analysis,  
import numpy
import sys
import locale
from scipy.optimize import fmin
import synergia
import s2_fish
from s2_solver_fftw import *
from s2_containers import *
from s2_deposit import *
import plotPotential

class electronFlockAnalysis:
#
    def __init__(self):
    
      self.data = [] # intented to be a list of numpy arrays, 8 doubles.
      self.nBin=50
      self.hist=numpy.zeros(self.nBin)
      self.histCoord=numpy.zeros(self.nBin)
      self.histFit=numpy.zeros(self.nBin)
      self.histErr=numpy.zeros(self.nBin)
      self.rMax=25.
      
    # with no unit chage... 
    def readFlock(self, fName):
      inputF=open(fName,'r')
      text=inputF.readlines()
      print " Number of lines ", len(text)
      n=0
      for line in text:
        if (n==0): # first line is table names.
	  n+=1 
	  continue
        else:
	  ee=numpy.zeros(9)
	  strs=line.split(' ')
#	  print " line ", line 
#	  print " splitted ", strs
	  k=0
          for num in strs:
	    if (len(num) == 0): continue;
	    if (num=='\n'): break
	    if (k == 0):
	       ee[8]=float(locale.atoi(num))
	    else:
	       ee[k-1]=locale.atof(num)
	    k+=1
	  self.data.append(ee)  
	n+=1
#	if (n > 10): break
      inputF.close()
#      print " e, fith one ", self.data[5]
      
    def average(self, k):
      if (len(self.data) == 0): return 0 
      a=0
      for ee in self.data:
        a+=ee[k]
      a/=len(self.data)
      return a
       
    def minmax(self, k):
      if (len(self.data) == 0): return 0 
      aMin=-1.e10
      aMax=1.0e10
      for ee in self.data:
        a=ee[k]
	if (a < aMin): aMin=a
	if (a > aMax): aMax=a
      return aMin, aMax 
    
    
    def studyPotential(self,sizes):
      gridnum=64
      griddim = (gridnum,gridnum,2*gridnum)
      rhoElectrons = Real_scalar_field(griddim,sizes,(0.0,0.0,0.0))
      indices=numpy.array([0,0,0],'i')
      offsets=numpy.array([0., 0., 0.], 'd')
      total_charge_per_cell_vol = 0.0
      h=(rhoElectrons.get_cell_size())
      weight0 = 1.0 / (h[0] * h[1] * h[2])
      rhoElectrons.get_points().zero_all();
      for ee in self.data:
      # Note the unit change in the electron dump file.
        location= numpy.array([ee[0]/1000., ee[2]/1000., ee[4]/1000.],'d')
        indices = rhoElectrons.get_leftmost_indices(location);
        offsets = rhoElectrons.get_leftmost_offsets(location);
        for i in range(2): 
          for j in range(2): 
            for k in range(2):
              weight = weight0 * (1 - i - (1 - 2 * i) * offsets[0]) * \
                                    (1 - j - (1 - 2 * j) * offsets[1]) * \
                                    (1 - k - (1 - 2 * k) * offsets[2])
	      newInd = numpy.array([indices[0] + i, indices[1] + j, indices[2] + k], 'i')		    
              rhoElectrons.get_points().add_to_point(newInd, weight);
              total_charge_per_cell_vol += weight;
        
      
      fftwh = Fftw_helper(griddim, False)
      phiElectron = solver_fftw_open(rhoElectrons,fftwh,0, True)
      plotPotential.plotPotentialVertE(phiElectron, 1.0, "ElectronPotential")
    
    def makeXYZVtkFile(self, sizes, filename, zoff):
    
      nx = 100
      ny = 100
      nz = 50 
      dx=sizes[0]/nx
      dy=sizes[1]/ny
      dz=sizes[2]/nz
      f = open(filename,"w")
      f.write("# vtk DataFile Version 2.0\n")
      f.write("File written by octave script write_vtk\n")
      f.write("ASCII\n\n")
      f.write("DATASET STRUCTURED_POINTS\n")
      line="DIMENSIONS %d %d %d \n" % (nx, ny, nz)
      f.write(line)
      line="ORIGIN    %g %g %g \n" % (-sizes[0]/2, -sizes[1]/2, -sizes[2]/2)
      f.write(line)
      line="SPACING   %g %g %g\n" % (dx, dy, dz)
      f.write(line)
      line="POINT_DATA %d\n" % (nx*ny*nz)
      f.write(line)
      f.write("SCALARS scalars unsigned_short \n")
      f.write("LOOKUP_TABLE default\n")
      vals=numpy.zeros((nx, ny, nz),'i') 
      for ee in self.data:
        iz=int((ee[4]/1000. + 0.5*dz + sizes[2]/2 + zoff)/dz)
#        print " ee[4] ", ee[4], " iz " , iz 
        if (iz < 0): continue
        if (iz > nz-1): continue
        iy=int((ee[2]/4. + 0.5*dy + sizes[1]/2)/dy) # 4 is a rescaling factor.. 
#        print " ee[2] ", ee[1], " iy " , iy 
        if (iy < 0): continue
        if (iy > ny-1): continue
        ix=int((ee[0]/4. + 0.5*dx + sizes[0]/2)/dx)
        if (ix < 0): continue
        if (ix > nx-1): continue
	vals[ix][iy][iz] += 1;
      for k in range(nz):
        for j in range(ny):
	  line=""
          for i in range(nx):
	    line += " %d"% vals[i][j][k] 
          line +="\n"
	  f.write(line)
      f.close()
      
    def makeXYByVtkFile(self, sizes, filename, zoff):
    
      nx = 100
      ny = 100
      nz = 100 
      dx=sizes[0]/nx
      dy=sizes[1]/ny
      dz=sizes[2]/nz
      f = open(filename,"w")
      f.write("# vtk DataFile Version 2.0\n")
      f.write("File written by octave script write_vtk\n")
      f.write("ASCII\n\n")
      f.write("DATASET STRUCTURED_POINTS\n")
      line="DIMENSIONS %d %d %d \n" % (nx, ny, nz)
      f.write(line)
      line="ORIGIN    %g %g %g \n" % (-sizes[0]/2, -sizes[1]/2, -sizes[2]/2)
      f.write(line)
      line="SPACING   %g %g %g\n" % (dx, dy, dz)
      f.write(line)
      line="POINT_DATA %d\n" % (nx*ny*nz)
      f.write(line)
      f.write("SCALARS scalars unsigned_short \n")
      f.write("LOOKUP_TABLE default\n")
      vals=numpy.zeros((nx, ny, nz),'i') 
      for ee in self.data:
        iz=int((ee[3] + 0.5*dz + sizes[2]/2 + zoff)/dz)
        print " ee[3] ", ee[3], " iz " , iz 
        if (iz < 0): continue
        if (iz > nz-1): continue
        iy=int((ee[2]/1000. + 0.5*dy + sizes[1]/2)/dy) 
        print " ee[2] ", ee[2], " iy " , iy 
        if (iy < 0): continue
        if (iy > ny-1): continue
        ix=int((ee[0]/1000. + 0.5*dx + sizes[0]/2)/dx)
        print " ee[0] ", ee[0], " ix " , ix 
        if (ix < 0): continue
        if (ix > nx-1): continue
	vals[ix][iy][iz] += 1;
      for k in range(nz):
        for j in range(ny):
	  line=""
          for i in range(nx):
	    line += " %d"% vals[i][j][k] 
          line +="\n"
	  f.write(line)
      f.close()
        
    def histoAndFitR(self):
      for kr in range(self.nBin):
        self.hist[kr]=0.
	self.histErr[kr]=0.
      k=0
      # fill the bins.. 
      for ee in self.data:
        r=numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
	kr=int(self.nBin*r/self.rMax)
	if (kr > (self.nBin-1)): continue
	self.hist[kr]+=1.
      # we want to fit 1/rdN/dr distribution. 
      for kr in range(self.nBin):
        rb = self.rMax*(float(kr) + 0.5)/float(self.nBin)
	self.histCoord[kr]=rb
	self.histErr[kr]=numpy.sqrt(self.hist[kr])/rb
	if (self.hist[kr] < 3):
	  self.histErr[kr]=1.7/rb  
	self.hist[kr]/=rb
	
      # Now we fit

      myFunc=self.fitFuncForFittingR
      xStart=[20., 0.05, 0.14, 0.001]
      resFit=fmin(myFunc, xStart, xtol = 0.001, ftol = 0.01, maxiter = 50000)
      aCst, aSlope, aQuad, aCubic = resFit
#      aCubic=0. 
      for kr in range(self.nBin):
         rb = self.rMax*(float(kr) + 0.5)/float(self.nBin)
	 rb-=self.rMax/2.0  
         self.histFit[kr]=aCst*(1.0 + aSlope*rb + aQuad*rb*rb + aCubic*rb*rb*rb)
      print " Number of calls chisq ", self.fitFuncForFittingR(resFit) 
      print " cst = ", aCst
      print " slope = ", aSlope
      print " quad = ", aQuad
      print " cubic = ", aCubic
	      
    def histoAndFitPhi(self, useTop):
      for kr in range(self.nBin):
        self.hist[kr]=0.
	self.histErr[kr]=0.
      k=0
      # fill the bins.. 
      for ee in self.data:
        phi=numpy.arctan2(ee[2], ee[0])
	if ((phi < 0.) & (useTop)): continue;
	if ((phi > 0.) & (not useTop)): continue;
	if (phi < 0.):
	  phi*=-1.0 
	kr=int(phi*self.nBin/numpy.pi)
	if (kr > (self.nBin-1)): continue
	self.hist[kr]+=1.
      # we want to fit 1/rdN/dr distribution. 
      for kr in range(self.nBin):
        rb = numpy.pi*(float(kr) + 0.5)/float(self.nBin)
	self.histCoord[kr]=rb
	self.histErr[kr]=numpy.sqrt(self.hist[kr])
	if (self.hist[kr] < 3):
	  self.histErr[kr]=1.7  
	
      # Now we fit

      myFunc=self.fitFuncForFittingPhi
      xStart=[20., 0.05, 0.14, 4000., 0.25]
      resFit=fmin(myFunc, xStart, xtol = 0.001, ftol = 0.01, maxiter = 50000)
      aCst, aSlope, aQuad, aPeak, aSigPeak = resFit
      for kr in range(self.nBin):
         rb = numpy.pi*(float(kr) + 0.5)/float(self.nBin)
	 rb-=numpy.pi/2.0  
         self.histFit[kr]=aCst*(1.0 + aSlope*rb + aQuad*rb*rb)
	 self.histFit[kr]+= aPeak*numpy.exp(-1.0*rb*rb/(2.0*aSigPeak*aSigPeak))
      print " Number of calls chisq ", self.fitFuncForFittingPhi(resFit) 
      print " cst = ", aCst
      print " slope = ", aSlope
      print " quad = ", aQuad
      print " Peak value = ", aPeak 
      print " Sigma Peak = ", aSigPeak 
      
    def pyPlotRHist(self):
      import pylab
      pylab.plot(self.histCoord, self.hist, 'bs')
      histPlus=self.hist+self.histErr
      pylab.plot(self.histCoord, histPlus, 'r.')
      histMinus=self.hist-self.histErr
      pylab.plot(self.histCoord, histMinus, 'r.')
      pylab.plot(self.histCoord, self.histFit, 'g-')
      pylab.show()
      
    def pyPlotPhiHist(self):  # same code... That's good.. we'll probably junk one of these. 
      import pylab
      pylab.plot(self.histCoord, self.hist, 'bs')
      histPlus=self.hist+self.histErr
      pylab.plot(self.histCoord, histPlus, 'r.')
      histMinus=self.hist-self.histErr
      pylab.plot(self.histCoord, histMinus, 'r.')
      pylab.plot(self.histCoord, self.histFit, 'g-')
      pylab.show()

    def fitFuncForFittingR(self, x):
   
      chisq=0.
      cst, slope, quad, cubic=x
      for kr in range(self.nBin):
        rb= self.rMax*(float(kr) + 0.5)/float(self.nBin)
	rb-=self.rMax/2.0  
        pred=cst*(1.0 + slope*rb + quad*rb*rb + cubic*rb*rb*rb)
#       pred=cst*(1.0 + quad*rb*rb + cubic*rb*rb*rb)
#       pred=cst*(1.0 + slope*rb + quad*rb*rb)
        dd = pred-self.hist[kr]
        chisq+=dd*dd/(self.histErr[kr]*self.histErr[kr])
      return chisq
    
    def fitFuncForFittingPhi(self, x):
   
      chisq=0.
      cst, slope, quad, peak, sigPk=x
      for kr in range(self.nBin):
        rb= numpy.pi*(float(kr) + 0.5)/float(self.nBin)
	rb-=numpy.pi/2.0  
        pred=cst*(1.0 + slope*rb + quad*rb*rb)
	pred+=peak*numpy.exp(-1.0*rb*rb/(2.0*sigPk*sigPk))
        dd = pred-self.hist[kr]
        chisq+=dd*dd/(self.histErr[kr]*self.histErr[kr])
      return chisq
       
    
 
