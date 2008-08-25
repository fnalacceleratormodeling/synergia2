#!/usr/bin/env bwpython
# Electron flock, Analysis, simple Main  
import numpy
import sys
import minuit
from electronFlockAnalysis import electronFlockAnalysis 

if ( __name__ == '__main__'):

#  fName="ElectronFlockBC_V3e11ppbMediumPNegLW3_Crossing_9.txt"
  fName="./results1/ElectronFlockDR_V3e11ppbMediumPNegLW13_Crossing_21.txt"
  myAnal=electronFlockAnalysis()
  myAnal.readFlock(fName)
  print " average y ", myAnal.average(2)
#  myAnal.histoAndFitR()
#  myAnal.pyPlotRHist()
#  myAnal.histoAndFitPhi(False)
#  myAnal.pyPlotPhiHist()
  maxBeanPipeInX=0.055
  maxBeanPipeInY=0.025 # Radius..
#  sizeNow = (2.0*maxBeanPipeInX+0.0005, 2.0*maxBeanPipeInY+0.0005, 32.)
  sizeNow = (250.*2.0*maxBeanPipeInX+0.0005, 250.0*2.0*maxBeanPipeInY+0.0005, 10.)
#  myAnal.studyPotential(sizeNow)
  sizeNowBy = (2.0*maxBeanPipeInX+0.0005, 2.0*maxBeanPipeInY+0.0005, 0.04)
  myAnal.makeXYByVtkFile(sizeNowBy, "./results1/tVTKXYbyLW13_C21.vtk", 0.)
