#!/usr/bin/env bwpython
# Electron flock Analysis,  
import numpy
import sys
import locale
import minuit

nBin=50
hist=numpy.zeros(nBin)
histCoord=numpy.zeros(nBin)
histFit=numpy.zeros(nBin)
histErr=numpy.zeros(nBin)
rMax=25.

class electronFlockAnalysis:
#
    def __init__(self):
    
      self.data = [] # intented to be a list of numpy arrays, 8 doubles.
      
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
        
    def histoAndFitR(self):
      global nBin
      global hist
      global histCoord
      global histFit
      global histErr
      global rMax
      for kr in range(nBin):
        hist[kr]=0.
	histErr[kr]=0.
      k=0
      # fill the bins.. 
      for ee in self.data:
        r=numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
	kr=int(nBin*r/rMax)
	if (kr > (nBin-1)): continue
	hist[kr]+=1.
      # we want to fit 1/rdN/dr distribution. 
      for kr in range(nBin):
        rb = rMax*(float(kr) + 0.5)/float(nBin)
	histCoord[kr]=rb
	histErr[kr]=numpy.sqrt(hist[kr])/rb
	if (hist[kr] < 3):
	  histErr[kr]=1.7/rb  
	hist[kr]/=rb
	
      # Now we fit
#      import fitRadial
# test 
# 
#      afRad=fitRadial()
#      print " testing fitRadial ", afRad.fitFuncForFittingR(4.55, 1.0, 0.)

#      myFunc=afRad.fitFuncForFittingR
      m = minuit.Minuit(fitFuncForFittingR, cst=40., slope=1.0e-9, quad=0.0047, cubic=0.0001)
      m.fixed["slope"]=True
#      m = minuit.Minuit(fitFuncForFittingR, cst=55., slope=0.01, quad=0.0025)
      m.maxcalls=5000
      m.migrad()
      aCst=m.values["cst"]
      aSlope=m.values["slope"]
      aQuad=m.values["quad"]
      aCubic=m.values["cubic"]
#      aCubic=0. 
      for kr in range(nBin):
         rb = rMax*(float(kr) + 0.5)/float(nBin)
         histFit[kr]=aCst*(1.0 + aSlope*rb + aQuad*rb*rb + aCubic*rb*rb*rb)
      print " Number of calls ", m.ncalls, " chisq ", m.fval 
      print " cst = ", aCst
      print " slope = ", aSlope
      print " quad = ", aQuad
      print " cubic = ", aCubic
	
      # definition of the chis-square radial fit. 
      
    def pyPlotRHist(self):
      import pylab
      global nBin
      global hist
      pylab.plot(histCoord, hist, 'bs')
      histPlus=hist+histErr
      pylab.plot(histCoord, histPlus, 'r.')
      histMinus=hist-histErr
      pylab.plot(histCoord, histMinus, 'r.')
      pylab.plot(histCoord, histFit, 'g-')
      pylab.show()

class fitRadial:

  def __init__(self):
      nBinDum=100
      
def fitFuncForFittingR(self, cst, slope, quad, cubic):
#def fitFuncForFittingR(cst, slope, quad):
   
    global nBin
    global hist
    global histErr
    global rMax
    chisq=0.
    for kr in range(nBin):
      rb= rMax*(float(kr) + 0.5)/float(nBin)   
      pred=cst*(1.0 + slope*rb + quad*rb*rb + cubic*rb*rb*rb)
#      pred=cst*(1.0 + slope*rb + quad*rb*rb)
      dd = pred-hist[kr]
      chisq+=dd*dd/(histErr[kr]*histErr[kr])
    return chisq
    
       
    
 
