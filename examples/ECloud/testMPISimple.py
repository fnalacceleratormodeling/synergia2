import numpy
from mpi4py import MPI


class testMPISimple:
#
    def __init__(self):
    
      self.data = [] # intended to hold phasespace vectors. 
      self.myRank=MPI.rank
      self.aXVal=0.
      self.aXValErr=0.
      self.aYVal=0.
      self.aYValErr=0.
      self.numP=10
      
    def generate(self):
      if (self.myRank == 2): # lazy processor.
        return 
      x=numpy.random.normal(5., 1., self.numP)
      y=numpy.random.normal(10., 1., self.numP)
      aM=0.
      aSq=0.
      for i in range(self.numP):
        ee=numpy.array([x[i], y[i]], 'd')
	self.data.append(ee)
      
    def analyze(self):  
      aXM=0.
      aXSq=0.
      aYM=0.
      aYSq=0.
      self.aXVal = 0.
      self.aYVal = 0.
      self.aXValErr=-1.0
      self.aYValErr=-1.0
      if (len(self.data) < 2): return
      for ee in self.data:
        aXM += ee[0]
        aYM += ee[1]
        aXSq += ee[0]*ee[0]
        aYSq += ee[1]*ee[1]
	
      self.aXVal = aXM/self.numP
      self.aYVal = aYM/self.numP
      self.aXValErr=numpy.sqrt((aXSq - self.numP*self.aXVal*self.aXVal)/(self.numP-1))
      self.aYValErr=numpy.sqrt((aYSq - self.numP*self.aYVal*self.aYVal)/(self.numP-1))
    
    # Simple send receive.. We assume all the lists are the same size.
    # very inefficient... May wait for Godo is data is length 0 !!  
    def spitOut(self):
 
      if (self.myRank == 0): 
        fName="TesMPIXYOut.txt"
        out=open(fName,'w')
        out.write("x y \n")
      for ee in self.data:
        if (self.myRank != 0):
	  MPI.COMM_WORLD.Send(ee,dest=0, tag=None)
        else:
	  for proc in range(1,MPI.size):
	    eeTmp=MPI.COMM_WORLD.Recv(source=proc)
	    line=" %10.5g " % eeTmp[0]
	    line+=" %10.5g \n" % eeTmp[1] 
	    out.write(line)
# on node 0	    
	  line=" %10.5g " % ee[0]
	  line+=" %10.5g \n" % ee[1] 
	  out.write(line)
#
      if (self.myRank == 0): 
	 out.close()
	 
    def spitOutFast(self):
# with test on empty, and send averages..
 
#      print " From node ", self.myRank, " spit out fast..." 
      sumDataTmp=numpy.array([float(len(self.data)), self.aXVal, self.aXValErr, \
                               self.aYVal, self.aYValErr],'f')
      sumAllTmp=[]
      if (self.myRank == 0):
        for proc in range(0,MPI.size):
          sumAllTmp.append(sumDataTmp) # temporary copy
      
      if (self.myRank != 0):
	MPI.COMM_WORLD.Send(sumDataTmp,dest=0, tag=None)
      else:
	for proc in range(1,MPI.size):
	  sumDataTmp2=MPI.COMM_WORLD.Recv(source=proc)
          sumAllTmp[proc]=sumDataTmp2
	  
#      print " From node ", self.myRank, " at phase 2 " 
      if (self.myRank == 0): 
        fName="TesMPIXYOut.txt"
        out=open(fName,'w')
        out.write("x y \n")
      if (self.myRank != 0):
	if (len(self.data) > 0): 
	  MPI.COMM_WORLD.Send(self.data,dest=0, tag=None)
      else:
	for proc in range(1,MPI.size):
	  if (sumAllTmp[proc][0] > 0.): 
	    dataTmp=MPI.COMM_WORLD.Recv(source=proc)
	    for ee in dataTmp:
	      line=" %10.5g " % ee[0]
	      line+=" %10.5g \n" % ee[1] 
	      out.write(line)
	  else:
	    print " Node ", proc, " has produced no data "     
# on node 0	    
	for ee in dataTmp:
	  line=" %10.5g " % ee[0]
	  line+=" %10.5g \n" % ee[1] 
	  out.write(line)
#
      if (self.myRank == 0): 
	 out.close()
	 # and handle the analysis
	 proc=0  
	 for asum in sumAllTmp:
	   print " for node  ", proc, " Average in X is " , asum[1], " +- ",  asum[2] 
	   print " .......  Average in Y is " , asum[3], " +- ",  asum[4] 
	   proc+=1
# Obsolete.. 	    
    def sumOut(self):
      aXVSum=0.
      if (self.myRank !=0):
	MPI.COMM_WORLD.Send([self.aXVal, 1, MPI.DOUBLE],dest=0, tag=None)
      else:
	for proc in range(1,MPI.size):
	  aVTmp=MPI.COMM_WORLD.Recv(source=proc)
	  print " Receive aXVal value from node ", proc, " = ", aVTmp
	  aXVSum+=aVTmp[0]
	aXVSum+=self.aXVal
        aXVSum/=MPI.size
	print " Collector node, final value ", aXVSum 
	  
