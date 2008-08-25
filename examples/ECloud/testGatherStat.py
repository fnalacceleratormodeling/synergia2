
import numpy
import sys
from mpi4py import MPI

class aFlock:
  def __init__(self):
    self.data = [] # stuff..
    self.lenDataAll=0 # For MPI use... 
    self.lenDataBPAll=0 # same
    self.nBad=0  # Number of electron lost in inf. loop or other propagation errors. 

  def addSome(self, numAv):
    nn=int(numpy.random.rand(1)*float(numAv))
    self.data=[]
    for k in range(nn):
      ee=numpy.zeros(5)
      self.data.append(ee)

    print " MPI Rank ", MPI.rank, " created ", len(self.data), " thingies " 

  def gatherNumInVaccum(self):
      if (MPI.size == 1):
        self.lenDataAll=len(self.data)
      self.lenDataAll=0
      if (MPI.rank !=0):
	MPI.COMM_WORLD.Send([len(self.data), 1, MPI.INT],dest=0, tag=None)
      else:
	for proc in range(1,MPI.size):
	  aVTmp=MPI.COMM_WORLD.Recv(source=proc)
          self.lenDataAll+=aVTmp[0]
#          print " MPI Rank ", MPI.rank, " received ", aVTmp[0] 
	self.lenDataAll+=len(self.data)
#      print " MPI Rank ", MPI.rank, " all collected " 
      self.lenDataAll=MPI.COMM_WORLD.Bcast(self.lenDataAll, root=0)
#      print " MPI Rank ", MPI.rank, " Broadcast val ", self.lenDataAll 
      
  def numInVaccum(self):
       return self.lenDataAll
        
aF=aFlock()
aF.addSome(45)
MPI.COMM_WORLD.Barrier()
aF.gatherNumInVaccum()
print " MPI Rank ", MPI.rank, " Gathered, total ", aF.numInVaccum()
aNumHere=0.
if (MPI.rank==0):
   aNumHere=99.
aNumHere=MPI.COMM_WORLD.Bcast(aNumHere, root=0)
print " MPI Rank ", MPI.rank, " aNumHere ", aNumHere

MPI.Finalize()

   
