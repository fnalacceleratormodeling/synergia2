#!/usr/bin/env python
import numpy
import math
import sys
import tables
from mpi4py import MPI

x = 0
xprime = 1
y = 2
yprime = 3
z = 4
zprime = 5

x_xprime = 0
x_y = 1
xprime_y = 2
x_yprime = 3
xprime_yprime = 4
y_yprime = 5
x_z = 6
xprime_z = 7
y_z = 8
yprime_z = 9
x_zprime = 10
xprime_zprime = 11
y_zprime = 12
yprime_zprime = 13
z_zprime = 14

xprime_x = 0
y_x = 1
y_xprime = 2
yprime_x = 3
yprime_xprime = 4
yprime_y = 5
z_x = 6
z_xprime = 7
z_y = 8
z_yprime = 9
zprime_x = 10
zprime_xprime = 11
zprime_y = 12
zprime_yprime = 13
zprime_z = 14

def get_diagnostics_old(bunch):
#    means = numpy.zeros([6],numpy.float64)
    means = numpy.zeros([6],'d')
    stds = numpy.zeros([6],numpy.float64)    
    for i in range(0,6):
        means[i] = numpy.average(bunch[i,:])
    mom2s = numpy.zeros([6,6],numpy.float64)
    corrs = numpy.zeros([6,6],numpy.float64)
    for i in range(0,6):
        for j in range(i,6):
            tmp = bunch[i,:]*bunch[j,:]
            mom2s[i,j] = numpy.average(tmp) - means[i]*means[j]
            mom2s[j,i] = mom2s[i,j]
    for i in range(0,6):
        stds[i] = math.sqrt(mom2s[i,i])
    for i in range(0,6):
        for j in range(i,6):
            corrs[i,j] = mom2s[i,j]/(stds[i]*stds[j])            
	    corrs[j,i] = corrs[i,j]
    return means,stds,mom2s,corrs


def get_diagnostics(bunch):
    MPI.COMM_WORLD.Barrier()
#    num_local=bunch.num_particles_local() it doesn't work
    num_local=numpy.shape(bunch)[1]
    Nptlocal=numpy.array(num_local,'d')
    data_local = numpy.zeros([6],'d')
    tmp=numpy.zeros([6,6], 'd')
      
    
    for i in range(0,6):
        data_local[i] = numpy.sum(bunch[i,:])
        
   
    for i in range(0,6):
        for j in range(i,6):
            tmp[i,j] = numpy.sum(bunch[i,:]*bunch[j,:])
	    tmp[j,i]=tmp[i,j]
               	
 
    means=numpy.zeros([6],'d')
    num_part=numpy.array(0.0, 'd') 
    mom2s = numpy.zeros([6,6], 'd')
    MPI.COMM_WORLD.Allreduce([Nptlocal,MPI.DOUBLE], [num_part, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce([data_local,MPI.DOUBLE], [means, MPI.DOUBLE], op=MPI.SUM)
    MPI.COMM_WORLD.Allreduce([tmp,MPI.DOUBLE], [mom2s, MPI.DOUBLE], op=MPI.SUM)
    
    means=means/num_part
    mom2s= mom2s/num_part-numpy.outer(means,means)	
    
    stds = numpy.zeros([6],numpy.float64) 
    corrs = numpy.zeros([6,6],numpy.float64)
    for i in range(0,6):
        stds[i] = math.sqrt(mom2s[i,i])
   
    for i in range(0,6):
        for j in range(i,6):
            corrs[i,j] = mom2s[i,j]/(stds[i]*stds[j])            
	    corrs[j,i] = corrs[i,j]
    
    return means,stds,mom2s,corrs


class Diagnostics_impact:
    def __init__(self,units):
        self.s = []
        self.means = []
        self.stds = []
        self.u = units
        # n.b. we are ignoring mom2, corrs!

    def add(self,s,bunch):	
        means,stds,mom2s,corrs = get_diagnostics(bunch.particles())
	for i in range(0,6):
	    means[i]=means[i]/self.u[i]
	    stds[i] = stds[i]/self.u[i]
        self.s.append(s)
	self.means.append(means)
	self.stds.append(stds)
       

    def get_coord_stds(self,bunch):
        print "about to get diagnostics"
        print type(bunch.particles())
        means,stds,mom2s,corrs = get_diagnostics(bunch.particles())
        print "got diagnostics"
        return (stds[x],stds[y],stds[z])
    
    def get_s(self):
        return numpy.array(self.s)    
    def get_means(self):
        return numpy.array(self.means)
    def get_stds(self):
        return numpy.array(self.stds)
    
    def write(self,filename_prefix):
        # same format as fort.24
	if MPI.COMM_WORLD.Get_rank() ==0:
	    fx = open(filename_prefix+"_x.dat","w")
	    fy = open(filename_prefix+"_y.dat","w")
	    fz = open(filename_prefix+"_z.dat","w")
	    for i in range(0,len(self.s)):
	        #fx.write("%g %g %g %g %g\n" % \
		    #(self.s[i],
		    #self.mean[xprime][i], self.std[xprime][i],
		    #self.mean[x][i], self.std[x][i]))
	        #fy.write("%g %g %g %g %g\n" % \
		    #(self.s[i],
		    #self.mean[yprime][i], self.std[yprime][i],
		    #self.mean[y][i], self.std[y][i])) 
	        #fz.write("%g %g %g %g %g\n" % \
		    #(self.s[i],
		    #self.mean[zprime][i], self.std[zprime][i],
		    #self.mean[z][i], self.std[z][i]))       
		    
		    
	        fx.write("%g %g %g %g %g\n" % \
                     (self.s[i],
                     self.means[i][x], self.stds[i][x],
                     self.means[i][xprime], self.stds[i][xprime]))
                fy.write("%g %g %g %g %g\n" % \
                     (self.s[i],
                     self.means[i][y], self.stds[i][y],
                     self.means[i][yprime], self.stds[i][yprime]))
                fz.write("%g %g %g %g %g\n" % \
                     (self.s[i],
                     self.means[i][z], self.stds[i][z],
                     self.means[i][zprime], self.stds[i][zprime])) 		    
	    fx.close()
	    fy.close()
	    fz.close()

    def write_hdf5(self,filename_prefix,compress_level=1):
        f = tables.openFile(filename_prefix+".h5",mode = "w")
        # n.b. filter (and compress_level) not (yet) used
        filter = tables.Filters(complevel=compress_level)
        root = f.root
        hdfarray = f.createArray(root,'s',numpy.array(self.s),"position")
        hdfarray = f.createArray(root,'mean',numpy.array(self.means),"centroid")
#        hdfarray = f.createArray(root,'mom2',numpy.array(self.mom2s),"second moments")
#        hdfarray = f.createArray(root,'corr',numpy.array(self.corrs),"correlation coefficients")
#        hdfarray = f.createArray(root,'diagmom4',numpy.array(self.diagmom4s),"fourth moments on diagonal")
        hdfarray = f.createArray(root,'std',numpy.array(self.stds),"standard deviation")
        #hdfarray = f.createArray(root,'emitx',numpy.array(self.emitxs),"x emittance")
        #hdfarray = f.createArray(root,'emity',numpy.array(self.emitys),"y emittance")
        #hdfarray = f.createArray(root,'emitz',numpy.array(self.emitzs),"z emittance")
        #hdfarray = f.createArray(root,'emitxy',numpy.array(self.emitxys),"x-y emittance")
        #hdfarray = f.createArray(root,'emitxyz',numpy.array(self.emitxyzs),"x-y-z emittance")
        hdfarray = f.createArray(root,'units',numpy.array(self.u),"units")
        f.close()