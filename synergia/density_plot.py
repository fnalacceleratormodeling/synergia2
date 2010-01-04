#!/usr/bin/env python

import hist2d
import numpy

have_pylab = False
try:
    import pylab
    have_pylab = True
except:
    have_pylab = False
    
import numpy

def density_plot(px_in,py_in,n):
    global have_pylab
    if have_pylab:
        px = numpy.array(px_in)
        py = numpy.array(py_in)
        xmin = numpy.min(px)
        xmax = numpy.max(px)
        ymin = numpy.min(py)
        ymax = numpy.max(py)
        X,Y = pylab.meshgrid(numpy.arange(xmin,xmax,(xmax-xmin)/n),
                             numpy.arange(ymin,ymax,(ymax-ymin)/n))
        freq = numpy.zeros([n,n],'d')
        hist2d.hist2d(px,xmin,xmax,n,py,ymin,ymax,n,len(px),freq)
        pylab.pcolor(X,Y,freq,shading='flat')
    else:
        print "density_plot: pylab not available"

def density_plot_alex(px_in,py_in,n,text=None):
    nh0=13.249
    nv0=8.211	
    global have_pylab
    if have_pylab:
        px = numpy.array(px_in)
        py = numpy.array(py_in)
        xmin1 = numpy.min(px)
        xmax1 = numpy.max(px)	
        ymin1 = numpy.min(py)
        ymax1 = numpy.max(py)
	print "xmax=",xmax1, "ymax=",ymax1 
	print "xmin=",xmin1, "ymin=",ymin1
	
	
	
	dx1=(xmax1-xmin1)/n
	dy1=(ymax1-ymin1)/n
	
	xmax = numpy.max((xmax1,nh0+dx1))
	ymax = numpy.max((ymax1,nv0+dy1))	
	
	xmin=13.21#xmin1
	ymin=8.156#ymin1
	
	dx=(xmax-xmin)/n
	dy=(ymax-ymin)/n
	print "dx,dy=",dx,dy
        X,Y = pylab.meshgrid(numpy.arange(xmin,xmax+dx-0.0000001,dx),
                            numpy.arange(ymin,ymax+dy-0.00000001,dy))
	    
	     
        freq = numpy.zeros([n,n],'d')
	
	
        hist2d.hist2d(px,xmin,xmax,n,py,ymin,ymax,n,len(py),freq)
	

        pylab.pcolor(X,Y,freq,shading='flat')
	
	x1 = numpy.array([nh0])
	y1 = numpy.array([nv0])
	pylab.plot(x1,y1,'ro',markersize=10)
	
	
	pylab.xlim( (xmin, xmax) )
	pylab.ylim( (ymin, ymax) )	
	pylab.xlabel(r'$\nu_H$',{'color'    : 'k', 'fontsize'   : 20 })
	pylab.ylabel(r'$\nu_V$',{'color'    : 'k', 'fontsize'   : 20 })
	
	xarray= numpy.arange(xmin,xmax+dx,5*(xmax-xmin)/n)
	xlabel=[]	
	for a in xarray:
	  a1=numpy.floor(a*1000)/1000	
	  xlabel.append(str(a1))		
	pylab.xticks(numpy.arange(xmin,xmax+dx,5*(xmax-xmin)/n), xlabel)
			
       # pylab.text(xmin+1*dx,ymax-3*dy, 'I=2.69\nE=4 \nppb=4.2e+11  ',{'color'    : 'r', 'fontsize'   : 20 })
	pylab.text(xmin+0.002,ymax-0.01, text,{'color'    : 'r', 'fontsize'   : 20 })
	
    
    else:
        print "density_plot: pylab not available"





if __name__ == '__main__':
    import tables
    import sys
    filename = sys.argv[1]
    xcol = int(sys.argv[2])
    ycol = int(sys.argv[3])
    n = int(sys.argv[4])
    
    f = tables.openFile(filename)
    p = f.root.particles
    
    density_plot(p[xcol,:],p[ycol,:],n)
    f.close()
    pylab.show()
