#!/usr/bin/env python

import hist2d
import Numeric

have_pylab = False
try:
    import pylab
    have_pylab = True
except:
    have_pylab = False
    
import MLab

def density_plot(px_in,py_in,n):
    global have_pylab
    if have_pylab:
        px = Numeric.array(px_in)
        py = Numeric.array(py_in)
        xmin = MLab.min(px)
        xmax = MLab.max(px)
        ymin = MLab.min(py)
        ymax = MLab.max(py)
        X,Y = pylab.meshgrid(Numeric.arange(xmin,xmax,(xmax-xmin)/n),
                             Numeric.arange(ymin,ymax,(ymax-ymin)/n))
        freq = Numeric.zeros([n,n],'d')
        hist2d.hist2d(px,xmin,xmax,n,py,ymin,ymax,n,len(px),freq)
        pylab.pcolor(X,Y,freq,shading='flat')
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
