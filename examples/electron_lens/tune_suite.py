#!/usr/bin/env python
import numpy as np

# suite of functions that can do different tune calculations

# calculate the simple x y z tunes just from fft
def basic_tunes(coords):
    # coords has shape (6,n)
    n = coords.shape[1]
    df = 1.0/n
    # get maximum amplitude location
    maxn = n/2
    xt = np.abs(np.fft.fft(coords[0,:]))
    # probably don't want 0 frequency
    locmax = np.argmax(xt[0:maxn])
    xtune = (locmax)*df
    yt = np.abs(np.fft.fft(coords[2,:]))
    locmax = np.argmax(yt[0:maxn])
    ytune = (locmax)*df
    zt = np.abs(np.fft.fft(coords[4,:]))
    locmax = np.argmax(zt[0:maxn])
    ztune = (locmax)*df
    return (xtune,ytune,ztune)

# calculate tunes using the CFT algorithm
def cft_tunes(coords, search_range=((0.0,0.5),(0.0,0.5),(0.0,0.5))):
    # coords has shape (6,n)
    n = coords.shape[1]
    #if n>100:
    #    print "warning! you have a lot of points.  This will be slow"

    # normal FFT precision is 1/n, CFT gives addition 1/n factor
    df = 1.0/n**2
    t = np.arange(n,dtype='d')

    # loop over x,y,z
    tunes = np.zeros(3)
    for pln in range(3):
        f = np.arange(search_range[pln][0],search_range[pln][1],df)
        nf = len(f)
        cft = np.zeros(nf,dtype='d')
        idxcoord = 2*pln

        for i in range(nf):
            expfact = np.exp(-2.0*np.pi*(1.0j)*f[i]*t)
            cft[i] = abs(np.dot(expfact, coords[idxcoord,:]))

        cftmaxloc = np.argmax(cft)
        tunes[pln] = f[cftmaxloc]
    
    return tuple(tunes)

# get interpolated tunes
def interp_tunes(coords, search_range=((0.0,0.5),(0.0,0.5),(0.0,0.5))):
    # coords has shape (6,n)
    n = coords.shape[1]
    # loop pver x. y, z
    tunes = np.zeros(3)
    maxn = n/2
    df = 1.0/n
    f = np.arange(n,dtype='d')/n
    for pln in range(0,6,2):
        xt = abs(np.fft.fft(coords[pln, :]))
        locmax = np.argmax(xt[0:maxn])
        xtp = xt[locmax]
        if xt[locmax-1] > xt[locmax+1]:
            dir=-1.0
            xtp2 = xt[locmax-1]
        else:
            dir=1.0
            xtp2 = xt[locmax+1]

        tune = f[locmax] + dir*(xtp2/(xtp+xtp2))/n
        tunes[pln/2] = tune

    return tuple(tunes)
        
    


# get the (fractional) tunes of a set of coordinates from a single track
def refined_tunes(coords):
    # coords has shape (6,n)
    n = coords.shape[1]
    df = 1.0/n

    if n <= 100:
        # really, if there are less than 100 points, just use
        # the cft tunes
        ctunes = cft_tunes(coords)
        return ctunes

    # first get the basic tunes
    btunes = basic_tunes(coords)

    if np.any(btunes == 0.0):
        # if some of the basic tunes are 0.0, then there was
        # some problem and it's probably not right. Get the CFT tunes
        # for the first 100 points.
        ctunes = cft_tunes(coords[:,:100])
        xrange = (ctunes[0]-2.0*df, ctunes[0]+2.0*df)
        yrange = (ctunes[1]-2.0*df, ctunes[1]+2.0*df)
        zrange = (ctunes[2]-2.0*df, ctunes[2]+2.0*df)
    else:
        xrange = (btunes[0]-2.0*df, btunes[0]+2.0*df)
        yrange = (btunes[1]-2.0*df, btunes[1]+2.0*df)
        zrange = (btunes[2]-2.0*df, btunes[2]+2.0*df)
        
    tunerange = (xrange,yrange,zrange)
    ctunes = cft_tunes(coords, tunerange)
    return ctunes
