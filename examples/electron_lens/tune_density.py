 # -*- coding: utf-8 -*-
import os
import sys
import tables
from mpi4py import MPI
import argparse

#from mpl_toolkits.mplot3d import Axes3D as axes3d
#from matplotlib import cm
#from matplotlib  import font_manager
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
### for Palatino and other serif fonts use:
##rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
import numpy
from numpy.fft import fft





if __name__ == "__main__":
    #if len(sys.argv) < 3:
    #    raise RuntimeError,"usage: tunes_density.py tracks-filename  output_label  "
    os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"

    parser = argparse.ArgumentParser()
    parser.add_argument("track_file", type=str, help="input file for tracks")
    parser.add_argument("output_label", type=str, help="label for track file")
    parser.add_argument("--initurn", default=96, help="initial turn for tune calculation", type=int)
    parser.add_argument("--endturn", default=2496, help="final turn for tune calculation", type=int)
    parser.add_argument("--peakonly", action='store_true', help="use only peak for tune spectrum")
    parser.add_argument("--tracklist", default=None, type=str)
    parser.add_argument("--npart", default=10000, help="number of particles")    
    args = parser.parse_args()
    #print args
    initurn=args.initurn
    endturn=args.endturn
    nmacropartices=args.npart
    output_label=args.output_label
    output_label += "_"+str(initurn)
    tracks_file = args.track_file
    peakonly = args.peakonly
    tracklist = args.tracklist
    
    com_size=MPI.COMM_WORLD.Get_size() 
    rank=MPI.COMM_WORLD.Get_rank()
    if (rank==0):
        print "Reading tracks from file ", tracks_file
        print "number of macroparticles considered=", nmacropartices
        print "start turn: ", initurn
        print "end turn: ", endturn
        print "output label: ", output_label
        print "peakonly: ", peakonly
        if tracklist:
            print "use tracks from file: ", tracklist
        else:
            print "using all tracks from file"

    #print "com_size=",com_size
    #print "rank=",rank

    #h5f = tables.openFile(tracks_file)
    h5f = tables.open_file(tracks_file)
    h5tracks = h5f.root.track_coords[:, 6, 0]
    h5tracksq = list(h5tracks)
    if tracklist:
        fp = open(tracklist)
        usetracks = [float(pid) for pid in fp.readlines()]
        fp.close()
        if rank == 0:
            print "read in ", len(usetracks), " particle IDs"
    else:
        usetracks = h5tracks

    ntracks = len(usetracks)

    trkindex = [0]*ntracks
    for i,trk in enumerate(usetracks):
        trkindex[i] =  h5tracksq.index(trk)

    nmacropartices = min(nmacropartices, ntracks)
    if rank == 0:
        print "number of macroparticles: ", nmacropartices

    lnpart=int(numpy.floor(1.0*nmacropartices)/com_size)
    rest=nmacropartices-lnpart*com_size  
    counts=[]
    offsets=[]
    for ip in range(com_size):
      if (ip<rest): 
        counts.append(lnpart+1)
        offsets.append(ip*(lnpart+1))
      else:
        counts.append(lnpart)
        offsets.append(ip*lnpart+rest)


    #print "rank=",rank,"  counts=",counts, "offsets=",offsets
    local_count=counts[rank]
    local_offset=offsets[rank]
    print "rank=",rank,"local  count=",local_count, "local_offset=",local_offset

    #ntracks = h5f.root.track_coords.shape[0]
    #if (ntracks != nmacropartices):
    #      h5f.close()
    #      raise RuntimeError,"ntracks read from file is not eqalt to input nmacropartices"
                 
    nturns=endturn-initurn
    lenspect = nturns/2
    print "length of spectrum: ", lenspect
    local_npart=local_count
      
    #print "local number of particles=",local_npart, "rank=",rank
    #print "number of turns=",nturns, "rank=",rank

    local_tune_density=numpy.zeros((lenspect,lenspect))
    fourierx_spct=numpy.zeros(lenspect)
    fouriery_spct=numpy.zeros(lenspect)

    #print "rank: ", rank, " local_npart: ", local_npart
    for idp in range(local_npart):
        print "rank: ", rank, " processing particle ", usetracks[local_offset+idp]
        #print "idp: ", idp, ", trkindex[idp]: ", trkindex[idp], ", endturn: ", endturn, ", initurn: ", initurn
        xcoords = h5f.root.track_coords[trkindex[idp], 0, initurn:endturn]
        ycoords = h5f.root.track_coords[trkindex[idp], 2, initurn:endturn]
        
        fourierx_spct=numpy.square(abs(fft(xcoords)[:lenspect]))
        if peakonly:
            peakbin = fourierx_spct.argmax()
            fourierx_spct=numpy.zeros(lenspect)
            fourierx_spct[peakbin] = 1.0
            #print "rank: ", rank, ", particle: ", idp, " x spectrum peak bin: ", peakbin
        normx= numpy.sum(fourierx_spct)
        fouriery_spct=numpy.square(abs(fft(ycoords)[:lenspect]))
        if peakonly:
            peakbin = fouriery_spct.argmax()
            fouriery_spct=numpy.zeros(lenspect)
            fouriery_spct[peakbin] = 1.0
            #print "y spectrum peak bin: ", peakbin
        normy= numpy.sum(fouriery_spct)
        if ((trkindex[idp]>=0) and (normx >1.e-16) and (normy >1.e-16)):
            local_tune_density += numpy.outer(fourierx_spct,fouriery_spct)/normx/normy

  
    #print "local total density checkx=",numpy.sum(local_tune_density),"  rank=",rank
    #print "**************************************************"


    tune_density=numpy.zeros((lenspect,lenspect))
    MPI.COMM_WORLD.Barrier()
    MPI.COMM_WORLD.Reduce(local_tune_density,tune_density , op=MPI.SUM, root=0)
    #tune_density /=1.*nmacropartices
  

    if (rank==0):
       print "final total density checkx=",numpy.sum(tune_density),"  rank=",rank
       print "************************"
       numpy.save("DOT_2D_"+output_label+".npy",tune_density)
       print "DOT saved"
   
