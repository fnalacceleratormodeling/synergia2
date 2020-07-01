from __future__ import print_function
# -*- coding: utf-8 -*-
import os
import sys
import tables
from matplotlib  import pyplot
from mpl_toolkits.mplot3d import Axes3D as axes3d
from matplotlib import cm
from matplotlib  import font_manager
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import numpy

booster=0
modes=0
oforodo=0
iota=0
ell=1


if booster:
  import booster_constants as beamconst
elif modes:
  import modes_constants as beamconst
elif oforodo:
 import oforodo_constants as beamconst
elif iota: 
   import iota_constants as beamconst
elif ell:    
    import electron_lens_constants as beamconst

    

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("usage: pace_charge_analyze.py file1 file2 ...")
    
    os.environ["HDF5_DISABLE_VERSION_CHECK"] = "2"

    BE_tune_shifts=0
    r0=1.535e-18 # proton clasical radius = e^2/(4 pi epsilon m_p c^2)
    Npb=beamconst.Npb
    Q_s=beamconst.Q_s
    orbit_length=beamconst.orbit_length
    gamma=beamconst.gamma
    beta=beamconst.beta
    w0=beamconst.w0 #2pi/L
   # zrms=beamconst.zrms
    #factor=Npb*r0/numpy.sqrt(2*numpy.pi)/(gamma*gamma*gamma*beta*beta*w0)
    #initial_num_macroparticles=beamconst.initial_num_macroparticles
    #factor=Npb*r0/numpy.sqrt(2*numpy.pi)/(gamma*gamma*gamma*beta*beta*w0)/initial_num_macroparticles

    initial_num_macroparticles=beamconst.initial_num_macroparticles
    factor=Npb*r0/numpy.sqrt(2*numpy.pi)/(gamma*gamma*gamma*beta*beta*w0)/initial_num_macroparticles


    yrange_low=beamconst.yrange_low
    yrange_up=beamconst.yrange_up
    xrange_low=beamconst.xrange_low
    xrange_up=beamconst.xrange_up
    Nplim_low=beamconst.Nplim_low
    Nplim_up=beamconst.Nplim_up
 
    print("Npb=",Npb)
# All bunch quantities read are in the bunch z lab  frame
    for filen in range(1,len(sys.argv)):
        spfile = sys.argv[filen]
        print("Reading file ", spfile)
        h5f = tables.open_file(spfile)
        if hasattr(h5f.root,'inc_tune_shift'):
              print("space_charge rectangular or 3d hockney")
              inc_tune_shift=h5f.root.inc_tune_shift[:,:,:]
              mean=h5f.root.mean[:,:]
              step_betas=h5f.root.step_betas[:,:]        
              s=h5f.root.s[:]
              repetition=h5f.root.repetition[:]
              std=h5f.root.std[:,:]
              macroparticles=h5f.root.num_macroparticles[:]
        h5f.close()
        print("File ", spfile," read")
        print("len(repetition)=",len(repetition))
         


        nkicks=0;
        for i in range(len(repetition)):    
          if (repetition[i]==0):
            nkicks += 1
          else:
            break

        
        print("number of spc kicks per turn = ", nkicks)
        nturns=int(len(s)/nkicks)
        print("number of turns=",nturns)
        print(" s len=",len(s))
        xtune_shift=[]
        ytune_shift=[]
        be_xtune_shift=[]
        be_ytune_shift=[]
        nregression=2
        for nt in range(nturns):
            xtune_shift.append(sum(inc_tune_shift[nt*nkicks:(nt+1)*nkicks,0,nregression]))
            ytune_shift.append(sum(inc_tune_shift[nt*nkicks:(nt+1)*nkicks,1,nregression]))
            xinteg=0. # integrate beta(s)/zrms*xrms*(xrms+yrms)
            yinteg=0.# integrate beta(s)/zrms*yrms*(xrms+yrms)
            if (BE_tune_shifts):
              for step in range(nt*nkicks, (nt+1)*nkicks):
                  tmp=(std[step,0]+std[step,1])*std[step,2]*beta
                  tmp /= macroparticles[step];
                  xinteg += step_betas[step,0]/(std[step,0]*tmp)
                  yinteg += step_betas[step,0]/(std[step,1]*tmp)
              be_xtune_shift.append(factor*xinteg/nkicks/Q_s)
              be_ytune_shift.append(factor*yinteg/nkicks/Q_s)


        labelname="file "+str(filen)
        be_labelname="BE file "+str(filen)

        pyplot.figure(1)
       
        pyplot.subplot(2,1,1)
        #pyplot.plot(xtune_shift,"o-",label=labelname)
        pyplot.plot(xtune_shift,"o-",label="horizontal")
        numpy.savetxt("xtune_shift.txt", xtune_shift, header="x tune shift from run directory "+os.getcwd())
        if (BE_tune_shifts):
           pyplot.plot(be_xtune_shift,"x-",label=be_labelname)
       # pyplot.xlim([xrange_low,xrange_up])
       # pyplot.ylim([yrange_low,yrange_up])
        #pyplot.ylim([0,1])
        #pyplot.xlabel('turn',fontsize=30)
        pyplot.xlabel('turn')
        #pyplot.ylabel('Hor $\Delta Q_{sc}/Q_{s}$',fontsize=30)
        #pyplot.ylabel(r'Hor $\Delta Q_{sc}/Q_{s}$',fontsize=30)
        #pyplot.ylabel(r'Hor $\Delta Q_{sc}$',fontsize=30)
        pyplot.ylabel(r'Hor $\Delta Q_{sc}$')
        #pyplot.xticks(fontsize=25)
        #pyplot.yticks(fontsize=25)
        pyplot.legend(loc=0)


        pyplot.subplot(2,1,2)
        #pyplot.plot(ytune_shift,"o-",label=labelname)
        pyplot.plot(ytune_shift,"o-",label="vertical")
        numpy.savetxt("ytune_shift.txt", ytune_shift, header="y tune shift from run directory "+os.getcwd())
        if (BE_tune_shifts):
           pyplot.plot(be_ytune_shift,"x-",label=be_labelname)
      #  pyplot.xlim([xrange_low,xrange_up])
      #  pyplot.ylim([yrange_low,yrange_up])
        #pyplot.ylim([0,1])
      
        #pyplot.xlabel('turn',fontsize=30)
        pyplot.xlabel('turn')
       # pyplot.ylabel('Ver $\Delta Q_{sc}/Q_{s}$',fontsize=30)
        #pyplot.ylabel(r'Ver $\Delta Q_{sc}/Q_{s}$',fontsize=30)
        #pyplot.ylabel(r'Ver $\Delta Q_{sc}$',fontsize=30)
        pyplot.ylabel(r'Ver $\Delta Q_{sc}$')
        #pyplot.xticks(fontsize=25)
        #pyplot.yticks(fontsize=25)
        pyplot.legend(loc=0)

    pyplot.tight_layout()
    pyplot.savefig("tuneshifts.png")
    pyplot.savefig("tuneshifts.svg")
    pyplot.show()



