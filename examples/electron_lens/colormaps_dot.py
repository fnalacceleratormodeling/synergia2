# -*- coding: utf-8 -*-
import os
import sys
import tables
#import matplotlib
from matplotlib  import pyplot
from mpl_toolkits.mplot3d import Axes3D as axes3d
from matplotlib import cm
import matplotlib.ticker as ticker 
from matplotlib  import font_manager
from matplotlib import rc
from matplotlib.colors import LogNorm
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
#rc('font',size=35)
rc('font',size=24)
#matplotlib.rcParams['xtick.major.size'] = 20
#matplotlib.rcParams['xtick.major.width'] = 4
#atplotlib.rcParams['xtick.minor.size'] = 10
#matplotlib.rcParams['xtick.minor.width'] = 2
import numpy
from tune_diagram import tune_diagram

if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise RuntimeError,"usage: colormap_dot.py input_file, coords(xy, zx, zy)  output_label"
    
    input_file=sys.argv[1]
    print "input file=",input_file
    axes=sys.argv[2]
    output_label=sys.argv[3]
    print "axes=",axes
    

    log_scale=0
    limit_z_axe=0
    # the resonance lines are getting drawn between 0 and 1 which
    # is not where the plot limits are
    show_resonance_lines=1
    show_cbar=1
    roder=2
    zminlim=0
    zmaxlim=50.0
    #zmaxlim=0.5e-5
    #zmaxlim=1.5e-5
   #zmaxlim=1
   
 
    #Q0x=0.723756 #hep_nodog
    #Q0y=0.829121 #hep_nodog
    #title="n6e10,HEP, nodog"
    
    #Q0x=0.73219  #hep
    #Q0y=0.828922  #hep
    #title="n7e10,HEP, no wake,lost"

    #Q0x=0.763345#Q0 nodog
    #Q0y=0.788129#Q0 nodog
    title="n6e10,HEP, no dog, ideal, ch=-17,-9"
 
    Q0x=3.74
    Q0y=3.74
    #Q0x=3.610729
    #Q0y=3.881919

    
    title=output_label

    npt=6
    Q0s=0.02
    
    
    if (axes=="xy"):
        bare_tunex=Q0x
        bare_tuney=Q0y       
        xlim_low=0.
        xlim_up=1
        ylim_low=0
        ylim_up=1
        #ylim_low=xlim_low      
        #ylim_up=xlim_up
        xlabel='$Q_x$'
        ylabel='$Q_y$'
    elif (axes=="zx"):
        bare_tunex=1-Q0s
        bare_tuney=Q0x
        xlim_low=0.5
        xlim_up=1
        ylim_low=0.5
        ylim_up=1.
        xlabel='$Q_s$'
        ylabel='$Q_x$'
    elif (axes=="zy"):    
        bare_tunex=1-Q0s        
        bare_tuney=Q0y
        xlim_low=0.5
        xlim_up=1
        ylim_low=0.5
        ylim_up=1.
        xlabel='$Q_s$'
        ylabel='$Q_y$'
    else:  raise RuntimeError,"unkown axes option input"    
        

    if log_scale:
       zminlim=0.001
       zmaxlim=30


    #xlim_low=0.45
    #xlim_up=1
    
    #ylim_low=xlim_low
    #ylim_up=xlim_up
    
    
    
    
    
    #bucket_length=5.64529
    #beta=0.720748
    #gamma=1.44259
    #lattice_length=474.204
    #number of measurements per turn=48

   

    #pctx1=0.2841
    #pcty1=0.30396 #y135dm0p11
    ##pctx1=0.2841
    ##pcty1=0.30396 #y150dm0p11
    #pctx2=0.2896
   #pcty2=0.30396#y119dm0p11

    #pctx3=0.30396 #x200dm0p11
    #pcty3=0.30396#x200dm0p11
    #pctx4=0.30396#x181dm0p11
    #pcty4=0.316078#x181dm0p11

    #pctx5=0.28304#x170dm0p11
    #pcty5=0.30396#x170dm0p11


    
     
   # print "input_file[-4:]=",input_file[-4:]
    if(input_file[-4:]==".dat"):
      wwx=[]
      wwy=[]
      rho=[] 
      foo = open(input_file, "r")
      for line in foo:     
        words=line.split()
        if (words[0] != "#"):
          wwx.append(float(words[0]))
          wwy.append(float(words[1]))
          rho.append(float(words[2])) 
      foo.close()
      ngrid=int(numpy.sqrt(len(wwx)))
      tune_density=numpy.zeros((ngrid,ngrid))     
      for iwx in range(ngrid):
        for iwy in range(ngrid):
          tune_density[iwx,iwy] = rho[iwx*ngrid+iwy]
        
    elif ( input_file[-4:]==".npy"):
      tune_density=numpy.load(input_file)
      print "tune_density shape=",tune_density.shape
      ngrid=tune_density.shape[0]


   
    print "nr grid points=",ngrid 
 
   


  
    deltax=xlim_up-xlim_low
    deltay=ylim_up-ylim_low
   
    #zmaxlim=0.00040
    print "ngrid: ", ngrid, ", xlim_low: ", xlim_low, ", xlim_up: ", xlim_up, ", ylim_low: ", ylim_low, ", ylim_up: ", ylim_up
    tune_density_cm=tune_density[numpy.meshgrid(numpy.arange(int(ngrid*xlim_low),int(ngrid*xlim_up)), numpy.arange(int(ngrid*ylim_low),int(ngrid*ylim_up)))]
    
################################################################
    dv=1./ngrid 
    #pyplot.figure(1)
    #pyplot.plot(numpy.arange(int(ngrid*xlim_low),int(ngrid*xlim_up))*dv,tune_density[iwc,int(ngrid*ylim_low):int(ngrid*ylim_up)]*(1.-rest)+tune_density[iwc+1,int(ngrid*ylim_low):int(ngrid*ylim_up)]*rest)
    #pyplot.plot(numpy.arange(int(ngrid*xlim_low),int(ngrid*xlim_up))*dv, tune_density[iwc,int(ngrid*ylim_low):int(ngrid*ylim_up)],'r-')
    #pyplot.plot(numpy.arange(int(ngrid*xlim_low),int(ngrid*xlim_up))*dv, tune_density[iwc+1,int(ngrid*ylim_low):int(ngrid*ylim_up)],'g-')
    
    
    

 
    pyplot.figure(2,figsize=(10.2,8))
   
    
    if log_scale:
      pyplot.imshow(tune_density_cm,origin='lower',extent=[xlim_low,npt*xlim_up,ylim_low,npt*ylim_up],norm=LogNorm())
    else:
      pyplot.imshow(tune_density_cm,origin='lower',extent=[xlim_low,npt*xlim_up,ylim_low,npt*ylim_up ])
   # pyplot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
      #pyplot.clim(zminlim,zmaxlim)


    
     
    #cbar= pyplot.colorbar(ticks=[0, 0.00015, 0.0003])
   # cbar.ax.set_yticklabels([0, 1.5, 3])
    

    #pyplot.xticks([xlim_low+deltax/5.,xlim_low+2*deltax/5.,xlim_low+3*deltax/5.,xlim_up-deltax/5.])
    #pyplot.yticks([ylim_low+deltay/5.,ylim_low+2*deltay/5.,ylim_low+3*deltay/5.,ylim_up-deltay/5.])
    xlow_show=3.0
    #xlow_show=3.5
    xup_show=4.0
    #xup_show=5
    ylow_show=3.0
    #ylow_show=3.5
    yup_show=4.0
    #yup_show=5
    pyplot.xlim(xlow_show, xup_show)
    pyplot.ylim(ylow_show, yup_show)
    #pyplot.xticks([2.5,3,3.5,4]) XXXX EGS
    #pyplot.yticks([2.5, 3,3.5,4]) XXX EGS
    pyplot.gca().xaxis.tick_top()
    
   # pyplot.xticks([0.26,0.29,0.32])
    #pyplot.yticks([0.26,0.29,0.32])
   # pyplot.gca().caxis.major.formatter.set_powerlimits((0,0))    
   # pyplot.plot([ch_tune-kmode*Qs,ch_tune-kmode*Qs],[ylim_low,ylim_up],'w--',linewidth=2)
   # pyplot.plot([ch_tune,ch_tune],[ylim_low,ylim_up],'w-',linewidth=2)
    #pyplot.plot([ch_tune-kmode*Qs-Qs,ch_tune-kmode*Qs-Qs],[ylim_low,ylim_up],'m--',linewidth=0.5)
    #pyplot.plot([ch_tune-kmode*Qs-2*Qs,ch_tune-kmode*Qs-2*Qs],[ylim_low,ylim_up],'m--',linewidth=0.5)
    #pyplot.plot([ch_tune-kmode*Qs+2*Qs,ch_tune-kmode*Qs+2*Qs],[ylim_low,ylim_up],'m--',linewidth=0.5)
   # plot_lines=True

   # if  not(tune_diff=="d1"):
      #pyplot.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs,ch_tune-kmode*Qs],'m-')
     # pyplot.plot([xlim_low,xlim_up],[ch_tune,ch_tune],'w-')
     # pyplot.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs-Qs,ch_tune-kmode*Qs-Qs],'m--',linewidth=0.5)
     # pyplot.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs-2*Qs,ch_tune-kmode*Qs-2*Qs],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up],[xlim_low,xlim_up],'m--',linewidth=0.5)
     # pyplot.plot([xlim_low+kmode*Qs,xlim_up],[xlim_low,xlim_up-kmode*Qs],'m--',linewidth=0.5)
     # pyplot.plot([xlim_low,xlim_up-kmode*Qs],[xlim_low+kmode*Qs,xlim_up],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low+kmode*Qs,xlim_up],[xlim_low,xlim_up-kmode*Qs],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up],[ch_tune,ch_tune],'m--',linewidth=0.5)
     # pyplot.plot([xlim_low,xlim_up],[ch_tune-2*Qs,ch_tune-2*Qs],'m--',linewidth=0.5)
     # pyplot.plot([ch_tune-kmode*Qs-Qs,ch_tune-kmode*Qs-Qs],[ylim_low,ylim_up],'m--',linewidth=0.5)
     
    #pyplot.plot([xlim_low,xlim_up],[bare_tuney-kmode*Qs,bare_tuney-kmode*Qs],'g-')
    #pyplot.plot([xlim_low,xlim_up],[ch_tune,ch_tune],'g-')
    pyplot.plot([bare_tunex],[bare_tuney],'wD',markersize=8)
    if show_resonance_lines:
        print "roder: ", roder
        tune_diagram(roder,['black', 'red', 'cyan', 'white', 'orange', 'green'], linew=1, xbias=3, ybias=3)

   # pyplot.axis([xlim_low,npt*xlim_up,ylim_low,npt*ylim_up])
    pyplot.axis([xlow_show, xup_show,ylow_show, yup_show])
    
    plot_lines=True
    if plot_lines:
       pyplot.plot([xlow_show, xup_show],[ylow_show, yup_show],'w--',linewidth=0.5) 
       pyplot.plot([3.5, 3.5],[ylow_show, yup_show],'w--',linewidth=0.5)
       pyplot.plot([xlow_show, xup_show],[3.5,3.5],'w--',linewidth=0.5)
       pyplot.plot([xlow_show, xup_show], [3.0, 3.0], 'w--', linewidth=0.5)
       pyplot.plot([3.0, 3.0], [ylow_show, yup_show], 'w--', linewidth=0.5)
       pyplot.plot([4.0, 4.0], [ylow_show, yup_show], 'w--', linewidth=0.5)
       pyplot.plot([4.5, 4.5], [ylow_show, yup_show], 'w--', linewidth=0.5)
       pyplot.plot([xlow_show, xup_show], [4.0, 4.0], 'w--', linewidth=0.5)
       pyplot.plot([xlow_show, xup_show], [4.5, 4.5], 'w--', linewidth=0.5)
    #  pyplot.plot([xlim_low,xlim_up],[xlim_low,xlim_up],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low+kmode*Qs,xlim_up],[xlim_low,xlim_up-kmode*Qs],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up-kmode*Qs],[xlim_low+kmode*Qs,xlim_up],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low+2.*kmode*Qs,xlim_up],[xlim_low,xlim_up-2.*kmode*Qs],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up-2.*kmode*Qs],[xlim_low+2.*kmode*Qs,xlim_up],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low+3.*kmode*Qs,xlim_up],[xlim_low,xlim_up-3.*kmode*Qs],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up-3.*kmode*Qs],[xlim_low+3.*kmode*Qs,xlim_up],'m--',linewidth=0.5)
      #pyplot.plot([ch_tune-kmode*Qs+2*Qs,ch_tune-kmode*Qs+2*Qs],[xlim_low,xlim_up],'m--',linewidth=0.5)      
      #pyplot.plot([ch_tune-kmode*Qs-Qs,ch_tune-kmode*Qs-Qs],[xlim_low,xlim_up],'m--',linewidth=0.5)
      #pyplot.plot([ch_tune-kmode*Qs-2*Qs,ch_tune-kmode*Qs-2*Qs],[xlim_low,xlim_up],'m--',linewidth=0.5)
      #pyplot.plot([ch_tune-kmode*Qs-3*Qs,ch_tune-kmode*Qs-3*Qs],[xlim_low,xlim_up],'m--',linewidth=0.5)
      #pyplot.plot([ch_tune-kmode*Qs-4*Qs,ch_tune-kmode*Qs-4*Qs],[xlim_low,xlim_up],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs+2*Qs,ch_tune-kmode*Qs+2*Qs],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs+Qs,ch_tune-kmode*Qs+Qs],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs-Qs,ch_tune-kmode*Qs-Qs],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs-2*Qs,ch_tune-kmode*Qs-2*Qs],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs-3*Qs,ch_tune-kmode*Qs-3*Qs],'m--',linewidth=0.5)
      #pyplot.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs-4*Qs,ch_tune-kmode*Qs-4*Qs],'m--',linewidth=0.5)

    #plot_particles=False
    #if plot_particles:
      ##pyplot.plot([pctx1],[pcty1],'bo',markersize=6)
      ##pyplot.plot([pctx2],[pcty2],'go',markersize=6)
      #pyplot.plot([pctx3],[pcty3],'ro',markersize=6)
      #pyplot.plot([pctx4],[pcty4],'ko',markersize=6)
      #pyplot.plot([pctx5],[pcty5],'mo',markersize=6)
   
    #zmaxlim=3.5 # for d0 0.5%
    #zmaxlim=5
    #zmaxlim=0.035
    pyplot.suptitle(title,fontsize=20)
    #pyplot.clim(zminlim,zmaxlim)
    if limit_z_axe:
        pyplot.clim(zminlim,zmaxlim)
    if log_scale:       
       pyplot.colorbar(ticks=[0.001,0.01,0.1, 1])
    
     # pyplot.colorbar(ticks=[10, 20, 30])#.set_label('arb. units', rotation=90)
      #pyplot.colorbar(ticks=[0 ,1, 2, 3])#.set_label('arb. units', rotation=90)
      #pyplot.colorbar(ticks=[0 ,1.5, 3])#.set_label('arb. units', rotation=90)
      #pyplot.colorbar(ticks=[0 ,0.03])#.set_label('arb. units', rotation=90)
    if show_cbar:  
        cbar= pyplot.colorbar()
   # cbar.set_label('arb. units', rotation=90)
    #pyplot.xlabel(xlabel,fontsize=30)
    pyplot.xlabel(xlabel,fontsize=20)
    #pyplot.ylabel(ylabel,fontsize=30)
    pyplot.ylabel(ylabel,fontsize=20)
################################################################       
    #tx_grid, ty_grid =numpy.meshgrid(numpy.arange(int(ngrid*xlim_low),int(ngrid*xlim_up)), numpy.arange(int(ngrid*ylim_low),int(ngrid*ylim_up)))    
    #txpoints= tx_grid/(1.*ngrid)
    #typoints= ty_grid/(1.*ngrid)
    #ax = axes3d(pyplot.figure(3))
    ##ax.set_zlim3d(zminlim,zmaxlim)
    #ax.plot_surface(txpoints, typoints,tune_density_cm, rstride=1, cstride=1, cmap='RdBu',
            #linewidth=0, antialiased=False)
   ## ax.plot_wireframe(txpoints, typoints,tune_density_cm, rstride=1, cstride=1)
    #ax.plot([bare_tunex],[bare_tuney],zs=0.0 , c='g',marker='o',markersize=6)
    #ax.plot([ch_tune-kmode*Qs,ch_tune-kmode*Qs],[ylim_low,ylim_up], zs=0.0, c='g',linewidth=1)
    #ax.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs,ch_tune-kmode*Qs],zs=0.0, c='g',linewidth=1)
    #ax.set_xlabel("tune x")
    #ax.set_ylabel("tune y")
   
    ##ax.view_init(elev=45, azim=50)
    #ax.view_init(elev=45, azim=-55)
    ##ax.view_init(elev=0, azim=-90)
    #if plot_lines:
      #ax.plot([xlim_low,xlim_up],[xlim_low,xlim_up], zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low+kmode*Qs,xlim_up],[xlim_low,xlim_up-kmode*Qs],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low,xlim_up-kmode*Qs],[xlim_low+kmode*Qs,xlim_up],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low+2.*kmode*Qs,xlim_up],[xlim_low,xlim_up-2.*kmode*Qs],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low,xlim_up-2.*kmode*Qs],[xlim_low+2.*kmode*Qs,xlim_up],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low+3.*kmode*Qs,xlim_up],[xlim_low,xlim_up-3.*kmode*Qs],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low,xlim_up-3.*kmode*Qs],[xlim_low+3.*kmode*Qs,xlim_up],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([ch_tune-kmode*Qs+2*Qs,ch_tune-kmode*Qs+2*Qs],[xlim_low,xlim_up],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([ch_tune-kmode*Qs+Qs,ch_tune-kmode*Qs+Qs],[xlim_low,xlim_up],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([ch_tune-kmode*Qs-Qs,ch_tune-kmode*Qs-Qs],[xlim_low,xlim_up],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([ch_tune-kmode*Qs-2*Qs,ch_tune-kmode*Qs-2*Qs],[xlim_low,xlim_up],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([ch_tune-kmode*Qs-3*Qs,ch_tune-kmode*Qs-3*Qs],[xlim_low,xlim_up],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([ch_tune-kmode*Qs-4*Qs,ch_tune-kmode*Qs-4*Qs],[xlim_low,xlim_up],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs+2*Qs,ch_tune-kmode*Qs+2*Qs],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs+Qs,ch_tune-kmode*Qs+Qs],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs-Qs,ch_tune-kmode*Qs-Qs],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs-2*Qs,ch_tune-kmode*Qs-2*Qs],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs-3*Qs,ch_tune-kmode*Qs-3*Qs],zs=0.0, c='y',linewidth=0.5)
      #ax.plot([xlim_low,xlim_up],[ch_tune-kmode*Qs-4*Qs,ch_tune-kmode*Qs-4*Qs],zs=0.0, c='y',linewidth=0.5)
#################################################################   
    x_tune_density=numpy.zeros(ngrid)
    total_weight=0.
    for iwx in range(ngrid):
       for iwy in range(ngrid):
          x_tune_density[iwx] +=  tune_density[iwx,iwy]
          total_weight += tune_density[iwx,iwy]

    

    print "total weight=",total_weight
    fo0 = open("DOT_X_"+output_label+".dat", "w")
    fo0.write("# tune_x   XDOT[:]    ")
    fo0.write("\n")   
    for it in range(int(ngrid*0.5)):
      fo0.write(str(it*dv))
      fo0.write("    ")
      fo0.write(str(x_tune_density[it]))
      fo0.write("\n") 
    fo0.close() 

#################################################################   
    y_tune_density=numpy.zeros(ngrid)
    total_weight=0.
    for iwy in range(ngrid):
       for iwx in range(ngrid):
          y_tune_density[iwy] +=  tune_density[iwx,iwy]
          total_weight += tune_density[iwx,iwy]

    

    print "total weight=",total_weight
    fo0 = open("DOT_Y_"+output_label+".dat", "w")
    fo0.write("# tune_xy  XDOT[:]    ")
    fo0.write("\n")   
    for it in range(int(ngrid*0.5)):
      fo0.write(str(it*dv))
      fo0.write("    ")
      fo0.write(str(y_tune_density[it]))
      fo0.write("\n") 
    fo0.close() 



################################################################ 
    pyplot.savefig("tuneplot-"+output_label+".png")
    pyplot.show()

  
