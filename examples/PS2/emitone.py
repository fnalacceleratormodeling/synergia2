#!/usr/bin/env python

from matplotlib  import pyplot
import tables
#import options
import sys
import numpy
    
if __name__ == "__main__":
    #opts = options.Options("emit")
    #opts.add("dir","point2","job directory")
    #opts.parse_argv(sys.argv)
    
    
    table =  tables.openFile("nosp/ps2_output.h5")
    #table =  tables.openFile("%s/tmp-diagnostics.h5" % opts.get("dir"),"r")
    class Empty: pass
    diags = Empty()
    for node in table.root:
        setattr(diags,node.name,node.read())
    table.close()
    pyplot.plot(diags.s/1000,diags.emitxyz/diags.emitxyz[0],label='I=0.')
   # pyplot.plot(diags.s/1000,diags.emity/diags.emity[0],label='I=0.')
   # pyplot.plot(diags.s[0:len(diags.s):60]/1000,diags.emitx[0:len(diags.s):60]/diags.emitx[0],label='I=0.')
   # pyplot.plot(diags.s[0:len(diags.s):60]/1000,diags.emitxy[0:len(diags.s):60]/diags.emitxy[0],label='I=0.')
    
    table =  tables.openFile("sp_k60/ps2_output.h5")
    #table =  tables.openFile("%s/tmp-diagnostics.h5" % opts.get("dir"),"r")
    class Empty: pass
    diags = Empty()
    for node in table.root:
        setattr(diags,node.name,node.read())
    table.close()
    pyplot.plot(diags.s/1000,diags.emitxyz/diags.emitxyz[0],label='I=4.2e+11, kicks=60')
   #pyplot.plot(diags.s/1000,diags.emity/diags.emity[0],label='I=4.2e+11, kicks=60')
   #pyplot.plot(diags.s[0:len(diags.s):60]/1000,diags.emitx[0:len(diags.s):60]/diags.emitx[0],label='I=4.2e+11, kicks=60')
   # pyplot.plot(diags.s[0:len(diags.s):60]/1000,diags.emitxy[0:len(diags.s):60]/diags.emitxy[0],label='I=4.2e+11, kicks=60')
   
    
    table =  tables.openFile("sp_k100/tmp_output-03.h5")
    #table =  tables.openFile("%s/tmp-diagnostics.h5" % opts.get("dir"),"r")
    class Empty: pass
    diags = Empty()
    for node in table.root:
        setattr(diags,node.name,node.read())
    table.close()
    pyplot.plot(diags.s/1000,diags.emitxyz/diags.emitxyz[0],label='I=4.2e+11, kicks=100')
   #pyplot.plot(diags.s/1000,diags.emity/diags.emity[0],label='I=4.2e+11, kicks=100')
   #pyplot.plot(diags.s[0:len(diags.s):100]/1000,diags.emitx[0:len(diags.s):100]/diags.emitx[0],label='I=4.2e+11, kicks=100')
   # pyplot.plot(diags.s[0:len(diags.s):100]/1000,diags.emitxy[0:len(diags.s):100]/diags.emitxy[0],label='I=4.2e+11, kicks=100')
    
    
    
  
    table =  tables.openFile("sp_k200/ps2_output.h5")
    class Empty: pass
    diags = Empty()
    for node in table.root:
        setattr(diags,node.name,node.read())
    table.close()
    pyplot.plot(diags.s/1000,diags.emitxyz/diags.emitxyz[0],label='I=4.2e+11, kicks=200')
  # pyplot.plot(diags.s/1000,diags.emity/diags.emity[0],label='I=4.2e+11, kicks=200')
   # pyplot.plot(diags.s[0:len(diags.s):200]/1000,diags.emitx[0:len(diags.s):200]/diags.emitx[0],label='I=4.2e+11, kicks=200')
   # pyplot.plot(diags.s[0:len(diags.s):200]/1000,diags.emitxy[0:len(diags.s):200]/diags.emitxy[0],label='I=4.2e+11, kicks=200')
    
    table =  tables.openFile("sp_k400/tmp_output-13.h5")
    class Empty: pass
    diags = Empty()
    for node in table.root:
        setattr(diags,node.name,node.read())
    table.close()
    pyplot.plot(diags.s/1000,diags.emitxyz/diags.emitxyz[0],label='I=4.2e+11, kicks=400')
   # pyplot.plot(diags.s/1000,diags.emity/diags.emity[0],label='I=4.2e+11, kicks=400')
   # pyplot.plot(diags.s[0:len(diags.s):400]/1000,diags.emitx[0:len(diags.s):400]/diags.emitx[0],label='I=4.2e+11, kicks=400')
   # pyplot.plot(diags.s[0:len(diags.s):400]/1000,diags.emitxy[0:len(diags.s):400]/diags.emitxy[0],label='I=4.2e+11, kicks=400')
    
    
    pyplot.legend(loc=0)
    pyplot.xlabel('s(m x 1.0e3)')
    pyplot.ylabel('emitXYZ/emitXYZ(s=0)')
    pyplot.title('Emittance XYZ, 1000 turns', fontsize=18)
    #pyplot.ylabel('emitY/emitY(s=0)')
    #pyplot.title('Emittance Y, 1000 turns', fontsize=18)
    #pyplot.ylabel('emitX/emitX(s=0)')
    #pyplot.title('Emittance X, 1000 turns', fontsize=18)
    #pyplot.ylabel('emitXY/emitXY(s=0)')
    #pyplot.title('Emittance XY, 1000 turns', fontsize=18)
    
    pyplot.show()
