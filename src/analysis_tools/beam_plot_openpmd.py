import sys
import argparse
#from openpmd_api import io
#import matplotlib.pyplot as plt
from dataclasses import dataclass
from enum import Enum

class coords(Enum):
    x = 0
    xp = 1 
    y = 2
    yp = 3
    cdt = 4 
    dpop = 5
    pz = 7 
    energy = 8
    t = 9 
    z = 10 

@dataclass
class Options:
    hist: bool
    hcoord: coords
    vcoord: coords
    inputfile: str
    outputfile: str
    num_contour: int  
    show: bool = True
    iteration: int = 0
    bins: int = 10
    minh: float =-sys.float_info.max
    maxh: float = sys.float_info.max
    minv: float = -sys.float_info.max
    maxv: float = sys.float_info.max
    contour: bool = False

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="OpenPMD-series-filename",
                    type=str)



