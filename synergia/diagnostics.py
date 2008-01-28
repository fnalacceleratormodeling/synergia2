#!/usr/bin/env python
import Numeric
import string

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

class Diagnostics_impact_orig:
    def __init__(self, dirname = "."):
        self.dirname = dirname
        self._read_s_mean_std()
        self._read_emit()
        self._read_corr()
        self._read_num_part()

    def _read_s_mean_std(self):
        # read the number of lines in fort.24
        num = len(open("%s/fort.24" % self.dirname ,"r").readlines())
        self.num = num
        self.s = Numeric.zeros(num,'d')
        self.mean = Numeric.zeros((num,6),'d')
        self.std = Numeric.zeros((num,6),'d')
        for offset in range(0,3):
            f = open("%s/fort.%d" % (self.dirname,24 + offset),"r")
            line = f.readline()
            index = 0
            while line:
                cols = string.split(line)
                try:
                    self.s[index] = float(cols[0])
                    self.mean[index,offset*2] = float(cols[1])
                    self.std[index,offset*2] = float(cols[2])
                    self.mean[index,offset*2+1] = float(cols[3])
                    self.std[index,offset*2+1] = float(cols[4])
                    index += 1
                except:
                    pass
                line = f.readline()
            f.close()

    def _read_emit(self):
        num = self.num
        self.emitx = Numeric.zeros(num, 'd')
        self.emity = Numeric.zeros(num, 'd')
        self.emitz = Numeric.zeros(num, 'd')
        self.emit4d = Numeric.zeros(num, 'd')
        self.emit6d = Numeric.zeros(num, 'd')
        f = open("%s/fort.31" % self.dirname, "r")
        line = f.readline()
        index = 0
        while line:
            cols = string.split(line)
            try:
                self.emitx[index] = float(cols[0])
                self.emity[index] = float(cols[1])
                self.emitz[index] = float(cols[2])
                self.emit4d[index] = float(cols[3])
                self.emit6d[index] = float(cols[4])
                index += 1
            except:
                pass
            line = f.readline()
        f.close()

    def _read_corr(self):
        num = self.num
        self.corr = Numeric.zeros((num,15),'d')
        f = open("%s/fort.32" % self.dirname, "r")
        line = f.readline()
        index = 0
        while line:
            cols = string.split(line)
            try:
                for i in range(0,15):
                    self.corr[index,i] = float(cols[i])
                index += 1
            except:
                pass
            line = f.readline()
        f.close()

    def _read_num_part(self):
        num = self.num
        self.num_part = Numeric.zeros(num,'d')
        f = open("%s/fort.28" % self.dirname, "r")
        line = f.readline()
        index = 0
        while line:
            try:
                cols = string.split(line)
                self.num_part[index] = float(cols[3])
                index += 1
            except:
                pass
            line = f.readline()
        f.close()
        
