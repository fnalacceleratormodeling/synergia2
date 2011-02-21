#!/usr/bin/env python

import sys
import numpy
from pygsl import sf
from math import exp, sqrt, pi
from matplotlib import pyplot
import tables

# spherical gaussian charge density
def gauss_density(x, y, z):
    global Q, sigma
    r2 = x * x + y * y + z * z
    return Q * exp(-r2 / (2 * sigma ** 2)) / (sigma * sqrt(2 * pi)) ** 3

# electric potential due to gauss_density
def gauss_field(x, y, z):
    global Q, sigma, epsilon0
    r = sqrt(x * x + y * y + z * z)
    if r == 0.0:
        return Q / (4.0 * pi * epsilon0) * sqrt(2.0) / (sqrt(pi) * sigma)
    return Q / (4.0 * pi * epsilon0 * r) * sf.erf(r / (sqrt(2.0) * sigma))[0];

c = 299792458.0
mu0 = 4 * pi * 1.0e-7
epsilon0 = 1.0 / (c * c * mu0)
e = 1.6021892e-19
Q = 1.7e11 * e

n_sigma = 8.0
sigma = 1.3e-3
n = 16
n2 = 2 * n
size = sigma * n_sigma
h = size / n
left = -size / 2.0 + h / 2.0

def calculate(coeff, plot=False):
    rho2 = numpy.zeros([n2, n2, n2], numpy.float64)
    for i in range(0, n):
        x = left + i * h
        for j in range(0, n):
            y = left + j * h
            for k in range(0, n):
                z = left + k * h
                rho2[i, j, k] = gauss_density(x, y, z)

    G2 = numpy.zeros([n2, n2, n2], numpy.float64)
    G000 = coeff / h
    for i in range(0, n + 1):
        dx = i * h
        mi = n2 - i
        if (mi == n2):
            mi = i
        for j in range(0, n + 1):
            dy = j * h
            mj = n2 - j
            if (mj == n2):
                mj = j
            for k in range(0, n + 1):
                dz = k * h
                mk = n2 - k
                if (mk == n2):
                    mk = k
                if((i == 0) and (j == 0) and (k == 0)):
                    G2[0, 0, 0] = G000
                else:
                    G2[i, j, k] = 1.0 / sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                # seven mirror images
                G2[mi, j, k] = G2[i, j, k]
                G2[mi, mj, k] = G2[i, j, k]
                G2[mi, mj, mk] = G2[i, j, k]
                G2[mi, j, mk] = G2[i, j, k]
                G2[i, mj, k] = G2[i, j, k]
                G2[i, mj, mk] = G2[i, j, k]
                G2[i, j, mk] = G2[i, j, k]

    rho2hat = numpy.fft.fftn(rho2)
    G2hat = numpy.fft.fftn(G2)
    phi2hat = numpy.zeros([n2, n2, n2], numpy.complex64)
    for i in range(0, n2):
        for j in range(0, n2):
            for k in range(0, n2):
                phi2hat[i, j, k] = rho2hat[i, j, k] * G2hat[i, j, k]

    phi2 = numpy.fft.ifftn(phi2hat)

    phi = numpy.zeros([n, n, n], numpy.float64)
    phie = numpy.zeros([n, n, n], numpy.float64)
    error = numpy.zeros([n, n, n], numpy.float64)

    for i in range(0, n):
        x = left + i * h
        for j in range(0, n):
            y = left + j * h
            for k in range(0, n):
                z = left + k * h
                phi[i, j, k] = phi2[i, j, k].real * h * h * h / (4 * pi * epsilon0)
                phie[i, j, k] = gauss_field(x, y, z)
                error[i, j, k] = (phie[i, j, k] - phi[i, j, k]) / phie[i, j, k]


    max_error = numpy.max(numpy.max(numpy.max(error)))
    min_error = numpy.min(numpy.min(numpy.min(error)))
    mean_error = numpy.mean(numpy.mean(numpy.mean(abs(error))))

    if plot:
        # comparison point
        p = n / 2
        pyplot.plot(phi[:, p, p], 'o', label='result')
        pyplot.plot(phie[:, p, p], label='expected')
        pyplot.show()

    return max_error, min_error, mean_error


if __name__ == "__main__":
    find_coeff = False
    if find_coeff:
        cs = []
        maxes = []
        mines = []
        meanes = []
        for i in range(0, 20):
            coeff = 2.65 + i * 0.05
            maxe, mine, meane = calculate(coeff)
            cs.append(coeff)
            maxes.append(maxe)
            mines.append(mine)
            meanes.append(meane)
        nmeanes = numpy.array(meanes, numpy.float64)
        pyplot.plot(cs, maxes, label='max(error)')
        pyplot.plot(cs, mines, label='min(error)')
        pyplot.plot(cs, nmeanes * 10, label='mean(abs(error))*10')
        pyplot.legend()
        pyplot.show()
    else:
        coeff = 2.8  # best fit value using find_coeff search above
        max_error, min_error, mean_error = calculate(coeff, True)
        print coeff, max_error, min_error, mean_error


