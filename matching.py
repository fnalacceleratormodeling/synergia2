#!/usr/bin/env python

from math import sqrt

def match_twiss_width(width,alpha,beta):
    """Calculate input parameters for a matched beam of given width
    using Courant-Snyder (Twiss) parameters. Returns
        (width_prime,r,emittance),
    where width_prime is the width in the conjugate coordinate and r
    is the correlation coefficient."""
    gamma = (1+alpha**2)/beta
    emittance = width**2/beta
    width_prime = sqrt(gamma*emittance)
    r = -alpha/sqrt(1+alpha**2)
    return (width_prime,r,emittance)

def match_twiss_emittance(emittance,alpha,beta):
    """Calculate input parameters for a matched beam of given width
    using Courant-Snyder (Twiss) parameters. Returns
        (width,width_prime,r),
    where width and width_prime are the width in the coordinate and
    its conjugate and r is the correlation coefficient."""
    gamma = (1+alpha**2)/beta
    width = sqrt(beta*emittance)
    width_prime = sqrt(gamma*emittance)
    r = -alpha/sqrt(1+alpha**2)
    return (width,width_prime,r)
