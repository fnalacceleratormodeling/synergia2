#!/usr/bin/env python

import physics_constants
import math
import numpy

class Beam_parameters:
    def __init__(self, mass_GeV, charge_e, kinetic_energy_GeV,
                 initial_phase_rad=0, scaling_frequency_Hz=1, transverse=0, adjust_zlength_to_freq=1):
        self.mass_GeV = mass_GeV
        self.charge_e = int(charge_e)
        self.kinetic_energy_GeV = kinetic_energy_GeV
        self.initial_phase_rad = initial_phase_rad
        self.scaling_frequency_Hz = scaling_frequency_Hz
        self.default_correlation_coeffs()
        self.sigma_x_m = None
        self.lambda_x_GeVoc = None
        self.sigma_y_m = None
        self.lambda_y_GeVoc = None
        self.sigma_z_m = None
        self.lambda_z_GeVoc = None
        self.transverse = transverse
        self.adjust_zlength_to_freq=adjust_zlength_to_freq
        self.x_params(1,1)
        self.y_params(1, 1)
        self.z_params(1, 1, 1)

    def x_params(self,sigma,lam,
                 r=None,mismatch=1,mismatch_p=1,offset=0,offset_p=0):
	"""The correlation coefficient r can also be set by the xpx argument
	to correlation_coeffs"""
        self.sigma_x_m = sigma
        self.lambda_x_GeVoc = lam
        if r != None: self.xpx = r
        self.mismatch_x = mismatch
        self.mismatch_px = mismatch_p
        self.offset_x_m = offset
        self.offset_px = offset_p
        
    def y_params(self,sigma,lam,
                 r=None,mismatch=1,mismatch_p=1,offset=0,offset_p=0):
	"""The correlation coefficient r can also be set by the ypy argument
	to correlation_coeffs"""
        self.sigma_y_m = sigma
        self.lambda_y_GeVoc = lam
        if r != None: self.ypy = r
        self.mismatch_py = mismatch_p
        self.offset_y_m = offset
        self.offset_py = offset_p
        
    def z_params(self,sigma,lam, z_length=None,
                 r=None,mismatch=1,mismatch_p=1,offset=0,offset_p=0,
                 num_peaks=1):
	"""The correlation coefficient r can also be set by the zpz argument
	to correlation_coeffs, zlength is used to make a  periodic bunch of length z_length 
    or can be the length of a transverse bunch"""
        self.sigma_z_m = sigma
        self.lambda_z_GeVoc = lam
        if r !=None : self.zpz = r
        self.mismatch_z = mismatch
        self.mismatch_pz = mismatch_p
        self.offset_z = offset
        self.offset_pz = offset_p
        self.num_zpeaks = num_peaks

        self.z_length=z_length
        if self.adjust_zlength_to_freq: 
            self.z_length = 2.0*math.pi*self.get_beta()*physics_constants.PH_MKS_c/self.get_omega()
        elif self.transverse:
            if z_length==None: raise RuntimeError, "please provide z_length for a transverse beam"
            
        



    def get_mass(self):
        return self.mass_GeV
    
    def get_kinetic_energy(self):
        return self.kinetic_energy_GeV
    
    def get_charge(self):
        return self.charge_e
    
    def get_z_length(self):
        return self.z_length
    
    def get_t_length(self): 
        if (self.z_length==None):
            t_length=None
        else:
           (Cxy, Cxpyp, Cz, Czp) = self.get_conversions()
           t_length=self.z_length*Cz
        return t_length
        
    
    def get_z_peaks(self):
        return self.num_zpeaks

    def get_omega(self):
        return self.scaling_frequency_Hz * 2* math.pi

    def get_gamma(self):
        return self.kinetic_energy_GeV/self.mass_GeV + 1.0

    def get_beta(self):
        """returns the velocity of the reference particle in natural units"""
        gamma = self.get_gamma()
        return math.sqrt(1.0 - 1.0/(gamma*gamma))

    def get_conversions(self):
        c = physics_constants.PH_MKS_c # m/s
        w = self.get_omega() # rad/s
        beta = self.get_beta()
        m = self.mass_GeV # GeV
        # Unit conversion: X^impact_i = C_i X^real_i
        Cxy   = w/c
        Cz    = w/(beta*c)
        Cxpyp = 1.0/m
        Czp   = 1.0/m # really a conversion from p_z to p_t, so includes an
                    # extra factor of beta, i.e., 1/(beta*m) * beta     
                    # !!  space charge and impedance kicks assume Cxpyp=Czp, be careful when/if you change this...
        return (Cxy, Cxpyp, Cz, Czp)
    
    def default_correlation_coeffs(self):
	"""Sets all correlation coefficients to zero."""
        self.xpx = 0.0
        self.xy = 0.0
        self.pxy = 0.0
        self.xpy = 0.0
        self.pxpy = 0.0
        self.ypy = 0.0
        self.xz = 0.0
        self.pxz = 0.0
        self.yz = 0.0
        self.pyz = 0.0
        self.xpz = 0.0
        self.pxpz = 0.0
        self.ypz = 0.0
        self.pypz = 0.0
        self.zpz = 0.0
	
    def correlation_coeffs(self, xpx=None, xy=None, pxy=None,
                           xpy=None, pxpy=None, ypy=None, xz=None,
                           pxz=None, yz=None, pyz=None, xpz=None,
                           pxpz=None, ypz=None, pypz=None, zpz=None):
        """Correlation coefficients for use the with gaussian covariance
        distribution"""
        if xpx: self.xpx = xpx
        if xy: self.xy = xy
        if pxy: self.pxy = pxy
        if xpy: self.xpy = xpy
        if pxpy: self.pxpy = pxpy
        if ypy: self.ypy = ypy
        if xz: self.xz = xz
        if pxz: self.pxz = pxz
        if yz: self.yz = yz
        if pyz: self.pyz = pyz
        if xpz: self.xpz = xpz
        if pxpz: self.pxpz = pxpz
        if ypz: self.ypz = ypz
        if pypz: self.pypz = pypz
        if zpz: self.zpz = zpz
    
    def get_covariances(self):
        c  = numpy.zeros((6,6),'d')
        (Cxy, Cxpyp, Cz, Czp) = self.get_conversions()
        # Unit conversion: X^impact_i = C_i X^real_i
        Cx = Cxy
        Cy = Cxy
        #Cz = Cz
        Cxp = Cxpyp
        Cyp = Cxpyp
        #Czp = Czp
        
        c[0,0] = self.sigma_x_m**2 * Cx**2
        c[0,1] = c[1,0] = self.sigma_x_m*self.lambda_x_GeVoc*self.xpx*Cx*Cxp
        c[1,1] = self.lambda_x_GeVoc**2 * Cxp**2
        c[0,2] = c[2,0] = self.sigma_x_m*self.sigma_y_m*self.xy*Cx*Cy
        c[1,2] = c[2,1] = self.lambda_x_GeVoc*self.sigma_y_m*self.pxy*Cxp*Cy
        c[2,2] = self.sigma_y_m**2 * Cy**2
        c[0,3] = c[3,0] = self.sigma_x_m*self.lambda_y_GeVoc*self.xpy*Cx*Cyp
        c[1,3] = c[3,1] = self.lambda_x_GeVoc*self.lambda_y_GeVoc*self.pxpy*Cxp*Cyp
        c[2,3] = c[3,2] = self.sigma_y_m*self.lambda_y_GeVoc*self.ypy*Cy*Cyp
        c[3,3] = self.lambda_y_GeVoc**2 * Cyp**2
        c[0,4] = c[4,0] = self.sigma_x_m*self.sigma_z_m*self.xz*Cx*Cz
        c[1,4] = c[4,1] = self.lambda_x_GeVoc*self.sigma_z_m*self.pxz*Cxp*Cz
        c[2,4] = c[4,2] = self.sigma_y_m*self.sigma_z_m*self.yz*Cy*Cz
        c[3,4] = c[4,3] = self.lambda_y_GeVoc*self.sigma_z_m*self.pyz*Cyp*Cz
        c[4,4] = self.sigma_z_m**2 * Cz**2
        c[0,5] = c[5,0] = self.sigma_x_m*self.lambda_z_GeVoc*self.xpz*Cx*Czp
        c[1,5] = c[5,1] = self.lambda_x_GeVoc*self.lambda_z_GeVoc*self.pxpz*Cxp*Czp
        c[2,5] = c[5,2] = self.sigma_y_m*self.lambda_z_GeVoc*self.ypz*Cy*Czp
        c[3,5] = c[5,3] = self.lambda_y_GeVoc*self.lambda_z_GeVoc*self.pypz*Cyp*Czp
        c[4,5] = c[5,4] = self.sigma_z_m*self.lambda_z_GeVoc*self.zpz*Cz*Czp
        c[5,5] = self.lambda_z_GeVoc**2 * Czp**2
        
        return c

    def get_means(self):
        (Cxy, Cxpyp, Cz, Czp) = self.get_conversions()
        # Unit conversion: X^impact_i = C_i X^real_i
        Cx = Cxy
        Cy = Cxy
        #Cz = Cz
        Cxp = Cxpyp
        Cyp = Cxpyp
        #Czp = Czp

        return numpy.array(
                [self.offset_x_m * Cx,
                 self.offset_px * Cxp,
                 self.offset_y_m * Cy,
                 self.offset_py * Cyp,
                 self.offset_z * Cz,
                 self.offset_pz * Czp],'d')
                 
    def get_distparam(self):
        (Cxy, Cxpyp, Cz, Czp) = self.get_conversions()
        # Unit conversion: X^impact_i = C_i X^real_i
        Cx = Cxy
        Cy = Cxy
        #Cz = Cz
        Cxp = Cxpyp
        Cyp = Cxpyp
        #Czp = Czp
        if self.transverse:
            dist_num = 667.0
        else:
            dist_num = 666.0
        param = [self.sigma_x_m**2 * Cx**2,
                 (self.sigma_x_m*self.lambda_x_GeVoc*self.xpx)*Cx*Cxp,
                 self.lambda_x_GeVoc**2 * Cxp**2,
                 (self.sigma_x_m*self.sigma_y_m*self.xy)*Cx*Cy,
                 (self.lambda_x_GeVoc*self.sigma_y_m*self.pxy)*Cxp*Cy,
                 self.sigma_y_m**2 * Cy**2,
                 (self.sigma_x_m*self.lambda_y_GeVoc*self.xpy)*Cx*Cyp,
                 (self.lambda_x_GeVoc*self.lambda_y_GeVoc*self.pxpy)*Cxp*Cyp,
                 (self.sigma_y_m*self.lambda_y_GeVoc*self.ypy)*Cy*Cyp,
                 self.lambda_y_GeVoc**2 * Cyp**2,
                 (self.sigma_x_m*self.sigma_z_m*self.xz)*Cx*Cz,
                 (self.lambda_x_GeVoc*self.sigma_z_m*self.pxz)*Cxp*Cz,
                 (self.sigma_y_m*self.sigma_z_m*self.yz)*Cy*Cz,
                 (self.lambda_y_GeVoc*self.sigma_z_m*self.pyz)*Cyp*Cz,
                 self.sigma_z_m**2 * Cz**2,
                 (self.sigma_x_m*self.lambda_z_GeVoc*self.xpz)*Cx*Czp,
                 (self.lambda_x_GeVoc*self.lambda_z_GeVoc*self.pxpz)*Cxp*Czp,
                 (self.sigma_y_m*self.lambda_z_GeVoc*self.ypz)*Cy*Czp,
                 (self.lambda_y_GeVoc*self.lambda_z_GeVoc*self.pypz)*Cyp*Czp,
                 (self.sigma_z_m*self.lambda_z_GeVoc*self.zpz)*Cz*Czp,
                 self.lambda_z_GeVoc**2 * Czp**2,
                 self.offset_x_m * Cx,
                 self.offset_px * Cxp,
                 self.offset_y_m * Cy,
                 self.offset_py * Cyp,
                 self.offset_z * Cz,
                 self.offset_pz * Czp,
                 self.num_zpeaks,
                 dist_num,
                 dist_num]
        return numpy.array(param,'d')
                  
    def get_nparam(self):
        return 30

    def get_dist_type(self):
        if self.transverse:
            return 667
        else:
            return 666

    def get_transverse(self):
        return self.transverse
